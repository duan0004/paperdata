#!/usr/bin/env python3
"""
NG15yr_production_HD.py — T2.4 production MCMC: HD-correlated GWB
==================================================================

相比 NG15yr_powerlaw_full.py 的改动：
  1. 信号模型：CURN → HD-correlated common process（utils.hd_orf）
  2. 步数：50_000 → 2_000_000（可 resume 继续到 5M+）
  3. 初值：随机先验 → 在官方后验中心 (-14.20, 3.25) 附近 ±0.5 抖动
  4. resume=True：支持中断后从 chain_1.txt 续跑
  5. 中间摘要：每 50_000 步保存一次 JSON

参考：
  - 数据：arXiv:2306.16213
  - 噪声：arXiv:2306.16214
  - HD/CURN 比较：arXiv:2306.16213 §VI
"""

import os, sys, json, logging, warnings, time, threading
import numpy as np

logging.getLogger("pint").setLevel(logging.WARNING)
logging.getLogger("pint.observatory").setLevel(logging.WARNING)
warnings.filterwarnings("ignore")

# ── 配置 ─────────────────────────────────────────────────────────────
N_STEPS          = 500_000        # 冷启动链，500k 步 (~7 h)
N_GWBFREQS       = 14
GAMMA_COMMON     = None           # 自由 γ
USE_EMP_DISTR    = True
USE_HD           = True           # ← 启用 HD 空间相关
THIN             = 10
ISAVE            = 10_000         # 进度打印间隔
SUMMARY_ROWS_INT = 5_000          # 每 5000 行（=50k 步）摘要

# 初值中心：冷启动远离 warm chain 后验中心，用于 Gelman-Rubin R̂
# warm chain 结果: log10_A=-14.20, γ=3.25 → 冷启动分散到低振幅+陡谱
LOG10A_INIT      = -15.50         # 比后验低 ~10σ (σ=0.12)
GAMMA_INIT       = 5.00           # 比后验高 ~5σ (σ=0.32)，且远离 13/3
INIT_JITTER      = 1.0            # 大 jitter，进一步去相关

# ── 路径 ──────────────────────────────────────────────────────────────
BASE_DIR    = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR    = os.path.join(BASE_DIR, "data", "NG15yr", "tutorials", "data")
FEATHER_DIR = os.path.join(DATA_DIR, "feathers")
WN_DICT     = os.path.join(DATA_DIR, "15yr_wn_dict.json")
EMP_DISTR   = os.path.join(DATA_DIR, "15yr_emp_distr.json")
OUT_DIR     = os.path.join(BASE_DIR, "results", "T2_NG15yr")
CHAIN_DIR   = os.path.join(OUT_DIR, "chains_production_HD_cold")

os.makedirs(CHAIN_DIR, exist_ok=True)

T_START = time.time()
def elapsed():
    s = time.time() - T_START
    h, rem = divmod(int(s), 3600); m, sec = divmod(rem, 60)
    return f"{h:02d}:{m:02d}:{sec:02d}"

print("=" * 70)
print("NG15yr_production_HD — T2.4 production MCMC（HD 空间相关）")
print("=" * 70)
print(f"  脉冲星：全部 67 颗")
print(f"  信号：HD-correlated common process（{N_GWBFREQS} 频率分量）")
print(f"  目标步数：{N_STEPS:,}（可 resume）")
print(f"  γ：{'自由拟合' if GAMMA_COMMON is None else GAMMA_COMMON}")
print(f"  初值中心：log10_A={LOG10A_INIT}, γ={GAMMA_INIT} ± {INIT_JITTER}")
print(f"  链目录：{CHAIN_DIR}")
print(f"  启动时间：{time.strftime('%Y-%m-%d %H:%M:%S')}")
print()
sys.stdout.flush()

# 检查路径
for path, name in [(FEATHER_DIR, "feather 目录"), (WN_DICT, "白噪声字典"),
                   (EMP_DISTR, "经验分布")]:
    if not os.path.exists(path):
        print(f"错误：{name} 不存在：{path}"); sys.exit(1)

# ── 加载脉冲星 ────────────────────────────────────────────────────────
print(f"[T={elapsed()}] 加载脉冲星 ..."); sys.stdout.flush()
from enterprise_extensions.load_feathers import load_feathers_from_folder
psrs = load_feathers_from_folder(FEATHER_DIR)
print(f"  [T={elapsed()}] {len(psrs)} 颗脉冲星已加载"); sys.stdout.flush()

with open(WN_DICT) as f:
    wn_params = json.load(f)
print(f"  白噪声参数：{len(wn_params)}")

tmin = np.min([p.toas.min() for p in psrs])
tmax = np.max([p.toas.max() for p in psrs])
Tspan = tmax - tmin
print(f"  Tspan = {Tspan / 3.15576e7:.2f} 年"); sys.stdout.flush()

# ── PTA 模型 ─────────────────────────────────────────────────────────
print(f"\n[T={elapsed()}] 建模（HD = {USE_HD}）..."); sys.stdout.flush()
from enterprise.signals import parameter, utils, signal_base, selections
from enterprise.signals import white_signals, gp_signals

selection = selections.Selection(selections.by_backend)
mn = white_signals.MeasurementNoise(
        efac=parameter.Constant(),
        log10_t2equad=parameter.Constant(),
        selection=selection)
ec = white_signals.EcorrKernelNoise(
        log10_ecorr=parameter.Constant(), selection=selection)

# Per-pulsar red noise
log10_A_rn = parameter.Uniform(-20, -11)
gamma_rn   = parameter.Uniform(0, 7)
pl_rn = utils.powerlaw(log10_A=log10_A_rn, gamma=gamma_rn)
rn = gp_signals.FourierBasisGP(spectrum=pl_rn, components=30, Tspan=Tspan)

# Common GWB process (HD or CURN)
log10_A_gw = parameter.Uniform(-18, -14)("log10_A_gw")
if GAMMA_COMMON is None:
    gamma_gw = parameter.Uniform(0, 7)("gamma_gw")
else:
    gamma_gw = parameter.Constant(GAMMA_COMMON)("gamma_gw")
cpl = utils.powerlaw(log10_A=log10_A_gw, gamma=gamma_gw)

if USE_HD:
    orf = utils.hd_orf()
    gw  = gp_signals.FourierBasisCommonGP(spectrum=cpl, orf=orf,
              components=N_GWBFREQS, Tspan=Tspan, name="gw")
else:
    gw  = gp_signals.FourierBasisGP(spectrum=cpl, components=N_GWBFREQS,
              Tspan=Tspan, name="gw")

tm  = gp_signals.MarginalizingTimingModel(use_svd=True)
s   = tm + mn + ec + rn + gw
pta = signal_base.PTA([s(p) for p in psrs])
pta.set_default_params(wn_params)

ndim = len(pta.params)
print(f"  PTA ndim = {ndim}")
sys.stdout.flush()

# 保存参数名
with open(os.path.join(CHAIN_DIR, "pars.txt"), "w") as f:
    for n in pta.param_names: f.write(n + "\n")

# ── 初值（围绕官方后验中心）──────────────────────────────────────────
print(f"\n[T={elapsed()}] 构造初值（围绕 log10_A={LOG10A_INIT}, γ={GAMMA_INIT}）...")
def make_x0():
    x0 = np.array([p.sample() for p in pta.params]).flatten()
    for i, name in enumerate(pta.param_names):
        if name == "log10_A_gw":
            x0[i] = LOG10A_INIT + INIT_JITTER * (np.random.rand() - 0.5)
        elif name == "gamma_gw":
            x0[i] = GAMMA_INIT + INIT_JITTER * (np.random.rand() - 0.5)
    return x0

x0 = make_x0()
lnL = pta.get_lnlikelihood(x0); lnPr = pta.get_lnprior(x0)
print(f"  init lnL = {lnL:.4f}, lnPr = {lnPr:.4f}")
if not np.isfinite(lnL):
    for k in range(50):
        x0 = make_x0()
        lnL = pta.get_lnlikelihood(x0)
        if np.isfinite(lnL):
            print(f"  resampled (try {k+1}): lnL = {lnL:.4f}"); break
    else:
        print("  错误：无法找到有限初值"); sys.exit(1)
sys.stdout.flush()

# ── 经验分布 ─────────────────────────────────────────────────────────
emp_distr_list = None
if USE_EMP_DISTR:
    try:
        TUT = os.path.dirname(DATA_DIR)
        sys.path.insert(0, TUT)
        from extra_functions import EmpiricalDistribution2D
        with open(EMP_DISTR) as f: raw = json.load(f)
        emp_distr_list = [EmpiricalDistribution2D(v) for v in raw.values()]
        print(f"  经验分布：{len(emp_distr_list)} 个")
    except Exception as e:
        print(f"  经验分布加载失败：{e}")

# ── 后台摘要线程 ─────────────────────────────────────────────────────
chain_file = os.path.join(CHAIN_DIR, "chain_1.txt")
_stop = threading.Event()

def watcher():
    last = 0
    names = pta.param_names
    iA = next((i for i,n in enumerate(names) if n=="log10_A_gw"), None)
    iG = next((i for i,n in enumerate(names) if n=="gamma_gw"), None)
    while not _stop.is_set():
        time.sleep(60)
        if not os.path.exists(chain_file): continue
        try:
            ch = np.loadtxt(chain_file)
            if ch.ndim < 2 or len(ch) < 50: continue
            n = len(ch)
            if n - last < SUMMARY_ROWS_INT: continue
            last = n
            step = n * THIN
            burn = max(50, n // 4)
            post = ch[burn:, :]
            summary = {
                "task": "T2.4-production-HD",
                "model": "HD" if USE_HD else "CURN",
                "n_pulsars": len(psrs),
                "n_steps_target": N_STEPS,
                "n_steps_estimated": int(step),
                "chain_rows": int(n),
                "burnin_rows": int(burn),
                "elapsed": elapsed(),
                "pct_complete": f"{100.0 * step / N_STEPS:.1f}%",
                "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
            }
            if iA is not None and len(post) > 10:
                a = post[:, iA]
                summary["log10_A_gw"] = {
                    "median": float(np.median(a)), "std": float(np.std(a)),
                    "p16": float(np.percentile(a, 16)),
                    "p84": float(np.percentile(a, 84))}
            if iG is not None and GAMMA_COMMON is None and len(post) > 10:
                g = post[:, iG]
                summary["gamma_gw"] = {
                    "median": float(np.median(g)), "std": float(np.std(g)),
                    "p16": float(np.percentile(g, 16)),
                    "p84": float(np.percentile(g, 84))}
            tag = f"step{int(step):08d}"
            sp = os.path.join(OUT_DIR, f"T2_summary_prodHD_cold_{tag}.json")
            with open(sp, "w") as f: json.dump(summary, f, indent=2, ensure_ascii=False)
            print(f"\n[T={elapsed()}] 摘要 @ {step:,} 步：")
            if "log10_A_gw" in summary:
                a = summary["log10_A_gw"]
                print(f"  log10_A_gw = {a['median']:.4f} ± {a['std']:.4f}")
            if "gamma_gw" in summary:
                g = summary["gamma_gw"]
                print(f"  gamma_gw   = {g['median']:.4f} ± {g['std']:.4f}")
            sys.stdout.flush()
        except Exception as e:
            print(f"\n[watcher 错误] {e}"); sys.stdout.flush()

threading.Thread(target=watcher, daemon=True).start()
print(f"[T={elapsed()}] 后台 watcher 已启动")

# ── MCMC ─────────────────────────────────────────────────────────────
# 检查是否 resume
existing = os.path.exists(chain_file) and os.path.getsize(chain_file) > 1000
RESUME   = existing
print(f"\n[T={elapsed()}] {'RESUME 已存链' if RESUME else '全新启动'}，目标 {N_STEPS:,} 步")
print(f"  预计速率 ~{60 if USE_HD else 170} 步/秒 → ETA {N_STEPS/(60 if USE_HD else 170)/3600:.1f} 小时")
sys.stdout.flush()

from PTMCMCSampler.PTMCMCSampler import PTSampler as ptmcmc
cov = np.diag(np.ones(ndim) * 0.01**2)

sampler = ptmcmc(ndim, pta.get_lnlikelihood, pta.get_lnprior, cov,
                 outDir=CHAIN_DIR, resume=RESUME)

if emp_distr_list is not None:
    try:
        from enterprise_extensions import sampler as ee_sampler
        jp = ee_sampler.JumpProposal(pta, empirical_distr=emp_distr_list)
        sampler.addProposalToCycle(jp.draw_from_empirical_distr, 10)
        print(f"  经验分布提议已注入（weight=10）")
    except Exception as e:
        print(f"  经验分布提议失败：{e}")

t0 = time.time()
sampler.sample(x0, N_STEPS, SCAMweight=30, AMweight=15, DEweight=50,
               thin=THIN, isave=ISAVE)
dt = time.time() - t0
print(f"\n[T={elapsed()}] sample() 返回，用时 {dt/3600:.2f} 小时，"
      f"速率 {N_STEPS/max(1,dt):.2f} 步/秒")
sys.stdout.flush()

_stop.set()

# ── 最终摘要 ─────────────────────────────────────────────────────────
chain = np.loadtxt(chain_file)
print(f"\nchain shape = {chain.shape}")
burn = max(100, len(chain) // 4)
names = pta.param_names
iA = next((i for i,n in enumerate(names) if n=="log10_A_gw"), None)
iG = next((i for i,n in enumerate(names) if n=="gamma_gw"), None)

final = {
    "task": "T2.4-production-HD",
    "model": "HD" if USE_HD else "CURN",
    "n_pulsars": len(psrs),
    "n_steps": N_STEPS,
    "Tspan_yr": float(Tspan / 3.15576e7),
    "chain_shape": list(chain.shape),
    "burnin_rows": int(burn),
    "total_elapsed": elapsed(),
    "completed_at": time.strftime("%Y-%m-%d %H:%M:%S"),
}
if iA is not None:
    a = chain[burn:, iA]
    final["log10_A_gw"] = {"median": float(np.median(a)), "std": float(np.std(a)),
                          "p16": float(np.percentile(a, 16)),
                          "p84": float(np.percentile(a, 84))}
    print(f"  log10_A_gw = {final['log10_A_gw']['median']:.4f} ± "
          f"{final['log10_A_gw']['std']:.4f}（官方：-14.20 ± 0.13）")
if iG is not None and GAMMA_COMMON is None:
    g = chain[burn:, iG]
    final["gamma_gw"] = {"median": float(np.median(g)), "std": float(np.std(g)),
                         "p16": float(np.percentile(g, 16)),
                         "p84": float(np.percentile(g, 84))}
    print(f"  gamma_gw   = {final['gamma_gw']['median']:.4f} ± "
          f"{final['gamma_gw']['std']:.4f}（官方：3.25 ± 0.35）")

fp = os.path.join(OUT_DIR, "T2_summary_prodHD_cold_final.json")
with open(fp, "w") as f: json.dump(final, f, indent=2, ensure_ascii=False)
print(f"\n最终摘要：{fp}")
print("=" * 70)
