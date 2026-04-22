#!/usr/bin/env python3
"""
check_mcmc_progress.py
监控 NG15yr 全量 MCMC 运行进度（T2.4）

用法：
    python3 check_mcmc_progress.py

功能：
  - 读取 chains_full/chain_1.txt 和 pars.txt
  - 报告当前链行数、接受率、log-likelihood
  - 报告 GWB log10_A 和 gamma 的运行估计（中位数 ± 1σ）
  - 报告进度百分比和估计剩余时间
"""

import os, sys, json, time
import numpy as np

# ── 路径 ──────────────────────────────────────────────────────────────
BASE_DIR   = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT_DIR    = os.path.join(BASE_DIR, "results", "T2_NG15yr")
CHAIN_DIR  = os.path.join(OUT_DIR, "chains_full")
CHAIN_FILE = os.path.join(CHAIN_DIR, "chain_1.txt")
PARS_FILE  = os.path.join(CHAIN_DIR, "pars.txt")
LOG_FILE   = os.path.join(OUT_DIR, "mcmc_full_log.txt")

N_STEPS_TARGET = 50_000
THIN           = 10      # PTMCMCSampler 默认每 10 步写一行

def fmt_time(seconds):
    """将秒数格式化为 HH:MM:SS"""
    h, rem = divmod(int(seconds), 3600)
    m, s   = divmod(rem, 60)
    return f"{h:02d}:{m:02d}:{s:02d}"

def read_jump_rates():
    """读取各跳跃提议的接受率"""
    jump_files = {
        "SCAM":     os.path.join(CHAIN_DIR, "covarianceJumpProposalSCAM_jump.txt"),
        "AM":       os.path.join(CHAIN_DIR, "covarianceJumpProposalAM_jump.txt"),
        "EmpDistr": os.path.join(CHAIN_DIR, "draw_from_empirical_distr_jump.txt"),
        "All":      os.path.join(CHAIN_DIR, "jumps.txt"),
    }
    rates = {}
    for name, path in jump_files.items():
        if os.path.exists(path):
            try:
                data = np.loadtxt(path)
                if data.ndim == 0:
                    rates[name] = float(data)
                elif len(data) >= 2:
                    # 格式：[n_accepted, n_proposed] 或直接是接受率
                    if data[-1] > 1.0:
                        # 看起来是累计计数
                        rates[name] = data[0] / max(1, data[1]) if len(data) >= 2 else float(data[-1])
                    else:
                        rates[name] = float(data[-1])
            except Exception:
                pass
    return rates

def main():
    print("=" * 65)
    print("NANOGrav 15yr 全量 MCMC 进度监控 (T2.4)")
    print(f"检查时间：{time.strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 65)

    # ── 检查链目录 ────────────────────────────────────────────────────
    if not os.path.isdir(CHAIN_DIR):
        print(f"\n  链目录不存在：{CHAIN_DIR}")
        print("  全量 MCMC 尚未启动，或目录路径有误。")
        return

    # ── 读取参数名 ────────────────────────────────────────────────────
    param_names = []
    if os.path.exists(PARS_FILE):
        with open(PARS_FILE) as f:
            param_names = [line.strip() for line in f if line.strip()]
        print(f"\n参数总数：{len(param_names)}")
    else:
        print(f"\n  参数名文件不存在：{PARS_FILE}")
        print("  MCMC 可能尚未完成初始化。")

    # ── 读取链文件 ────────────────────────────────────────────────────
    if not os.path.exists(CHAIN_FILE):
        print(f"\n  链文件不存在：{CHAIN_FILE}")
        print("  MCMC 可能尚未开始写入数据。")

        # 尝试显示日志末尾
        if os.path.exists(LOG_FILE):
            print(f"\n--- 日志末尾（{LOG_FILE}）---")
            with open(LOG_FILE) as f:
                lines = f.readlines()
            for line in lines[-20:]:
                print(line, end="")
        return

    # 读取链
    try:
        chain = np.loadtxt(CHAIN_FILE)
    except Exception as e:
        print(f"\n  读取链文件失败：{e}")
        return

    if chain.ndim == 1:
        chain = chain.reshape(1, -1)

    n_rows   = len(chain)
    n_cols   = chain.shape[1] if chain.ndim == 2 else 0
    n_steps_done = n_rows * THIN  # 估计已完成的步数

    print(f"\n链文件大小：{os.path.getsize(CHAIN_FILE) / 1024:.1f} KB")
    print(f"链行数（已写入）：{n_rows:,}")
    print(f"链列数（参数数）：{n_cols}")
    print(f"估计已完成步数：~{n_steps_done:,}  (假设 thin={THIN})")
    print(f"目标步数：{N_STEPS_TARGET:,}")

    pct = 100.0 * n_steps_done / N_STEPS_TARGET
    bar_len = 40
    filled = int(bar_len * pct / 100)
    bar = "#" * filled + "-" * (bar_len - filled)
    print(f"\n进度：[{bar}] {pct:.1f}%")

    # ── 文件修改时间估算速度 ──────────────────────────────────────────
    mtime = os.path.getmtime(CHAIN_FILE)
    age_s = time.time() - mtime
    print(f"\n链文件最后修改：{fmt_time(int(age_s))} 前")

    if os.path.exists(LOG_FILE):
        log_mtime = os.path.getmtime(LOG_FILE)
        log_age_s = time.time() - log_mtime
        print(f"日志文件最后修改：{fmt_time(int(log_age_s))} 前")

        # 从日志中找速度（"步/s" 字段）
        try:
            with open(LOG_FILE) as f:
                lines = f.readlines()
            speed_lines = [l for l in lines if "步/s" in l]
            if speed_lines:
                last_speed_line = speed_lines[-1]
                # 解析格式：速度：X.X 步/s
                import re
                m = re.search(r"速度：([\d.]+)\s*步/s", last_speed_line)
                if m:
                    steps_per_s = float(m.group(1))
                    remaining_steps = max(0, N_STEPS_TARGET - n_steps_done)
                    eta_s = remaining_steps / max(1e-9, steps_per_s)
                    print(f"\n当前采样速度（来自日志）：{steps_per_s:.2f} 步/s")
                    print(f"预计剩余时间：{fmt_time(eta_s)}")
                    print(f"预计完成时刻：{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time() + eta_s))}")
        except Exception:
            pass

    # ── log-likelihood（最后列 - 2）────────────────────────────────
    # PTMCMCSampler chain 格式：[params..., log10_post, log10_like, accept_rate, ...]
    if n_rows > 5 and n_cols > len(param_names):
        # 最后几列是 sampler 元数据
        lnL_col  = len(param_names)      # lnpost 或 lnlike 在 param 列之后
        recent   = chain[-min(500, n_rows):, :]
        try:
            lnpost_vals = recent[:, lnL_col]
            print(f"\nlog-posterior（最近 {len(recent)} 行）：")
            print(f"  最新值：{lnpost_vals[-1]:.2f}")
            print(f"  均值：  {np.mean(lnpost_vals):.2f}")
            print(f"  最大值：{np.max(lnpost_vals):.2f}")
        except Exception as e:
            print(f"\n（无法提取 log-posterior：{e}）")

    # ── 接受率 ────────────────────────────────────────────────────────
    jump_rates = read_jump_rates()
    if jump_rates:
        print(f"\n跳跃接受率：")
        for k, v in jump_rates.items():
            print(f"  {k:12s}：{v:.4f}  ({v*100:.1f}%)")
    else:
        print("\n（跳跃接受率文件尚未生成）")

    # ── GWB 参数估计 ──────────────────────────────────────────────────
    if len(param_names) > 0 and n_rows >= 20:
        gwb_A_idx   = [i for i, n in enumerate(param_names)
                       if "gw" in n.lower() and "log10_a" in n.lower()]
        gwb_gam_idx = [i for i, n in enumerate(param_names)
                       if "gw" in n.lower() and "gamma" in n.lower()]

        burnin = max(5, n_rows // 4)
        post_burn = chain[burnin:, :]

        print(f"\nGWB 参数估计（burn-in = {burnin} 行，后验 {len(post_burn)} 行）：")

        if gwb_A_idx:
            col = gwb_A_idx[0]
            samp = post_burn[:, col]
            med  = np.median(samp)
            std  = np.std(samp)
            p16, p84 = np.percentile(samp, [16, 84])
            print(f"  log10_A_gw  = {med:.4f}  (±1σ: {std:.4f})")
            print(f"              68% CI: [{p16:.4f}, {p84:.4f}]")
            print(f"              参考值：-14.5 ± 0.5 (NG15yr CURN)")

        if gwb_gam_idx:
            col = gwb_gam_idx[0]
            samp = post_burn[:, col]
            med  = np.median(samp)
            std  = np.std(samp)
            p16, p84 = np.percentile(samp, [16, 84])
            print(f"  gamma_gw    = {med:.4f}  (±1σ: {std:.4f})")
            print(f"              68% CI: [{p16:.4f}, {p84:.4f}]")
            print(f"              SMBHB 预期：{13/3:.4f}")
    elif n_rows < 20:
        print(f"\n  链行数 ({n_rows}) 不足 20 行，暂不计算 GWB 参数统计")

    # ── 最近日志 ──────────────────────────────────────────────────────
    if os.path.exists(LOG_FILE):
        print(f"\n--- 日志末尾 20 行（{LOG_FILE}）---")
        with open(LOG_FILE) as f:
            lines = f.readlines()
        for line in lines[-20:]:
            print(line, end="")

    # ── 中间摘要 ──────────────────────────────────────────────────────
    summary_files = sorted([
        f for f in os.listdir(OUT_DIR)
        if f.startswith("T2_summary_full_step") and f.endswith(".json")
    ])
    if summary_files:
        latest_summary = os.path.join(OUT_DIR, summary_files[-1])
        print(f"\n最新中间摘要：{summary_files[-1]}")
        with open(latest_summary) as f:
            s = json.load(f)
        for k, v in s.items():
            if k not in ("pulsar_names",):
                print(f"  {k}: {v}")

    print("\n" + "=" * 65)
    print("监控完成。再次运行此脚本以更新进度。")
    print("=" * 65)


if __name__ == "__main__":
    main()
