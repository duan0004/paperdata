# Plan D Pilot — Decision Memo (2026-04-21)

**Status**: INVESTIGATE — neither hard gate passed, but the result is
**scientifically informative** and reframes (does not kill) Plan D.

## Numerical findings

From `hd_30f_fs.core` post-burn, 9,188 samples, 14 low-k bins:

| metric | value | interpretation |
|---|---|---|
| max \|off-diagonal correlation\| | **0.227** | modest, not dominant |
| mean \|off-diagonal correlation\| | **0.030** | very small |
| #\|corr\|≥0.1 pairs (i<j) | 7 / 91 | sparse |
| top eigenvalue of correlation matrix | **1.483** | 48% above unit |
| MP noise upper λ₊ (p=14, N=9188) | 1.080 | |
| #signal modes (λ > 1.2·λ₊) | **1** | one mode, not ≥3 |
| effective rank r_eff | 13.56 / 14 | near-isotropic |
| Gaussianity (n_bins with \|skew\|<0.5, \|kurt-3\|<2) | 8 / 14 | 6 high-k bins rail-bound |

Low-k subspace scan (bins 0..n−1 for n ∈ {4,6,8,10,12,14}): **max \|corr\| and
λ_max saturate at n = 6**. The single enhancement mode lives in bins 0–5.
Adding more bins only dilutes the mean correlation.

## What this means for Plan D

### The CAR-closes-lnB-gap narrative is **falsified**.

Rough budget for lnZ shift from restoring bin covariance:
Δ ln Z_CAR ≈ ½ log det(Σ/diag(Σ)) = ½ Σ log(λ_k) ≈ **−0.4 nat** (full 14-bin)
or up to ~1 nat for narrow-feature templates that weight the one
enhanced mode.

The observed 2-dex (≈5 nat) absolute-Bayes-factor gap between our
refit and NG15yr Table 3 **cannot** be explained by bin-independence
bias. It must therefore come from:

1. **Template approximation** (semi-analytic vs. PTArcade full integration)
   — most likely dominant
2. **Fixed-γ SMBHB baseline** (quantified: −0.76 nat, so not dominant either)
3. Residual from the ceffyl-KDE log10_ρ → template conversion itself

This **contradicts §IV.C's current claim** that the bin-independence
bias is "plausibly a larger source of the absolute-Bayes-factor gap
than the template approximation itself." §IV.C must be edited.

### The ceffyl approximation is nearly exact.

Mean \|off-diag corr\| = **0.03**. For the bulk of the 14-D posterior,
ceffyl is essentially unbiased. This is a **useful community result**
in its own right — no prior work has quantified this.

### Plan D's original premise is weakened, but not dead.

CAR still produces a **new, publishable observable** (mode-1 projection
coefficient as a discriminator) with ~0.5–1 nat of discriminating power
per template — small but nontrivial, and novel.

## Revised verdict

**Hold 2026-05-05 kill-switch; pivot the Plan D narrative, do not abort.**

### Revised PRL story (if we continue Plan D)

Title idea: *"Covariance-Aware Refit of the NANOGrav 15-Year Stochastic
Gravitational-Wave Signal: Ceffyl is Nearly Exact, Template Choice is
Not"*

Key claims (in descending novelty):
1. We quantify the bin-to-bin posterior covariance of the NG15yr HD
   14-bin free spectrum: max \|corr\| = 0.23, mean = 0.03 → ceffyl is
   validated at ≲1-nat accuracy.
2. The 2-dex lnB gap is therefore attributed to template fidelity, not
   likelihood construction. Full PTArcade integration for SMBHB / SIGW
   closes this gap to within X nat.
3. A 12-model unified semi-analytic library, calibrated against
   PTArcade on 3 reference cases, reproduces NG15yr Table 3 Bayes
   factors to Y accuracy.
4. The one discriminating mode in Σ_bin distinguishes peaked vs.
   power-law spectra at Z nat.

This is still PRL-plausible but **leans more heavily on points 2–3**,
which require completing the template library and PTArcade integration.
Point 3 is just "do the original Phase 1 / Phase 2b-2 work".

### Kill-switch re-evaluation at 2026-05-05

- **GO** condition: PTArcade integration for SMBHB + SIGW-Gauss (Phase
  2b-2) closes lnB gap to <1 nat → points 2–3 substantiated → PRL
  remains viable.
- **NO-GO** condition: PTArcade integration still leaves >1.5 nat
  residual → downgrade to Plan C (PRD only). Story becomes:
  "comprehensive template-library + systematics study", no PRL-grade
  novelty.

### Soft Plan D alternative (lnB(N_bins) truncation)

Run the auxiliary pilot — compute lnZ for each of the 12 templates as
a function of N_bins ∈ {8, 10, 12, 14, 16, 18}. If wrong-shape templates
show systematic drift (e.g., SMBHB vs data show lnB increasing with
N_bins because data deviates from power law in noise-dominated bins)
while right-shape templates are flat, this is a new template-
discrimination observable **that's completely orthogonal to CAR**.

Schedule this as a parallel pilot in Week 2.

## Immediate action items (revised)

1. ✅ Document this memo (D1 complete, with revised findings).
2. Update §IV.C of the paper to reflect the measured bin correlation
   strengths — flag that bin-independence bias is **not** the dominant
   source of the lnB gap. (Deferred to when PTArcade numbers land.)
3. **Pivot the remaining D-week tasks**:
   - D2/D3 (template projection onto eigenmodes) — still useful, but
     expected yield is ~0.5-nat discriminator, not a 5-nat one.
   - Add **D2b**: lnB(N_bins) auxiliary pilot over SMBHB / SIGW / CS.
   - D4 (third chain) — still worth doing for R̂ tightening.
   - D5/D6 (templates) — unchanged, these are must-haves regardless.
4. Hold the 2026-05-05 gate, with PTArcade integration as the new go/no-go
   test (rather than eigenmode count).

## Scientific bottom line

We set out to confirm CAR as a 5-nat effect. We found it is a ~0.5-nat
effect. This is useful information — it tells us where the bias budget
actually lives (templates, not likelihoods). The paper trajectory
adjusts; the project does not abort.
