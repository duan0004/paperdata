# PRL CAR Covariance Null Calibration

**Date**: 2026-04-21 13:51:20  
**Input**: `results/T2_NG15yr/covariance/corr_bin.npy`  
**Null**: independent Gaussian samples with the same sample count and bin count.  
**Simulations**: 5000  

## Observed Statistics

| statistic | observed | null p-value |
|---|---:|---:|
| max abs off-diagonal correlation | 0.227106 | 0.00019996 |
| mean abs off-diagonal correlation | 0.029916 | 0.00019996 |
| lambda max | 1.482627 | 0.00019996 |
| entropy effective rank | 13.781766 | 0.00019996 |
| participation effective rank | 13.555656 | 0.00019996 |
| signal modes above 1.2 MP upper edge | 1 | 0.00019996 |

## Null Percentiles

| statistic | p50 | p95 | p99 | p99.9 | max |
|---|---:|---:|---:|---:|---:|
| max abs off-diagonal correlation | 0.027899 | 0.035972 | 0.040569 | 0.046070 | 0.049159 |
| mean abs off-diagonal correlation | 0.008326 | 0.009421 | 0.009902 | 0.010447 | 0.010999 |
| lambda max | 1.064490 | 1.078452 | 1.084591 | 1.092497 | 1.104194 |
| entropy effective rank | 13.990141 | 13.992341 | 13.993147 | 13.993962 | 13.994566 |
| participation effective rank | 13.980313 | 13.984715 | 13.986305 | 13.987969 | 13.989136 |

## Interpretation

The observed correlation structure is not consistent with pure independent
Gaussian sampling noise: the leading correlation mode is real.  However, the
mean absolute off-diagonal correlation remains small and both effective-rank
definitions stay close to the full 14-bin rank.  This supports the PRL claim that
bin covariance exists but is subdominant compared with template/density
provenance for the leading official-density evidence reproduction.
