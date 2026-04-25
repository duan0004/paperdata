# Local Experiment Reproducibility Index

Generated: 2026-04-25T10:13:44

This index records the local P1-P4 run products. Commands are intended to be rerun from the project root.

## Commands

### P1

```bash
python3 code/pta_identifiability_ppc.py
python3 code/bayes_factors_ptarcade_ceffyl_density_extended.py
python3 code/prl_evidence_ti_qmc_crosscheck.py
```

### P2

```bash
python3 code/ng15_2mpz_single_model_resource_gate.py  # matrix over models/counts as recorded in results/lss_tomography/
python3 code/ng15_2mpz_hypermodel_smoke_gate.py --tag prod_nside64_namaster_reference_gate_20260424 --model-set pair-low --components-gw 14 --components-rn 30 --max-pulsars 24 --nsteps 100 --clean
python3 code/ng15_2mpz_null_reproduction_gate.py --tag prod_nside64_namaster_reference_gate_20260424 --components-gw 14 --components-rn 30 --epsilon-max 2 --grid-points 5 --include-red-noise --output-suffix eps2_g5_rn30
python3 code/lss_orf_null_tests.py --tag prod_nside64_weight_nrand36 --nside 64 --nulls 128 --seed 20260425
python3 code/lss_tomography_production_gate.py --nside 64 --tag prod_nside64_local_nulls_20260425 --nulls 128 --seed 20260425 --skip-download
```

### P3

```bash
python3 code/pta_full_loader_noise_tests.py
python3 code/manifest_v2_baseline_gate.py --tag 2026-04-25-local-refresh
python3 code/pta_direct_combination_manifest_audit.py --inpta-band NB --tag 2026-04-25-local-refresh
python3 code/pta_direct_combination_manifest_audit.py --inpta-band WB --tag 2026-04-25-local-refresh
python3 code/inpta_loader_failure_diagnostic.py
python3 code/pta_timing_baseline_pilot.py --arrays NG15_TIMING --max-per-array 8 --components-gw 4 --output-suffix ng15_timing_8_gw4
python3 code/pta_timing_baseline_pilot.py --arrays NG15_TIMING EPTA PPTA --max-per-array 4 --components-gw 4 --output-suffix 3array_4each_gw4
python3 code/pta_timing_baseline_pilot.py --arrays NG15_TIMING EPTA PPTA --max-per-array 4 --components-gw 8 --components-rn 5 --include-red-noise --output-suffix 3array_4each_gw8_rn5
```

### P4

```bash
python3 code/prl_decisive_evidence_figure.py
python3 code/prl_bridge_evidence_figure.py
```

## Key Files And Checksums

| file | bytes | sha256 |
|---|---:|---|
| `theory/Local_runnable_followup_experiment_plan_2026-04-25.md` | 14861 | `e701a88655ce5a6705221905a4d10422e9a97f02386691dd4cb5b5845e550c06` |
| `results/T2_NG15yr/bayes_factors/ptarcade_ceffyl_density_extended.json` | 7300 | `4fdf2922111d809743787a6a0074bfa8b871fcbdc06cb5c9fb5283729a0fafff` |
| `results/T2_NG15yr/bayes_factors/prl_evidence_ti_qmc_crosscheck.json` | 41295 | `8ebe348453c3b74fc94884fa549c1d4e07d8e6ac33f38e88fcc13eef2c64caeb` |
| `results/T2_NG15yr/bayes_factors/prl_evidence_ti_qmc_crosscheck.md` | 2030 | `e3ce80db9cdd7cd64b8da67768ef152f616b64f078ce18b5122518b0e295e761` |
| `theory/PRL_evidence_uncertainty_budget_2026-04-25.md` | 3241 | `137b531aa58db575edef4855d633e06e52f3816f26b390e990788d20950dfdc0` |
| `results/lss_tomography/lss_tomography_production_gate_prod_nside64_local_nulls_20260425.json` | 19036 | `8038dde511db0e6a6cabff87d193fcc241451d1acebec2b866e10ff199cfe19f` |
| `results/lss_tomography/lss_tomography_production_gate_prod_nside64_local_nulls_20260425.md` | 4890 | `1da1cd88f7c0cb24d2805c5cda54d865aafcf08f93bcd74ccb1fd521cb66bf0e` |
| `results/lss_tomography/lss_orf_null_tests_prod_nside64_weight_nrand36.md` | 2724 | `a33fe7c9d815130f50e14cb80b2a463c81a241a2e2f8f8cbbc9add1e493db53c` |
| `results/lss_tomography/ng15_2mpz_hypermodel_smoke_gate_pair-low_prod_nside64_namaster_reference_gate_20260424_gw14_rn30_maxpsr24_n100.json` | 3329 | `f8e9b38951d8162c0cc9f34f7f894dadeed8ce23d68209623acf9d00137d55f1` |
| `results/lss_tomography/ng15_2mpz_null_reproduction_gate_prod_nside64_namaster_reference_gate_20260424_eps2_g5_rn30.json` | 15146 | `519ebbc31b5d015118b8beaa998d6c41ddfcb1ce93b04a9091293cada2e5216b` |
| `results/lss_tomography/ng15_2mpz_null_reproduction_gate_prod_nside64_namaster_reference_gate_20260424_eps3_g7_rn30.json` | 24365 | `b17cbd10951795b56c47fd8078f9c305ddc018d154f15c0be583a471a539e8d5` |
| `results/5pta_timing/full_loader_noise_tests.md` | 2207 | `83b342fc1dc856e81c868a3b4fad1b3b5daad3017a38770eb7a991b3c69c4fbf` |
| `results/5pta_timing/full_loader_noise_summary.csv` | 622 | `f03b76af8fcd01a3e4384b3cc13b599077d124ac7f7b49bc55ffed9468bf5b2b` |
| `theory/PTA_manifest_v2_exact_baseline_gate_2026-04-25-local-refresh.md` | 7159 | `7b2906706f335af48eea138133669b58740200a64d8d1ff67875fdc4cd2a08d5` |
| `results/5pta_timing/manifest_v2_exact_baseline_gate_2026-04-25-local-refresh.json` | 9733 | `86ed190663cea5ce907587f3c26e09a356b8f050d99cfb8d39ec90e309e60495` |
| `theory/5PTA_direct_combination_manifest_audit_nb_2026-04-25-local-refresh.md` | 1928 | `2cdc369ba3666fa2832683a986a27c5dbd4ec1da209f11249bd0f054d2b24c06` |
| `theory/5PTA_direct_combination_manifest_audit_wb_2026-04-25-local-refresh.md` | 1928 | `e56bfb1c5fe89062e6eed13b471a95fa2b92dd3e1d84703cf74bb243680fcfee` |
| `results/5pta_timing/inpta_loader_failure_diagnostic.md` | 1611 | `10bdc92697cf75f73399f8dad009cd4376239fc013866e01212a845da586c062` |
| `theory/5PTA_duplicate_policy_worklist_2026-04-25.md` | 9606 | `f2c14a82d5b873af016bdb9b378f436a613ce06223d76a6d669b59975b2021a2` |
| `results/5pta_timing/baseline_powerlaw_pilot_ng15_timing_8_gw4.json` | 4464 | `4b8181cdd3ae6afba62938f952dbf988c20db51a437c44aafb246bd20f362094` |
| `results/5pta_timing/baseline_powerlaw_pilot_3array_4each_gw4.json` | 6120 | `525b96f520f30b3e617fa9628e6e1f89e1f07cff4597436f3f78214c8279e77f` |
| `results/5pta_timing/baseline_powerlaw_pilot_3array_4each_gw8_rn5.json` | 6952 | `f8f36ea3d457201091ff4c7501024a301cd68bb5daa63173a2e60ab42c177e8f` |
| `results/T2_NG15yr/figures/prl_decisive_evidence_figure.pdf` | 19845 | `41777c8a3a4163af89935dbb654791c7bec215779b1077582dd13b5964c490a6` |
| `results/T2_NG15yr/figures/prl_decisive_evidence_figure.png` | 150157 | `06a474e077419ee6638bd85e5462faa2cb8b7a4e82c6233baaa1f577b8cea642` |
| `results/T2_NG15yr/figures/prl_decisive_evidence_figure_data.json` | 4126 | `688ee2f18f5176a405ce23d168d46f29372d69ad92f93cbf27257530a447642a` |
| `results/T2_NG15yr/figures/prl_bridge_evidence_figure.pdf` | 21117 | `ab6b8030e12abd74db21005548d5743628dad627608a4ba1fd03a5b21ee4b59a` |
| `results/T2_NG15yr/figures/prl_bridge_evidence_figure.png` | 91022 | `0f10dea3b9cdb35a819c072e5050a5e1536db857d256a1420d8bfd998b6a2a16` |

## Gate Interpretation

- P1 products can support PRL stability and uncertainty-budget wording.
- P2 products are LSS geometry/resource/mechanics gates only. They do not reproduce the NG15+2MPZ published null and must not be reported as LSS evidence.
- P3 products are manifest, loader, policy, and mechanical timing gates only. They do not constitute timing-level family evidence.
- P4 products are synchronized figures and local archive metadata.

