# Manifest v2 and Exact Timing-Level Baseline Gate

Generated: 2026-04-25T10:00:22

This file separates local ingestion status from public availability.  A `BLOCKED`
status below means not yet present locally or not yet harmonized with this
project's exact reproduction manifest; it does not mean the public input is absent.

## Public Five-PTA Reference Target

- Reference source: `arXiv:2512.08666 source lines 338-350, 359-428`.
- Direct-combination priority: `NG15 > EPTA > PPTA > MPTA; InPTA omitted from priority order because its pulsars are observed by other PTAs in the reference analysis.`.
- White-noise policy summary: `NG15/PPTA/MPTA and seven EPTA backends use EFAC/EQUAD/ECORR; InPTA and the remaining EPTA backends use EFAC/EQUAD.`.

| array | reference count | local count |
|---|---:|---:|
| `EPTA` | 25 | 25 |
| `InPTA_DR1` | 14 | 14 |
| `MPTA_4p5yr` | 83 | 83 |
| `NG15` | 68 | 68 |
| `PPTA` | 32 | 32 |

## Local Inventory

| product | count | status |
|---|---:|---|
| `NG15` | 67 | `PRESENT_LOCALLY` |
| `NG15_TIMING` | 68 | `PRESENT_LOCALLY` |
| `EPTA` | 25 | `PRESENT_LOCALLY` |
| `PPTA` | 32 | `PRESENT_LOCALLY` |
| `InPTA_DR1_NB` | 14 | `PRESENT_LOCALLY` |
| `InPTA_DR1_WB` | 14 | `PRESENT_LOCALLY` |
| `MPTA_4p5yr` | 83 | `PRESENT_LOCALLY` |

## Reduced Development Namespace

The current three-array timing diagnostics use an `NG15-priority deduplicated
reduced 3-array development namespace`, not a full public five-PTA direct
combination.

- Retained pulsars: `87`.
- Skipped duplicate occurrences: `37`.
- Manifest CSV: `results/5pta_timing/manifest_v2_ng15_priority_dedup_3array_development_2026-04-25-local-refresh.csv`.
- Unique-only retained pulsars: `60`.
- Unique-only skipped duplicate occurrences: `64`.
- Unique-only CSV: `results/5pta_timing/manifest_v2_unique_only_3array_development_2026-04-25-local-refresh.csv`.

## Duplicate Sources Across NG15/EPTA/PPTA

| pulsar | arrays | occurrences |
|---|---|---:|
| `J0030+0451` | `EPTA, NG15, PPTA` | 3 |
| `J0437-4715` | `NG15, PPTA` | 2 |
| `J0613-0200` | `EPTA, NG15, PPTA` | 3 |
| `J0900-3144` | `EPTA, PPTA` | 2 |
| `J1012+5307` | `EPTA, NG15` | 2 |
| `J1022+1001` | `EPTA, NG15, PPTA` | 3 |
| `J1024-0719` | `EPTA, NG15, PPTA` | 3 |
| `J1455-3330` | `EPTA, NG15` | 2 |
| `J1600-3053` | `EPTA, NG15, PPTA` | 3 |
| `J1640+2224` | `EPTA, NG15` | 2 |
| `J1643-1224` | `NG15, PPTA` | 2 |
| `J1713+0747` | `EPTA, NG15, PPTA` | 3 |
| `J1730-2304` | `EPTA, NG15, PPTA` | 3 |
| `J1738+0333` | `EPTA, NG15` | 2 |
| `J1741+1351` | `NG15, PPTA` | 2 |
| `J1744-1134` | `EPTA, NG15, PPTA` | 3 |
| `J1751-2857` | `EPTA, NG15` | 2 |
| `J1832-0836` | `NG15, PPTA` | 2 |
| `J1843-1113` | `EPTA, NG15` | 2 |
| `J1857+0943` | `EPTA, PPTA` | 2 |
| `J1909-3744` | `EPTA, NG15, PPTA` | 3 |
| `J1910+1256` | `EPTA, NG15` | 2 |
| `J1911+1347` | `EPTA, NG15` | 2 |
| `J1918-0642` | `EPTA, NG15` | 2 |
| `J2124-3358` | `EPTA, NG15, PPTA` | 3 |
| `J2145-0750` | `NG15, PPTA` | 2 |
| `J2322+2057` | `EPTA, NG15` | 2 |

## Local Count Discrepancy Audit

| product | status | detail | missing/incomplete |
|---|---|---|---|
| `NG15_feathers` | `LOCAL_COUNT_67_REFERENCE_68` | Local feather directory contains 67 files; ENTERPRISE loader returns 67 pulsars. | `J0614-3329` |
| `NG15_full_timing_reference_tree` | `REFERENCE_TREE_HAS_68_NAMES_LOADER_SMOKE_PASS` | Reference-project NANOGrav 15-year timing residual tree has 68 unique pulsars. The feather-missing target J0614-3329 has par=True and tim=True. NG15_TIMING loader smoke status=PASS, available=68, loaded=68, failures=0. | `official NG15 noise-policy wiring remains required before timing-level evidence` |
| `PPTA_DR3_toas_and_parameters_all` | `METADATA_HAS_32_PAR_31_TIM` | Metadata lists 32 standard .par files, 32 single-pulsar-noise .par files, and 31 .tim files. | `J1741+1351` |
| `PPTA_DR3_github_analysis_codes_data_all` | `GITHUB_SOURCE_HAS_32_PAR_32_TIM` | GitHub source has 32 standard .par files and 32 .tim files at commit fdbe6eb1c86d4c6cf2f1f518711e44ad1a9fd3fa. | `` |
| `MPTA_4p5yr_legacy_DOI_partim` | `LOCAL_ARCHIVE_AND_EXTRACTED_TREE_HAVE_74_PAIRS` | legacy local archive previously used by this project: extracted tree has 74 .par names and 74 .tim names (74 matched pairs); archive has 74 .par and 74 .tim entries. | `` |
| `MPTA_4p5yr_DataCentral_partim_20260425` | `LOCAL_ARCHIVE_AND_EXTRACTED_TREE_HAVE_83_PAIRS` | refreshed DataCentral partim download: extracted tree has 83 .par names and 83 .tim names (83 matched pairs); archive has 83 .par and 83 .tim entries. | `` |
| `MPTA_4p5yr_Anisotropy_supplement_zip` | `NO_TIMING_INPUTS_FOUND` | Supplement zip contains 10 entries, including 9 MP4 files, and 0 timing-like or manifest-like entries. | `does not provide timing inputs or duplicate/noise-policy metadata` |
| `MPTA_4p5yr_archives_tar` | `RAW_DLY_ARCHIVE_NO_PAR_TIM` | archives.tar.gz contains 10100 entries, 86 directories, 85 pulsar directories under data_august23_32ch, 10014 .dly files, and 0 .par / 0 .tim entries. | `raw .dly archive does not provide baseline-ready .par/.tim inputs` |

## Exact Baseline Gate

| gate | status | detail |
|---|---|---|
| NG15 timing products | `PRESENT_LOCALLY` | NANOGrav 15-year reference par/tim target count 68; reference target count 68. The separate feather development loader has 67 pulsars and must not be mixed into exact timing-level evidence. |
| EPTA DR2 timing products | `PRESENT_LOCALLY` | EPTA DR2 par/tim count 25; reference target count 25. |
| PPTA DR3 timing products | `PRESENT_LOCALLY` | PPTA DR3 local par/tim pair count 32; reference target count 32.  Count now matches when using the GitHub reference source. |
| InPTA DR1 exact target inputs | `PRESENT_BUT_POLICY_BLOCKED` | Local tree contains 14 NB and 14 WB InPTA DR1 par/tim pairs; reference target count is 14.  Read-only diagnostics show NB+tempo2 loads 13/14 and fails on J0751+1807_NB, WB+PINT loads 13/14 and fails on J0613-0200_WB, NB+PINT fails all 14 on TCB/TDB handling, and WB+tempo2 fails all 14 with empty-TOA loader errors.  Exact baseline reproduction therefore remains blocked until the published NB/WB, timing-package, and TCB-to-TDB preprocessing policy is identified. |
| MPTA exact target inputs | `PRESENT_BUT_POLICY_BLOCKED` | Local MPTA par/tim pair count is 83; reference target count is 83. Production baseline use also requires exact noise-policy and duplicate-combination harmonization. |
| duplicate pulsar handling | `BLOCKED_BY_POLICY_HARMONIZATION` | Current reduced 3-array diagnostics use priority deduplication or unique-only dropping.  A production baseline must reproduce the published direct-combination policy for duplicated pulsars before family evidence. |
| timing-level baseline reproduction | `NOT_STARTED_FOR_EXACT_TARGET` | No exact public five-PTA baseline likelihood reproduction has passed locally; family evidence remains on hold. |

## Decision

- Do not extend the current reduced free-red chains as a production route.
- First reproduce the selected public timing-level baseline with the same duplicate-handling and noise policy.
- Only after that gate passes should minimal family evidence be run on timing-level likelihoods.
