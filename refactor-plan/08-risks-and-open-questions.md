GEMS3K commit: 54b5f819ba83589daa9cd9aa0a9ff699adb28b1d (branch develop_optima, 2026-11-03)

# 08 — Risks and open questions

## 1. Blocker found in this session (updated after access was granted)

**gems-benchmark is now accessible, but the freeze still cannot run outside the owner's machine.** `freeze.sh`
calls `mode_compare` from the owner's local `GEMS3K/debug-optima-vs-reaktoro/` (in no repository), and
`freeze.sh`, `build-gems3k.sh` and `recheck.py` hard-code `/home/dmiron/...` paths; the recorded standard freeze is
in the git-ignored `.claude/docs/freeze/`. What ran here (03 §7): trace self-check OK, CTest 92/96, `ci.sh` 114 of
132 rows exact (18 rows differ for environment reasons). Every Stage A step depends on the freeze (05 §0), so
**A0 must run on the owner's machine first, or `mode_compare` must be published (Q21)**.

## 2. Risks

| # | Risk | Likelihood / impact | Mitigation |
|---|---|---|---|
| R1 | A "pure move" changes FP results through inlining or a changed translation unit | medium / high | A-rule 3: keep hot callers and callees together; instruction counts + freeze per step; `objdump` identity for A10a |
| R2 | Removing the hard-coded `/usr/local` paths picks up another Eigen/spdlog (Eigen 5 rounds the pivoted QR differently) | low / high | A1 fails on Eigen ≠ 3.4.x unless an explicit option is set; compile-flag diff script |
| R3 | A thread-local flag is read on a path outside the call chain (A6) | medium / high | per-flag read/write list in each PR; debug mirror asserts for one release; DECIDE golden test |
| R4 | GEMSGUI breaks: it subclasses TMultiBase, reads MULTI, defines `BASE_PARAM::read/write`, binds fields by pointer | medium / high | 02 §2.4–§2.5 kept; `static_assert` on the 38-member prefix; GEMSGUI build in GEMS3K CI (owner, Q13) |
| R5 | Stage A tempts "small fixes" that change answers (FitVar copy, OK test, nDC loop, model caches) | high / medium | all are listed as B1/B2 with owner decisions; Stage A copies today's behaviour exactly |
| R6 | Freezes taken on different machines differ (glibc picks `__ieee754_log_fma` by CPU; aarch64/Apple clang contract FMA by default) | medium / high | one machine class for all Stage A freezes; B10/B12 make it explicit. **Measured (03 §7):** disabling glibc's FMA/AVX dispatch removes all 5 CI findings of `f_GEOTHERM` (incl. ITG 687 → 688); gcc 14.3 vs 16.2 changes nothing |
| R7 | The activity call is not idempotent (lnGmo, F0 smoothing via FitVar[3]); any change of call count changes results | certain / high | Stage A never changes the number or order of activity calls; B4 verifies bit-identity on the Optima path first |
| R8 | The Optima fork's feature macros (`OPTIMA_LU_PARTIAL_PIVOT`, `…SNAP_LEAVING`, `…STALL_ESCAPE`, `…REJECT_WORSE`) silently do nothing against upstream Optima | medium / medium | print the macros in EFF (A13); `gems3k_build_features()` (A1) |
| R9 | `NativeOptimaFinish/Rescue` reset the iteration budget to −1 instead of restoring it; nesting a cap around them would lose the outer budget | low today / medium after refactor | A6b keeps the behaviour and adds a test; B1 restores the previous value |
| R10 | `TNodeArray` warm results depend on the OpenMP thread count | certain with OpenMP / medium | A11 names the policy; B11 changes it (scored) |
| R11 | Plan scale: 25–35 Stage A PRs, each with two freeze pairs (~hours each) | certain / schedule | batch small S steps; freeze jobs in parallel (`FREEZE_JOBS`); no build during freezes |
| R12 | The substitute corpus under-represents the owner's large systems; B priorities may shift | medium / medium | B0 and each B step re-measure on the owner's corpus before merging |
| R13 | dn/db users (if any) get values from the primary Optima attempt, not the accepted answer | unknown / medium | B1; document now (A14) |
| R14 | ThermoFun does not cover every T-P method used in the GEMS databases (known: `f_PitzerTHE` water data) | medium / high (blocks B14) | coverage part of the campaign (10 §3 c); grid kept until the list is empty |
| R15 | ThermoFun direct calculation is slower than grid interpolation for transport with many temperatures | medium / medium | measure (10 §3 d); exact (T, P) memo or on-demand table (12 §5) |
| R16 | Moving the models into ThermoSol/ThermoKin targets changes inlining and so FP results | medium / medium | A-rule 3 checks; B15 scored even when identity is expected |
| R17 | Two copies of the phase-property code (TMultiBase vs TActivity) differ in a detail; merging changes results | medium / high | A16 merges only token-identical code |

## 3. Open questions (each with the default this plan assumes)

| # | Question | Assumed default |
|---|---|---|
| Q1 | Can the planning/implementation sessions get read access to gems-benchmark (freeze tools, corpus, CTest)? | Granted during this session (03 §7). The freeze still needs Q21. |
| Q2 | May the plan files be pushed to a branch `refactor-plan` in GEMS3K? | Not pushed. They are in the untracked folder `refactor-plan/` and in the zip. |
| Q3 | Single declaration of `BASE_PARAM`: X-macro (`.def` file) or a generator script (Python → header)? | X-macro (no build-time tool; works with qmake and MSVC). |
| Q4 | Keep the qmake files (`gems3k.pri`, `node-gem.pro`, `gemcalc.pro`) and the MinGW Makefiles? | Keep in Stage A, generated from the CMake lists by a script; drop in 5.0 if nobody uses them. |
| Q5 | Removal date for deprecated API and internal headers | Deprecate in the release that closes Stage A (4.5); remove in 5.0, not before 2027-12 and not before xgems and GEMSGUI have shipped without them. |
| Q6 | Should GEMSGUI generate its `o_pax*` objects, `vis_od.ini.json` rows and `ResetTail` from the same `.def`? | Yes, as a later GEMSGUI PR after A3 (optional). |
| Q7 | Fix `TNodeArray`'s success test (Optima OKs counted as errors)? | Yes, in B1; GEMSGUI's override checked in the same change. |
| Q8 | May the plan propose changes to the Optima fork (thread-local trace state, EFF of feature macros)? | Yes, as separate PRs in ../optima, owner-approved; none needed for Stage A results. |
| Q9 | `native_cert_dualfree_resolved` (process-wide, resolved once): is "first call in the process" the intended scope? | Yes; keep it with `std::call_once` (A8). |
| Q10 | The probe `GEMS3K_PHSTAB_PROBE` opens its file twice (two statics). One handle changes line order in the probe file. Accept? | Accept (probe output, not freeze). |
| Q11 | FP contract for releases: pin `-ffp-contract=off` on all platforms? | Yes, but only in Stage B (B10), scored; Stage A changes no flags. |
| Q12 | GPU work: start only after a CPU batch experiment shows ≥ 70 % of node-iterations on the common path? | Yes. |
| Q13 | Add GEMSGUI and xgems builds to GEMS3K CI (as consumers)? | Yes for xgems (public); GEMSGUI if CI time allows. |
| Q14 | Keep `TActivity`/`TKinetics`? | **Answered by the owner: yes** — they are the seeds of ThermoSol and ThermoKin (11); not deprecated; the TMultiBase copies go (A16). |
| Q15 | Is the stale-DATABR-on-ERR behaviour (trap 6) relied on by anyone, or may B change it to "NaN-fill on ERR"? | Keep (contract); a B change only with consumer agreement. |
| Q16 | Line endings: pin per file with `.gitattributes` (keep CRLF for `ipm_chemical.cpp`, `s_kinmet.cpp`) or convert both to LF? | Pin as they are (no content change in Stage A); convert in a separate whitespace-only PR later if the owner wants. |
| Q17 | Application class for quality management (de-RSE/GI guideline) | Highest class (library used by other groups and tools). |
| Q18 | Is the dn/db feature (`optima_want_sensitivity`) used by any consumer today? | Not found in xgems or GEMSGUI; treat as internal until a getter is published (07 §3). |
| Q19 | Two environment variables change results (`GEMS3K_OPTIMA_NET_RESUME`, `GEMS3K_FILLOUT_BUDGET`). Move them into `BASE_PARAM`? | Keep the env names; read once into RunConfig and print in EFF (A8); adding `pa_` fields is a later owner decision (BASE_PARAM contract). |
| Q20 | Established terms: the plan keeps "GEMS3K IPM solver", "warm start", "cold start". New terms introduced: RunContext, WarmState, PhaseModel, DiagnosticsSink, SolverWorkspace. OK? | Yes; no established term was changed. |
| Q21 | Publish `mode_compare` (the freeze driver, now only in the owner's local `GEMS3K/debug-optima-vs-reaktoro/`) and the standard freeze, and make the gems-benchmark tool paths configurable? | Yes, as a gems-benchmark PR before A1; until then every freeze runs on the owner's machine. |
| Q22 | Which CTest failures are known? Brief: "44 tests, 4 known failures (2026-10-17)". Here: 96 tests, failing `optima_regression.aop.restart`, `proposed.aop`, `dimreduce.gate` (+ `ci.baseline`, invalid with Optima). | Assume these three are the known Optima-side failures; the owner confirms on the A0 machine. |
| Q23 | Freeze environment: which CPU, glibc, ThermoFun and compiler recorded the CI baseline and the standard freeze? 18 of 132 CI rows differ here for environment reasons. | Record them in every freeze header from A0 on; one environment for all Stage A freezes. |
| Q24 | GEMSGUI's own T-P methods (`s_supcrt*`, `m_dcomp*`, `m_reac*`, `ms_mtparm`, ~9 000 lines): replace with ThermoFun too? | Yes, in a separate GEMSGUI track behind the same comparison (10 §6). (Asked; no answer yet.) |
| Q25 | Where do ThermoSol and ThermoKin live while they are made? | Inside GEMS3K first (`thermosol/`, `thermokin/` targets and namespaces), own repositories (thermohub) when the API and the v2 record format are stable (B17). (Asked; no answer yet.) |
| Q26 | Tolerances of the ThermoFun comparison (10 §3) | \|ΔG0\| ≤ 1 J/mol; relative 1e-6 for V0, S0, Cp0, H0 and water properties; freeze tolerances for equilibria. |
| Q27 | Names | ThermoSol and ThermoKin, as given by the owner. |
| Q28 | Sorption models (TSorpMod stub + sorption code in TMultiBase): in ThermoSol, or a library of their own? | In ThermoSol (they are phase models). |
| Q29 | Does ThermoSol own the concentration-scale conversions (molality, mole fraction, surface density), or does GEMS3K pass concentrations? | ThermoSol owns them (per-phase evaluation, 12 §4). |
| Q30 | Sorption and surface complexation (sites, densities, electrostatic models) in the first headless builder (13 S2), or later? | In S2, with GEMSGUI's arrays. |
| Q31 | Should `system.yaml` also hold processes and sweeps (titrations, T series)? | No; a separate `run.yaml` later. |
| Q32 | Where is the `system.yaml` schema documented? | `Docs/system-yaml.md` plus a JSON Schema file that the reader uses. |
| Q33 | Should GEMSGUI write `system.yaml` too? | Yes: first through the S3 converter, then directly. |
