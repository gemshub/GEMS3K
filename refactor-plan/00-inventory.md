GEMS3K commit: 54b5f819ba83589daa9cd9aa0a9ff699adb28b1d (branch develop_optima, 2026-11-03)

# 00 — Inventory

## 1. Repositories used (read only; nothing changed, committed or pushed)

| Repository | Branch / head | Role here | Status |
|---|---|---|---|
| gemshub/GEMS3K | develop_optima / 54b5f81 | subject | read in full (see 09 for depth) |
| gemshub/optima (fork) | optima_develop / 145aeed | solver behind AOP/SOP/HOP/SHP/ROP | adapter-facing classes read |
| gemshub/xgems | develop_optima / c1ea07f | C++/Python consumer | API use listed (02) |
| gemshub/GEMSGUI | develop_gems3kv5 / 1fa83316 | GUI consumer, subclasses TMultiBase/TNode/TNodeArray | API + BASE_PARAM I/O listed (02) |
| thermohub/thermofun | v2 / 824cd14 | thermodynamic data engine | v2 rules read (legacy-and-v2.md, autodiff.md) |
| orchestra_cpp | main / f88bb1b | design study (07) | read |
| reaktoro | main | design study | layout + activity model API read |
| solvers/{cea, MAGEMin, phreeqc3, thermochimica, Equilipy, qpax} | default | design study | layout read |
| gemshub/gems-benchmark | develop / dc220de (main 5651c52 has no freeze tools) | freeze tools, corpus, CTest, recheck.py | not accessible during the planning; granted later in the session. CTest and the CI check ran (03 §7); the freeze cannot run here (08 §1, Q21) |

## 2. GEMS3K library sources (`GEMS3K/GEMS3K/`, 77 files, 74 035 lines)

Grouped by what the code does (today the file names often do not say this). "CRLF" marks the two files with
Windows line endings (trap 9). Hidden state is listed in §4.

### 2.1 Public node API and node arrays (consumers call these)

| File | Lines | Holds |
|---|---|---|
| node.h | 1636 | `class TNode`: GEM_init/GEM_run/GEM_read_dbr/GEM_write_dbr, ~150 getters/setters, name/index maps, Set_pH/Eh_target, GEM_trace_regimes, GEM_rescue_step, layer declarations l.1460–1474 |
| node.cpp | 2477 | GEM_run layer stack (l.297–1663), 25 thread-locals (l.203–228), packDataBr/unpackDataBr (l.1730–1950), GEM_from_MT ×6, GEM_to_MT, CalcTime |
| node2.cpp | 1353 | name lookup, index maps, DC_G0…DC_U0 grid accessors, Set_DC_G0, Get_Psat |
| node_copy.cpp | 762 | GEM_init (file and string forms), read/write DBR, datach/databr copy, ThermoFun engine load, `load_all_thermodynamic_from_thermo` |
| node_activities.cpp | 488 | TNode front end of `TActivity` (setTemperature, setSpeciation, updateActivityCoefficients, getIntegralPhaseProperties, …) — **seed of ThermoSol** (11) |
| node_kinetics.cpp | 893 | `TKinetics` (KM_* functions; the TMultiBase copy is in ipm_chemical4.cpp) — **seed of ThermoKin** (11) |
| node_format.cpp | 293 | VTK output of DATABR |
| node_trace.cpp | 208 | GEM_trace_regimes (snapshot + restore of MULTI and DATABR) |
| nodearray.h / .cpp / nodearray2.cpp | 794 / 810 / 849 | TNodeArray: OpenMP node loop, work-node copy in/out, transport helpers, logs, VTK; `static TNodeArray* na` |
| activities.h / kinetics.h | 400 / 212 | `ACTIVITY`, `KINETICS` structs (copies of MULTI fields) |

### 2.2 Data model and input/output

| File | Lines | Holds |
|---|---|---|
| ms_multi.h | 2105 | `MULTI` (C struct, ~400 raw pointers), `BASE_PARAM` (l.63–539, ~200 fields, 11 "reserved"), `class TMultiBase` (≈300 methods), resolvers l.587–675, enum `MULTI_DYNAMIC_FIELDS` l.2056–2098 |
| datach.h / databr.h | 142 / 266 | `DATACH` (system definition, T/P grids), `DATABR` (node data), status codes (databr.h l.139–148) |
| datach_api.h / .cpp | 188 / 848 | free functions `dbr_dch_api::*`: grid lookup, string/stream I/O, realloc/free/reset |
| datach_formats.cpp | 1036 | DATACH/DATABR text, JSON, binary formats (templated on the I/O policy) |
| ms_multi_format.cpp | 1718 | IPM file read/write (descriptor tables, 164-row count ×2), second copy of physical constants |
| ms_multi_file.cpp | 1387 | set_def (MULTI defaults), binary to_file/from_file, text dump, solmod JSON dump |
| ms_multi_copy.cpp | 1242 | copyMULTIData (copies scalars by memory layout, counts 39/53), multi_realloc (~660 lines of new[]), multi_kill |
| ms_multi_diff.cpp | 481 | `pa_p_` defaults, TMultiBase ctor, 23 alloc_* virtuals, testMulti, GUI hook stubs, thermodynamic loading from the grid |
| io_template.h/.cpp, io_keyvalue.*, io_nlohmann.*, io_simdjson.* | 399 / 491 / 416 / 616 | policy-based I/O (`TPrintArrays<TIO>`, `TReadArrays<TIO>`, `outField` descriptors) |
| gdatastream.h / .cpp | 132 / 302 | binary stream with byte-order handling |
| gems3k_impex.h / .cpp | 222 / 186 | `.lst` file parsing (`GEMS3KGenerator`), file modes |
| jsonconfig.h / .cpp | 246 / 470 | `GemsSettings` singleton (`gems3k-config.json`), logger set-up |

### 2.3 Solver kernel (GEMS3K IPM solver)

| File | Lines | Holds |
|---|---|---|
| ipm_simplex.cpp | 1069 | `CalculateEquilibriumState` (native entry), Scale/RescaleSystem, InitalizeGEM_IPM_Data, MultiConstInit, GEM_IPM_Init, simplex (Karpov–Chudnenko) |
| ipm_main.cpp | 5814 | trace writers (RUN/SET/EFF/KEY/CERT/DECIDE), GibbsEnergyMinimization, GEM_IPM, MBR, InteriorPointsMethod, linear solves (normal equations, Householder QR, Eigen QR, JAMA), Dikin criterion, step size, golden section, fill-out, certification, determinacy checks |
| ipm_chemical.cpp (CRLF) | 2166 | primal chemical potentials, GX, TotalGibbsEnergy, PSSC, StabilityIndexes, PhaseInsertionCeiling, PhaseSelect, Set_DC_limits, ConvertDCC |
| ipm_chemical2.cpp | 1443 | concentrations, sorption (Gouy–Chapman, IS_EtaCalc), phase sums |
| ipm_chemical3.cpp | 1212 | CalculateActivityCoefficients dispatch, SolModCreate (model factory), smoothing factor, PhaseSpecificGamma |
| ipm_chemical4.cpp | 859 | kinetics/metastability (TKinMet driving), RunKineticsStep |
| ipm_level.cpp | 386 | leveled start, column-generation seed, TwoPhaseSimplexMinSum (LP), NativeTpdPhase use |
| num_methods.h / .cpp | 180 / 504 | interpolation (LagranInterp), small numerics |
| tnt*.h, jama_*.h | 1029 / 583 | TNT arrays, JAMA Cholesky/LU (vendored, old) |

### 2.4 Optima adapter (built with `USE_OPTIMA_SOLVER`)

| File | Lines | Holds |
|---|---|---|
| ipm_optima.h | 115 | control-condition types, `RedoxFreeCheck`, adapter declarations |
| ipm_optima.cpp | 7195 | Optima entry `CalculateEquilibriumStateOptima` (l.2464–6870, one function of ~4 400 lines with ~40 retry/repair tiers), HOP (l.6874), NativeOptimaFinish/Rescue (l.7085/7157), pre-solve (l.1745), PotentialSpaceFinish (l.977), LP seeds, twin labelling, pH/Eh targets |

### 2.5 Activity, kinetics and sorption models

| File | Lines | Holds |
|---|---|---|
| s_solmod.h | 2456 | `SolutionData`, `TSolMod` base and ~30 subclass declarations |
| s_solmod.cpp | 498 | TSolMod base (aliases MULTI arrays), multisite helpers |
| s_solmod2.cpp | 6230 | fluid EoS: TCGFcalc, TCORKcalc, TPR78calc, TPRSVcalc, TSRKcalc, TSTPcalc |
| s_solmod3.cpp | 3679 | solid/liquid solutions: Ideal, Margules, Regular, Subregular, Guggenheim, Redlich–Kister, Van Laar, NRTL, Wilson, Berman, CEF, MBW |
| s_solmod4.cpp | 4811 | aqueous: Limiting law, Debye–Hückel, Davies, Helgeson, Karpov, Shvarov, ELVIS |
| s_solmod5.cpp | 2835 | EUNIQUAC, Pitzer, SIT |
| s_solmod6.cpp | 187 | TModOther |
| s_activity.cpp / 2 / 3 | 1003 / 486 / 1142 | `TActivity`: where GEMS links the phase models (SolModCreate…DarkenProp, concentrations, primal potentials); same code as in TMultiBase — **seed of ThermoSol** (11 §2) |
| s_kinmet.h / s_kinmet.cpp (CRLF) | 623 / 1460 | TKinMet kinetic rate models |
| s_sorpmod.h / .cpp | 416 / 142 | TSorpMod — unfinished stub (sorption lives in TMultiBase, ipm_chemical2.cpp) |

### 2.6 Support

| File | Lines | Holds |
|---|---|---|
| m_const_base.h | 379 | constants, char-coded enums (DC classes, model codes, phase codes) |
| verror.h / v_detail.h / v_detail.cpp | 88 / 232 / 72 | `TError`, `Error`/`ErrorIf`, cancel flag, `gems3k` logger |
| v_service.h / .cpp | 205 / 293 | string helpers |
| gems3k_timing.h | 92 | `GEMS3K_TIMING_FILE` per-part timer |
| gems3k_version.h | 7 | **tracked file overwritten by every CMake configure** |
| lazy_file_sink.h | 76 | spdlog sink opened on first write |

## 3. Build, tests, tools, docs

- **CMake**: root `CMakeLists.txt` (190 lines; C++20; GNU flags `-pthread -Wall -Wno-…`; default Release;
  hard-coded `/usr/local/include`, `/usr/local/lib`, `/usr/local/include/eigen3`; configure writes the tracked
  `gems3k_version.h`; `BUILD_GEM2MT` points to a folder that does not exist). `GEMS3K/CMakeLists.txt` (101):
  `file(GLOB *.cpp)` → one library, built twice (shared + static). Options `USE_OPTIMA_SOLVER` (OFF),
  `USE_THERMOFUN` (ON), `ENABLE_BENCHMARK_DIAGNOSTICS`, `BUILD_TOOLS`, `BUILD_SOLMOD(_PYTHON)`, `BUILD_TESTING` (OFF).
- **Second build system**: qmake `gems3k.pri`, `node-gem/*.pro`, `gemcalc/gemcalc.pro`, `Makefile.MinGW`,
  `Makefile.lib` — source lists kept by hand.
- **Tests**: `tests/pstall_default_test.cpp` (1 test). Its data folder `nodearray-gem-build/TestRW` is not in the
  repository, so it cannot run from a clean clone. All real regression testing is in gems-benchmark (private).
- **Tools**: `tools/` kva2json, optima_test, thermofun_test, thread_test, Cu-Pourbaix data; `node-gem/` (transport
  example), `gemcalc/` (single-node calculator), `tsolmod4rkt/` (TSolMod outside TMultiBase + Python module).
- **Docs**: `Docs/GLOSSARY.md`, `SOLVER-ADVANCED-SETTINGS.md`, `SOLVER-ADVANCED-USAGE.md`, `SOLVER-OPTIONS-GUIDE.md`,
  `SOLVER-PATH-GRAPH.md/.html`, Doxygen config.
- **Build check in this session**: conda env (gcc 16.2, Eigen 3.4.0, spdlog 1.17, fmt 12, ThermoFun 0.7.3), Optima
  fork built and installed to a private prefix, GEMS3K built with `USE_OPTIMA_SOLVER=ON` and `BUILD_SOLMOD=OFF`
  (tsolmod4rkt did not find `spdlog/spdlog.h` because of the hard-coded include paths). 36 warnings at `-Wall`.

## 4. Mechanisms that hide state (complete list found; removal plan in 05, steps A3–A9)

| # | Mechanism | Where | What it hides | Removal step |
|---|---|---|---|---|
| H1 | ~25 thread-local flags between GEM_run layers | node.cpp l.203–228 (+ externs in ipm_optima.cpp l.74–113) | re-entry guards, budgets, mode switches | A6 → `RunContext` passed down |
| H2 | thread-local outputs of solver helpers | ipm_level.cpp (`g_cgConverged`, `g_cgPseudo`, `g_cgGap`, `g_cgDual`), ipm_chemical2.cpp `g_aqueousAbsent`, ipm_main `g_ipmNoProgArmed/Fired`, `g_ipmLevelUsed/Off`, `g_tpdStrongSearch`, ipm_optima `g_ps0GaveUp`, `g_coldLvlUsed`, `g_lvlShortcut`, `g_lvlCgStart`, `g_finishFromSuccess`, `g_optimaQuickTry` | results returned through globals | A6 → return structs (F.21) |
| H3 | function-static thread_local caches | node.cpp `bigWarned` l.267, `warned` l.1179; ipm_optima `warnedPairs` l.759, `traceXPrev/tracePhPrev/traceCumIt` l.2608–2610; ipm_main `s_ipmKktFallbacks` l.94; SolveIpmAugmented scratch vectors | warn-once and scratch state | A7 (scratch → workspace), A8 (warn-once → logger state) |
| H4 | process-wide mutable statics | ipm_main `native_cert_dualfree_resolved` l.316 (data race), `warned` sets with mutex l.2672/2775; v_detail `g_calculationCancel`; nodearray `TNodeArray::na`; jsonconfig `GemsSettings data`; Optima `trace::state()` | cross-thread state | A8 |
| H5 | function-local memo statics in activity models | s_solmod2.cpp l.1422/1483/1525/1578/1758 (TCGFcalc), s_solmod4.cpp l.4678 (TELVIS::trapzd) | wrong-instance cache, data race | A9 (move into the instance) |
| H6 | environment variables read inside the solver | `GEMS3K_OPTIMA_NET_RESUME` (ms_multi.h l.634, changes results), `GEMS3K_FILLOUT_BUDGET` (ipm_main l.4579, changes results, read per call), probes (§3 of 09) | hidden inputs | A8 (results: into BASE_PARAM or an explicit `DiagnosticsConfig`; probes kept) |
| H7 | TSolMod objects alias MULTI arrays | s_solmod.cpp ctor; node_trace.cpp restores MULTI by value copy because the models point into it | ownership through raw pointers | A10, then ThermoSol `PhaseState`/`PhaseProps` (11 §4, B15) |
| H8 | non-idempotent activity call | `CalculateActivityCoefficients(LINK_UX_MODE)`: lnGmo accumulation and F0 smoothing via `FitVar[3]` (ipm_chemical l.417, ipm_chemical3 SmoothingFactor) | call history changes results | documented in A; the smoothing stays a GEMS3K solver feature outside ThermoSol (11 §2); changed only in B |
| H9 | MULTI shared across nodes of a TNodeArray worker | nodearray.cpp l.142–170; ipm_chemical4 l.80–82 comment | previous node's answer seeds the next | A11 (make the seeding explicit and test it) |
| H10 | getters that mutate | `TotalGibbsEnergy()` (X/XF/XFA), `Get_GibbsEnergy()`, `Get_PMc`, `packDataBr()` → `ReportAbsentPhaseIndexes()` | state changes on read | A5 |
| H11 | member flags used as hidden parameters | `optima_disable_smoothing`, `optima_hop_leg`, `optima_want_sensitivity`, `hop_split`, `trace_growth_stuck`, `prev_scfact`, `phase_floor_parked`, `zeroed_floor_amount`, `warm_floor_seed`, `prevBulkMol_`, `warmInputChanged_` | per-call parameters stored on the object | A6 (call parameters), A11 (warm-start state object) |
| H12 | memory-layout copies by count | ms_multi_copy.cpp `copyValues(&pm.N, …, 39)`, `(&pm.TC, …, 53)`; GEMSGUI `readArray(&PC,10)`, `(&DG,28)` | struct layout as contract | A2 (static_assert), A3 (named copy) |
| H13 | shared settings object | `pa_standalone` shared_ptr copied between nodes; `base_param()` returns non-const | silent settings sharing | A4 (const view + explicit LocalSettings) |
| H14 | stale-on-failure outputs | GEM_run_single catch without packDataBr; TNodeArray copies back | old values look new | A5 (documented, flagged) — behaviour kept in A |

## 5. Libraries around GEMS3K (target, see 10 and 11)

| Library | Role | Source today | Plan |
|---|---|---|---|
| ThermoFun (thermohub/thermofun, v2) | standard-state properties at T, P; water properties | GEMS3K grid interpolation (DATACH arrays + `LagranInterp`) or ThermoFun; GEMSGUI computes the grids (`s_supcrt*`, `m_dcomp*`, `m_reac*`, `ms_mtparm`) | the only source after the comparison campaign (10) |
| ThermoSol (new) | phase properties: activity and mixing models, excess properties, concentrations | `TActivity` + TSolMod family + TSorpMod + tsolmod4rkt | 11 §2 |
| ThermoKin (new) | kinetics: rate laws, metastability bounds, surface areas | `TKinetics` + TKinMet | 11 §3 |
