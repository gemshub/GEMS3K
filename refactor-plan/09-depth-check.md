GEMS3K commit: 54b5f819ba83589daa9cd9aa0a9ff699adb28b1d (branch develop_optima, 2026-11-03)

# 09 — Depth self-check

This file checks that the reading is deep enough to plan from. It was written before `04-target-architecture.md`,
as the task asks. Line numbers refer to the commit above. "Verified" means the line was read in this session.

Result: **all 14 traps found and verified, 0 items "not found", 0 guessed.** The pass rule (at most two "not found"
or guessed items) is met on the first pass.

What was read: every header in `GEMS3K/GEMS3K/` in full; every solver source in full (`node*.cpp`, `ipm_*.cpp` incl.
the 7 195-line `ipm_optima.cpp`, `ms_multi_*.cpp`, `nodearray*.cpp`, `s_solmod.cpp`); the I/O, data and support
files by function index plus the parts that touch traps, hidden state or the API; the activity models
(`s_solmod2-6.cpp`, 17 742 lines of formulas) and `s_kinmet.cpp`/`s_activity*.cpp` by structure (class list,
constructors, aliasing, statics, fixed arrays, allocation) and by targeted reads. This split is recorded in
`NOTES.md` per file. The model formulas themselves were not checked line by line; the plan does not change them.

## 1. Traps of brief §3

| # | Trap | Where (file : function : line) | How the code makes it happen |
|---|------|------------------------------|------------------------------|
| 1 | `BASE_PARAM` is a cross-repo contract; 7 edit sites, 2 hard-coded counts | `ms_multi.h` struct `BASE_PARAM` l.63–539, enum `MULTI_DYNAMIC_FIELDS` l.2056–2098; `ms_multi_format.cpp` descriptor table `MULTI_dynamic_fields` l.72–254, writer `to_text_file_gemipm` l.261 (`TPrintArrays<TIO> prar( 164, ...)` **l.272**), reader `from_text_file_gemipm` l.861 (`TReadArrays<TIO> rddar( 164, ...)` **l.1014**, `case f_pa_...` labels); `ipm_main.cpp` `native_trace_run_header` SET line l.166–222; `ms_multi_diff.cpp` `pa_p_` l.6–23 (positional aggregate of the first 37 fields); GEMSGUI `Modules/m_param.cpp` `BASE_PARAM::read/write` l.236–247 (`readArray(&PC,10)`, `readArray(&DG,28)`) | The count 164 is typed twice by hand. A script check in this session found 164 descriptor rows = 164 enum entries, same order, every enum read; so the counts are right today. The table, the enum, the writer and the reader are four separate lists that only agree by discipline. GEMSGUI reads the first 38 members by memory position and defines `BASE_PARAM::read/write` itself (declared in GEMS3K `ms_multi.h` l.537–538). GEMSGUI also binds 42 newer fields by name (`o_pax*` data objects, `Modules/v_mod.h`, `Resources/data/vis_od.ini.json`). A new field (after `tprn`) needs 8 edits in GEMS3K (the brief's 7 sites, with the two counts counted separately); a field shown in GEMSGUI needs up to 5 more there (enum entry, `SetPtr`, `vis_od.ini.json`, `ResetTail`, defaults). |
| 2 | No runtime setter for `pa_`; a field before `<END_DIM>` is never read | `ms_multi_format.cpp` writer emits `<END_DIM>` after the 8 static fields l.341; reader static pass `rdar( 8, MULTI_static_fields, ...)` l.905; `io_template.h` `TReadArrays::findNext` l.273–287; `io_keyvalue.cpp` `KeyValueRead::has_next` l.123–148; `ms_multi.h` `pa_standalone` l.1138, `base_param()` l.1180, `LocalSettings` l.1188–1202 | The static pass knows only 8 names. `findNext` skips every unknown name silently, and `has_next` reads the next token as a name (a value token is also eaten as a "name") until it meets `END_DIM`. So a `pa_` field placed there is dropped with no message. Settings live in `TMultiBase::pa_standalone` (a `shared_ptr` shared by node copies). `TNode` has no setter; `base_param()` is `const` but returns a non-const pointer, so code with `multi_ptr()` can still write settings. |
| 3 | Rescaling by `ScFact = pa_DG/ΣB`; values stored mid-solve are in internal units | `ipm_simplex.cpp` `SystemTotalMolesIC` l.636, `ScaleSystemToInternal` l.652–743, `RescaleSystemFromInternal` l.747–839, call in `CalculateEquilibriumState` l.538–542; Optima entry `ipm_optima.cpp` l.2594–2597 and exits l.6808, l.6821, l.6830 | Scale-in multiplies MBX, FitVar[0], VXc, HXc, HX_, FX, Yw, DUL (<1e6, recorded), DLL, Y, X, XY, XU, B, BFC, XF, YF, XFs, FVOL, FWGT, XFA, YFA, PUL, PLL, BF and the sorption arrays. Not scaled: U, MU, W, Falp, lnGam, kinetic outputs, titrant amounts, and every `std::vector` cache on `TMultiBase`; `zeroed_floor_amount`/`warm_floor_seed` are scaled by hand by the caller; C is not scaled on entry but is divided on exit (correct only because C is recomputed during the solve). Anything a mechanism stores during the solve is in internal units unless it converts by hand. |
| 4 | Two copies of G0 | `node2.cpp` `DC_G0` l.358 (reads `CSD->G0` grid, LagranInterp order 6); `node.h` `Get_solver_G0` l.954 (returns `pm.G0`), comment l.950; `ms_multi_diff.cpp` `load_all_thermodynamic_from_grid` (`pm.G0[j] = ConvertGj_toUniformStandardState( Go+Gg+Ge, j, k )`); `node_copy.cpp` `load_all_thermodynamic_from_thermo` l.617 | The grid path writes `pm.G0` = converted value (with `Guns` and the standard-state term), so even without ThermoFun `DC_G0()` ≠ `pm.G0`. The ThermoFun path writes `pm.G0` and `tpp_G` only and never `CSD->G0`, so `DC_G0()` returns the stale exported grid. xgems calls `DC_G0` (2 sites). |
| 5 | `GEM_run()` only recalculates for a `NEED_*` code | `node.cpp` `GEM_run_single` l.1438–1456 | A non-NEED `NodeStatusCH` returns the old status at once. The check moved from `GEM_run` into the innermost layer; the outer layers (rescue, cap, dry redo, tries) still run their set-up before it. |
| 6 | On a FAIL row every reported value is stale | `node.cpp` `GEM_run_single` catch block l.1646–1663; `GEM_run` l.314–319 rethrows after `LogError`; `nodearray.cpp` `RunGEM` l.529–546 | `catch(TError&)` sets `ERR_*` and never calls `packDataBr()`, so DATABR keeps the previous answer. `TNodeArray::RunGEM` then copies the stale work-node DATABR back into the node array (`MoveWorkNodeToArray`) on failure too. |
| 7 | `GEM_run()` is a stack of wrappers with thread-local flags | `node.cpp` `GEM_run` l.297, `GEM_run_rescue` l.464–783, `GEM_run_capped` l.790–904, `GEM_run_dry_redo` l.908–937, `GEM_run_tries` l.940–1185, `GEM_run_single_legacy_retry` l.1192, `GEM_run_aop_cgseed` l.1233, `GEM_run_optima_cold_retry` l.1249, `GEM_run_cold_retry` l.1315, `GEM_run_single` l.1410; thread-locals l.203–228; `BeginCalculation` l.245 called from `GEM_run_tries` l.944 | The named layers exist already (brief §10.1 is done), and each calls the next directly. But they still talk through about 25 thread-local flags (re-entry guards `g_optimaCapCall`, `g_optimaDryRedoCall`; budgets `g_optimaIterBudget`, `g_optimaCapAllowance`; switches `g_optimaFixAccept`, `g_optimaLSEscapeOff`, `g_optimaLineSearchRetry`, `g_optimaCgSeedArmed`, `g_ipmLevelOff`, `g_optimaDryOff`). The Optima objective's convergence hook writes `g_optimaIterBudget` (ipm_optima.cpp l.3881–3950). `BeginCalculation` resets per attempt, not per `GEM_run`. |
| 8 | `DC_name_to_xCH()` returns the first match | `node2.cpp` `DC_name_to_xCH` l.123–133; `DC_name_to_xCH_map` l.135 | Linear scan, returns the first index, logs a warning on miss. The `_map` version loops `k < CSD->nPHb` but indexes DCH arrays `nDCinPH[k]`/`PHNL[k]` (mixes DBR count with DCH index space). |
| 9 | Line endings vary per file | `file` on `GEMS3K/GEMS3K/*`: **`ipm_chemical.cpp` and `s_kinmet.cpp` are CRLF**; all others LF | Scripted edits and patch tools must keep each file's line ending; a mixed-ending commit makes a whole-file diff. |
| 10 | Only one wall-clock decision (`pa_OptimaMaxSeconds`) | `ipm_optima.cpp` `OptimaReducedPreSolve` `PreStallWatch` l.1868–1885; full solve `StallWatch` l.3695–3736, hook test l.4061 (`CalculationCancelled() || stallWatch->overBudget()`) | Both watches use `std::chrono::steady_clock` against `pa_OptimaMaxSeconds` (default 0 = off). Every other clock use is reporting only (`pm.t_start/t_end` with `clock()` in ipm_main/ipm_simplex/ipm_optima; `GEMS3K_BENCHMARK_DIAGNOSTICS` timers ipm_main l.5100–5306). Second non-deterministic stop: `CalculationCancelled()` (`v_detail.cpp` l.23–25, process-wide atomic), driven by the caller. |
| 11 | Optima zeroes the bound error whatever its sign; absolute mass-balance cutoff | Optima `Optima/ResidualErrors.cpp` l.97–109 (`ex(ju).fill(0.0)`; basic variable `== xlower || == xupper` → `ex = 0`); GEMS3K re-check `ipm_optima.cpp` l.5858–5960 (KKT sign test per species, tol `pa_GAS`); `CheckMassBalanceResiduals` l.6170 with absolute cutoff min(DHBM·1e10, 1e-2) mol = 1e-3 mol at DHBM = 1e-13 (comments l.6381, l.6462), relative per-IC test added after | Optima reports "converged" for a species held at a bound with the wrong gradient sign; GEMS3K catches it after the solve. The absolute cutoff cannot see a trace-element loss of ~1e-13 mol; the CERT record and the added relative test do. |
| 12 | PSSC + `GX()` truncation make the exact zeros; `PhaseInsertionCeiling`; `PC >= 2` gate | `ipm_chemical.cpp` (CRLF) `GX` l.747, truncation `if( pm.X[i] < pm.DcMinM ) pm.X[i]=0.;` **l.765–766**; `PhaseSelectionSpeciationCleanup` l.1390–1849 (insertion l.1523–1539, `Y<DcMinM → 0`); `PhaseInsertionCeiling` l.1363; gate `ipm_main.cpp` l.3197 `pa_p->PC >= 2 && nCNud <= 0`; `Restore_Y_YF_Vectors` ipm_main l.5440 | Exact zeros come from three places: the GX line-search truncation, PSSC elimination/cleanup, and `Restore_Y_YF_Vectors` for whole phases. The insertion amount is `min(pm.DFYsM, PhaseInsertionCeiling)` once per solve. `PC = 1` takes the old `PhaseSelect` path (l.3283). |
| 13 | Every mechanism prints SET and EFF; resolvers are single-source | `ipm_main.cpp` `native_trace_run_header` l.124 (RUN l.147, SET l.166, EFF l.260–281); resolvers `ms_multi.h` `optima_dimreduce_passes` l.594, `optima_presolve_pass_budget` l.604, `optima_net_resume_mode` l.631, `optima_multisite_phase_count` l.649, `optima_earlystability_at` l.669; DECIDE writer `native_trace_decide` (node_trace.cpp) | Solver and trace call the same inline resolver, so EFF shows what ran. Exception found: `optima_net_resume_mode()` reads env `GEMS3K_OPTIMA_NET_RESUME` once per process, so an environment variable changes solver behaviour (EFF prints "env"). |
| 14 | Diagnostics env-gated and zero-cost when unset; spdlog channels | `ipm_main.cpp` `GEMS3K_NATIVE_TRACE_FILE` l.61, `GEMS3K_IPM_PROBE` l.84 (call site l.4252); `ipm_optima.cpp` `GEMS3K_PHSTAB_PROBE` l.1663 and l.1700, `GEMS3K_EARLYTREND_PROBE`(`_N`) l.3814–3826, `GEMS3K_OPTIMA_TRACE_FILE`, `GEMS3K_PRESOLVE_RESID_PROBE`, `GEMS3K_PRESOLVE_HANDOVER_PROBE`; `gems3k_timing.h` `GEMS3K_TIMING_FILE` l.45; Optima `Trace.cpp` `GEMS3K_OPTIMA_STEPTRACE_FILE` (fallback `REAKTORO_TRACE_FILE`); loggers `gems3k` v_detail.cpp l.7, `ipm` ipm_main l.49, `tnode` node.cpp l.47, `solmod` s_solmod.cpp l.36, `kinmet` s_kinmet.cpp l.32 | Each probe is a function-static `FILE*` opened once by a `getenv` lambda (`static FILE* dbg = [](){...}()`), so the cost when unset is one branch. Two separate statics open the same `GEMS3K_PHSTAB_PROBE` path twice. Some `getenv` calls run per call, not once: `GEMS3K_OPTIMA_TRACE_FILE` (per pre-solve pass and per solve), `GEMS3K_FILLOUT_BUDGET` (ipm_main l.4579, **changes results**, every `FilloutBudgetValue` call). |

## 2. Floating-point order points (brief §4)

All verified in code; Stage A must not change any of them:

1. Sparse `MassBalanceResiduals` (ipm_main l.4388): `C[ii] -= A·Y` in species order over the compressed index
   (`Build_compressed_xAN`, drops |a| ≤ 1e-12).
2. Normal-equation assembly `AA += a·a·W` in `MakeAndSolveSystemOfLinearEquations` (l.5097), species-major.
3. Hand-written Householder QR in `SolveIpmAugmented` (l.4855, `pa_IpmAugmentedKKT = 2`), row loops r = k..M.
4. Eigen `ColPivHouseholderQR` for La ≥ `pa_IpmQrSortPivot` (100) with stable row sort — Eigen-version sensitive.
5. JAMA Cholesky → LU fallback (MBR and IPM fallback).
6. `GX()` sum per phase/species (ipm_chemical l.747) and the golden-section `OptimizeStepSize` (ipm_main l.4448).
7. `DikinsCriterion` (l.5318) and `DC_DualChemicalPotential` (ascending IC index over the compressed index).
8. `StabilityIndexes` sums (ipm_chemical l.1854).
9. Smoothing blend `F0 = Fold + dF0·SmoothingFactor()` (ipm_chemical l.417) — history-dependent.
10. Optima FD Hessian columns (ipm_optima l.3367–3410), `SymEigFloorInPlace` (Jacobi), pre-solve copy of both.
11. Activity models write into MULTI arrays through aliases (TSolMod ctor, s_solmod.cpp) — any change of call
    count changes the smoothing history (comment in `memoize_objective_if_linesearch`, ipm_optima l.120).

## 3. Five non-obvious facts (with evidence)

1. **`copyMULTIData` copies MULTI scalars by memory layout, and the counts are already off by one.**
   `ms_multi_copy.cpp` l.37–50: `copyValues( &pm.N, &otherPM->N, 39 )` and `copyValues( &pm.TC, &otherPM->TC, 53 )`.
   Checked against `ms_multi.h` l.683–760: the long-int block N…SolveCallCount has **40** fields (the last is not
   copied); the double block TC…YFk has **52** scalars and is followed by `FitVar[5]`, so the copy takes FitVar[0]
   but not FitVar[1..4] (FitVar[3] is the F0 smoothing state). Same failure class as trap 1, in a place nobody lists.
2. **Warm modes in `TNodeArray` depend on the OpenMP thread count.** `nodearray.cpp` `CalcIPM_List` l.142–170: each
   thread copies `calcNode` once and reuses the copy for all nodes of its chunk, so node *i−1*'s MULTI state
   (pm.U, pm.Y, FitVar, member caches) seeds node *i*. SOP/SHP decide "warm usable" from `pm.U` (ipm_optima
   l.2642–2650; HOP l.6886–6897), not from the DBR. With the default static schedule the chunk borders move with the
   thread count. Also `CalcIPM_Node` (l.194, l.259) counts only `OK_GEM_AIA || OK_GEM_SIA` as success, so every
   Optima-mode OK is reported as an error; and `static TNodeArray* na` (l.37) points at the last array built.
3. **Activity models hold process-wide memo caches keyed only on (T, ρ).** `s_solmod2.cpp` `TCGFcalc::DIntegral`
   l.1422 (`static double TOld, roOld; static double a,b,c,d,e;`), also l.1483, l.1525, l.1578, and `FWCA` l.1758
   (`static double TOld, roOld, F;` returned when `essentiallyEqual(T,TOld) && essentiallyEqual(ro,roOld)`). Two
   CG-fluid phases, or two nodes in two threads, at the same T and ρ can get each other's value; under OpenMP it is
   also a data race. `TELVIS::trapzd` (s_solmod4 l.4678) keeps `static double s` across calls by design.
4. **Optima's trace state is process-wide and always counted, and GEMS3K DECIDE records read it.** Optima
   `Trace.cpp` `static State s;` (not thread_local); `LU.cpp` l.142/147 and `BacktrackSearch.cpp` l.138 increment
   it on every solve, trace or not. GEMS3K `CalculateEquilibriumStateOptima` (`LuRecord`/`SnapRecord`, l.2489–2505)
   reports `lupartial`/`snapleaving` as differences of these counters, so in a parallel run a `# dec` line can
   count another thread's events.
5. **dn/db is taken from the primary Optima attempt, not from the accepted answer.** `ipm_optima.cpp` l.4515–4527:
   `optima_dndb` is filled from `sensitivity.xc` right after the first `timed_solve`; the comment says a retry
   "solves a different problem". When a retry tier produces the answer, the stored dn/db belongs to another point.
   Any extension that uses dn/db (07) must recompute it at the final state.

Further facts recorded in `NOTES.md` (each with evidence): ThermoFun's cache key holds a reference to the P
lvalue, hence `static thread_local double funP` in `load_all_thermodynamic_from_thermo` (node_copy.cpp);
`from_text_file_gemipm` loops `ii < dCH->nDC` over the nIC-sized `ccIC` (out-of-bounds read, ms_multi_format
~l.946); `NativeOptimaFinish`/`NativeOptimaRescue` reset `g_optimaIterBudget` to −1 in a destructor instead of
restoring it (l.7109–7110, 7169–7170); `TotalGibbsEnergy()` writes X/XF/XFA, so the getter `Get_GibbsEnergy()`
mutates state; `packDataBr()` runs `ReportAbsentPhaseIndexes()` (up to ~600 activity evaluations) on every call;
`RescaleToSize` is a no-op; the only test needs a data folder that is not in the repository; two environment
variables change results (`GEMS3K_OPTIMA_NET_RESUME`, `GEMS3K_FILLOUT_BUDGET`).

## 4. MULTI fields and layer calls of the seven key functions

Layer order (outermost first; each calls the next directly):
`TNode::GEM_run` → `GEM_run_rescue` → `GEM_run_capped` → `GEM_run_dry_redo` → `GEM_run_tries` →
`GEM_run_single` → {`TMultiBase::CalculateEquilibriumState` (AIA/SIA) | `CalculateEquilibriumStateOptima`
(AOP/SOP/ROP) | `CalculateEquilibriumStateHOP` (HOP/SHP)} → `packDataBr`.
Side entries from `GEM_run_tries`: `GEM_run_single_legacy_retry`, `GEM_run_aop_cgseed`, `GEM_run_optima_cold_retry`,
`GEM_run_cold_retry` (each re-enters `GEM_run_single`). Native path inside the kernel:
`CalculateEquilibriumState` → `InitalizeGEM_IPM_Data` → `ExcludeRedundantDCs` → `RunKineticsStep` →
`ScaleSystemToInternal` → `GibbsEnergyMinimization` (`GEM_IPM_Init` → `GEM_IPM_InitialApproximation`/simplex →
`GEM_IPM`: MBR → `InteriorPointsMethod` → PSSC → MBR) → `RescaleSystemFromInternal`.

| Function | Reads (main) | Writes (main) | Called by |
|---|---|---|---|
| `TNode::GEM_run(bool)` node.cpp l.297 | DATABR `NodeStatusCH`, `bIC` (copy `bBefore`), targets (`pHtarget`, `Ehtarget`) | DATABR `bIC` (restoreEnteredComposition), control conditions (ApplyDbrTargets), `g_errorsHandledByCaller`, timing flush | caller (xgems, GEMSGUI, TNodeArray::RunGEM) |
| `TNode::GEM_run_single` node.cpp l.1410 | `NodeStatusCH`, DATABR via `unpackDataBr` → MULTI `kTau, kdT, TCc, Tc, Pc, VXc, DUL, DLL, B, Aalp, Falp, XF` (+ warm: `IT, MBX, IC, Eh, Y, PUL, PLL, YF, FVOL, FWGT, BF, XFA, YFA, C, U, Gamma, lnGam`) | MULTI `pNP` (0/1), `kTau/kdT/ITau/pKMM`, `hop_split`; after solve DATABR via `packDataBr`; `NodeStatusCH` → OK/BAD/ERR; `IterDone`; `CalcTime` | `GEM_run_tries` and the four side entries |
| PSSC `PhaseSelectionSpeciationCleanup` ipm_chemical l.1390 | `Y, Gamma, GamFs (K2==0), DUL, DLL, L1, PHC, DCC, DFYsM, DSM, DcMinM, pa DF/DFM/GAS/PRD`, `B`, `A` (ceiling) | `Y` (insert / eliminate / cleanup / exact zeros), `XY` (restore copy), `GamFs`, `K2++`, `Falp`, `NMU`, `EMU` (via StabilityIndexes), `F` (PrimalChemicalPotentials), scratch `YFk/logXw/logYFk/Yw/aqsTail` | `GEM_IPM` (native, gate `PC >= 2 && nCNud <= 0`, ipm_main l.3197) |
| `GX(LM)` ipm_chemical l.747 | `Y, MU, G, L1, DcMinM, lowPosNum, LO, PHC` | `X` (= Y + LM·MU, truncated below DcMinM), `XF, XFA` (TotalPhasesAmounts), returns FX; scratch `logXw, logYFk` | `InteriorPointsMethod` (entry, golden-section probes, final step), `AutoInitialApproximation` |
| `ScaleSystemToInternal(ScFact)` ipm_simplex l.652 | `ScFact` (from `SystemTotalMolesIC`: `B`, `N`, `E`, pa `DG`) | `SizeFactor, MBX, FitVar[0], VXc, HXc, HX_, FX, Yw, DUL, DLL, Y, X, XY, XU, B, BFC, XFs, XF, YF, FVOL, FWGT, XFA, YFA, PUL, PLL, BF`, sorption `XetaA/B/D, XFTS`; records `DUL/DLL/PUL_pre/postScale_` | `CalculateEquilibriumState` l.538–542, `CalculateEquilibriumStateOptima` l.2594–2597 (both when `pa_DG > 1e-5`) |
| `TNode::packDataBr()` node.cpp l.1730 (virtual; GEMSGUI overrides) | MULTI `pNP, TCc, Pc, IT, VXc, FX (sentinel check), HXc, IC, pH, pe, Eh, MBX, XF, Aalp, Falp, FVOL, FWGT, XFA, PUL, PLL, BF, X, Gamma, DUL, DLL, B, C, U, BFC`, `phase_floor_parked` | DATABR `NodeStatusCH` (mode reset), `TK, P, IterDone, Vs, Gs, Hs, IC, pH, pe, Eh, Ms, xPH, aPH, omPH, vPS, mPS, xPA, amru, amrl, bPS, xDC, gam, dul, dll, bIC (from pm.B), rMB, uIC, bSP`; also MULTI `Falp`/`omPH` via `ReportAbsentPhaseIndexes()` (activity evaluations) | `GEM_run_single` (success paths), rescue pool snapshots, `TNodeArray::SaveToNode` |
| Optima entry `CalculateEquilibriumStateOptima` ipm_optima l.2464 | `N, L, Ls, FI, FIs, L1, PHC, A, B, G0, fDQF, DLL, DUL, LO, SB, SM, SF, ICC, pNP, stkey, sMod, T, XwMinM, DSM, PhMinM, DCCW, U` (warm test), ~60 pa_ fields | `Y, X, XF, XFA, YF, YFA, U (= −ye), F, lnGam, Gamma (activity), F0, B (titrants), C, pH/pe/Eh (CalculateConcentrations), FX (TotalGibbsEnergy), FI1, FI1s, ITF=0, ITG, MK, PZ, Ec, t_start/t_end/t_elap_sec, errorCode/errorBuf`; members `prev_scfact, phase_floor_parked, zeroed_floor_amount, optima_dndb(_rows/_cols), optima_control_conditions, phfTwin/phfFired`; thread-locals `g_optimaIterBudget, g_optimaPreSolveEmpty, g_ps0GaveUp, g_coldLvlUsed, g_lvlShortcut, g_lvlCgStart, g_tpdStrongSearch, g_optimaQuickTry`; Optima `trace::state()` | `GEM_run_single` (AOP/SOP/ROP), `CalculateEquilibriumStateHOP` (HOP/SHP, `runKinetics=false`), `NativeOptimaFinish`, `NativeOptimaRescue` |

## 5. Verdict

0 of 14 traps not found; 0 guessed; field lists verified against code. The reading is deep enough to write
`04-target-architecture.md`. Known limit: the activity-model formulas were read by structure, not line by line;
the plan treats them as a black box with a fixed interface and fixed floating-point order.
