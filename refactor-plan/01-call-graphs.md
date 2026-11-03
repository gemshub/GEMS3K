GEMS3K commit: 54b5f819ba83589daa9cd9aa0a9ff699adb28b1d (branch develop_optima, 2026-11-03)

# 01 — Call graphs

Arrows mean "calls". `[TL: x]` marks a thread-local flag set or read on that edge; `[M: x]` a TMultiBase member
flag; `[E: x]` an environment variable. File and line refer to the commit above.

## 1. Initialisation

```
TNode::GEM_init(lst)                                   node_copy.cpp l.107
├─ clear_ThermoEngine
├─ GEMS3KGenerator(lst)                                 gems3k_impex.cpp  (.lst parsing, file mode)
├─ binary:  dbr_dch_api::datach_from_file  → TMultiBase::read_multi           (GemDataStream)
│  text/json: dbr_dch_api::read_dch_format_stream (datach_formats.cpp)
│             [f_thermofun] load_ThermoEngine(...)      node_copy.cpp l.590
│             TMultiBase::read_ipm_format_stream → from_text_file_gemipm<TIO>  ms_multi_format.cpp l.861
│                 ├─ set_def()  (virtual; GEMSGUI overrides)
│                 ├─ static pass: TReadArrays(8 fields) … until <END_DIM>        (trap 2)
│                 ├─ multi_realloc(PAalp, PSigm) (virtual) → 23 alloc_* (virtual)
│                 └─ dynamic pass: TReadArrays(164 fields) switch(f_…)          (trap 1)
├─ Pai/Tai from DATACH grid
├─ read_dbr_format_file(first DBR) → databr_from_text_file / databr_from_file
└─ init_into_gems3k()  (virtual; GEMSGUI: empty)
Errors: TError/std::exception → return 1 + ipmlog_error; unknown → −1.
Second form: GEM_init(dch_json, ipm_json, dbr_json, thermo_json) l.234 (strings, used by xgems/Python).
```

## 2. One equilibrium call: the layer stack (node.cpp)

```
TNode::GEM_run(uPrimalSol)                        l.297
├─ WarnIfSystemTooLarge                          [static thread_local bigWarned]
├─ ApplyDbrTargets(pHtarget/Ehtarget)            → TMultiBase::SetControlCondition_pH/_Eh
├─ restoreEnteredComposition; bBefore = bIC
├─ g_errorsHandledByCaller++                      [TL]
└─ GEM_run_rescue                                 l.464   (pa_FailRescue; six modes only)
   ├─ snapshots D0 / pool (DATABR copies)         [snapshot pattern ×12]
   ├─ GEM_run_capped  ← first attempt and each rescue step (1 warm, 2 nudges, 3 walk, 4 LocalSettings+looser,
   │                                                         5 other solver, 6 best balanced as BAD)
   └─ native_trace_quiet(true/false) around rescue steps
GEM_run_capped                                    l.790   (pa_OptimaCallCap; AOP/SOP without kinetics)
   [TL: g_optimaCapCall (re-entry guard), g_optimaIterBudget, g_optimaCapAllowance]
   └─ GEM_run_dry_redo (nudged retries, last try without cap)
GEM_run_dry_redo                                  l.908   (pa_OptimaDryRetry)
   [TL: g_optimaDryRedoCall, g_optimaDryFailed, g_optimaDryOff; adjusts g_optimaIterBudget]
   └─ GEM_run_tries
GEM_run_tries                                     l.940
   ├─ TMultiBase::BeginCalculation()  (resets g_optimaTargetWarned, g_optimaPreSolveEmpty, g_lvlFullCalls — per ATTEMPT)
   ├─ quick first try (pa_OptimaFastFirst, AOP)   LocalSettings + [TL: g_optimaFixAccept, g_optimaQuickTry]
   │    fast-cut accept → AbsentPhaseWorstDrive, PhaseHessianTwins, NativeTpdPhase [TL: g_tpdStrongSearch]
   ├─ GEM_run_single  (main attempt)
   ├─ line-search retry / stall-escape retry      [TL: g_optimaLineSearchRetry, g_optimaLSEscapeOff, g_optimaEscRetryNext]
   ├─ IpmLevelStart retry (natlvlretry)           [TL: g_ipmLevelOff, g_ipmLevelUsed]; DATABR snapshot lvlBefore
   ├─ GEM_run_cold_retry          l.1315 (only when pa_FailRescue <= 0)
   ├─ GEM_run_optima_cold_retry   l.1249 (SOP/SHP)
   ├─ GEM_run_aop_cgseed          l.1233 (AOP)    [TL: g_optimaCgSeedArmed]
   ├─ GEM_run_single_legacy_retry l.1192
   ├─ WarnCollapsedTwins; aqueous-gone note       [TL: g_aqueousAbsent]
   └─ (each side entry calls GEM_run_single again)
GEM_run_single                                    l.1410
   ├─ not NEED_* → return old status              (trap 5)
   ├─ check_TP; pNP = warm ? 1 : 0; unpackDataBr(uPrimalSol)        DATABR → MULTI
   ├─ native_trace_run_header (RUN/SET/EFF)       [E: GEMS3K_NATIVE_TRACE_FILE]
   ├─ kinetics set-up (kTau, kdT, ITau, pKMM); hop_split reset      [M: hop_split]
   ├─ dispatch by NodeStatusCH:
   │   AIA/SIA  → [TL: g_ipmNoProgArmed] CalculateEquilibriumState()           §3
   │              catch E06IPM/E07IPM → NativeOptimaRescue()                   §5
   │              ipmnoprog retry cold; NativeOptimaFinish()                   §5
   │   AOP/SOP/ROP → CalculateEquilibriumStateOptima(…, referenceMode = ROP)   §4
   │   HOP/SHP  → CalculateEquilibriumStateHOP(…, warm = SHP)                  §5
   ├─ packDataBr()  (virtual; GEMSGUI overrides)  MULTI → DATABR; ReportAbsentPhaseIndexes() (activity calls)
   ├─ IterDone; testMulti() → OK/BAD by NodeStatusCH (as rewritten by packDataBr)
   ├─ native_trace_run_result (KEY/CERT)
   └─ catch(TError) → ERR_*, NO packDataBr (trap 6); catch(std::exception) → T_ERROR_GEM
```

## 3. GEMS3K IPM solver (native kernel)

```
TMultiBase::CalculateEquilibriumState(NumIterFIA, NumIterIPM, runKinetics)   ipm_simplex.cpp l.509
├─ InitalizeGEM_IPM_Data()                       l.847
│   ├─ MultiConstInit() → multiConstInit_PN() (virtual: PZ = pa DW)
│   ├─ initalizeGEM_IPM_Data_GUI() (virtual hook)
│   ├─ pTPD < 2 → DC_LoadThermodynamicData() → load_all_thermodynamic_from_thermo | _from_grid   (trap 4)
│   ├─ Alloc_internal() → Build_compressed_xAN()  (new[] every call)
│   ├─ Alloc_uDD; RescaleToSize(true) (no-op)
│   └─ warm: loadData(false) (virtual), X = Y, TotalPhasesAmounts, CalculateConcentrations, Gamma
├─ [M: trace_growth_stuck = false; phase_floor_parked.clear(); warm_floor_seed ⇄ zeroed_floor_amount]
├─ ExcludeRedundantDCs()  + RAII RedundantRestore
├─ RunKineticsStep()  (ipm_chemical4.cpp l.57) → TKinMet models   (if runKinetics)
├─ pa_DG > 1e-5: SystemTotalMolesIC(); ScaleSystemToInternal(ScFact)                     (trap 3)
├─ [M: prev_scfact, warmInputChanged_, prevBulkMol_]
├─ GibbsEnergyMinimization()                       ipm_main.cpp
│   ├─ GEM_IPM_Init()  (ipm_simplex l.944) → GEM_IPM_Init_gui1/2 (virtual), CalculateActivityCoefficients(LINK_TP_MODE)
│   ├─ cold: GEM_IPM_InitialApproximation → AutoInitialApproximation() → SolveSimplex (SPOS/START/NEW/WORK/FIN)
│   │        or leveled start (pa_IpmLevelStart) → ColumnGenerationSeed (ipm_level.cpp) [TL: g_cg*]
│   └─ GEM_IPM(rLoop)
│       ├─ MassBalanceRefinement (MBR): MakeAndSolveSystemOfLinearEquations(N, true) (JAMA Cholesky→LU),
│       │   DikinsCriterion(true), StepSizeEstimate(true), Y += LM·MU; strict guard l.4067 (brief §8)
│       ├─ InteriorPointsMethod(status)            l.4104
│       │   loop IT < pa_IIM:
│       │     PrimalChemicalPotentials(F, Y, YF, YFA)
│       │     WeightMultipliers(false)
│       │     MakeAndSolveSystemOfLinearEquations(N, false)
│       │         └─ SolveIpmAugmented(N): hand Householder QR | Eigen ColPivHouseholderQR (La ≥ 100) | → normal eq. (JAMA)
│       │     u-divergence (PLLG); DikinsCriterion → PCI; StepSizeEstimate → LM
│       │     OptimizeStepSize(LM): golden section → GX(·) ≥ 2× ; FX = GX(LM1)          (trap 12: truncation l.765)
│       │     PD ≥ 2: CalculateActivityCoefficients(LINK_UX_MODE)  (EVERY iteration)      §6
│       │     [E: GEMS3K_IPM_PROBE] IPMIT line
│       │     PCI ≤ DXM → converged; stall window; [TL: g_ipmNoProgArmed → g_ipmNoProgFired] → E06IPM
│       │   Restore_Y_YF_Vectors()  (whole-phase zeros)
│       ├─ PSSC: PhaseSelectionSpeciationCleanup(kfr, kur, cleanup)  gate PC ≥ 2 && nCNud ≤ 0 (ipm_main l.3197)
│       │     StabilityIndexes(), PhaseInsertionCeiling(), DC_RaiseZeroedOff(); loop back to MBR/IPM (K2 < 6)
│       │   PC == 1: KarpovsPhaseStabilityCriteria / PhaseSelect / SpeciationCleanup (old path l.3283)
│       ├─ fill-out (FilloutBudgetValue [E: GEMS3K_FILLOUT_BUDGET per call]), certification, determinacy,
│       │   trace-growth check (sets [M: trace_growth_stuck] for NativeOptimaFinish), twin split
│       └─ MBR again; CheckMassBalanceResiduals
├─ RescaleSystemFromInternal(ScFact)
└─ return t_elap_sec (clock())
```

## 4. Optima path (ipm_optima.cpp)

```
CalculateEquilibriumStateOptima(NumIterFIA, NumIterIPM, referenceMode, runKinetics)     l.2464
├─ RAII guards: SmoothingGuard [M: optima_disable_smoothing = true], FloorRecordGuard (PhaseHessianEnd),
│               LuRecord/SnapRecord (Optima trace::state() counters), RedundantRestore, RedoxTargetRemove, PrevScFactGuard
├─ InitalizeGEM_IPM_Data(); ExcludeRedundantDCs(); RunKineticsStep() (not on the HOP leg)
├─ CheckRedoxFree() → SetControlCondition_Eh (pa_OptimaRedoxFree ≥ 2)
├─ ScaleSystemToInternal (pa_DG > 1e-5); SubFloorElementCheck
├─ GEM_IPM_Init_gui1 (virtual); CalculateActivityCoefficients(LINK_TP_MODE)
├─ seed:
│   warm (pNP=1 and any pm.U ≠ 0): Y as is
│   cold: ColumnGenerationSeed [TL: g_optimaCgSeedArmed] | leveled start (pa_OptimaLevelStart) [TL: g_tpdStrongSearch,
│         g_lvlShortcut, g_lvlCgStart] | LPFeasibilitySeed (TwoPhaseSimplexMinSum) | cold level (pa_OptimaColdLevel,
│         twin copies apart) [TL: g_coldLvlUsed] | DetectSolventCollapseAndReseed | DetectPhaseCollapseAndReseed | uniform
├─ dim-reduction pre-solve OptimaReducedPreSolve (cold; 2 attempts: tol 10, then rank −3) [TL: g_ps0GaveUp, g_optimaPreSolveEmpty]
│     each pass: Optima::Problem/State/Solver, own objective lambda (copy of the full one), PreStallWatch (wall clock, trap 10)
├─ [E: GEMS3K_PRESOLVE_HANDOVER_PROBE] dump
├─ Problem build: Aex = A (+ titrant columns), be = B, bounds from DLL/DUL/dcFloor, c = B and bec = I if
│     [M: optima_want_sensitivity]; objective lambda (§6 per evaluation; FD Hessian columns for opts.ibasicvars;
│     PhaseHessianFloor blocks with SymEigFloorInPlace); memoize_objective_if_linesearch
├─ State: x = max(Y, xlower); ye = −U (warm / after pre-solve) | g_cgDual | LPGibbsDual (cold)
├─ Options: maxiters max(2000, IIM), tolerance OptimaTol, line search / LU / snap (fork macros), quick-try fix-and-accept
├─ StallWatch + convergence.check hook: quick-try watchers [TL: g_optimaQuickTry], budget [TL: g_optimaIterBudget],
│     early-stability trend/cap, dry-out, frozen cut, tail stop, CalculationCancelled, wall clock (trap 10)
├─ PRIMARY SOLVE timed_solve(solver, problem, state[, sensitivity]) → dn/db read-out (optima_dndb) from THIS attempt
├─ retry/repair tiers (one function; DECIDE names): earlystop → dryout/tailpin/tailaccept/dryback/drycheck/drywet →
│   waterskip/waterkeep → dryretry → extinction/extallow (refills budget from g_optimaCapAllowance) → phasesel →
│   twinreadmit → netunpin → finishfirst (PotentialSpaceFinish) → warmfailfast/stallnet → netresume/netresolve/netfallback
│   ([E: GEMS3K_OPTIMA_NET_RESUME] via optima_net_resume_mode()) → tracerespread
├─ CalculateConcentrations; collapse check; KKT sign check (trap 11) l.5858; dual-free; titrant check
├─ OptimaZeroAbsent (optimarepair, zeroabsent); TPD accept (tpdaccept/-repair/-rejected/tpdveto/-dry); aqsanity; floorspecies
├─ CheckMassBalanceResiduals (absolute cutoff, trap 11) + relative per-IC test; verdict: MK = 2 soft fail, W21IPM
├─ StrandedElementCheck
├─ catch TError → Rescale, rethrow; catch std::exception (Optima throws runtime_error) → Rescale, E91IPM
├─ RescaleSystemFromInternal; OptimaTwinLabel
└─ FX = TotalGibbsEnergy() (save/restore X, XF, XFA); FI1/FI1s (GEMSGUI record size)
```

## 5. Hybrid and finish paths

```
CalculateEquilibriumStateHOP(…, warmNative)                                   l.6874
├─ warm requested but all pm.U == 0 → cold native (warn once per system)
├─ native leg: pNP = warm; CalculateEquilibriumState()  (§3)
│     catch E06IPM with [TL: g_ipmNoProgFired] → Optima warm from last state, then cold
│     warm failed → native cold; both failed → Optima cold
├─ HopLegGuard [M: optima_hop_leg = nativeOk]; HopLegRecorder (destructor writes [M: hop_split] + DECIDE hop-legsplit)
└─ Optima leg: CalculateEquilibriumStateOptima(…, runKinetics = false)
      on throw: restore Y, X, U, FX, B, controls; activity + concentrations; MK = 2 (BAD_GEM_HOP)

NativeOptimaFinish (pa_NativeOptimaFinish > 0 && [M: trace_growth_stuck])     l.7085
└─ snapshot 13 arrays; Budget RAII [TL: g_optimaIterBudget = cap; ~ → −1]; [M: optima_hop_leg]; Optima warm;
   accept if clean and G lower by > 1e-8, else restore everything
NativeOptimaRescue (after E06IPM/E07IPM)                                       l.7157
└─ two attempts (warm with hop leg, then cold), same Budget pattern
```

## 6. One activity call

```
TMultiBase::CalculateActivityCoefficients(LinkMode)           ipm_chemical3.cpp
├─ LINK_TP_MODE: SolModCreate per phase (factory switch on sMod codes, ~30 classes) → TSolMod(SolutionData*)
│     (SolutionData points INTO MULTI arrays: x, lnGamma, lnGamConf…, aIPc, aGEX …) ; PTparam()
├─ LINK_UX_MODE (per evaluation):
│     CalculateConcentrations(X, XF, XFA) → CalculateConcentrationsInPhase (all phases), IS_EtaCalc, phase_bcs
│     per phase k: TSolMod::MixMod() (writes lnGamma via alias) ; PhaseSpecificGamma ; DC_PrimalChemicalPotentialUpdate
│        → F0 = Fold + dF0·SmoothingFactor()  (unless [M: optima_disable_smoothing])   (history: FitVar[3], lnGmo)
│     sorption: Gouy–Chapman / electrostatic terms; calculateActivityCoefficients_scripts (virtual, GEMSGUI)
└─ TKinMet / TSorpMod not involved (kinetics in RunKineticsStep)
Callers: InteriorPointsMethod (every iteration, PD ≥ 2), GEM_IPM_Init, Optima objective (every evaluation and every
FD column), PotentialSpaceFinish (every FD column), pre-solve, TPD searches, ReportAbsentPhaseIndexes (from packDataBr),
HOP restore, NativeOptimaFinish restore.
```

## 7. Node arrays (TNodeArray)

```
TNodeArray::CalcIPM_List(mode, start, end)            nodearray.cpp l.142 (OpenMP) / l.234 (serial)
└─ per thread: TNode wrkNode(*calcNode)  (copy once; reused for every node of the chunk)
   └─ CalcIPM_Node(ii) → SmartMode, NeedGEMS → RunGEM(wrkNode, ii, Mode, C1)       l.529
        ├─ CopyWorkNodeFromArray(ii)  DATABR array → work DATABR (field by field)
        ├─ NodeStatusCH = |Mode|; CalcNodeServer → wrkNode->GEM_run(uPrimalSol = Mode<0 || SIA)    §2
        └─ MoveWorkNodeToArray(ii)    work DATABR → array (also after a failure)
   success test: OK_GEM_AIA || OK_GEM_SIA only
TNodeArray::RunGEM(Mode, nNodes, …)  l.291 uses static TNodeArray::na
GEMSGUI: TNodeArrayGUI overrides CalcIPM_List/CalcIPM_Node.
```

## 8. Data ownership today

```
TNode ── owns ─▶ DATACH* CSD, DATABR* CNode, TMultiBase* (multi), TActivity, TKinetics, ThermoEngine (shared_ptr)
TMultiBase ── owns ─▶ MULTI pm (~400 raw new[] arrays, multi_realloc/multi_kill), TSolMod* phSolMod[], TKinMet*
               phKinMet[], TSorpMod* phSorpMod[], AA/BB/arrL/arrAN, shared_ptr<BASE_PARAM> pa_standalone (shared
               between copies), std::vector caches (prev_*, *_floor_*, optima_*)
TSolMod ── aliases ─▶ MULTI arrays (x, lnGamma, …) ; owns aIP, multisite y/mn/mns
TMultiBase::node1 ─▶ back-pointer to TNode (used for names, ThermoFun, T/P)
TNodeArray ── owns ─▶ calcNode (TNode), arrays of DATABR (NodT0, NodT1), grid; static na
Optima::Problem/State ── own ─▶ Eigen copies of A, b, bounds (rebuilt every solve); objective lambda captures `this`
```
