GEMS3K commit: 54b5f819ba83589daa9cd9aa0a9ff699adb28b1d (branch develop_optima, 2026-11-03)

# Running notes — GEMS3K refactor planning

GEMS3K commit: `54b5f819ba83589daa9cd9aa0a9ff699adb28b1d` (develop_optima, 2026-11-03)

Purpose: survive a context reset. Raw facts with file:line. Polished into 00..09 later.

## Session facts
- Clones under /home/user: GEMS3K (full history, 1580 commits), optima (optima_develop 145aeed, full, 1046 commits),
  GEMSGUI (develop_gems3kv5 1fa8331), xgems (develop_optima c1ea07f), thermofun (v2 824cd14),
  orchestra_cpp (main f88bb1b), reaktoro (main f587235, depth 50), solvers/{cea 2f79a64, MAGEMin 1298456,
  phreeqc3 992c7ec, thermochimica 1c076b2, Equilipy 7a29c5e, qpax a5721f8} (depth 50).
- **gems-benchmark: NOT reachable** (private; add_repo "no access"; anonymous clone asks for credentials).
  Consequence: no freeze tools, no Resources corpus, no CTest suite, no recheck.py in this session.
  Substitute corpus for profiling: 34 `*-dat.lst` projects in GEMS3K/{node-gem-build,gemcalc-build,tools-build,
  tsolmod4rkt-build,tools/Cu-Pourbaix} and xgems/{demos/resources,tests/gems3k}.
- Machine: 4 cores, 15 GB, g++ 13.3 (system), clang++, cmake 3.28, valgrind; no perf; no conda (micromamba 2.10
  downloaded to /home/user/hub-tools, env at /home/user/hub-tools/env from conda-forge).
- micro.mamba.pm is blocked by the proxy; github releases and conda.anaconda.org work.

## Per-file notes
(appended while reading)

### ms_multi.h (2105 lines, LF) — read in full
- `BASE_PARAM` (l.63-539): 10 shorts PC,PD,PRD,PSM,DP,DW,DT,PLLG,PE,IIM; 28 doubles DG..DKIN; `char* tprn` (l.110);
  then ~70 trailing fields with default initializers (PSTALL l.118 ... OptimaTwinLabel l.535). Several "reserved, no
  effect" kept for order: OptimaMaxStepRatio, OptimaPhaseCompaction, MbTrendPhaseDecay, MbPivotSplit,
  OptimaReadmitSeed, LpDualFillout, IpmLoopTweaks, OptimaFDDiagFloor, OptimaLSWindow, GAR, GAH. Macros
  GEMS3K_DEFAULT_OPTIMA_{TPDACCEPT,CGSEED,COLDRETRY,FINISH} override compiled defaults (build-time defaults!).
  `write(GemDataStream&)/read()` binary I/O (positional).
- `MULTI_DYNAMIC_FIELDS` enum l.2056-2098: the per-field I/O index, BASE_PARAM fields from f_pa_DB.. then new
  f_pa_PSTALL..f_pa_OptimaTwinLabel (order matches struct trailing order). `MULTI_STATIC_FIELDS` l.2051.
- `struct OptimaQuickTryWatch` + `extern thread_local OptimaQuickTryWatch g_optimaQuickTry` (defined node.cpp) l.544-559
  = HIDDEN STATE #1 (thread-local quick-try watch).
- Auto-resolved settings (single source, trap 13): `optima_call_cap()` l.587, `optima_dimreduce_passes()` l.594
  (AUTO gate nDC>=200 -> 8 passes), `optima_presolve_pass_budget()` l.604, `optima_earlystability_at()` l.669
  (cap 200 cold / 25 warm if multisite), `optima_net_resume_mode()` l.631 = function-local static reading env
  GEMS3K_OPTIMA_NET_RESUME once (HIDDEN STATE #2: process-wide static from env), `optima_multisite_phase_count()`.
  `kTotalGibbsEnergyUnset = 7777777.` l.573 (FX seed sentinel; packDataBr checks).
- `MULTI` typedef struct l.678-1093: plain C struct of raw pointers; counters N,L,FI,FIs,Ls,...; K2,PZ,pNP,pESU,
  pIPN,pBAL,pTPD,ITaia,MK,W1, SolveCallCount; doubles T,P,..,FX,PCI,CondNum,SolveTimeMs,DXM,RT,Yw,...; ~150 arrays
  (A, B, U, X, Y, XY, XF, YF, XFA, YFA, G, G0, F, F0, lnGam, Gamma, Wx, W, MU, EMU, NMU, Fx, DUL, DLL, fDQF, ...);
  name arrays SB, SM, SF ...; class codes RLC, RSC, ICC, DCC, PHC, DCCW; cutoffs scaled XwMinM..DFYsM,
  SizeFactor. `t_start,t_end clock_t` and `t_elap_sec` (wall-clock read).
- `PhStabViolation`, `PhStabCensus` structs.
- `extern const BASE_PARAM pa_p_` (defaults; defined ms_multi_diff.cpp).
- `class TMultiBase` l.1134-1981: members PAalp_, PSigm_, `shared_ptr<BASE_PARAM> pa_standalone` (settings, shared
  between nodes!), DUL/DLL/PUL pre/post-scale vectors (scaling sentinel record), kinDegenWarned_, prevBulkMol_,
  prevTC_, prevP_, prevDulPositive_, warmInputChanged_; `friend class TNode`; `static shared_ptr<spdlog::logger>
  ipm_logger`; `TNode* node1` public back-pointer; `LocalSettings` RAII (copy-on-write of pa_standalone for retries);
  `GetPM()`; many public methods: TotalGibbsEnergy, Cert* report functions, ReportAbsentPhaseIndexes,
  NativeTpdPhase, PhaseHessianTwins, WarmTwinSplitCheck, ColumnGenerationSeed, AbsentPhaseWorstDrive,
  CalculateEquilibriumState(NumIterFIA, NumIterIPM, runKinetics), InitalizeGEM_IPM_Data, DC_LoadThermodynamicData;
  #ifdef USE_OPTIMA_SOLVER: optima_want_sensitivity, optima_dndb (+rows/cols) — **dn/db already exists on Optima path**;
  insBudgetTried (declared outside ifdef "to keep member position"; ABI-sensitive layout under ifdef!);
  GetControlConditions/RestoreControlConditions; CalculateEquilibriumStateOptima(…, referenceMode, runKinetics),
  CalculateEquilibriumStateHOP(…, warmNative), NativeOptimaFinish, NativeOptimaRescue, SetControlCondition_pH/Eh,
  ClearControlConditions, OptimaMaxMassBalanceResidual, HasAqueousPhase, CheckRedoxFree, PhaseHessianBegin/End,
  WarnCollapsedTwins, OptimaTwinLabel, phfTwin, phfFired, DetectSolventCollapseAndReseed,
  DetectPhaseCollapseAndReseed, LPFeasibilitySeed, PotentialSpaceFinish, LPGibbsDual, WorstPhaseStabilityViolation,
  OptimaReducedPreSolve. Then pTSolMod, CheckMassBalanceResiduals, ConvertGj_toUniformStandardState,
  PhaseSpecificGamma, HelmholtzEnergy, InternalEnergy.
  protected: optima_control_conditions, optima_disable_smoothing (flag read by SmoothingFactor = HIDDEN COUPLING),
  optima_hop_leg (flag read by dimreduce gate = HIDDEN COUPLING), prev_scfact, HopLegSplit hop_split, `MULTI pm;
  MULTI* pmp;` sizeN, AA, BB, arrL, arrAN (compressed A), phSolMod (TSolMod**), phSorpMod, phKinMet; uDD divergence
  arrays U_mean,U_M2,U_CVo,U_CV,ICNud; ipmKktRescues; SolveIpmAugmented; virtual GUI hooks (get_PAalp_PSigm,
  STEP_POINT, alloc_* x22, loadData, testTSyst, calculateActivityCoefficients_scripts, initalizeGEM_IPM_Data_GUI,
  multiConstInit_PN, GEM_IPM_Init_gui1/2) — **GEMSGUI subclasses TMultiBase** (TMulti in GUI overrides these).
  Methods grouped by source file: ipm_chemical.cpp (XmaxSAT_IPM2, DC_DualChemicalPotential, Set_DC_limits,
  TotalPhasesAmounts, DC_PrimalChemicalPotentialUpdate, DC_PrimalChemicalPotential, PrimalChemicalPotentials,
  KarpovCriterionDC, KarpovsPhaseStabilityCriteria, StabilityIndexes, DC_GibbsEnergyContribution, GX, ConvertDCC,
  getXvolume); ipm_chemical2.cpp (GasParcP, phase_bcs, phase_bfc, bfc_mass, CalculateConcentrationsInPhase,
  CalculateConcentrations, IS_EtaCalc, GouyChapman, SurfaceActivityCoeff); ipm_chemical3.cpp (SmoothingFactor,
  SetSmoothingFactor, CalculateActivityCoefficients, SolModCreate, SolModParPT, SolModActCoeff, SolModExcessProp,
  SolModIdealProp, SolModStandProp, SolModDarkenProp); ipm_chemical4.cpp (CalculateKinMet, RunKineticsStep, KM_*);
  ipm_main.cpp (GibbsEnergyMinimization, GEM_IPM, MassBalanceRefinement, MassBalanceReproject,
  EnergyDeterminacyCheck, ExcludeRedundantDCs/RestoreRedundantDCs, StrandedElementCheck, TraceGrowthCheck,
  trace_growth_stuck/note, phase_floor_parked, zeroed_floor_amount, warm_floor_seed, InteriorPointsMethod,
  AutoInitialApproximation, MassBalanceResiduals, OptimizeStepSize, DC_ZeroOff, DC_RaiseZeroedOff,
  ApplyFilloutBudget, FilloutBudgetValue, PhaseInsertionCeiling, RaiseDC_Value, MetastabilityLagrangeMultiplier,
  WeightMultipliers, MakeAndSolveSystemOfLinearEquations, DikinsCriterion, StepSizeEstimate, Restore_Y_YF_Vectors,
  RescaleToSize, SpeciationCleanup, PhaseSelectionSpeciationCleanup, PhaseSelect, GEM_IPM_InitialApproximation);
  ipm_simplex.cpp (SolveSimplex, SPOS, START, NEW, WORK, FIN, SystemTotalMolesIC, ScaleSystemToInternal,
  RescaleSystemFromInternal, MultiConstInit, GEM_IPM_Init). public tail: elementsOfInterest, ICIsNumericalTrace,
  ICIsDefaultSeed, `static void BeginCalculation()` (resets per-calculation thread/static state; GEMSGUI must call it
  = HIDDEN STATE contract), RefreshFromAmounts, PhaseFloorParked, SubFloorElementCheck.
- Free functions (trace): native_trace_file() (env GEMS3K_NATIVE_TRACE_FILE, opened once, append, process lifetime),
  native_trace_quiet(bool) nests, ipm_probe_file() (GEMS3K_IPM_PROBE), native_trace_run_header(pm,pa,mode),
  native_trace_decide(fmt,...), first_note_for_system(systemKey, noteKey) (thread-level memo), native_trace_run_result,
  native_cert_dualfree_reset/set/get (global three-valued state).
- Enums THERM_POTENTIALS, NUM_POTENTIALS, volume_code.
- Layout note: MULTI member order matters to binary I/O? (check ms_multi_file.cpp). TMultiBase layout matters to
  GEMSGUI ABI ("trailing member, so no existing offset moves").

### node.h (1636 lines) — read in full
- includes datach_api.h, activities.h, v_detail.h, v_service.h, ms_multi.h, ThermoFun/ThermoFun.h (unless
  NO_THERMOFUN). `class TNode { friend class TNodeArray; ...}`.
- private: dbr_file_name, ipmlog_error, current_input_set_name, current_output_set_name, load_thermodynamic_data
  (bool flag, cleared by Set_mLook/Set_PMc/Get_PMc(!)/Set_DMc — a *getter* Get_PMc mutates load flag l.1283).
- protected: `MULTI* pmm; DATACH* CSD; DATABR* CNode;` unique_ptr<ThermoFun::ThermoEngine> thermo_engine;
  thermo_json_string; `TMultiBase* multi_base` + `shared_ptr<TMultiBase> internal_multi`; shared_ptr<TActivity> atp;
  shared_ptr<TKinetics> kip; ACTIVITY* AiP; CalcTime, PrecLoops, NumIterFIA, NumIterIPM; floorParkedAmount
  (per-DBR-phase, read by Ph_Present); multi_ptr().
- public static loggers `node_logger`, `ipmlog_file`.
- API (coupled codes): GEM_init(lst) / GEM_init(dch,ipm,dbr,fun json strings); datach_to/from_string,
  databr_to/from_string, gemipm_to_string; system_id; input/output_set_name; code_error_IPM, description_error_IPM;
  [USE_OPTIMA_SOLVER] Set_pH_target, Set_Eh_target, Clear_ControlConditions, keepApiTarget,
  Get_ControlCondition_titrant; GEM_restore_MT (2 overloads), GEM_from_MT (6 overloads: 8,8a,8b,8c,8d,8e),
  GEM_from_MT_time, GEM_set_MT (obsolete), GEM_step_MT (writes pmm->ITau directly); GEM_read_dbr (file / json
  string); **GEM_run(bool uPrimalSol)** l.518; GEM_rescue_step(); SyncMultiWithDataBr(); struct TraceRegime +
  GEM_trace_regimes(factors, traceRel, tol, mode, ofInterest); GEM_set_elements_of_interest /
  GEM_elements_of_interest; [OPTIMA] GEM_run_ROP_or_AIA (inline, writes NodeStatusCH and re-enters GEM_run);
  GEM_CalcTime, GEM_Iterations, GEM_IterationsHOP; GEM_write_dbr (file / string); GEM_print_ipm;
  write_ThermoFun_format_stream (virtual); GEM_to_MT; GEM_to_MT_time; pCSD(), pCNode(), pActiv(), otherPMM();
  cTC,cTK,cP,cMs,cVs; setNodeHandle; Access_GEM_IMP_init; get_sizeTSolMod; get_ptrTSolMod(void*);
  name/index maps IC_name_to_xCH, DC_name_to_xCH (first match, trap 8), DC_name_to_xCH_map, xCH_to_DC_name,
  xCH_to_PH_name, xCH_to_IC_name, xCH_to_ccPH, Get_LsMod, Get_LsMdc, Ph_name_to_xCH, IC/DC/Ph_name_to_xDB,
  DC_xCH_to_xDB_map, IC/DC/Ph_xCH_to_xDB, check_IC/DC/Phase_xCH, IC/DC/Ph_xDB_to_xCH, Phx_to_DCx, PhtoDC_DCH,
  DCtoPh_DCH, PhtoDC_DBR, DCtoPh_DBR; `virtual packDataBr()` l.870, `unpackDataBr(bool)`; check_TP, check_grid_T/P,
  grid_P, check_grid_TP, gridTP, get_Ppa_sat, get_grid_index_Ppa_sat, Get_Psat; Set_DC_G0 (writes CSD grid),
  Set_aPH; DC_G0 (reads CSD->G0 grid, trap 4) l.945; Get_solver_G0 (pm.G0) l.954; Get_solver_L; Get_CertRank;
  DC_V0, DC_H0, DC_S0, DC_Cp0, DC_A0, DC_U0; DensArrayH2Ow, EpsArrayH2Ow, EpsH2Ow, DenH2Ow, VisH2Ow, EpsH2Og,
  DenH2Og; Ph_Volume, Ph_GibbsEnergy, Ph_Enthalpy, Ph_Entropy, Ph_HeatCapacityCp, Ph_Moles, Ph_Mole (obsolete),
  Ph_Present(xBR, relFloor), Ph_Mass, Ph_SatInd, Ph_BC (may allocate with new[] — caller frees!); Get_mIC, Get_pH,
  Get_GibbsEnergy (TotalGibbsEnergy, mutates pm), Get_pe, Get_Eh, Get_IC; Set_TK, Set_P, Get_P, Get_TK, Set_bIC,
  Get_uIC, Get_bIC, Set_dll, Set_dul, Set_nDC, Get_nDC, Get_muDC, Get_aDC, Get_cDC, Get_gDC, DCmm, ICmm, DCaJI;
  MULTI-direct: Set_IC_b (pmm->B, internal units!), Set_mLook, Set_PMc, Get_PMc, Get_sMod, Set_DMc, Get_DMc,
  IC_b, DC_n (pmm->X), DC_a, GEMSFIT: Get_IPc_IPx_DCc_indices, Get_NPar_NPcoef_MaxOrd_NComp_NP_DC, Get_aIPc,
  Get_aIPx, Get_aDCc, Set_aIPc, Set_aDCc, Set_Tk, Set_Pb; DC_c, DC_g, DC_lng, DC_mu, DC_mu0; databr_to_vtk (virtual);
  ipmLogError, clearipmLogError; load_all_thermodynamic_from_thermo, clear_ThermoEngine, load_ThermoEngine;
  databr_free. TActivity API: InitCopyActivities, setTemperature, setPressure, updateStandard* (5),
  updateThermoData, setSpeciation, updateConcentrations, initActivityCoefficients, updateActivityCoefficients,
  getIntegralPhaseProperties, updateChemicalPotentials, updateTotalGibbsEnergy, updateActivities,
  updateChemicalData.
- protected: rescueStep; ApplyDbrTargets(mode), dbrTarget_pH/Eh, dbrTargetWarned, static modeUsesTargets;
  titrEntered_, titrSolved_, restoreEnteredComposition, rememberTitration; apiTarget_pH/Eh;
  **GEM_run layers declared** l.1460-1474: GEM_run_single, WarnIfSystemTooLarge, GEM_run_rescue (pa_FailRescue),
  GEM_run_capped (pa_OptimaCallCap), GEM_run_dry_redo (pa_OptimaDryRetry), GEM_run_tries (pa_OptimaFastFirst),
  GEM_run_cold_retry(maxNudges) (pa_ColdRetryNudges), NudgedAnswerBalances, GEM_run_optima_cold_retry,
  GEM_run_aop_cgseed, GEM_run_single_legacy_retry. → the "named layers" split (brief §10.1) looks partly DONE in
  this head; check node.cpp whether they still re-enter GEM_run() with thread-local flags.
  allocMemory, freeMemory, databr_realloc, read_dbr_format_file, write_dbr_format_file(checkTP), databr_*_to_vtk,
  datach_copy, databr_copy, CheckMtparam, virtual init_into_gems3k, WarnSolidsCodedAqueous.

### datach.h (142) — read in full
- `DATACH` typedef struct: dims nIC,nDC,nPH,nPS,nDCs,nTp,nPp,iGrd,nAalp,nICb,nDCb,nPHb,nPSb,mLook; arrays nDCinPH, xic,
  xdc, xph; Ttol, Ptol, dRes1/2; TKval, Pval, Psat, A [nIC][nDC]... actually indexed A[xic + xdc*nIC] (DC-major),
  ICmm, DCmm, DD; grids denW, denWg, epsW, epsWg, G0, V0, S0, H0, Cp0, A0, U0 [nDC][nPp][nTp]; name lists
  std::vector<std::string> ICNL, DCNL, PHNL (C++ members inside a C-style typedef struct); ccIC, ccDC, ccPH.
- enums DATACH_STATIC_FIELDS, DATACH_DYNAMIC_FIELDS (I/O field indexes, positional contract for files).

### databr.h (266) — read in full
- `DATABR` struct: NodeHandle, NodeTypeHY, NodeTypeMT, NodeStatusFMT, NodeStatusCH, IterDone; TK, P, Vs, Vi, Ms, Mi,
  Gs, Hs, Hi, IC, pH, pe, Eh, Tm, dt; **#ifndef NO_NODEARRAYLEVEL**: Dif, Vt, vp, eps, Km, Kf, S, Tr, h, rho, al, at,
  av, hDl, hDt, hDv, nto (struct layout depends on a macro: ABI hazard between builds); arrays bIC, rMB, uIC, xDC,
  gam, dul, dll, aPH, xPH, vPS, mPS, bPS, xPA, bSP, amru, amrl, omPH; pHtarget, Ehtarget (trailing).
- `enum NODECODECH` l.137-196: 0..9 native + T_ERROR_GEM=9; AOP 10-13, SOP 14-17, ROP 18-21, HOP 22-25, SHP 26-29.
  Pattern NEED=base, OK=base+1, BAD=+2, ERR=+3.
- NODECODEFMT, NODETYPE, NODEINEX, DATABR_FIELDS (f_pHtarget/f_Ehtarget after VTK-only fields "so indices do not change").

### s_solmod.h (2456) — read in full
- Enums fluid_mix_rules (MR_* chars), dc_class_codes, tp_codes (char-coded closed sets; string/char dispatch).
- `struct SolutionData` l.83-137: the hand-over record MULTI -> model: counts (NSpecies, NParams, NPcoefs, MaxOrder,
  NPperDC, NSublat, NMoiet, NDQFpDC), Mod_Code, Mix_Code, and ~20 RAW POINTERS INTO MULTI arrays (DC_Codes, TP_Code,
  arIPx, arSM, arIPc, arDCc, arMoiSN, arSitFr, arSitFj, arGEX, arDQFc, arPparc, arWx, arlnGam, arlnDQFt, arlnRcpt,
  arlnExet, arlnCnft, arCTermt, arVol, aphVOL (pm.FVOL[k]), aphXF (pm.XF[k])), T_k, P_bar.
  => models do not own their inputs/outputs: they alias MULTI. Any MULTI realloc invalidates every TSolMod.
- `class TSolMod` l.140-309: protected ModCode, MixCode, DC_Codes, PhaseName, DC_Names, NComp, NPar, NPcoef, MaxOrd,
  NP_DC, NSub, NMoi, NDQFpc, aIPx, R_CONST, Tk, Pbar, aIPc, aIP, aDCc, aGEX, aPparc, aMoiSN, aSitFR, DQFcf, x, aVol,
  phVOL, phXF; results Gex..Uex, Gid..Uid, Gdq..Udq, Grs..Urs; lnGamConf, lnGamRecip, lnGamEx, lnGamDQF, CTerm,
  lnGamma (pointers into MULTI); multisite arrays y, mn, mns, fjs, aSitFj; alloc/init/free_multisite, IdealMixing,
  ideal_conf_entropy, return_sitefr, retrieve_sitefr; `static constexpr bool kReportEosMixtureVolume = false`;
  set_phVOL_from_molar. public: `static shared_ptr<spdlog::logger> solmod_logger`; ctors TSolMod(SolutionData*),
  TSolMod(NSpecies, Mod_Code, T, P); virtual PureSpecies, PTparam, MixMod, ExcessProp(double*), IdealProp(double*),
  StandardProp(double*), Set_Felect_bc; UpdatePT; testSizes; Get_lnGamma; to_text_file, to_json_file,
  to_json_stream. Virtual interface = the natural plug-in seam for activity models (already polymorphic).
- Subclasses (one per model; ~30): TIdeal; EOSPARAM (helper) + TCGFcalc (Churakov-Gottschalk EoS; ctor
  (sd, aphWGT, arX)); TPRSVcalc, TSRKcalc, TPR78calc (cubic EoS, near-duplicate classes: same members, AB(),
  FugacityPT, FugacityPure, Cardano, MixParam, FugacityMix, FugacitySpec, ResidualFunct, MixingWaals/Const/Temp —
  DUPLICATED LOGIC x3); TCORKcalc, TSTPcalc (fluids, Cardano again); TVanLaar, TRegular, TRedlichKister, TNRTL,
  TWilson; TBerman (multisite, reciprocal; ctor (sd, G0) — takes pm.G0 pointer); TCEFmod (Calphad CEF, ctor (sd,
  G0)); TMBWmod (modified Bragg-Williams, "TEST MODEL", ctor (sd,G0)); aqueous: TSIT, TPitzer (ctor (sd, arM, arZ,
  dW, eW): molalities, charges, water density/eps pointers; getIc/getIa/getIn linear search in hot path; p_sum),
  TEUNIQUAC (fixed-size debug arrays gammaDH[200], gammaC[200], gammaR[200] — overflow if NComp>200), TELVIS
  (#ifdef ELVIS_SPEED block declares `U` twice: dead code that cannot compile; debug arrays [200]), THelgeson,
  TDavies, TLimitingLaw, TDebyeHueckel, TKarpov, TShvarov (EDH variants: duplicated IonicStrength/BgammaTP/IonsizeTP/
  Gfunction/GShok2 in THelgeson, TKarpov, TShvarov = DUPLICATED LOGIC x3); TModOther (hard-coded minerals Amphibole1
  .. Talc), TMargules, TSubregular, TGuggenheim; TSCM_NEM (surface complexation, EDLmod, SCM_conf_entropy).
- All models: raw new/delete arrays (alloc_internal/free_internal), `long int` return codes, no exceptions.

### s_kinmet.h (623) — read in full
- `MAXDCNAME_ = 16, MAXPHNAME_ = 16, MAXSYMB_ = 4` redeclared (duplicates of m_const_base values; "see v_mod.h").
- enum kinmet_controls (char codes; note duplicate char values in one enum: KM_MOD_PAL_='P', KM_PRO_NUPR_='P',
  KM_NUCL_SSMP_='P', KM_LNK_PVOL_='P' — codes are positional per field, enum is not a type), affin_term_op_codes.
- `struct KinMetData` hand-over record: codes, NComp, nlPh, nlPc, nPRk, nSkr, nrpC, naptC, nAscC, nFaceC, T, P, kTau,
  kdT, IS, pH, pe, Eh, nPh, mPh, vPh, sAPh, LaPh, OmPh, sFact, sSA, sgw, sgg, rX0, hX0, sVp, sGP, nPul, nPll; raw
  pointers into MULTI: arlPhc, arfeSAr, arrpCon, arapCon, arAscp, SM_, arDCC_, arPhXC_, arocPRk_, arxSKr_, arym_
  (molalities all species), arla_ (lg activities), arxp_, armp_, arvp_, arasp_ (phase arrays), arnx_, arnxul_ (DUL
  output, direct write), arnxll_ (DLL output), arWx_, arVol_.
- `struct TKinReact` (per parallel reaction): constants ko, Ko, Ap, Ea, bI, bpH, bpe, bEh, pPR, qPR, mPR, uPR, OmEff,
  GamN; results arf, cat, aft, k, K, rPR, rmol.
- `class TKinMet`: static kinmet_logger; ~60 state fields; arPRt; spcfu/spcfl split coefficients; kTot,gTot,rTot,vTot;
  particle geometry; virtual PTparam, RateMod, RateInit, SplitMod, SplitInit, SetMetCon, UpdateFSA, GetModFSA,
  UpdatePT, UpdateTime, GetMdTime, testSizes, PRrateCon, CorrSpecSurfArea. Derived: TMWReaKin (SSReaKinInit/Mod),
  TUptakeKin (uptake, Thien-Kulik-Curti), TIonExKin, TAdsorpKin, TNucleKin (TBD stubs).
- writes DUL/DLL of MULTI directly through arnxul_/arnxll_ (kinetics changes the solver's bounds in place).

### s_sorpmod.h (416) — read in full
- "Stub declaration of new versions of SCMs": SorptionSiteData, TSurfSiteMod, SorptionData, TSorpMod (virtual
  SorptionSpecies, PTparam, IsothermMod, ElstatMod; UpdatePT; testSizes), TNEMcalc (declares PureSpecies, MixMod,
  ExcessProp, IdealProp, SRKCalcFugPure(!) — copy-paste of the SRK interface). Check s_sorpmod.cpp: likely mostly
  unused (old sorption is in ipm_chemical2.cpp SurfaceActivityCoeff / GouyChapman). Candidate DEAD CODE.

### ipm_optima.h (115) — read in full
- `struct EqControlCondition`: name, target, tolerance(-1 → defaultToleranceFn), std::function defaultToleranceFn,
  stoich vector<(IC row, coef)>, std::function fixedGradientFn(target), achievedValueFn(achievedMuj),
  needsAqueousPhase, slot (L+k), outputs titrantAmount (REAL moles), achievedValue, targetMet. Plain C++ so a
  non-Optima build can hold the (empty) list. Closures capture indices and read G0[]/T at call time (hidden aliasing).

### nodearray.h (794) — read in full
- ProcessProgressFunction, PhaseDataLogFunction (std::function); enum PTCODE; struct LOCATION; struct
  TestModeGEMParam {mode, useSIA, step, cdv, cez}.
- `class TNodeArray : std::enable_shared_from_this<TNodeArray>` — **PRIVATE inheritance** (class default), so
  shared_ptr cannot see the base: weak_this is never set and `getptr()` (shared_from_this) throws bad_weak_ptr.
  (NON-OBVIOUS FACT candidate; verify nobody calls getptr()).
- protected: internal_Node (shared_ptr<TNode>), calcNode (TNode*), NodT0/NodT1 (DATABR** arrays: previous/current
  time), anNodes, sizeN/M/K, size, grid, tcNode, iaNode, warnedOutside; allocMemory/freeMemory; getGrid,
  isLocationInNode; setNodeArray (2), checkNodeArray, CopyWorkNodeFromArray, MoveWorkNodeToArray, NeedGEMS,
  SmartMode, ErrorGEMsMessage, virtual CalcIPM_Node, allocNewDBR (raw new DATABR), virtual CalcNodeServer,
  pVisor_Message; InitNodeArray; protected ctors TNodeArray(nNod), (N,M,K).
- public: **`static TNodeArray* na`** (GLOBAL singleton pointer used by all node0_*/node1_* macros) l.205; getptr;
  getCalcNode; pCSD; gridTP; DCaJI; ICmm; DCmm; Ph_xDB_to_xCH; static create(nNod)/(N,M,K); Set_DC_G0; LinkToNode;
  SaveToNode (calls calcNode->packDataBr()); genGEMS3KInputFiles; GEMS3k_write_dbr (3); GEMS3k_read_dbr (2); iNode,
  indN/M/K; nNodes, SizeN/M/K, pNodT0/1, piaNode, ptcNode; RunGEM (node, ndx, Mode, arr), RunGEM(N,M,K...),
  RunGEM(Mode, nNodes, nodeArray, nodeFlags, retCodes) ("optional openmp parallelization"); CalcIPM_One,
  CalcIPM_List (virtual); GEM_init (lst files; json strings); GEM_print_ipm; get_mPH, get_vPH, get_bPH;
  CopyNodeFromTo; logDiffsIC, logProfileAqIC, logProfileTotIC, logProfilePhMol, logProfilePhVol, logProfileAqDC;
  SetGrid, FindNodeFromLocation, GetNodeSizes, GetNodeLocation, GetSize, GetNodeMass, warnParticleOutsideOnce,
  resetOutsideWarning, MoveParticleMass, databr_to_vtk.
- ~120 macros node0_*/node1_* dereferencing TNodeArray::na (global). Layout = array of structs (DATABR* per cell),
  each DATABR holding its own heap arrays: worst case for batching/GPU (pointer chasing, no contiguous SoA).

### activities.h (400) — read in full
- `ACTIVITY` struct: copy of the MULTI fields needed for activity evaluation (TK,P,RT,ln5551,IC,pH,pe,Eh,Pc,Tc,
  LsMod..PhLin, denW/denWg/epsW/epsWg[5], PMc, DMc, MoiSN, SitFr, lPhc, DQFc, sorption arrays, H0,S0,Cp0, tpp_G, G0, G,
  lnAct, lnGam, lnGmo, lnDQFt.., Gamma, lnGmf, lnGmM, lnScalT.., VPh..UPh, sMod, dcMod, SM, SF, ICC, DCC, PHC, DCCW,
  MM, FWGT, YOF, N,NR,L,Ls,LO,FI,FIs,K2,L1, A, Awt, XF, YF, XFA, YFA, X, Y, Fx, Wx, F, F0, U, Falp, GamFs, fDQF, MU,
  EMU, NMU, Y_la, Y_m, Pparc, EZ, Vol, FVOL, BF, BFC, cutoffs, ICmin, aqsTail, lowPosNum, logXw, logYFk, YFk).
- `class TActivity` (2014-2018 "alternative Activity API", DK/AL): act, csd, cnd, cno, AA/BB/arrL/arrAN, phSolMod,
  phSorpMod, load, IT, IIM, AG, DGC, FitVar; methods mirror TMultiBase one-to-one (DC_DualChemicalPotential,
  TotalPhasesAmounts, DC_PrimalChemicalPotentialUpdate, DC_PrimalChemicalPotential, PrimalChemicalPotentials,
  KarpovCriterionDC, KarpovsPhaseStabilityCriteria, StabilityIndexes, DC_GibbsEnergyContribution, GX, ConvertDCC,
  getXvolume, phase_bcs, phase_bfc, bfc_mass, CalculateConcentrationsInPhase, CalculateConcentrations,
  CalculateActivityCoefficients, SolModCreate.. SolModDarkenProp) + generic setTemperature.. updateChemicalData.
  => DUPLICATED LOGIC: s_activity*.cpp re-implements ipm_chemical*.cpp. Not on the solve path (verify).

### kinetics.h (212) — read in full
- `KINETICS` struct (copy of MULTI kinetics fields) and `class TKinetics` (csd, cnd, cno, phKinMet, load) with KM_*
  methods mirroring TMultiBase's ipm_chemical4.cpp = DUPLICATED LOGIC (node_kinetics.cpp). header guard named
  Kinetics_H_ but file comment says activities.h (copy-paste).

### m_const_base.h (379) — read in full
- Constants: S_OFF/S_ON/S_REM/A_NUL, MST=6, DFCN=6, MAXICNAME=6, MAXSYMB=4, TDBVERSION=64, MAXDCNAME=16, MAXPHNAME=16,
  EQ_RKLEN=58, MPP_* (mixed-phase property columns), MIXPHPROPS=5. Enums (all char-coded): DC_LIMITS, SolDCodes
  (DCCW: U,I,S,W), solmod_switches (SM_* model codes 'I','B','$','#',... + MR_* mixing rules — **one enum mixes
  index values 0..11 and char codes**), IC_CLASSES, DC_CLASSES, SORPTION_DC macro, PH_CLASSES, sorption_control,
  ph_kinmet_controls (duplicate of s_kinmet.h kinmet_controls with no trailing underscore), SPPA_UNITS.
  => DUPLICATED enums (s_solmod.h re-declares MR_*_, s_kinmet.h KM_*_). Strings/chars dispatch everywhere.

### verror.h (88) — read in full
- `extern shared_ptr<spdlog::logger> gems_logger`; `_GEMIPM_version_stamp`; `struct TError {mess,title}` (virtual
  dtor), `TFatalError : TError`; `[[noreturn]] Error(title,msg)` (throws TError), `LogError`, `ErrorIf`.
- **`extern thread_local int g_errorsHandledByCaller`** (HIDDEN STATE: Error() logs at debug instead of error while > 0).
- `SetCalculationCancel(bool)` / `CalculationCancelled()` — cross-thread cancel flag (added in head commit 2ccbef2):
  global (check whether atomic and process-wide: one GUI cancel stops all threads' calculations?).

### v_detail.h (232) / v_service.h (205) — read in full
- v_detail: is<T>, approximatelyEqual/essentiallyEqual/definitelyGreater/LessThan, approximatelyZero, noZero (eps =
  numeric_limits::epsilon, absolute!), isAllZero, FLOAT_EMPTY/DOUBLE_EMPTY = numeric_limits::min() (sentinel
  values), IsDoubleEmpty, FLOAT/DOUBLE_INF*/NAN sentinels (DOUBLE_NAN = 1.797693e308 — finite "nan"), InfMinus/Plus,
  Nan templates (Nan = max()), is_minusinf/is_plusinf/is_nan, floating_point_to_string (precision 15!),
  string_to_floating_point, string2value, value2string. Note: text output precision 15 digits → round trip of a
  double through text is lossy (17 needed). Relevant to "freeze" via files and to any I/O refactor.
- v_service: char_array_to_string, strip, replace, replaceall, replace_json_inf_nan, split, regexp_extract_string,
  extract_string_json, extract_int_json, f_getline, u_makepath, replace_all, u_create_directory, u_splitpath,
  u_getpath, ROUND, fillValue, copyValues (overloads incl. float<->double, short->long), ltrim/rtrim/trim,
  name_for_message, to_string(arr,size).

### io_template.h (347), io_keyvalue.h (254), io_nlohmann.h (248), io_simdjson.h (283) — read in full
- Already a policy-based I/O design: namespace io_formats; `struct outField {name, alws, readed, indexation,
  comment}` = the field descriptor table; `TRWArrays` (findFld, setAlws/NoAlws, reset, testRead);
  `TPrintArrays<TIO>` (writeField, addField, writeArray (several), writeArrayF, write_sel_array);
  `TReadArrays<TIO>` (findNext → returns field index for a `switch`, readNext, readArray overloads).
  TIO drivers: KeyValueWrite/KeyValueRead (`<name> values`, 40 per line), NlohmannJsonWrite/Read (#ifdef
  USE_NLOHMANNJSON), SimdJsonWrite/Read (#ifndef USE_NLOHMANNJSON; pimpl SimdJsonImpl). JSON choice is compile-time.
- Writing doubles: specialisations writeValue(double)/add_value(double) in .cpp (check precision there).
- Field order in files = order of `outField` table rows + `MULTI_DYNAMIC_FIELDS` enum; reader dispatch by
  enum index in a switch (ms_multi_format.cpp). This is where trap 1's "descriptor row / writeField / read case /
  two counts" live.

### jsonconfig.h (246), gdatastream.h (132), num_methods.h (180), gems3k_impex.h (222), gems3k_timing.h (92),
### gems3k_version.h (7), lazy_file_sink.h (76), datach_api.h (188) — read in full
- jsonconfig.h: home_dir(); JsonConfigSection (nlohmann or simdjson impl), TJsonConfig (section/value/
  value_or_default/value_must_exist), `class GemsSettings : TJsonConfig` with STATIC members settings_file_name,
  data_logger_directory, logger_section_name, default_gems3k_loggers, gems3k_logger_pattern, log_file_size,
  log_file_count, log_thermodynamic; gems3k_clear_loggers/update_loggers/update_log_level; `gemsSettings()` global
  singleton accessor. (GLOBAL STATE: logger config.)
- gdatastream.h: GemDataStream (binary, endian swap) operator<</>> per type; writeArray/readArray. `long` is 8 bytes on
  Linux, 4 on Windows → binary files not portable across platforms (operator<<(long)).
- num_methods.h: getNpoints, getStep, LagranInterp (3 overloads float/double), LagranInterp1D, quot/prod2/prod3
  derivative helpers, GoldenSection(+Two) via function pointer.
- gems3k_impex.h: class GEMS3KGenerator {enum IOModes f_key_value, f_binary, f_json, f_thermofun, f_kv_thermofun;
  static default_type_f; ext(), gen_*_name, ctor reading *-dat.lst; get_*_path; load_dbr_lst_file ...}.
- gems3k_timing.h: namespace gems3k_timing: enum Part (Call, Native, OptimaCall, OSolve, Obj, ObjFxx, FdMix, FdPure,
  FdBlock, ObjRestore, ActCoef, Psf, Tpd, Wpsv, PreSolve, Level); env GEMS3K_TIMING_FILE (function-static cached);
  `inline thread_local Acc acc`; RAII Scope (steady_clock); flush(requested, status, iterations) writes a TIMING line
  per call. **Existing per-part profiler — use it for 03-profile.** Zero cost when unset (one cached bool).
- gems3k_version.h is GENERATED INTO THE SOURCE TREE by root CMakeLists.txt l.150 configure_file(... GEMS3K/
  gems3k_version.h): every configure rewrites a tracked file (hash, branch, ThermoFun version). The committed copy
  is stale (hash 3652086, branch develop_optima2, ThermoFun 0.6.0). NON-OBVIOUS FACT. (Restore after build.)
- lazy_file_sink.h: gems3k_detail::lazy_rotating_file_sink (opens the log file at first message).
- datach_api.h: namespace dbr_dch_api: datach/databr to/from text (TIO templates), to/from string, format streams,
  binary to/from file, gridTP, is_Psat_grid, grid_P, check_TP, change_TP, check_grid_T/P/TP, datach_reset/realloc/
  free, databr_reset/realloc/free_internal. (Free functions over C structs: good seam for the data module.)

### tnt*.h / jama*.h (vendored NIST TNT + JAMA, ~1600 lines) — skimmed (third-party, unchanged)
- Array1D/Array2D (ref-counted i_refvec), JAMA::Cholesky, JAMA::LU. Used ONLY by ipm_main.cpp: l.37-38 include;
  l.4926 `JAMA::LU<double> lu(KKT)` (augmented KKT), l.5228 `JAMA::Cholesky<double> chol(A)` and l.5251 JAMA::LU
  (normal equations in MakeAndSolveSystemOfLinearEquations). Eigen used in ipm_main.cpp (l.34 <Eigen/QR>,
  l.4948-4952 `thread_local Eigen::MatrixXd S; VectorXd sb, su; ColPivHouseholderQR qr` = thread-local scratch,
  HIDDEN STATE but harmless reuse) and in ipm_optima.cpp. Two LA libraries in one solver; JAMA loops fix their own
  FP order (Stage A must keep them; Stage B may replace by Eigen LLT/PartialPivLU = FP-order change).

### node.cpp (2477) — part 1 (l.1-1680)
- l.47 static `TNode::node_logger = spdlog::stdout_color_mt("tnode")` (static-init order dependency on spdlog
  registry); l.49 ipmlog_file created once in ctor via `std::once_flag once_flag` (GLOBAL, l.60) → lazy file sink
  "ipmlog.txt" in GemsSettings dir; also installs Optima::setWarningHandler (process-global).
- ctor TNode(): allocMemory() = new DATACH, new DATABR, `internal_multi.reset(new TMultiBase(this))`, set_def,
  pmm = GetPM(), atp (TActivity), kip (TKinetics) **always allocated though not used on the solve path**.
  copy ctor: datach_copy, databr_copy, copyMULTI, control conditions, Pai/Tai/Fdev.
- **THREAD-LOCAL FLAGS (node.cpp l.203-228)**: g_optimaLineSearchRetry, g_optimaLSEscapeOff, g_optimaEscRetryNext,
  g_optimaCgSeedArmed, g_optimaTargetWarned, g_optimaIterBudget (-1), g_optimaFixAccept, g_optimaQuickTry,
  g_optimaPreSolveEmpty, g_lvlFullCalls, g_optimaCapCall, g_optimaCapAllowance, g_optimaDryRedoCall,
  g_optimaDryFailed, g_optimaDryOff; extern from ipm_main.cpp: g_aqueousAbsent, g_ipmLevelUsed, g_ipmNoProgArmed,
  g_ipmNoProgFired, g_ipmLevelOff; also `extern thread_local bool g_tpdStrongSearch` (l.1000, declared locally in a
  block!). Plus function-local `static thread_local bool bigWarned` (l.267), `warned` (l.1179).
- `TMultiBase::BeginCalculation()` (l.245, static, defined in node.cpp!) resets only g_optimaTargetWarned,
  g_optimaPreSolveEmpty, g_lvlFullCalls. It is called from GEM_run_tries() (l.944) i.e. once per *attempt* (every
  rescue/cap/dry re-entry), not once per GEM_run(). NON-OBVIOUS FACT: "per GEM_run()" comments (l.207, l.222) are
  inaccurate; the reset scope is the tries layer.
- GEM_run(bool) l.297: WarnIfSystemTooLarge → ApplyDbrTargets(requested) → restoreEnteredComposition → bBefore copy →
  timing Scope(Call) → `g_errorsHandledByCaller++` → GEM_run_rescue() → (catch TError: LogError, rethrow) →
  log ipmlog_error if failed → rememberTitration → gems3k_timing::flush.
  NOTE: GEM_run() itself no longer checks NEED_* — GEM_run_single() does (l.1438-1456: returns the old status for a
  non-NEED code, trap 5).
- Layer stack (outermost first): GEM_run → GEM_run_rescue (pa_FailRescue, l.464-783) → GEM_run_capped
  (pa_OptimaCallCap, l.790-904) → GEM_run_dry_redo (pa_OptimaDryRetry, l.908-937) → GEM_run_tries (l.940-1185:
  BeginCalculation; pa_OptimaFastFirst quick try (AOP only); lsretry/escretry; IpmLevelStart retry (natlvlretry);
  GEM_run_cold_retry (only when FailRescue<=0); GEM_run_optima_cold_retry (SOP/SHP); cgseed (AOP);
  WarnCollapsedTwins; aqgone note) → GEM_run_single (l.1410). Side entries: GEM_run_single_legacy_retry (l.1192),
  GEM_run_aop_cgseed (l.1233), GEM_run_optima_cold_retry (l.1249), GEM_run_cold_retry (l.1315).
  The layers call each other directly (brief §10.1 split is DONE), but they still pass state through the
  thread-local flags above (re-entry guards g_optimaCapCall, g_optimaDryRedoCall; budget g_optimaIterBudget,
  g_optimaCapAllowance; mode switches g_optimaFixAccept, g_optimaLSEscapeOff, g_optimaLineSearchRetry,
  g_optimaCgSeedArmed, g_ipmLevelOff, g_optimaDryOff).
- Rescue (pa_FailRescue) l.464: only six modes; snapshot D0; steps (1) warm from prev answer, (2) nudges
  bIC[i]*(1 ± k·1e-15) with alternating sign per IC (i odd: -), (3) walk bPrev→b0 RescueWalkSteps steps with ≤3
  halvings, (4) LocalSettings + looser settings (DK*10, OptimaTol 1e-6, DcFloor 1e-12, Tau 1e-12) then finish,
  (5) other solver warm from best pool entry, (6) best balanced attempt as BAD. Pool of DATABR snapshots.
  Traces suppressed during rescue (native_trace_quiet(true)) → DECIDE rescue lines still written? (decide may ignore
  quiet; check node_trace.cpp). Sums ITF/ITG/IterDone/CalcTime across attempts.
- SNAPSHOT PATTERN duplicated ≥12 times: `DATABR* d=new DATABR; databr_reset(d,1); databr_realloc(CSD,d);
  {DATABR* live=CNode; CNode=d; databr_copy(live); CNode=live;}` — swaps the member pointer CNode to reuse
  databr_copy() (copies INTO CNode). Allocation per call: GEM_run_tries allocates lvlBefore on EVERY native/HOP/SHP
  call when pa_IpmLevelStart>0 (default 2), and lsBefore on every Optima call with OptimaLineSearch>0 (default 1.5)
  → heap traffic per solve (profile).
- GEM_run_capped: cap only for AOP/SOP without kinetics; nudges + warm finish; last try without cap.
- GEM_run_tries quick try: LocalSettings + bp->DG=1e-5, OptimaDcFloor 1e-16, Tau 1e-16, LS 0, StallWindow 0,
  DimReduce -1, FDHessianDelay 0; g_optimaFixAccept; g_optimaQuickTry reset; fastcutaccept uses
  AbsentPhaseWorstDrive + PhaseHessianTwins + NativeTpdPhase with g_tpdStrongSearch.
- GEM_run_single l.1410: CalcTime=0..; check_TP; NEED_GEM_SIA/SOP → pNP=1, unpackDataBr(uPrimalSol);
  AIA/AOP/ROP/HOP/SHP → pNP=0 (dt>0 forces uPrimalSol=true), unpackDataBr; else return status (trap 5);
  native_trace_run_header(*pmm, pa, mode); kinetics setup kTau/kdT/ITau/pKMM; hop_split reset; dispatch:
  HOP/SHP → NoProgArm + CalculateEquilibriumStateHOP(.., SHP); AOP/SOP/ROP → CalculateEquilibriumStateOptima(..,
  ROP); else native: NoProgArm(NativeOptimaFinish>0) + CalculateEquilibriumState() with catch TError E06IPM/E07IPM
  → NativeOptimaRescue; ipmnoprog retry cold; then NativeOptimaFinish. Without USE_OPTIMA_SOLVER: warn + fall back.
  Then packDataBr(); IterDone; testMulti() → BAD_/OK_ by NodeStatusCH (as rewritten by packDataBr!);
  native_trace_run_result. catch(TError) → ERR_* WITHOUT packDataBr (trap 6, l.1646-1663); std::exception →
  T_ERROR_GEM.
  The status mapping chain (NEED→OK/BAD/ERR) is repeated 3 times as if/else ladders (l.1597-1662) though the
  enum is regular (NEED+1/+2/+3): candidate for one function.

### node.cpp part 2 (l.1680-2477)
- GEM_CalcTime (reads pmm->K2, ITF, ITG), GEM_Iterations (PrecLoops, NumIterFIA, NumIterIPM), GEM_IterationsHOP.
- **packDataBr()** l.1730-1828 (virtual): mode reset l.1738-1758 (ROP/HOP/SHP kept; AOP/SOP by pmm->pNP; ELSE →
  NEED_GEM_AIA/SIA by pNP: trap "new mode must be listed here" CONFIRMED). Reads MULTI: pNP, TCc, Pc, IT (IterDone =
  pmm->IT), VXc, FX (checks kTotalGibbsEnergyUnset), HXc, IC, pH, pe, Eh, MBX; calls
  `ReportAbsentPhaseIndexes()` inside try{}catch(...){} (swallows all errors; MUTATES pm.Falp then copies it);
  XF, Aalp, Falp, PhaseFloorParked→floorParkedAmount; FVOL, FWGT, XFA, PUL, PLL (nPSb), BF (bPS), X, Gamma, DUL,
  DLL (nDCb), B, C, U, BFC (nICb). Writes DATABR: NodeStatusCH, TK, P, IterDone, Vs, Gs, Hs, IC, pH, pe, Eh, Ms, xPH,
  aPH, omPH, vPS, mPS, xPA, amru, amrl, bPS, xDC, gam, dul, dll, bIC, rMB, uIC, bSP. Note **bIC written back from
  pmm->B** (titration results flow out this way; rescaling must be undone before this point — trap 3).
- **unpackDataBr(uPrimalSol)** l.1838-1950: CheckMtparam() (T/P change → reload thermo flag); kTau, kdT, TCc, Tc, Pc,
  VXc; DUL/DLL (+ check DUL<DLL → Error); B (+ check B < pa->DB for all but the LAST IC → Error; assumes charge is
  last); Aalp, Falp (omPH), XF; if uPrimalSol: IT=IterDone, MBX, IC, Eh, **Y** = xDC, PUL/PLL, XF=YF=Ph_Mole(ii),
  FVOL, FWGT, BF, XFA=YFA, C, U, Gamma, lnGam via PhaseSpecificGamma(.., 1).
- GEM_from_MT_time / GEM_set_MT (duplicate bodies); stray free-function declaration `double GEM_to_MT_time(...)` at
  l.1983 (dead; the member is defined elsewhere); GEM_restore_MT x2, GEM_to_MT (p_aPH = aPH*Ph_Mass(ii)),
  GEM_from_MT x6 (8, 8a convolves xDC into bIC, 8b copies xDC/gam else SIA→AIA, 8c, 8d, 8e) — 6 near-duplicate
  bodies; dead `useSimplex=false` variables. SIA→AIA downgrade only for native SIA, not SOP/SHP (inconsistency).

### node2.cpp (1353, LF) — read in full
- get_Ppa_sat, get_grid_index_Ppa_sat, Get_Psat (LagranInterp1D order 3); get_ptrTSolMod, Access_GEM_IMP_init
  (calls GEM_IPM_Init), get_sizeTSolMod.
- Name lookup: IC_name_to_xCH / **DC_name_to_xCH (first match, l.123-133 — trap 8)** / Ph_name_to_xCH: linear scans
  that log `gems_logger->warn` on miss; DC_name_to_xCH_map (l.135) loops k < CSD->nPHb but indexes DCH arrays
  nDCinPH[k]/PHNL[k] — mixes DBR count with DCH index space (wrong when nPHb < nPH); DC_xCH_to_xDB_map.
  IC/DC/Ph_xCH_to_xDB: linear search over xic/xdc/xph + warn (O(n) per call; used inside loops in Ph_* below).
- check_IC/DC/Phase_xCH; Phx_to_DCx, DCtoPh_DCH, PhtoDC_DCH, DCtoPh_DBR, PhtoDC_DBR.
- Set_DC_G0 (writes CSD->G0 grid; load_thermodynamic_data=false), Set_aPH, **DC_G0 (CSD->G0 grid, LagranInterp
  order 6) — trap 4**, DC_V0 (5), DC_H0 (5), DC_S0 (4), DC_Cp0 (3), DC_A0 (5), DC_U0 (5), EpsArrayH2Ow,
  DensArrayH2Ow (falls back to Psat grid index), EpsH2Ow, DenH2Ow, EpsH2Og, DenH2Og (`VisH2Ow` declared, check def).
- Ph_Volume (single-component: DC_V0 × xDC), Ph_GibbsEnergy (new[]/delete[] each call; Σ bPS·uIC·RT·mol),
  Ph_Enthalpy/Ph_Entropy/Ph_HeatCapacityCp (standard-state sums, no mixing), Ph_Moles, Ph_Present (floorParkedAmount
  exact-equality test), Ph_Mole, Ph_Mass, Ph_SatInd (omPH), Ph_BC (new[] if null — leak by contract).
- Get_muDC (Σ A·uIC), **Get_aDC uses DC_G0 (grid), not pm.G0 → with ThermoFun it compares a dual potential from
  pm.G0 with a grid G0 (same trap-4 shape)**; DC_mu0 also returns DC_G0 (grid). Get_cDC (assumes aqueous phase is
  DBR phase 0 and its DCs are first in DCH), Get_mIC (uses pmm->LO — a MULTI index — as a DBR DC index in DCaJI).
  Get_pH/pe/Eh/IC; DC_a (10^Y_la); GEMSFIT accessors Get_IPc_IPx_DCc_indices (ip_DCc uses LsMdc[k] not
  LsMdc[k*3] — inconsistent with Get_NPar..., Set_aDCc/Get_aDCc also use LsMdc[index_phase]), Set/Get_aIPc,
  Get_aIPx, Set/Get_aDCc, Set_Tk, Set_Pb; DC_c, DC_mu (Fx/RT), DC_mu0.

### node_trace.cpp (208) — read in full
- GEM_trace_regimes(): picks trace ICs (b <= traceRel·Σb or elementsOfInterest), snapshot DATABR (same pointer
  swap) + full MULTI snapshot (`new TMultiBase(this)`, copyMULTI) restored with copyMULTIData(.., false) because
  TSolMod objects point into live MULTI arrays (OWNERSHIP FACT: restoring MULTI by reallocation would dangle every
  TSolMod). Re-enters GEM_run() (full layer stack) 1+factors times; classifies LINEAR/SATURATED/BOUNDARY/NONLINEAR/
  STRANDED/FAILED. GEM_set_elements_of_interest / GEM_elements_of_interest (multi_base->elementsOfInterest).

### ipm_main.cpp (5814, LF) — part 1 (l.1-1370)
- l.49 static TMultiBase::ipm_logger = stdout_color_mt("ipm"); `#define uDDtrace false`.
- Trace infrastructure lives HERE (not in a trace module): native_trace_file() l.57 (function-local static FILE*,
  env GEMS3K_NATIVE_TRACE_FILE, fopen "a", never closed); `static thread_local int native_trace_quiet_depth` l.72
  (quiet only suppresses RUN/BULK/SET/EFF and KEY/CERT; DECIDE always written); ipm_probe_file() l.80
  (GEMS3K_IPM_PROBE); `static thread_local long s_ipmKktFallbacks, s_ipmQrSortPivot` (first-only DECIDE);
  first_note_for_system() l.98 (thread_local set, cleared at 256 entries); native_trace_decide() l.105 (vfprintf +
  fflush per record: one syscall per DECIDE when tracing on).
- native_trace_run_header() l.124-283: RUN (mode name switch on literal 1,5,10,14,18,22,26), BULK, **SET (one
  printf with ~125 format specs in BASE_PARAM order: trap-1 site #7)**, EFF (dimreduce, earlystability, netresume,
  presolve budget, callcap — calls the same resolver functions as the solver: trap 13 single source).
- **`static int native_cert_dualfree_resolved` l.316 — process-wide, NOT thread_local** (data race when nodes run
  in parallel threads; CERT line of one thread can carry another thread's value). NON-OBVIOUS FACT.
- Cert instruments (report-only): CertSymEigMin (cyclic Jacobi), CertPrimalPotentials (mirror of
  PrimalChemicalPotentials with locals), cert_species_scorable, cert_dc_floor, cert_box_state, CertKktMax,
  CertDualFreeDirs (MGS 1e-8), CertCondFromGram, CertJacobiScale, CertRank, CertCurvMin (FD Hessian per phase;
  save/restore 11 arrays by copy; comment: **CalculateActivityCoefficients(LINK_UX_MODE) is NOT idempotent — it
  accumulates pm.lnGmo and blends pm.F0 through smoothing FitVar[3]**, and TSolMod objects hold their own internal
  composition that no array copy restores → re-run at base X to "re-seed"), CertStabTPD (successive substitution
  from ≤64 starts, 300 it each).
- PhaseModelAlone(k); **ReportAbsentPhaseIndexes() l.1079 — called from packDataBr() on EVERY call** (not trace
  gated): for absent non-ideal phases with Falp==-1 or all Gamma==1, up to kEvalBudget=600 activity evaluations
  (each restores ~11 L-sized arrays by copy), writes pm.Falp (→ DATABR omPH). Hidden per-call cost + a "report"
  that changes an output field. DECIDE phindex-ownx.
- thread_locals defined l.1221-1226: g_tpdStrongSearch, g_ipmLevelUsed, g_ipmNoProgArmed, g_ipmNoProgFired,
  g_ipmLevelOff.
- RefreshFromAmounts() l.1227 (X=Y; TotalPhasesAmounts ×2; CalculateActivityCoefficients; CalculateConcentrations;
  FX=TotalGibbsEnergy with X/XF/XFA restored by copy — TotalGibbsEnergy has side effects).
- PhaseHessianTwins() (sort stoichiometry rows per phase; O(FIs²)), WarmTwinSplitCheck() (pa_IpmTwinSplit; DECIDE
  twinsplit), NativeTpdPhase() l.1330 (timing Scope Tpd; same restore-by-copy evalY closure; ideal shortcut; alone
  → SolModActCoeff only).
- (ipm_main part 1b, l.1370-2367) NativeTpdPhase cont.: starts (corners, ideal, centroid, grid if
  g_tpdStrongSearch & n∈{3,4}, 1 % grid for binaries), pa_TpdStrongBasins selection, follow() successive
  substitution, golden-section refine for binaries; std::function recursion `rec`.
- native_trace_run_result() l.1521: KEY (FNV-1a XOR-folded assemblage hash; names "~2" for repeats), CERT
  (mb_rel/abs/chg, kkt_max, curv_min, dual_dirs, dual_resolved from the process-static, stab_ss ...), RANK.
  Cert work (CertCurvMin, CertStabTPD perturb & restore MULTI + TSolMod state) runs ONLY with the trace open →
  risk: trace on/off could change a later warm solve if a restore is inexact (documented as bit-identical).
- native_trace_mb_of/native_trace_mb, native_trace_stage (STAGE records), native_trace_mbr_exit (MBRX):
  comment l.1798-1805 documents the cause of "cold MBR returns iRet=0 on a state its own test rejects": guard
  `if( pa_p->DW && ( WhereCalledFrom == 0L || pm.pNP ) )` false for a cold call's 2nd MBR (K2>=1, pNP=0) →
  explains brief §8 (cold path answers failing own mass balance; only warm calls check).
- MassBalanceReproject(amt, keepPartial) l.1828-2050: residual C = B - A·amt (row loop over j: FP order fixed);
  sort species by amount; MGS rank-revealing pivots (1e-8); dense Gaussian elimination with partial pivoting
  (hand-written); clamp negatives; ≤8 passes; accept only if worst rel strictly falls; sync X,Y; TotalPhasesAmounts
  ×2; CalculateConcentrations; DECIDE mbreproject / -reverted / -partialkept / -partialdiscarded; warn text.
- ICIsNumericalTrace, ICIsDefaultSeed, SubFloorElementCheck (warn; toReal = 1/SizeFactor).
- EnergyDeterminacyCheck() l.2164 (pa_DeterminacyWarn; env GEMS3K_DETERMINACY_PROBE → stderr DETPROBE lines):
  independentSubset (MGS twice), K0 assembly, hand-written LU with partial pivoting, per-phase compliance.
- DUPLICATED DENSE LINEAR ALGEBRA: modified Gram-Schmidt rank tests appear in CertDualFreeDirs, CertRank (×2),
  MassBalanceReproject, EnergyDeterminacyCheck (independentSubset), ipm_optima.cpp free-dual search; hand-written
  LU/Gauss in MassBalanceReproject and EnergyDeterminacyCheck; Jacobi eig in CertSymEigMin (copy of
  SymEigFloorInPlace in ipm_optima.cpp); plus JAMA and Eigen. → one small dense-LA utility (Stage A: move with
  identical loop order; Stage B: replace by Eigen).

## Build status (session)
- Optima (optima_develop 145aeed) Release → /home/user/hub-install/optima (conda gcc 16.2, Eigen 3.4.0).
- GEMS3K Release -g -fno-omit-frame-pointer, USE_OPTIMA_SOLVER=ON, USE_THERMOFUN=ON (ThermoFun 0.7.3 conda),
  BUILD_SOLMOD=OFF → /home/user/hub-install/gems3k. 36 warnings at -Wall (18× -Wtemplate-id-cdtor, sign-compare,
  unused-but-set). tsolmod4rkt target FAILS outside /usr/local: root CMakeLists.txt l.101-102, l.167 hard-code
  include_directories("/usr/local/include"), ("/usr/local/include/eigen3") — SOLMOD_OBJECT does not get the
  spdlog include from the imported target (finding for build-system step).
- Configure rewrites tracked GEMS3K/gems3k_version.h (restored with git checkout after each configure).
- Profiling driver: scratchpad/prof/prof_driver.cpp (GEM_init + N× GEM_run in a mode; warm modes seed with one
  cold call then sweep bIC×(1+1e-3·sin(r+1))). Corpus: scratchpad/corpus (21 projects copied).

### ipm_main.cpp part 2 (l.2367-2926)
- EnergyDeterminacyCheck cont. (per-phase compliance solve; first_note_for_system; DECIDE undetermined).
- TraceGrowthCheck() l.2438 (pa_GAS; sets trace_growth_stuck/trace_growth_note → used by NativeOptimaFinish;
  DECIDE tracegrowth).
- ExcludeRedundantDCs() l.2515-2693: identical stoichiometry/class/G0,H0,S0,Cp0,Vol within 1e-12; (a) in-phase
  (skips IPx members, LO, multisite, sorption), (b) single-species phases (keep rule pure/name/first); hold sets
  DLL=DUL=0, RLC=BOTH_LIM, moves X/Y amount to the kept copy; **`static std::mutex warnedMutex; static
  std::set<std::string> warned`** (process-wide warn-once memo, grows unbounded) — same in StrandedElementCheck.
  RestoreRedundantDCs restores RLC/DLL/DUL. DECIDE redundant-dc(-suspicious/-constrained).
- StrandedElementCheck() l.2715 (env GEMS3K_STRANDED_PROBE → stderr).
- **GibbsEnergyMinimization()** l.2797: Reset_uDD; reset CondNum/CondNumDiag/SolveTimeMs/CondNumTimeMs/
  SolveCallCount; CALL trace line; `FORCED_AIA:` label; GEM_IPM_Init(); warm: IT=max(30,ITaia),
  Set_DC_limits(true); IAstatus = GEM_IPM_InitialApproximation(); INIT trace line; warm with an uncarried IC →
  pNP=0, MK=0, `goto FORCED_AIA` (DECIDE warmnocarrier); if IAstatus==false: GEM_IPM(-1); ITaia=IT (cold);
  U_r=U·RT; GasParcP(); IT=ITG; warm zero check (only if warmInputChanged_) …
- GibbsEnergyMinimization cont. (l.2926-3049): warm zero check (test amount 1e-20; share>1e-6 → MK=2; DECIDE
  warmzero); WarmTwinSplitCheck → MK=2; **MK==2: warm → pNP=0, MK=0, goto FORCED_AIA (cold restart); cold with
  nCNud<=0 → Error(errorCode); cold with u-divergence → ITG=0, goto FORCED_AIA** (retry inside the kernel — a 2nd
  retry ladder below the TNode layers); FitVar[0]=bfc_mass(); ANSWER trace line; TraceGrowthCheck (if !MK&&!PZ);
  final mass-balance check of pm.X (warn only) + MassBalanceReproject(pm.X) if fails (pa_MbReproject);
  EnergyDeterminacyCheck; StrandedElementCheck; testMulti() if MK||PZ.
- **GEM_IPM(rLoop)** l.3058-3476 (= the IPM driver; goto-based state machine with labels `mEFD:` and local
  `FORCED_AIA:` — same label name as in GibbsEnergyMinimization but a different target): W1=K2=Ec=MK=PZ=0;
  insBudgetTried.assign(FI,0); `mEFD:` STAGE MBR enter → **MassBalanceRefinement(pm.K2)** → STAGE exit; eRet 1..5:
  warm → MK=2 + goto FORCED_AIA(finish), cold → Error(); STAGE IPM → **InteriorPointsMethod(status)** → eRet 0:
  CalculateActivityCoefficients(LINK_PP_MODE); 1,2: warm→MK=2, cold→Error (E06IPM); 3,4: warm→MK=2 else Error;
  5 (u divergence): MK=2 always.
  PSSC gate `if( pa_p->PC >= 2 && nCNud <= 0 )` (trap 12: PC=1 → old PhaseSelect path l.3283; PC=0 → none):
  cleanupStatus = (PC>2 ? 0 : 1); **PhaseSelectionSpeciationCleanup(k_miss, k_unst, cleanupStatus)**; then
  MassBalanceReproject(pm.Y, true) if k_miss<0; rcode 1 → PZ=0; 0 → PZ=1, `goto mEFD` (new loop: MBR+IPM again,
  K2 counts loops); -1 → W08IPM, warm→MK=2, cold→PZ=2.
  PC==1: PhaseSelect + optional SpeciationCleanup (PRD≠0) + CalculateActivityCoefficients(LINK_PP_MODE).
  PC==3: XmaxSAT_IPM2(). Then (nCNud<=0) **second MassBalanceRefinement(pm.K2)**; OK → X=Y, TotalPhasesAmounts,
  CalculateConcentrations, PD>=2 → CalculateActivityCoefficients(LINK_UX_MODE); errors → warm MK=2 / cold Error.
  t_end=clock(); t_elap_sec (CPU clock, report only); `FORCED_AIA:` FI1/FI1s count (threshold
  min(PhMinM,1e-22)); **G[i]=G0[i] for all i** (G reset at every exit); t_elap_sec again.
- **GEM_IPM_InitialApproximation()** l.3486-3732: B[i] < pa->DB → raise to DB (W20IPM, PZ=3 then testMulti, PZ=0);
  **RescaleToSize(false)** (cutoffs/insertions scaled: SizeFactor); AllPhasesPure → pNP=0; COLD: ITaia=0,
  TotalPhasesAmounts, PCI=1, logCDvalues reset, lnGam=lnGmf (else 0), Gamma=1, G=G0+lnGam,
  CalculateConcentrations + IC=0 ("important for LPP AIA reproducibility"), zero Xeta/Xpsi/XcapD for sorption;
  IT=0, ITF+=1; **AutoInitialApproximation()** (simplex); pa_IpmLevelStart: if all mixed phases "checkable"
  (no aq/gas/fluid/sorption/ion-ex) → g_tpdStrongSearch=true; ColumnGenerationSeed(nCg, 1e-6) → Y=nCg,
  lvlLeaveOut (=IpmLevelStart>=2 && no mixed member), g_ipmLevelUsed=true (DECIDE natlvl);
  Y→X, MassBalanceResiduals, CalculateConcentrations; AllPhasesPure → return true (solved by LP);
  yLpFill snapshot if FilloutBudgetValue()>0; **DC_RaiseZeroedOff(lvlLeaveOut ? Ls : 0, L)** (class fill-out
  amounts); ApplyFilloutBudget; X=Y; PD==2 → FitVar[4]=-1 (no smoothing first time),
  CalculateActivityCoefficients(LINK_UX_MODE), FitVar[4]=1. WARM: warm_floor_seed (Optima zeroed floor amounts
  restored as start; DECIDE warmfloorseed), TotalPhasesAmounts, X=Y, PCI=1, logCD reset, PD>=2 → same
  FitVar[4]=-1/act/1 dance; pNP<=-1 → DC_RaiseZeroedOff(0,L). FitVar[4] is used as a FLAG (-1) and a value (1).
- **MassBalanceRefinement(WhereCalledFrom)** l.3755-…: TotalPhasesAmounts(Y); jK=MetastabilityLagrangeMultiplier()
  (≥0 → E17IPM ret 5); bestY tracking (ResidualMetric lambda over all N); loop IT1<pa->DP (ITF++ per iteration):
  NR=N (N-1 if aqueous gone: YF[0]<DSM && YFA[0]<XwMinM — assumes the aqueous phase is phase 0); 
  MassBalanceResiduals(N,L,A,Y,B,C); convergence test 3 variants (MbClassRule>0 / DT==0 relative / DT≠0 combined
  abs+rel; abs cutoff 10^-|DT|); best snapshot; converged → X=Y, TotalPhasesAmounts, PD==2 →
  CalculateActivityCoefficients(LINK_UX_MODE), return iRet; else **WeightMultipliers(true)**;
  **MakeAndSolveSystemOfLinearEquations(N, true)** (sRet==1 → iRet=1 E05IPM, continues); **DikinsCriterion(N,true)**;
  **LM=StepSizeEstimate(true)**; LM<1e-6 → iRet=3 E03IPM break; LM=min(LM,1); **Y[J] += LM*MU[J]** (FP order:
  plain loop); PSTALL stall detector (maxDeltaY, oscillation counter ≥10) → `iRet = 0; break;` (unconditional
  reset; the degenerate_cause test is commented out → `degenerate_cause` set-but-unused warning).
- MBR exit: best-residual restore (only on non-ideal exits); **strict guard `if( pa_p->DW && ( WhereCalledFrom
  == 0L || pm.pNP ) )` l.4067 → iRet=2 E04IPM**; else lenient: X=Y, TotalPhasesAmounts, PD==2 →
  CalculateActivityCoefficients(LINK_UX_MODE); return iRet (often 0) = the cold 2nd-MBR exemption (brief §8).
- **InteriorPointsMethod(status)** l.4104-4378 (the GEMS3K IPM main descent): stall-window constants; resets
  ipmKktRescues, s_ipmKktFallbacks, s_ipmQrSortPivot; lnGmo=lnGam; **FX=GX(0)**; YFA=XFA; loop IT1<pa->IIM
  (IT++, ITG++ per iteration): NR (aqueous gone → N-1, assumes aqueous = phase 0); **PrimalChemicalPotentials(F,
  Y, YF, YFA)**; Uc[J][0]=U if PCI<=10·DXM; **WeightMultipliers(false)**; **MakeAndSolveSystemOfLinearEquations(N,
  false)** (ret 1 → E07IPM return 1); u-divergence detector (PLLG): Increment_uDD, Check_uDD (W14IPM, return 5);
  **PCI=DikinsCriterion(N,false)**; **LM=StepSizeEstimate(false)**; **LM1=OptimizeStepSize(LM)** (golden section,
  0.382/0.618, ≥2 GX evaluations, stop on |FX1-FX2| <= 1e-13·max); **FX=GX(LM1)** (GX copies Y→X!); logCDvalues
  shift; aqueous elimination if XF[0]<DSM or X[LO]<=XwMinM (XF[0]=XFA[0]=0 — phase 0 assumption again);
  **PD>=2 → status=CalculateActivityCoefficients(LINK_UX_MODE)** EVERY iteration (warm & status → return 3);
  GEMS3K_IPM_PROBE IPMIT line (O(N·L) per iteration when set); **PCI <= DXM → CONVERGED** (PCI set to
  DXM·0.999999 "for smoothing"); pa_IpmStallWindow noise-stall accept (vectors with erase(begin()) per iteration,
  ring buffer (W+1)·L of X copies); pa_IpmNoProgress (g_ipmNoProgArmed) → g_ipmNoProgFired, break → falls to
  E06IPM (return 2) by design; u-divergence conditional convergence (PZ=5); **Restore_Y_YF_Vectors()**.
  Max iterations → E06IPM return 2.
- MassBalanceResiduals(N,L,A,Y,B,C) l.4388: **sparse loop over compressed A (arrL/arrAN): C[ii] -= A[jj*N+ii]*Y[jj]
  in species order** (FP order fixed by species order and by the compressed index; Build_compressed_xAN).
  CheckMassBalanceResiduals (cutoff min(DHBM·1e10, 1e-2), E02IPM).
- OptimizeStepSize(LM) l.4448: golden section with gotos SH1/SH2/SH3/OCT; returns midpoint (A+B)/2.
- DC_ZeroOff, DC_RaiseZeroedOff(jStart,jEnd,k) (class-dependent insertion amounts DFYaqM/DFYwM/DFYidM/DFYhM/DFYrM/
  DFYsM/DFYcM — the "fill-out"), FilloutBudgetValue() (reads env GEMS3K_FILLOUT_BUDGET on EVERY call, getenv in
  solver path), ApplyFilloutBudget.
- MetastabilityLagrangeMultiplier() l.4641 (pa_DKIN; clamps Y into [DLL+E, DUL-E]; returns j of a negative Y).
- **WeightMultipliers(square)** l.4674: W_j by RLC (NO/LOWER: Y-DLL; UPPER: DUL-Y; BOTH: min), squared (MBR) with
  ±1.34e120 clamp, or max(.,0) (IPM); bad code → E16IPM Error.
- `#define a(j,i) (*(pm.A+(i)+(j)*Na))` macro (A stored species-major: column j contiguous over ICs).
- GEMS3K_BENCHMARK_DIAGNOSTICS block (ENABLE_BENCHMARK_DIAGNOSTICS CMake option): DiagRatioConditionProxy,
  PowerIterationMaxEig, InverseIterationMinEig (steady_clock timing into pm.SolveTimeMs/CondNumTimeMs).
- **SolveIpmAugmented(N)** l.4855 (pa_IpmAugmentedKKT=1/2, default 2, main IPM loop only): `thread_local`
  scratch vectors act, dreg, Q, b, Rd, v (reuse → no per-iteration allocation after warm-up); active set
  Y > min(lowPosNum, DcMinM); dreg = Σ W a² per row; zero-row rescue 1e-12·dmax (DECIDE ipmkkt-zerorow, first only);
  mode 1: dense KKT (La+N)² with JAMA::LU (slow); **mode 2: hand-written Householder QR of the (La+N)×N stack
  [W^½ A_act; D^½]** column-major (FP order: row loops r=k..M), back-substitution into pm.U; rank test 1e-14·rmax;
  **pa_IpmQrSortPivot (default 100): La >= 100 → rows stable_sorted by size + Eigen ColPivHouseholderQR
  (thread_local Eigen::MatrixXd S etc., threshold 1e-14)** — Eigen version sensitivity lives here (brief §4).
  Returns 2 → caller falls back to normal equations (DECIDE ipmkkt-fallback, first only).
- **MakeAndSolveSystemOfLinearEquations(N, initAppr)** l.5097: augmented path for !initAppr; else Alloc_A_B(N)
  (reallocates only on size change?), zero AA; assemble AA += a(jj,ii)·a(jj,kk)·W[jj] over active species via
  compressed index (FP order: species-major accumulation); BB = C (MBR) or Σ F·a·W (IPM normal equations);
  MBR: Jacobi scaling Dscale (std::vector allocated EVERY call); TNT Array2D(N,N,AA) wraps AA (copy?) ...
- MakeAndSolve cont.: **JAMA::Cholesky chol(A); is_spd → solve; else JAMA::LU (isNonsingular else return 1)**;
  MBR: Uefd = B·Dscale (unscale); IPM: U = B. (Two LA libraries: hand QR / Eigen QR for IPM, JAMA for MBR and the
  fallback.) Array2D(N,N,AA) — TNT Array2D from a raw pointer (row pointers built per call: allocation per solve).
- **DikinsCriterion(N, initAppr)** l.5318: MU_j = W_j·(Σ a·U (− F_j for IPM)); PCI = Σ|qMu| (IPM) or 1/sqrt(Σ qMu²)
  (MBR); clamp 1.34e120. DC_DualChemicalPotential(U, A+J·N, N, J) per species (FP order: i loop).
- **StepSizeEstimate(initAppr)** l.5377 (first-limiting species; MBR ×0.95 or PCI; IPM 1/sqrt(PCI) cap 1e10).
- **Restore_Y_YF_Vectors()** l.5440 (phase with XF<=DSM or sorption with XFA<ScMin → Y=0, lnGam=0, YF=YFA=0 —
  this is where exact zeros of whole phases are written in the IPM loop; trap 12 context).
- **RescaleToSize()** l.5474: SizeFactor is a local constant 1. → pm.SizeFactor=1 and every *M cutoff = pa value
  (the "rescale to size" is a no-op; the real scaling is ScaleSystemToInternal in ipm_simplex.cpp). DEAD LOGIC.
- Alloc_A_B (realloc on size change), Free_A_B, **Build_compressed_xAN()** (|a|>1e-12 index; raw new[]),
  Free_compressed_xAN, Free_internal, Alloc_internal, setErrorMessage/addErrorMessage (C buffers errorCode[100],
  errorBuf[1024]), uDD: Alloc_uDD, Free_uDD, Reset_uDD, Increment_uDD (3-point moving mean via Uc), Check_uDD.
- PSSC, GX, TotalGibbsEnergy, SpeciationCleanup, PhaseInsertionCeiling, PhaseSelect, RaiseDC_Value, StabilityIndexes
  are in ipm_chemical.cpp; CalculateEquilibriumState, Scale/Rescale, AutoInitialApproximation, InitalizeGEM_IPM_Data,
  MultiConstInit, GEM_IPM_Init in ipm_simplex.cpp; testMulti, DC_LoadThermodynamicData,
  load_all_thermodynamic_from_grid in ms_multi_diff.cpp; ColumnGenerationSeed, AbsentPhaseWorstDrive in ipm_level.cpp.
  => file names do not say what is inside (ms_multi_diff holds thermo loading; ipm_simplex holds the solve entry).

### ipm_simplex.cpp (1069) — read in full
- AutoInitialApproximation() l.39: U=0; raw new[] DN, DU(L+N), B1, AA1(T), STR(T), NMB(L+1) per call (try/catch
  delete); EPS=pa->EPS, GZ=1/EPS; sparse pack of -G and A (|v|>1e-19, ROW-major walk `A+i*N+j`); SolveSimplex;
  Y=DU; TotalPhasesAmounts; **FX=GX(0)**; MassBalanceResiduals.
- Karpov-Chudnenko 1992 simplex: SPOS, START, NEW (new double[M+1] per call — allocation INSIDE the simplex
  iteration), WORK (same), FIN; SolveSimplex(M,N,T,GZ,EPS,...) dense (M+1)² tableau, IT=200 max; after 200
  iterations with EPS<=1e-6 it silently goes to FIN (no error) — only EPS>1e-6 raises E01IPM. gotos MK3/MK4/FINISH.
- **CalculateEquilibriumState(NumIterFIA, NumIterIPM, runKinetics)** l.509 (= native entry): timing Scope(Native);
  InitalizeGEM_IPM_Data(); trace_growth_stuck=false; phase_floor_parked.clear(); warm_floor_seed.swap(
  zeroed_floor_amount); ExcludeRedundantDCs() + RAII RedundantRestore; t_start=clock(); ITF=ITG=0;
  RunKineticsStep() (if runKinetics); **DG>1e-5 → ScFact=SystemTotalMolesIC(); ScaleSystemToInternal(ScFact)**;
  prev_scfact; warm_floor_seed*=ScFact; warmInputChanged_ (bulk/T/P/DUL>0 change, relative 1e-12);
  switch(pm.tMin): A_TV_/U_SV_ bodies commented out (dead), H_PS_/_S_PH_/_S_UV_ nothing, default
  **GibbsEnergyMinimization()**; catch(TError) → RescaleSystemFromInternal, ITF/ITG, rethrow; normal →
  RescaleSystemFromInternal; NumIterFIA=ITF, NumIterIPM=ITG; return t_elap_sec (CPU clock()).
- SystemTotalMolesIC(): TMols=Σ_{i<N-E} B; SMols=pa->DG; ScFact=DG/TMols.
- **ScaleSystemToInternal(ScFact)** l.652 (trap 3): SizeFactor=ScFact; MBX, FitVar[0], VXc, HXc, HX_, FX, Yw ×ScFact;
  per j: DUL (if <1e6) with pre/post record, DLL (always) with record, Y, X, XY, XU; per i: B, BFC (**C not scaled**,
  commented); per k<FI: XFs, XF, YF, FVOL, FWGT; per k<FIs: XFA, YFA, PUL (<1e6, record), PLL (no record), BF;
  sorption XetaA/B/D, XFTS. NOT scaled: U, MU, W, Falp, lnGam, IC/pH (intensive, fine), kinetic outputs PrT/PkT,
  emRd, control-condition titrants, **zeroed_floor_amount/warm_floor_seed (scaled manually by the caller)**,
  any std::vector cache in TMultiBase (prevBulkMol_ stores B/ScFact explicitly).
- **RescaleSystemFromInternal(ScFact)** l.747: inverse; DUL/DLL/PUL exact restore when untouched since scale-in
  (record match), else divide; **C /= ScFact (asymmetric: C was never multiplied; correct only because C is
  recomputed in internal units during the solve)**; clears records; SizeFactor=1.
- **InitalizeGEM_IPM_Data()** l.847: MultiConstInit(); initalizeGEM_IPM_Data_GUI() (virtual hook); pTPD<2 →
  DC_LoadThermodynamicData(); **Alloc_internal() → Build_compressed_xAN() reallocates arrL/arrAN on EVERY call**
  (allocation per solve); Alloc_uDD; MBX=Σ B·Awt/1000; RescaleToSize(true) (no-op scale 1); warm: loadData(false)
  (virtual), X=Y, TotalPhasesAmounts, CalculateConcentrations, lnGmo=lnGam, Gamma=PhaseSpecificGamma.
- MultiConstInit() l.913: FI1.., ITF=ITG=0, PD=pa->PD, Ec=K2=MK=0, W1, is, js, next, ln5551, lowPosNum=1.66e-24,
  logXw=-16, logYFk=-9, DXM=pa->DK, **FX=kTotalGibbsEnergyUnset**, pH trash check, YMET, PCI=1, FitVar[4]=1,
  multiConstInit_PN() (virtual).
- GEM_IPM_Init() l.944: Uefd=0; cold: zero X,Y,lnGam,lnGmo,MU,XU,EMU,NMU,W,F,F0; Gamma=1; mixed phases:
  GEM_IPM_Init_gui1() (virtual), warm & E & LO: CalculateConcentrations, IS_EtaCalc, GouyChapman for sorption;
  **CalculateActivityCoefficients(LINK_TP_MODE)**; else PD=0, pNP=0, pIPN=1; GEM_IPM_Init_gui2() (virtual).
- pTSolMod(xPH), HelmholtzEnergy(x) / InternalEnergy(TC,P) (unused stubs for other potentials; A_P/U_TP commented).

### ipm_chemical.cpp (2166, **CRLF**) — part 1 (l.1-1119)
- XmaxSAT_IPM2() (PC==3; writes DUL of surface species from site capacities; SATX/MASDJ/Nfsp).
- **DC_DualChemicalPotential(U, AL, N, j)** l.204: Σ over the compressed index (arrL/arrAN) of U[ii]*AL[ii]
  (ii<N). NOTE: the compressed index drops |a| <= 1e-12 (Build_compressed_xAN) — sparse paths ignore tiny
  stoichiometry entries that dense paths (Cert*, TraceGrowthCheck, WarmTwinSplit, EnergyDeterminacy, Optima
  adapter) include → two semantics for "A" in one solver. FP order: ascending IC index.
- Set_DC_limits(InitState) (PLIM; converts PLL/PUL and DLL/DUL from units RFSC/RSC (mol, molal, mole/mass/volume
  fraction) into moles; E11IPM/E12IPM). Called warm in GibbsEnergyMinimization → writes DUL/DLL inside the scaled
  region (why ScaleSystem records pre/post values).
- TotalPhasesAmounts(X, XF, XFA) l.385 (sum per phase; carrier DC_ASYM_CARRIER → XFA) — FP order species-wise.
- **DC_PrimalChemicalPotentialUpdate(j,k)** l.417: F0 from lnGmM (+ sorption Coulomb terms psiA/psiB·Ez·FRT, TLM/
  BSM/CCM variants, surface density terms log(FactSurT), FactSur/(1+FactSur)) + lnGam; **smoothing: if sMod
  SGM_MODE != SM_IDEAL and !optima_disable_smoothing and |dF0|>=1e-5 → F0 = Fold + dF0·SmoothingFactor()**
  (history-dependent: Fold = pm.F0[j] from the previous evaluation) — the hidden-coupling flag
  optima_disable_smoothing is read here and in SmoothingFactor().
- DC_PrimalChemicalPotential(G, logY, logYF, asTail, logYw, DCCW): F by DCCW class (single, asym species,
  asym carrier with 1 - 1/(1-asTail), symmetric); 7777777 sentinel.
- **PrimalChemicalPotentials(F, Y, YF, YFA)** l.635: per phase skip tests (L1==1 & YF<PhMinM; YF<=DSM; aqueous
  with Y[LO]<=XwMinM); YF>=1e6 → E13IPM Error; **writes MULTI scalars as scratch: pm.YFk, pm.logXw, pm.logYFk,
  pm.Yw, pm.aqsTail** (NonLogTerm = 1 - YFk/Yf unless Pitzer); F[j] = DC_PrimalChemicalPotential(pm.G[j], log Y …).
- DC_GibbsEnergyContribution (unused: inlined in GX).
- **GX(LM)** l.747: LM<lowPosNum → X=Y; else **X = Y + LM·MU and `if( pm.X[i] < pm.DcMinM ) pm.X[i]=0.;`
  l.765-766 — the exact-zero truncation (trap 12)**; TotalPhasesAmounts; FX = Σ_k Σ_j x·(G + log x - logXw/logYFk)
  in phase/species order (FP order point); writes pm.logXw/pm.logYFk scratch. Called by InteriorPointsMethod (FX at
  entry, golden-section probes, final step), AutoInitialApproximation.
- **TotalGibbsEnergy()** l.838: mirror of GX(0) with Gj=G0+fDQF+F0 (path-independent G); copies Y→X, TotalPhasesAmounts
  (side effects; called by TNode::Get_GibbsEnergy(), so a *getter* mutates X/XF/XFA). DUPLICATE of GX/CertPrimal.
- ConvertGj_toUniformStandardState(g0, j, k): g0/RT + ln5551 (aqueous/sorption), + log P (gas/fluid), + YOF if
  testTSyst() (virtual hook).
- ConvertDCC() (DCC → DCCW; remaps DC_SUR_GROUP/COMPLEX → DC_SSC_A0 and IPAIR → DC_WSC_A0 IN PLACE in pm.DCC:
  input data mutated); getXvolume(); KarpovCriterionDC(...).

### ipm_chemical.cpp part 2 (l.1119-2166)
- KarpovsPhaseStabilityCriteria() (PC==1 path; EMU/NMU; uses GamFs/Y_la for removed phases; Falp = Σ − 1;
  cold absent mixture → -1).
- SpeciationCleanup(threshold, cutoff) (PC==1 + PRD path; MjuPrimal F vs MjuDual Fx; Yj/exp(ΔMju)).
- **PhaseInsertionCeiling(j)** l.1363 = min_i b_i/a(j,i) over ordinary ICs (trap 12 budget function).
- **PhaseSelectionSpeciationCleanup(kfr, kur, CleanupStatus)** l.1390-1849 (PSSC): MjuDiffCutoff = GAS (if >1e-6)
  else 1e-3; AmountThreshold = 10^-max(|PRD|,4); K2==0 → GamFs = Gamma (copy of γ before the first PSSC);
  XY = Y; PrimalChemicalPotentials; **StabilityIndexes()**; K2++; PSSC trace lines (enter / per phase / exit);
  pass 1 over phases: skip sorption & any phase with a species having DUL<1e6 or DLL>0 (KinConstrPh);
  stable (−DFM<logSI<DF) & amount<DSM → zero (ZERO_NEGLIGIBLE); logSI>=DF & amount<DSM → INSERT: single-species
  phase insertion amount = **min(pm.DFYsM, PhaseInsertionCeiling)**, budget-sized only once per solve
  (insBudgetTried[k]); mixed phase → DC_RaiseZeroedOff(jb, jb+L1, k) + Y[jb]=insAmt (!: overwrites the first
  end-member's raised amount); kfr=k, MassBalanceViolation; logSI<=−DFM & present → ELIMINATE (Y=0) — exact
  zeros #2; pass 2 (CleanupStatus && PRD): kinetic clamps, then cleanup Y/exp(ΔMju) with overflow guard ±608/609,
  pure-phase zeroing below AmountThreshold unless composition-limited (DECIDE cleanupzero skipped), over-stable
  pure phase re-insertion at DFYsM; Y<DcMinM → 0 (exact zeros #3). Exit: status 1/0/-1 (K2<6 → loop again, else
  -1 and Y restored from XY).
- **StabilityIndexes()** l.1854: lnPc, Xw from phase 0 if aqueous; per species ln a_dual = lg_to_ln·Y_la − class
  corrections (lnFmol−lnXw, lnPc, fDQF) clamped [−608,609]; x_est = exp(.)/γ (γ from Gamma or saved GamFs);
  StabIndex = Σ x_est; NMU, EMU written; Falp = log10(Σ) (|logSI|<log10(DSM) → 0); cold absent mixture without γ →
  Falp = −1 (the "placeholder" that ReportAbsentPhaseIndexes later replaces). FP order: species order.
- PhaseSelect(kfr, kur, RaiseStatus) l.1977 (PC==1 Selekt2; DC_ZeroOff/DC_RaiseZeroedOff loops; K2>4 → -1),
  RaiseDC_Value(j) (class insertion amounts).
- Duplicated per-class switch tables: DCC → DCCW (ConvertDCC), DCC → insertion amount (DC_RaiseZeroedOff,
  RaiseDC_Value: same table twice), DCC → concentration scale (StabilityIndexes, Get_cDC, DC_c,
  ConvertGj_toUniformStandardState, CalculateConcentrationsInPhase?) = 6+ switch tables over the same char codes.

## PROFILE OBSERVATIONS (survey, running)
- lb_test (LimBrine) SIA sweep (bIC×(1+1e-3 sin r)): r=0 OK_GEM_AIA (warm→cold fallback inside kernel, 1087 it),
  r=3/r=4 ERR_GEM_SIA after the full rescue ladder: 14.6 s / 8.3 s per call; GEMS3K_TIMING: native 0.024 s (12
  attempts), optima 14.60 s (5 calls), osolve 14.60, obj 10.98 s (134 649 objective calls), objfxx 10.09, **fdmix
  9.15 s (1.54 M finite-difference Hessian columns of mixed-phase species), act 9.30 s (1.81 M
  CalculateActivityCoefficients calls)**, objrestore 0.70 s. → in escalated calls ~85-90 % of time is activity-model
  evaluation for FD Hessian columns. Normal native AIA call: 2.7 ms.

### ipm_chemical2.cpp (1443, LF) — read in full
- `thread_local bool g_aqueousAbsent` defined l.34 (set by CalculateConcentrations; read by GEM_run_tries and
  native_trace_run_result).
- phase_bcs (sparse BF = Σ A·X with |X|<1e-16 skip), phase_bfc (solids → BFC, |X|<1e-12 skip; different
  thresholds), bfc_mass.
- CalculateConcentrationsInPhase(X,XF,XFA,Factor,MMC,Dsur,jb,je,k) l.109: per species Fx = dual potential
  (DC_DualChemicalPotential with NR), Wx, VL, Y_la (dual log activity by class), Y_m (molality), Y_w, **pH from the H+
  dual potential**, IC (ionic strength) increments, FVOL, IC_m/IC_wm. Assumes aqueous = phase 0 (`pm.PHC[0] ==
  PH_AQUEL`, `!k`).
- **CalculateConcentrations(X, XF, XFA)** l.306: resets BFC, Wx, VL, VXc; per phase FWGT/FVOL; single-component
  phases; aqueous (LO && !k) resets IC, IC_m/lm/wm; absent phases → pH=pe=Eh=0 and g_aqueousAbsent=true;
  `pm.Wx[j] = 0.0` inside a loop over jj (l.387: uses j not jj — writes the same element repeatedly; latent bug);
  phase_bcs → BF; aqueous: MMC, Factor, **pe = ln_to_lg·U[N-1], Eh = 0.000086·U[N-1]·T** (hard-coded constant;
  charge assumed last IC); sorption: **XFA[k] clamped to ScMinM in place (mutates its input)**, XFTS, logYFk;
  CalculateConcentrationsInPhase; VXc += FVOL; IC_lm, IC_wm normalisation.
- IS_EtaCalc() (surface charge per plane XetaA/B/D), GouyChapman(jb,je,k) (diffuse-layer potentials, many
  truncations ±0.4..1.1 V with status codes 60-76), SurfaceActivityCoeff(jb,je,..,k) (SACT isotherms: Langmuir,
  QCA, Frumkin, Pivovarov, virial, BET stub, obsolete SAT_COMP/NCOMP, SAT_SOLV; smoothing with lnGmo ±0.2).
  => the classic sorption ACTIVITY MODELS live inside TMultiBase (solver class), not in TSolMod/TSorpMod:
  plug-in boundary violated; s_sorpmod.* is the unfinished replacement.

### ipm_chemical3.cpp (1212) — read in full
- SmoothingFactor() (optima_disable_smoothing → 1; FitVar[4]<0 → 1; FitVar[3] or FitVar[4]); SetSmoothingFactor(mode)
  (AG, DGC, IIM; three variants; writes FitVar[3]) — smoothing state lives in MULTI.FitVar (shared scratch).
- PhaseSpecificGamma(j, jb, je, k, DirFlag) (lnGam ↔ Gamma with asymmetry NonLogTerm; Pitzer excluded).
- **CalculateActivityCoefficients(LinkMode)** l.291 (timing Scope ActCoef):
  LINK_TP_MODE: FitVar[3]=1; per mixed phase (and single-species gas): index offsets ipb/jpb/jdb/jmb/jsb walked
  cumulatively (offsets recomputed at every call by prefix sums over LsMod/LsMdc/LsMdc2); Wx floor, fDQF=0,
  lnGmo=lnGam; SolModCreate + SolModParPT(k, code).
  LINK_PP_MODE: SolModExcessProp/IdealProp/StandProp/DarkenProp → GPh..UPh.
  **LINK_UX_MODE: SetSmoothingFactor; CalculateConcentrations(X, XF, XFA) (full pass over ALL phases and species,
  incl. dual potentials); lnGam=0, Gamma=1 for all L; IS_EtaCalc (+ zero psi)**; then per phase k: pmpXFk = Σ X;
  **XF[k] overwritten by pmpXFk if XF[k] < DSM (mutation)**; dispatch on PHC + sMod char code (SGM_MODE ==
  SM_STNGAM): aqueous (IC > ICmin) → SolModActCoeff; gas/fluid → SolModActCoeff (6 ifs); condensed → SolModActCoeff;
  sorption → GouyChapman + SurfaceActivityCoeff; calculateActivityCoefficients_scripts() (virtual, GEMSGUI);
  END_LOOP (UX): lnGmo = lnGam + lnGamG; Gamma = PhaseSpecificGamma; **F0 = DC_PrimalChemicalPotentialUpdate(j,k)
  (smoothing)**; **G = G0 + fDQF + F0**.
  => ONE CALL EVALUATES EVERY PHASE. The Optima FD Hessian (fdmix) perturbs one species per column and pays a full
  CalculateActivityCoefficients each time (profile: 1.8 M calls). Per-phase evaluation is only exact if each
  TSolMod is a pure function of (its phase composition, T, P, IC) — must be certified per model (some keep state:
  PTparam caches, fluid EoS density start values?, lnGmo-based smoothing in SACT).
- **SolModCreate(...)** l.597: fills SolutionData with POINTERS into MULTI arrays (IPx, PMc, DMc, Wx, lnGam, lnDQFt,
  lnRcpt, lnExet, lnCnft, CTerms, FVOL+k, XF+k, DCC, MoiSN, SitFr, fDQF, Pparc, dcMod, DQFc, Y_m, EZ, Vol, SM,
  denW, epsW, FWGT+k, X+jb, G0+jb); reuses phSolMod[k] if testSizes() (then only UpdatePT); else `new` model by
  char code switch (≈27 cases; SM_AQELVIS commented out → phSolMod[k]=nullptr but SM_AQELVIS IS listed in
  SolModActCoeff/ParPT → ErrorIf "Invalid index of phase" if a project uses ELVIS: DEAD/BROKEN path);
  debug level → writes solmod_<k>.json files into the CWD.
- SolModParPT / SolModActCoeff / SolModExcessProp / SolModIdealProp / SolModStandProp: each repeats the same list of
  ≈27 model codes in a switch → **adding one activity model = edit ≥8 places** (enum m_const_base.h, s_solmod.h class,
  s_solmodN.cpp, SolModCreate, 5 wrapper switches, CalculateActivityCoefficients dispatch, PhaseModelAlone,
  ExcludeRedundantDCs? , optima_smod_is_multisite, maybe GEMSGUI). SolModDarkenProp writes zeros (stub).
  SolModActCoeff at debug level appends solmod_act_coef_<k>.txt per call (file I/O in hot path when debug).
- Alloc_TSolMod/Free_TSolMod, Alloc_TSorpMod/Free_TSorpMod (raw new[] of pointers).

### ipm_chemical4.cpp (859) — read in full
- RunKineticsStep() l.57: degenerate-box report (snapshot DUL/DLL, DECIDE degenbox, warn-once via kinDegenWarned_);
  **comment l.80-82: "a cold call has not [copied the amounts], so pm.X still holds the previous solve's answer,
  which on a shared solver can be another node's"** → MULTI state is shared between the nodes of a TNodeArray that
  use one calcNode; history (X, F0 smoothing state, lnGmo, TSolMod internal state, prevBulkMol_, warm_floor_seed,
  zeroed_floor_amount, phase_floor_parked) crosses node boundaries. Cold call: X:=Y, CalculateConcentrations
  (DECIDE kincold), restore after; pH from H+ molality·γ (DECIDE kinph); pKMM<2 → CalculateKinMet(LINK_TP_MODE)
  + (LINK_IN_MODE) or (LINK_PP_MODE).
- CalculateKinMet(LinkMode) (prefix-sum offsets over LsKin/LsUpt/LsPhl; only PH_LIQUID/SINCOND/SINDIS phases).
- KM_Create (KinMetData with ~25 pointers into MULTI incl. **DUL+jb/DLL+jb as direct outputs**; factory switch on
  kmod[0]; TUptakeKin gets IC_m, emRd, emDf); KM_ParPT (dynamic_cast), KM_InitTime/UpdateTime/UpdateFSA/ReturnFSA/
  InitRates/CalcRates/InitSplit/CalcSplit/SetAMRs/CalcUptake/InitUptake — 11 wrappers each with the same switch and
  **C-style downcasts `(TMWReaKin*)phKinMet[k]`** (virtual calls would do; a wrong code is UB). ReturnFSA writes
  pm.Aalp[k], PfFact, PrT, PkT, PvT, PUL/PLL. Alloc/Free_TKinMet.

### ipm_level.cpp (386) — read in full
- thread_locals g_cgConverged, g_cgPseudo, g_cgGap, g_cgDual (OUTPUTS of ColumnGenerationSeed passed to callers
  through thread-local globals); extern g_tpdStrongSearch.
- **TwoPhaseSimplexMinSum(nRows, nCols, std::function a(i,j), b, xOut, cost, yOut)** — free function: dense
  two-phase simplex with Bland's rule; element access through std::function in every inner loop (slow); verifies
  dual sign. = THIRD LP solver in the code (Karpov SolveSimplex in ipm_simplex.cpp; this one; LPFeasibilitySeed /
  LPGibbsDual in ipm_optima.cpp probably use this one — check).
- ColumnGenerationSeed(nOut, tol) (timing Level): species LP priced G0+fDQF → dual → NativeTpdPhase per condensed
  mixed phase → pseudo-compound columns, ≤20 rounds; saves/restores pm.U; spreads pseudo amounts over end-members;
  DECIDE cgseed / tpdscreen; mb residual check.
- AbsentPhaseWorstDrive(worstPhase) (pure: reduced cost; mixed: NativeTpdPhase with g_tpdStrongSearch).

### ipm_optima.cpp part 1 (l.1-1639) — Optima adapter, helpers
- Includes Optima/Optima.hpp, Sensitivity.hpp, Trace.hpp. `timed_solve` template wraps solver.solve in Scope OSolve.
- 29 thread_local mentions. extern from node.cpp: g_optimaCgSeedArmed, g_optimaTargetWarned, g_optimaIterBudget(-1 none),
  g_optimaCapAllowance, g_optimaDryFailed/Off/RedoCall, g_optimaPreSolveEmpty, g_lvlFullCalls, g_optimaFixAccept,
  g_optimaLineSearchRetry, g_optimaLSEscapeOff, g_optimaEscRetryNext; from ipm_main: g_ipmNoProgFired, g_tpdStrongSearch;
  from ipm_level: g_cgConverged, g_cgPseudo, g_cgGap, g_cgDual(vector); from ipm_chemical2: g_aqueousAbsent.
  static thread_local here: g_ps0GaveUp, g_coldLvlUsed, g_lvlShortcut, g_lvlCgStart; defines g_finishFromSuccess (l.102).
  Function-static thread_local: warnedPairs set (l.759), traceXPrev/tracePhPrev/traceCumIt (l.2608-2610).
- memoize_objective_if_linesearch (l.120): comment states GEMS3K objective is NOT pure (CalculateActivityCoefficients
  LINK_UX_MODE accumulates lnGmo, blends F0 via FitVar[3]); memoises f at bit-identical x so line-search re-evaluation is
  stable. => evidence for trap "activity not idempotent".
- optima_trace_env(): GEMS3K_OPTIMA_STEPTRACE_<name> with REAKTORO_TRACE_<name> fallback. optima_trace_hessian_dump_now:
  function-static env caches; optima_trace_write_matrix CSV.
- apply_optima_lu / _snap / _linesearch: gated by fork macros OPTIMA_LU_PARTIAL_PIVOT, OPTIMA_BACKTRACK_SNAP_LEAVING,
  OPTIMA_LINESEARCH_STALL_ESCAPE, OPTIMA_LINESEARCH_REJECT_WORSE (silently no-op on upstream Optima => build-config trap).
- SymEigFloorInPlace (l.244, Jacobi eigen floor) duplicates CertSymEigMin (ipm_main). kLnFmol=log(1000/18.01528) hard-coded.
  optima_default_box_moles(pm, DG) (l.321).
- SetControlCondition_pH/_Eh (l.339/378): look up "H","Zz","H+" by name via node1 (TNode back-pointer); lambdas capture
  `this`; Eh uses 0.000086*T (approx k/e, not F/RT constant from m_const_base). ClearControlConditions l.418.
- OptimaMaxMassBalanceResidual l.423, HasAqueousPhase l.440, CheckRedoxFree l.451 (recursive std::function, 200000 cap),
  OptimaTwinLabel l.575 (swaps 12 DC arrays + 9 phase arrays then RefreshFromAmounts), PhaseHessianBegin/End l.664/672,
  WarnCollapsedTwins l.710 (static thread_local warnedPairs).
- ic_index_by_name l.793 "works without a node (GEMSGUI solver has none)". DetectSolventCollapseAndReseed l.805 (names "H","O").
  DetectPhaseCollapseAndReseed l.845 gates on LPGibbsDual reduced costs (DECIDE reseedR).
- PotentialSpaceFinish l.977-1497: Newton on fixed phase set; FD Hessian via evalF (one CalculateActivityCoefficients per
  column); hand-written dense LU factorDense/solveFactored + 2 refinement steps (another dense-LA duplicate); LM damping;
  exact-penalty merit; Armijo 30 halvings; <=6 tabu set changes; 60 it; n<=900; TPD entry of absent phase;
  pa_OptimaFinishStall; DECIDE finish; restores state on failure.
- LPFeasibilitySeed l.1499, LPGibbsDual l.1539 -> TwoPhaseSimplexMinSum (third LP solver).
- WorstPhaseStabilityViolation l.1559: TotalPhasesAmounts(Y), StabilityIndexes, kinetic exemption, clamped detection by
  NMU, structural floor rule, `static const kRelPresenceEps 1e-12`.
- Function index: OptimaReducedPreSolve l.1745, CalculateEquilibriumStateOptima l.2464 (to ~6870, ~4400 lines one
  function!), CalculateEquilibriumStateHOP l.6874, NativeOptimaFinish l.7085, NativeOptimaRescue l.7157.

## Survey complete (21 projects incl. ls_test twice: two -dat.lst; 5 calls/mode; GEMS3K_TIMING files in prof/survey)
- All exit 0. Failures only: lb_test SIA ok=1/5 (23 s), SOP ok=3/5 (175 s), SHP ok=3/5 (161 s); mines23 SIA ok=3/5.
- Slow cold Optima: CemGEMS-keyvalue / Export-condensed AOP ~555-608 ms/call (it=494/call) vs AIA 21-22 ms;
  mines23 AOP 1848 ms vs AIA 591 ms. CASHNK1 AOP 62 ms vs AIA 1.8 ms.
- Typical small systems: 0.1-3 ms/call any mode. Warm < cold in IPM except process/Kaolinite (SIA 949 it vs AIA 775).

### ipm_optima.cpp part 2 (l.1639-4100)
- Wpsv end: GEMS3K_PHSTAB_PROBE found at l.1663 and l.1700 (function-static FILE* opened once via getenv lambda,
  `static FILE* dbg = [](){...}()`, two separate statics for the same env var, two fopen of same path "a").
- OptimaReducedPreSolve l.1745-2462: active-set pre-solve; LPGibbsDual pricing (tol>0 threshold, <0 rank rule);
  Problem(Aex,be,xlower,xupper) copied from pm.A/pm.B per pass; objective lambda is a near copy of the full-path
  objective (ideal-mixing Hessian, FD columns for ibasicvars, PhaseHessianFloor block) -> duplication (~200 lines);
  sign: state.ye = -pm.U (Optima dual has opposite sign); dead code: `const int mode = 3; const bool capOn = true;`
  makes mode==2 branch unreachable (l.2231-2245). getenv("GEMS3K_OPTIMA_TRACE_FILE") per pass,
  GEMS3K_PRESOLVE_RESID_PROBE. Writes pm.X,pm.Y,pm.U; discard() restores Y,U only (X, XF, F, lnGam not restored ->
  caller re-establishes with reestablish lambda).
- TRAP 10 evidence: PreStallWatch.maxSeconds = pa_OptimaMaxSeconds, steady_clock (l.1868-1885); StallWatch same for
  full solve (l.3695-3736), checked in convergence hook l.4061 `CalculationCancelled() || stallWatch->overBudget()`.
  Comment l.5550: "Gated on stalled, not timedOut". The other clock uses are reporting only (pm.t_start/t_end
  clock(), ipm_main l.3445/3470/3649, ipm_simplex 528/615/627) and GEMS3K_BENCHMARK_DIAGNOSTICS timers
  (ipm_main l.5100-5306). CalculationCancelled() (v_detail.cpp l.25, static std::atomic<bool>) is a second
  non-deterministic stop, caller-driven (also node.cpp l.567).
- More mutable statics: s_solmod4.cpp l.4678 `static double s;` (process-wide, in an activity model!); ipm_main
  l.2672/2775 static mutex+set warned (warn-once). s_ipmKktFallbacks is static thread_local (ipm_main l.94),
  reset l.4125.
- CalculateEquilibriumStateOptima l.2464: RAII guards SmoothingGuard (optima_disable_smoothing), FloorRecordGuard
  (PhaseHessianEnd), LuRecord/SnapRecord (fork macros), RedundantRestore, RedoxTargetRemove, PrevScFactGuard.
  Order: InitalizeGEM_IPM_Data -> ExcludeRedundantDCs -> clock -> RunKineticsStep (if runKinetics; false on HOP
  leg) -> CheckRedoxFree (pa_OptimaRedoxFree>=2 adds Eh target pe 4) -> ScaleSystemToInternal if DG>1e-5 ->
  SubFloorElementCheck -> GEM_IPM_Init_gui1 (GUI hook) -> CalculateActivityCoefficients(LINK_TP_MODE) ->
  warm usable test (any pm.U != 0 and finite) -> cold seed: ColumnGenerationSeed (pa_OptimaCgSeed &&
  g_optimaCgSeedArmed) / leveled start (pa_OptimaLevelStart; g_tpdStrongSearch toggled) / LPFeasibilitySeed /
  cold-level (pa_OptimaColdLevel, twin copies start apart) / DetectSolventCollapseAndReseed /
  DetectPhaseCollapseAndReseed / uniform seed -> X=Y, activity UX -> dim-reduce pre-solve (two attempts, tol 10
  then rank -3; g_optimaPreSolveEmpty) -> GEMS3K_PRESOLVE_HANDOVER_PROBE dump -> Problem build.
- MULTI written by the Optima entry (first part): Y, X, XF, XFA, YF, YFA, U (dual), F, lnGam (via activity), t_start,
  t_end, t_elap_sec, ITF, ITG, Ec, MK, PZ, errorMessage; members prev_scfact, phase_floor_parked,
  zeroed_floor_amount, optima_disable_smoothing, optima_control_conditions, phfTwin/phfFired.
  Reads: N, L, Ls, FI, FIs, L1, PHC, A, B, G0, DLL, DUL, LO, SB, SM, SF, ICC, pNP, stkey, sMod, T, XwMinM.
- dn/db SENSITIVITY (for 07): optima_want_sensitivity -> dims.c = N; problem.c = pm.B; problem.bec = identity;
  "Sensitivity::xc is dn/db"; with dims.c = 0 Optima skips the sensitivity solve. (read-out later in function)
- Control conditions (pH/Eh titrants): extra slots L..L+R, fixedGradientFn(target), Aex column = stoich, bounds
  +-2*default box.
- Full objective lambda l.3204-3520: captures `this`, &fixedGrad (by reference! problem must not outlive scope),
  shared_ptr state (phLast, phDec, phMax, nzSeen, fdSuppress, traceLastF). Per call: X=x, TotalPhasesAmounts,
  CalculateActivityCoefficients(LINK_UX_MODE), PrimalChemicalPotentials; f = sum X F - tau sum ln X (pure);
  Hessian: aqueous molality block, mole-fraction (pa_OptimaMoleFracHessian), FloorCurvature 1/2/3, log-barrier,
  FD columns for opts.ibasicvars (one activity call each; GEMS3K_TIMING FdMix/FdPure), PhaseHessianFloor block
  (FD + symmetrise + SymEigFloorInPlace); restore (activity again). Trace HESS/HESSBASIS lines and CSV dump
  (std::ofstream) when GEMS3K_OPTIMA_STEPTRACE_HESSDUMP style env set. ipm_logger->info inside objective (nzSeen).
  => profile: FD Hessian = N activity calls per Hessian = dominant cost in escalated calls (lb_test).
- State seeding: x = max(Y, xlower); warm floor snap via prevScFactIn/ScFact; dual ye = -U on warm or after
  pre-solve; cold: g_cgDual or LPGibbsDual if water share >= 0.2 (names "H").
- Options: maxiters = max(2000, IIM), tolerance OptimaTol; referenceMode (ROP) keeps library defaults.
  pa_OptimaLSWarmTrials, apply_optima_snap (warm, not HOP leg), apply_min_max_fix_and_accept for quick try
  (g_optimaFixAccept), maxiters=1 when g_lvlShortcut.
- StallWatch (l.3687-3735) + convergence.check hook (l.3880-4100+): quick-try watchers use thread_local
  g_optimaQuickTry (OptimaQuickTryWatch: cstep/cf0, step/f0/preMin/postMin/jumpAt, errPrev/frozenRun,
  rises/fPrev/riseSteps); hook mutates g_optimaIterBudget (global budget shared across retries in GEM_run);
  early-stability trend/cap (optima_earlystability_at resolver); GEMS3K_EARLYTREND_PROBE (+_N) at l.3814-3826
  (function-static FILE*); dry-retry stop (water < 1e-3 cap for 100 it); pa_OptimaFrozenCut; tail stop.
  Hidden-state: the hook's results flow back via shared_ptr flags and thread_locals, not via Result.

### ipm_optima.cpp part 3 (l.4100-7195), read by sections + DECIDE index
- The retry ladder INSIDE CalculateEquilibriumStateOptima (one function, l.2464-6870, ~4400 lines), in order:
  tail stop (tailstop) / stall (optimastall) l.4177-4244; tiered Hessian cheap attempt (pa_OptimaFDHessianDelay)
  l.4339; start-state phase-stability scan (presolvephstab) l.4399; PRIMARY SOLVE l.4458-4466
  (`timed_solve(solver, problem, state[, sensitivity])`); early-probe net (earlystop) l.4479; dry-out (dryout,
  tailpin, tailaccept, dryback, drycheck, drywet) l.4568-4806; water skip/keep l.4969-4991; dry retry l.5069;
  phase-extinction tier (extinction, extallow: refills g_optimaIterBudget from g_optimaCapAllowance l.5208);
  phase-selection loop (phasesel) l.5314; twin readmit l.5435; net restore (netunpin) l.5468; finish-first
  l.5533; stall net (warmfailfast, stallnet) l.5545; early-stability net (netresume/netresolve/netfallback,
  mode from optima_net_resume_mode(), env GEMS3K_OPTIMA_NET_RESUME overrides!) l.5596-5710; trace respread
  l.5814; CalculateConcentrations l.5855; collapse check l.5891; KKT check l.5858-5960; tail KKT; dual-free
  (dualfree) l.6087; titrant; OptimaZeroAbsent (optimarepair, zeroabsent) l.6214-6374; TPD accept (tpdaccept*,
  tpdveto) l.6458-6525; aqueous sanity (aqsanity) l.6547; floorspecies l.6758; verdict l.6770-6800
  (pm.MK = 2 soft fail, W21IPM message; StrandedElementCheck); catch TError / std::exception -> Rescale, rethrow
  as E91IPM; RescaleSystemFromInternal; OptimaTwinLabel; pm.FX = TotalGibbsEnergy() with X/XF/XFA save-restore
  (TotalGibbsEnergy mutates X/XF/XFA!); FI1/FI1s count for GEMSGUI record size.
- dn/db READ-OUT l.4515-4527: only from the PRIMARY solve; optima_dndb (row-major (L+R) x N), optima_dndb_rows/cols
  members on TMultiBase. A retry solves a different problem, so dn/db is that of the primary attempt even when the
  answer came from a retry => correctness hazard for 07 (sensitivity must be recomputed at the accepted answer).
- TRAP 11 evidence: post-solve KKT sign check l.5858-5960: gradJ = F[j] (- tau/X barrier) - sum U[i] A[i,j];
  degenerate box -> 0; barrier-held -> max(-g,0); at lower -> max(-g,0); at upper -> max(g,0); interior |g|;
  tol = pa_GAS. Absolute cutoffs: comment l.6381 "absolute cutoff ... does not see a trace-IC loss of ~1e-13 mol";
  l.6462 massBalanceBadIC absolute min(DHBM*1e10, 1e-2) mol, plus a relative per-IC test |C_i|<=B_i*DHBM added.
  Optima side: confirm in ../optima (bound error zeroed).
- CalculateEquilibriumStateHOP l.6874-7083: native leg CalculateEquilibriumState (pm.pNP warm/cold) with fallbacks:
  E06IPM + g_ipmNoProgFired -> Optima warm from last state; warm fail -> cold native; both fail -> cold Optima.
  Optima leg runKinetics=false; on throw restores Y,X,U,FX,B,controls and recomputes (activity, concentrations),
  MK=2, W21IPM, BAD_GEM_HOP. HopLegGuard (optima_hop_leg), HopLegRecorder destructor writes hop_split + DECIDE
  hop-legsplit. => member flags optima_hop_leg / hop_split are cross-call hidden state on TMultiBase.
- NativeOptimaFinish l.7085: gated pa_NativeOptimaFinish>0 && trace_growth_stuck; snapshots 13 arrays + FX;
  Budget RAII sets g_optimaIterBudget = cap and ~Budget sets -1 (NOT the previous value: nesting would lose an outer
  budget; fragile). Accept only if MK==PZ==0 and Gopt < Gnat - 1e-8; else restore all (DECIDE nativefinish).
  NativeOptimaRescue l.7157: two attempts (warm with hop leg, then cold), same Budget pattern. Callers: node.cpp
  GEM_run_single l.1530-1579.

### ms_multi_format.cpp (1718) — IPM file read/write (trap 1, 2)
- Duplicate physical constants l.38-49: R_CONSTANT=8.31451 (CODATA 1986), F=96485.309, k=1.380658e-23, NA=6.0221367e23,
  H2O_mol_to_kg etc. -> check against m_const_base.h (another copy?).
- MULTI_static_fields (8 rows) l.61-70 and MULTI_dynamic_fields (164 rows; stale comment "//80") l.72-254.
  Checked by script: 164 descriptor rows == 164 enum entries f_sMod..f_pa_OptimaTwinLabel in ms_multi.h l.2056-2098,
  same order, 0 mismatches; every enum is read (172 `case` labels incl. static); 131 written via writeField (the
  rest are arrays via writeArray). TRAP 1 hard-coded counts: `TPrintArrays<TIO> prar( 164, ...)` l.272 and
  `TReadArrays<TIO> rddar( 164, ...)` l.1014; static `prar1( 8, ...)` / `rdar( 8, ...)` l.271/905.
  Template to_text_file_gemipm (l.261) / from_text_file_gemipm (l.861), explicit instantiations l.1706-1714 for
  NlohmannJson OR SimdJson (USE_NLOHMANNJSON) and KeyValue. Write is gated per field by
  `!brief_mode || pa_p->X != pa_p_.X` (pa_p_ = defaults; essentiallyEqual for doubles).
- TRAP 2 mechanism: the writer emits `<END_DIM>` after the 8 static fields (l.341). The reader's static pass uses
  rdar with only the 8 static descriptors; TReadArrays::findNext (io_template.h l.273) loops has_next() and drops
  every unknown name silently; KeyValueRead::has_next (io_keyvalue.cpp l.123-148) reads the next token as a name
  (so a value token is also consumed as a "name") and returns false at END_DIM. A pa_ field put before END_DIM is
  thus skipped with no message. Same silent skip of unknown tags in the dynamic block.
- No runtime setter for pa_ on TNode: settings live in TMultiBase::pa_standalone (shared_ptr<BASE_PARAM>, ms_multi.h
  l.1138), shared between nodes made by copy; `virtual BASE_PARAM* base_param() const` (l.1180) hands out a
  NON-const pointer from a const method (anyone with multi_ptr() can write). LocalSettings RAII (l.1188-1202) makes a
  private copy for a retry (used by the rescue layer). 11 fields are "reserved, no effect" kept only for position.
- BUG (non-obvious fact): from_text_file_gemipm, EZ block (~l.946-956): `for( ii=0; ii<dCH->nDC; ii++ ) if(
  dCH->ccIC[ii] == IC_CHARGE ) break;` loops over nDC on the nIC-sized ccIC array (out-of-bounds read when nDC > nIC,
  always true), then iZ is used as IC column. Works only because the charge IC is found before nIC is passed.
- from_text_file_gemipm also: set_def(); sizes from DATACH; pm.LO = last DC of phase 'a'; pm.PG = size of g/p/f
  phase; pm.E reset if last IC is not charge (fix "09.12.2025"); multi_realloc; copies A, Awt, SB, ICC, MM, ...;
  restores dCH->DCmm when zero (WRITES into DATACH while reading MULTI).

### ms_multi_diff.cpp (481) + ThermoFun load (node_copy.cpp l.617) — trap 4
- pa_p_ (const BASE_PARAM defaults) l.6-23: POSITIONAL aggregate init of the first 37 fields (PC..tprn); the other
  ~200 fields take default member initializers. => BASE_PARAM must stay an aggregate; any X-macro must keep the
  member order; pa_p_ is the "default" against which brief_mode writes (ms_multi_format).
- TMultiBase ctor: pa_standalone = copy of pa_p_; null-inits AA, BB, arrL, arrAN, U_mean/U_M2/U_CVo/U_CV, ICNud,
  phSolMod, phSorpMod, phKinMet, Guns, Vuns, tpp_G/S/Vm; set_def().
- 23 alloc_* functions (IPx, PMc, DMc, MoiSN, SitFr, DQFc, PhLin, lPhc, xSMd, IsoPc/Sc/Ct, EImc, mCDc, xSKrC, ocPRkC,
  feSArC, rpConC, apConC, AscpC, UMpcC, xICuC): each `delete[]; new[]` of a raw MULTI pointer -> manual memory.
- testMulti(): returns 1 if pm.MK || pm.PZ (warns at PSM>=2 to TNode::ipmlog_file). GUI hooks as empty virtuals:
  calculateActivityCoefficients_scripts, initalizeGEM_IPM_Data_GUI, GEM_IPM_Init_gui1/gui2, multiConstInit_PN
  (pm.PZ = pa DW!). GEMSGUI overrides these (to check).
- DC_LoadThermodynamicData(aNa): node1 or aNa; tries na->load_all_thermodynamic_from_thermo(TK,P) (ThermoFun), else
  load_all_thermodynamic_from_grid; pm.pTPD = 2.
- load_all_thermodynamic_from_grid: T,P, Psat at P<1e-5, RT = R_CONSTANT*T, denW/epsW/denWg/epsWg [5] by grid or
  LagranInterp (orders 6/5/4/3 hard-coded per property); per DC: Go from dCH->G0 grid (interp order 6), tpp_G[j] =
  Go, pm.G0[j] = ConvertGj_toUniformStandardState(Go + Guns[j]), Vol (and A[xVol] for VOL_CONSTR), S0, H0, Cp0, A0,
  U0. TRAP 4: DC_G0() returns the raw grid value (CSD->G0, before Guns and the standard-state conversion); the
  solver uses pm.G0 (converted). With ThermoFun, load_all_thermodynamic_from_thermo writes pmm->G0/tpp_G from
  thermo_engine and does NOT touch CSD->G0, so DC_G0() returns the stale exported grid. node.h l.950 documents it;
  Get_solver_G0 (node.h l.954) returns pm.G0.
- NON-OBVIOUS: load_all_thermodynamic_from_thermo uses `static thread_local double funP, funPg, funPge;` because
  ThermoFun's property-cache key stores a REFERENCE to the P lvalue (std::tuple<double,double,double&,string>); a
  local would leave dangling keys and a later lookup could return another substance's data. ThermoFun also writes
  Psat back into funP at P=0. => hidden state forced by a dependency; ThermoFun v2 should fix the key.
- Aqueous phase assumed at CSD->ccPH[0] (index 0) in the ThermoFun path; water key = DCNL[pmm->LO].
- USE_THERMO_LOG blocks (GemsSettings::log_thermodynamic) write CSV logs (thermodynamic-log*.csv).

### nodearray.cpp / nodearray2.cpp (1659) — parallel model
- `static TNodeArray* na` (nodearray.cpp l.37, nodearray.h l.205): process-wide pointer to the LAST constructed
  array (set in both ctors l.51/65); RunGEM(Mode,nNodes,...) l.291 uses `na->` -> two arrays in one process clash.
- OpenMP (useOMP) CalcIPM_List l.142: `#pragma omp parallel shared(...) private( wrkNode, n )` AND an inner
  `TNode wrkNode( *calcNode )` -> the private() copy is default-constructed per thread and unused (waste), plus an
  unused outer copy. Each thread copies calcNode once and REUSES that work node for all nodes of its chunk: the
  MULTI state (Y, U, lnGam, FitVar, and TMultiBase members) of node ii-1 seeds node ii. Default static schedule
  => the chunk boundaries, and so which node seeds which, depend on the thread count => SOP/SHP (warm, which use
  pm.U/pm.Y of the work node, not the DBR) can give thread-count-dependent iteration counts / answers.
  CalcNodeServer l.548: uPrimalSol only for Mode<0 or NEED_GEM_SIA (not SOP/SHP).
- CalcIPM_Node l.194/259: success test `RetCode==OK_GEM_AIA || RetCode==OK_GEM_SIA` only -> OK_GEM_AOP/SOP/HOP/SHP/ROP
  are reported as errors by TNodeArray (BUG for Optima modes; check callers in xgems / GEMSGUI / node-gem).
- CopyWorkNodeFromArray l.53-130: field-by-field copy DATABR (scalars "//6", "//13", "//19" counts in comments),
  NO_NODEARRAYLEVEL ifdef for transport fields; `wrkNode->pCNode()->aPH = nullptr` when nAalp == 0 (drops pointer).
  MoveWorkNodeToArray l.131 mirror. Transport helpers (MoveParticleMass, grids, VTK, log profiles).
- RunGEM(wrkNode, iNode, Mode, nodeArray) l.529: copy in, NodeStatusCH = |Mode|, GEM_run(uPrimalSol), copy out
  (always, also on failure -> stale values copied back, trap 6 at array level).

### ms_multi_copy.cpp (1242) — copyMULTIData layout contract (NON-OBVIOUS)
- copyMULTIData l.37: `pa_standalone = otherMulti.pa_standalone;` (settings SHARED by shared_ptr between copies,
  e.g. every OMP work node shares one BASE_PARAM; LocalSettings makes a private copy for a retry);
  `copyValues( &pm.N, &otherPM->N, 39 )` and `copyValues( &pm.TC, &otherPM->TC, 53 )` copy MULTI scalars by
  memory layout. Checked: the long int block N..SolveCallCount has 40 fields (ms_multi.h l.683-728) -> the last
  (SolveCallCount) is not copied; the double block TC..YFk has 52 scalars followed by FitVar[5] -> the copy takes
  FitVar[0] but not FitVar[1..4] (FitVar[3] is the F0 smoothing state). Adding a scalar to MULTI shifts these
  silently (same class of bug as the BASE_PARAM counts, trap 1).
- None of the newer TMultiBase members (prev_scfact, optima_control_conditions, phase_floor_parked,
  zeroed_floor_amount, warm floor seed, hop_split, trace_growth_*, DUL_preScale_, prevBulkMol_, kinDegenWarned_,
  phfTwin...) are copied: a copied node starts with them empty.
- multi_realloc l.355-1019 (~660 lines of new[] by flags), multi_kill l.1020 (delete[]).

### Activity models s_solmod*.cpp (18 240 lines) — structural scan + targeted reads
- Class per file: s_solmod2 (6230): TCGFcalc, TCORKcalc, TPR78calc, TPRSVcalc, TSRKcalc, TSTPcalc (fluid EoS);
  s_solmod3 (3679): TBerman, TCEFmod, TGuggenheim, TIdeal, TMBWmod, TMargules, TNRTL, TRedlichKister, TRegular,
  TSubregular, TVanLaar, TWilson; s_solmod4 (4811): TDavies, TDebyeHueckel, TELVIS, THelgeson, TKarpov,
  TLimitingLaw, TShvarov; s_solmod5 (2835): TEUNIQUAC, TPitzer, TSIT; s_solmod6 (187): TModOther.
- TSolMod base (s_solmod.cpp 498): ctor(SolutionData*) ALIASES MULTI arrays ("Direct access to index list and
  parameter coeff arrays!"): x = arWx, lnGamma = arlnGam, lnGamConf/Recip/Ex/DQF, CTerm, aIPx, aIPc, aGEX, aPparc,
  aMoiSN, aSitFR, aDCc, DC_Codes, phVOL, phXF, aVol, DQFcf -> the model writes straight into MULTI; ctor zeroes
  lnGamConf..CTerm in MULTI. Owns only aIP (new[NPar]) and multisite y/mn/mns (new[] nests). R_CONST(8.31451)
  hard-coded per instance (a 4th copy of R). solmod_logger = spdlog::stdout_color_mt("solmod") at static init
  (throws if a logger "solmod" already exists in the registry -> static-init-order hazard for embedders).
- DATA-RACE / WRONG-CACHE (NON-OBVIOUS): TCGFcalc (Churakov-Gottschalk fluid EoS) uses function-local mutable
  statics as memo caches: DIntegral l.1422 `static double TOld,roOld; static double a,b,c,d,e;`, also l.1483,
  l.1525, l.1578 (+KOLD), FWCA l.1758 `static double TOld,roOld,F;` returning the cached F when
  essentiallyEqual(T,TOld) && essentiallyEqual(ro,roOld). The cache key ignores the instance and its parameters,
  so two CG fluid phases, or two nodes in two threads, at the same T and ro can get each other's value; under
  OpenMP it is also a data race. TELVIS::trapzd (s_solmod4 l.4678) `static double s;` is the Numerical Recipes
  stateful integrator (state kept ACROSS calls by design) -> not reentrant. s_solmod5 l.1485/1493 static tables
  ak1/ak2 are non-const but read-only (should be static const).
- Fixed-size scratch: Rho = new double[NComp][11], cfh/cfc new double*[10] (s_solmod2 l.4090, 5076-5078);
  McI_PT_array new double[13] (s_solmod5 l.539). TEUNIQUAC/TELVIS [200] arrays (noted earlier).

### Support files (structural scan)
- jsonconfig.cpp: GemsSettings singleton `static GemsSettings data(settings_file_name)` (l.304), file name
  "gems3k-config.json" (l.286) read on first use; sets logger levels/sinks (rotating file sink, lazy file sink),
  log_thermodynamic flag (USE_THERMO_LOG CSV dumps). home_dir() via getenv HOME / HOMEDRIVE+HOMEPATH, static
  aHomeDir. => a JSON file in the working directory changes logging = hidden input (not a result change).
  gems3k_update_loggers(use_stdout, logfile, level) = public API for embedders.
- v_detail.cpp: gems_logger = stdout_color_mt("gems3k") at static init; thread_local g_errorsHandledByCaller (Error()
  logs at debug instead of error while > 0); static atomic g_calculationCancel + SetCalculationCancel /
  CalculationCancelled; Error() always throws TError; ErrorIf(bool, const std::string&, const std::string&) builds two
  std::string from literals on EVERY call, also when false (hot-path cost if used in loops -> Stage B check).
- No other mutable statics/thread_locals/getenv in s_kinmet, s_activity*, s_sorpmod, io_*, gdatastream,
  gems3k_impex, v_service, num_methods, datach_*, ms_multi_file, node_format, node_kinetics (grep).
- node_copy.cpp GEM_init(lst) l.107: GEMS3KGenerator parses the .lst; binary (GemDataStream datach_from_file +
  read_multi) or text/json (read_dch_format_stream; f_thermofun -> load_ThermoEngine; read_ipm_format_stream);
  Pai/Tai from grid; Fdev1/Fdev2 hard-coded 1e-6; read first DBR; init_into_gems3k(). Returns 0 OK, 1 on TError /
  std::exception (message in ipmlog_error), -1 on unknown exception. Second GEM_init(dch_json, ipm_json, dbr_json,
  thermo) l.234 for in-memory strings (xgems/Python).

### Build files and tests
- GEMS3K/CMakeLists.txt: `file(GLOB SOURCE_FILES *.cpp)` -> one monolithic library; shared AND static are both
  built (each compiles every file: 2x build time). USE_OPTIMA_SOLVER and ENABLE_BENCHMARK_DIAGNOSTICS options are
  declared in the subfolder; USE_THERMOFUN adds USE_SPDLOG_PRECOMPILED PUBLIC. CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS.
  Headers installed by pattern (all *.h public, incl. internal ones such as tnt/jama, gems3k_timing.h).
- Root CMakeLists.txt: C++20 required (l.11); GNU flags `-pthread -Wall -Wno-misleading-indentation
  -Wno-ignored-attributes -Wno-variadic-macros -Wno-deprecated` (no -march, no -ffp-contract); default build type
  Release (-O3 -DNDEBUG for GCC). gnu++20 => GCC default -ffp-contract=fast: no FMA on baseline x86-64, but on
  aarch64 (FMA in base ISA) or with -march=native the compiler fuses a*b+c => a different FP order from the same
  source (cross-platform freeze differences; Stage A must not add -march; Stage B could pin -ffp-contract=off).
  Hard-coded include_directories("/usr/local/include") + link_directories("/usr/local/lib") l.101-102 and
  ("/usr/local/include/eigen3") l.167 (leaks system headers before conda ones). configure_file writes the TRACKED
  GEMS3K/gems3k_version.h (l.149-151). BUILD_GEM2MT -> add_subdirectory(nodearray-gem) which does NOT exist.
- tests/: ONE test (pstall_default_test.cpp, `pa_PSTALL` default reset on re-read). Its data dir
  `${CMAKE_SOURCE_DIR}/nodearray-gem-build/TestRW` is NOT in the repository -> the test cannot run from a clean
  clone. BUILD_TESTING default OFF. Real regression testing lives in the private gems-benchmark (not accessible).
- Folders: Docs (GLOSSARY.md, SOLVER-ADVANCED-SETTINGS.md, SOLVER-ADVANCED-USAGE.md, SOLVER-OPTIONS-GUIDE.md,
  SOLVER-PATH-GRAPH.md/.html, doxy, solmod_params.ods); node-gem (main.cpp, qmake .pro/.pri, Makefile.MinGW);
  gemcalc (main.cpp, .pro); tools (kva2json, optima_test, thermofun_test, thread_test, Cu-Pourbaix data);
  tsolmod4rkt (solmodfactory*, solmodengine*, python binding, demos). qmake .pro/.pri files duplicate the CMake
  source lists (a second build system to keep in sync).

### ../optima (fork, branch optima_develop, head 145aeed "Refactor Recorder::iter and Trace classes"), 11 886 lines
- Interface used by the adapter: Optima::Dims{x, be, c}, Problem(dims){f (ObjectiveFunction std::function), Aex,
  be, xlower, xupper, c, bec}, State{x, ye}, Options{maxiters, convergence.tolerance (default 1e-8),
  convergence.check (ConvergenceCheck = std::function<bool(ConvergenceCheckArgs{dims, F, E, uo, u, result})>),
  linesearch.*, backtracksearch.apply_min_max_fix_and_accept, output.*}, Solver{setOptions, solve(problem,
  state[, sensitivity])}, Result{succeeded, iterations}, Sensitivity{xc, pc, yec, ...}, trace::{file(), state(),
  beginSolve()}. Objective callback writes res.f, res.fx, res.fxx (dense L x L), res.diagfxx, res.succeeded and
  reads opts.eval.fxx, opts.ibasicvars. Ownership: Problem/State own Eigen copies of A, b, bounds (copied from MULTI
  every solve); the objective lambda reads/writes MULTI through `this` (side effects on pm.X, XF, XFA, F, lnGam, ...).
  Sign convention: Optima ye = -GEMS3K U.
- TRAP 11 Optima side: ResidualErrors.cpp l.97-109: `ex(ju).fill(0.0)` (unstable variables) and
  `if(xbs[i] == xbslower[i] || xbs[i] == xbsupper[i]) ex[jbs[i]] = 0.0;` -> a basic variable exactly on a bound gets
  zero optimality error whatever the gradient sign; the comment notes Eigen 5 `select` over index views read freed
  memory (hence a plain loop). The GEMS3K side re-checks the sign (ipm_optima.cpp l.5858-5960). Mass balance
  absolute cutoff: GEMS3K CheckMassBalanceResiduals used on the Optima path (l.6170) with absolute
  min(DHBM*1e10, 1e-2) mol = 1e-3 mol at DHBM=1e-13 (comment l.6462) -> trace ICs invisible; relative per-IC test
  added after.
- Optima global state (NON-OBVIOUS): Trace.cpp `static State s;` (process-wide, NOT thread_local) and
  `static FILE* fp` (GEMS3K_OPTIMA_STEPTRACE_FILE, fallback REAKTORO_TRACE_FILE); counters are incremented
  unconditionally inside the solver (LU.cpp l.142/147 lupartial/lufallbacks, BacktrackSearch.cpp l.138 snaps), so
  with two threads solving at once they race, and the GEMS3K DECIDE records lupartial/snapleaving (computed as
  differences of these counters by LuRecord/SnapRecord in CalculateEquilibriumStateOptima) can count another
  thread's events. Also Exception.cpp l.35 static WarningState, LinearSolverNullspace.cpp l.44 static FILE*,
  l.299 static long nprobe, BacktrackSearch.cpp l.148 static FILE* bpf.

### Consumers (compatibility contract, for 02)
- xgems (develop_optima, c1ea07f; C++20; CI ubuntu/windows/macos): includes <GEMS3K/node.h> only. TNode methods used
  (count): pCNode 30, check_Phase_xCH 24, pCSD 21, Ph_Mass 14, check_DC_xCH 12, cTK 10, cP 10, Ph_Mole 10,
  DC_xDB_to_xCH 10, Ph_Volume 8, Ph_HeatCapacityCp/Entropy/Enthalpy 5, Ph_Present 3, Set_bIC/TK/P, Ph_SatInd,
  Ph_Moles, Ph_GibbsEnergy, Get_cDC, GEM_write_dbr, GEM_run, GEM_read_dbr, GEM_init, DC_name_to_xCH_map,
  DC_name_to_xCH, DC_G0, DC_Cp0/V0/U0/S0/H0/A0, phEnthalpies, check_TP, check_IC_xCH, cVs, cMs, Set_pH_target,
  Set_Eh_target, Set_dul, Set_dll, Set_DC_G0, Ph_name_to_xCH, IC_name_to_xCH, Get_nDC, Get_muDC, Get_gDC, Get_aDC,
  Get_TK, Get_P, Get_ControlCondition_titrant, Clear_ControlConditions, GEM_trace_regimes (TraceRegime struct),
  GEM_rescue_step, DCmm, DCaJI; free function gemsSettings().gems3k_update_loggers. DATABR fields: xDC, bIC,
  NodeStatusCH, dul, dll, Ms, pe, pH, Vs, IterDone, IC, Hs, Gs, Eh; DATACH: ccPH, nDCinPH, nPS, PHNL, ICNL, DCNL, nPH,
  nIC, nDC, ICmm, DCmm. Status codes NEED/OK of AIA, SIA, AOP, SOP, HOP, SHP, ROP. No BASE_PARAM / TMultiBase /
  TSolMod / TNodeArray use. Python binding (pybind11) mirrors ChemicalEngine (pyxGEMS.hpp 2856, pyxdGEMS.hpp 1968).
  NOTE xgems DC_G0 -> trap 4 applies to users of xgems with ThermoFun.
- GEMSGUI (develop_gems3kv5, 1fa83316; C++20, Qt6; CI linux/macos/windows; links installed GEMS3K::GEMS3K[-static]):
  * class TMulti : public TMultiBase, public TCModule (Modules/Submods/ms_multi_new.h l.95) overrides: base_param()
    (returns the project's BASE_PARAM), set_def, multi_realloc, testMulti, DC_LoadThermodynamicData,
    get_PAalp_PSigm, STEP_POINT, alloc_IPx/PMc/DMc/MoiSN/SitFr/DQFc/PhLin/lPhc/... (all 23 alloc_* virtuals),
    plus the empty GUI hooks (GEM_IPM_Init_gui1/2, initalizeGEM_IPM_Data_GUI, calculateActivityCoefficients_scripts,
    multiConstInit_PN) -> these virtuals and their signatures are API. TNodeGUI : TNode overrides packDataBr(),
    write_ThermoFun_format_stream, init_into_gems3k(). TNodeArrayGUI : TNodeArray overrides CalcIPM_List,
    CalcIPM_Node, pVisor_Message.
  * BASE_PARAM binary contract: GEMSGUI DEFINES BASE_PARAM::read/write (declared in GEMS3K ms_multi.h l.537-538) in
    Modules/m_param.cpp l.236-247: `readArray(&PC, 10)` (10 long ints from PC), `readArray(&DG, 28)` (28 doubles
    from DG), 4 chars -> the first 38 members' memory order and types are an on-disk format. SPP_SETTING::read adds
    ver, DCpct[352], &NP 9, Pi 19. v_object.cpp l.765-790 reads the SPPpar record at a fixed 758-byte baseline
    ("BASE_PARAM's on-disk size differs ... corrupting every object read after this one"). GEMSGUI defaults:
    `SPP_SETTING pa_` (l.49) and `BASE_PARAM dfBase[5]` (l.93) are POSITIONAL aggregate initializers;
    ResetTail (l.195-235) copies ~40 newer fields by name; 42 data objects o_pax* (Modules/v_mod.h enum, positional)
    bound by SetPtr(&pa.p.Field) (m_param.cpp l.401-440) and described in Resources/data/vis_od.ini.json.
    => adding a BASE_PARAM field visible in GEMSGUI touches: GEMS3K member, MULTI_DYNAMIC_FIELDS enum, descriptor,
    writeField, read case, 2 counts, SET line, pa_p_ (if in the positional part) + GEMSGUI v_mod.h enum, SetPtr,
    vis_od.ini.json, ResetTail, dfBase/pa_ (if positional) = 9 + 5 sites.
  * GEMSGUI uses NEED/OK/BAD/ERR of AIA, SIA, AOP, SOP, HOP, SHP (+ NEED_GEM_PIA in comments only).
- In-repo consumers: node-gem (GEM_from_MT, GEM_from_MT_time, GEM_to_MT, GEM_step_MT, GEM_restore_MT, GEM_init,
  GEM_read_dbr, GEM_run, IC_name_to_xDB, Ph_SatInd, Ph_name_to_xDB, pCNode, pCSD); gemcalc (DC_A0/Cp0/G0/H0/S0/U0/
  V0, DC_a/c/g/mu/n, name and index maps, GEM_print_ipm/read_dbr/write_dbr/run/init, Get_aDC/cDC/gDC/muDC/nDC,
  Ph_Enthalpy/Mass/SatInd/Volume, cP, cTK); tools (Set/Get pH/Eh targets, Clear_ControlConditions,
  code_error_IPM, description_error_IPM); tsolmod4rkt includes s_solmod.h, datach_api.h, io_* , jsonconfig.h,
  m_const_base.h, num_methods.h, v_service.h, verror.h and builds TSolMod objects from its own SolutionData
  (solmodfactory) -> TSolMod ctor(SolutionData*) and SolutionData layout are API for tsolmod4rkt + its Python module.

### ORCHESTRA C++ (../orchestra_cpp, f88bb1b, 11 083 lines, Apache-2.0, Meeussen; literal Java->C++ translation)
- Expander (Expander.h/.cpp 362): text pre-processor; "objects" = parameterised text templates, expanded
  recursively (define class, append/insert, include, delete, expand references); "static sharedBasket" (global).
- Parser -> ExpressionNode tree (virtual double evaluate(), constant(), optimize(parser) = constant folding);
  ExpressionGraph (189): expressions in 5 evaluation stages (stage1/2 append, 3-5 prepend) = ordered evaluation of
  derived variables from unknowns. MemoryNode/Node hold the variable values per node; NodeType holds the variable
  layout and "static" values shared by all nodes of a type; Calculator::clone and Node::clone for copies (per node
  state is explicit -> thread-friendly design).
- Calculator (774) solves the UnEq group (unknowns + equations) by Newton with a FORWARD-DIFFERENCE Jacobian
  (UnEqGroup::calculateJacobian: offset each unknown, recompute all residuals: n+1 full evaluations per Jacobian),
  step control with nodeBetween() interpolation from the last successful node.
- Fit for GEMS3K activity models: the IDEA (model = data: a declarative expression graph per phase, evaluated in
  stages, values per node in an explicit memory block) fits closed-form excess models (ideal, Margules, regular,
  Redlich-Kister, Van Laar, Guggenheim, Debye-Hueckel/Davies/limiting law, SIT, much of Pitzer). It does not fit
  iterative EoS models (CG, PRSV/PR78/SRK cubic roots, STP), ELVIS (numerical integration), EUNIQUAC with its
  internal tables, or CEF multisite without extensions. ORCHESTRA itself interprets the tree (virtual call per node)
  and uses FD derivatives, so copying it as is would make the activity call (already the dominant cost in
  escalated Optima calls) slower. Useful form: describe models declaratively, then generate/compile C++ with
  analytic or automatic derivatives (reverse mode gives the whole d ln gamma / dn block in ~one evaluation vs N FD
  columns). Recommendation for 04/07: adopt the "model as data + explicit per-node memory" idea; reject runtime
  interpretation and the text Expander (GEMS3K already has JSON input).

## Guidance search (2026-11-03) — verified sources
- C++ Core Guidelines: isocpp.github.io is BLOCKED by the session network policy; rule titles verified from the
  official source file https://raw.githubusercontent.com/isocpp/CppCoreGuidelines/master/CppCoreGuidelines.md
  (repo github.com/isocpp/CppCoreGuidelines): P.11 Encapsulate messy constructs; I.2 Avoid non-const global
  variables; I.3 Avoid singletons; I.4 Make interfaces precisely and strongly typed; I.11 Never transfer ownership
  by a raw pointer or reference; I.22 Avoid complex initialization of global objects; I.25 Prefer empty abstract
  classes as interfaces to class hierarchies; I.27 For stable library ABI, consider the Pimpl idiom; F.8 Prefer pure
  functions; F.20 prefer return values to output parameters; F.21 return a struct for multiple out values; C.9
  Minimize exposure of members; C.20 rule of zero; C.21 rule of five; C.41 A constructor should create a fully
  initialized object; C.131 Avoid trivial getters and setters; C.149 unique_ptr/shared_ptr; R.1 RAII; R.3 raw
  pointer is non-owning; R.11 Avoid calling new and delete explicitly; R.20 unique_ptr/shared_ptr for ownership;
  Con.1 immutable by default; Con.2 const member functions; ES.20 Always initialize an object; ES.23 {}-init;
  ES.45 Avoid magic constants; ES.46 Avoid lossy conversions; Enum.3 class enums; Per.1 Don't optimize without
  reason; Per.2 Don't optimize prematurely; Per.6 No performance claims without measurements; Per.7 Design to
  enable optimization; Per.16 compact data structures; Per.19 Access memory predictably; CP.1 Assume multi-threaded;
  CP.2 Avoid data races; CP.3 Minimize explicit sharing of writable data; E.2 throw to signal failure; NR.5 Don't use
  two-phase initialization; SL.con.1 std::array/vector instead of C arrays; T.1 templates to raise abstraction.
- UK Analysis Function: AQuA Book (HM Treasury) https://www.gov.uk/guidance/the-aqua-book (proportionality,
  assurance throughout, verification & validation, RIGOUR); "Quality assurance of code for analysis and research"
  (Duck Book, ONS) https://best-practice-and-impact.github.io/qa-of-code-guidance/intro.html, checklists
  https://best-practice-and-impact.github.io/qa-of-code-guidance/checklists.html (version listed 2023.2).
- xSDK community policies https://xsdk.info/policies (versions differ: 0.4.0 M1 CMake options...M16 install;
  0.6.0 M14 64-bit, M16 debug build; check current); BSSw PSIP https://bssw-psip.github.io and
  https://betterscientificsoftware.github.io/PSIP-Tools (incremental goals, Progress Tracking Cards).
- Wilson et al. "Best Practices for Scientific Computing" PLOS Biology 12(1) e1001745 (2014)
  doi:10.1371/journal.pbio.1001745; "Good enough practices in scientific computing" PLoS Comput Biol (2017)
  doi:10.1371/journal.pcbi.1005510.
- de-RSE "Foundational Competencies and Responsibilities of a Research Software Engineer" arXiv:2311.11457 (v3 Aug
  2024), https://de-rse.org/learn-and-teach/competencies/; GI/de-RSE template guideline (quality management,
  testing) https://www.se-rwth.de/publications/GI-and-DE-RSE-template-guideline-for-the-efficient-development-of-research-software.pdf
  (date inconsistent in the document). US-RSE: referenced the German paper at US-RSE24; no US-RSE guideline found.
- Agner Fog optimization manuals https://www.agner.org/optimize/ (keep data used together, sequential access,
  automatic vs explicit vectorization). Data-oriented design: Acton CppCon 2014
  https://isocpp.org/blog/2015/01/cppcon-2014-data-oriented-design-and-c-mike-acton ; SoA vs AoS.
- std::execution (P2300R10) adopted into C++26 WD at St. Louis July 2024 (Sutter trip report
  https://herbsutter.com/2024/07/02/trip-report-summer-iso-c-standards-meeting-st-louis-mo-usa/); std::simd P1928
  targets C++26 (https://isocpp.org/files/papers/P1928R15.pdf) - final merge not confirmed by search -> treat as
  "C++26, compiler support immature" in 07.
- Kokkos Kernels batched dense (small matrices, Serial/Team/TeamVector, nestable in parallel regions)
  https://kokkos.org/kokkos-kernels/docs/API/batched/dense-index.html ; batched gesv on GPUs (P3HPC 2025)
  https://sc25.supercomputing.org/proceedings/workshops/workshop_pages/ws_p3hpc104.html ; batched sparse
  https://www.osti.gov/biblio/1962268 ; SYCL batched iterative (Intel) https://ar5iv.labs.arxiv.org/html/2308.08417.
- Reaktoro: ChemicalSystem / ChemicalState / EquilibriumSolver / EquilibriumSpecs / EquilibriumConditions /
  EquilibriumSensitivity (dndw) https://reaktoro.org/api/classReaktoro_1_1EquilibriumSensitivity.html ; local
  ../reaktoro: one folder per layer (Common, Core, Equilibrium, Kinetics, Models/ActivityModels, Math, Water,
  Serialization, Singletons, Transport, Extensions), each file with .py.cxx binding and .test.cxx next to it;
  ActivityModel = Model<ActivityProps(ActivityModelArgs)>, ActivityModelGenerator = Fn<ActivityModel(SpeciesList)>,
  chain() to compose; autodiff::real throughout.
- MAGEMin (Riel et al. 2022, G-cubed doi:10.1029/2022GC010427): parallel C library + Julia interface; LP +
  extended Partitioning Gibbs Energy + gradient local minimisation; parallel over single-point calculations (CPU).
- Cantera: ThermoPhase objects shared between MultiPhase and vcs_MultiPhaseEquil; equilibrate() front end chooses
  ChemEquil / VCS (https://cantera.org/dev/cxx/da/df6/classCantera_1_1vcs__MultiPhaseEquil.html).
- PHREEQC/IPhreeqc (Charlton & Parkhurst 2011, C&G): all globals/statics of PHREEQC moved into the IPhreeqc class,
  C functions became methods, small string-based API, in-memory data exchange
  (https://water.usgs.gov/water-resources/software/PHREEQC/IPhreeqc.pdf) - the same move GEMS3K needs for its
  thread_locals.
- OSQP embedded code generation (Banjac, Stellato et al. 2017) https://web.stanford.edu/~boyd/papers/osqp_embedded.html
  (parametric problem -> generated library-free C, factorisation caching, warm start); FORCES (code generation,
  structure exploitation). Clarabel, CEA, IPOPT: no search hit; local ../solvers/cea has CMake + AGENTS.md;
  IPOPT TNLP known pattern (user implements a callback class) - cite as "not checked online".
- ThermoFun v2 (../thermofun 824cd14): docs/explanation/legacy-and-v2.md: v2 format named methods/parameters with
  units (YAML/JSON); legacy positional JSON kept; schedule 1.0 (v2 documented, old calls via Compat/, deprecated)
  -> 1.x (legacy read logs a message naming the converter) -> 2.0 (Compat/ removed, legacy reader stays). Uses
  autodiff::real two-pass (T and P seeded) for derivatives (docs/explanation/autodiff.md). => model for the
  BASE_PARAM / IPM-file migration and for AD in activity models.
- Local ../solvers: Equilipy 7a29c5e (Python), MAGEMin 1298456 (C, Makefile), cea 2f79a64 (CMake, presets),
  phreeqc3 992c7ec (CMake + autotools), qpax a5721f8 (Python/JAX QP), thermochimica 1c076b2 (Fortran, CMake).

## PROFILE: GEMS3K_TIMING per-part summary over the survey (prof/timing_summary.txt; warm modes n=6 incl. seed)
- Activity-coefficient calls ("act", inclusive time of CalculateActivityCoefficients) = 15-63 % of call time in
  EVERY mode and project; native AIA/SIA 16-46 % (CASHNK1 44-46 %), Optima 14-63 %.
- Optima calls: objective + Hessian ("obj") 13-74 %; FD Hessian columns ("fdmix") dominate escalated calls
  (lb_test SIA/SOP/SHP: fdmix 59-62 %, act 60-63 %, 0.46-3.7 M activity calls per call).
- Dimension-reduction pre-solve dominates large cold Optima calls: CemGEMS-keyvalue / Export-condensed AOP
  presolve 85 % of 555-608 ms (16 409 activity calls per call), SOP 69-70 %; mines23 AOP 57 % (+ level 22 %).
- PotentialSpaceFinish + TPD dominate CASHNK1 AOP: psf 89 %, tpd 79 % of 62 ms (AIA 1.75 ms).
- Native IPM modes: "native" ~99 %; inner split needs callgrind (running: prof/cg).
- Native modes that fall into Optima through the rescue/finish layers: Thermo-time-in and test01 AIA 80-83 %
  optima (7.9-14.4 ms vs 1.7 ms SIA), mines23 AIA 72 % optima (591 ms).
- Activity call cost ~18 us (CemGEMS: 0.30 s / 16 409); small systems 2-6 us.
- Variation: per-call wall time over 5 repeated identical cold calls varies ~5-10 % (CemGEMS AOP 0.528/0.558/
  0.585 s) -> only instruction counts (callgrind) are deterministic; wall time needs >=5 repeats, median.

## PROFILE: callgrind pass 1 (whole program; prof/cg, summary prof/cg_summary.txt; shares = inclusive Ir / GEM_run Ir)
- Build profiled: hub-install/gems3k (Release -O3, no -g line info for the lib, symbols present), Optima fork,
  conda gcc 16.2, Eigen 3.4.0. Driver prof_driver (5 calls; warm modes 1 seed + 5).
- Native AIA/SIA (sw_test, CalciteIC, CASHNK1, CemGEMS): activity 16-46 %; GX incl. golden-section
  OptimizeStepSize 19-37 % (GX itself is ~90 % of golden: every probe = full X update + TotalPhasesAmounts + log
  sums); linear solve 10-55 % (CemGEMS AIA 55 %: SolveIpmAugmented QR 54.8 %, Eigen ColPivHouseholderQR 37.7 %
  because La >= pa_IpmQrSortPivot = 100); PrimalChemicalPotentials 2-3 %; PSSC + StabilityIndexes < 0.5 %;
  simplex AutoInitialApproximation 0.2-14 % (sw_test AIA 13.7 %).
- Optima AOP/SOP: Optima::Solver::solve 9-98 %; activity 15-60 %; PrimalChemicalPotentials 1-10 %; pre-solve
  84 % (CemGEMS AOP: 5.7 G Ir for one call vs 0.2 G Ir AIA = 28x); PotentialSpaceFinish 85-91 % with NativeTpdPhase
  73-80 % (CASHNK1 AOP/SOP).
- Diagnostics cost when unset: native_trace* 0.0-0.1 %, CERT 0.0 % -> zero-cost claim holds.
- malloc/free/new and printf shares in pass 1 include GEM_init and the driver -> re-measured with
  --toggle-collect='TNode::GEM_run(bool)' (prof/cgt, running).
- Static init (dl-init: spdlog loggers, ThermoFun, iostreams) = 30 % of the whole CalciteIC AIA program (23.6 M Ir),
  ~0.5 GEM_run worth: matters for short-lived tools, not for transport loops.

## PROFILE: callgrind pass 2 (--toggle-collect='TNode::GEM_run(bool)', self Ir by category, % of GEM_run Ir)
sw_test_AIA 31.3M: ipmQR 18.2, log 17.6, GX+TotalPhasesAmounts 15.4, other 13.0, simplex 11.6, conc+gamma 11.5,
  chemPot 6.3, alloc 3.2, exp 1.4, models 0.8
sw_test_AOP 46.4M: conc+gamma 25.0, Optima 15.7, Eigen FullPivLU 15.0, other 13.7, chemPot 9.9, log 8.4, alloc 4.5
CalciteIC_SIA 50.1M: ipmQR 23.1, log 20.6, GX+TPA 19.4, other 13.4, conc+gamma 12.5, chemPot 6.5, alloc 1.2
CalciteIC_SOP 28.2M: conc+gamma 26.2, other 16.6, Optima 15.1, chemPot 9.8, FullPivLU 9.5, log 8.8, alloc 5.8
CemGEMS_AIA 201.5M: Eigen QR/gemv 41.8, ipmQR 11.4, GX+TPA 11.0, log 9.5, other 8.2, conc+gamma 6.5, simplex 2.0
CemGEMS_AOP 5709M: conc+gamma 26.9, FullPivLU 15.2, chemPot 13.7, other 11.5, log 7.6, Optima 7.3, memset 4.8
  (dense res.fxx.setZero per Hessian), pow 4.4, exp 3.5, models 2.8
CASHNK1_AIA 86.7M: models 21.9 (TBerman dGref_dysigma/ExcessPart/KronDelta), log 19.8, GX+TPA 16.5, ipmQR 11.0
CASHNK1_HOP 135.5M: models 26.7, log 15.2, GX+TPA 10.7, conc+gamma 10.1, ipmQR 7.1, alloc 5.1
=> The activity-model formulas are cheap (0.8-2.8 % self) except multisite Berman (22-27 %). The "activity call"
   cost is mostly the wrapper around the models: CalculateConcentrations(+InPhase), PhaseSpecificGamma, phase_bcs,
   IS_EtaCalc for ALL phases on every call (also for one FD column that changes one phase), plus log/exp/pow.
=> glibc dispatches __ieee754_log_fma (FMA variant chosen at run time by CPU features): libm results can differ
   between CPUs with and without FMA -> freeze comparisons across machines need the same CPU class (08).
=> printf/strtod inside GEM_run: 0 %. Allocation: 0.8-5.8 % self.

## gems-benchmark (access granted after the first hand-over; develop dc220de, main 5651c52)
- Clone /home/user/gems-benchmark (shallow), develop checked out DETACHED (no local branch). Default branch `main`
  lacks the freeze tools; `develop` has GUIDE.md, tools/ (freeze.sh, freeze_diff.py, recheck.py, ci.sh, ...),
  tests/ (96 CTest tests registered here; GUIDE says 42), regression/, Resources/ (259 .lst; corpora gems3k 40,
  gems3k-fail 15, gems3k-native-fail 5, gems3k-psina 28, gems3k-proposed T1-T12, gems4r-corium, solvus, under_review).
- FREEZE NOT RUNNABLE HERE: freeze.sh l.59-62 hard-codes BM=/home/dmiron/git/hub/gems-benchmark and
  CMP=${FREEZE_CMP:-/home/dmiron/git/hub/GEMS3K/debug-optima-vs-reaktoro/mode_compare}; `mode_compare` is not in
  GEMS3K (code search: 0 hits) nor in gems-benchmark. GUIDE "Known limits": the freeze cannot run in CI.
  build-gems3k.sh also hard-codes SRC/CONDA paths of the owner's machine. No recorded freeze file in the repo.
- Built tools+tests (build-solvus/, ignored by .gitignore) against hub-install/gems3k (Optima ON, ThermoFun ON,
  BUILD_SOLMOD OFF, conda gcc 16.2, Eigen 3.4.0). FOUND: include/difftest/optima_regression.h l.608, l.643 use
  ::getpid() without <unistd.h> -> compile error with gcc 16.2 / this libstdc++; worked around with
  CMAKE_CXX_FLAGS=-include unistd.h (no source change).
- trace_selfcheck.sh: OK (RUN 3, BULK 3, SET 3, EFF 3, KEY 3, DECIDE 21 records).
- CTest (Optima build, -j3, 91 s): 92/96 pass. FAIL: optima_regression.aop.restart (1 check:
  earlynet_f_CalcDolo "cost is the unarmed solve plus the wasted probe" armed=18 unarmed=18 set=-4),
  proposed.aop (hoknob_T2b degeneracy penalty ratio 1.0x, ratchet >= 20; crossing_T11 "aop is penalised
  approaching the boundary"), dimreduce.gate (T8_aq101 AUTO -> OFF below the gate; T8_aq201 AUTO -> ON above the
  gate; wrong-side cost 1.0x), ci.baseline (invalid here: baseline is for a no-Optima library; 34 findings).
  Brief: "44 tests, 4 known failures (2026-10-17)" - which 4 not stated; owner to confirm.
- tools/ci.sh (stock GEMS3K, no Optima, built from /home/user/GEMS3K head): first configure FAILED - GEMS3K's
  default BUILD_SOLMOD_PYTHON=ON needs pybind11 even with BUILD_SOLMOD=OFF (cmake/modules/GEMS3KFindDeps.cmake l.10)
  -> plan A1 item. With --gems3k-cmake -DBUILD_SOLMOD_PYTHON=OFF: 132 rows, 28 findings (26 ANSWER: Vs rel
  5e-8..9e-5 on 9 projects; phase amounts e.g. CSHK 4.15e-6 -> 5.07e-6 mol (18 %), Gibbsite 1.08e-4 rel,
  feldspars ~1e-4 rel; 2 COST: f_GEOTHERM ITG 687 -> 688, o_Solvus series2 ITG 372 -> 398). No fixture
  mismatch (input files identical, FNV hashes in "# fix" lines).
- Baseline header: recorded 2026-10-20T14:24:17Z, gems3k 5.2.0, "built by gcc 14.3 ... floating-point answers
  depend on this - the record is not portable across toolchains". This session: gcc 16.2 -> discriminating run
  with gcc 14.3 (conda-forge 14.3.0) on the same head started (build-ci-gcc14).
- Each configure of GEMS3K rewrote the tracked GEMS3K/gems3k_version.h (restored with git checkout every time;
  the stop hook flagged it as an uncommitted change) -> confirms A1 item.
