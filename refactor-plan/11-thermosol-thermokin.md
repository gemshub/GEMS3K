GEMS3K commit: 54b5f819ba83589daa9cd9aa0a9ff699adb28b1d (branch develop_optima, 2026-11-03)

# 11 — ThermoSol (phase models) and ThermoKin (kinetics): libraries outside the solver

**Goal (owner):** the activity and mixing models leave GEMS3K and form a library, **ThermoSol**, built like ThermoFun
v2 but for phase properties. **ThermoFun v2 is the guideline** (how to organise records, methods and formats).
**`TActivity` shows what moves**: it is where GEMS links the phase models today. In the same way **`TKinetics` becomes
ThermoKin**, the kinetics library. GEMS3K stays the equilibrium solver that uses ThermoFun (10), ThermoSol and
ThermoKin.

```
            ThermoFun (standard states, T-P)        thermohub/thermofun
                 │ G0…U0, water properties
     ┌───────────┴───────────┐
 ThermoSol (phase models)   ThermoKin (kinetics)     inside GEMS3K first (thermosol/, thermokin/), own repos later
     │ lnγ, excess props,       │ rates, bounds, surface areas
     │ d lnγ/dn (later)         │
     └───────────┬────────────┘
             GEMS3K kernel (GEMS3K IPM solver, Optima adapter, run layers)
```

## 1. ThermoFun v2 principles, applied to phase models

From thermofun `docs/explanation/data-model.md`, `legacy-and-v2.md`, `autodiff.md` and `ThermoFun/Methods/`:

| ThermoFun v2 principle | ThermoSol / ThermoKin form | Replaces in GEMS3K |
|---|---|---|
| Records come first; formats are only ways in and out | a **phase record**: species list, model id + implementation, named parameters with units, T-P coefficients, sources. A **kinetic record**: process, rate law, parameter sets, affinity terms, surface-area model | positional arrays `sMod[8]`, `LsMod`, `LsMdc`, `PMc`, `DMc`, `IPx`, `MoiSN`, `DQFc`, `LsKin`, `xSKrC`, `ocPRkC`, `feSArC`, `rpConC`, `apConC`, `UMpcC` in MULTI and the IPM file |
| Values are SI inside; a reader converts units | parameters carry units in the v2 file; legacy arrays are read with GEMS's implicit units | implicit units per model |
| Missing is not zero (`notdefined` status) | a parameter not given is absent, the model reports `notdefined` | zero-filled arrays |
| A record does not change; the engine binds it to method objects once and recomputes only what depends on a change | `Engine::bind(record)` once per system; T-P dependent parameters (`PTparam`) recomputed only when T or P changes; per-phase results recomputed only when that phase's inputs change | `SolModCreate` + `PTparam` + `MixMod` on every call; all phases on every call (03 §3: 25–27 % of Optima-run instructions are this work) |
| A method belongs to a record, not to a class | the record names its model; no switch on a character code | `switch (sMod[SPHAS_TYP])` in `SolModCreate`/`SolModActCoeff`/… and `KM_Create` (≥ 8 edit sites per new model, 09) |
| Water passes its properties on | the solvent's density and dielectric constant (with derivatives) are an explicit input of the aqueous models, from ThermoFun | `new TPitzer( &sd, aM, aZ, pm.denW, pm.epsW )` raw pointers (ipm_chemical3.cpp l.740–761) |
| Constants are data (`Methods/ParameterSets`) | model constants (Debye–Hückel A/B coefficients, Pitzer universal terms, …) in parameter sets | constants in code (`s_solmod4/5.cpp`) |
| `MethodRegistry`: one factory per method id and implementation | `ModelRegistry`: one line registers a phase model or a rate law | the switch tables |
| `Methods/Symbolic` (`SymbolicExpression`) | **symbolic excess models**: an excess Gibbs function given as an expression in the record, with automatic derivatives (the ORCHESTRA "model as data" idea, 04 §7) | hard-coded formulas for every new simple model |
| autodiff two-pass (T, P seeded) | derivatives of lnγ with respect to composition (d lnγ/dn), T and P | finite differences in the Optima adapter (03: fdmix up to 62 %) |
| `Compat/` + staged deprecation 1.0 → 1.x → 2.0 | `Compat` reads and writes the legacy GEMS arrays byte for byte; v2 YAML/JSON is the documented format | — |
| Statuses and errors are not part of the calculation | `PhaseProps::status`, set after the calculation | `ErrorIf` inside model code |

## 2. What moves into ThermoSol (from `TActivity`)

`TActivity` (activities.h, s_activity.cpp, s_activity2.cpp, s_activity3.cpp) and its TNode front end
(node_activities.cpp) already separate "phase properties at given T, P and composition" from the solver.

| Part of `TActivity` | Goes to | Notes |
|---|---|---|
| `ACTIVITY` struct (activities.h) — the phase-model subset of MULTI | `thermosol::PhaseState` (inputs) and `thermosol::PhaseProps` (outputs) | inputs: T, P, amounts per species, phase totals, solvent properties, ionic strength; outputs: lnγ, Gex/Hex/Sex/Cpex/Vex/Aex/Uex, ideal and standard terms, Darken terms, site fractions, optional derivative blocks |
| `SolModCreate`, `SolModParPT`, `SolModActCoeff`, `SolModExcessProp`, `SolModIdealProp`, `SolModStandProp`, `SolModDarkenProp`, `Alloc/Free_TSolMod`, `Alloc/Free_TSorpMod` (s_activity3.cpp) | `thermosol::Engine` (bind, PT update, evaluate) | the registry of §1 replaces the switches |
| `CalculateConcentrations`, `CalculateConcentrationsInPhase`, `phase_bcs`, `phase_bfc`, `bfc_mass` (s_activity2.cpp) | `thermosol` (concentration scales are part of the phase model: molality, mole fraction, surface density) | per phase → only the changed phase is recomputed (12 §4) |
| `DC_PrimalChemicalPotential(Update)`, `PrimalChemicalPotentials`, `TotalPhasesAmounts`, `ConvertDCC`, `ConvertGj_toUniformStandardState`, `PhaseSpecificGamma`, `SetSmoothingFactor` | split: the chemical potential of a species in its phase → ThermoSol; the GEMS-specific smoothing and units → GEMS3K adapter | the F0 smoothing history (H8) stays a GEMS3K solver feature, outside ThermoSol |
| TSolMod family (s_solmod.h, s_solmod*.cpp: ~30 models, 18 240 lines), `TSorpMod` (stub) and the sorption code in TMultiBase (ipm_chemical2.cpp) | `thermosol/models/` (fluid, solution, aqueous, sorption) | memo statics removed (H5, A9) |
| tsolmod4rkt (`SolModFactory`, `SolModEngine`, Python module) | ThermoSol loaders and Python binding | tsolmod4rkt is merged, not kept beside |
| TNode API `setTemperature`, `setPressure`, `updateStandard*`, `setSpeciation`, `updateConcentrations`, `init/updateActivityCoefficients`, `getIntegralPhaseProperties`, `updateChemicalPotentials`, `updateTotalGibbsEnergy`, `updateActivities`, `updateChemicalData` (node_activities.cpp) | kept on TNode as the front end, forwarding to ThermoSol | public API (02); now documented |
| **Stays in GEMS3K** (solver, not phase properties): `GX`, `DC_GibbsEnergyContribution`, `StabilityIndexes`, `KarpovCriterionDC`, `KarpovsPhaseStabilityCriteria`, `DC_DualChemicalPotential`, `CheckMassBalanceResiduals`, `Build_compressed_xAN`, `Alloc_A_B` | `gems3k_kernel` | — |

## 3. What moves into ThermoKin (from `TKinetics` and `TKinMet`)

| Part | Goes to | Notes |
|---|---|---|
| `KINETICS` struct and `TKinetics` (kinetics.h, node_kinetics.cpp): `KM_Create`, `KM_ParPT`, `KM_InitTime`, `KM_UpdateTime`, `KM_UpdateFSA`, `KM_ReturnFSA`, `KM_InitRates`, `KM_CalcRates`, `KM_InitSplit`, `KM_CalcSplit`, `KM_SetAMRs`, `KM_InitUptake`, `KM_CalcUptake`, `Alloc/Free_TKinMet` | `thermokin::Engine` | the duplicate in TMultiBase (ipm_chemical4.cpp, `RunKineticsStep` and its KM_* copies) goes; GEMS3K calls ThermoKin |
| `TKinMet` base and `TMWReaKin`, `TUptakeKin`, `TIonExKin`, `TAdsorpKin`, `TNucleKin` (s_kinmet.h l.226–592, s_kinmet.cpp, CRLF): `PTparam`, `RateMod`, `RateInit`, `SplitMod`, `SplitInit`, `SetMetCon`, `UpdateFSA`, `GetModFSA`, `UpdatePT`, `UpdateTime`, `PRrateCon`, `CorrSpecSurfArea`, uptake tables | `thermokin/models/` | `KinMetData` (s_kinmet.h l.97) → kinetic record |
| **Stays in GEMS3K**: applying the result as metastability bounds `DLL`/`DUL` (`Set_DC_limits`), the time step loop of the caller | `gems3k_kernel` | ThermoKin returns rates and bounds; it does not run the solver |

## 4. API sketch (new code, `namespace thermosol` / `thermokin`)

```cpp
namespace thermosol {
// PhaseModel calculates the properties of one phase.
// It keeps no state between calls. The caller owns the workspace.
class PhaseModel {
public:
    virtual ~PhaseModel() = default;
    virtual void updatePT(double T_K, double P_Pa, const SolventProps& water) = 0;  // only when T, P change
    virtual PhaseProps evaluate(const PhaseState& in, Workspace& ws) const = 0;   // thread-safe
    virtual bool hasCompositionDerivatives() const { return false; }              // d lnγ/dn block (B16)
};
// Engine binds the phase records of one system to models, once.
class Engine {
public:
    Engine(const Database& phases, std::shared_ptr<ThermoFun::ThermoEngine> standardStates);
    void bind(const SystemView& system);                       // species and phases of one GEMS system
    void evaluateAll(const AmountsView&, std::span<PhaseProps> out, Workspace&);           // one node
    void evaluateBatch(const BatchView&, BatchProps& out, BatchWorkspace&);                // many nodes (SoA)
};
}
namespace thermokin {
class RateModel {
public:
    virtual ~RateModel() = default;
    virtual void updatePT(double T_K, double P_Pa) = 0;
    virtual Rates rates(const KineticState& in, Workspace& ws) const = 0;   // rates, bounds, surface areas
};
}
```

Inside GEMS3K, `TActivity` (and `TKinetics`) become the adapters that fill `PhaseState` from MULTI and copy
`PhaseProps` back, so MULTI, GEMSGUI and the TNode API keep working (02).

## 5. Steps (also listed in 05 and 06)

| Step | Stage | What | FP class | Check |
|---|---|---|---|---|
| A16 | A | One phase-property implementation: TMultiBase's `CalculateActivityCoefficients` path and `TActivity` share code where it is token-identical; where they differ, keep both with a cross-reference (A12 rule) | none | freeze row for row |
| A17 | A | `ModelRegistry` replaces the `SolModCreate` and `KM_Create` switches (same constructors, same order of calls) | none | freeze; one test per registered model |
| B15 | B | `thermosol/` and `thermokin/` folders, CMake targets and namespaces inside GEMS3K; `PhaseModel`/`RateModel` interfaces; `Compat` readers/writers for the legacy arrays (byte-identical round trip); v2 YAML/JSON records | none intended; moving code across libraries changes inlining (A-rule 3), so scored | freeze (expect identical), instruction counts, IPM round trip |
| B16 | B | derivative blocks: analytic or AD `d lnγ/dn`, d/dT, d/dP per model family (= 06 B5) | FP change | scored (06 §0) |
| B17 | B | move to own repositories (proposed `thermohub/thermosol`, `thermohub/thermokin`), merge tsolmod4rkt, Python modules, a Reaktoro adapter (`ActivityModelGenerator` wrapper) | none | consumers' CI (xgems, GEMSGUI) |
| — | — | `TActivity`/`TKinetics` are **not** deprecated (this reverses 05 A12/A14 and 08 Q14 of the first version): they are the front ends | — | — |

## 6. Order and dependencies

A9 (model statics) and A10 (files into modules) come first; A16/A17 need A2 (golden tests) and A4 (settings
access). B15 needs 10 §2 (`ThermoSource`, so ThermoSol gets standard states and water properties from one place).
B16 is where the speed comes from (03: FD columns and activity calls dominate escalated Optima calls). B17 only
after the v2 record format has been used on the whole corpus for one release.

## 7. Open points (08)

Q25 (where ThermoSol/ThermoKin live first), Q27 (names), Q28 (sorption models: ThermoSol or a third library), Q29
(does ThermoSol own the concentration-scale conversions, or does GEMS3K pass concentrations).
