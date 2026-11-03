GEMS3K commit: 54b5f819ba83589daa9cd9aa0a9ff699adb28b1d (branch develop_optima, 2026-11-03)

# 10 — Standard-state T-P properties from ThermoFun only

**Goal (owner):** GEMS3K takes all standard-state T-P properties from ThermoFun. Its own T-P path is first made
switchable, then compared in full with ThermoFun, and removed when the comparison shows that ThermoFun reproduces it.
The same applies to GEMSGUI, which today computes the grids that GEMS3K interpolates (§6).

## 1. Today: two paths into the same MULTI arrays

| | Grid path (own) | ThermoFun path |
|---|---|---|
| Entry | `TMultiBase::load_all_thermodynamic_from_grid` (ms_multi_diff.cpp) | `TNode::load_all_thermodynamic_from_thermo` (node_copy.cpp l.617) |
| Chosen by | `DC_LoadThermodynamicData` (ms_multi_diff.cpp l.293): ThermoFun if `thermo_engine` is set, else grid | `.lst` file mode `f_thermofun` → `load_ThermoEngine` (node_copy.cpp l.590) |
| Source data | DATACH lookup arrays (datach.h l.50–108): `G0, V0, S0, H0, Cp0, A0, U0` [nDC][nPp][nTp]; `denW, denWg, epsW, epsWg` [5][nPp][nTp]; `Psat` [nTp]; grid `TKval`, `Pval`; `Ttol`, `Ptol`; `mLook` | ThermoFun JSON database; `thermoPropertiesSubstance`, `propertiesSolvent`, `electroPropertiesSolvent` (ThermoEngine.h l.148–182) |
| Method | `LagranInterp` with a fixed order per property (G0 6, V0 5, S0 4, H0 5, Cp0 3, A0 5, U0 5; water 5–6), or the exact node when T and P are on the grid (`check_grid_TP`) | direct calculation at (T, P), ThermoFun property cache |
| Writes | `pm.G0` (after `ConvertGj_toUniformStandardState`), `tpp_G`, `Vol` (+ `A[xVol]` for `VOL_CONSTR`), `S0`, `H0`, `Cp0`, `A0`, `U0`, `denW/epsW/denWg/epsWg[5]`, `P` (Psat when P < 1e-5) | the same MULTI arrays; **not** `CSD->G0` (trap 4) |
| Reload gate | `CheckMtparam` (node_copy.cpp l.563): reload when T or P moved more than `Ttol`/`Ptol`; `pm.pTPD`, `load_thermodynamic_data` | same gate |
| Other users | `TNode::DC_G0 … DC_U0`, `Get_Psat` (node2.cpp, 22 `LagranInterp` calls; trap 4: `DC_G0` returns the grid value); `TActivity::updateStandard*` (node_activities.cpp) call `DC_G0` etc. | — |
| Known issues | interpolation error between nodes (not measured yet, §3 a) | `static thread_local double funP, funPg, funPge` because ThermoFun keeps a reference to P in its cache key (09 §3); `f_PitzerTHE` gets no water density/dielectric (gems-benchmark GUIDE "Known limits") |

**Stays in GEMS3K (not a T-P method):** `ConvertGj_toUniformStandardState` (G0/RT + ln 55.51 for aqueous and
sorption species, + ln P for gases and fluids, + `YOF` surface term), the `Guns`/`Vuns` increments, `fDQF`/DQF terms,
and the molar-volume row of A for `VOL_CONSTR`. They are corrections that GEMS applies to standard-state values,
whatever their source.

## 2. Target: one interface, three implementations

```cpp
// ThermoSource gives the standard-state properties at one T and P.
// The solver does not know where they come from.
class ThermoSource {
public:
    virtual ~ThermoSource() = default;
    // Values for each species in DATACH order (SI: J/mol, J/(mol K), m3/mol).
    virtual StandardProps standardProperties(double T_K, double P_Pa) = 0;
    // Water: density and dielectric constant with T and P derivatives, liquid and vapour.
    virtual SolventProps solventProperties(double T_K, double P_Pa) = 0;
    // Saturation pressure of water at T (Pa).
    virtual double psat(double T_K) = 0;
};
```

- `ThermoFunSource` — wraps `ThermoFun::ThermoEngine` (target default). Optional derivatives d/dT, d/dP come from
  ThermoFun's `ThermoScalar` (`ddt`, `ddp`, autodiff two-pass, thermofun `docs/explanation/autodiff.md`).
- `GridSource` — today's `load_all_thermodynamic_from_grid` and `DC_*` interpolation, moved without change (legacy;
  removed in §4).
- `CompareSource` — calls both and writes one line per species and property with the two values and their
  difference (env `GEMS3K_THERMO_COMPARE_FILE`, zero cost when unset), and returns the grid values, so that results
  do not change while the comparison runs.

`DC_LoadThermodynamicData`, `TNode::DC_G0 … DC_U0`, `Get_Psat` and `TActivity::updateStandard*` call the source.
The reload gate (`CheckMtparam`) stays in front of the source. Water properties go to the phase models as an explicit
input (11 §2), not as `pm.denW`/`pm.epsW` raw pointers.

## 3. Comparison campaign — the gate for removing the grid path

All four parts must pass before step B14 (removal). Results are kept in gems-benchmark (new `tools/compare/thermofun-gems`
outputs, recorded in `recheck.py`).

| Part | What | How | Pass rule (owner decision, 08 Q26; proposed) |
|---|---|---|---|
| a. Property level | every species of every corpus project; at the grid nodes and half-way between nodes (this also measures the grid's own interpolation error) | existing gems-benchmark `tools/compare/thermofun-gems/` (`tfprops.cpp`, `tfcheck.cpp`, `compare.py`, `summary.py`) driven over all `-dch` files; `CompareSource` for the values GEMS3K actually uses | \|ΔG0\| ≤ 1 J/mol; relative 1e-6 for V0, S0, Cp0, H0; water density and dielectric constant with all 5 derivative entries relative 1e-6; or a documented reason per exception |
| b. Equilibrium level | every corpus project in all six modes | the existing `f_`/`j_` twins (11 `f_`, 17 `j_` in `Resources/gems3k`); a ThermoFun twin for every grid-only project, made with GEMSGUI's ThermoFun export (`v_thermofun.cpp`); STD and RAW freezes of both, `freeze_diff.py`; regression diagrams and Pourbaix maps (`regression/run_regression.sh`) | no lost answer; G and composition within the freeze tolerances, or explained by the (a) differences |
| c. Coverage | every T-P method used in the GEMS databases | list the methods used by GEMSGUI (`s_supcrt*.cpp`, `m_dcomp*.cpp`, `m_reacdc.cpp`, `m_reac2.cpp`) and map each to a ThermoFun method (`ThermoFun/Methods/Solutes` (HKFgems, ADgems, Anderson91, HollandPowell98), `Water` (HGK-JN gems, Bradley–Pitzer, Fernández 1997, Helgeson–Kirkham, Sverjensky 2014, …), `Fluids`, `Minerals`, `Reactions`, `Symbolic`) | the list of methods without a ThermoFun equivalent is empty; `f_PitzerTHE` solved first |
| d. Cost | time to load all species at a new (T, P) | callgrind on node-gem / GEMS4R-like T sweeps: grid interpolation versus ThermoFun direct | if ThermoFun is slower by more than the owner accepts, use the (T, P) memo or the on-demand table of 12 §5 |

## 4. Switch and removal schedule (ThermoFun `legacy-and-v2.md` pattern)

1. **A15 (Stage A, no FP change):** `ThermoSource` with `GridSource` and `ThermoFunSource` wrapping today's two
   functions; the choice stays exactly as today (`thermo_engine` set → ThermoFun). Freeze row for row.
2. **B13 (Stage B, scored):** run the campaign (§3). When it passes, ThermoFun is the default for every project
   that has ThermoFun data; GEMSGUI exports ThermoFun data for every project (§6).
3. **B14 (deprecation, then removal):**
   - release N: GEMS3K writers stop writing grid arrays (DCH files with `nTp = nPp = 1` placeholder or without
     grid fields; format decision in 02 §6); readers still accept grids and log one `info` line naming the
     converter (`kva2json`-like tool that adds ThermoFun data);
   - release N + 12 months (not before 5.0): remove `GridSource`, the `DC_*` grid interpolation and the T-P use of
     `LagranInterp`; `num_methods` keeps the generic interpolation only if something else still uses it;
   - the DATACH grid fields stay in the struct until 5.0 (02 §2.2 contract), then go with a migration note.

## 5. Requests to ThermoFun (thermohub/thermofun, owner's repository)

1. Pressure **by value** in `thermoPropertiesSubstance`, `propertiesSolvent`, `electroPropertiesSolvent` (today
   `double& P`, ThermoEngine.h l.148–182); this removes the `static thread_local funP` workaround and its hidden
   state (09 §3).
2. A batch call: all species of one system at one (T, P), returning SoA arrays (one call per T-P change instead of
   one per species; prepares 12 §5 and the GPU batch of 07 §1).
3. Thread safety of one `ThermoEngine` shared by many nodes (today each `TNode` holds its own `shared_ptr`).
4. Vapour-state water: done (GEMS3K a9c65d9 reads steam properties in the vapour state).

## 6. GEMSGUI track (separate pull requests in GEMSGUI)

GEMSGUI computes the grids with its own T-P methods: `Modules/s_supcrt.cpp` and `s_supcrt2.cpp` (2 957 lines,
SUPCRT/HKF and water), `m_dcomp.cpp`, `m_dcomp2.cpp`, `m_reacdc.cpp`, `m_reac2.cpp`, `Submods/ms_mtparm.cpp`,
`s_tpwork.cpp` (about 9 000 lines in total). It already exports ThermoFun data (`Modules/v_thermofun.cpp`, 963
lines). Plan: (1) export ThermoFun data for every project (needed for §3 b); (2) compute the values that GEMSGUI shows
through ThermoFun, with the comparison of §3 a in GEMSGUI's own tests; (3) remove the GEMSGUI T-P methods in the
same schedule as §4. Owner decision 08 Q24.

## 7. What changes for users

- Projects need ThermoFun data after B14 (02 §6 migration: converter, 12-month overlap).
- Results move within the tolerances of §3 (interpolation error disappears); this is a data change, scored as such,
  never mixed with a solver change in one freeze pair.
- Transport codes that relied on fast grid interpolation get the memo/table of 12 §5 if §3 d shows a cost.
