GEMS3K commit: 54b5f819ba83589daa9cd9aa0a9ff699adb28b1d (branch develop_optima, 2026-11-03)

# 12 — Performance techniques (C++ optimisation catalogue)

This file collects the techniques that can make GEMS3K faster, and says for each one what it attacks (measured
share from 03), how it changes floating-point results, and where it is in the plan. 06 (Stage B) has the ordered
steps; 07 has the extensions; this file is the reference behind them.

**FP class** of a technique:
- **I — bit-identical**: same operations in the same order (allowed in Stage A, still checked by the freeze);
- **F — FP-order change**: the same mathematics, different rounding (Stage B, scored against the nudge spread);
- **R — result change**: a different algorithm or an approximation, answers move within a tolerance (Stage B,
  scored on G and composition, opt-in where it is an approximation).

## 0. Where the time is (from 03, for reference)

| Area | Share | Measured on |
|---|---|---|
| failing calls through the rescue ladder | ×1 000–13 000 a normal call | lb_test warm modes |
| finite-difference Hessian columns (`fdmix`) | 13–62 % of time in Optima modes | GEMS3K_TIMING parts |
| work around the activity models (concentrations, γ conversion, phase sums) | 25–27 % of instructions (Optima), 6–13 % (native) | callgrind pass 2 |
| activity-model formulas themselves | 0.8–2.8 % (Berman multisite 22–27 %) | callgrind pass 2 |
| linear algebra | 10–53 % | callgrind pass 2 |
| GX + phase sums, log | 11–19 % + 10–21 % (native) | callgrind pass 2 |
| dimension-reduction pre-solve | 57–85 % of large cold Optima calls | GEMS3K_TIMING parts |
| allocation, memset | 0.8–5.8 %, memset 4.8 % | callgrind pass 2 |

## 1. C++ basics

| Technique | Attacks | Class | Step |
|---|---|---|---|
| One workspace per solver object instead of `new[]`/`std::vector` per call (compressed index, simplex buffers, TNT `Array2D`, DATABR snapshots, Hessian `setZero` only on touched blocks) | allocation 1–6 %, memset 4.8 % | I (if every buffer is written before it is read) | A7, B9 |
| Cache `log(X)`, `log(XF)` of the current iterate; reuse them in GX, chemical potentials and the line search | log 10–21 % | I if the cached value is the same expression; F otherwise | B8 |
| No virtual call per species: one virtual call per phase that loops inside | activity wrapper | I | B15 (PhaseModel API) |
| Move invariant work out of loops (T-P terms, model constants) | activity models | I | B15 |
| LTO; PGO | 0–5 % (not measured) | I in principle (no reassociation without `-ffast-math`), but inlining changes contraction on FMA targets → treat as F and score | B10 |
| Explicit FP contract (`-ffp-contract=off`) and no `-march` in release builds; an optional tuned build | reproducibility | F (one time) | B10 |
| `ErrorIf(bool, const std::string&, …)` builds strings even when false → overloads with `const char*` | small, hot paths | I | A12 |

## 2. Vectorisation

- **What vectorises without changing results:** element-wise loops with no cross-element dependency, e.g.
  `X[j] = Y[j] + LM*MU[j]` (GX), weight updates, scaling, copies. GCC/Clang do this at `-O3` for SoA data without
  reassociation. Check with `-fopt-info-vec` and keep the FP class I.
- **What changes results:** reductions (GX sum, mass-balance residuals, dot products in the dual potentials, phase
  sums) need reassociation to vectorise (`-ffast-math`/`-fassociative-math` or explicit SIMD); vector math libraries
  (`libmvec`, SVML, SLEEF) round `log`/`exp`/`pow` differently from scalar glibc → class F, Stage B only.
- **Data layout:** MULTI is already structure-of-arrays (one array per quantity, species contiguous), which suits
  SIMD (Agner Fog: keep data used together close, access sequentially). The species → phase index (`L1`) gives
  contiguous blocks per phase. Keep it.
- **Tools:** Eigen 3.4 (already used; vectorised kernels), `std::simd` (C++26 target, not yet on all consumer
  compilers, 04 §7), xsimd or Highway as an interim choice. No hand-written intrinsics in solver code.
- **Where it pays:** batched evaluation over many nodes (§7) more than within one small node (N ≈ 5–30, L ≈ 20–500).

## 3. Solving reduced systems

| Technique | Today | Proposed | Class |
|---|---|---|---|
| Species reduction | `OptimaReducedPreSolve` (ipm_optima.cpp l.1745): active set from LP seeds, readmission by reduced cost | keep; reuse the factorisation and objective state between passes; start the second attempt from the first's active set (B6) | R (path), F |
| Redundant species | `ExcludeRedundantDCs` (held at zero per call) | keep | I |
| Sparse stoichiometry | compressed A index (`Build_compressed_xAN`, \|a\| > 1e-12) rebuilt every solve | build once per system, rebuild only when A changes (B9) | I |
| Active phase set | PSSC inserts/removes phases; all species stay in the linear systems | factorise only over present phases plus candidate phases; keep absent phases out of the KKT system until a stability test asks for them | R |
| Elements with fixed potential | not used | eliminate ICs whose dual is fixed by a pure phase or a control condition (pH/Eh target) before the solve | F/R |
| Species at the floor | stay in the factorisation with tiny weights | take species with `X < DcMinM` (exact zeros, trap 12) out of the factorisation; put them back on the stability test | F |
| Optima canonical form | `Echelonizer`/`EchelonizerExtended::updateWithPriorityWeights`, `updateOrdering`; nullspace and rangespace solvers | reuse the canonical form across warm calls of the same system (A fixed) instead of rebuilding it | F |

## 4. Computing only what changes

| What changes | What is recomputed today | What is enough | Class | Step |
|---|---|---|---|---|
| one species in an FD Hessian column | concentrations, γ conversion and models of **all** phases | the phase that holds the species (and phases coupled by electrostatics) | I where the activity call is idempotent (Optima path: smoothing off); F otherwise | B4 |
| the bulk composition b only (transport, warm modes) | everything from scratch: compressed index, thermodynamic data reload test, factorisation structure | A and its index are unchanged: reuse index, symbolic factorisation, ordering; reload nothing | I/F | B9, B7 |
| T or P | all standard states (grid or ThermoFun), all model PT parameters | ThermoFun v2 rule: "the engine calculates again only what depends on a change"; `CheckMtparam` already gates the reload by `Ttol`/`Ptol`; ThermoSol recomputes `PTparam` only on a T-P change (11 §1) | I | A15, B15 |
| a few weights in the IPM normal equations | full rebuild and factorisation of `AA` every iteration | low-rank updates do not apply (all weights change each IPM iteration); use them only in the line-search retries and in MBR polishing where few weights move | F | B7 |
| step length λ in the golden-section search | full `GX(λ)`: X update for all species, phase sums, all logs | `GX(λ)` along a fixed direction: phase sums are linear in λ, only the log terms need X; precompute per-phase direction sums | F | B8 |
| one phase's amount (PSSC insertion/elimination) | full `TotalPhasesAmounts`, chemical potentials of all species | update that phase's totals and potentials | I/F | B8 |

## 5. Memoisation

| Memo | Today | Proposal | Keys and rules |
|---|---|---|---|
| objective value at an identical x (line search) | `memoize_objective_if_linesearch` (ipm_optima.cpp l.120) | keep | bit-identical x (correct: the objective is not pure, H8) |
| standard-state properties at (T, P) | ThermoFun property cache (with the reference-key defect, 10 §5); grid interpolation (to be removed) | an exact-value cache per (T, P) in `ThermoFunSource`; for transport with many temperatures, an **on-demand table**: ThermoFun values computed at the T-P points that occur, kept in a deterministic map | bit-exact (T, P) keys → no FP change versus ThermoFun direct; an interpolating table would be class R and needs its own comparison (10 §3) |
| phase-model PT parameters | `SolModParPT` on T-P change only (via `pm.pTPD`) | the same inside ThermoSol records (11 §1) | per instance, keyed on (T, P) and the record |
| FD Hessian across line-search trials and retries starting from the same point | recomputed | reuse when x is bit-identical | bit-identical x |
| warm-start state per node | `WarmState` on TMultiBase (A11) | per node, never shared between nodes (B11) | node id |
| **Rules** | the H5 counter-example: `TCGFcalc` caches keyed on (T, ρ) only, process-wide, shared by instances and threads (09 §3 fact 3) | every memo: per instance; keyed on **all** inputs; bit-exact keys; deterministic eviction (count-based, never time-based — trap 10); one memo per thread or a lock | — |

## 6. On-demand learning (smart equilibrium)

**Idea** (Leal, Kyas, Kulik, Saar, "Accelerating reactive transport modeling: on-demand machine learning algorithm
for chemical equilibrium calculations", Transport in Porous Media, 2020; Reaktoro `SmartEquilibriumSolver`,
`EquilibriumPredictor`, `EquilibriumSensitivity`, `SmartKineticsSolver` in `../reaktoro/Reaktoro/Equilibrium` and
`Kinetics`; references from memory, not re-checked online in this session):

1. Store solved states with their sensitivities (∂n/∂b, ∂n/∂T, ∂n/∂P).
2. For a new input, find a stored state close to it (nearest neighbour; Reaktoro clusters by the phase assemblage).
3. Predict: n ≈ n₀ + ∂n/∂b·(b − b₀) + ∂n/∂T·(T − T₀) + ∂n/∂P·(P − P₀).
4. Accept the prediction if a cheap test passes (Reaktoro: chemical potentials within a tolerance; for GEMS3K: the
   CERT record (mass balance, charge) and the stability test of the absent phases); otherwise run the full solve
   and store the result.

**For GEMS3K:**

- Needs dn/db (and dn/dT, dn/dP) **at the accepted answer** (07 §3; today dn/db comes from the primary Optima
  attempt only, 09 §3 fact 5 → B1) and exact mixing blocks for the native path (B16).
- Plug-in point: a run layer before `runSingle` (04 §3.2), active only when the caller asks (new option on
  `TNodeArray` or on `GEM_run`; no new NEED code unless the owner wants one).
- Reproducibility: a deterministic store (insertion order = call order, fixed search order, no time-based
  eviction), so a run repeats exactly; results differ from full solves within the acceptance tolerance → class R,
  opt-in, scored like a solver change and reported per mode.
- A cheaper first step with no acceptance risk: **warm start from the nearest stored state** (class R only through
  the iteration path; the answer is a full solve).
- Expected gain: in transport sweeps, most calls are close to an earlier one (not measured here; Reaktoro reports
  large speed-ups for reactive transport — to be measured on the owner's transport cases).

## 7. GPU and batched execution

Summary of 07 §1: GEMS3K calls are small (0.1–20 ms) and branch-heavy, so a GPU helps only with batches of nodes
that take the same path. **The first GPU candidate is batched ThermoSol evaluation** (`Engine::evaluateBatch`,
11 §4): the same phase model on many nodes, SoA over nodes, no branching, and it is the largest share in Optima runs
(§0). Next: batched standard states from ThermoFun (10 §5 batch call), then the common path of the native solver
(MBR + IPM iterations with a fixed phase set), with nodes that need phase selection or retries leaving the batch.
Backend: Kokkos (Kokkos Kernels batched dense solvers) or SYCL; decided only after the CPU batch experiment of
07 §1.

## 8. Parallel CPU over nodes

After A6–A9 and B2 (no shared mutable state), `TNodeArray` runs one work node per thread with
`WarmSeedPolicy::OwnPreviousState` so results do not depend on the thread count (B11). Parallelism inside one node is
not worth it at these sizes.

## 9. Order (summary)

Stage A (class I only): A7 workspace, A15 ThermoSource, A16/A17 one phase-property implementation and model
registry. Stage B, in order of measured share: B3 (failure tail) → B4 (only the changed phase) → B5/B16 (derivative
blocks) → B6 (pre-solve reuse) → B7 (linear algebra reuse) → B8 (GX/log caching) → B9 (allocation) → B10 (flags,
LTO) → B11 (parallel nodes) → on-demand learning (opt-in) → GPU batch (after the CPU batch test).
