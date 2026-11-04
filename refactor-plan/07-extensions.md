GEMS3K commit: 54b5f819ba83589daa9cd9aa0a9ff699adb28b1d (branch develop_optima, 2026-11-03)

# 07 — Extensions: GPU, methods from other sectors, dn/db sensitivities

For each extension: the plug-in point in the target architecture (04), what Stage A/B must prepare, and the
evidence. None of these is part of Stage A.

## 1. GPU and batched execution

**What it would do.** Solve many nodes at once (transport: 10³–10⁶ nodes per step, same chemical system,
different bulk composition, T, P). The work per node is small (0.1–20 ms on one CPU core, 03 §1), so a GPU helps
only with **batches** of nodes that run the same code path.

**Evidence that limits the idea.**
- The solve path is branch-heavy: phase selection, retry tiers, rescue ladder (01 §2–§4). Nodes in one batch
  diverge after a few iterations. Escalated calls cost ×1 000–13 000 (03 §1); one such node stalls its batch.
- The hot arithmetic per iteration is small dense algebra (N ≈ 5–30 ICs, L ≈ 20–500 species): QR/LU of size N
  or (L + N), dot products, logs (03 §3). That is the size range of batched small solvers
  ([Kokkos Kernels batched dense](https://kokkos.org/kokkos-kernels/docs/API/batched/dense-index.html): Serial/Team
  routines nestable in a parallel region; batched gesv on GPUs,
  [P3HPC 2025](https://sc25.supercomputing.org/proceedings/workshops/workshop_pages/ws_p3hpc104.html)).
- The activity models are the other hot part, written per phase in virtual C++ classes with aliases into MULTI
  (H7) and process-wide statics (H5): not portable to a device as they are.

**Plug-in point.** `gems3k_kernel` gets a batch entry beside the single-node entry:
`solveNativeBatch(BatchView, BatchContext&)` where `BatchView` is structure-of-arrays over nodes (one array per
MULTI field, nodes contiguous: data-oriented design). The batch kernel runs the **common path only** (MBR + IPM
iterations + GX line search with a fixed phase set); a node that needs phase selection, a retry or the Optima path
leaves the batch and is solved by the normal CPU path. Backend choice: Kokkos (C++17, CUDA/HIP/SYCL/OpenMP back
ends) for the batch kernel; `std::execution`/`std::simd` (C++26) are not yet usable on all consumer compilers
(04 §7).

**Preparation in A/B.** A6 (RunContext: no globals), A7 (workspace), A9/B2 (model statics removed), A11/B11
(warm seeds per node, deterministic), A10 (kernel separated from node), B4/B5 (activity block interface with
explicit inputs/outputs — the device needs it). First experiment after B11: CPU-only batch (SoA, OpenMP SIMD) on
the corpus to measure divergence (share of nodes leaving the batch per iteration) before any GPU code.

**Decision proposed.** Do not start GPU work before the CPU batch experiment shows ≥ 70 % of node-iterations on
the common path for the owner's transport cases (08 Q12).

## 2. Methods from other sectors

| Method (source) | What it is | Fit for GEMS3K | Plug-in point | Prepare |
|---|---|---|---|---|
| Declarative models + generated code with derivatives (ORCHESTRA idea; OSQP/FORCES code generation; Reaktoro/ThermoFun autodiff) | describe a mixing model as data, generate C++ with analytic/AD derivatives | **Good** for closed-form excess models; removes FD Hessian columns (B5) | `models/` registration table (04 §3.9) | A10, B4 |
| Levelling / Partitioning Gibbs Energy (MAGEMin, [Riel et al. 2022](https://doi.org/10.1029/2022GC010427)) | LP over pseudo-compounds then local minimisation | GEMS3K already has the leveled start and column generation (ipm_level.cpp, `pa_IpmLevelStart`, `pa_OptimaLevelStart`); **extend**, do not add a third copy | `kernel/level.cpp` | A6e (outputs as structs), A12 (one LP solver) |
| Warm-started parametric QP/IPM (OSQP factorisation caching; [Banjac, Stellato et al.](https://web.stanford.edu/~boyd/papers/osqp_embedded.html)) | reuse factorisation when only b changes | **Good** for transport sweeps (SIA/SOP: same A, small Δb) | `numerics/` workspace; Optima fork options | A7, B7 |
| Smart equilibrium / on-demand learning (Reaktoro `EquilibriumPredictor`, `SmartEquilibriumSolver`) | predict the new state from a stored state and its sensitivity; accept if a cheap check passes | **Promising** for transport; needs exact dn/db at the accepted point (§3) and a certified acceptance test (CERT exists) | `run/` layer before `runSingle` | §3, B1 (dn/db at final state), B11 |
| Homotopy / continuation (common in process simulation) | follow the solution from an easy problem to the target | Already present as the rescue "walk" (pa_RescueWalkSteps); **keep**, make it a separate layer | `run/rescue` | A10d |
| Interval or certified checks (verified computing) | certify phase stability | CERT/TPD exist; **no change** | `kernel/certify` | — |

## 3. dn/db sensitivities

**State today.** With `optima_want_sensitivity = true` the Optima adapter sets `dims.c = N`, `problem.c = B`,
`problem.bec = I` and reads `Sensitivity::xc` into `optima_dndb` (row-major (L+R)×N, members
`optima_dndb_rows/_cols`) — **from the primary solve only** (ipm_optima.cpp l.4515–4527; 09 §3 fact 5). The
native path has none. No public `TNode` getter exposes it (to check with the owner; xgems does not use it).

**Evidence of value.** Reaktoro exposes `EquilibriumSensitivity` (dn/db, dn/dT, dn/dP;
[API](https://reaktoro.org/api/classReaktoro_1_1EquilibriumSensitivity.html)); it is used for implicit coupling in
transport (Newton on the transport–chemistry system) and for smart prediction (§2).

**Plug-in point.** A `Sensitivity` result on `SolveOutcome` (04 §3.3), filled by every solver leg at the
**accepted** state:
- Optima leg: one extra linear solve with Optima's sensitivity solver at the final state (the KKT matrix is
  available after convergence) — recompute when a retry tier produced the answer (B1).
- Native leg: from the converged IPM normal equations: with the active set fixed, dn/db solves
  [H Aᵀ; A 0]·[dn; dy] = [0; I] restricted to non-zero species — the same structure as `SolveIpmAugmented`; needs
  the exact Hessian of G (the mixing blocks of B5; with FD columns only an approximation).
- At phase boundaries dn/db is discontinuous; report the active set with the matrix and flag "near boundary" when
  a phase amount is within `pa_DcMin`-scale of zero (CERT already computes margins).

**Preparation.** B1 (dn/db at the accepted point), B5 (exact mixing blocks), A10c (Optima tiers as functions, so
"which tier produced the answer" is a return value), 02 §6 migration rules for the new public getter
(`TNode::GEM_sensitivity_dn_db()` returning a matrix view and the active set).

**Verification.** Finite-difference check on the corpus: perturb each b_i by ±1e-7 relative, solve warm, compare
with dn/db·Δb on rows with an unchanged phase set; report the relative error per project (expected < 1e-5).

See also: `10-thermofun-tp.md` (standard states from ThermoFun only), `11-thermosol-thermokin.md` (phase-model
and kinetics libraries; the batched `Engine::evaluateBatch` is the first GPU candidate) and
`12-performance-techniques.md` (techniques catalogue, on-demand learning §6). Systems made in code or YAML without
GEMSGUI (`13-headless-system-definition.md`) make batch, GPU and on-demand-learning tests easy to set up.

## 4. Order (after Stage B core steps)

1. dn/db at the accepted point for the Optima leg (B1) → public getter (small).
2. Mixing-block derivatives (B5) → native dn/db.
3. CPU batch experiment (after B11) → GPU decision (08 Q12).
4. Smart prediction layer (needs 1, 2 and the CERT acceptance test).
