GEMS3K commit: 54b5f819ba83589daa9cd9aa0a9ff699adb28b1d (branch develop_optima, 2026-11-03)

# 06 — Stage B: speed and behaviour changes (scored like solver changes)

Stage B starts only after Stage A is complete and released (the Stage A tests and diagnostics are the
instruments of Stage B). Every B step can change floating-point order or behaviour, so it is scored like a solver
change, never by row-for-row identity.

## 0. Scoring protocol for every B step

1. **Two freeze pairs** (brief §5): STD vs STD and RAW vs RAW, all six modes, the full corpus (72–87 projects).
2. **Order of verdict** (`freeze_diff.py`): a lost answer outranks any saving; G, Vs, Ms must not move on a row that
   converged in both; per-mode scorecard; exit 1 on a lost or changed answer.
3. **Answers are scored on G and composition**, not on OK/FAIL alone (brief §4: two assemblages 1e-6 RT apart can
   differ in phases; a 29-point composition error cost 2.3e-4 RT). Eh is not scored.
4. **Iteration counts are judged against the known spread, not one draw.** Before a B step that changes FP order,
   run the 1e-15 relative nudge set (brief §4: the GEMS3K IPM solver's count is a lottery on 20 of 42 projects;
   Optima-mode status lotteries on `10THb_351`, `case-275`, `case-276`) on base and change: report the median and
   the range per project; a change "saves" only if the median moves outside the base range.
5. **Cost is reported in two units**: instructions inside `GEM_run` (callgrind, deterministic) on the 8-case set of
   03 §3 plus 5 owner-chosen projects, and wall time (median of ≥ 5 runs, one machine class).
6. **All six modes always** (brief §8: one optimisation that was 8.8 % cheaper in one mode cost +27 % over six,
   SHP +364 %, and lost 2 answers).
7. Each step lists the **expected gain from 03** and states it as a hypothesis; a step that does not reach half of
   its expected gain is reconsidered before merging (Core Guidelines Per.6).

Step card fields: *Change* · *Basis (03)* · *Expected gain* · *Risk* · *Scoring notes* · *Size* · *Depends on*.

## 1. Order

```
B0 nudge baseline ─▶ B1 behaviour fixes ─▶ B2 thread safety ─┬─▶ B3 escalation cost
                                                              ├─▶ B4 activity-call scope
                                                              ├─▶ B5 analytic/AD Hessian blocks (largest)
                                                              ├─▶ B6 pre-solve cost
                                                              ├─▶ B7 linear algebra
                                                              ├─▶ B8 GX / golden section
                                                              └─▶ B9 allocation
B10 compile/FP policy (owner decision) ; B11 node batches (after B2) ; B12 libm reproducibility (optional)
```

## 2. Steps

### B0 — Nudge baseline
- **Change:** none in code. Record, per project and mode, the spread of iterations and status under 10 nudges of
  1e-15 (the existing experiment), for STD and RAW, on the Stage A release.
- **Basis:** brief §4. **Gain:** none; it is the yardstick for B3–B10.
- **Size:** owner's machine time. **Depends on:** A done.

### B1 — Behaviour fixes found in this review (each needs an owner decision, 08)
- **Change** (one PR each):
  - `TNodeArray::CalcIPM_Node`: accept every `OK_*` code (today only AIA/SIA; Optima OKs are reported as errors).
  - `copyMULTIData`: copy the scalar blocks by name (today 39 of 40 long ints and FitVar[0] only; 09 §3 fact 1).
  - `from_text_file_gemipm`: loop `ii < dCH->nIC` over `ccIC` (today `nDC`, an out-of-bounds read).
  - dn/db: recompute the sensitivity at the accepted answer when a retry tier produced it (09 §3 fact 5), or flag
    the stored matrix as "from the primary attempt".
  - `NativeOptimaFinish`/`Rescue` `Budget`: restore the previous `g_optimaIterBudget` (now a `RunContext` field)
    instead of −1.
  - `DC_name_to_xCH_map`: loop over the DCH phase count, or map the DBR index (node2.cpp l.135).
- **Basis:** 09. **Gain:** correctness; node-array users see OK for Optima modes.
- **Risk:** GEMSGUI/node-gem behaviour changes (B1 first item); FitVar copy can change warm results of copied nodes.
- **Scoring:** full protocol; node-array freeze rows (new) for the first two items.
- **Size:** 6 × S. **Depends on:** A done.

### B2 — Thread safety of models and solver (pre-condition for B11)
- **Change:** if A9 changed the freeze, do it here (scored); make the Optima trace counters thread-local in the fork;
  remove the remaining process-wide scratch; TSan clean run of `tools/thread_test.cpp` on 4 threads × 6 modes.
- **Basis:** 09 §3 facts 3, 4. **Gain:** correctness under threads. **Size:** S–M. **Depends on:** A6–A9.

### B3 — Escalation cost (the tail)
- **Change:** (a) one `IterationBudget` per `GEM_run` across all layers and tiers (today budgets are per layer and the
  rescue ladder multiplies them); (b) stop a warm call early when the stall/FD signature of lb_test appears
  (0.46–3.7 M activity calls per call) and go to the cold path once, not through every tier; (c) keep the FD
  Hessian of the last accepted iterate across line-search retries and attempts that start from the same point.
- **Basis:** 03 §1 finding 3 (lb_test warm: 4.6–35 s per failing call vs 2.7 ms normal; fdmix 59–62 %).
- **Expected gain:** −50 % or more of total corpus time where failures exist; no change on clean rows.
- **Risk:** lost answers (some rescues succeed late). **Scoring:** lost answers decide; report rescues won/lost
  per tier from the `# dec` lines.
- **Size:** M. **Depends on:** A6 (RunContext budgets), B0.

### B4 — Scope of one activity call
- **Change:** per-phase dirty flags in `CalculateActivityCoefficients(LINK_UX_MODE)`: an FD column that changes one
  species of phase k recomputes concentrations, γ conversion and the model of phase k only (plus phases coupled
  by electrostatics, which already disable `fdPureSkip`).
- **Basis:** 03 §3: concentrations + γ conversion = 25–27 % of instructions in Optima runs; act = 15–63 % of time.
- **Expected gain:** 15–25 % of Optima-mode instructions on multi-phase systems.
- **Risk:** the activity call is not idempotent (H8: lnGmo accumulation, F0 smoothing via FitVar[3]). Smoothing is
  off on the Optima path (`optima_disable_smoothing`), so skipping unchanged phases should be bit-identical there;
  on the native path keep full calls. Verify bit-identity first; if bit-identical on the corpus, the step can be
  merged with the Stage A rule (owner's choice).
- **Size:** M. **Depends on:** A10 (models behind PhaseModel), A7.

### B5 — Analytic or automatic derivatives for the mixing blocks
- **Change:** `PhaseModel` returns `d lnγ / d n` for its phase (analytic where written, forward-mode AD as in
  ThermoFun/Reaktoro (`autodiff::real`) elsewhere); the Optima objective uses the block instead of FD columns
  (`fdmix`, `fdblock`), the PotentialSpaceFinish Newton too. Start with Ideal, Margules, Regular, Redlich–Kister,
  Van Laar, Guggenheim (closed form), then Berman/CEF (multisite, the expensive one in CASHNK1), then aqueous
  (DH/Davies/Helgeson, Pitzer, SIT). Models without a block keep FD columns.
- **Basis:** 03 §2 (fdmix 59–62 % in escalations; one activity call per FD column); 07 §3.
- **Expected gain:** up to the FD share of the call (20–60 % in Optima modes); fewer escalations if the Hessian is
  more exact (hypothesis).
- **Risk:** largest FP change of the plan; Hessian quality changes Optima paths (status lotteries). Needs B0.
- **Size:** L (per model family M). **Depends on:** A10, B4.

### B6 — Pre-solve cost on large cold Optima calls
- **Change:** reuse the factorisation and the objective state between passes; start the second attempt from the
  first attempt's active set instead of from scratch; skip the pre-solve when the LP seed's support already
  covers ≥ x % of species (threshold measured); measure `pa_OptimaPreSolveFirstIters` choices on the corpus.
- **Basis:** 03 §2: presolve 57–85 % of large cold Optima calls (CemGEMS AOP 555 ms, 16 409 activity calls).
- **Expected gain:** −30…−60 % on such calls. **Risk:** changed Optima paths; status lotteries.
- **Size:** M. **Depends on:** A10c, B0.

### B7 — Linear algebra
- **Change:** Optima: `FullPivLU` (9.5–15 % of Optima-run instructions) → partial pivoting where the fork's
  `OPTIMA_LU_PARTIAL_PIVOT` path shows equal answers (the macro exists; `pa_OptimaPartialPivLU`); native: reuse
  QR workspace (A7), blocked Householder for La ≥ 100, check the `pa_IpmQrSortPivot` threshold (Eigen QR is 42 %
  of CemGEMS AIA); JAMA (MBR, fallback) → Eigen LDLT/PartialPivLU.
- **Basis:** 03 §3. **Expected gain:** 5–25 % depending on size. **Risk:** pivot order = FP order (brief §4: Eigen 5
  vs 3.4 already differ). **Size:** M. **Depends on:** A7, A10b, B0.

### B8 — GX and the golden-section search
- **Change:** cache per-phase sums and logarithms for the base point Y; evaluate GX(LM) on the step direction
  incrementally; stop the golden section on the existing criterion with fewer full evaluations.
- **Basis:** 03 §3: GX + TotalPhasesAmounts 11–19 % plus most of the 10–21 % log share in native runs.
- **Expected gain:** 10–20 % of native-mode instructions. **Risk:** truncation below DcMinM (trap 12) must happen at
  the same place; FP order of sums changes. **Size:** M. **Depends on:** A10b, B0.

### B9 — Allocation
- **Change:** build the compressed A index once per system and on change only (today every solve); allocate
  DATABR snapshots lazily (today every call in `GEM_run_tries`, 3 KB–MBs); simplex NEW/WORK buffers and TNT Array2D
  from the workspace; `res.fxx.setZero()` only on the touched blocks (memset 4.8 % in CemGEMS AOP).
- **Basis:** 03 §3 allocation 0.8–5.8 %, memset 4.8 %. **Expected gain:** 2–8 %. **Risk:** low (no FP change
  expected; verify bit-identity first, as in B4). **Size:** S–M. **Depends on:** A7.

### B10 — Compile and FP policy (owner decision, 08 Q11)
- **Change:** state the FP contract in CMake: `-ffp-contract=off` (or `on`) explicitly for all platforms, so x86-64,
  aarch64 and Apple clang builds compute the same; decide on `-march` for release builds (portable default,
  optional tuned build); try LTO.
- **Basis:** NOTES "Build files" (gnu++20 default contraction; FMA in the aarch64 base ISA; xgems/GEMSGUI CI build
  on macOS). **Gain:** reproducibility across platforms; LTO +0–5 %. **Risk:** every row's iteration count can move
  (lottery). **Size:** S. **Depends on:** B0.

### B11 — Parallel node batches
- **Change:** `TNodeArray` with one work node per thread **and** `WarmSeedPolicy::OwnPreviousState` (each node
  seeds from its own previous answer, kept in its DATABR/WarmState), so results do not depend on the thread count;
  deterministic schedule; then the batched interfaces of 07 §1.
- **Basis:** 09 §3 fact 2 (today warm modes depend on the thread count). **Gain:** near-linear scaling over nodes for
  transport users; deterministic results. **Risk:** different warm seeds → different iteration counts (scored).
- **Size:** M. **Depends on:** B2, A11.

### B12 — Cross-CPU reproducibility of libm (optional)
- **Change:** offer a build option that uses a correctly rounded `log/exp/pow` (for example CORE-MATH) instead of
  glibc's CPU-dispatched variants, for freeze machines.
- **Basis:** 03 §4 (`__ieee754_log_fma` chosen at run time). **Gain:** freezes comparable across machines.
  **Cost:** some speed (measure). **Size:** S. **Depends on:** B0.

## 3. Measured basis summary

| Step | Measured share it attacks (03) | Where measured |
|---|---|---|
| B3 | failing warm calls ×1 000–13 000 the normal call; fdmix 59–62 % | survey, timing parts (lb_test) |
| B4 | concentrations + γ conversion 25–27 % of instructions (Optima) | callgrind pass 2 |
| B5 | FD columns: fdmix 13–62 % of time in Optima modes | timing parts |
| B6 | pre-solve 57–85 % of large cold Optima calls | timing parts (CemGEMS, mines23) |
| B7 | FullPivLU 9.5–15 %; Eigen QR 42 % (CemGEMS AIA); hand QR 11–23 % | callgrind pass 2 |
| B8 | GX + phase sums 11–19 %, log 10–21 % | callgrind pass 2 |
| B9 | allocation 0.8–5.8 %, memset 4.8 % | callgrind pass 2 |
| B10–B12 | reproducibility, not speed | build analysis, libm dispatch |
