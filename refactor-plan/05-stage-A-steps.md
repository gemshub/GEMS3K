GEMS3K commit: 54b5f819ba83589daa9cd9aa0a9ff699adb28b1d (branch develop_optima, 2026-11-03)

# 05 — Stage A: behaviour-preserving steps

## 0. Rules for every Stage A step

**A-rule 1 — Same answers.** The STD freeze and the RAW freeze of the step must match the freezes of the step
before **row for row**: status, G, Vs, Ms, iterations, and every `# dec` line, in all six modes (AIA, SIA, AOP,
SOP, HOP, SHP), plus `# set`, `# eff`, `# key`. One differing row fails the step. There is no tolerance in Stage A.

**A-rule 2 — Same floating-point program.** No change of operation order, no `-ffast-math`, no FMA contraction
change, no `-march`, no Eigen change (stays 3.4.0), no compiler change, no change of compile flags or
definitions per translation unit. Check: `compile_commands.json` of the step and of the base give the same flags
for every source file that exists in both (script `tools/check_flags.py`, added in A1).

**A-rule 3 — Same inlining scope for hot code.** A step that moves a hot function (03 §3: GX, TotalPhasesAmounts,
SolveIpmAugmented, CalculateConcentrations*, PhaseSpecificGamma, chemical potentials, the Optima objective) to
another translation unit must keep the callers that inline it in the same unit, or show that the freeze still
matches. Report the instruction count inside `GEM_run` (callgrind `--toggle-collect='TNode::GEM_run(bool)'`) on
the 8-case set of 03 §3 before and after; a change above 2 % on any case needs a note in the step record.

**A-rule 4 — Diagnostics unchanged.** RUN/SET/EFF/KEY/CERT/DECIDE records byte for byte on the golden inputs
(A2); probes and environment variable names unchanged; logger names unchanged.

**A-rule 5 — Process.** One step = one branch = one pull request; build in a separate build folder; no build while
a freeze runs (brief §5); check `ldd <tool> | grep -i -E 'gems3k|optima'` before a freeze; check
`CMAKE_HOME_DIRECTORY` of a reused build folder; no commit or push without the owner's request; no
`Co-Authored-By: Claude` trailer (owner rule).

**Step card format** (PSIP style): *Moves* · *Why* · *Risk* · *Verification* · *Size* (S < 200 changed lines,
M 200–1 000, L > 1 000 mostly moves) · *Depends on* · *Other repos*.

## 1. Order of work

```
A0 baseline ─▶ A1 build ─▶ A2 contract tests ─┬─▶ A3 BASE_PARAM single declaration ─▶ A4 settings access
                                               ├─▶ A5 getters / stale outputs
                                               └─▶ A6 RunContext (a–e) ─▶ A7 workspace ─▶ A8 process statics, env
A9 model statics (after A2) ; A10 file moves and splits (after A6, A7, A13) ; A11 warm state (after A6)
A12 de-duplication (after A10) ; A13 diagnostics component (after A2) ; A14 tests and docs (continuous)
```

Expected length: 25–35 pull requests. Steps A6a–e, A10a–d and A12 are the long ones.

## 2. Steps

### A0 — Baseline and harness (owner's machine)
- **Moves:** nothing in code. Reproduce the recorded freeze of the committed head (STD and RAW) with the current
  harness; record callgrind instruction counts for the 8-case set (03 §3) and for 5 owner-chosen corpus projects.
- **Why:** brief §5: "reproduce the recorded freeze of the committed head before trusting a new harness". This
  session could not (gems-benchmark private).
- **Risk:** a non-reproducible baseline (machine class, CPU FMA dispatch of libm, 03 §4). Mitigation: one machine
  class for all Stage A freezes.
- **Verification:** freeze matches the record row for row; if not, stop and find why before A1.
- **Status from this session (2026-11-03, 03 §7):** trace self-check OK; CTest 92/96; `ci.sh` reproduces 114
  of 132 recorded native rows exactly, and the other 18 rows differ for environment reasons (same result with
  gcc 14.3 and 16.2, at the baseline's code and at head; glibc FMA dispatch alone moves one project). The freeze
  could not run (driver `mode_compare` not published). **Before A1:** publish `mode_compare` (and the recorded
  standard freeze) in gems-benchmark, make `freeze.sh`, `build-gems3k.sh` and `recheck.py` take their paths from
  the environment instead of `/home/dmiron/...`, and record in each freeze header the CPU flags, glibc,
  compiler and conda package versions (08 Q21).
- **Size:** none in GEMS3K (S in gems-benchmark). **Depends on:** —. **Other repos:** gems-benchmark (tool
  paths, `mode_compare`, the `<unistd.h>` include of `include/difftest/optima_regression.h`).

### A1 — Build system: explicit targets, no hidden paths
- **Moves:** `GEMS3K/CMakeLists.txt` from `file(GLOB)` to explicit lists grouped by the future modules (04 §2),
  as `OBJECT` libraries feeding the unchanged `GEMS3K` and `GEMS3K-static` targets; remove
  `include_directories("/usr/local/include")`, `link_directories("/usr/local/lib")`,
  `include_directories("/usr/local/include/eigen3")` (root l.101–102, l.167) in favour of `find_package` targets;
  write `gems3k_version.h` into the build tree (stop rewriting the tracked file, root l.149–151); add
  `gems3k_build_features()` (01/02 §1: USE_OPTIMA_SOLVER, NO_NODEARRAYLEVEL) checked by `GEM_init`; drop the
  dangling `BUILD_GEM2MT` subdirectory; make `BUILD_SOLMOD_PYTHON` follow `BUILD_SOLMOD` (today a default
  configure fails without pybind11 even with `BUILD_SOLMOD=OFF`, 03 §7); add `tools/check_flags.py` (A-rule 2); add `.gitattributes` that pins
  `ipm_chemical.cpp` and `s_kinmet.cpp` to CRLF and the rest to LF (trap 9; no content change).
- **Why:** reproducible builds, module boundaries visible, tracked file no longer dirty after configure; the
  hard-coded paths broke the tsolmod4rkt build in this session.
- **Risk:** removing `/usr/local/include` can change which Eigen/spdlog headers are found (an Eigen 5 in
  `/usr/local` would change FP results, brief §4). Mitigation: CMake prints the Eigen version found and fails on
  anything but 3.4.x unless `GEMS3K_ALLOW_EIGEN5=ON`.
- **Verification:** A-rule 2 script (same flags per file); identical object files for unchanged sources
  (`objdump -d` hash of each `.o`, expected equal except the version string unit); freeze (A-rule 1).
- **Size:** M. **Depends on:** A0. **Other repos:** gems-benchmark build script (new build folder name only);
  xgems/GEMSGUI: none (`find_package(GEMS3K)` unchanged).

### A2 — Contract tests and compile-time checks
- **Moves:** new `tests/` targets with data inside the repository (copy the `TestRW` inputs that the existing test
  needs — today missing): (1) `static_assert` on `offsetof`/type of the first 38 `BASE_PARAM` members (02 §2.5);
  (2) `static_assert(std::size(MULTI_dynamic_fields) == f_pa_OptimaTwinLabel + 1)` and the same for the static
  table; (3) asserts that document the MULTI scalar blocks used by `copyMULTIData` (40 long ints from `N`,
  52 doubles from `TC`, then `FitVar`) — documenting today's 39/53 copy, not fixing it; (4) IPM/DCH/DBR write →
  read → write byte-identical round trip on 6 corpus projects; (5) trace golden test: RUN/SET/EFF/KEY/CERT/DECIDE
  of 3 small projects × 6 modes compared byte for byte; (6) a test that every `NEED_*` code has a branch in
  `packDataBr` and in the `TNodeArray` success test (the latter is expected to fail today → marked
  `WILL_FAIL` with a link to B1).
- **Why:** turns traps 1, 9, 13, 14 and the H12 layout contract into failing tests instead of reviewer memory
  (Wilson 2014 "plan for mistakes").
- **Risk:** none for results (tests only). Golden files must be produced on the A0 machine class.
- **Verification:** tests pass on the base commit (except the marked one); freeze unchanged.
- **Size:** M. **Depends on:** A1. **Other repos:** none.

### A3 — `BASE_PARAM` from one declaration
- **Moves:** new `data/base_param_fields.def` with one line per field after the 38-member positional prefix:
  `GEMS3K_PA_FIELD(type, name, default, "pa_tag", "comment", flags)`. It generates: the members (same order), the
  `MULTI_DYNAMIC_FIELDS` entries, the descriptor rows, the `writeField` calls (with the same `brief_mode` /
  `pa_p_` comparison), the reader `case` labels, both counts (`std::size`), and the SET line pieces in
  `native_trace_run_header`. The prefix (PC…tprn) stays hand-written with A2's `static_assert`s. The
  `MULTI_dynamic_fields` order stays identical (it is the writing order).
- **Why:** trap 1: seven edit sites and two hand counts become one line; "three fields once shipped unreadable".
- **Risk:** a different order of `writeField` calls or of SET fields (diagnostics change); a changed default.
  Mitigation: generate, then compare the generated tables with today's by a script (names, order, defaults,
  comments).
- **Verification:** A2 round-trip and golden tests; every corpus IPM file written in full and in brief mode
  before/after is byte-identical; `sizeof(BASE_PARAM)` and every `offsetof` identical (static_assert); freeze.
- **Size:** M (mostly generated). **Depends on:** A2. **Other repos:** GEMSGUI optional later PR: generate
  `ResetTail`, the `o_pax*` enum entries and `vis_od.ini.json` rows from the same `.def` (08 Q6).

### A4 — Settings access
- **Moves:** add `const BASE_PARAM& TMultiBase::settings() const` (calls `base_param()` — the GEMSGUI override
  stays the source); replace internal reads `base_param()->X` by `settings().X` (mechanical, ~1 000 sites); writes
  only through `LocalSettings`. A clang-tidy check or a grep in CI rejects new non-const uses inside GEMS3K.
- **Why:** H13: today `base_param() const` returns a writable pointer shared by node copies.
- **Risk:** low; a hidden write found by the compiler becomes visible (each one is listed in the PR).
- **Verification:** compiles; freeze.
- **Size:** L (mechanical). **Depends on:** A3. **Other repos:** none.

### A5 — Getters that change state; stale outputs
- **Moves:** add pure functions `gibbsEnergyOf(const AmountsView&)` (same formula and order as
  `TotalGibbsEnergy()`, but writing scratch into a local buffer); keep `TotalGibbsEnergy()`,
  `Get_GibbsEnergy()`, `Get_PMc()` unchanged (they forward and still write X/XF/XFA where they do today); move
  `ReportAbsentPhaseIndexes()` out of the `try{}catch(...)` inside `packDataBr` into a named
  `refreshAbsentPhaseIndexes()` called at the same point; document trap 6 (stale DATABR on ERR) in `node.h` with an
  ASD-STE100 note:
  *"If the status is ERR_GEM_*, the values in DATABR are from the previous call. Do not use them."*
- **Why:** H10, H14: make side effects visible without changing them.
- **Risk:** the pure variant must not be used by the solver in Stage A (it would skip the X/XF/XFA writes). Rule:
  only new callers use it.
- **Verification:** freeze; golden trace.
- **Size:** S. **Depends on:** A2. **Other repos:** xgems/GEMSGUI may switch to the pure getter later (optional).

### A6 — `RunContext` replaces the thread-local flags (five sub-steps)
- **Moves:** introduce `RunContext` (04 §3.1), created in `GEM_run` and passed by reference through the layers and
  into `CalculateEquilibriumState*` (new overloads; old signatures forward with a context they create).
  - **A6a** re-entry guards: `g_optimaCapCall`, `g_optimaDryRedoCall`.
  - **A6b** budgets: `g_optimaIterBudget`, `g_optimaCapAllowance` (incl. the Optima convergence hook, the
    extinction tier refill, and the `Budget` RAII in NativeOptimaFinish/Rescue — keep its "reset to −1" behaviour,
    08 R9).
  - **A6c** attempt switches: `g_optimaFixAccept`, `g_optimaQuickTry`, `g_optimaLSEscapeOff`,
    `g_optimaLineSearchRetry`, `g_optimaEscRetryNext`, `g_optimaCgSeedArmed`, `g_ipmLevelOff`, `g_optimaDryOff`,
    `g_tpdStrongSearch`.
  - **A6d** attempt memory and notes: `g_optimaPreSolveEmpty`, `g_lvlFullCalls`, `g_optimaTargetWarned`,
    `g_optimaDryFailed`, `g_ipmNoProgArmed/Fired`, `g_ipmLevelUsed`, `g_aqueousAbsent`, `g_finishFromSuccess`;
    `BeginCalculation()` → `RunContext::beginAttempt()` (same per-attempt scope).
  - **A6e** helper outputs → return structs: `ColumnGenerationSeed` (`g_cgConverged`, `g_cgPseudo`, `g_cgGap`,
    `g_cgDual`), Optima seeds (`g_coldLvlUsed`, `g_lvlShortcut`, `g_lvlCgStart`, `g_ps0GaveUp`).
  Member flags used as parameters (`optima_hop_leg`, `optima_disable_smoothing` guard, `referenceMode`,
  `runKinetics`) → `OptimaLeg` argument (H11), with the member kept as a mirror where code outside the call
  chain reads it (smoothing).
- **Why:** H1, H2, H11; brief §3.7 ("make each a named layer calling the next directly" — done — but the flags
  remain); C++ Core Guidelines I.2, CP.3, F.21; PHREEQC→IPhreeqc precedent.
- **Risk:** a flag is read on a path not in the call chain (e.g. through `TNodeArray`, a GEMSGUI override, or a
  second `GEM_run` from inside a layer). Mitigation: per sub-step, `grep` list of every read/write in the PR; the
  DECIDE golden test catches a changed decision; keep the `thread_local` as a debug-build mirror that asserts it
  equals the context value for one release (removed in A14).
- **Verification:** freeze per sub-step (`# dec` lines are the sensitive part); golden trace; TSan build of the
  thread test (`tools/thread_test.cpp`) — expected: fewer reports than before, none new.
- **Size:** 5 × M. **Depends on:** A2 (A4 helps). **Other repos:** none (public signatures unchanged).

### A7 — Solver workspace
- **Moves:** `SolverWorkspace` member of `TMultiBase` holding the scratch now in function-static `thread_local`s
  (`SolveIpmAugmented` vectors `act, dreg, Q, b, Rd, v` and Eigen matrices, ipm_main), the compressed A index
  (`arrL/arrAN`), simplex `new[]` buffers, the `Dscale` vector; same sizes, same fill, same order, still rebuilt
  where it is rebuilt today.
- **Why:** H3; prepares B allocation work and node batching.
- **Risk:** a scratch buffer that today carries values between calls (e.g. a `thread_local` vector that is not
  cleared) — moving it per object changes the carry-over between nodes. Mitigation: for each buffer, prove "written
  before read" in the PR, or keep it per thread.
- **Verification:** freeze; instruction count (A-rule 3).
- **Size:** M. **Depends on:** A6. **Other repos:** none.

### A8 — Process-wide statics and environment reads
- **Moves:** `native_cert_dualfree_resolved` (ipm_main l.316, process-wide `static int`, data race) → resolved
  once per `RunConfig` with the same value semantics (read the resolver first: if it depends on the first
  system's settings, keep "first call in the process" semantics and make it `std::call_once`, 08 Q9);
  warn-once sets (ipm_main l.2672, l.2775; node.cpp `bigWarned`, `warned`; ipm_optima `warnedPairs`) → a
  `WarnOnce` registry in `DiagnosticsSink` with the same scope (process or thread) as today; `GEMS3K_FILLOUT_BUDGET`
  (read per `FilloutBudgetValue` call) → read once per call into `RunConfig`, printed in EFF;
  `GEMS3K_OPTIMA_NET_RESUME` unchanged (already once per process, already in EFF); `TNodeArray::na` kept for the
  API but `RunGEM(Mode, nNodes, …)` uses `this`; Optima fork: `trace::state()` → `thread_local` (separate PR in
  ../optima, owner's decision, 08 Q8).
- **Why:** H4, H6; CP.2.
- **Risk:** a warning printed once per thread instead of once per process (logs only). Changing
  `native_cert_dualfree_resolved` semantics could change CERT values in a multi-project process.
- **Verification:** freeze (one project per process, as the freeze runs); golden trace; TSan.
- **Size:** M. **Depends on:** A6. **Other repos:** optima (one PR, optional).

### A9 — Activity-model statics into the instance
- **Moves:** `TCGFcalc` memo statics (s_solmod2 l.1422, l.1483, l.1525, l.1578, l.1758) → members of the
  `TCGFcalc` instance; `TELVIS::trapzd` `static double s` (s_solmod4 l.4678) → member; read-only tables
  (s_solmod5 l.1485/1493) → `static const`.
- **Why:** H5: wrong-instance cache and a data race (09 §3 fact 3).
- **Risk:** **this can change results** where two CG-fluid phases (or two threads) shared a cache hit. In
  single-thread freezes a hit needs two instances with equal T and ρ (unlikely). Rule: if the freeze shows any
  difference, the step moves to Stage B as a scored bug fix (B2) and Stage A keeps the statics.
- **Verification:** freeze; a new unit test with two CG phases at the same (T, ρ).
- **Size:** S. **Depends on:** A2. **Other repos:** tsolmod4rkt (same files; no API change).

### A10 — Files into module folders; split the large functions
- **A10a** `git mv` files into the folders of 04 §2 (no content change except include paths; installed header names
  unchanged). Size L (moves). Verification: identical object code per function (`objdump` per symbol) except
  debug info; freeze.
- **A10b** split `ipm_main.cpp` (5 814 lines) by topic: `kernel/ipm_loop.cpp` (MBR, InteriorPointsMethod,
  step size, golden section), `kernel/linsolve.cpp` (normal equations, SolveIpmAugmented, Eigen QR, JAMA),
  `kernel/certify.cpp`, `kernel/fillout.cpp`, `diag/native_trace.cpp`. Hot callers and callees stay together
  (A-rule 3).
- **A10c** split `CalculateEquilibriumStateOptima` into the ~15 tier functions of 04 §3.4 (extract function,
  statements unchanged and in the same order; shared locals become an `OptimaRun` struct). The pre-solve and the
  full solve get one shared `makeObjective()` **only if** the two lambda bodies are token-identical after the
  extraction; where they differ (reduced index maps), keep two bodies with a comment.
- **A10d** move the GEM_run layers to `run/` (`pipeline.cpp`, one file per layer).
- **Why:** goal 4 (one job per module); review size; today file names hide content.
- **Risk:** inlining changes (A-rule 3); a lambda capture by reference to a local that now lives in another frame
  (today the objective captures `&fixedGrad`, 04 §5).
- **Verification:** freeze per sub-step; instruction counts; ASan/UBSan build of the 8-case set.
- **Size:** 4 × L. **Depends on:** A6, A7, A13. **Other repos:** none (installed headers keep names);
  gems-benchmark scripts that grep source paths, if any (08 Q1).

### A11 — Warm-start state as one object; node-array seeding made explicit
- **Moves:** `WarmState` (U, Y, prev_scfact, zeroed_floor_amount/warm_floor_seed, prevBulkMol_, prevTC_, prevP_,
  prevDulPositive_, warmInputChanged_, phase_floor_parked) on `TMultiBase`. `TNodeArray` work-node seeding: keep
  today's behaviour (the previous node of the same thread seeds the next) but name it
  (`WarmSeedPolicy::PreviousNodeOfThisWorker`) and add a test that reproduces it.
- **Why:** H9, H11; 09 §3 fact 2. Makes the B decision (seed from the node's own previous state) a one-line change.
- **Risk:** low if fields keep their write points.
- **Verification:** freeze; a TNodeArray test (serial) comparing per-node results before/after.
- **Size:** M. **Depends on:** A6. **Other repos:** GEMSGUI review (`TNodeArrayGUI` overrides `CalcIPM_Node`).

### A12 — De-duplication without numerical change
- **Moves:** one `DatabrSnapshot` RAII helper for the 12 copies of the snapshot pattern (node.cpp); one
  `mapStatus(need, outcome)` for the three NEED→OK/BAD/ERR ladders (node.cpp l.1597–1662); one body for the six
  `GEM_from_MT` variants (node.cpp) and the duplicate `GEM_from_MT_time`/`GEM_set_MT`; one table per character code
  only where tables are identical; dense kernels only where bit-identical (compare `SymEigFloorInPlace` with
  `CertSymEigMin`, the two dense LU/MGS copies; merge identical ones, keep the others with a cross-reference).
  `TActivity`/`TKinetics` duplicates of TMultiBase code: mark `[[deprecated]]` (02 §6), no merge in A.
- **Why:** goal 3.
- **Risk:** two "identical" kernels that differ in one loop order (FP). Rule: merge only after a token diff shows
  identity.
- **Verification:** freeze; golden trace.
- **Size:** M. **Depends on:** A10. **Other repos:** none.

### A13 — Diagnostics component
- **Moves:** trace writers (ipm_main RUN/SET/EFF/KEY/CERT, node_trace, DECIDE writer) and probes to `diag/`;
  `DiagnosticsSink` owns the `FILE*` handles (today one function-static per probe, two for `GEMS3K_PHSTAB_PROBE`)
  — same open time (first use), same append mode; `native_trace_decide()` forwards. Logging levels as brief
  §10.4 (default warn, decisions info, numbers debug) only where already so; a level change is a separate PR with
  the owner's approval.
- **Why:** goal 5; H3/H4 for probe handles.
- **Risk:** output order changes when two handles of one file become one (the double PHSTAB handle): the probe
  file content order can change. Accept (probe, not freeze) or keep two handles (08 Q10).
- **Verification:** golden trace byte for byte; probe files byte for byte on the golden inputs.
- **Size:** M. **Depends on:** A2. **Other repos:** gems-benchmark (reads traces; unchanged).

### A14 — Tests, documentation, clean-up
- **Moves:** unit tests for resolvers (`optima_dimreduce_passes` etc.), status mapping, snapshot helper,
  name maps (trap 8 first-match behaviour), scaling round trip (trap 3: scale → rescale restores every scaled
  array; lists the not-scaled ones); module `README.md` per folder in ASD-STE100 style; remove the debug mirrors
  of A6; deprecation notes for internal headers (02 §5); update `Docs/` (GLOSSARY keeps the fixed terms).
- **Why:** BSSw/Duck Book: tested, documented modules.
- **Risk:** none for results.
- **Verification:** CTest green; freeze.
- **Size:** M (continuous). **Depends on:** each step adds its tests. **Other repos:** none.

## 3. Hidden-state removal map (00 §4 → steps)

| Mechanism | Step | Removed in A? |
|---|---|---|
| H1 layer thread-locals | A6a–d | yes |
| H2 thread-local helper outputs | A6e | yes |
| H3 function-static thread-local caches/scratch | A7 (scratch), A8 (warn-once) | yes |
| H4 process-wide statics | A8 (+ optima PR) | yes, except warn-once registry (documented exception) and `TNodeArray::na` (API) |
| H5 model memo statics | A9 | yes, or moved to B2 if the freeze changes |
| H6 env reads that change results | A8 | moved into RunConfig and EFF; names kept |
| H7 TSolMod aliases MULTI | A10 + B (PhaseModel blocks) | wrapped in A; removed in B |
| H8 non-idempotent activity call | — | documented in A; changed only in B (it is FP order) |
| H9 MULTI shared across nodes | A11 | made explicit; behaviour kept |
| H10 mutating getters | A5 | pure variants added; old kept |
| H11 member flags as parameters | A6 (OptimaLeg), A11 (WarmState) | yes |
| H12 layout copies by count | A2 (documented), A3 (BASE_PARAM) | MULTI copy fixed in B1 |
| H13 shared writable settings | A4 | yes (const access) |
| H14 stale outputs on failure | A5 (documented) | behaviour kept (contract, 02 §3) |
