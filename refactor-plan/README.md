GEMS3K commit: 54b5f819ba83589daa9cd9aa0a9ff699adb28b1d (branch develop_optima, 2026-11-03)

# GEMS3K restructuring plan — summary

**What this is.** A plan to make GEMS3K easier to change and then faster, in two strictly separated stages. No
code was changed, committed or pushed. The detail is in files 00–09; `NOTES.md` holds the reading notes.

## Main findings

1. **The code works, but much of its state is hidden.** About 25 thread-local flags carry decisions between the
   retry layers of `GEM_run`. Some activity models keep process-wide caches keyed only on T and density, so two
   phases or two threads can get each other's values. Two environment variables change results.
2. **Some contracts exist only by discipline.** A new `BASE_PARAM` field needs 8 edits in GEMS3K
   (the brief's 7 sites, with the count 164 typed twice by hand) and up to 5 more in GEMSGUI. `copyMULTIData`
   copies fields by memory layout. Its counts are already off by one, so `FitVar[1..4]` are not copied.
3. **There are real bugs to fix later.** `TNodeArray` counts Optima "OK" codes as errors. Its warm-mode results
   depend on the thread count. dn/db is taken from the first Optima attempt, not from the accepted answer. The IPM
   reader has an out-of-bounds loop. These fixes change answers, so they belong to Stage B, after an owner decision.
4. **Where the time goes** (21 public projects, measured). Small systems take 0.1–3 ms per call in any mode.
   Cold Optima runs 25–35× slower than the GEMS3K IPM solver on mid-size systems. Failing warm calls cost
   1 000–13 000× a normal call. In Optima runs, a quarter of the work is bookkeeping around the activity models
   (all phases are recomputed for every finite-difference column). The model formulas themselves are cheap.
5. **All 14 traps of the brief were found in the code.** Each has its file, function and line (09).

## Target structure (04)

Eleven modules in folders and CMake targets: core, data, io, thermo, models, numerics, kernel (GEMS3K IPM solver),
optima (adapter), run (retry layers), node, nodearray, plus diag (traces). The installed library and headers
keep their names. A `RunContext` object replaces the thread-local flags. `BASE_PARAM` gets one declaration
file. Activity models move behind a `PhaseModel` interface. The diagnostics stay byte for byte.

## Order of work

- **Stage A (05, behaviour kept, 25–35 small pull requests).** Steps in order:
  - A0 baseline freeze on the owner's machine;
  - A1 build cleanup;
  - A2 contract tests;
  - A3 `BASE_PARAM` from one declaration;
  - A4–A5 settings access and getters;
  - A6 `RunContext`;
  - A7–A9 workspace, statics and environment variables;
  - A10 files into modules and splitting the 4 400-line Optima function;
  - A11 warm-start state;
  - A12 de-duplication;
  - A13 diagnostics component;
  - A14 tests and docs.

  Every step must reproduce the STD and RAW freezes row for row, including `# dec`.
- **Stage B (06, scored like a solver change).** Steps in order:
  - B0 nudge baseline;
  - B1 bug fixes;
  - B2 thread safety;
  - B3 failure cost;
  - B4 activity-call scope;
  - B5 analytic Hessian blocks;
  - B6 pre-solve;
  - B7 linear algebra;
  - B8 line search;
  - B9 allocations;
  - B10–B12 floating-point policy, parallel node batches and libm.
- **Extensions (07).** dn/db at the accepted answer comes first. No GPU work starts until a CPU batch test shows
  that most nodes stay on the common path.

## Five biggest risks (08)

1. **The freeze runs only on the owner's machine.** gems-benchmark became accessible: CTest passes 92 of 96, and
   the CI check reproduces 114 of 132 recorded rows (the rest differ because of the machine, not the code). But
   the freeze driver is in no repository, so the owner must run the baseline (A0) or publish the driver.
2. **Moving hot code can change the floating-point results** through inlining. Each step checks the freeze and
   the instruction counts.
3. **A hidden flag may be read on an unexpected path** while the thread-locals are replaced. This is handled by
   per-flag lists, debug mirrors and the DECIDE golden test.
4. **GEMSGUI depends on internals:** it subclasses `TMultiBase`, reads `MULTI` and stores `BASE_PARAM` by
   position. Those parts stay fixed, with compile-time checks.
5. **The activity call is not repeatable** (it keeps a history), and freezes differ between CPU types. So Stage A
   never changes the number of calls, and all freezes use one machine class.

Open questions with assumed defaults: 08 §3 (23 items).
