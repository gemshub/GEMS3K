GEMS3K commit: 54b5f819ba83589daa9cd9aa0a9ff699adb28b1d (branch develop_optima, 2026-11-03)

# 03 — Profile

## 0. Limits of this profile (read first)

- **Not the owner's corpus.** gems-benchmark was not accessible when this profile was made (it was granted later
  in the session; the checks then run are in §7). The freeze itself still cannot run here (§7). The profile uses a substitute corpus of
  **21 projects** found in the public repositories (GEMS3K `node-gem-build`, `tools-build`, `tools`,
  `tsolmod4rkt-build`; xgems `demos/resources`, `tests/gems3k`). Larger systems from the owner's corpus may move
  the balance.
- **The committed freeze was not reproduced** (the freeze driver and the recorded freezes are on the owner's
  machine only, §7). Every number in §1–§3 is a measurement of this build on this machine.
- Machine: 4 cores, 15 GB, Linux 6.18 container. No `perf`. Deterministic instruction counts from
  `valgrind --tool=callgrind`; wall time from `std::chrono::steady_clock` and from the existing
  `GEMS3K_TIMING_FILE` per-part timer (`gems3k_timing.h`).
- Build: conda gcc 16.2, Eigen 3.4.0, spdlog 1.17, fmt 12, ThermoFun 0.7.3, Optima fork `optima_develop` 145aeed,
  GEMS3K `USE_OPTIMA_SOLVER=ON`, `USE_THERMOFUN=ON`, Release (`-O3 -DNDEBUG`), `BUILD_SOLMOD=OFF`.
- Driver: `prof_driver` (source in §6): `GEM_init`, then N calls of `GEM_run`. Warm modes first make one cold
  call (seed), then sweep `bIC` by ×(1 + 1e-3·sin(r+1)) like a transport step.

## 1. Mode survey (wall time, 5 calls per mode, every project)

All 126 runs ended without a crash. Per call (ms), selected projects (full table: `NOTES.md` "Survey complete"):

| Project | AIA | SIA | AOP | SOP | HOP | SHP | Remark |
|---|---|---|---|---|---|---|---|
| gems_files | 0.24 | 0.09 | 0.69 | 0.30 | 0.63 | 0.26 | smallest |
| CalciteIC | 1.66 | 1.44 | 2.49 | 0.86 | 1.64 | 1.03 | typical small |
| sw_test | 0.66 | 0.48 | 1.37 | 1.50 | 1.35 | 1.13 | typical small |
| CASHNK1 | 1.79 | 0.46 | **62.4** | 1.30 | 2.87 | 2.88 | AOP 35× AIA |
| CemGEMS-keyvalue | 22.1 | 4.67 | **555** | 23.7 | 51.1 | 41.3 | AOP 25× AIA |
| mines23 | 591 | 424 (3/5 ok) | **1848** | 193 | 550 | 309 | largest |
| lb_test | 2.72 | **4596 (1/5 ok)** | 4.78 | **35052 (3/5 ok)** | 3.67 | **32180 (3/5 ok)** | escalation |

Findings:

1. Small systems cost 0.1–3 ms per call in every mode. The modes differ by at most ~3×.
2. **Cold Optima (AOP) is 25–35× the GEMS3K IPM solver on mid-size systems** (CemGEMS, CASHNK1) and 3× on mines23.
3. **Escalation dominates the tail.** lb_test warm calls that fail go through the whole rescue ladder:
   4.6–35 s per call versus 2.7 ms for a normal call (×1 000–13 000). Mean cost on such a corpus is set by the
   failures, not by the normal path.
4. Native modes can fall into Optima through the finish/rescue layers: Thermo-time-in and test01 AIA spend
   80–83 % in Optima (7.9–14.4 ms, versus 1.7 ms SIA); mines23 AIA spends 72 % in Optima.

## 2. Where the time goes (GEMS3K_TIMING parts, wall time, share of the call)

From 693 `TIMING` lines (all projects × modes, `prof/timing_summary.txt`):

- **Activity-coefficient calls** (`act`, inclusive `CalculateActivityCoefficients`) take **15–63 % of call time in
  every mode and project**: native AIA/SIA 16–46 %, Optima 14–63 %. Cost per call ≈ 2–6 µs (small systems),
  ≈ 18 µs (CemGEMS).
- Optima objective + Hessian (`obj`) 13–74 %. In escalated calls the finite-difference Hessian columns (`fdmix`)
  dominate: lb_test SIA/SOP/SHP fdmix 59–62 %, act 60–63 %, **0.46–3.7 million activity calls per call**.
- **Dimension-reduction pre-solve** (`presolve`) dominates large cold Optima calls: CemGEMS AOP 85 %
  (16 409 activity calls per call), SOP 69–70 %, mines23 AOP 57 %.
- **PotentialSpaceFinish + TPD search** dominate CASHNK1 AOP (psf 89 %, tpd 79 %).
- Wall time of five identical cold calls varies by 5–10 % (CemGEMS AOP 0.528 / 0.558 / 0.585 s). Use instruction
  counts for decisions; for wall time take ≥ 5 repeats and the median.

## 3. Where the instructions go (callgrind, collected inside `TNode::GEM_run` only)

Self instruction counts by category, % of all instructions inside `GEM_run` (`--toggle-collect=TNode::GEM_run(bool)`):

| Run (calls) | M Ir | Linear algebra | GX + phase sums | log/exp/pow | Concentrations + γ conversion | Chem. potentials | Model formulas | Optima own | Simplex | Allocation | Other |
|---|---|---|---|---|---|---|---|---|---|---|---|
| sw_test AIA (5) | 31 | 18.2 (hand QR) | 15.4 | 19.2 | 11.5 | 6.3 | 0.8 | — | 11.6 | 3.2 | 13.0 |
| CalciteIC SIA (5) | 50 | 23.1 (hand QR) | 19.4 | 22.1 | 12.5 | 6.5 | 0.9 | — | 0.2 | 1.2 | 13.4 |
| CemGEMS AIA (1) | 202 | 53.2 (Eigen QR 41.8 + hand QR 11.4) | 11.0 | 11.6 | 6.5 | 4.4 | 0.8 | — | 2.0 | 1.5 | 8.2 |
| CASHNK1 AIA (5) | 87 | 11.0 | 16.5 | 20.9 | 8.1 | 4.6 | **21.9** (Berman) | — | 0.4 | 3.7 | 12.3 |
| sw_test AOP (5) | 46 | 16.1 (FullPivLU 15.0) | 0.8 | 11.5 | **25.0** | 9.9 | 1.6 | 15.7 | — | 4.5 | 13.7 |
| CalciteIC SOP (5) | 28 | 10.5 (FullPivLU 9.5) | 0.8 | 12.0 | **26.2** | 9.8 | 1.6 | 15.1 | — | 5.8 | 16.6 |
| CemGEMS AOP (1) | 5 709 | 15.6 (FullPivLU 15.2) | 1.1 | 15.5 | **26.9** | 13.7 | 2.8 | 7.3 | — | 0.8 (+ memset 4.8) | 11.5 |
| CASHNK1 HOP (5) | 136 | 8.6 | 10.7 | 16.6 | 10.1 | 5.0 | **26.7** (Berman) | 2.3 | 0.2 | 5.1 | 14.0 |

Attribution to the five areas the task names:

1. **Activity models**: the model formulas are cheap (0.8–2.8 % self) except the multisite Berman model
   (22–27 %). The expensive part is the **work around the models**: on every activity call
   `CalculateConcentrations`/`…InPhase`, `PhaseSpecificGamma`, `phase_bcs` and `IS_EtaCalc` run for **all**
   phases (6.5–27 %), plus the logarithms. An FD Hessian column changes one species of one phase but pays for all
   phases. With chemical potentials and log/exp/pow, the "activity" share reaches 30–55 % of instructions in
   Optima runs.
2. **Linear algebra**: native — the hand-written Householder QR (`SolveIpmAugmented`, 11–23 %) and, from
   La ≥ 100 (`pa_IpmQrSortPivot`), Eigen `ColPivHouseholderQR` (CemGEMS 42 %); Optima — Eigen `FullPivLU` inside
   Optima, 9.5–15 %. The JAMA paths are cold (fallback, MBR).
3. **Phase selection**: PSSC + StabilityIndexes < 0.5 % in native runs. In Optima runs the phase-selection
   *machinery* (pre-solve, PotentialSpaceFinish, TPD) is what costs (§2), through repeated solves and activity
   calls, not through its own arithmetic.
4. **I/O and logging inside `GEM_run`**: printf/strtod 0.0 %; `native_trace*` 0.0–0.1 % with tracing off — the
   zero-cost claim of brief §3.14 holds. (Outside `GEM_run`, static initialisation of loggers, iostreams and
   ThermoFun is 30 % of a short CalciteIC program.)
5. **Allocation**: malloc/free/new 0.8–5.8 % self (sources: `Build_compressed_xAN` every solve, DATABR snapshots
   in `GEM_run_tries`, `new[]` inside simplex NEW/WORK, TNT `Array2D` per normal-equation solve, `std::vector`
   temporaries in the Optima adapter). Plus `memset` 4.8 % in CemGEMS AOP: `res.fxx.setZero()` of a dense L×L
   Hessian on every Hessian evaluation.

Also: `GX` with `TotalPhasesAmounts` (11–19 %) is called ≥ 2 times per golden-section step on the full vector.

## 4. What this means for the plan

- Stage A (behaviour-preserving) will not make GEMS3K faster, and must not try. The profile shows that no Stage A
  step touches a hot loop except by moving it; moving code must keep inlining (same translation unit or LTO off/on
  unchanged), so A steps check instruction counts as well as the freeze (05 §0).
- Stage B priorities by measured share (06): (1) escalation cost (rescue ladder, FD columns: up to ×13 000);
  (2) work around the activity models per call (25–27 % in Optima runs) and the FD Hessian (fdmix 59–62 % in
  escalations); (3) pre-solve in large cold Optima calls (57–85 %); (4) linear algebra (10–53 %); (5) GX/golden
  section (11–19 %); (6) allocation (1–6 %).
- Floating-point note found while profiling: glibc chooses `__ieee754_log_fma` at run time from the CPU
  features. A freeze taken on a CPU without FMA can differ from one taken with FMA even with the same binary.
  Freeze pairs must come from one machine class (08 R6).

## 5. Not measured (and why)

- The owner's 72–87-project corpus (not accessible).
- Cache misses and branch mispredictions (`perf` not available; `cachegrind` possible, not run: the hot loops are
  small and the instruction mix already explains the cost).
- Multi-threaded scaling of `TNodeArray` (no OpenMP build was made; the thread-count dependence of warm modes in
  09 §3 fact 2 follows from the code, not from a measurement).
- GPU: nothing to measure yet (07).

## 6. How to repeat (scripts used; outside the repository)

Driver (`prof_driver.cpp`, built against the installed GEMS3K with `-O2 -g -DUSE_OPTIMA_SOLVER -DUSE_SPDLOG_PRECOMPILED`):

```cpp
// Usage: prof_driver <path-dat.lst> <AIA|SIA|AOP|SOP|HOP|SHP> <nrep> [sweep_rel]
#include "node.h"
#include <chrono>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>
static long modeCode(const std::string& m)
{
    if (m == "AIA") return NEED_GEM_AIA;  if (m == "SIA") return NEED_GEM_SIA;
    if (m == "AOP") return NEED_GEM_AOP;  if (m == "SOP") return NEED_GEM_SOP;
    if (m == "HOP") return NEED_GEM_HOP;  if (m == "SHP") return NEED_GEM_SHP;
    return -1;
}
int main(int argc, char* argv[])
{
    if (argc < 4) return 2;
    const std::string lst = argv[1], mname = argv[2];
    const long mode = modeCode(mname);
    const int nrep = std::atoi(argv[3]);
    const double sweep = (argc > 4) ? std::atof(argv[4]) : 1e-3;
    TNode node;
    if (node.GEM_init(lst.c_str())) return 1;
    DATACH* dch = node.pCSD();
    DATABR* dbr = node.pCNode();
    std::vector<double> b0(dbr->bIC, dbr->bIC + dch->nICb);
    const bool warm = (mode == NEED_GEM_SIA || mode == NEED_GEM_SOP || mode == NEED_GEM_SHP);
    const long coldMode = (mode == NEED_GEM_SIA) ? NEED_GEM_AIA : (mode == NEED_GEM_SOP) ? NEED_GEM_AOP
                         : (mode == NEED_GEM_SHP) ? NEED_GEM_HOP : mode;
    if (warm) { dbr->NodeStatusCH = coldMode; node.GEM_run(false); }
    for (int r = 0; r < nrep; r++) {
        const double f = warm ? (1. + sweep * std::sin((double)(r + 1))) : 1.;
        for (long i = 0; i < dch->nICb; i++) dbr->bIC[i] = b0[(size_t)i] * f;
        dbr->NodeStatusCH = mode;
        auto a = std::chrono::steady_clock::now();
        long st = node.GEM_run(false);
        auto b = std::chrono::steady_clock::now();
        std::printf("CALL r=%d mode=%s status=%ld it=%ld sec=%.6f G=%.12e\n", r, mname.c_str(), st,
                    dbr->IterDone, std::chrono::duration<double>(b - a).count(), node.Get_GibbsEnergy());
    }
    return 0;
}
```

Commands: `GEMS3K_TIMING_FILE=<out>.timing prof_driver <lst> <MODE> 5` per project and mode;
`valgrind --tool=callgrind --toggle-collect='TNode::GEM_run(bool)' prof_driver <lst> <MODE> 5`;
`callgrind_annotate --auto=no --threshold=100 <file>`. Note: `Get_GibbsEnergy()` mutates X/XF/XFA (09); the
driver calls it after the timed region only.

## 7. gems-benchmark checks (run after access was granted, 2026-11-03)

gems-benchmark `develop` dc220de (the default branch `main` 5651c52 has no freeze tools). Built against the same
GEMS3K install as §0 (Optima on, ThermoFun on). Results:

| Check | Result here | Notes |
|---|---|---|
| `tools/trace_selfcheck.sh` | **OK** — RUN, BULK, SET, EFF, KEY, DECIDE all emitted | the install writes every record the freeze reads |
| CTest, 96 registered tests (GUIDE says 42; feature-watch cases added) | **92 pass, 4 fail** in 91 s | `optima_regression.aop.restart` (1 check: `earlynet_f_CalcDolo`, armed = unarmed = 18 it), `proposed.aop` (ratchets `hoknob_T2b` penalty 1.0× < 20×, `crossing_T11`), `dimreduce.gate` (AUTO on the wrong side of the gate for T8_aq101/T8_aq201), `ci.baseline` (not valid with an Optima build) — 08 Q22 |
| `tools/ci.sh` (stock GEMS3K without Optima, 132 recorded native rows) at head 54b5f81 | **28 findings** (26 answer, 2 cost) on 9 projects; 114 rows reproduce exactly | baseline recorded 2026-10-20 with gcc 14.3 |
| same, built with gcc 14.3 instead of 16.2 | **the same 28 findings, identical to the last digit** | the compiler version does not move answers here |
| same, at GEMS3K cc60bb7 (the code the baseline was recorded from, 2 min later) | **the same 28 findings** | the 26 commits since did not change any CI row here |
| same, with glibc's FMA/AVX dispatch disabled (`GLIBC_TUNABLES=glibc.cpu.hwcaps=-AVX2,-FMA,…`) | **23 findings**: all 5 of `f_GEOTHERM` (incl. ITG 687 → 688) disappear | libm CPU dispatch alone changes answers and iteration counts (R6) |
| `tools/freeze.sh` | **cannot run here** | needs `mode_compare` from the owner's local `GEMS3K/debug-optima-vs-reaktoro/` (in no repository); hard-coded `/home/dmiron/...` paths (freeze.sh l.59–62, build-gems3k.sh, recheck.py l.22–28) |
| `tools/recheck.py --tier fast` | 2 of 51 claims hold — **not meaningful here** | reads the standard freeze from `.claude/docs/freeze/` (git-ignored, owner's machine only) |

Build problems found on the way (A1 and a gems-benchmark fix): a default GEMS3K configure needs pybind11 even with
`BUILD_SOLMOD=OFF` (`BUILD_SOLMOD_PYTHON` defaults ON); every configure rewrites the tracked `gems3k_version.h`;
gems-benchmark `include/difftest/optima_regression.h` l.608/l.643 call `::getpid()` without `<unistd.h>` (fails with
gcc 16.2; worked around here with `-include unistd.h`, no source change).

Consequence for the plan: the remaining 23 findings are environment-dependent (most likely the ThermoFun or glibc
version; 6 of the 9 projects are ThermoFun `f_` projects). Freeze pairs for Stage A must come from **one fixed
environment** (machine class, glibc, conda package versions), recorded in the freeze header (05 A0, 08 R6).
