# GEMS3K 5.0.0 — options guide

What is new in this version, how to choose a solver, how to run it from code, how to change a setting in the
project file, what every setting does, and the known limitations. What runs, in what order, in each mode, and
what is tried when a calculation fails, is in [Solver paths](SOLVER-PATH-GRAPH.md)
([drawn version](SOLVER-PATH-GRAPH.html)).

## The two solvers

GEMS3K computes chemical equilibrium by minimising the Gibbs energy of the system. It has two solvers, and
you choose one for each calculation (see [The calculation modes](#the-calculation-modes)).

**The GEMS3K IPM solver** is the solver GEMS3K has always had: an interior-point method (IPM) followed by a
refinement of the mass balance. It is fast on almost every system and is the default (modes AIA and SIA).
Its settings are the `pa_...` entries of the project's `-ipm` file:
[long-standing settings](#gems3k-ipm-solver--long-standing-settings),
[original settings](#original-settings-not-listed-above) and
[settings added in this version](#gems3k-ipm-solver--settings-added-in-this-version).

**The Optima solver** is a second solver, based on the open-source Optima optimisation library
([github.com/gemshub/optima](https://github.com/gemshub/optima)) with changes made for GEMS3K. It is more
robust for solid solutions, melts and other systems where several mixing phases compete (modes AOP and SOP),
and it can also check the answer of the GEMS3K IPM solver (modes HOP and SHP). It can also fix the pH or Eh
([Fixing the pH or Eh](#fixing-the-ph-or-eh-optima-modes-only)). Its settings are also `pa_...` entries of the
`-ipm` file, with names starting `pa_Optima`: [main settings](#optima-solver--main-settings) and
[supplemental settings](#optima-solver--supplemental-settings). The Optima library has further options of its
own; they are not set in GEMS3K project files, and GEMS3K chooses them. They are described in the Optima
library's own options guide (`docs/OPTIONS-GUIDE.md` in its repository).

Some general settings, such as `pa_PE`, `pa_DHB` and `pa_DG`, apply to both solvers. All settings are in
[Settings](#settings), which also shows how to add one to a project file.

## What is new

**A second equilibrium solver.** Besides the GEMS3K IPM solver (interior-point method with mass-balance
refinement), GEMS3K can now compute equilibrium with the Optima solver. Optima is more robust for solid
solutions, melts and other systems where several mixing phases compete. You choose the solver for each
calculation by the value you put in `NodeStatusCH` before calling `GEM_run()`. The GEMS3K IPM modes remain the
default.

**Fewer failed calculations.** When a calculation does not converge, GEMS3K now retries automatically —
from a slightly changed start, from scratch instead of from the previous result, or by finishing the
calculation on the phases already found — before it reports a failure.

## The calculation modes

Each mode name says how the calculation starts and which solver does the work.

| Value | Mode | Stands for | What it does |
|---|---|---|---|
| 1 | AIA | **A**utomatic **I**nitial **A**pproximation | GEMS3K IPM solver. Builds its own starting point from scratch. |
| 5 | SIA | **S**mart **I**nitial **A**pproximation | GEMS3K IPM solver. Starts from the previous result held in the node. |
| 10 | AOP | **A**utomatic initial approximation with **OP**tima | Optima solver, starting from scratch. |
| 14 | SOP | **S**mart initial approximation with **OP**tima | Optima solver, starting from the previous result. |
| 22 | HOP | **H**ybrid with **OP**tima | GEMS3K IPM solver from scratch, then Optima refines and checks its answer. |
| 26 | SHP | **S**mart **H**ybrid with o**P**tima | As HOP, but the GEMS3K IPM part starts from the previous result. |

In HOP the GEMS3K IPM answer is kept, and the result is flagged, if the Optima part fails, so HOP does not
give a worse answer than the GEMS3K IPM solver alone.

## Which mode to use

| Situation | Recommended mode | Why |
|---|---|---|
| Most systems, single calculations | **AIA** | The GEMS3K IPM solver is the fastest on almost every system it can solve, often by a large margin. |
| Systems with many species (hundreds or more) | **AIA**, or **HOP** if you want a checked answer | Optima alone is much slower on large systems; HOP gets Optima's check for little more than the GEMS3K IPM solver. |
| Strongly non-ideal systems: solid solutions, melts, miscibility gaps, phase diagrams | **AOP**, or **HOP** | Optima finds the correct set of stable phases far more often where several mixing phases compete. |
| The GEMS3K IPM solver fails, or warns that its answer does not satisfy the mass balance | **AOP** or **HOP** | Optima solves most of these cases. |
| Sequences of small changes: sweeps in temperature or composition, titrations, transport steps | **SOP** after a first point in AOP (or **SIA** after AIA) | Starting from the previous result is many times cheaper than starting from scratch. |
| Sequences that cross phase boundaries | **AOP** or **SOP**; avoid **SIA** | SIA can keep the phases of the previous point after they stop being stable. |
| Systems with very little water (for example nearly dry cement) | **AIA** or **HOP** | Optima can fail when water is close to running out; the GEMS3K IPM solver handles these. |

The steps each mode runs, and the retries it makes before reporting a failure, are shown in
[Solver paths](SOLVER-PATH-GRAPH.md).

## Using it from code

The project files are read once, in `GEM_init()`. There is no call that changes a setting afterwards, so
to try another value, edit the `-ipm` file (see "Settings" below) and initialise a new node.

```cpp
#include "node.h"

TNode node;
if (node.GEM_init("Cu-dat.lst"))               // reads the project, including its pa_ settings
    return 1;

// Choose the solver for this call: the mode is the value of NodeStatusCH.
node.pCNode()->NodeStatusCH = NEED_GEM_AOP;    // Optima, cold start
long status = node.GEM_run(false);
if (status != OK_GEM_AOP)
    std::cerr << "not fully trustworthy, status " << status << "\n";
std::cout << "pH " << node.Get_pH() << "  Eh " << node.Get_Eh() << "\n";

// A second, similar calculation: start from the previous answer (warm start).
node.pCNode()->NodeStatusCH = NEED_GEM_SOP;    // set it again before EVERY call
status = node.GEM_run(false);

// Fix the pH and let the solver add the acid or base needed (Optima modes only).
node.Set_pH_target(7.0);
node.pCNode()->NodeStatusCH = NEED_GEM_SOP;
status = node.GEM_run(false);
```

Two things to remember. `GEM_run()` does nothing unless `NodeStatusCH` was set to a `NEED_GEM_...` value
just before it, so set it before every call. And the same node can be used for any mode; for a comparison
of two settings, make two nodes from two copies of the project, one per `-ipm` file.

A complete working example is `tools/optima_test.cpp`.

## Fixing the pH or Eh (Optima modes only)

In the Optima modes AOP and SOP you can ask for a given pH, a given Eh (in volts), or both. The solver adds
whatever acid, base or electron donor the target needs, and reports how much it added. This is not available
in AIA, SIA, HOP or SHP, and only in a build with Optima.

```cpp
node.Clear_ControlConditions();           // targets stay set until cleared
node.Set_pH_target(7.0);                  // pH 7
node.Set_Eh_target(0.2);                  // Eh 0.2 V (optional; may be used alone)

node.pCNode()->NodeStatusCH = NEED_GEM_AOP;   // or NEED_GEM_SOP; set before every call
long status = node.GEM_run(false);

if (status == OK_GEM_AOP)
    std::cout << "pH " << node.Get_pH() << "  Eh " << node.Get_Eh()
              << "  H+ removed (mol) " << node.Get_ControlCondition_titrant("pH")
              << "  electrons removed (mol) " << node.Get_ControlCondition_titrant("Eh") << "\n";
```

- The tolerance is an optional second argument, `Set_pH_target(7.0, 1e-3)`. Without it, it is derived from
  `pa_GAS`.
- The amount added is part of the answer: the bulk composition returned by the node already contains it.
  To start the next point from the original composition, write the original bulk composition back first.
- Targets stay active for every later call until `Clear_ControlConditions()`.
- A target may not be reachable (for example Eh in a system with no redox couple). Check
  the status, and do not rely on the Eh in such a system.

## Settings

Settings are read from the project's `-ipm` file (`pa_...` entries). Settings that are not in the file
keep their default. The table lists the settings that are new or changed in this version.

**Adding a setting that is not in the file.** Files exported by earlier versions do not contain the new
settings. To use a value other than the default, add the line yourself, next to the other `pa_` entries.

In a JSON file (`…-ipm.json`), add it inside the `"ipm"` block, for example right after `"pa_PLLG"`:

```json
    "pa_PLLG": 10000,
    "pa_OptimaTol": 1e-09,
    "pa_OptimaMaxSeconds": 60,
    "tMin": 0,
```

In a text file (`…-ipm.dat`), add it after the `<END_DIM>` line, for example right after `<pa_PLLG>`:

```
<pa_PLLG>  30000
<pa_OptimaTol>  1e-09
```

Keep each setting only once in the file: if a name appears twice, the later value is used. A setting
placed before `<END_DIM>` in a text file is not read.

### Optima solver — main settings

These are the settings most users may want to adjust.

| Setting | Default | What it does | When it can be useful | Other values |
|---|---|---|---|---|
| `pa_OptimaTol` | 1e-8 | How precisely Optima must converge. | Raise it slightly to trade accuracy for speed on large systems; lower it only if an answer looks insufficiently converged. | Any positive number; smaller is stricter. |
| `pa_OptimaMaxSeconds` | 0 (off) | A time limit in seconds for one calculation, including retries. Results then depend on the speed of the computer. | Long batch or transport runs where one stuck calculation must not hold up the rest. | Any positive number of seconds. |
| `pa_OptimaDimReduce` | 0 (automatic) | For large systems, first solves with only the species that matter, then adds the rest. Automatic means on from 200 species. | Switch it on for a system just under 200 species that is slow; switch it off if a large system gives a different answer with it. | A positive number = on, with at most that many passes; a negative number = off. |
| `pa_OptimaColdRetry` | 2 | If a calculation started from the previous result (SOP, SHP) fails, start again from scratch. 2 gives up on the first attempt quickly. | Use 1 if the warm start usually succeeds but needs many iterations; 0 only when testing. | 0 = no retry; 1 = retry only after the first attempt used its full budget. |
| `pa_OptimaFinish` | 1 (on) | If Optima stops just short of converging, finish the job using the phases it has found. | Leave on. Turn off only to see what Optima alone returns. | 0 = off. |
| `pa_OptimaTpdAccept` | 1e-6 | Accepts a not-quite-converged answer when a direct check shows that no missing phase should be present (within this tolerance). | Leave on. Systems where Optima gets very close but does not meet its own tolerance. | 0 = off. |
| `pa_OptimaAcceptRepair` | 0 (off) | With the setting above, also corrects a small error in the element totals before accepting. | Systems with very little water, or other cases where Optima ends close to the answer with slightly unbalanced element totals. | 1 = on. |
| `pa_OptimaStallWindow` | 500 | Stops a calculation that has not improved for this many iterations, so a retry can start sooner. | Lower it to save time on systems that often get stuck; raise it if slow but steady systems are stopped too early. | 0 = never stop early. |
| `pa_OptimaZeroAbsent` | 2 | How absent species are reported: as a tiny amount, with the element totals corrected if needed. | Use 1 if you need absent phases to read exactly zero in the output. | 0 = tiny amounts as they are; 1 = absent phases set to exactly zero. |
| `pa_OptimaLineSearch` | 1.5 | When a step makes the result much worse (by more than this factor), try a shorter step. | Leave on. A negative value can save time on easy systems that rarely need it. | 0 = off; a negative value = use it only when retrying a failed calculation. |

### Optima solver — supplemental settings

Fine-tuning of the solver. The defaults suit almost all systems; change these only when
investigating a difficult case.

| Setting | Default | What it does | When it can be useful | Other values |
|---|---|---|---|---|
| `pa_OptimaLSStallEscape` | 10 | If shorter steps keep getting nowhere this many times in a row, take one full step to break out. | Systems where Optima gets stuck at the same error for many iterations. | 0 = off. |
| `pa_OptimaLSRejectWorse` | 0 (off) | If a shorter step did not help, use the full step instead. | A few difficult systems where the shorter steps keep making tiny progress; it can slow others. | 1 = on. |
| `pa_OptimaCgSeed` | 1e-6 | A better starting point for a retry from scratch, which considers which phases are likely to be stable. | Leave on. Systems with many competing phases where the first start from scratch fails. | 0 = off. |
| `pa_OptimaDimReduceTol` | 10 | How generous the first guess of "species that matter" is (larger includes more). | Raise it if the reduced first solve misses species that end up present. | A negative value −m = start with the m×N cheapest species. |
| `pa_OptimaPreSolveFirstIters` | 6000 | Iteration limit for each pass of the reduced first solve, so a pass that goes nowhere is abandoned. | Large systems where the reduced first solve takes long without converging. | 0 = no limit. |
| `pa_OptimaEarlyStabilityAt` | 0 (automatic) | Lets Optima notice early that a phase should disappear, instead of shrinking it slowly. | Systems where a phase disappears slowly over hundreds of iterations. | A positive number = check at that iteration; a negative number −N = check after a phase has shrunk N times in a row. |
| `pa_OptimaDcFloor` | 0 (automatic) | The smallest amount Optima can represent for an absent species. Automatic derives it from the mass-balance tolerance `pa_DHB`. | Systems with elements present only in traces, where the automatic value is too coarse. | Any positive amount in moles. |
| `pa_PhaseHessianFloor` | 0.01 | Keeps the curvature of mixed phases well-behaved. | Solid solutions or melts that tend to split into two phases (miscibility gaps). | 0 = off. |
| `pa_OptimaFDHessian` | 1 (on) | Computes a more exact curvature for non-ideal phases, at extra cost. | Turn off to speed up systems with weakly non-ideal phases; keep on for strongly non-ideal ones. | 0 = off (faster, sometimes less robust). |
| `pa_OptimaFDHessianDelay` | 0 | Uses the cheaper approximate curvature for this many iterations first. | Large non-ideal systems where the early iterations are expensive. | Any positive number of iterations. |
| `pa_OptimaMoleFracHessian` | 0 (off) | A more exact curvature for ideal mixed phases. | Try it if a system with ideal solid solutions converges slowly; it helps some systems and slows others. | 1 = on. |
| `pa_LogBarrierTau` | 1e-16 | A tiny numerical push that keeps pure minerals away from exactly zero. | Rarely needs changing. | Any small positive number. |

### GEMS3K IPM solver — long-standing settings

These settings existed before this version and keep their meaning. Most projects carry them
in their file already.

| Setting | Default | What it does | When it can be useful | Other values |
|---|---|---|---|---|
| `pa_DK` | 1e-6 | How precisely the GEMS3K IPM solver must converge. Smaller gives a more precise answer at the cost of more iterations. | Solid solutions whose end-members are nearly alike: 1e-7 gives a noticeably more accurate composition. | Any small positive number. |
| `pa_DHB` | 1e-13 | How closely the amount of each element in the answer must match the amount put in, as a fraction of that amount. | Loosen it slightly if a calculation fails only because of this check; tighten it for systems with trace elements. | Usually between 1e-9 (looser) and 1e-15 (stricter). |
| `pa_DT` | 0 | How the check above is applied. 0 applies it the same way to every element. | Systems mixing major and trace elements where the relative check is too strict for the major ones. | A value of −6 or less = for major elements, use a fixed amount instead of a fraction. |
| `pa_IIM` | 7000 | Maximum number of iterations of the main calculation before it gives up. | Raise it for large or difficult systems that run out of iterations. | Up to 9999. |
| `pa_DP` | 130 | Maximum number of iterations of the step that makes the element totals add up. | Raise it if that step runs out of iterations on a difficult system. | Any positive number. |
| `pa_DW` | 1 (on) | Treats running out of the iterations above as an error. | Turn off only to inspect a result that would otherwise be rejected. | 0 = do not treat it as an error. |
| `pa_PE` | 1 (on) | Requires the answer to be electrically neutral. | Leave on for any system with ions; off only for systems with no charged species. | 0 = off. |
| `pa_PC` | 2 | How the solver decides which phases are present. 2 is the current method, which also adds back phases that were wrongly left out. | Leave at 2. The old method is kept for comparison with old results. | 1 = the old method. |
| `pa_DF` | 0.01 | How clearly a missing phase must be stable before it is added to the answer. | Lower it if phases that should appear are missed; raise it if phases flicker in and out along a sweep. | Any small positive number; smaller adds phases more readily. |
| `pa_DFM` | 0.01 | How clearly a present phase must be unstable before it is removed from the answer. | Adjust together with `pa_DF` when phases flicker in and out along a sweep. | Any small positive number. |
| `pa_PD` | 2 | How often activity coefficients (the corrections for non-ideal behaviour) are recalculated. 2 = at every iteration. | Leave at 2. Lower values can save time on nearly ideal systems. | 0 = once at the start; 1 = only while fixing the element totals; 3 = only in the main calculation. |
| `pa_AG` | 1 | Damping of the non-ideal corrections between iterations. | Lower it if a strongly non-ideal system oscillates and fails to converge. | From −1 to 1. |
| `pa_DGC` | 0 | A second damping setting that works together with the one above. | Use 0.001 for systems with sorption. | From −1 to 1. |
| `pa_PRD` | −5 | A final clean-up that removes tiny amounts of species left over at the end of the calculation. −5 means amounts below 1e-5 are cleaned up. | Make it more negative if real trace amounts are being removed; 0 to see the raw result. | 0 = off; −6 or less = clean up only smaller amounts. |
| `pa_GAS` | 0.001 | How strict the final clean-up is. A larger value keeps more trace amounts in the answer. | Systems with very little water (0.002 together with `pa_OptimaAcceptRepair`), or when trace phases matter. | Any small positive number. |
| `pa_DS` | 1e-20 | Smallest amount of a phase, in moles, that is still reported as present. | Raise it to hide meaningless trace phases in the output. | Any small positive amount. |
| `pa_XwMin` | 1e-13 | If the amount of water falls below this (moles), the aqueous solution is removed from the answer. | Drying or nearly dry systems, to control when the aqueous solution disappears. | Any small positive amount. |
| `pa_ScMin` | 1e-13 | If the amount of a solid that carries a sorption surface falls below this (moles), the sorption phase is removed. | Sorption systems where the sorbent dissolves almost completely. | Any small positive amount. |
| `pa_PhMin` | 1e-20 | If a solution phase other than the aqueous one falls below this amount (moles), it is removed with all its species. | Systems with many solid solutions or melts present in trace amounts. | Any small positive amount. |
| `pa_DcMin` | 1e-33 | If a species inside a mixed phase falls below this amount (moles), it is removed. | Rarely needs changing. | Any small positive amount. |
| `pa_DB` | 1e-17 | Smallest amount of an element, in moles, that the solver works with in the bulk composition (charge excluded). | Systems with elements at extremely low amounts. | Any small positive amount. |
| `pa_ICmin` | 1e-5 | Below this ionic strength (molal), activity coefficients of aqueous species are taken as 1, i.e. the solution is treated as ideal. | Very dilute waters, if the ideal treatment starts too early or too late. | Any small positive number. |
| `pa_DG` | 1000 | The solver internally scales the whole system to this total number of moles. Results are reported in the original amounts. | Very small or very large systems; switch off only to compare with unscaled results. | Below 1e-4 = no scaling. |
| `pa_EPS` | 1e-10 | How precisely the automatic starting point (used by AIA) is computed. | Loosen it if the starting point cannot be found; tighten it for systems with trace elements. | From 1e-6 (looser) to 1e-14 (stricter). |
| `pa_DFYw`, `pa_DFYaq`, `pa_DFYid`, `pa_DFYr`, `pa_DFYh`, `pa_DFYc` | 1e-5 | Small amounts (moles) given to water, aqueous species, species of ideal and non-ideal mixed phases, and pure phases that are zero in the automatic starting point, so the main calculation can work with them. | Lower them for systems with elements in traces, so the start does not disturb the element totals. | Any small positive amount. |
| `pa_DFYs` | 1e-6 | Amount (moles) at which a pure phase is added back when the solver decides it was wrongly left out. It is never more than the bulk composition can supply. | Rarely needs changing. | Any small positive amount. |
| `pa_PLLG` | 30000 | Checks whether the calculation is drifting off in the wrong direction and stops it if so. 1 to 1000 is the useful range for this check. | Lower it to catch failing calculations earlier; 0 if the check stops calculations that would succeed. | 0 = no check; 30000 or more also allows full diagnostic tracing. |
| `pa_PSM` | 1 | How many diagnostic messages are written. | Use 2 or 3 when investigating a failed calculation; 0 for large batch runs. | 0 = none (no log file); 2 = also warnings; 3 = detailed trace. |
| `pa_DNS` | 12.05 | Standard density of surface sites (per nm²), used for the activity of surface species in sorption models. | Sorption models that use a different standard site density. | Any positive number. |
| `pa_IEPS` | 0.001 | How precisely the surface terms of sorption models are computed. | Tighten it for sorption systems that converge poorly. | From 0.01 to 1e-6. |
| `pa_DKIN` | 1e-10 | Tolerance on species amounts that are held between an upper and a lower limit (for example in kinetic or metastability calculations). | Kinetic or metastability calculations with very small allowed ranges. | Any small positive amount. |

### Original settings not listed above

These settings come from earlier versions and are unchanged. Most are written by GEMS when it exports a
project, and are best left as exported.

| Setting | Typical value | What it does | When it can be useful | Other values |
|---|---|---|---|---|
| `pKin` | 1 | Lets the solver apply the limits on species amounts that come from the kinetics and metastability constraints of the project. | Set 0 to ignore all such limits and compute the full equilibrium. | 0 = ignore the limits. |
| `PV` | 1 | Adds the volume as a balance constraint, used for indifferent equilibria at saturated vapour pressure. | Leave as exported; relevant only for systems with a gas or fluid at saturation pressure. | 0 = off. |
| `PAalp` | `+` | Uses the specific surface area of phases in the calculation. | Set `-` to ignore surface areas, for example to compare with a run that has none. | `+` = use, `-` = ignore. |
| `PSigm` | `+` | Uses the specific surface free energy of phases in the calculation. | Set `-` to ignore it. | `+` = use, `-` = ignore. |
| `tMin` | 0 | Reserved for choosing which thermodynamic potential is minimised. Only the Gibbs energy is in use. | Do not change. | — |
| `PSOL` | 0 | Reserved (number of species in liquid hydrocarbon phases). | Do not change. | — |
| `pa_GAR`, `pa_GAH` | 1, 1000 | Reserved, no effect. | Do not change. | — |

Also in the file, and not settings to tune: `Lads`, `FIa` and `FIat` (sizes of the sorption part of the
system), `sMod`, `LsMod`, `IPxPH` and `PMc` (which mixing model each phase uses, and its parameters), and
`ID_key` (the name of the system). They describe the project's contents. Changing them by hand corrupts
the project; change them only by exporting the project again from GEMS.

### GEMS3K IPM solver — settings added in this version

| Setting | Default | What it does | When it can be useful | Other values |
|---|---|---|---|---|
| `pa_ColdRetryNudges` | 4 | If a calculation from scratch (AIA) fails, retry up to this many times with a microscopically changed composition, then finish at the exact one. | Leave on. Systems where AIA fails occasionally. | 0 = no retry. |
| `pa_IpmStallWindow` | 30 | Stops once the energy and composition have clearly stopped changing for this many iterations. | Leave on. Saves iterations on systems that converge slowly at the end. | 0 = off. |
| `pa_IpmAugmentedKKT` | 2 | How each step's equations are solved. 2 is the most accurate. | Leave at 2; other values only to compare with older results. | 0 = the original method; 1 = an intermediate method. |
| `pa_MbReproject` | 1 (on) | A final touch-up that makes the element totals add up exactly. Also used after phase selection, where a partly successful touch-up is now kept because it only seeds the next iterations. | Leave on. Turn off only to see the raw result. | 0 = off. |
| `pa_PSTALL` | 1 (on) | Lets the mass-balance step give up early when it stops improving. | Leave on. Saves time on systems where that step stalls. | 0 = off. |
| `pa_MbClassRule` | 0 (off) | Checks the totals of major elements in absolute terms and of trace elements in relative terms; the value is the trace/major ratio. | Systems with trace elements at very low amounts. | A positive ratio = on. |
| `pa_FilloutBudget` | 0 (off) | Limits how much the starting guess may disturb the element totals, as a fraction of each element's amount. | Systems with elements present only in traces. | A positive fraction = on. |
| `pa_DeterminacyWarn` | 0.01 | Warns when the amount of a phase in the answer is fixed by the energy only to worse than this relative uncertainty, i.e. when small amounts are not reliable. | Leave on; it tells you when reported small amounts should not be relied on. | 0 = no warning. |
| `pa_StabTPD` | 1 (on) | Reports, in the diagnostic trace only, whether a mixed phase left out of the answer should have been present. Does not change the answer. | Investigating whether a solid solution or melt was wrongly left out. | 0 = off. |

**Removed settings.** These names are still accepted in project files, so existing files load unchanged,
but they no longer have any effect: `pa_LpDualFillout`, `pa_OptimaLSWindow`, `pa_IpmLoopTweaks`,
`pa_OptimaFDDiagFloor`, `pa_MbPivotSplit`, `pa_OptimaPhaseCompaction`, `pa_OptimaReadmitSeed`,
`pa_MbTrendPhaseDecay`, `pa_OptimaMaxStepRatio`, `pa_GAR`, `pa_GAH`.

## Known limitations

- **Optima is slower per calculation** than the GEMS3K IPM solver, often by one to two orders of magnitude on
  large systems. Use it where it improves the answer, not by default. This will be improved in future releases. 
- **Some strongly non-ideal systems** may still be difficult to solve, even with Optima.
- **SIA along a sequence** can keep phases that are no longer stable when a phase boundary is crossed.
- **Systems with very little water:** close to the point where water runs out, AOP and SOP can fail and SHP
  can return a wrong result. Use AIA or HOP there. Slightly further from that point, AOP and SOP work if the
  project sets `pa_OptimaAcceptRepair = 1` and `pa_GAS = 0.002`.
- **Redox not fixed by the system.** Without a redox couple, the energy and amounts are well defined but Eh
  is not, and it can differ between runs and modes. Compare solvers on energy and phases, not Eh.
- **Eh and pH can be unstable** under tiny input changes while the phases and energy stay the same.
- **No free water solution** (very dry systems): pH, Eh and ionic strength are not meaningful and can differ
  widely between solvers. The phases are still correct.
- **Trace elements:** their mass balance can be relatively off in the Optima modes while the energy and the
  main elements are exact.
- **Solid solutions near a miscibility gap** and duplicated twin phases cost many more iterations, and the
  twin labels can swap between modes. The energy is the same.
- **GEMS3K IPM solver iteration counts** vary on many systems with tiny input changes. The answer does not.
- **Trace phases:** the GEMS3K IPM solver can leave out a slightly supersaturated phase that Optima finds.
- **Large systems** can cost tens of seconds per point in the Optima modes, and some fail where the GEMS3K IPM
  solver succeeds.

## Compatibility

Project files from earlier versions load unchanged. The order of settings in the project files is the
same, so files written by this version can be read by programs that use the same format.

## Technical note: logging and tracing

Nothing is traced by default. The only output is the normal log written through spdlog (messages and warnings
on the screen, and the `ipmlog.txt` file for the GEMS3K IPM solver). Everything below is switched on by the
user and costs nothing while it is off.

**Log level (spdlog).** The log is written by named loggers: `gems3k`, `ipm`, `tnode`, `kinmet`, `solmod`
(and `chemicalfun`, `thermofun` when ThermoFun is used). Their level is set in the file `gems3k-config.json`, in
the section `log`. Levels are `trace`, `debug`, `info` (default), `warning`, `error`, `critical` and `off`:

```json
{
  "log": {
    "level": "info",
    "logs-directory": "logs",
    "module_level": { "ipm": "debug", "tnode": "warning" }
  }
}
```

**What to use.** The default is `info`, which is right for most uses: warnings and errors are shown, and the
detailed messages stay off. You need no `gems3k-config.json` for this. Use `warning` when running many
calculations in a loop (a sweep, or a coupled transport code), `debug` for one logger (usually `ipm`) when
investigating a failed calculation, and `error` with `pa_PSM = 0` for very large batches.

The file is optional: without it every logger runs at `info`. The repository contains one example,
`tools-build/gems3k-config.json`, meant for development. It sets `gems3k` to `debug` and `thermofun` and
`chemicalfun` to `trace`, and also writes a log file (`"file"` section with `path`, `size`, `count`), so do
not copy it for normal use.

`level` applies to all loggers, and `module_level` overrides it for one logger. `logs-directory` is where
the log files go. The level can also be changed in code with `gems3k_update_log_level()`. The setting
`pa_PSM` in the project file controls what the GEMS3K IPM solver writes to `ipmlog.txt` (0 = nothing,
2 = also warnings, 3 = detailed trace, see [Settings](#settings)).

**Event trace (environment variable).** Set `GEMS3K_NATIVE_TRACE_FILE=<path>` before the program starts. The
solver then appends a record for each call to that file: the settings in force, the bulk composition, the
decisions the solver took, and the outcome (phases present with amounts, pH, pe, ionic strength). If the
variable is not set, nothing is written. The file can be tens of megabytes on a long run, so write it to
scratch space and delete it afterwards.

**Further switches for investigating one problem.** Each is an environment variable that takes a file path,
is read once, and writes only when set:

| Variable | What it writes |
|---|---|
| `GEMS3K_IPM_PROBE` | Per-iteration progress of the GEMS3K IPM solver |
| `GEMS3K_PHSTAB_PROBE` | Per-phase stability verdicts of the Optima path |
| `GEMS3K_EARLYTREND_PROBE` | What the early stability check would see |
| `GEMS3K_OPTIMA_TRACE_FILE` | Optima's own iteration table for the reduced problem |

A few more instruments (condition-number estimates, per-phase timing) exist only in a build made with
`-DENABLE_BENCHMARK_DIAGNOSTICS=ON`; they add run time and are for benchmarking.
