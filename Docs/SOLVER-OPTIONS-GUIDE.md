# GEMS3K 5.2.0 — options guide

This guide tells you what is new in this version. It tells you how to choose a solver and how to run it from
code. It tells you how to change a setting in the project file. It describes each setting and the known
limits.

Words with a special meaning are in the [Glossary](GLOSSARY.md). To see what runs in each mode, in which order, and what GEMS3K tries when a calculation fails, read
[Solver paths](SOLVER-PATH-GRAPH.md) ([drawn version](SOLVER-PATH-GRAPH.html)).

## The two solvers

GEMS3K calculates chemical equilibrium. It finds the lowest Gibbs energy of the system. GEMS3K has two
solvers. You choose one for each calculation (see [The calculation modes](#the-calculation-modes)).

**The GEMS3K IPM solver** is the solver that GEMS3K always had. It uses an interior-point method (IPM). A
refinement of the mass balance follows. This solver is fast on almost all systems. It is the default (modes
AIA and SIA). Its settings are the `pa_...` entries in the `-ipm` file of the project:
- [long-standing settings](#gems3k-ipm-solver--long-standing-settings)
- [original settings](#original-settings-not-listed-above)
- [settings added in this version](#gems3k-ipm-solver--user-setting-added-in-this-version)

**The Optima solver** is a second solver. It is based on the open-source Optima optimisation library
([github.com/gemshub/optima](https://github.com/gemshub/optima)). The GEMS3K team made changes for GEMS3K.
- Optima is more robust for solid solutions, melts and other systems where several mixing phases compete
  (modes AOP and SOP).
- Optima can also check the result of the GEMS3K IPM solver (modes HOP and SHP).
- Optima can fix the pH or the Eh. See [Fixing the pH or Eh](#fixing-the-ph-or-eh-optima-modes-only).

The settings of Optima are also `pa_...` entries in the `-ipm` file. Their names start with `pa_Optima`. See
[user settings](#optima-solver--user-settings). Settings for tests are in
[Advanced settings](SOLVER-ADVANCED-SETTINGS.md).

The Optima library has more options. You cannot set them in the GEMS3K project files. GEMS3K chooses them.
The options guide of the Optima library describes them (`docs/OPTIONS-GUIDE.md` in its repository).

Some general settings apply to both solvers. Examples are `pa_PE`, `pa_DHB` and `pa_DG`. See
[Settings](#settings). That section also shows how to add a setting to a project file.

## What is new

**A second equilibrium solver.** GEMS3K can now calculate equilibrium with the Optima solver. Optima is more
robust for solid solutions, melts and other systems where several mixing phases compete. You choose the
solver for each calculation. Put a value in `NodeStatusCH` before you call `GEM_run()`. The GEMS3K IPM modes
are still the default.

**Fewer failed calculations.** A calculation can fail to converge. Before GEMS3K reports a failure, it now
tries again automatically. It uses one of these methods:
- a start that is changed by a very small amount
- a cold start instead of a warm start
- a finish on the phases that it already found

**A rescue for failed calls.** A call can still fail. Then GEMS3K tries it again in steps (`pa_FailRescue`).
If no step gives a good result, you get the best attempt that keeps the element totals. Its status is BAD.
`GEM_rescue_step()` tells you which step gave the result. See [Reading the result](#reading-the-result).

**Faster calculations.**
- The GEMS3K IPM solver stops when its energy does not improve. Optima then finishes the result.
- A cold call without aqueous, gas or sorption phases starts from a better estimate.
- Optima uses faster linear algebra on large systems. It does fewer steps that are not necessary.

The results do not change. These changes are on by default. The switches are in
[Advanced settings](SOLVER-ADVANCED-SETTINGS.md).

## The calculation modes

The name of each mode tells you how the calculation starts and which solver does the work.

| Value | Mode | Stands for | What it does |
|---|---|---|---|
| 1 | AIA | **A**utomatic **I**nitial **A**pproximation | GEMS3K IPM solver. It makes its own start from scratch. |
| 5 | SIA | **S**mart **I**nitial **A**pproximation | GEMS3K IPM solver. It starts from the previous result in the node. |
| 10 | AOP | **A**utomatic initial approximation with **OP**tima | Optima solver. It starts from scratch. |
| 14 | SOP | **S**mart initial approximation with **OP**tima | Optima solver. It starts from the previous result. |
| 22 | HOP | **H**ybrid with **OP**tima | GEMS3K IPM solver from scratch. Then Optima refines and checks the result. |
| 26 | SHP | **S**mart **H**ybrid with o**P**tima | The same as HOP. The GEMS3K IPM part starts from the previous result. |

In HOP, the Optima part can fail. Then GEMS3K keeps the GEMS3K IPM result and marks it. For this reason, HOP
does not give a worse result than the GEMS3K IPM solver alone.

### Older programs that use the numbers 10 and 11

In some older versions, status 10 meant "an external solver succeeded". Status 11 meant "it failed". In this
version, 10 requests an Optima cold start (AOP). Status 11 means that this solve succeeded.

An old program can read 11 as a failure. That program then treats a success of Optima as an error. Use the
named constants from `databr.h` (`NEED_GEM_AOP`, `OK_GEM_AOP` and so on). Do not use the raw numbers. Check
all old code that compares a status with 10 or 11.

## Which mode to use

| Situation | Mode | Reason |
|---|---|---|
| Most systems, single calculations | **AIA** | The GEMS3K IPM solver is the fastest on almost all systems that it can solve. It is often much faster. |
| Systems with many species (hundreds or more) | **AIA**, or **HOP** if you want a checked result | Optima alone is much slower on large systems. HOP gives the check of Optima. It costs a little more than the GEMS3K IPM solver. |
| Strongly non-ideal systems: solid solutions, melts, miscibility gaps, phase diagrams | **AOP**, or **HOP** | Optima finds the correct set of stable phases more often when several mixing phases compete. |
| The GEMS3K IPM solver fails, or warns that its result does not satisfy the mass balance | **AOP** or **HOP** | Optima solves most of these cases. |
| Series of small changes: sweeps in temperature or composition, titrations, transport steps | **SOP** after a first point in AOP (or **SIA** after AIA) | A warm start costs much less than a cold start. |
| Series that cross phase boundaries | **AOP** or **SOP**. Do not use **SIA**. | SIA can keep the phases of the previous point after they are not stable. |
| Systems with very little water (for example, nearly dry cement) | **AIA** or **HOP** | Optima can fail when the water is almost gone. The GEMS3K IPM solver handles these systems. |

[Solver paths](SOLVER-PATH-GRAPH.md) shows the steps that each mode runs. It also shows the tries that each
mode makes before it reports a failure.

## Use from code

GEMS3K reads the project files one time, in `GEM_init()`. No call changes a setting after that. To try
another value, edit the `-ipm` file (see [Settings](#settings)). Then initialise a new node.

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

// A second, similar calculation: start from the previous result (warm start).
node.pCNode()->NodeStatusCH = NEED_GEM_SOP;    // set it again before EVERY call
status = node.GEM_run(false);

// Fix the pH and let the solver add the acid or base that it needs (Optima modes only).
node.Set_pH_target(7.0);
node.pCNode()->NodeStatusCH = NEED_GEM_SOP;
status = node.GEM_run(false);
```

Remember these two points:
- `GEM_run()` does nothing unless `NodeStatusCH` has a `NEED_GEM_...` value. Set the value just before each
  call.
- You can use the same node for all modes. To compare two settings, make two nodes from two copies of the
  project. Use one `-ipm` file for each node.

A complete example is `tools/optima_test.cpp`.

### Reading the result

- `OK_GEM_...`: the result passed all checks.
- `BAD_GEM_...`: there is a result, but it did not pass all checks. Read the log.
- `ERR_GEM_...`: there is no result. The values in the node are from before the call.
- `GEM_rescue_step()` gives the rescue step that made the result:
  - 0: the call itself
  - 1: a warm restart from the previous result
  - 2: a composition changed by a very small amount
  - 3: a walk from the previous result
  - 4: looser settings
  - 5: the other solver
  - 6: the best attempt (the status is always BAD)

```cpp
long status = node.GEM_run(false);
if (node.GEM_rescue_step() == 6)
    std::cerr << "no good result; best balanced attempt returned\n";
```

### Which phases are present

The Optima solver never sets an amount to exactly zero. A phase that does not form stays at a very small floor
amount. The floor is about the size of `pa_DHB` for each species (for example, 1e-13 mol). The GEMS3K IPM
solver removes such phases. For this reason, one result can show more phases with an amount above zero in an
Optima mode than in a GEMS3K IPM mode.

- **Old rule, `pa_DS`.** A phase is present when its amount is more than `pa_DS` (default 1e-20 mol). The
  GEMS3K IPM solver uses this rule to remove phases. Graphical front ends use it to show the phases that are
  present. The value is much lower than the floor of Optima. With Optima, each floor phase passes this rule.
- **New rule, `Ph_Present(k)`.** A phase is present when it holds more than 1e-9 of the total amount of all
  phases. A second argument changes this share. The amounts do not change. With this rule, the Optima modes
  and the GEMS3K IPM modes show the same phases.

```cpp
for (long k = 0; k < node.pCSD()->nPHb; k++)
    if (node.Ph_Present(k))
        std::cout << "phase " << k << ": " << node.Ph_Moles(k) << " mol\n";
```

## Fixing the pH or Eh (Optima modes only)

In the Optima modes AOP, SOP, HOP and SHP, you can request a pH, an Eh (in volts), or both. The solver adds
the acid, base or electron donor that the target needs. It reports how much it added.

- In HOP and SHP, the Optima part meets the target. The GEMS3K IPM part runs first and ignores the target.
- This function is not available in AIA or SIA.
- This function needs a build with Optima.

```cpp
node.Clear_ControlConditions();           // targets stay set until cleared
node.Set_pH_target(7.0);                  // pH 7
node.Set_Eh_target(0.2);                  // Eh 0.2 V (optional; can be used alone)

node.pCNode()->NodeStatusCH = NEED_GEM_AOP;   // or SOP, HOP, SHP; set before every call
long status = node.GEM_run(false);

if (status == OK_GEM_AOP)
    std::cout << "pH " << node.Get_pH() << "  Eh " << node.Get_Eh()
              << "  H+ removed (mol) " << node.Get_ControlCondition_titrant("pH")
              << "  electrons removed (mol) " << node.Get_ControlCondition_titrant("Eh") << "\n";
```

- The tolerance is an optional second argument: `Set_pH_target(7.0, 1e-3)`. If you do not give it, GEMS3K
  derives it from `pa_GAS`.
- The amount that the solver added is part of the result. The bulk composition that the node returns
  includes it. To start the next point from the original composition, first write the original bulk
  composition back.
- A target stays active for each later call until you use `Clear_ControlConditions()`.
- A target can be impossible to reach. An example is an Eh in a system with no redox couple. Check the
  status. Do not trust the Eh in such a system.
- A target can use all the water (for example, an extreme pH). GEMS3K then reports the result as not fully
  trustworthy (`BAD_GEM_...`).
  - A failed SOP or SHP call is solved again with a cold start. It can return `OK_GEM_AOP` (see
    `pa_OptimaColdRetry`).
  - A failed HOP call returns the GEMS3K IPM result without the target.
- An `OK` status means that the solver met the target in its own model. An extreme target can still give a
  result that has no physical meaning (for example, a pH far outside 0 to 14). Judge the target. Do not
  only look at the status.

## Settings

GEMS3K reads the settings from the `pa_...` entries in the `-ipm` file of the project. A setting that is not
in the file keeps its default value. The tables list the settings that a user can change. Settings for tests
are in [Advanced settings](SOLVER-ADVANCED-SETTINGS.md).

**How to add a setting that is not in the file.** Files that earlier versions exported do not contain the new
settings. To use a value other than the default, add the line to the file. Put it next to the other `pa_`
entries.

In a JSON file (`…-ipm.json`), add the setting in the `"ipm"` block. For example, add it after `"pa_PLLG"`:

```json
    "pa_PLLG": 10000,
    "pa_OptimaTol": 1e-09,
    "pa_OptimaMaxSeconds": 60,
    "tMin": 0,
```

In a text file (`…-ipm.dat`), add the setting after the `<END_DIM>` line. For example, add it after
`<pa_PLLG>`:

```
<pa_PLLG>  30000
<pa_OptimaTol>  1e-09
```

Put each setting in the file one time only. If a name is in the file two times, GEMS3K uses the later value.
GEMS3K does not read a setting that is before `<END_DIM>` in a text file.

### Optima solver — user settings

These are the Optima settings that a user can change. Advanced settings has all other Optima settings
([Advanced settings](SOLVER-ADVANCED-SETTINGS.md)).

| Setting | Default | What it does | When it can help | Other values |
|---|---|---|---|---|
| `pa_OptimaTol` | 1e-8 | Sets how precisely Optima must converge. | Make it a little larger to get more speed and less accuracy on large systems. Make it smaller only if a result is not converged enough. | Any positive number. A smaller number is stricter. |
| `pa_OptimaMaxSeconds` | 0 (off) | Sets a time limit in seconds for one calculation with its tries. The results then depend on the speed of the computer. | Long batch runs or transport runs. One stuck calculation must not stop the other calculations. | Any positive number of seconds. |
| `pa_OptimaDimReduce` | 0 (automatic) | For large systems, the solver first solves with only the important species. Then it adds the other species. Automatic means that it is on from 200 species. | Switch it on for a slow system with fewer than 200 species. Switch it off if a large system gives a different result with it. | A positive number = on, with at most that number of passes. A negative number = off. |
| `pa_OptimaColdRetry` | 2 (fail fast) | A warm-start calculation (SOP, SHP) can fail. Then the solver makes a new cold start. The value 2 stops the first attempt early. | Use 1 if the warm start usually succeeds but needs many iterations. Use 0 only for tests. | 0 = no new try. 1 = a new try only after the first attempt used its full budget. 2 = the first attempt stops early. |
| `pa_OptimaAcceptRepair` | 0 (off) | Optima can accept a result that is almost converged. Then the solver first corrects a small error in the element totals. | Systems with very little water. Other cases where Optima ends near the result with element totals that are not exactly balanced. | 1 = on. |
| `pa_OptimaZeroAbsent` | 2 | Sets how absent species are shown. They are shown as a very small amount. If necessary, the solver corrects the element totals. | Use 1 if absent phases must show exactly zero in the output. | 0 = very small amounts without a change. 1 = absent phases set to exactly zero. |
| `pa_OptimaRedoxFree` | 1 (write to the log) | Checks before the calculation if the bulk fixes the redox level. This is the case if each element has one oxidation state and no spare oxidant or reductant exists. If the bulk does not fix the level, the solver writes this in the log. | Systems with no redox buffer. In these systems, the input does not fix Eh and the species that depend on it. | 0 = off. 2 = also aim the calculation at pe 4 when you give no pH target or Eh target. |
| `pa_OptimaCallCap` | 0 (no cap) | Sets the largest number of iterations for one calculation (AOP, SOP). After that number, the solver repeats the calculation with amounts that are changed by a very small amount. | Batch runs where a calculation sometimes runs for a very long time. | -1 = chosen from the size of the system. A positive number = that number of iterations. |
| `pa_OptimaDryOut` | 0 (off) | The water can run out during a calculation. Then the solver takes the water phase out and solves the solids alone. It compares this result with a solve that keeps the water. It keeps the result with the lower energy. If the dry solve fails, the calculation continues with the water back. The solver does not use this if an element can only be in the water. It needs `pa_OptimaStallWindow` to be on. | Use 1 for systems where the water is used up by design, for example cement hydration with little water. It is off by default. If the water only becomes low during a calculation, the dry result can be wrong. | 0 = off. 1 = on. |

### GEMS3K IPM solver — long-standing settings

These settings were in earlier versions. They have the same meaning. Most projects have them in the file.

| Setting | Default | What it does | When it can help | Other values |
|---|---|---|---|---|
| `pa_DK` | 1e-6 | Sets how precisely the GEMS3K IPM solver must converge. A smaller value gives a more precise result. It needs more iterations. | Solid solutions with end-members that are almost the same. The value 1e-7 gives a clearly more accurate composition. | Any small positive number. |
| `pa_DHB` | 1e-13 | Sets how closely the amount of each element in the result must agree with the input amount. The value is a fraction of that amount. | Make it a little looser if a calculation fails only because of this check. Make it stricter for systems with trace elements. | Usually from 1e-9 (looser) to 1e-15 (stricter). |
| `pa_DT` | 0 | Sets how the check above is used. The value 0 uses it in the same way for each element. | Systems with major elements and trace elements. The relative check can be too strict for the major elements. | A value of −6 or less = for major elements, use a fixed amount instead of a fraction. |
| `pa_IIM` | 7000 | Sets the largest number of iterations of the main calculation. After that number, the calculation stops. | Make it larger for large or difficult systems that use all their iterations. | Up to 32767. A larger value becomes 32767, with a warning. |
| `pa_DP` | 130 | Sets the largest number of iterations of the step that makes the element totals agree. | Make it larger if that step uses all its iterations on a difficult system. | Any positive number. |
| `pa_DW` | 1 (on) | Treats the use of all the iterations above as an error. | Switch it off only to look at a result that GEMS3K would reject. | 0 = not an error. |
| `pa_PE` | 1 (on) | Requires that the result has no electric charge. | Keep it on for each system with ions. Switch it off only for systems with no charged species. | 0 = off. |
| `pa_PC` | 2 | Sets how the solver decides which phases are present. The value 2 is the current method. It also adds phases that the solver wrongly left out. | Keep the value 2. The old method is only for the comparison with old results. | 1 = the old method. |
| `pa_DF` | 0.01 | Sets how clearly stable a missing phase must be before the solver adds it. | Make it smaller if phases that must appear are missing. Make it larger if phases go in and out along a sweep. | Any small positive number. A smaller number adds phases more easily. |
| `pa_DFM` | 0.01 | Sets how clearly unstable a present phase must be before the solver removes it. | Change it together with `pa_DF` if phases go in and out along a sweep. | Any small positive number. |
| `pa_PD` | 2 | Sets how often the solver calculates the activity coefficients again. These coefficients correct for non-ideal behaviour. The value 2 means at each iteration. | Keep the value 2. Lower values can save time on systems that are almost ideal. | 0 = one time at the start. 1 = only when it fixes the element totals. 3 = only in the main calculation. |
| `pa_AG` | 1 | Sets the damping of the non-ideal corrections between iterations. | Make it smaller if a strongly non-ideal system oscillates and does not converge. | From −1 to 1. |
| `pa_DGC` | 0 | A second damping setting. It works together with the setting above. | Use 0.001 for systems with sorption. | From −1 to 1. |
| `pa_PRD` | −5 | A final clean-up. It removes very small amounts of species that stay at the end of the calculation. The value −5 means that it cleans up amounts below 1e-5. | Make it more negative if real trace amounts are removed. Use 0 to see the raw result. | 0 = off. −6 or less = clean up only smaller amounts. |
| `pa_GAS` | 0.001 | Sets how strict the final clean-up is. A larger value keeps more trace amounts in the result. | Systems with very little water (0.002 together with `pa_OptimaAcceptRepair`). Cases where trace phases are important. | Any small positive number. |
| `pa_DS` | 1e-20 | Sets the smallest amount of a phase, in moles, that is still shown as present (the old rule, see [Which phases are present](#which-phases-are-present)). | Make it larger to hide trace phases that have no meaning. With Optima, use `Ph_Present()`. | Any small positive amount. |
| `pa_XwMin` | 1e-13 | The amount of water can become less than this value (moles). Then the solver removes the aqueous solution from the result. | Drying systems or almost dry systems. Use it to control when the aqueous solution goes away. | Any small positive amount. |
| `pa_ScMin` | 1e-13 | The amount of a solid with a sorption surface can become less than this value (moles). Then the solver removes the sorption phase. | Sorption systems where the sorbent almost fully dissolves. | Any small positive amount. |
| `pa_PhMin` | 1e-20 | A solution phase other than the aqueous phase can become less than this amount (moles). Then the solver removes it with all its species. | Systems with many solid solutions or melts in trace amounts. | Any small positive amount. |
| `pa_DcMin` | 1e-33 | A species in a mixed phase can become less than this amount (moles). Then the solver removes it. | You almost never need to change it. | Any small positive amount. |
| `pa_DB` | 1e-17 | Sets the smallest amount of an element, in moles, that the solver uses in the bulk composition (not the charge). | Systems with elements at extremely low amounts. | Any small positive amount. |
| `pa_ICmin` | 1e-5 | Below this ionic strength (molal), the solver sets the activity coefficients of aqueous species to 1. The solution is then ideal. | Very dilute waters, if the ideal treatment starts too early or too late. | Any small positive number. |
| `pa_DG` | 1000 | The solver scales the whole system internally to this total number of moles. It shows the results in the original amounts. | Very small or very large systems. Switch it off only to compare with results that are not scaled. | Below 1e-4 = no scaling. |
| `pa_EPS` | 1e-10 | Sets how precisely the solver calculates the automatic start (used by AIA). | Make it looser if the solver cannot find the start. Make it stricter for systems with trace elements. | From 1e-6 (looser) to 1e-14 (stricter). |
| `pa_DFYw`, `pa_DFYaq`, `pa_DFYid`, `pa_DFYr`, `pa_DFYh`, `pa_DFYc` | 1e-5 | Small amounts (moles). The solver gives them to water, to aqueous species, to species of ideal and non-ideal mixed phases, and to pure phases that are zero in the automatic start. The main calculation can then use them. | Make them smaller for systems with elements in traces. The start then does not disturb the element totals. | Any small positive amount. |
| `pa_DFYs` | 1e-6 | Sets the amount (moles) at which the solver adds a pure phase again when it decides that the phase was wrongly left out. The amount is never more than the bulk composition can supply. | You almost never need to change it. | Any small positive amount. |
| `pa_PLLG` | 30000 | The solver checks if the calculation goes in the wrong direction. It stops the calculation if it does. The useful range for this check is 1 to 1000. | Make it smaller to find failing calculations sooner. Use 0 if the check stops calculations that would succeed. | 0 = no check. 30000 or more also allows full diagnostic tracing. |
| `pa_PSM` | 1 | Sets how many diagnostic messages the solver writes. | Use 2 or 3 to find the cause of a failed calculation. Use 0 for large batch runs. | 0 = none (no log file). 2 = also warnings. 3 = detailed trace. |
| `pa_DNS` | 12.05 | Sets the standard density of surface sites (per nm²). The sorption models use it for the activity of surface species. | Sorption models that use a different standard site density. | Any positive number. |
| `pa_IEPS` | 0.001 | Sets how precisely the solver calculates the surface terms of sorption models. | Make it stricter for sorption systems that converge badly. | From 0.01 to 1e-6. |
| `pa_DKIN` | 1e-10 | Sets the tolerance on species amounts that stay between an upper limit and a lower limit (for example, in kinetic or metastability calculations). | Kinetic or metastability calculations with very small allowed ranges. | Any small positive amount. |

### Original settings not listed above

These settings are from earlier versions. They did not change. GEMS writes most of them when it exports a
project. Keep them as exported.

| Setting | Typical value | What it does | When it can help | Other values |
|---|---|---|---|---|
| `pKin` | 1 | Lets the solver use the limits on species amounts from the kinetics and metastability constraints of the project. | Use 0 to ignore all such limits and calculate the full equilibrium. | 0 = ignore the limits. |
| `PV` | 1 | Adds the volume as a balance constraint. It is for indifferent equilibria at saturated vapour pressure. | Keep it as exported. It is only for systems with a gas or fluid at saturation pressure. | 0 = off. |
| `PAalp` | `+` | Uses the specific surface area of phases in the calculation. | Use `-` to ignore surface areas, for example to compare with a run that has none. | `+` = use. `-` = ignore. |
| `PSigm` | `+` | Uses the specific surface free energy of phases in the calculation. | Use `-` to ignore it. | `+` = use. `-` = ignore. |
| `tMin` | 0 | Reserved. It is for the choice of the thermodynamic potential that the solver minimises. Only the Gibbs energy is in use. | Do not change it. | — |
| `PSOL` | 0 | Reserved (number of species in liquid hydrocarbon phases). | Do not change it. | — |
| `pa_GAR`, `pa_GAH` | 1, 1000 | Reserved. They have no effect. | Do not change them. | — |

The file also has `Lads`, `FIa` and `FIat` (sizes of the sorption part of the system). It has `sMod`,
`LsMod`, `IPxPH` and `PMc` (the mixing model of each phase and its parameters). It has `ID_key` (the name of
the system). Do not tune these entries. They describe the contents of the project. If you change them by
hand, the project is not correct. To change them, export the project again from GEMS.

### GEMS3K IPM solver — user setting added in this version

| Setting | Default | What it does | When it can help | Other values |
|---|---|---|---|---|
| `pa_DeterminacyWarn` | 0.01 | Writes a warning when the energy fixes the amount of a phase in the result only to worse than this relative uncertainty. Small amounts are then not reliable. | Keep it on. It tells you when you must not trust the reported small amounts. | 0 = no warning. |

### Both solvers — rescue of failed calls

| Setting | Default | What it does | When it can help | Other values |
|---|---|---|---|---|
| `pa_FailRescue` | 2 (on) | A call can end with an error or a BAD result. Then GEMS3K tries it again in steps. The steps are: a warm restart from the previous result, compositions that are changed by a very small amount, a walk from the previous result, looser settings, and the other solver. GEMS3K finishes each successful step at the exact composition with the normal settings. With the value 2, you get the best attempt that keeps the element totals if no step succeeds. Its status is BAD. If there is a pH target or an Eh target, GEMS3K keeps only an attempt that meets the target. | Keep the value 2. This is important for transport codes, which need a result for each node. It costs more iterations only on calls that fail. | 0 = off. 1 = steps 1 to 5 only. |

## Known limits

- **Optima is slower for each calculation** than the GEMS3K IPM solver. It is typically 2 times slower. It
  can be up to 10 times slower on large systems with many solid solutions. Use it where it gives a better
  result. Do not use it by default.
- **Some strongly non-ideal systems** can be difficult to solve, also with Optima.
- **SIA in a series** can keep phases that are not stable when the series crosses a phase boundary.
- **Systems with very little water.** Near the point where the water runs out, the modes can give different
  trace phases. Optima can need many more iterations.
  - For systems where the design uses up the water, set `pa_OptimaDryOut = 1`.
  - If AOP or SOP fails there, set `pa_OptimaAcceptRepair = 1` and `pa_GAS = 0.002`. You can also use AIA or
    HOP.
- **Redox not fixed by the system.** A system can have no redox couple. Then the energy and the amounts are
  well defined, but the Eh is not. The Eh can be different between runs and between modes. Compare solvers
  on energy and phases, not on Eh. To get a defined Eh, set `pa_OptimaRedoxFree = 2` (it aims at pe 4) or
  give an Eh target.
- **Eh and pH can be unstable** when the input changes by a very small amount. The phases and the energy do
  not change.
- **No free water solution** (very dry systems). The pH, the Eh and the ionic strength have no meaning. They
  can be very different between solvers. The phases are correct.
- **Trace elements.** In a few cases, the relative mass balance is not good in the Optima modes. The energy
  and the main elements are exact.
- **Solid solutions near a miscibility gap and twin phases that are duplicates** need many more iterations.
  The twin labels can change between modes. The energy is the same.
- **Iteration counts of the GEMS3K IPM solver** change on many systems when the input changes by a very small
  amount. The result does not change.
- **Trace phases.** The GEMS3K IPM solver can leave out a phase that is a little supersaturated. Optima finds
  it.
- **Large systems** can need tens of seconds for each point in the Optima modes.
- **Trace elements below the species floors.** The floor amounts of the species of a trace element can add up
  to more than the amount of the element. An Optima call can then stay BAD. The rescue cannot solve such a
  call. Use AIA or HOP.
- **Rescued calls cost more.** A call that the rescue saves (`pa_FailRescue`) can use several times the
  iterations of a normal call. A call that no step can save uses the most iterations.

## Compatibility

Project files from earlier versions load without change. The order of the settings in the project files is
the same. Programs that use the same format can read the files that this version writes.

## Technical note: logging and tracing

By default, GEMS3K does not trace anything. The only output is the normal log through spdlog. This is the
messages and warnings on the screen, and the `ipmlog.txt` file of the GEMS3K IPM solver. The user switches on
all items below. They cost nothing while they are off.

**Log level (spdlog).** Named loggers write the log: `gems3k`, `ipm`, `tnode`, `kinmet` and `solmod`. When
you use ThermoFun, the loggers `chemicalfun` and `thermofun` also write. Set the level in the file
`gems3k-config.json`, in the section `log`. The levels are `trace`, `debug`, `info` (default), `warning`,
`error`, `critical` and `off`:

```json
{
  "log": {
    "level": "info",
    "logs-directory": "logs",
    "module_level": { "ipm": "debug", "tnode": "warning" }
  }
}
```

**Which level to use.** The default level `info` is correct for most uses. It shows warnings and errors. It
does not show the detailed messages. You do not need a `gems3k-config.json` file for it.
- Use `warning` when you run many calculations in a loop (a sweep, or a coupled transport code).
- Use `debug` for one logger (usually `ipm`) to find the cause of a failed calculation.
- Use `error` with `pa_PSM = 0` for very large batches.

The file is optional. If it is not there, each logger runs at `info`. The repository has one example,
`tools-build/gems3k-config.json`. It is for development. It sets `gems3k` to `debug`. It sets `thermofun` and
`chemicalfun` to `trace`. It also writes a log file (section `"file"` with `path`, `size` and `count`). Do
not copy it for normal use.

`level` applies to all loggers. `module_level` replaces it for one logger. `logs-directory` is the folder for
the log files. You can also change the level in code with `gems3k_update_log_level()`. The setting `pa_PSM`
in the project file controls what the GEMS3K IPM solver writes to `ipmlog.txt` (0 = nothing, 2 = also
warnings, 3 = detailed trace, see [Settings](#settings)).

**Event trace (environment variable).** Set `GEMS3K_NATIVE_TRACE_FILE=<path>` before the program starts.
The solver then adds a record for each call to that file. The record has these items:
- the settings that are in use
- the bulk composition
- the decisions of the solver
- the outcome (phases that are present with amounts, pH, pe, ionic strength)

If the variable is not set, the solver writes nothing. On a long run, the file can be tens of megabytes. Write
it to scratch space. Delete it after use.

**More switches to find the cause of one problem.** Each switch is an environment variable that has a file
path. GEMS3K reads it one time. It writes only when you set the variable.

| Variable | What it writes |
|---|---|
| `GEMS3K_IPM_PROBE` | The progress of the GEMS3K IPM solver at each iteration. |
| `GEMS3K_PHSTAB_PROBE` | The stability result for each phase on the Optima path. |
| `GEMS3K_EARLYTREND_PROBE` | What the early stability check would see. |
| `GEMS3K_OPTIMA_TRACE_FILE` | The iteration table of Optima for the reduced problem. |
| `GEMS3K_OPTIMA_STEPTRACE_FILE` | A debug record of each Optima step. |
| `GEMS3K_TIMING_FILE` | The time of each call, divided by solver part. |

A build with `-DENABLE_BENCHMARK_DIAGNOSTICS=ON` has a few more tools (condition-number estimates and timing
for each phase). They add run time. They are for benchmarks.
