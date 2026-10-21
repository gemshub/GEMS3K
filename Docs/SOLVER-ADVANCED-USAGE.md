# GEMS3K 5.2.0 — advanced usage

This page is for batch runs, reactive transport, grids and other special cases. For the solvers, the modes and the
user settings, see the [options guide](SOLVER-OPTIONS-GUIDE.md). For the settings that are for tests, see
[Advanced settings](SOLVER-ADVANCED-SETTINGS.md). The [Glossary](GLOSSARY.md) explains the special words.

To change a setting, edit the `-ipm` file of the project (see [Settings](SOLVER-OPTIONS-GUIDE.md#settings)).

## Suggested settings for common cases

Most projects need no change. The table gives the settings that help in some common cases. They are suggestions
for one project or one run. They do not change the defaults.

| Case | Settings and mode | Why |
|---|---|---|
| Reactive transport and batch runs: each node must have a result | Keep the defaults `pa_FailRescue = 2` and `pa_RescueBudget = 100000`. Solve the first point cold and the next points warm (SIA, SOP). `pa_PSM = 0` for large runs. | The rescue tries a failed node again in steps. The limit stops a node that no step can save. `pa_PSM = 0` writes no log file. |
| Some calls run for a very long time | `pa_OptimaCallCap = -1` | Sets a cap from the size of the system. The results do not depend on the computer. `pa_OptimaMaxSeconds` also limits each solve, but then the results depend on the speed of the computer. |
| Phase diagrams and Pourbaix maps where a failed point is acceptable | `pa_FailRescue = 0`. Use AOP or SOP, not SIA, along a row. | A point that cannot be solved, for example a target outside the water window of the model, then costs much less. You get the error instead of the best attempt. SIA can keep phases after a phase boundary. |
| Water used up by design (for example cement hydration with little water) | `pa_OptimaDryOut = 1`. If AOP or SOP still fails: `pa_OptimaAcceptRepair = 1` and `pa_GAS = 0.002`, or use AIA or HOP. | Optima solves the solids without the water phase when the water runs out. |
| No redox couple in the system | `pa_OptimaRedoxFree = 2`, or give an Eh target. | The bulk composition does not fix the Eh. The setting aims at pe 4. |
| Solid solutions with end-members that are almost the same | `pa_DK = 1e-7` | Gives a clearly more accurate composition with the GEMS3K IPM solver. |
| Sorption models | `pa_DGC = 0.001` | Damps the non-ideal corrections between iterations. |
| Absent phases must show exactly zero | `pa_OptimaZeroAbsent = 1` | Optima keeps an absent species at a very small floor amount by default. |
| Large systems (hundreds of species) | AIA, or HOP for a checked result. With Optima: `pa_OptimaTol` a little larger. | Optima alone is much slower on large systems. |
| Trace elements | A stricter `pa_DHB`. Give each element more than the floor amounts of its species need (see [Known limits](SOLVER-OPTIONS-GUIDE.md#known-limits)). | Keeps the mass balance of small amounts exact. |
| Find the cause of a failed calculation | `pa_PSM = 2` or `3` | Writes warnings or a detailed trace to `ipmlog.txt`. |
