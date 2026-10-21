# GEMS3K 5.2.0 — advanced settings

These settings are for tests. Keep them at the default value. For the user settings, see the
[options guide](SOLVER-OPTIONS-GUIDE.md#settings).

Each setting here switches a correctness fix on or off, switches a speed tuning on or off, or sets an
internal limit. They are separate so that you can test each mechanism alone. A project file keeps one of them
only if its value is not the default. GEMS3K reads all of them.

The [Glossary](GLOSSARY.md) explains the special words.

| Setting | Default | What it does |
|---|---|---|
| `pa_OptimaFinish` | 1 (on) | Optima can stop just short of convergence. Then a final step finishes the job on the phases that Optima found. |
| `pa_OptimaTpdAccept` | 1e-6 (on) | Accepts a result that is not fully converged when a direct check shows that no missing phase must be present (within this tolerance). |
| `pa_OptimaStallWindow` | 200 (on) | Stops a calculation that does not improve for this number of iterations. A new try can then start sooner. |
| `pa_OptimaLineSearch` | 1.5 | A step can make the result worse by more than this factor. Then the solver tries a shorter step. |
| `pa_OptimaLSStallEscape` | 10 (on) | Shorter steps can fail this number of times in a row. Then the solver takes one full step to get out. |
| `pa_OptimaLSRejectWorse` | 0 (off) | If a shorter step does not help, the solver uses the full step. |
| `pa_OptimaCgSeed` | 1e-6 (on) | A better start for a new try with a cold start. It uses the phases that are likely to be stable. |
| `pa_OptimaDimReduceTol` | 10 | Sets how many species the first guess of "species that matter" includes. A larger value includes more species. |
| `pa_OptimaPreSolveFirstIters` | 6000 | Sets the iteration limit for each pass of the reduced first solve. The solver stops a pass that does not progress. |
| `pa_OptimaEarlyStabilityAt` | 0 (automatic) | Lets Optima find early that a phase must disappear. Without it, the phase becomes smaller slowly. |
| `pa_OptimaDcFloor` | 0 (automatic) | Sets the smallest amount that Optima can represent for an absent species. |
| `pa_PhaseHessianFloor` | 0.01 (on) | Keeps the curvature of mixed phases in a good range. |
| `pa_OptimaFDHessian` | 1 (on) | Calculates a more exact curvature for non-ideal phases. This costs more time. |
| `pa_OptimaFDHessianDelay` | 0 | Uses the cheaper approximate curvature for this number of first iterations. |
| `pa_OptimaMoleFracHessian` | 0 (off) | Uses a more exact curvature for ideal mixed phases. |
| `pa_OptimaFastRises` | 10 (on) | The quick first try of a cold-start calculation stops after this number of steps where the energy goes up by more than 1 %. The full calculation then runs. |
| `pa_OptimaFastCollapse` | 1000 (on) | The quick first try stops at once when its energy grows to more than this number of times its start value. |
| `pa_OptimaFrozenCut` | 20 (on) | Stops a calculation when its energy and error do not change at all for this number of iterations. The finishing step and the new tries then take over. |
| `pa_OptimaTailStop` | 100 (on) | Stops a calculation when its energy does not change for this number of iterations and the remaining error is only on trace amounts. The solver checks the stop point for phase stability. If the check passes, the solver takes the point as the result. |
| `pa_OptimaTailStopWarm` | 10 (on) | The same, for a warm calculation (SOP, SHP). |
| `pa_OptimaTailDryOnly` | 2 (on) | Sets where the tail stop acts. |
| `pa_OptimaWaterKeep` | 2 (on) | The value 1: the solver keeps a converged result when the new try with more water fails. |
| `pa_OptimaSnapLeaving` | 2 (on) | In a warm restart, a species that goes to zero can limit the step many times in a row. After this number of steps, the solver puts the species on its floor. The species then does not become smaller a little at a time. |
| `pa_OptimaFastFirst` | 60 (on) | A cold-start calculation (AOP) first tries a quick solve that takes full steps. If the quick solve does not end well, the normal calculation runs from the same start. |
| `pa_OptimaFastFrozen` | 3 (on) | The quick first try stops when this number of steps in a row change almost nothing. |
| `pa_OptimaFastDiverge` | 10 (on) | The quick first try stops when its energy jumps far above the start and its error does not improve by this factor. |
| `pa_OptimaColdLevel` | 1 (on) | If the reduced first solve is not used, a cold-start calculation starts from a simple energy-minimising estimate. Two copies of one solid solution start at different compositions. |
| `pa_OptimaLevelStart` | 1 (on) | For systems with only solids and other mixed condensed phases, a cold-start calculation first finds which mixed phases form and what their compositions are. It then starts from that result. |
| `pa_OptimaPreSolveSeed` | 1 (on) | The first pass of the reduced first solve starts from the simple energy-minimising estimate. If a pass needs more than 300 steps, the solver drops the reduced first solve for that calculation. |
| `pa_OptimaEmptyPreSolveSkip` | 1 (on) | The reduced first solve can find no usable set of species in a calculation. Then the new tries of that calculation do not use it. |
| `pa_OptimaCapNudges` | 3 (used only with `pa_OptimaCallCap`) | Sets how many times the solver repeats a calculation that hit its cap. Each repeat changes the amounts by a very small amount. |
| `pa_OptimaFinishFirst` | 1 (on) | After a stuck calculation fails, the solver first runs the finishing step on the point where the calculation stopped. It keeps that result if no phase tends to form. |
| `pa_OptimaCollapseCheck` | 1 (on) | When the solver checks a result, it judges a phase that disappeared by its pure species. A phase that must come back is then not accepted as correctly absent. |
| `pa_OptimaTpdVeto` | 1 (on) | The direct check that accepts a result that is not fully converged (`pa_OptimaTpdAccept`) cannot overrule the phase stability check. This is for a mixed phase that the direct check cannot test, for example a water phase that disappeared. |
| `pa_OptimaTwinReadmit` | 1 (on) | The solver can switch off a mixed phase because a second copy of it took over. After the calculation, the solver tries that phase again if the phase stability check says that it must form. |
| `pa_OptimaFloorCurvature` | 3 (on) | Uses the exact curvature of the logarithm term for each species in a mixed phase. The numerical estimate is too low for very small amounts. |
| `pa_OptimaEmptySolvent` | 1 (on) | If the water in a cold-start calculation is at or below its minimum amount, the solver adds water. This is the case when the dissolved amounts are larger than the water. |
| `pa_OptimaDryRetry` | 0 (off) | A calculation can end with almost no liquid water while the bulk can form much more. Then the solver solves again from a start that first puts hydrogen and oxygen into water. |
| `pa_OptimaFastAcceptCut` | 1 (on) | The quick first try (`pa_OptimaFastFirst`) can use all its iterations and still end OK. |
| `pa_OptimaLevelFinish` | 1 (on) | The leveled start of a cold-start calculation can hold mixed phases. Then the finishing step first runs on that start. |
| `pa_OptimaFinishStall` | 15 (on) | The finishing step stops after this number of steps without progress. |
| `pa_OptimaLSWarmTrials` | 5 (on) | On a warm calculation (SOP, SHP), the line search tries at most this number of step lengths. |
| `pa_OptimaPartialPivLU` | 100 (on) | Optima uses a faster method for its linear equations when there are at least this number of equations. |
| `pa_OptimaFDPureSkip` | 1 (on) | The solver sets the curvature of a pure phase to exactly zero. It does not do an activity calculation. |
| `pa_OptimaTraceRespread` | 1 (on) | After the calculation, an element that only the floor amounts of its species hold gets its real amount back. |
| `pa_TpdStrongBasins` | 4 (on) | Both solvers: the thorough check of whether a mixed phase can form tries each start composition one time. |
| `pa_LogBarrierTau` | 1e-16 | A very small numerical push that keeps pure minerals away from exactly zero. |
| `pa_ColdRetryNudges` | 4 (on) | If a cold-start calculation (AIA) fails, the solver tries again up to this number of times. Each try uses a composition that is changed by a very small amount. The solver then finishes at the exact composition. |
| `pa_IpmStallWindow` | 30 (on) | Stops when the energy and the composition do not change for this number of iterations. |
| `pa_IpmAugmentedKKT` | 2 | Sets how the solver solves the equations of each step. |
| `pa_MbReproject` | 1 (on) | A final touch-up that makes the element totals agree exactly. |
| `pa_PSTALL` | 1 (on) | Lets the mass-balance step stop early when it does not improve. |
| `pa_MbClassRule` | 0 (off) | Checks the totals of major elements in absolute terms. It checks the totals of trace elements in relative terms. The value is the trace/major ratio. |
| `pa_FilloutBudget` | 0 (off) | Limits how much the start guess can change the element totals. The limit is a fraction of the amount of each element. |
| `pa_StabTPD` | 1 (on) | Reports, in the diagnostic trace only, if a mixed phase that is not in the result must be there. |
| `pa_IpmNoProgress` | 500 (on) | Needs a build with Optima. The GEMS3K IPM solver stops if the energy does not improve for this number of iterations. |
| `pa_IpmLevelStart` | 2 (on) | A cold-start calculation (AIA) on a system without aqueous, gas or sorption phases starts from the leveled solution. This is a linear calculation with the mixed phases at their best compositions. |
| `pa_IpmQrSortPivot` | 100 (on) | A step with at least this number of active species uses a faster and more accurate method for its equations. |
| `pa_RescueWalkSteps` | 8 | Sets the number of equal steps in the walk from the previous result to the new composition. |
| `pa_RescueBudget` | 100000 (on) | The rescue (`pa_FailRescue`) skips steps 1 to 4 after the call has used this number of iterations. The count includes the first attempt and both solvers. Step 5, the other solver, still runs one time. 0 = no limit. |
| `pa_OptimaStallNetFactor` | 5 (on) | After a stall, Optima solves one more time with the stall check off. This last solve gets at most this number times the iterations of the solve that stalled, and at least 2000. 0 = the normal limit. |
| `pa_NativeOptimaFinish` | 500 (on) | Needs a build with Optima. A solid can be left out or be in a trace amount and still tend to form. Then Optima finishes the result. The solver keeps the new result only if it is better. The calculation can stop at its iteration limit, or it can fail to solve its equations. Then Optima finishes it from the point where it stopped. |

## Removed settings

Project files can still contain these names, so existing files load without change. The names have no effect:
`pa_LpDualFillout`, `pa_OptimaLSWindow`, `pa_IpmLoopTweaks`, `pa_OptimaFDDiagFloor`, `pa_MbPivotSplit`,
`pa_OptimaPhaseCompaction`, `pa_OptimaReadmitSeed`, `pa_MbTrendPhaseDecay`, `pa_OptimaMaxStepRatio`,
`pa_GAR`, `pa_GAH`.
