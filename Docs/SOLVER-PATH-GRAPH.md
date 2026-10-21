# GEMS3K 5.2.0 — solver paths

This page shows what happens when you call `GEM_run()` in each calculation mode. It shows which solver runs
first. It shows what the solver tries when a calculation fails. It shows when the solver stops.

The [Glossary](GLOSSARY.md) explains the special words. A drawn version of this page is in [SOLVER-PATH-GRAPH.html](SOLVER-PATH-GRAPH.html). To choose a mode and to
see what each setting does, read the [options guide](SOLVER-OPTIONS-GUIDE.md). Settings for tests are in
[Advanced settings](SOLVER-ADVANCED-SETTINGS.md).

## Summary

- You choose the mode for each call. Put a value in `NodeStatusCH`. GEMS3K does not change your choice.
  One exception: in a build without Optima, GEMS3K uses the GEMS3K IPM solver.
- Each mode tries many times before it reports a failure. Each new try starts from the same inputs.
- If no try gives a good result, you get the best attempt that keeps the element totals. Its status is BAD.
  See [The rescue](#the-rescue-pa_failrescue).
- With `pa_FailRescue` = 0, you get the **first** failed result without change.
- GEMS3K counts all tries. The iteration count that you read after the call includes every try.
- If kinetics is on (the time step is more than zero), GEMS3K does not use the retry ladder in
  [After a failed call](#after-a-failed-call-the-retry-ladder). The tries **inside** a solver and the rescue
  still run.

## The modes

| Mode | First | Then | If it fails |
|---|---|---|---|
| **AIA** (1) | GEMS3K IPM solver from scratch | — | The solver tries again from some compositions that are changed by a very small amount. Then it goes back to the exact composition. |
| **SIA** (5) | GEMS3K IPM solver from the previous result | — | The solver starts again from scratch in the same call. The status shows AIA. |
| **AOP** (10) | Optima from scratch | — | The solver tries again with a plainer line search. Then it tries once from a better first guess. |
| **SOP** (14) | Optima from the previous result | — | The solver tries again with a plainer line search. Then it does a full AOP from scratch. |
| **HOP** (22) | GEMS3K IPM solver from scratch | Optima, from the result of the first part | If Optima fails, GEMS3K keeps the GEMS3K IPM result and marks it. If the GEMS3K IPM solver fails, Optima runs alone from scratch. |
| **SHP** (26) | GEMS3K IPM solver from the previous result | Optima, from the result of the first part | The same as HOP. If the warm GEMS3K IPM start fails, the solver first starts again from scratch. |

After these tries, each mode can go to the rescue. The rescue is on by default (`pa_FailRescue`).

ROP (18) is a sixth Optima mode. It uses the default choices of the Optima library. It is only for the
comparison of results during development. Do not use it for normal work.

## What each call does first

In all modes, `GEM_run()` does these steps:

1. It reads the temperature, the pressure and the bulk composition from the node.
2. It checks that the temperature and the pressure are in the range of the data tables.
3. It loads the previous result in the warm modes (SIA and SOP).
4. SHP uses the previous result for its GEMS3K IPM part if a usable result exists. If no usable result exists,
   or if the warm attempt fails, SHP makes a cold start.
5. All other modes prepare a new start.
6. If kinetics is on, it does one kinetics time step.
7. It holds back, for this call only, each species that is an exact copy of another species.
8. If `pa_DG` is more than 1e-5, it scales all amounts to a fixed total. This keeps the numbers at a good
   size. At the end, it scales the amounts back.

Then it runs the solver for the chosen mode.

## Inside the GEMS3K IPM solver

The GEMS3K IPM solver runs in AIA and SIA. It also runs as the first part of HOP and SHP.

1. **First guess.**
   - In a cold start, a simple linear calculation finds a set of start amounts.
   - In a warm start, the solver uses the previous result.
   - If the system has only pure phases, this step gives the final result.
   - If the system has no aqueous, gas or sorption phases, the solver starts from the leveled solution. This
     is the linear calculation with the mixed phases at their best compositions (`pa_IpmLevelStart`).
     If this call fails, GEMS3K solves it again from the usual start.
2. **Mass balance fix.** The solver changes the amounts until each element agrees with the bulk composition.
3. **Main descent.** The solver lowers the Gibbs energy step by step and keeps the mass balance. It stops
   when the steps are small enough.
   - If the energy does not improve for 500 iterations (`pa_IpmNoProgress`), the descent stops. Optima then
     finishes the result from that point.
   - If Optima fails, GEMS3K solves the call again without this stop.
4. **Phase check.** This step runs when `pa_PC` = 2 (the default). It removes phases that must not be there.
   It adds phases that must be there. It removes very small amounts.
   - It never adds a phase in a larger amount than the bulk composition can supply.
   - If the phase check changes the phases, steps 2 to 4 run again. The number of passes has a limit.
5. **Second mass balance fix.** This step uses the final set of phases.
6. **Final checks.** If an element does not agree with the bulk composition, the solver does one repair
   (`pa_MbReproject`, on by default). If the repair does not help, the solver gives the result and writes a
   warning in the log.
7. **Optima finish.** This step needs a build with Optima (`pa_NativeOptimaFinish`).
   - A solid can be left out of the result but still tend to form. Then Optima finishes the result.
     GEMS3K keeps the new result only if it is better.
   - The descent can stop at its iteration limit. It can also stop when it cannot solve its equations.
     Then Optima finishes the result from that point.

**Restarts in the solver.**
- If a warm start (SIA) fails in one of these steps, the solver makes a new cold start in the same call.
  The status shows AIA.
- If a cold start fails because the numbers grow without limit, the solver starts again one more
  time. It uses a simpler phase check. Then it stops.

## Inside the Optima solver

The Optima solver runs in AOP and SOP. It also runs as the second part of HOP and SHP.

1. **Quick first try.** This step is only in AOP (`pa_OptimaFastFirst`, on by default). It is a short, plain
   solve that takes full steps.
   - It stops early if the energy keeps going up or becomes very large.
   - It stops early if its steps do not change anything.
   - It stops early if its error does not improve after a jump in energy.
   - If it does not end OK, the normal solve below starts from the same inputs. GEMS3K counts the iterations
     of both.
   - If it uses all its iterations but ends OK, GEMS3K keeps the result when no missing phase can form
     (`pa_OptimaFastAcceptCut`).
2. **First guess.**
   - Cold start: a linear calculation puts the mass in as few species as it can. If this leaves almost no
     water, or leaves a mixed phase empty, GEMS3K corrects the guess first. On a second try (see below), GEMS3K
     uses a better guess that also checks the mixed phases.
   - Leveled start (`pa_OptimaColdLevel`, `pa_OptimaLevelStart`, on by default). It is only for the quick try
     and the first full attempt. The first guess also uses the mixed phases.
     - If the system has only solids and other mixed condensed phases, GEMS3K first finds which mixed phases
       form and what their compositions are. If no mixed phase forms, Optima only confirms this result.
     - If mixed phases form, the finishing step first runs on that start. Optima then starts from the new
       point (`pa_OptimaLevelFinish`).
     - If the smaller problem (step 3) is not used, the start is the result of the linear calculation and its
       element potentials. Two copies of one solid solution start at different compositions.
   - Warm start: GEMS3K uses the previous result. If the node has never been solved, GEMS3K writes a
     warning and makes a cold start.
   - In HOP and SHP: the start is the result of the GEMS3K IPM solver.
3. **Smaller problem first.** This step is for systems with 200 species or more, in AOP and SOP. The solver
   solves a smaller problem with only the species that are likely to matter. Then it adds each species that
   it needs. If this step does not settle, GEMS3K drops it and solves the full problem.
4. **Main solve.** Optima lowers the Gibbs energy. It stops in these cases:
   - It stops if it makes no progress for 200 steps (`pa_OptimaStallWindow`). A new try can then start sooner.
   - It stops if its energy and its error do not change at all for 20 steps (`pa_OptimaFrozenCut`). The
     finishing step and the new tries then take over. This stop is not in HOP and SHP.
   - On a warm call, the line search tries at most 5 step lengths (`pa_OptimaLSWarmTrials`).
   - On a warm call, an error can stay only on a solid solution that is almost gone. Then the solve stops and
     holds that solid at zero. The solver then finishes. The hold stays only if no absent solid tends to form.
     If a solid does tend to form, the solve continues from the point where it stopped
     (`pa_OptimaTailDryOnly` = 2).
   - For a system with a multisite solid solution, the solver first looks for a phase that must disappear.
     This look is short: 200 steps in a cold start, 25 steps in a warm start.
5. **Tries in the call.** These tries run in this order. They run only while the result is not good.
   1. Put the water back, if the result has almost no water.
   2. Put back a mixed phase that became empty.
   3. If two identical mixed phases compete, remove the one that disappears.
   4. Check each phase. Remove a phase that must not be there. Bring back a phase that must be there. Solve
      again. Repeat until nothing changes.
   5. As a last step, solve again without the early stop, or without the short first look.
   6. After a stuck solve fails, run the finishing step on the point where the solve stopped. Keep its result
      if no phase tends to form. Run the full new solves only if the result is not kept
      (`pa_OptimaFinishFirst`).
6. **Finish.** If Optima stops just short of the result, a final step finishes the job on the phases that
   Optima found (`pa_OptimaFinish`, on by default). The finishing step stops early if its error, its energy
   and its step size do not change for 15 steps (`pa_OptimaFinishStall`).
7. **Final checks.** The solver checks these points:
   - Does each element agree with the bulk composition?
   - Does a missing phase need to form?
   - Does the result meet the pH target and the Eh target?

   If a check fails, GEMS3K reports the result as not fully trustworthy. Before the checks, an element that
   only the floor amounts of its species hold gets its real amount back (`pa_OptimaTraceRespread`).
8. **Accept by direct check.** Optima can fail its own stopping test. If the mass balance holds, and a direct
   test shows that no missing phase must form, GEMS3K accepts the result (`pa_OptimaTpdAccept`, on by default).
9. **Touch-up.** An accepted result can have a small error in the element totals. The solver tries a small
   repair. It keeps the repair only if the repair helps (`pa_OptimaZeroAbsent` = 2, the default). With the
   value 1, the species of absent phases are set to exactly zero.

**Cap on one call.** The cap is for AOP and SOP (`pa_OptimaCallCap`, off by default). The whole call, with all
its tries, can use only a set number of iterations.
- If a call fails and the cap is used up, GEMS3K repeats the call. It changes the bulk composition by a very
  small amount. It does this up to `pa_OptimaCapNudges` times. Each repeat has a new cap.
- GEMS3K finishes a result from a repeat at the exact composition.
- If no repeat works, GEMS3K runs the call one more time without the cap. For this reason, the cap does not
  lose a result.

**Dry-out path.** This path is optional (`pa_OptimaDryOut`, off by default). It works as follows:
- The water can run out in the main solve. Then the solve stops and GEMS3K takes the water phase out.
- The solver solves the solids alone.
- GEMS3K compares this result with a solve that keeps the water. It keeps the result with the lower energy.
- If the dry solve fails, the usual solve continues with the water back.
- GEMS3K does not use this path if an element can only be in the water.

Three settings act after the water is out:
- The tail stop (`pa_OptimaTailStop`, `pa_OptimaTailStopWarm`, `pa_OptimaTailDryOnly`) ends a solve whose
  energy does not change and whose remaining error is only on trace amounts. A phase stability check follows.
- `pa_OptimaWaterKeep` keeps a dry result that converged or that the tail stop ended. GEMS3K does not try it
  again with more water.
- With `pa_OptimaTailDryOnly` = 2 (the default), the tail stop also acts on warm calls. Step 4 shows this.

Only the Optima modes can fix the pH or the Eh. See the options guide.

## The hybrid modes (HOP and SHP)

1. The GEMS3K IPM solver runs. In HOP it makes a cold start. In SHP it makes a warm start.
   - If the node of SHP has never been solved, the solver makes a cold start and writes a warning.
   - If the solver stops because the energy does not improve (`pa_IpmNoProgress`), Optima starts from that
     point. It first starts warm. If that fails, it makes a cold start.
   - If the solver starts from the leveled solution and fails, GEMS3K solves it again from the usual start.
2. GEMS3K saves the result of this part.
3. Optima runs and starts from the saved result.
4. The outcome is one of these:
   - Optima succeeds. You get the result of Optima.
   - Optima fails. GEMS3K puts back the saved GEMS3K IPM result. The status is BAD, for example
     `BAD_GEM_HOP`. For this reason, HOP is never worse than the GEMS3K IPM solver alone.
   - The GEMS3K IPM solver failed in step 1. In SHP, the solver first runs again with a cold start. If there is
     still no result, Optima runs alone with a cold start, as in AOP.

The iteration count of a hybrid call is the sum of the two parts. `GEM_IterationsHOP()` gives the two parts
separately.

## After a failed call: the retry ladder

When a solver returns a failure, `GEM_run()` can try again from the same inputs. The steps run in the order
of the table. The ladder stops at the first success. If a step does not succeed, it gives back the result that
existed before that step.

| Step | Applies to | What it does | Setting (default) |
|---|---|---|---|
| 1 | AOP, SOP, HOP, SHP | Solves again with the anti-stall step of the line search switched off. | `pa_OptimaLineSearch` (1.5), `pa_OptimaLSStallEscape` (10) |
| 2 | AIA only, after a failure with no result | Solves from up to 4 compositions that are changed by a very small amount. The first that works is the warm start for the exact composition. | `pa_ColdRetryNudges` (4) |
| 3 | SOP, SHP | Solves again with a cold start as AOP. If it succeeds, the status shows AOP. | `pa_OptimaColdRetry` (2) |
| 4 | AOP, and step 3 | Makes one more cold start, from the better first guess. | `pa_OptimaCgSeed` (1e-6) |

If `pa_OptimaColdRetry` = 2 (the default), a warm Optima call that goes badly stops early. It goes straight
to step 3. This is usually less expensive than the end of the warm try.

## The rescue (`pa_FailRescue`)

A call can still end with an error or a BAD result. Then the rescue tries the call again in steps. It stops at
the first step that succeeds. GEMS3K finishes each successful step at the exact composition with the normal
settings.

| Step | What it does |
|---|---|
| 1 | Solves warm from the previous result of the node. |
| 2 | Solves at compositions that are changed by a very small amount. |
| 3 | Goes in 8 steps from the previous result to the new composition (`pa_RescueWalkSteps`). |
| 4 | Solves with looser settings. |
| 5 | Solves with the other solver. |
| 6 | Gives the best attempt that keeps the element totals, with a BAD status (only if `pa_FailRescue` = 2). |

- If a kinetic call has no result from any step, the rescue skips steps 4 and 5.
- If there is a pH target or an Eh target, GEMS3K keeps only an attempt that meets the target.
- `GEM_rescue_step()` tells you which step gave the result. 0 means the call itself.

## How to read the result

- **OK** (`OK_GEM_...`): the result passed all checks.
- **BAD** (`BAD_GEM_...`): there is a result, but it did not pass all checks. Read the log. In HOP and SHP,
  this usually means that Optima failed and you have the GEMS3K IPM result.
- **ERR** (`ERR_GEM_...`): there is no result. All values that you read from the node are from before the call.
- The status can name a different mode from the mode that you asked for.
  - SIA changes to AIA after a new cold start.
  - SOP and SHP change to AOP after step 3 of the ladder.
- After a rescue, `GEM_rescue_step()` is not 0. If the value is 6, the result is the best attempt and the
  status is BAD.

## How to choose a mode

The full table is in the [options guide](SOLVER-OPTIONS-GUIDE.md#which-mode-to-use). In short:

- **One calculation:** use AIA. It is the fastest on almost all systems.
- **A checked result:** use HOP. It costs a little more than AIA and is never worse.
- **Solid solutions, melts and phase diagrams:** use AOP or HOP.
- **A series of small steps** (sweeps, titrations, transport): use AIA or AOP for the first point. Use SIA,
  SOP or SHP for the next points. A warm start costs much less.
- **A series that crosses phase boundaries:** use SOP or SHP, not SIA. SIA can keep a phase after the phase
  is not stable.
- **Very little water:** use AIA or HOP.

## Without Optima

In a build without Optima, these modes change:

| Mode you ask for | Mode that runs |
|---|---|
| AOP, ROP, HOP | AIA |
| SOP, SHP | SIA |

GEMS3K writes a warning in the log. The status then names the GEMS3K IPM mode.
