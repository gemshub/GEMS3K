# GEMS3K 5.0.0 — solver paths

What happens, step by step, when you call `GEM_run()` in each calculation mode: which solver runs first,
what it tries when something goes wrong, and when it gives up. A drawn version of the same page is
[SOLVER-PATH-GRAPH.html](SOLVER-PATH-GRAPH.html). How to pick a mode and what each setting does is in the
[options guide](SOLVER-OPTIONS-GUIDE.md).

## The short version

- You choose the mode for each call by putting a value in `NodeStatusCH`. GEMS3K never changes your choice
  on its own, except to fall back to the GEMS3K IPM solver in a build without Optima.
- Every mode tries hard before it reports a failure. Each extra try starts from the same inputs, and if it
  does not succeed, you get the **first** failed result back unchanged, not a worse one.
- Extra tries are counted. The iteration count you read after the call includes every try.
- With kinetics switched on (a time step greater than zero), the extra tries described in
  [After a failed call](#after-a-failed-call-the-retry-ladder) are skipped. The tries **inside** a solver
  still run.

## The modes at a glance

| Mode | First | Then | If it fails |
|---|---|---|---|
| **AIA** (1) | GEMS3K IPM solver from scratch | — | Re-solves from a few microscopically changed compositions, then returns to the exact one |
| **SIA** (5) | GEMS3K IPM solver from the previous result | — | Starts again from scratch inside the same call (reported as AIA) |
| **AOP** (10) | Optima from scratch | — | One more try with a plainer line search; then one try from a better first guess |
| **SOP** (14) | Optima from the previous result | — | One more try with a plainer line search; then a full AOP from scratch |
| **HOP** (22) | GEMS3K IPM solver from scratch | Optima, starting from that answer | Optima fails: the GEMS3K IPM answer is kept and flagged. GEMS3K IPM fails: Optima runs alone from scratch |
| **SHP** (26) | GEMS3K IPM solver from the previous result | Optima, starting from that answer | As HOP; and if the warm GEMS3K IPM start fails, it reruns from scratch first |

A sixth Optima mode, ROP (18), uses the Optima library's own default choices. It is kept only to compare
results during development and is not meant for normal use.

## What every call does first

Whatever the mode, `GEM_run()`:

1. Reads the temperature, pressure and bulk composition from the node, and checks that the temperature and
   pressure lie inside the data tables.
2. Loads the previous result (warm modes SIA, SOP). SHP uses a usable previous result for its GEMS3K IPM leg; if none is available or the warm attempt fails, it starts from scratch. All other modes prepare a fresh start.
3. If kinetics is switched on, runs one kinetics time step.
4. Holds back species that are exact copies of other species, for this call only.
5. If `pa_DG` is above 1e-5, scales all amounts to a fixed total so the numbers stay well sized, and scales
   them back at the end.

Then it runs the solver for the chosen mode.

## Inside the GEMS3K IPM solver (AIA, SIA, and the first part of HOP and SHP)

1. **First guess.** From scratch, a simple linear calculation finds a starting set of amounts. From the
   previous result, it reuses that result. If the system has only pure phases, this step already gives the
   answer.
2. **Mass balance fix.** Adjusts the amounts so that every element adds up to the bulk composition.
3. **Main descent.** Lowers the Gibbs energy step by step, keeping the mass balance, until the steps become
   small enough.
4. **Phase check** (with `pa_PC` = 2, the default). Removes phases that should not be there, adds back
   phases that should, and cleans up tiny amounts. It never adds a phase in a larger amount than the bulk
   composition can supply. If it changed the phases, steps 2–4 run again, up to a fixed number of passes.
5. **Second mass balance fix**, on the final set of phases.
6. **Final checks.** If an element still does not add up, it tries one repair (`pa_MbReproject`, on by
   default). If the repair does not help, the answer is still returned, with a warning in the log.

**Built-in restarts.** If a warm start (SIA) fails anywhere in these steps, the solver starts again from
scratch in the same call, and the status says AIA. If a start from scratch fails because the numbers start
to run away, it restarts once more with a simpler phase check before giving up.

## Inside the Optima solver (AOP, SOP, and the second part of HOP and SHP)

1. **Quick first try** (AOP only, `pa_OptimaFastFirst`, on by default). A short, plain solve that takes full
   steps. It gives up early if its energy keeps rising or blows up, if its steps stop changing anything, or if
   its error stops improving after a jump in energy. If it does not end OK, the normal solve below starts from
   the same inputs, and the iterations of both are counted.
2. **First guess.**
   - From scratch: a linear calculation that puts the mass into as few species as possible. If that leaves
     almost no water, or empties a mixed phase, the guess is corrected first. On a second try (see below) a
     better guess is used that also looks at mixed phases.
   - Leveled start (`pa_OptimaColdLevel`, `pa_OptimaLevelStart`, on by default; the quick try and the first
     full attempt only): the first guess also takes the mixed phases into account. In a system with only
     solids and other mixed condensed phases, the mixed phases that form, and their compositions, are worked
     out first; if none forms, Optima only confirms that answer. Where the smaller problem (next step) is not
     used, the start is the linear-calculation answer and its element potentials, with two copies of one solid
     solution starting at different compositions.
   - From the previous result: that result is reused. If the node has never been solved, a warning is
     written and the call starts from scratch instead.
   - In HOP and SHP: the GEMS3K IPM answer is the starting point.
3. **Smaller problem first** (systems with 200 species or more, AOP and SOP). Solves a smaller problem with
   only the species likely to matter, then adds back any species that turns out to be needed. If this does
   not settle, it is dropped and the full problem is solved.
4. **Main solve.** Optima lowers the Gibbs energy. It stops early if it has made no progress for 200 steps
   (`pa_OptimaStallWindow`), so a retry can start sooner. It also stops, so the finishing step and the retries
   can take over, when its energy and error have not changed at all for 20 steps (`pa_OptimaFrozenCut`; not
   in HOP and SHP). On systems with a multisite solid solution it
   takes a short first look (200 steps from scratch, 25 from a previous result) to spot a phase that should
   disappear.
5. **Tries inside the call**, in this order, only while the answer is not yet good:
   1. Put the water back, if the result left almost none.
   2. Put back a mixed phase that collapsed to nothing.
   3. If two identical mixed phases compete, remove the one that is disappearing.
   4. Check every phase: remove one that should not be there, or bring back one that should, and solve
      again. Repeated until nothing changes.
   5. Last resorts: solve once more without the early stop, or without the short first look.
   6. After a stuck solve has failed, the finishing step is tried first on where it stopped, and its answer
      is kept if no phase wants to form; the full re-solves run only otherwise (`pa_OptimaFinishFirst`).
6. **Finish.** If Optima stopped just short, a final step finishes the job on the phases it found
   (`pa_OptimaFinish`, on by default).
7. **Final checks.** Does every element add up? Should a missing phase form? Were the pH and Eh targets
   met? An answer that fails a check is reported as not fully trustworthy.
8. **Accept by direct check.** If Optima's own stopping test was not met but the mass balance holds and a
   direct test shows that no missing phase should form, the answer is accepted (`pa_OptimaTpdAccept`, on by
   default).
9. **Touch-up.** If an element does not quite add up in an accepted answer, a small repair is tried and
   kept only if it helps (`pa_OptimaZeroAbsent` = 2, the default). With the value 1, species of absent
   phases are instead set to exactly zero.

**Cap on one call** (AOP and SOP, `pa_OptimaCallCap`, off by default). The whole call, with all its tries,
may use only a set number of iterations. A call that fails with the cap used up is repeated with the bulk
composition changed by a tiny amount, up to `pa_OptimaCapNudges` times, each with a fresh cap; an answer found
this way is finished at the exact composition. If none works, the call is run once more without the cap, so the
cap never loses an answer.

**Dry-out path** (optional, `pa_OptimaDryOut`, off by default). When the water runs out during the main
solve, the solve stops, the water phase is taken out and the solids are solved alone. That answer is then
compared with a solve that keeps the water, and the one with the lower energy is kept. If the dry solve
fails, the ordinary solve goes on with the water back. It is not used when an element can only be held by the
water. Two further settings act after the water has been taken out: the tail stop (`pa_OptimaTailStop`,
`pa_OptimaTailStopWarm`, `pa_OptimaTailDryOnly`) ends a solve whose energy no longer changes while the
remaining error sits only on trace amounts, after a phase stability check; and `pa_OptimaWaterKeep` keeps a
converged or tail-stopped dry answer instead of retrying it with more water.

Only the Optima modes can fix the pH or Eh (see the options guide).

## The hybrid modes in detail (HOP, SHP)

1. The GEMS3K IPM solver runs: from scratch in HOP, from the previous result in SHP. If SHP's node has never
   been solved, it runs from scratch and writes a warning.
2. Its answer is saved.
3. Optima runs, starting from that answer.
4. The outcome:
   - Optima succeeds: you get Optima's answer.
   - Optima fails: the saved GEMS3K IPM answer is put back, and the status is BAD (for example
     `BAD_GEM_HOP`). So HOP is never worse than the GEMS3K IPM solver alone.
   - The GEMS3K IPM solver failed in step 1: in SHP it first reruns from scratch. If there is still no
     answer, Optima runs alone from scratch, as in AOP.

The iteration count of a hybrid call is the sum of both parts. `GEM_IterationsHOP()` gives the two parts
separately.

## After a failed call: the retry ladder

When the solver returns a failure, `GEM_run()` can try again from the same inputs. The steps run in this
order and stop at the first success. A step that does not succeed hands back the result as it was before
that step.

| Step | Applies to | What it does | Setting (default) |
|---|---|---|---|
| 1 | AOP, SOP, HOP, SHP | Solve once more with the line search's anti-stall step switched off | `pa_OptimaLineSearch` (1.5), `pa_OptimaLSStallEscape` (10) |
| 2 | AIA only, after a failure with no result | Solve from up to 4 microscopically changed compositions; the first that works is used as a warm start for the exact composition | `pa_ColdRetryNudges` (4) |
| 3 | SOP, SHP | Solve again from scratch as AOP. Reported as AOP if it succeeds | `pa_OptimaColdRetry` (2) |
| 4 | AOP, and step 3 | Solve from scratch once more, from the better first guess | `pa_OptimaCgSeed` (1e-6) |

With `pa_OptimaColdRetry` = 2 (the default), a warm Optima call that is going badly gives up early and goes
straight to step 3, which is usually cheaper than finishing the warm try.

## Reading the result

- **OK** (`OK_GEM_...`): the answer passed all checks.
- **BAD** (`BAD_GEM_...`): there is an answer, but it did not pass every check. Look at the log. In HOP and
  SHP this usually means Optima failed and you have the GEMS3K IPM answer.
- **ERR** (`ERR_GEM_...`): no answer. All values read from the node are left over from before the call.
- The status can name a different mode than the one you asked for: SIA → AIA after a restart from scratch,
  SOP or SHP → AOP after step 3 of the ladder.

## Which mode to choose

The full table is in the [options guide](SOLVER-OPTIONS-GUIDE.md#which-mode-to-use). In short:

- **One calculation:** AIA. It is the fastest by far on almost every system.
- **A checked answer:** HOP. It costs little more than AIA and is never worse.
- **Solid solutions, melts, phase diagrams:** AOP or HOP.
- **A series of small steps** (sweeps, titrations, transport): start with AIA or AOP, then use SIA or SOP
  (or SHP) for the following points. Starting from the previous result is much cheaper.
- **A series that crosses phase boundaries:** prefer SOP or SHP over SIA. SIA can keep a phase after it
  stops being stable.
- **Very little water:** AIA or HOP.

## Without Optima

In a build without Optima, AOP and ROP run as AIA, SOP as SIA, HOP as AIA and SHP as SIA, with a warning in
the log. The status then names the GEMS3K IPM mode.
