# GEMS3K 5.2.0 — glossary

These words have one meaning in the [options guide](SOLVER-OPTIONS-GUIDE.md), the
[solver paths](SOLVER-PATH-GRAPH.md) and the [advanced settings](SOLVER-ADVANCED-SETTINGS.md).

| Word | Meaning |
|---|---|
| Activity coefficient | A number that corrects the amount of a species for non-ideal behaviour. |
| Aqueous phase | The water solution. |
| Bulk composition | The total amount of each element in the system. |
| Cold start | A start that does not use a previous result. AIA, AOP and HOP always make a cold start. |
| Control condition, target | A pH or Eh that you ask for. The solver adds the acid, base or electron donor that it needs. |
| Curvature | How fast the Gibbs energy changes its slope when an amount changes. The solver uses it to choose a step. |
| Eh, pe, pH | Measures of the redox level and of the acidity of the aqueous phase. |
| Element potential | The change of the Gibbs energy when you add one unit of an element. Also called the dual. |
| Finish, finishing step | A final calculation on a fixed set of phases. It starts from the point where the solver stopped. |
| Floor | The smallest amount that Optima gives to a species that is absent. |
| GEMS3K IPM solver | The solver that GEMS3K has always had: an interior-point method, then a correction of the mass balance. Modes AIA and SIA. |
| Gibbs energy | The quantity that the solver makes as small as it can. |
| Iteration | One step of a solver. |
| Leveled solution | The result of a linear calculation where the mixed phases are at their best compositions. |
| Leveled start | A cold start from the leveled solution. |
| Line search | A method that tries shorter steps when a full step makes the result worse. |
| Mass balance | The condition that the amount of each element in the result agrees with the bulk composition. |
| Miscibility gap | A range of composition where a mixed phase splits into two phases. |
| Mixed phase | A phase with more than one species whose amounts change its properties. Examples are a solid solution and a melt. |
| Native | Used only in setting names, for example `pa_NativeOptimaFinish`. It means the GEMS3K IPM solver. |
| Node | The object in GEMS3K that holds one system and its result. |
| Optima solver | The second solver, based on the Optima library. Modes AOP and SOP. |
| Phase check | A test that finds which phases must be added to the result or removed from it. |
| Quick first try | A short Optima solve with full steps that a cold start (AOP) tries first. If it does not end well, the normal solve runs. |
| Rescue | The steps that GEMS3K tries when a call still fails after the retry ladder. |
| Retry ladder | A fixed list of new tries that GEMS3K runs, in order, after a call fails. |
| Solid solution | A solid phase whose end-members mix in any ratio. |
| Sorption | The binding of species to a solid surface. |
| Stall | A state where a solver runs but its error or energy does not improve. |
| Status | The value that `GEM_run()` returns. It starts with OK, BAD or ERR. |
| Tail stop | A stop of a solve whose energy does not change and whose error is only on trace amounts. |
| Trace amount | A very small amount of a species or an element. |
| Twin phases | Two phases that have the same end-members. |
| Warm start | A start from the previous result in the node. SIA, SOP and SHP make a warm start. |
