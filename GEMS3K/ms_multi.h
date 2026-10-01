//-------------------------------------------------------------------
// $Id$
//
/// \file ms_multi.h
/// Declaration of TMulti class, configuration, and related functions
/// based on the IPM work data structure MULTI
//
/// \struct MULTI ms_multi.h
/// Contains chemical thermodynamic work data for GEM IPM-3 algorithm
//
// Copyright (c) 1995-2013 S.Dmytriyeva, D.Kulik, T.Wagner
// <GEMS Development Team, mailto:gems2.support@psi.ch>
//
// This file is part of the GEMS3K code for thermodynamic modelling
// by Gibbs energy minimization <http://gems.web.psi.ch/GEMS3K/>
//
// GEMS3K is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation, either version 3 of
// the License, or (at your option) any later version.

// GEMS3K is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with GEMS3K code. If not, see <http://www.gnu.org/licenses/>.
//-------------------------------------------------------------------
//
#ifndef MS_MULTI_BASE_H
#define MS_MULTI_BASE_H

#include "m_const_base.h"
#include "verror.h"
#include "datach.h"
// TSolMod header
#include "s_solmod.h"
// TsorpMod and TKinMet
#include "s_sorpmod.h"
#include "s_kinmet.h"
#include "gems3k_impex.h"
#include "ipm_optima.h"

#include <cstdlib>
#include <cstdio>
#include <vector>

class GemDataStream;
class TProfil;
class TNode;

const int  QPSIZE = 180, // earlier 20, 40 SD oct 2005
           QDSIZE = 60;

// Physical constants - see m_param.cpp or ms_param.cpp
extern const double R_CONSTANT, NA_CONSTANT, F_CONSTANT,
    e_CONSTANT,k_CONSTANT, cal_to_J, C_to_K, lg_to_ln, ln_to_lg, H2O_mol_to_kg, Min_phys_amount;

/// The `{ }` at the end of each field's comment below is the value a STANDALONE GEMS3K run gets
/// when the project file does not pin that field - i.e. the entry in `pa_p_` (`ms_multi_diff.cpp`),
/// which is the only default this library has. Re-derived field by field on 2026-09-21 with
/// `gems-benchmark/tools/basefield_doc_check.py`: 26 of 36 braces named something else - 24 a
/// different NUMBER (`pa_DG` documented 1e5 against a compiled 1000. being the one that was
/// reported), plus `PE` and `IEPS`, which stated a value set and a settable range where every
/// other field states its default, and are now written the same way as the rest. The checker is
/// a registry claim (`inv-basefield-doc-defaults`), so this cannot rot again silently.
/// A brace naming a value no build produces is the same defect shape as a freeze's
/// `# set` line naming a configuration no row was produced at - it is read as the answer to
/// "what runs if I change nothing", and it was wrong. Note most corpus projects pin most of these
/// in their own -ipm.json, so the compiled default is what a project file's ABSENCE selects, not
/// what a typical run uses. Comment only: no default was changed.
struct BASE_PARAM /// Flags and thresholds for numeric modules
{
   short
           PC,   ///< Mode of PhaseSelect() operation ( 0 1 2 ... ) { 2 }
           PD,   ///< abs(PD): Mode of execution of CalculateActivityCoefficients() functions { 2 }.
                 ///< Modes: 0-invoke, 1-at MBR only, 2-every MBR it, every IPM it. 3-not MBR, every IPM it.
                 ///< if PD < 0 then use test qd_real accuracy mode
           PRD,  ///< Since r1583/r409: Disable (0) or activate (-5 or less) the SpeciationCleanup() procedure { -5 }
                 ///< TWO EFFECTS, and only the first is what this line used to say (corrected 2026-09-09c,
                 ///< plan v5 s100.5/s103.4 - the earlier note that PRD=0 "disables nothing" was itself wrong):
                 ///<  1. It gates PSSC's speciation-cleanup loop, `if( CleanupStatus && pa_p->PRD )`
                 ///<     (ipm_chemical.cpp) - so 0 really does disable that half, as stated. CleanupStatus
                 ///<     itself comes from pa_PC (1 at PC==2, 0 at PC>2, ipm_main.cpp), so BOTH must be set.
                 ///<  2. UNDOCUMENTED UNTIL NOW: it also sets PSSC's AmountThreshold = 10^-|PRD|, and the
                 ///<     "at least 4 decades" floor beside it is guarded by noZero(), so it does NOT apply
                 ///<     at PRD = 0. AmountThreshold is then 1.0 MOL rather than 1e-4 or 1e-5.
                 ///<     With the cleanup loop off, that threshold's one remaining reachable use is the
                 ///<     ELIMINATE / ELIMINATE_MBVIOL split: both branches delete the phase, and the
                 ///<     threshold decides only whether MassBalanceViolation is raised. So at PRD = 0 a
                 ///<     phase carrying up to a WHOLE MOLE is removed without flagging the mass balance,
                 ///<     i.e. without forcing the further IPM loop that would repair it; at -5 the same
                 ///<     removal above 1e-5 mol is flagged. Two corpus projects ship PRD = 0 (both CASHNK),
                 ///<     where it is latent - no ELIMINATE record fires on their one state point.
                 ///<     Clamping the threshold at PRD = 0 would move those projects' answers, so it is a
                 ///<     gated change and not a comment fix; this comment is the part that is free.
           PSM,  ///< Level of diagnostic messages: 0- disabled (no ipmlog file); 1- errors; 2- also warnings 3- uDD trace { 1 }
           DP,   ///< Maximum allowed number of iterations in the MassBalanceRefinement() procedure { 130 }
           DW,   ///< Since r1583: Activate (1) or disable (0) error condition when DP was exceeded { 1 }
           DT,   ///< Since r1583/r409: DHB is relative for all (0) or absolute (-6 or less ) cutoff for major ICs { 0 }
           PLLG, ///< IPM tolerance for detecting divergence in dual solution { 30000; 1 to 1000 is the
                 ///< documented working range and 0 disables the detection, but the shipped default is the
                 ///< |PLLG| >= 30000 case, which InteriorPointsMethod() reads as "allow complete tracing" }
           PE,   ///< Flag for using electroneutrality condition in GEM IPM calculations ( 0 or 1 ) { 1 }
           IIM   ///< Maximum allowed number of iterations in the MainIPM_Descent() procedure up to 9999 { 7000 }
           ;
         double DG,   ///< Standart total moles { 1000. }
           DHB,  ///< Maximum allowed relative mass balance residual for Independent Components ( 1e-9 to 1e-15 ) { 1e-13 }
           DS,   ///< Cutoff minimum mole amount of stable Phase present in the IPM primal solution { 1e-20 }
           DK,   ///< IPM-2 convergence threshold for the Dikin criterion { 1e-6 }. NOTE the shipped default
                 ///< sits ON the lower end of the interval this comment used to call the settable range
                 ///< (1e-6 < DK < 1e-4), and plan v5 s74 measured eight projects whose DK is BELOW their own
                 ///< numerical noise floor, where termination becomes a waiting time rather than a test
           DF,   ///< Threshold for the application of the Karpov phase stability criterion: (Fa > DF) for a lost stable phase { 0.01 }
           DFM,  ///< Threshold for Karpov stability criterion f_a for insertion of a phase (Fa < -DFM) for a present unstable phase { 0.01 }
           DFYw, ///< Insertion mole amount for water-solvent { 1e-5 }
           DFYaq,///< Insertion mole amount for aqueous species { 1e-5 }
           DFYid,///< Insertion mole amount for ideal solution components { 1e-5 }
           DFYr, ///< Insertion mole amount for major solution components { 1e-5 }
           DFYh, ///< Insertion mole amount for minor solution components { 1e-5 }
           DFYc, ///< Insertion mole amount for single-component phase { 1e-5 }
           DFYs, ///< Insertion mole amount used in PhaseSelect() for a condensed phase component  { 1e-6 }
           DB,   ///< Minimum amount of Independent Component in the bulk system composition (except charge "Zz") (moles) (1e-17)
           AG,   ///< Smoothing parameter for non-ideal increments to primal chemical potentials between IPM descent iterations { 1. }
           DGC,  ///< Exponent in the sigmoidal smoothing function, or minimal smoothing factor in new functions { 0. }
           GAR,  ///< Initial activity coefficient value for major (M) species in a solution phase before LPP approximation { 1 }
           GAH,  ///< Initial activity coefficient value for minor (J) species in a solution phase before LPP approximation { 1000 }
           GAS,  ///< Since r1583/r409: threshold for primal-dual chem.pot.difference (mol/mol) used in SpeciationCleanup() { 1e-3 }.
                 ///< before: Obsolete IPM-2 balance accuracy control ratio DHBM[i]/b[i], for minor ICs { 1e-3 }
           DNS,  ///< Standard surface density (nm-2) for calculating activity of surface species (12.05)
           XwMin,///< Cutoff mole amount for elimination of water-solvent { 1e-13 }
           ScMin,///< Cutoff mole amount for elimination of solid sorbent { 1e-13 }
           DcMin,///< Cutoff mole amount for elimination of solution- or surface species { 1e-33 }
           PhMin,///< Cutoff mole amount for elimination of  non-electrolyte solution phase with all its components { 1e-20 }
           ICmin,///< Minimal effective ionic strength (molal), below which the activity coefficients for aqueous species are set to 1. { 1e-5 }
           EPS,  ///< Precision criterion of the SolveSimplex() procedure to obtain the AIA ( 1e-6 to 1e-14 ) { 1e-10 }
           IEPS, ///< Convergence parameter of SACT calculation in sorption/surface complexation models { 1e-3; settable 0.01 to 1e-6 }
           DKIN; ///< Tolerance on the amount of DC with two-side metastability constraints  { 1e-10 }
    char *tprn;       ///< internal

    // Enable (1, default) or disable (0) stall detection in MassBalanceRefinement(): with it
    // disabled, a stalled MBR run no longer exits early (iRet reset to 0) but keeps iterating
    // until DP is exhausted, surfacing as a real "Maximum allowed number of MBR iterations
    // exceeded" failure instead. GEMS3K-internal only: deliberately placed after every other
    // field (with a default member initializer) so that none of GEMSGUI's positional
    // BASE_PARAM/SPP_SETTING aggregate initializers or its binary project-file (de)serialization
    // need to know this field exists at all - it always keeps its default of 1 there. Only
    // GEMS3K's own keyword-based ipm-dat I/O (ms_multi_format.cpp, "pa_PSTALL") can override it.
    short PSTALL = 1;

    // GEMS3K+Optima solver (ipm_optima.cpp, CalculateEquilibriumStateOptima(),
    // USE_OPTIMA_SOLVER builds only) tuning, GEMS3K-internal only - same
    // trailing-field/default-member-initializer placement as PSTALL above,
    // for the same reason (no GEMSGUI project-file impact; keyword-only I/O
    // via ms_multi_format.cpp). No existing pa_p field measures either
    // quantity, so unlike IIM/DHB/DW (reused as-is for this solver), these
    // get their own fields rather than borrowing an unrelated one.
    //
    // Optima's own KKT optimality-error convergence tolerance - deliberately
    // NOT the same quantity as pa_DK (GEMS3K's native Dikin-criterion
    // threshold measures a structurally different residual; reusing DK's
    // typical project values verbatim was confirmed to let Optima report a
    // false "converged" well short of the true stationary point). Default
    // matches Optima's own library default (Optima::ConvergenceOptions::
    // tolerance, Optima/ConvergenceOptions.hpp).
    double OptimaTol = 1.0e-8;

    // Logarithmic-barrier penalty weight added to the Optima objective for
    // pure single-species phases (species whose chemical potential does not
    // depend on composition - e.g. simple mineral phases), ported from
    // Reaktoro's own EquilibriumSetup.cpp (updateGibbsEnergy()/
    // updateGradX()/updateHessX() - read from source, not guessed). Default
    // matches Reaktoro's own default (EquilibriumOptions::epsilon *
    // logarithm_barrier_factor = 1e-16*1.0). IMPORTANT, confirmed
    // 2026-08-23: at this default magnitude the term does NOT fix the
    // known failure case (Resources/gems3k/j_Flowline_G_series1_..., where
    // plain Optima converges to a wrong-but-KKT-stationary point) - tested
    // directly, same wrong answer as with the term absent entirely. A much
    // larger tau (1e-2) does reach the right answer on that case but then
    // Optima's own convergence check never cleanly triggers (hits its
    // iteration cap, reports failure despite the iterate being correct);
    // 1e-4 was no better than 1e-16. The mechanism is real (ported
    // faithfully from Reaktoro's source) but Reaktoro's own robustness on
    // this class of problem is evidently NOT explained by this term at its
    // own default value - the actual mechanism is still unidentified (see
    // GEMS3K's CLAUDE.md, 2026-08-23, for the full investigation and what
    // was ruled out). Left present and exposed via ipm-dat because it's a
    // real, sourced piece of Reaktoro's own objective and doesn't regress
    // anything at its default value - not because it's known to fix the
    // open problem.
    double LogBarrierTau = 1.0e-16;

    // Relative trust-region cap on Optima's per-iteration Newton step,
    // passed through to the modified Optima::BacktrackSearchOptions::
    // max_step_ratio (Optima/BacktrackSearch.cpp/.hpp, local checkout
    // /home/dmiron/git/hub/optima, NOT the vendored/conda-packaged Optima -
    // this option does not exist in stock Optima and only takes effect
    // when GEMS3K is built against the modified local checkout, see
    // debug-optima-vs-reaktoro/README.md for the build recipe). Disabled
    // (0., no cap, byte-identical to stock Optima's own BacktrackSearch
    // behavior) by default - kept OFF deliberately, not merely un-tuned.
    // Tried first (2026-08-24) as the fix for the aqueous-solvent-
    // collapses-to-floor failure (GEMS3K's CLAUDE.md, 2026-08-23) and
    // directly disproven: on Resources/gems3k/j_Flowline_G_series1_...,
    // every tested value (2, 3, 5, 10) made the SAME case actively WORSE
    // (a clean-but-wrong convergence turned into an outright KKT-residual
    // blowup, ~1e15-1e16) rather than better - root-cause tracing (this
    // same CLAUDE.md entry) found the real defect was an AIA cold-start
    // seed placing the solvent below its own solutes' total mass, not a
    // step-size/globalization problem this cap could ever have addressed;
    // capping the step size just slowed the same wrong trajectory down
    // without changing its direction. Left in place (Optima source and
    // this field both) as available infrastructure for a genuinely
    // step-size-related failure mode, should one turn up on a different
    // system - re-validate on its own merits before enabling, don't
    // assume the j_Flowline finding above generalizes either way.
    double OptimaMaxStepRatio = 0.0;

    // Eigenvalue floor, as a fraction of the block's own largest |eigenvalue|,
    // applied to the exact (finite-differenced, symmetrised) curvature block of
    // every NON-AQUEOUS multicomponent phase in the Optima solver. 0 disables
    // the exact block entirely, restoring the previous behaviour (analytic
    // ideal-mixing curvature plus FD columns for Optima's basic set only).
    //
    // Why this exists: inside a miscibility gap the true curvature of a
    // solution phase is INDEFINITE - the negative eigenvalue along the
    // unmixing direction is what a spinodal is - so Newton needs a
    // positive-definite model to get a descent direction at all, and the floor
    // is what sets how far it steps along that soft direction. Without the
    // exact block, the ideal-mixing model is O(1/X) stiff where the truth goes
    // to zero and Newton contracts at 1 - H/B -> 1 (measured on the
    // sanidine-albite solvus: rate 0.9424 at 600 C, 0.9946 at 640, 1.0000 at
    // 650, taking FULL Newton steps throughout). With the exact block but no
    // floor, the step is unbounded and annihilates the phase.
    //
    // The default was chosen on that benchmark plus the 25-project suite, and
    // the response is monotone in the floor over most of its range (a larger
    // floor means shorter steps means more iterations) but NOT everywhere - see
    // GEMS3K's CLAUDE.md, 2026-08-26. Same trailing-field placement rules as
    // PSTALL above; NOTE that adding a field here also requires bumping BOTH
    // hardcoded counts in ms_multi_format.cpp (prar/rddar), which is what kept
    // OptimaTol/LogBarrierTau/OptimaMaxStepRatio unreadable until 2026-08-26.
    double PhaseHessianFloor = 0.01;

    /// Stall/freeze limit for the Optima solver (AOP/SOP/ROP): abandon a solve
    /// when, over a WINDOW of this many Newton iterations, the best-so-far
    /// optimality error has not fallen by a meaningful relative amount
    /// (1e-8, pa_OptimaTol's own default). 0 disables it.
    ///
    /// The test is CUMULATIVE over a window, not per-step. It was per-step
    /// ("this many CONSECUTIVE iterations with no improvement at all") until
    /// 2026-09-03, and that had a measured false positive which is the whole
    /// reason the compiled default below is 0 - see the paragraph on 588 C.
    /// With the cumulative test that specific false positive is gone: the same
    /// point now converges at 1075 with the window armed at 500, where the
    /// per-step rule killed it at 705. Plan v5 section 54 has the measurement.
    ///
    /// NOTE this is NOT the same rule as the pre-solve stall watch in
    /// OptimaReducedPreSolve(), which shares this one field but adds a second
    /// clause ("has the live error even MOVED within the window"). The two
    /// watches genuinely need different rules and the projects that decide it
    /// pull in opposite directions: f_TestSUP98's reduced pre-solve converges
    /// through a 661-iteration excursion with best-so-far exactly frozen and
    /// the live error swinging, so only the range clause saves it; f_/j_CASHNK's
    /// full solve has best-so-far exactly frozen for 19 consecutive windows
    /// while the live error runs a perfect limit cycle, so a range clause there
    /// would read the oscillation as movement and cost them their rescue. Both
    /// measured. Do not "unify" the two tests without re-measuring both.
    ///
    /// TRIAGE, NOT A FIX. It converts a hang into a fast, honest failure and
    /// lets the existing retry tiers start sooner; the underlying step-length
    /// failure is untouched (GEMS3K's CLAUDE.md / Docs/gems3k-optima-plan-v5.md
    /// section 15-16, 2026-08-27). Every remaining AOP failure in the whole
    /// benchmark corpus - f_/j_TestPNTDB, f_/j_TestSUP98, and the two largest
    /// gems3k-psina projects - is a FROZEN ITERATE: the objective and the
    /// residual are bit-identical from iteration 0, so those runs burn their
    /// entire budget (40 minutes at 1392 species) and return nothing.
    ///
    /// Keyed on the best-so-far error rather than on the objective or on the
    /// iterate displacement, because both of those were measured and DO NOT
    /// discriminate: f is frozen to 6 significant figures on f_CASHNK (which
    /// progresses) as well as on complex_1 (which does not), and complex_1's
    /// iterate actually moves MORE per window than f_CASHNK's.
    ///
    /// DEFAULT 500 since 2026-09-03. It was OFF from 2026-08-27, per user
    /// direction ("for cases known to have converged don't add limit of
    /// iterations or time"), on the reasoning that a project which converges
    /// today must not acquire a new way to fail. That direction was then
    /// refined: add the limit if it converges AND takes fewer iterations AND is
    /// faster, provided that holds for all cases OR a smart fallback is
    /// possible. Both halves are now answered:
    ///
    ///   - Where the watch fires it IS a strict improvement, by 10x:
    ///     f_/j_CASHNK converge in 1001 iterations / ~0.5 s with it and 10001 /
    ///     ~5 s without, same G/pH/Vs to every digit.
    ///   - It does NOT hold for all cases - it is a measured no-op on every
    ///     other project in all three corpora (Resources/gems3k 0 of 81 rows
    ///     changed, gems3k-fail 0 of 48, gems3k-psina never fires, the 301-point
    ///     solvus sweep identical). The watch simply never fires on a project
    ///     that is not deadlocked.
    ///   - So the fallback is what decides it, and it exists:
    ///     CalculateEquilibriumStateOptima()'s last retry tier re-solves once
    ///     with this window disarmed when the watch fired AND every other tier
    ///     failed. The worst arming can now do is spend one extra solve on a run
    ///     that was already failing. See plan v5 section 56.
    ///
    /// Set it to 0 to disable, or to a smaller value per project - but note
    /// three independent measurements say do not go below 500 (plan v5 41.2c,
    /// 49.5, 54.5), and that a window short enough to false-positive is now
    /// recoverable rather than fatal, not harmless.
    ///
    /// That principle was reached the hard way. A default of 500 looked safe:
    /// the longest no-improvement run among converging projects appeared to be
    /// mid_1's 141 (f_GEOTHERM 77, f_Solvus 40, CSHSnplus 30). But the 301-point
    /// solvus temperature sweep - which the 5-degree-sampled `solvus.aop` CTest
    /// does NOT cover - has a point at 588 C that stagnates for 500-705
    /// iterations and then RECOVERS, converging at 1075. At 500 it regressed
    /// from converged to failed.
    ///
    /// UPDATE 2026-09-03: that specific false positive is REMOVED by the
    /// cumulative test described above, and measured to be so - at 588 C the
    /// per-step rule fails at 705 while the cumulative one converges at 1075,
    /// and the whole 301-point sweep is byte-identical with the field armed at
    /// 500 or left at 0. The lesson that motivated the reversal - that a
    /// per-step rule needs the window to clear the worst RECOVERING stagnation
    /// anywhere in the corpus, which nobody can bound in advance - is what the
    /// cumulative test dissolves: it asks whether the window made PROGRESS,
    /// which is scale-free, rather than how long a plateau lasted. The default
    /// nevertheless stays 0, because flipping it changes behaviour on every
    /// project and that is a decision for the project owner, not a measurement.
    /// See plan v5 section 54 for the gates run at the flipped default.
    ///
    /// Currently enabled (500) in: f_/j_CASHNK, where the no-improvement run is
    /// 9781 and the limit hands control to the phase-extinction retry early
    /// (j_CASHNK 10001 -> 563 iterations, same G to every digit); and
    /// f_/j_TestPNTDB, where it is what makes them converge at all - the primary
    /// solve is frozen from iteration 0, and cutting it off reaches a retry that
    /// was previously unreachable.
    long int OptimaStallWindow = 500;

    /// Wall-clock budget in SECONDS for one Optima (AOP/SOP/ROP) solve,
    /// including its retries. 0 (the default) disables it.
    ///
    /// DEFAULT OFF, AND DELIBERATELY SO. Unlike OptimaStallWindow above - which
    /// counts iterations and is therefore bit-reproducible - a wall-clock limit
    /// makes the result depend on the machine, the build and what else is
    /// running. Enabling it globally would make GEMS3K non-deterministic. It is
    /// meant to be set PER PROJECT, in that project's own -ipm.json, for a
    /// system already known to be pathological, as a guard rather than a
    /// tolerance.
    ///
    /// Was set for f_/j_TestSUP98 at ~2x their native solve time, per user
    /// direction 2026-08-27: those two were the only projects in the corpus
    /// where AOP neither converged nor stalled - it progressed slowly and
    /// indefinitely (>30 min against native's ~0.1 s), so the iteration-based
    /// stall detector could not catch them. That setting carried its own
    /// instruction to remove it once the underlying cause was fixed, because it
    /// was a stop-gap rather than a finding.
    ///
    /// REMOVED FROM BOTH FIXTURES 2026-09-02, and the instruction is discharged:
    /// pa_OptimaDimReduce's size gate (section 35) solves both projects on the
    /// compiled default, 2103 it / 80 s and 1695 it / 68 s, so there is no
    /// longer an indefinite run for the budget to guard against. NO PROJECT IN
    /// ANY CORPUS SETS THIS FIELD NOW. It stays as a general per-project guard
    /// for a future pathological system - do NOT enable it globally, for the
    /// determinism reason above.
    ///
    /// Granularity is one Newton iteration - the check runs in the same
    /// per-iteration convergence hook as the stall detector, so a single
    /// iteration longer than the budget cannot be interrupted.
    double OptimaMaxSeconds = 0.0;

    /// Whether the Optima path computes the finite-difference "PartiallyExact"
    /// Hessian columns (Reaktoro's own default strategy). 1 = yes (current
    /// behaviour), 0 = skip them and rely on the analytic ideal-mixing block
    /// plus the regularised exact per-phase block (pa_PhaseHessianFloor).
    ///
    /// EXISTS BECAUSE IT LOOKS REDUNDANT AND EXPENSIVE, but that is not yet
    /// validated enough to change the default. The FD loop runs once per
    /// Optima-reported BASIC variable per iteration - basic variables number
    /// pm.N, the IC count - and each pass is a full
    /// CalculateActivityCoefficients(LINK_UX_MODE) + PrimalChemicalPotentials
    /// over all L species. So f_TestSUP98 pays 82 full activity evaluations over
    /// 923 species EVERY Newton iteration, and the cost grows as N x L, which
    /// matches the measured ~n^2.32 per-iteration scaling.
    ///
    /// Measured 2026-08-27 with it disabled - FEWER iterations and lower cost
    /// everywhere tried, with G bit-identical in every case:
    ///   f_Flowline  93 -> 78 it,   19 -> 12 ms
    ///   f_Solvus   145 -> 67 it,   52 -> 20 ms   (2.6x)
    ///   o_Solvus    96 -> 65 it,   42 -> 23 ms
    ///   mid_1     4184 -> 2620 it, 53 -> 30 s
    ///   f_GEOTHERM 1765 -> 1063 it, 34 -> 4.6 s  (7.4x)
    /// and on the 301-point solvus sweep Tc, worst limb error and out-of-tolerance
    /// count are all IDENTICAL at 2.2x lower cost - including the near-critical
    /// band, which is the case the FD columns were introduced for.
    ///
    /// BUT IT IS NOT REDUNDANT - so it stays ON by default. Measured per-mode
    /// with an adequate budget, disabling it breaks exactly TWO projects:
    ///   f_CASHNK    OK 10003 it -> FAIL 5116
    ///   j_TiQ_PRSV  OK 20 it    -> FAIL 7001
    /// Everything else is unaffected or faster without it, including both
    /// Pitzer systems (j_PitzerTHE 0.96x; f_/j_GEOTHERM 5-7x FASTER) and the
    /// Van Laar Solvus family (0.46-0.69x). So neither "Pitzer aqueous" nor
    /// "has a non-aqueous solution phase" predicts the requirement.
    ///
    /// Treat the requirement as a numerical property of the trajectory, NOT of
    /// the activity model: f_CASHNK needs it and j_CASHNK - same chemistry,
    /// same models, differing only in the thermodynamic-data path - does not.
    /// Do not infer a per-model rule from two cases that disagree across export
    /// formats of one system.
    ///
    /// What is solid is the cost, and that it is wasted on most systems.
    /// pa_PhaseHessianFloor's exact block already covers present end-members of
    /// non-aqueous multicomponent phases, and pure phases have zero true
    /// curvature (their chemical potential is composition-independent), so
    /// their FD columns are noise bought with a full activity evaluation each.
    /// Narrowing this loop to the columns nothing else covers is the follow-up
    /// - see Docs/gems3k-optima-plan-v5.md section 21.
    long int OptimaFDHessian = 1;

    /// Form of the ideal-mixing Hessian block for NON-AQUEOUS multicomponent
    /// (solution) phases in the Optima path. GEMS3K-only, keyword ipm-dat I/O.
    ///   0 (default) - diag(1/X[j]), no off-diagonal. Formally NOT the
    ///                 derivative of the F[] the same code computes.
    ///   1           - the full ideal mole-fraction Jacobian
    ///                 d ln x_j/d n_i = delta_ij/X[j] - 1/Xf, i.e. Leal, Kulik,
    ///                 Smith & Saar (2017) Eq. 80 (Docs/literature/), what
    ///                 Reaktoro's own non-aqueous approxfuncs assembles, and the
    ///                 exact derivative of DC_SYMMETRIC's F = G + ln n_j - ln nSum.
    ///
    /// DEFAULT OFF DESPITE BEING THE CORRECT DERIVATIVE, because it measures
    /// worse on this corpus overall. The extra term is a rank-1
    /// -(1/Xf)*ones*ones^T per block: identically zero along the unmixing
    /// direction, and the whole curvature along the phase-SCALING direction,
    /// where the true block is singular by construction. So it governs only how
    /// freely a phase's TOTAL amount moves.
    ///
    /// Measured 2026-08-27 (AOP):
    ///   f_CASHNK              10003 -> 643 it, 15.6x, same pH/Eh/Vs/Ms
    ///   f_Solvus 400 C          145 -> 101 it ; j_Solvus 108 -> 100
    ///   301-point solvus sweep    0 -> 13 convergence failures,
    ///                             worst limb err 2.7e-3 -> 2.7e-1,
    ///                             iterations 49k -> 187k
    ///   j_CASHNK                563 -> 1600 it ; f_/j_CalcDolo ~+10%
    /// Enable per project for phase-extinction-limited cases (a vestigial phase
    /// decaying slowly toward its floor); leave off otherwise.
    ///
    /// 2026-09-06 - THE COUPLING IS NOW APPLIED ONLY OVER THE END-MEMBERS THAT
    /// ARE PRESENT (ipm_optima.cpp, both objectives; see the long comment at the
    /// main one). The failures above were a hole in pa_PhaseHessianFloor's
    /// regularisation, which requires two present end-members and is therefore
    /// skipped exactly when a phase has been driven to its floor - which this
    /// term is what does. Measured after the fix, same corpus:
    ///   61-point solvus sweep    4 lost temperatures -> 0, worst limb error
    ///                            unchanged at 2.009e-4, Tc unchanged at 655
    ///   f_CASHNK                 1001 it, unchanged (the 10x win is preserved)
    ///   large_f_TestPNTDB        7916 -> 501 it ; mb_10TH_THEREDA 1051 -> 102
    ///   ctest at the flipped default   3 tests failing -> 1
    /// So the field no longer loses an answer anywhere in the corpus (freeze at
    /// the flipped default: 0 status changes and 0 G changes over 258 rows, all
    /// six modes). It is still default-off because the flip is a COST wash
    /// (+2.2 % corpus iterations) and it destabilises one project: t_Solvus640
    /// AOP goes from a 1.3x jitter spread (88-113) to 17.9x (88-1572) - which it
    /// already did before this fix (8.1x, 124-1001), so that is the field's
    /// doing, not the fix's. plan v5 section 85.
    long int OptimaMoleFracHessian = 0;

    /// PSSC-EQUIVALENT PHASE COMPACTION for the Optima path: number of Newton
    /// iterations to spend on a cheap CLASSIFICATION pass before the real
    /// solve. 0 (default) = off, behaviour unchanged.
    ///
    /// *** READ THIS FIRST - 2026-09-07c, plan v5 section 95.4-95.5. ***
    /// EVERY MEASUREMENT BELOW WAS TAKEN AT pa_OptimaDimReduce = -1, i.e. with
    /// the species-level dimension-reduction pre-solve explicitly DISABLED, and
    /// that is not a configuration this solver ships any more.
    ///
    /// pa_OptimaDimReduce is THREE-VALUED and 0 means AUTO, engaging 8 passes at
    /// or above kDimReduceAutoMinDC = 200 species (ipm_optima.cpp). Every project
    /// large enough to have an interesting absent-phase count is therefore already
    /// running it, and it removes the same species this field would pin - earlier,
    /// and from the whole solve rather than from one classification pass. So AT
    /// SHIPPED SETTINGS THIS FIELD IS A STRICT +1-ITERATION NO-OP. Measured on
    /// 07PSIna_G_mid_1, the very project sized below:
    ///
    ///                              probe=0    probe=50   probe=100
    ///   pa_OptimaDimReduce = -1      4184        1243        1846
    ///   pa_OptimaDimReduce = 0/AUTO   269         270         270
    ///
    /// The classification itself is not what fails - the probe pins 117 of 122
    /// phases, the full correct set, and 120 of 141 on T14_ball120. It simply has
    /// nothing left to remove. The same +1 no-op holds on all five T14_ball* rungs
    /// at probe 20/50/100/200, with G identical to nine digits.
    ///
    /// So the 3.4x win below is REAL and REPRODUCIBLE, and it is 4.6x WORSE than
    /// doing nothing at today's defaults. Before reviving this field, check
    /// whether pa_OptimaDimReduce still auto-engages on your project; if it does,
    /// this field cannot help you. gems-benchmark's tools/recheck.py pins the 269
    /// as `veto-compaction-subsumed` so that a change to that gate resurfaces here
    /// rather than silently making this comment true again.
    ///
    /// WHY A CLASSIFICATION PASS AT ALL. Native drops absent phases from the
    /// active set on every call (PhaseSelectionSpeciationCleanup(), pa_PC=2);
    /// the Optima path has always carried every absent phase as a live
    /// box-constrained unknown for the whole solve. Measured (plan-v5 section
    /// 20, and see the precondition above - these figures are at
    /// pa_OptimaDimReduce = -1): the ABSOLUTE COUNT of absent phases tracks this
    /// solver's cost far better than problem size - 07PSIna_G_mid_1 has 117 of
    /// 122 phases absent and runs 4184 iterations / 49.5 s against native's
    /// 124 / 19.6 ms.
    /// Dropping them would cut the unknown count from 265 to ~84, and at the
    /// measured ~n^2.32 per-iteration scaling that alone is ~10x.
    ///
    /// WHY IT CANNOT BE DONE POST-SOLVE, which is the whole reason this is a
    /// separate knob from the phase-selection repair loop that always runs: a
    /// stability index is a function of the DUAL potentials, so nothing can
    /// classify phases before a solve has produced one - and a pass that runs
    /// only after the full solve cannot make that full solve cheaper. A short,
    /// deliberately unconverged probe solve is the cheapest thing that produces
    /// a usable dual.
    ///
    /// SAFETY: a phase pinned here is fixed only in problem.xlower/xupper, not
    /// in pm.DUL/pm.DLL, so it is NOT exempt from the phase-assemblage
    /// stability check - a phase wrongly pinned reports "absent but stable" and
    /// is READMITTED by the phase-selection loop, with its original box
    /// restored and a seed taken from the bulk composition. Misclassification
    /// from a short probe is therefore self-correcting, not silent.
    ///
    /// Suggested starting value if trying this on a many-absent-phase project:
    /// 50-100. Too short and the probe classifies from a still-meaningless
    /// dual (costing readmission loops); too long and the probe costs what it
    /// was meant to save.
    long int OptimaPhaseCompaction = 0;

    /// Iterations of the CHEAP Hessian to attempt in the Optima path before
    /// falling back to the finite-difference PartiallyExact columns.
    /// GEMS3K-only, keyword ipm-dat I/O.  0 (default) = off, i.e.
    /// pa_OptimaFDHessian decides from the very first iteration, which is the
    /// behaviour shipped before this field existed.
    ///
    /// N > 0 makes ONE attempt of up to N iterations with the FD loop
    /// suppressed.  If it converges, its state is adopted and the primary
    /// solve re-runs from it STILL CHEAP, confirming an already-converged
    /// point in a handful of iterations.  If it does not, that trajectory is
    /// DISCARDED and the primary solve restarts from the original seed with
    /// the FD columns on.  Every retry after the primary solve uses FD
    /// regardless - reaching a retry is itself an escalation trigger.
    ///
    /// This is Reaktoro's own published arrangement: Leal, Kulik, Smith & Saar
    /// (2017) state that the exact-Newton method "is only used if the
    /// quasi-Newton method ... fails or takes more than a certain number of
    /// iterations to converge".  We otherwise run the exact columns always.
    ///
    /// TWO SHAPES THAT LOOK RIGHT AND ARE MEASURABLY WRONG - do not "simplify"
    /// back to either (plan v5 section 30.3 has the numbers):
    ///   * Switching FD on IN PLACE at iteration N, carrying straight on,
    ///     fails the 301-point solvus sweep at T = 560 C - a point that needs
    ///     556 FD iterations from the start and hits the 7001 cap under the
    ///     cheap Hessian - at EVERY delay tried (300/1000/1500/3000).  The
    ///     early cheap iterates land somewhere the exact columns cannot
    ///     recover from, so the trajectory must be discarded, not corrected.
    ///   * Re-running the primary solve with FD ON after a SUCCESSFUL cheap
    ///     attempt, to make every reported result FD-validated, returns
    ///     BAD_GEM_AOP on 142 of 301 solvus points with a KKT residual of ~0.5
    ///     at H+, at phase amounts identical to the passing run's to six
    ///     figures.  The FD columns evaluated AT a converged point are
    ///     near-singular.
    ///
    /// The trigger MUST be iteration count, never convergence: the cheap
    /// Hessian alone CONVERGES on the Solvus family, so a success-gated
    /// escalation would never fire.  (It converges to the RIGHT answer there
    /// because pa_PhaseHessianFloor's exact per-phase block is a separate
    /// mechanism that always runs - that block, not this loop, is what
    /// supplies the unmixing curvature.)
    ///
    /// Why it pays: with pa_OptimaFDHessian = 0 the corpus is 1.5-7.4x faster
    /// and bit-identical in G, and exactly two projects break (f_CASHNK,
    /// j_TiQ_PRSV - see pa_OptimaFDHessian).  A delay large enough to cover
    /// the slowest cheap-Hessian convergence (f_/j_GEOTHERM need ~1100) keeps
    /// that speed for everything that converges cheaply while still reaching
    /// the exact columns for the two that need them.  Measured at N = 1500:
    /// 11x total wall time on the 17-project subset and 2.4x on the 301-point
    /// solvus sweep, every G identical and the limb accuracy unchanged.
    ///
    /// MEASURED PER-PROJECT, NOT A DEFAULT (full 25-project sweep, 2026-08-28):
    /// enable on f_/j_GEOTHERM (8.5-10.5x), j_10TH_G_seawater (FAIL -> OK, 10000
    /// it -> 229, at native's G bit-for-bit), the Solvus family and j_PitzerTHE
    /// (2-4x).  Do NOT enable on f_/j_CASHNK, j_TiQ_PRSV or f_/j_TestPNTDB.
    /// f_CASHNK regresses 16x (643 -> 11504 it) and corpus total wall time is
    /// unchanged, so a blanket default buys nothing on average.
    ///
    /// That regression is NOT just the wasted attempt, and it is the thing to fix
    /// before reconsidering a default: discarding the attempt resets Optima::State
    /// but NOT the MULTI scratch state this objective mutates (pm.X, pm.F,
    /// pm.Gamma, ...), so the restarted solve is not on the trajectory the
    /// undelayed run takes - on f_CASHNK it misses the stall-detector hand-off at
    /// ~642 iterations and runs to its 10000-iteration maxiters cap instead.
    ///
    /// KNOWN INTERACTION, check it before raising this much further: the stall
    /// detector (pa_OptimaStallWindow; default 0 = off, see that field) can fire during the cheap
    /// attempt, which is treated as a non-convergence and so restarts with FD
    /// - the safe direction, but it makes a very large N wasteful.

    /// RE-GATED 2026-09-06 to COLD STARTS ONLY (pm.pNP == 0), for the same reason
    /// OptimaReducedPreSolve() is: a warm call already carries the consistent
    /// (primal, dual) pair this probe exists to manufacture, so there it is pure
    /// waste. Measured on a corpus freeze with the delay on before gating it - the
    /// warm modes were taxed systematically, about twenty SOP/SHP rows going from 1
    /// iteration to 2, for SHP +42 % and SOP +2 % overall, while AOP (the cold mode
    /// it is for) was -5 % and HOP, whose Optima leg is warm, was +7 %.
    ///
    /// THE DEFAULT FLIP WAS RE-TESTED 2026-09-06 UNDER THE NEW DECISION RULE AND
    /// REJECTED AGAIN - but for a sharper reason than the cost ratchets it trips.
    /// On the 121-step Cu-Pourbaix titration AOP loses 8 converged steps where native
    /// keeps all 121: that is LOSING ANSWERS on a sweep, which is the workflow AOP is
    /// kept for, not merely being slower on some projects. The corpus freeze is
    /// otherwise a wash - ZERO status changes in any of the six modes, AOP -5 % - with
    /// the wins concentrated on the Solvus family and f_GEOTHERM (1.6-2.3x) and one
    /// large loss, 10TH_G_00001 AOP 107 -> 1107. Set it per project on
    /// f_/j_GEOTHERM, j_10TH_G_seawater, the Solvus family and j_PitzerTHE; do NOT
    /// set it on f_/j_TestPNTDB, 10TH_G_00001 or Cu-Pourbaix_G_pHtitr.
    long int OptimaFDHessianDelay = 0;

    /// Lower bound on species amounts in the Optima path ("dcFloor"), in the
    /// same internally rescaled units as pm.X.  GEMS3K-only, keyword ipm-dat
    /// I/O.  0 (default) = derive it from pa_DHB, exactly as before this field
    /// existed; > 0 = use this value instead.
    ///
    /// WHY IT EXISTS - pa_DHB is overloaded, and that was measured rather than
    /// argued.  It is simultaneously (a) native's RELATIVE mass-balance
    /// tolerance in MassBalanceRefinement(), (b) this solver's box lower bound,
    /// and (c) via dcFloor*10, the phase presence threshold.  So lowering it to
    /// get a smaller species floor necessarily tightens native's convergence
    /// test as well: on 07PSIna_G_simple_0_0_0_150, pa_DHB = 1e-20 and below
    /// makes NATIVE fail outright (inside MBR, it=30/FIA=30/IPM=0) while the
    /// Optima path still returns a state.  There was no way to test a low floor
    /// without breaking the reference.
    ///
    /// Leal, Kulik, Smith & Saar (2017) publish tau = the species floor at
    /// <= 1e-25, twelve orders below our 1e-13..1e-15, and call 1e-8 known-bad
    /// "for trace elements or redox reactions with very low amounts of species"
    /// - which is exactly the O2(aq)/H2(aq) pair this branch has chased since
    /// the first Tier A session.  This field is what makes that testable.
    ///
    /// Measured so far (07PSIna_G_simple_0_0_0_150, whose Eh is wrong by 0.52 V
    /// purely because its only redox species sit at the floor): lowering the
    /// floor moves AOP's Eh 0.4968 -> 0.0395 and then stops, while native's own
    /// Eh moves the other way (-0.0261 -> -0.2104) when pa_DHB is lowered.  The
    /// two do not meet, which is consistent with that system's Eh being
    /// genuinely underdetermined rather than merely floor-limited.  So this is
    /// an instrument for the question, not a known fix.
    double OptimaDcFloor = 0.;

    /// Per-IC-class mass-balance convergence rule in MassBalanceRefinement().
    /// DEFAULT 0 = OFF, and off means the existing code path is taken unchanged.
    ///
    /// Kulik et al. 2013 (GEM IPM3), Appendix 2.2 prescribes a per-IC-CLASS test:
    /// a RELATIVE threshold for minor/trace independent components and an ABSOLUTE
    /// threshold for major ones ("for instance, H and O in aquatic systems").
    /// The implementation has never done that - with pa_DT != 0 it applies BOTH
    /// tests to EVERY IC (`||`), which is strictly stricter than either and is
    /// never a per-class selection, under a comment stating the paper's intent.
    ///
    /// When > 0, this value is the TRACE/MAJOR ratio threshold: IC i is treated as
    /// trace when B[i] < MbClassRule * max_k B[k], major otherwise. A ratio rather
    /// than an absolute amount, because pm.B[] is internally rescaled to pa_DG total
    /// moles, so an absolute classification would not be scale-invariant.
    /// Trace ICs then get the relative test, major ICs the absolute test.
    ///
    /// The major absolute cutoff is 10^-|DT| when |DT| >= 2 (pa_DT is the field that
    /// already encodes it, and naming it explicitly is the recommended usage);
    /// otherwise it falls back to DHBM*1e5, so the switch is usable standalone.
    /// DHBM-scaled cutoffs are an existing idiom here - cf. min(DHBM*1e10, 1e-2)
    /// in CheckMassBalanceResiduals().
    ///
    /// MEASURED 2026-08-29 over all 50 loadable benchmark projects with
    /// gems-benchmark's tools/mb_class_probe: 7 projects would flip FAIL -> PASS
    /// (07PSIna_G_ironsi, 10TH, Al-species, CSHSnplus, FeNaCl_FyGt_{Precip,
    /// Precip_HighpH, TransitionZone}) and 0 would flip PASS -> FAIL. Those 7 are
    /// seven of the eight projects whose native SIA cannot re-solve its own
    /// converged answer. The mechanism is NOT trace ICs held to a relative bar -
    /// it is MAJOR ICs (H, O, K) held to a relative bar where the paper says
    /// absolute: a major IC at ~111 mol against DHB = 1e-13 demands ~1e-11 mol
    /// absolute accuracy, near the arithmetic floor.
    ///
    /// DEFAULT OFF because this RELAXES a native convergence test on a path every
    /// GEM_run() takes. Turn it on per project to test; do not make it a default
    /// without its own full-suite gate.
    double MbClassRule = 0.;

    /// Trend-based vanishing-phase detection in the Optima path's phase-extinction
    /// retry. DEFAULT 0 = OFF; off takes the pre-existing twin-only path unchanged.
    ///
    /// Leal 2014 §2.3.3 defines an unstable phase with TWO clauses: below a threshold
    /// **and decreasing since the last iterate**. The extinction tier here uses an
    /// exactness argument instead (two phases with identical stoichiometry and G0 are
    /// interchangeable) precisely because both magnitude-only detectors fail on
    /// f_CASHNK - the redundant phase is at ~1e-5 and still falling when the solver
    /// gives up, and its logSI is ~0 because it shares a tangent plane with its twin.
    /// A trend test is the third option neither considers, and unlike the twin proof
    /// it does not require a twin to exist at all.
    ///
    /// When > 0, a phase is deactivated if its total fell monotonically for at least
    /// `MbTrendPhaseDecay` consecutive objective evaluations AND is now below 1e-2 of
    /// its own peak. Tracking lives in the objective callback, the only place that
    /// sees every iterate.
    ///
    /// MEASURED 2026-08-29: on j_CASHNK the trend criterion ALONE reproduces the twin
    /// criterion's answer exactly (G = -6.028481770e+03), i.e. it SUBSUMES it on the
    /// only case we have. Zero false positives on Solvus (both), Flowline, CASH+ and
    /// TiQ_PRSV; CTest 7/7 green with it enabled. Safety is partly structural - the
    /// tier only runs after a solve has already failed, so a converging project never
    /// reaches it.
    ///
    /// DEFAULT OFF for two honest reasons. (1) It reintroduces TUNED CONSTANTS (the
    /// consecutive-decrease count, the 1e-2 peak ratio) where the twin criterion was
    /// exact and project-derived - and this codebase already carries more hand-tuned
    /// thresholds than is healthy. (2) Its general case - a decaying phase with no
    /// twin - has NO test system yet; that is T3 in
    /// Docs/literature/PROPOSED-TEST-SYSTEMS.md. Until T3 exists its value cannot be
    /// properly assessed, only its harmlessness.
    long int MbTrendPhaseDecay = 0;

    /// Iteration budget for the FIRST Optima attempt only, so the existing
    /// phase-selection repair loop gets its turn EARLY instead of after the
    /// primary solve has run to convergence. 0 = off (unchanged behaviour).
    ///
    /// Motivation, measured on Resources/gems3k/f_Solvus_G_Test1_0_0_1000_400_0
    /// (2026-09-02): warm-restarting with Al,Na,K scaled by 0.998 makes
    /// Plagioclase dissolve slowly - ~6958 consecutive monotone decreases - and
    /// costs 7007 Optima iterations, against 71 at a 0.995 scale.  The repair
    /// loop DOES diagnose it correctly ("present but unstable, logSI gap 0.63")
    /// and removing the phase up front costs only 58 iterations - a 121x
    /// oracle win - but the check runs only after the primary solve has already
    /// spent ~6180 of them.  So the criterion is right and only its TIMING is
    /// wrong; this caps the first attempt so it fires sooner.
    /// NOT the same as pa_OptimaPhaseCompaction, which pins phases below a
    /// PRESENCE threshold - measured harmful here (7059/14168/21377 iterations
    /// at probe 50/100/200), because this phase is present-but-dissolving.
    ///
    /// NOW USABLE - the probe-then-full-budget safety net this comment used to
    /// prescribe was implemented 2026-09-03 (plan v5 section 58) and measured.
    /// Before it this field was unusable, and for the reason recorded here: the cap
    /// applies to the first attempt of EVERY call, including calls whose assemblage
    /// is already correct and which simply need their budget, so the same project's
    /// own cold solve (823 iterations) FAILED outright at a cap of 500 - status=13,
    /// and a G wrong in the 7th digit - because the repair loop had no stability
    /// violation to repair, only a budget shortfall. Re-confirmed on the current
    /// tree, by disabling the net, before it was adopted.
    ///
    /// The net treats the capped attempt as a PROBE: if the cap is what stopped it
    /// and nothing downstream recovered, re-solve ONCE from the original state at
    /// the full budget. So the only cost of a probe that finds nothing is the probe
    /// itself, and the worst case of setting this field is one wasted solve on a run
    /// that was otherwise fine - never a lost answer. Same one-way-bet argument, and
    /// the same shape, as pa_OptimaStallWindow's own net.
    ///
    /// MEASURED on Resources/gems3k/f_Solvus_G_Test1_0_0_1000_400_0 (2026-09-03):
    ///
    ///                       cap=0      cap=200      cap=500
    ///   decay case (warm)   7005 it    209 it       509 it     <- 33.5x at 200
    ///   plain cold solve     823 it   1024 it      1324 it     <- = cap + 824
    ///
    /// G, pH and Vs are identical to every digit in all six runs, and the decay
    /// case's Plagioclase lands on 5.355388e-16 in all three. The win is the repair
    /// loop firing at iteration 200 ("deactivating phase 1 - present but unstable,
    /// logSI gap 0.833") instead of after the primary solve has already spent ~6180
    /// iterations converging on an assemblage it then has to correct.
    ///
    /// NEGATIVE VALUES select a TREND trigger instead of a fixed cap: -N ends the
    /// first attempt as soon as some NON-SOLVENT multicomponent phase has fallen
    /// monotonically for N consecutive objective evaluations. Built 2026-09-03
    /// (plan v5 section 59) to make the probe FREE - it is spent only on a run that
    /// shows the signature, where an ordinary cold solve pays the positive form's
    /// cap on every call. It reuses pa_MbTrendPhaseDecay's own counters, and it is
    /// that criterion's first consumer a CONVERGING run can reach (its other one
    /// sits inside `if( !result.succeeded )`).
    ///
    /// It IS free, and where it acts it beats the cap:
    ///   f_Solvus_G_Test1 decay case  7005 -> 95 it  (73.7x; the cap gives 209)
    ///   f_Solvus_G_Test1 cold solve   823 -> 823 it (UNCHANGED; the cap gives 1024)
    ///   j_CASHNK                     1001 -> 56 it  (17.9x, G/pH/Vs identical)
    /// and 40 of the 43 projects in Resources/gems3k + gems3k-fail are untouched.
    ///
    /// *** BUT IT CAN RETURN A WRONG ANSWER, AND THE SAFETY NET DOES NOT CATCH IT.
    /// Measured on gems3k-fail/CSHSnplus_G_CSH1_5_bufs: OK 1298 it -> BAD 197 it
    /// with G worse by 0.8 RT and pH off by 0.15. The net only fires when the early
    /// look found NOTHING to repair; here it found something and acted on it, using
    /// a dual that was not yet accurate enough for the stability index to be right.
    /// So the "worst case is one wasted probe" argument covers the CAP form and NOT
    /// this one - looking early means acting on an inaccurate dual, and the repair
    /// loop's criterion is exact only to the extent that dual is.
    /// (Also measured: CASH+CsSr 106 -> 168 it at an identical answer - a true decay
    /// signature on a phase the repair loop then correctly declined to act on, so a
    /// magnitude clause would not have prevented it either.)
    ///
    /// A DEFAULT OF -50 WAS TRIED AND REJECTED on exactly that: `ctest` goes 8/8 ->
    /// 7/8, `fail_CSHSnplus` reporting "absent but stable - should be present".
    /// Note a 43-project corpus sweep did NOT catch it and the CTest suite did -
    /// the sweep silently skipped the two projects whose -dat.lst basename differs
    /// from their directory name.
    ///
    /// STILL DEFAULT 0, therefore, for BOTH forms. The cap form is safe but costs
    /// ~1.24x on an ordinary cold solve; the trend form is free but can be wrong.
    /// Against the refined decision rule (add it if it converges AND takes fewer
    /// iterations AND is faster, for all cases or with a smart fallback) neither
    /// form passes: one fails the strict-improvement clause, the other fails the
    /// fallback clause.
    ///
    /// What a defaultable version needs is a way to tell a trustworthy stability
    /// verdict from a premature one - e.g. require the dual to have settled before
    /// the early look is allowed to ACT, rather than gating only when it is taken.
    /// 2026-09-06 - THE TREND FORM IS NOW GATED ON THE DUAL HAVING SETTLED
    /// (ipm_optima.cpp, kEarlyTrendDualSettled; plan v5 section 87). The repair
    /// loop's verdict is a stability index computed from the dual, so acting on a
    /// dual that is still moving is how the trend form went wrong. Measured at the
    /// moment the trigger fires: max|dw|/max|w| = 2.97e-15 where acting is RIGHT
    /// (j_CASHNK) against 1.33e-04 where it is WRONG (CSHSnplus) - eleven orders
    /// apart, so the threshold (pa_OptimaTol) is not a tuned knob.
    /// With it, CSHSnplus goes BAD 197 -> OK 2119 at the un-armed answer to every
    /// digit, j_CASHNK keeps 1001 -> 56, and ctest passes 11/11 at the CHANGED
    /// default - which section 59.6's flip did not.
    ///
    /// The flip is STILL rejected, on a regression only the benchmark freeze sees:
    /// j_CASHNK's HOP row goes OK 4587 -> FAIL 5227. That is this field's own,
    /// pre-existing (the gate can only delay a trigger, never create one) and was
    /// simply never measured, because section 59.6 stopped at ctest. Cause: on the
    /// WARM leg the trigger fires on a phase that fell 0.09% - 2.212 to 2.210 over
    /// 342 evaluations - i.e. SETTLING, not decaying. The dual is accurate from
    /// iteration 1 on a warm start, so this gate is a no-op exactly there.
    /// The obvious fix - a magnitude clause - is blocked by the measurement
    /// recorded at kEarlyTrendDropRatio: on the decay case the dissolving phase is
    /// still at 14% of its peak when the look should happen, so any 10x clause
    /// makes the "early" look not early. The separation is 14% / 91% / 99.91%,
    /// which is a tuned constant with ~2x margin, not the 11 orders above.
    ///
    /// 2026-09-09 - BOTH FORMS FLIPPED CORPUS-WIDE AND SCORED (plan v5 section
    /// 101), 412 rows x 3 corpora x 6 modes each against the frozen standard.
    /// Both exit 1, AND THEY FAIL ON OPPOSITE ROWS - each form fixes exactly the
    /// row the other breaks:
    ///
    ///   row                baseline      cap = 200      trend = -50
    ///   CSHSnplus  AOP     OK 1298       FAIL           OK 2119
    ///   j_CASHNK   HOP     OK 4586       OK 293         FAIL 5810
    ///
    /// Cost: cap 0.884x overall (SHP 0.610x), trend 1.061x (SHP +29.6 %). Native
    /// and SIA are untouched by either, so this field cannot affect the
    /// native-only beta. Note the cap form REMOVES the j_CASHNK HOP regression
    /// that had kept the trend form unflipped for a month - so that blocker was
    /// the TRIGGER's, not the field's.
    ///
    /// 2026-09-09 - THE CAP FORM IS NOW GATED ON THE DUAL TOO, and is no longer a
    /// maxiters budget at all (plan v5 section 102; ipm_optima.cpp, alongside the
    /// trend form's kEarlyTrendDualSettled). The two forms had asymmetric gates
    /// for a structural reason and not a deliberate one: a budget is spent inside
    /// Optima's own stepping loop, so the capped attempt never entered the
    /// convergence hook where both iterates of u = (x, p, w) - and hence the
    /// dual's movement - are available. A budget cannot be conditional. The cap
    /// is therefore evaluated in that hook now, ONCE, at iteration N: "if the dual
    /// has settled (max|dw|/max|w| <= pa_OptimaTol) end the first attempt here;
    /// otherwise leave this call alone". N stops being a budget and becomes the
    /// single point at which the early look is offered.
    ///
    /// ASKED ONCE, NOT WAITED FOR - and that is a measured decision, not a
    /// simplification. The first version deferred instead ("stop at the first
    /// iteration >= N with the dual settled") and MEASURED WORSE THAN DOING
    /// NOTHING: f_Solvus_G_Test1 AOP 823 -> 1530, CSHSnplus AOP 1298 -> 2119.
    /// The reason is structural. On a converging run the dual settles when the run
    /// is nearly over - f_Solvus's settles at iteration 707 of 823 - so a deferred
    /// look is not early, and the probe it spends is very nearly a whole solve.
    /// Deferring turns the cap's honest bounded cost ("waste N iterations") into
    /// an unbounded one, on exactly the ordinary runs that have nothing to repair.
    /// Asking once at N keeps the bound: the look happens where the dual is
    /// ALREADY settled at N - a warm leg, or a run whose dual is determined early,
    /// which is where every one of the cap's wins lives - and both rows above
    /// return to their baselines exactly.
    ///
    /// Measured with the cap armed, so every sample is from an evaluation <= 200:
    /// the smallest dual movement anywhere inside CSHSnplus's capped window is
    /// 1.18e-03 against the 1e-8 threshold - five orders above it - and at the
    /// earliest sample the dual is moving by 100 % of its own magnitude. The gate
    /// blocks there by a wide margin, not a close one. On j_CASHNK it is 2.97e-15,
    /// twelve orders below, so the gate is a no-op and the HOP win is untouched.
    ///
    /// Also corrected the same day, against what section 101.4 and the -08d
    /// handoff both asserted: the cap's loss on CSHSnplus is NOT a budget
    /// shortfall the safety net should have absorbed. The net fires, the
    /// full-budget re-solve converges in 727 iterations, and the row still fails
    /// on the ASSEMBLAGE - the early look deactivated CASH+Sn on a logSI gap of
    /// 0.044 that the converged state contradicts by two orders (9.01). It is a
    /// wrong decision taken on noise, which is exactly what a dual gate is for and
    /// exactly what a bigger budget is not.
    long int OptimaEarlyStabilityAt = 0;

    /// Species-level dimension reduction for the Optima path (AOP/SOP only).
    /// 0 = OFF (unchanged behaviour); > 0 = ON, and the value is the maximum
    /// number of readmission ("pricing") passes the reduced pre-solve may make.
    ///
    /// WHY. Native's linear system is N x N over independent components; the
    /// Optima path hands Optima dims.x = L + R - species PLUS components. On
    /// f_/j_TestSUP98 that is ~1005 x 1005 against native's 82 x 82, and the
    /// factorisation cost measured on this corpus scales as ~n^2.32, so the
    /// size wall is a DIMENSION problem: pinning an absent species with a
    /// degenerate box (pa_OptimaPhaseCompaction) leaves the system exactly as
    /// large and cannot help, whereas OMITTING it can. Measured with
    /// tools/dimension_ceiling: a median 59% of species are absent at native's
    /// own converged answer, 68% on the TestPNTDB/TestSUP98 pair (13-14x per
    /// iteration) and 84-92% on the Resources/gems3k-psina giants (71-322x).
    /// The absent species are mostly AQUEOUS (359 of f_TestPNTDB's 471), which
    /// is why this is species-level and not phase-level.
    ///
    /// HOW. A pre-solve, before the ordinary full-dimension solve:
    ///   - start from the LP-feasibility seed's own support WIDENED by pricing
    ///     every species against the linearised-Gibbs LP's dual, keeping those
    ///     below OptimaDimReduceTol (see that field - the widening is what makes
    ///     this work at all on a large system). Both halves are LPs over
    ///     pm.A/pm.B alone and so are INDEPENDENT of any solver state, which is
    ///     required: deriving the set from a failing solver's own iterate is
    ///     measured-unsafe, since on j_10TH_G_seawater 4 of the 60 floor-pinned
    ///     species are genuinely present, dolomite among them at 1.0e-2 mol, so
    ///     that set would omit dolomite and turn an honest failure into a silent
    ///     wrong answer;
    ///   - solve over the active species only, with every omitted species held
    ///     at its own lower bound and its contribution moved into the RHS
    ///     (be[i] -= sum over omitted j of A[i,j]*xlower[j]);
    ///   - PRICE the omitted columns on the resulting dual - readmit every one
    ///     whose reduced gradient s_j = F[j] - sum_i U[i]*A[i,j] is negative,
    ///     i.e. exactly Optima's own is_lower_unstable test - and repeat.
    /// Readmission is MONOTONE (a species is never dropped again), which is
    /// what bounds the loop and keeps it free of the active-set cycling this
    /// solver has repeatedly shown when a set is allowed to shrink again.
    /// At a fixed point every omitted species is at its lower bound with a
    /// non-negative reduced gradient, which is the full problem's own KKT
    /// condition - so the reduced answer is a solution of the FULL problem,
    /// not an approximation to it.
    ///
    /// SAFETY. The pre-solve never returns an answer on its own: it writes its
    /// primal into pm.Y[] and its dual into pm.U[], and the ordinary
    /// full-dimension solve then runs warm-started from both. Warm-starting a
    /// correct (x,y) pair costs O(1) Optima iterations (measured: 3 on
    /// f_TestPNTDB, section 27), so the full solve becomes a cheap verification
    /// that also keeps every existing retry, KKT, mass-balance and
    /// phase-stability check running at full dimension and full index. If the
    /// pre-solve fails or gains nothing, it is discarded and the full solve
    /// runs exactly as it would have.
    ///
    /// NOT applied when control conditions are active (R > 0) or in ROP
    /// (reaktoroMode), which is a faithful port of Reaktoro's own mechanism.
    ///
    /// THREE-VALUED, and 0 is AUTO rather than off (changed 2026-09-02, see
    /// section 35 of Docs/gems3k-optima-plan-v5.md):
    ///   > 0  ON, and the value is the readmission-pass limit;
    ///   = 0  AUTO - on with kDimReduceAutoPasses passes when the system has at
    ///        least kDimReduceAutoMinDC species, off below that
    ///        (both in ipm_optima.cpp, at the call site);
    ///   < 0  explicitly OFF whatever the size.
    /// The size gate is there because the corpus sweep in section 33.1 measured
    /// the payoff as tracking the fraction of species the reduction can OMIT,
    /// which is strongly size-dependent: 85-96 % of species stay active below
    /// ~35 species (nothing to omit, so the pre-solve is pure overhead) against
    /// 49-56 % from 122 species up. Every measured loss is below the gate and
    /// every large win above it - see the plan section for the table, and for
    /// the honest caveat that nothing in the corpus sits between 154 and 265
    /// species, so the gate's exact value is interpolated rather than measured.
    long int OptimaDimReduce = 0;

    /// Threshold (in RT units) for OptimaDimReduce's INITIAL active set, and
    /// the whole reason that feature works on a large system at all. Read
    /// whenever the reduction actually runs, which since the size gate means
    /// OptimaDimReduce > 0 OR OptimaDimReduce == 0 (AUTO) on a system at/above
    /// kDimReduceAutoMinDC species - NOT only when the field is positive.
    ///
    /// The initial set was originally the LP-FEASIBILITY seed's own support -
    /// safe (independent of any solve, feasible by construction) but a VERTEX,
    /// so only N species. Measured 2026-09-02, that sparsity is what broke the
    /// feature on the two largest projects: 07PSIna_G_complex_1's 69-species
    /// pass 0 was itself unsolvable, and f_TestPNTDB's 49-species pass 0 priced
    /// 427 of its 641 omitted columns back in one step.
    ///
    /// Instead, price EVERY species against the dual of the LINEARISED-GIBBS LP
    /// (LPGibbsDual(), min sum G0[j]*n_j s.t. A n = b, n >= 0 - the same simplex,
    /// so still wholly independent of any solve) and activate every species with
    ///     s_j = G0[j] - sum_i y_i A[i,j]  <  OptimaDimReduceTol.
    /// s_j is how far species j is from being stable under a first-order model,
    /// in RT units, so the threshold is a physical quantity rather than a tuning
    /// index: a solution species can lower its own chemical potential by at most
    /// the mixing-entropy bonus -ln(x_j), i.e. ~35 RT at x_j ~ 1e-15, and native's
    /// own converged dual puts the largest s among genuinely PRESENT species at
    /// 32.5-34.1 on all three projects measured. A tolerance well below that is
    /// deliberately NOT trying to capture every present species - anything missed
    /// simply prices back in on the next pass, which is why the initial set can
    /// affect cost without ever affecting correctness.
    ///
    /// Measured (total pre-solve + verification iterations / wall):
    ///   07PSIna_G_mid_1   (265 sp): 2455 it / 3.85 s at the old support rule;
    ///                               261 it / 0.36 s at 5; 269 it / 0.45 s at 10;
    ///                               978 it / 2.05 s at 20; 5396 it / 27.2 s at 50
    ///   f_TestPNTDB       (690 sp): pre-solve discarded at the old rule (~500 s);
    ///                               discarded at 5; 501 it / 2.2 s at 10;
    ///                               747 it / 5.2 s at 20
    ///   07PSIna_G_complex_1 (1392): pass 0 unsolvable at the old rule, and the
    ///                               project has never converged at all;
    ///                               1602 it / 206 s at 10, native's G to 9 s.f.
    /// Hence the default 10. Too small and pass 0's dual is poor enough that the
    /// readmitted set is unsolvable (f_TestPNTDB at 5); too large and the single
    /// big pass costs more than the reduction saves (mid_1 at 50). The resulting
    /// initial sets sat at 3.4-4.6 x N on all three, which suggested a
    /// dimensionless rank rule as an alternative - now implemented and measured,
    /// see the sign overload below.
    ///
    /// SIGN OVERLOAD (section 33.3): a POSITIVE value is the RT threshold
    /// described above; a NEGATIVE value -m is the dimensionless RANK rule
    /// "admit the cheapest-priced species until the active set reaches m x N".
    /// Same LP, same prices, only the selection differs, so it inherits the
    /// independence argument unchanged.
    ///
    /// MEASURED ON ALL ELEVEN 200+ SPECIES PROJECTS, 2026-09-02 (section 39) -
    /// total iterations, identical G/pH/Vs in every converging arm:
    ///
    ///   project              sp     N   tol 10   -2      -2.5   -3     -4
    ///   07PSIna_G_mid_1      265   19     269    119     184    168     -
    ///   f_TestPNTDB          690   42     501    396     433    491     -
    ///   j_TestPNTDB          690   42     495   FAIL     425    448     -
    ///   f_TestSUP98          923   82    2103   1494    FAIL   1690     -
    ///   j_TestSUP98          923   82    1695   1476      -      -      -
    ///   complex_1 (25 C)    1392   59    1601   1216      -    2031     -
    ///   complex_1 (80 C)    1392   59    FAIL   FAIL      -    FAIL     -
    ///   edt_2               1566   61    FAIL    936    1224    788    982
    ///   vcomplex  (25 C)    1566   61    3492    667      -      -      -
    ///
    /// TWO CONCLUSIONS, and the second is why the default did not move:
    ///
    /// 1. The rank rule is worth setting PER PROJECT. It is the only thing that
    ///    solves 07PSIna_G_edt_2 at all (the default's pass 0 settles on 218
    ///    species, readmits 249, and the resulting 467-species pass never
    ///    finishes), and it cuts 07PSIna_G_vcomplex 5.2x. Both reach native's G
    ///    to 8-9 significant figures.
    ///
    /// 2. Its response is NOT monotone, so no value is adoptable as a default:
    ///    j_TestPNTDB FAILS at -2 and works at -2.5/-3, while f_TestSUP98 works
    ///    at -2, FAILS at -2.5, and works at -3. f_TestSUP98's failure is
    ///    LOGGED ("pass 1 did not converge on 529 of 923 species"), reproduces
    ///    uncontended with a bit-identical pass 0, and shows that SET SIZE DOES
    ///    NOT PREDICT SOLVABILITY on that project: 164 works, 205 fails, 246
    ///    works. It is which species, not how many. This corrects section
    ///    33.3's "single sharp optimum at 2 N ... monotone-with-a-cliff, almost
    ///    unique among this branch's knobs" - true on the two projects it was
    ///    measured on, false on eight.
    ///
    /// Note j_TestPNTDB's -2 failure is a COST regression only (9798 iterations
    /// / ~400 s against 495 / 3.3 s) - the pre-solve is discarded and the full
    /// solve reaches the same G - and that its f_ twin is FASTEST at the same
    /// setting. The two exports differ only in thermodynamic data, and here
    /// that difference decides whether an 84-species reduced problem is
    /// solvable at all; never pool an f_/j_ pair.
    ///
    /// The way to make the rank rule safe is NOT a better m: it is to retry the
    /// pre-solve once with a different initial set when it is discarded (the
    /// two rules fail on DISJOINT projects, so "rank first, shipped threshold
    /// as fallback" would have solved all eleven). Not built - it is a change
    /// to the retry chain, which this branch has repeatedly found reshuffles
    /// outcomes chaotically, and it needs the full corpus gate.
    double OptimaDimReduceTol = 10.;

    /// Appendix A of Leal et al. (2017), Eqs. 130-136 - a pivot/non-pivot
    /// SPLIT of native MBR's Schur-complement reduction. 0 = off (the naive
    /// reduction, unchanged), non-zero = on. NATIVE path only; nothing in the
    /// Optima path reads it.
    ///
    /// MakeAndSolveSystemOfLinearEquations()'s initAppr branch assembles
    ///     A[i,k] = sum_j a(j,i) a(j,k) W[j]
    /// which is the paper's own Eq. 82/83 reduction A D^-1 A^T with
    /// D_jj = 1/W[j] - the same equation, not an analogy. The paper warns that
    /// this reduction "often fails because of round-off errors ... since no
    /// pivoting was performed to avoid division by small numbers in the
    /// evaluation of D^-1".
    ///
    /// Eq. 133 splits the species by comparing each |D_jj| against the
    /// infinity-norm of the matching column of C (here C_(i,j) = a(j,i)):
    ///     non-pivot  <=>  1/W[j] < max_i |a(j,i)|  <=>  W[j]*max_i|a(j,i)| > 1
    /// Under MBR's QUADRATIC weight W[j] = (Y[j]-DLL[j])^2 the non-pivot set is
    /// therefore the ABUNDANT species - water and the major solutes - so this
    /// keeps water in the joint solve rather than eliminating it through a
    /// division, which is at least adjacent to the water H:O = 2:1 near-
    /// dependence traced on 2026-08-21. Eqs. 135/136 eliminate only the pivot
    /// block and solve the non-pivot unknowns jointly with the duals, giving a
    /// system of dimension N + |I_n|; when I_n is empty it degenerates exactly
    /// to today's assembly, and the implementation falls through to it.
    ///
    /// DISTINCT from the Jacobi preconditioner in the same function (ported
    /// 2026-09-02): Jacobi RESCALES A D^-1 A^T after assembly, Appendix A
    /// refuses to FORM part of it. Both are applied when this is on - comparing
    /// an unpreconditioned Appendix A against a preconditioned baseline would
    /// measure the loss of the preconditioner, not the gain of the split.
    ///
    /// HONEST BOUND: it stops amplification THROUGH D^-1; it cannot repair the
    /// genuine near-singularity of A D^-1 A^T itself, which is what water's
    /// fixed H:O = 2:1 stoichiometry produces. Plausible partial help, not a
    /// claimed cure - and the gate for it is "no regression", since native
    /// already solves essentially the whole corpus.
    ///
    /// The augmented system is symmetric but INDEFINITE (a saddle point: the
    /// non-pivot diagonal block is -1/W[j]), so it is solved by LU directly
    /// rather than attempting Cholesky first as the naive path does.

    long int MbPivotSplit = 0;

    /// Report a species that is correctly ABSENT as EXACTLY ZERO instead of at
    /// the Optima box's numerical lower bound. 0 = off (the default), 1 = on.
    ///
    /// WHY. Native and the Optima path disagree about what "absent" looks like
    /// in the answer, and the gap is twenty orders of magnitude. Native's
    /// line-search objective GX() truncates any trial amount below pa_DcMin
    /// (1e-33 on the psina projects) to exactly 0 - measured 1 113 007 times in
    /// one solve of 07PSIna_G_complex_1_0_1_80_0, and responsible for 609 of the
    /// 542 zeros in that answer BEFORE PhaseSelectionSpeciationCleanup() runs
    /// once (see the comment at that truncation, and plan-v5 section 60). A
    /// box-constrained method cannot do that: every species is held above
    /// dcFloor (pa_DHB, 1e-13 there) and can only approach it. So an absent
    /// species is reported at 1e-13 where native reports 0, and - the part that
    /// costs something - it remains a live unknown, with a chemical potential,
    /// a share of the convergence test, and the ability to throttle the shared
    /// step length (measured: 96-97% of all step throttling on two systems came
    /// from species sitting a hair above their floor).
    ///
    /// WHAT THIS DOES. Purely a post-solve reporting step: after every existing
    /// trustworthiness check has run and PASSED at full dimension, each species
    /// that is (a) not kinetically required to be present (DLL <= 0), (b) still
    /// sitting on the numerical floor rather than on a real constraint, and
    /// (c) correctly there by the solver's own KKT test (non-negative reduced
    /// gradient) - or forcibly excluded by a degenerate box - is set to exactly
    /// 0, and the derived output is recomputed from that state.
    ///
    /// SELF-GATING, which is what makes it safe: the zeroed state is put back
    /// through CheckMassBalanceResiduals() and is KEPT ONLY IF IT PASSES the
    /// same per-IC test the un-zeroed state just passed. It therefore cannot
    /// turn an accepted answer into a rejected one, and it cannot fire at all on
    /// a solve that was already going to be reported as BAD. The check is not a
    /// formality: dropping ~1e-13 from each of several hundred species perturbs
    /// the balance by an amount that is utterly negligible against a major IC
    /// and NOT negligible against a trace one - which is exactly the trace-IC
    /// sensitivity this branch has documented repeatedly.
    ///
    /// PHASE-LEVEL AND REBALANCED since 2026-09-15 (owner decisions, plan v5 s123.7). The self-gate above could
    /// not see a trace IC (CheckMassBalanceResiduals()' absolute cutoff is 1e-3 mol at pa_DHB = 1e-13), and the
    /// per-species zeroing mostly hit DISSOLVED species of the present aqueous phase, whose small amounts are
    /// real. Now: (1) a phase is zeroed as a whole only if every member passes the tests above - species of a
    /// present phase are never zeroed; (2) if that leaves any IC or the charge row past both its residual
    /// before zeroing and its tolerance (B_i*DHBM; charge: DHBM times the total charge carried), the removed
    /// amount is put back onto present carriers by MassBalanceReproject(); (3) if that fails, the zeroing is
    /// undone. DECIDE `zeroabsent phases= species= rebalanced= reverted= of=`.
    ///
    /// This does NOT reproduce native's mechanism, only its reported semantics.
    /// Native's species genuinely leave the problem mid-solve and stop
    /// obstructing it; these leave only the answer. Closing that half is the
    /// dimension reduction (pa_OptimaDimReduce), which omits species from the
    /// vector outright - but its final full-dimension verification pass puts
    /// them all back on the floor before the answer is written, which is the
    /// gap this field closes from the other end.
    ///
    /// WAS DEFAULT 1 from 2026-09-03 to 2026-09-15 (now 2 - see VALUE 2 below), set on the project owner's decision: matching
    /// native's reported semantics is worth having, and an absent species read
    /// as 0 rather than as 1e-13 is the more useful answer. Measured before the
    /// flip across Resources/gems3k + gems3k-fail, both arms per project:
    /// 43 projects, ZERO rows differing in status, iterations, G, pH, Eh, Vs or
    /// Ms; it fires on 41 of them (3 to 561 species); the mass-balance
    /// self-gate never once had to revert it; and the 301-point solvus sweep -
    /// which scores end-member MOLE FRACTIONS through a critical point, not a
    /// status - reproduces its record exactly while firing at every one of the
    /// 301 points. See plan-v5 section 61.
    ///
    /// One consequence to know when reading output: a consumer that divides by
    /// a species amount now meets 0 where it previously met 1e-13. Native has
    /// always written exact zeros, so anything reading native's output already
    /// copes; something written against this path's output specifically may
    /// not. Set to 0 to restore the floor-valued reporting.
    ///
    /// VALUE 2 - DEFAULT since 2026-09-15 (owner: "if a tiny amount is numerically better then don't use 0, only
    /// if 0 helps the solver"; plan v5 s123.8). Keep the amounts Optima returned - absent species stay at the
    /// floor, as Reaktoro reports them - and repair only a FAILING mass balance. Measured before building:
    /// (i) the zeros of value 1 never help an Optima leg, whose warm entry clamps x back to the floor, and they
    /// cost a native warm leg (07PSIna_G_simple_1 SHP: ITG 33 re-inserting from zeros, 0 from floor values);
    /// (ii) what value 1 bought over value 0 on 12 of 49 AOP projects was the MassBalanceReproject() its
    /// rebalance runs, which removed residuals the SOLVE had left, on major ICs too (f_/j_TestSUP98 H 111 mol at
    /// 18x its tolerance, f_CASHNK_G_Chen04C-3T Si 1 mol at 95x, 07PSIna_G_mid S at 280x) - not the zeros.
    /// So value 2: no zeroing; if any ordinary IC fails |C_i| <= B_i*DHBM on the accepted answer, project it
    /// back onto A.x = b (MassBalanceReproject) and KEEP the projection only if every ordinary IC then passes
    /// and the charge row is no worse than max(its residual before, DHBM x total charge carried); otherwise
    /// the answer is restored exactly as Optima returned it. DECIDE `optimarepair relbefore= relafter=
    /// chgbefore= chgafter= kept=`, written only when the repair was attempted.
    /// Values: 0 = floor values, nothing else (the RAW profile's value); 1 = zero absent phases + rebalance;
    /// 2 = floor values + bounded repair. Any other non-zero value behaves as 1.
    /// REACH: the fifteen T8/T14 exports (T8_aq*, T8ax2_nIC*, T14_ball*) pin 1 in their -ipm.json and keep
    /// value 1 until their project files change.
    long int OptimaZeroAbsent = 2;

    /// SEED A READMITTED SPECIES AT ITS PREDICTED AMOUNT instead of at the
    /// numerical floor, in OptimaReducedPreSolve()'s pricing loop. GEMS3K-only,
    /// keyword ipm-dat I/O. 0 (default) = off, readmission at the floor exactly
    /// as before. A POSITIVE value is a cap on the growth exponent, in RT units
    /// - see "the cap" below for why the knob is the exponent and not the
    /// amount.
    ///
    /// THE FORMULA IS White, Johnson & Dantzig (1958) Eq. 12a, the original RAND
    /// paper (Docs/literature/, and FINDINGS-White1958-RAND-vs-GEMS3K.md section 3).
    /// White derives it to compute trace species that were left out of the
    /// problem entirely, from the converged multipliers alone - his own example
    /// recovers NH3 at 2e-5 from a 10-species calculation that never contained
    /// it. In his notation x_i = xbar * exp[-c_i + sum_j a_ij pi_j]; here the
    /// multipliers are pm.U[] and the same statement is one line, because for a
    /// species whose chemical potential carries ln(x_j) with unit coefficient
    /// (DC_SYMMETRIC and DC_ASYM_SPECIES - see DC_PrimalChemicalPotential())
    /// the reduced gradient s_j = F[j] - sum_i U[i]*a(j,i) is the ONLY term
    /// that moves when x_j does, so setting s_j to zero gives
    ///
    ///     x_j^predicted = x_j^current * exp( -s_j ).
    ///
    /// Both inputs are already computed by the readmission loop, which prices
    /// exactly this s_j and readmits on its sign. So this is a use of a number
    /// the loop has in hand, not a new evaluation.
    ///
    /// WHY THIS IS NOT THE VARIANT ALREADY MEASURED AND REJECTED. That loop's
    /// own comment records interior seeding at "bulk-composition bound x 1e-6"
    /// costing 07PSIna_G_mid_1's pass 1 303 -> 2263 iterations with pass 2 then
    /// failing, and instructs that it not be re-tried "without a genuinely new
    /// hypothesis". 1e-6 is an arbitrary fraction; this is the amount the
    /// thermodynamics predicts. The recorded failure is evidence against
    /// arbitrary seeding, and does not bear on seeding at the predicted value.
    ///
    /// RESTRICTED BY SPECIES CLASS, and the exclusion is not fussiness.
    /// DC_SINGLE (a pure phase) has F = G with no x-dependence at all, so
    /// exp(-s_j) is not a prediction of anything - a pure phase's amount is set
    /// by the mass balance, not by its own potential, and s_j merely says
    /// whether it wants to exist. DC_ASYM_CARRIER (the solvent) carries logYF,
    /// which is itself a function of x_j, so the unit-coefficient argument does
    /// not hold for it either; it is in practice never omitted, since it is
    /// always in the LP seed's support. Both are left at the floor.
    ///
    /// THE CAP, and why the knob is the exponent. s_j is a difference of two
    /// potentials and on a large formula unit both are of order 1e3 (the
    /// standing seawater/Loeweite case: G0/RT = -7915.88 against sum(U*a) =
    /// -8107.22), so exp(-s_j) can overflow outright. Capping the exponent
    /// bounds the growth per pass in the units the quantity is actually
    /// expressed in - pa_OptimaReadmitSeed = 20 permits at most e^20 = 4.9e8 x
    /// the floor, i.e. ~5e-5 mol from a 1e-13 floor. Two further clamps apply
    /// unconditionally: the species' own box, and its stoichiometric ceiling
    /// min_i b_i/a(j,i) over the ICs it consumes - the same bound
    /// DetectPhaseCollapseAndReseed() already computes, and the only one of the
    /// three that is a statement about the system rather than about numerics.
    ///
    /// WHAT IT CAN AND CANNOT AFFECT. A readmitted species becomes a free
    /// variable inside its own box; this sets only where inside that box the
    /// NEXT pass starts from. It cannot change the reduced problem's solution,
    /// only the path to it - and the pre-solve's guarantee is unchanged either
    /// way, since a pass is handed over only at a fixed point and discarded
    /// otherwise. The risk it does carry is that a large seed leaves the
    /// starting point further from satisfying Aex*x = be than the floor did,
    /// which is the mechanism behind the 1e-6 failure above; that is what the
    /// cap is for, and why this ships off.
    ///
    /// MEASURED 2026-09-04, cold AOP, every 200+ species project in the corpus
    /// that the reduction is reachable on (both arms run in full; the six
    /// baselines all reproduce their recorded values exactly, so the A/B is
    /// clean). The answer is IDENTICAL in every row - same G to ten digits,
    /// same pH and Vs - so this only ever changes cost:
    ///
    ///   project            species   off   cap 20    change
    ///   07PSIna_G_mid_1       265     269    261      -3%
    ///   f_TestPNTDB           690     501    496      -1%
    ///   j_TestPNTDB           690     495    471      -5%
    ///   f_TestSUP98           923    2103   2008      -4.5%
    ///   j_TestSUP98           923    1695  10221     +503%   <- pass 1 discarded
    ///   07PSIna_G_complex_1  1392    1601   1671      +4.4%
    ///   07PSIna_G_edt_2      1566   11215    808     -93%   <- 13.9x, pass 1 no
    ///                                                       longer discarded
    ///
    /// THE MECHANISM WORKS EXACTLY AS THE FORMULA PREDICTS, which is worth
    /// separating from the verdict. On 07PSIna_G_mid_1 all 56 readmitted
    /// species are seeded and the passes that consume them get cheaper - pass 1
    /// 20 -> 15 iterations, pass 2 4 -> 1, i.e. the readmitted species start
    /// essentially AT the answer, which is precisely Eq. 12a's claim. On
    /// f_TestSUP98 pass 1 settles immediately (0 readmitted). The predicted
    /// amount is right.
    ///
    /// WHAT IT DOES NOT CONTROL is whether the next pass converges at all, and
    /// that dominates. j_TestSUP98's pass 1 is discarded WITH seeding and
    /// converges without it; 07PSIna_G_edt_2's is discarded WITHOUT seeding and
    /// converges with it. Same mechanism, opposite outcomes, and nothing
    /// observed predicts which - note in particular that f_TestSUP98 and
    /// j_TestSUP98 are the SAME chemistry differing only in the thermodynamic-
    /// data path, and they move in opposite directions (-4.5% against +503%).
    ///
    /// AND THE RESPONSE IS NON-MONOTONE IN THE CAP, this solver's now-familiar
    /// signature for a knob that should not be tuned. 07PSIna_G_mid_1 over
    /// caps 5/10/15/20/30/40/60 gives 276/275/262/261/268/268/268 against 269
    /// off - amplitude +-3%, no trend, settling once the cap stops binding and
    /// the stoichiometric ceiling takes over. Worse, f_TestPNTDB is 501 off,
    /// 7874 at cap 10 (pass 1 discarded) and 496 at cap 20: a SMALLER seed is
    /// the one that breaks it. That is the right way round for this mechanism -
    /// a truncated exponent is itself an arbitrary value, which is the very
    /// thing the 1e-6 variant failed for - so if this is used at all, use a cap
    /// large enough that the stoichiometric ceiling binds instead (>= 20 on
    /// this corpus), never a small one.
    ///
    /// THE DOWNSIDE IS BOUNDED BUT NOT SMALL. Both cliff cases were rescued by
    /// the pre-solve's own initial-set fallback (section 41) - no answer is
    /// lost, and G is unchanged - but the rescue costs a whole wasted pass.
    /// Under the standing decision rule that is a fallback without a strict
    /// improvement, so DEFAULT 0. The knob is BIMODAL rather than marginal -
    /// one 13.9x win, one 6x loss, five moves of +-5% - and both extremes are
    /// the same event in opposite directions. Per project: SET IT (= 20) on
    /// 07PSIna_G_edt_2, where it is the largest single-project win this field
    /// offers; do NOT set it on j_TestSUP98 or 07PSIna_G_complex_1; elsewhere
    /// do not set it blind, because given the f_/j_TestSUP98 split,
    /// establishing that it helps a project costs as much as it saves.
    /// See plan-v5 section 63 and
    /// Docs/literature/FINDINGS-White1958-RAND-vs-GEMS3K.md section 3.
    ///
    /// Blast radius, by construction: OptimaReducedPreSolve() runs only on a
    /// COLD AOP/SOP solve, only at or above the dimension-reduction size gate,
    /// never in ROP, never on a warm start (so never in HOP or SOP), and never
    /// with control conditions active.
    double OptimaReadmitSeed = 0.;

    /// Window (in IPM iterations) for the noise-stall test. 0 = off (the shipped
    /// default). When > 0, IPM convergence is accepted once, over a whole window,
    /// ALL FOUR of these hold:
    ///   - the energy pm.FX is flat to kIpmStallFXTol,
    ///   - the composition is flat: sum(X) and max(X) to kIpmStallCompTol, AND
    ///     every species BOTH to kIpmStallCompTol of the system total AND to
    ///     kIpmStallSpRel of its own largest amount over the window (a species
    ///     below kIpmStallSpNegl of the total is exempt from the second clause),
    ///   - pm.PCI increased on kIpmStallIncLo..kIpmStallIncHi of the steps, i.e. it
    ///     is BOUNCING rather than moving.
    ///
    /// WHY FOUR, AND WHY NONE OF THEM IS pa_DK. The IPM loop terminates on
    /// pm.PCI <= pm.DXM, and on part of this corpus that criterion stops carrying
    /// signal long before it is met: once the composition settles, pm.PCI is a
    /// difference of nearly-equal numbers and wanders in a band whose median sits
    /// ABOVE pm.DXM, so termination becomes a waiting time for a lucky draw. On
    /// f_Kaolinite the answer is final to 1e-11 by iteration 25 of 62-363 and
    /// P(PCI <= DXM) = 0.0063/iteration, predicting a mean of 184.8 against an
    /// observed 183.8 (plan v5 section 74).
    ///
    /// Every single signal was measured UNSAFE on its own:
    ///   energy flat alone     -> 115 % energy error on f_/j_TestSUP98   (74.9)
    ///   criterion flat alone  -> 2.9e-05 error on f_GEOTHERM            (80.2)
    ///   criterion bouncing    -> 11 % error, fires on 25 of 42          (80.2)
    ///   mass balance          -> UNUSABLE inside this loop: MBR establishes it and
    ///                            IPM then drifts from it monotonically by design,
    ///                            so it is never satisfied here          (80.5)
    /// An earlier version instead guarded on "pm.PCI is within 30x of pm.DXM",
    /// which works but ties the accepted accuracy to the tolerance being replaced -
    /// and flipping its default broke proposed.aop's "native G is
    /// settings-independent" assertion. The conjunction above removes that
    /// reference entirely, and that assertion now passes.
    ///
    /// Replayed on all 42 gems3k + gems3k-fail traces: 82.8 % of native's IPM
    /// iterations saved, worst |FX - FXfinal|/|FX| = 5.1e-10, firing on exactly the
    /// 9 noise-tail projects. The envelope is FLAT - window 15..60, comp tolerance
    /// 1e-6..1e-4 and bounce band 0.25..0.75 all give the same result - with one
    /// sharp cliff: an energy tolerance of 1e-7 costs 2e-03, so 1e-9 carries 100x
    /// margin. That flatness is what distinguishes this from the tuned knobs on
    /// this branch, whose responses are non-monotone.
    ///
    /// THE COMPOSITION TEST IS PER SPECIES, AND IT NEEDS BOTH NORMALISERS. An
    /// earlier version used only the aggregate sum(X)/max(X), which are blind to a
    /// REDISTRIBUTION at nearly constant total - a closing miscibility gap or an
    /// appearing phase - and flipping the default then failed solvus.critical and
    /// proposed.aop's crossing_T11. Measured while fixing it (plan v5 section 82):
    ///   per-PHASE on pm.XF[k]  -> fixes the BETWEEN-phase half only. The 301-point
    ///     solvus sweep goes from 5 non-convergences and 7 out-of-tolerance points
    ///     to 0 and 1. Motion WITHIN a phase (a solid solution's end-member
    ///     fractions) is invisible to it, so it is per species and not per phase.
    ///   total-relative alone   -> passes both solvus tests, fails crossing_T11.
    ///     T11's vestigial gas phase is 1.5e-6 OF THE TOTAL while moving 75 % of
    ///     ITSELF, so no total-relative threshold separates it from settled
    ///     rounding noise: at 1e-6 it misses the boundary, at 1e-8 it stops firing
    ///     at all (0.4 % of iterations saved against 79.7 %).
    ///   species-relative alone -> passes crossing_T11, fails solvus.native.
    /// The two clauses catch different things - large ABSOLUTE motion in a big
    /// species, and large RELATIVE motion in a small one - so both are required.
    /// The negligibility exemption is what makes the relative clause usable at all:
    /// without it a species resting at the 1e-13 floor wiggles by 100 % of itself
    /// and blocks acceptance for ever. Normalise on the window's MAXIMUM, not its
    /// newest value, or a species on its way OUT exempts itself as it vanishes.
    ///
    /// STILL DEFAULT OFF, but no longer because it fails. At the changed default the
    /// suite is 11/11 - the first time this field has passed its own gate - saving
    /// 79.2 % of native's IPM iterations at kIpmStallSpRel = 1e-2, and 78.3 % at the
    /// shipped 1e-3. Flipping it is an owner decision because it moves every recorded
    /// native iteration count in the benchmark corpus, not because it is unsafe.
    ///
    /// WHY 1e-3 AND NOT 1e-2. Measured envelope, combined test, 42 projects:
    ///   1e-4  passes, 0.0 % broad benefit (fires on 1 project)
    ///   1e-3  passes, 6.3 % broad         (fires on 3)   <- shipped, 30x margin
    ///   1e-2  passes, 9.9 % broad         (fires on 6, 3 become deterministic)
    ///   3e-2  FAILS
    /// (broad = excluding FeNaCl_FyGt_Precip_HighpH, which alone is 77 % of the
    /// corpus's native iterations because it sits at its 25000 cap; it goes to 111.)
    /// Benefit is monotone in this constant right up to a cliff, so 1e-2 sits 3x
    /// below a failure and 1e-3 sits 30x below one. The failure above the cliff is a
    /// wrong ANSWER at a phase boundary and the cost below it is only lost savings,
    /// so the asymmetry says err low. That one-decade-to-a-cliff shape is exactly
    /// what this branch elsewhere calls a tuned knob - unlike the window and the
    /// bounce band, whose envelopes are flat - so do not raise it without re-running
    /// the full suite at the changed default.
    short IpmStallWindow = 30;

    /// Repair an unsatisfied mass balance on the answer native is about to return,
    /// instead of adjusting what residual is acceptable. NATIVE path only.
    /// 0 = off (default), 1 = on.
    ///
    /// White, Johnson & Dantzig (1958), the note under their Table III: project Y
    /// back onto A.Y = b with one m x m solve over a chosen set of m species,
    /// dy_p = {(a_pj)^-1}(db_j). This is a third remedy for the ten projects whose
    /// native SIA cannot re-solve their own converged answer, and it differs in kind
    /// from the two already tracked - pa_DT (an absolute floor) and pa_MbClassRule
    /// both change the TEST; this changes the STATE.
    ///
    /// WHITE'S OWN PIVOT RULE - "the m most abundant species" - IS SINGULAR ON
    /// AQUEOUS CHEMISTRY AND IS NOT WHAT THIS IMPLEMENTS. Measured 2026-09-06 on the
    /// three projects whose mass-balance warning fires on every run: det(Ap) is
    /// EXACTLY zero on all three, because the most abundant species in an aqueous
    /// system are precisely the ones most likely to be exact stoichiometric sums of
    /// each other. The null combinations, named by the measurement:
    ///     H2O(l) - OH- - H+          (07PSIna_G_iron, Al-species)
    ///     Cl- + Na+ - NaCl@          (Al-species)
    ///     Na+ + OH- - NaOH(aq)       (FeNaCl_FyGt_Precip)
    /// White's own test case was a 10-species ideal gas mixture over 3 elements,
    /// where this does not arise.
    ///
    /// So the pivot set here is RANK-REVEALING: take species in decreasing amount and
    /// keep one only if it raises the rank. Measured against the two alternatives on
    /// the same three projects (cond / residual after / max |dy|/X):
    ///     minimum-norm, weight X     2e16-5e22 / 1e-11..9e-11 (INCOMPLETE) / 1e-8..9e-3
    ///     minimum-norm, weight sqrt(X) 3e10-6e12 / 1.4e-14      / 1e-3..0.40
    ///     rank-revealing pivot        8-48      / 0..1.4e-14    / 1e-3..0.14
    /// The rank-revealing basis is four orders better conditioned, removes the
    /// residual exactly, and touches only N species instead of 15-20. Note the
    /// minimum-norm variant weighted by amount - the physically appealing one, since
    /// it spreads the correction in proportion to what is there - fails: A.diag(X^2).A^T
    /// is the same near-singular matrix MBR assembles (water's H:O = 2:1), so the
    /// solve is too inaccurate to remove the residual it was computed for.
    ///
    /// Honest cost: an IC whose residual sits on a trace element can only be repaired
    /// by a trace species, so the largest RELATIVE move lands there - 14 % of Fe+2 on
    /// 07PSIna_G_iron, 3.9 % of H2@ on Al-species. In absolute terms those are ~1e-10
    /// mol and the energy change is ~1e-12 relative, but a speciation report will show
    /// them. That is why this is default-off rather than unconditional.
    short MbReproject = 1;

    /// pa_DeterminacyWarn: relative-uncertainty threshold for the "answer not fully
    /// determined by the energy" warning (TMultiBase::EnergyDeterminacyCheck(), native
    /// path, converged calls only). A present phase whose amount the minimised Gibbs
    /// energy fixes only to worse than this - sqrt(2 eps_G c_k)/n_k, eps_G the rounding
    /// floor of G, c_k the phase's least-energy compliance under A dn = 0 - is listed in
    /// one spdlog warning and a DECIDE `undetermined` trace record. 0 (or negative)
    /// switches the check off entirely, at zero cost. Read-only either way: it never
    /// changes the answer (identical on 71/71 runnable corpus projects, 2026-09-13).
    ///
    /// Default 1e-2. At that value it fires on 24 of 72 corpus projects, always on trace
    /// solids of 1e-9..1e-18 mol; validated against observed movement on
    /// 07PSIna_G_vcomplex_2_0_1_80_0 (predicted TiO2(am_hyd) +-20 % vs a 14 % trajectory
    /// move; CaSiO3(cr) +-2.6e-6 vs a 2.2e-6 jitter spread). Trailing member: GEMSGUI
    /// serialises BASE_PARAM positionally.
    double DeterminacyWarn = 1e-2;

    /// pa_ColdRetryNudges: recovery of a failed cold native call. When TNode::GEM_run() is asked for
    /// NEED_GEM_AIA (no kinetics) and ends in ERR_GEM_AIA, it re-solves cold at up to this many nudged bulk
    /// compositions, bIC[i] x (1 + s_i k 1e-15) with s_i alternating in sign over the ICs and
    /// k = +1, -1, +2, -2, ...; the first nudge that returns OK_GEM_AIA becomes the warm start of an SIA
    /// solve at the EXACT bIC, and that solve's outcome is returned as OK_GEM_AIA or BAD_GEM_AIA. The answer
    /// is an equilibrium at the requested composition, not the nudged one. If no attempt recovers, the
    /// original ERR_GEM_AIA is returned with its DATABR restored. ITF/ITG and IterDone count every attempt;
    /// each attempt leaves a DECIDE `coldretry` record, and the inner solves write no RUN/KEY records of
    /// their own. 0 = off.
    ///
    /// Basis (gems-benchmark tools/dilute_probe, HANDOFF-2026-09-14c s4): on T-cement's water sweep
    /// (50 steps x 5 nudges, native AIA) 4 cold solves fail; an SIA at the exact bIC started from a draw
    /// one 1e-15 nudge away that converged recovers all 4 in <= 2 iterations, a warm retry from the failed
    /// state recovers none, and a start from a solve with the trace ICs raised to 1e-6 mol recovers 2.
    /// Trailing member: GEMSGUI serialises BASE_PARAM positionally.
    long int ColdRetryNudges = 4;

    /// pa_OptimaPreSolveFirstIters: Optima iteration budget for each pass of the dimension-reduction
    /// pre-solve's FIRST attempt (OptimaReducedPreSolve() under the configured pa_OptimaDimReduceTol rule).
    /// A pass is otherwise allowed max(2000, pa_IIM) iterations; this lowers that to min(max(2000, pa_IIM),
    /// value) and can never raise it. The fallback attempt keeps the full budget. A pass that reaches the
    /// budget without converging is discarded exactly as before and the fallback attempt runs, so the worst
    /// case is the existing discard/retry path. 0 = off. Resolved by optima_presolve_pass_budget(), which
    /// the EFF trace line calls too; each `dimreducepass` DECIDE record carries the `budget=` it ran at.
    ///
    /// Basis (plan v5 s122.5, s122.7; gems-benchmark Docs/presolve-jitter-2026-09-12.tsv): over 5 projects
    /// x 9 draws of a 1e-15 bIC nudge (110 passes) the worst CONVERGED first-attempt pass costs 3868
    /// iterations (07PSIna_G_vcomplex_0_0_1_25_0) and every first-attempt pass that does not converge runs
    /// to the 10000 ceiling (T-cement on all 9 draws, T14_ball000 on 3). The second attempt does not
    /// separate (T-cement converges once in 8917), hence first attempt only. 6000 = 1.55x the worst
    /// converged pass; 5000 (1.29x) sits inside the +37 % by which one project's worst pass moved between
    /// its single census draw and the nine-draw ladder. Trailing member: GEMSGUI serialises BASE_PARAM
    /// positionally.
    long int OptimaPreSolveFirstIters = 6000;

    /// pa_LpDualFillout: size the species the LP zeroed from the LP's OWN DUAL instead of from the
    /// per-class constants pa_DFYaq/DFYw/DFYid/DFYh/DFYr/DFYc. Native cold (AIA) path only, at
    /// DC_RaiseZeroedOff()'s call site in GEM_IPM_InitialApproximation(); HOP/SHP's native leg
    /// inherits it. 0 = OFF and is the shipped default.
    ///
    /// x_j = X_k * exp( a_j^T u_LP - G_j ) in RT units - Karpov 1997 Eq. 12, the same formula
    /// pa_OptimaReadmitSeed uses at the Optima pre-solve's re-admission site. X_k is the LP's own
    /// amount of the species' phase.
    ///
    /// MEASURED AND REJECTED AS A DEFAULT, 2026-09-25 (plan v5 137.3; FABLE Phase 3 WP5, CLOSED).
    /// Kept switchable on owner instruction the same day - "might come back in the future" - NOT
    /// because the measurements were inconclusive. They were not:
    ///   value 1  loses THREE answers (f_/j_GEOTHERM at the iteration cap, f_Solvus_G_test3) and
    ///            costs +14.5 % ITF / +14.2 % ITG over the 72 projects OK in both arms. The amount
    ///            it writes is closer to the converged one than the class constant on 457 of 12 579
    ///            species. Variance is extreme: 34.0x worse on j_TiQ_PRSV and 0.13x on
    ///            FeNaCl_FyGt_TransitionZone.
    ///   value 2  the composition ceiling dominates the class floor: loses TWENTY answers. Offered
    ///            only so 137.4's veto stays reproducible - do not use it on real work.
    ///   value 3  value 1 applied only within 8 RT of the leveling hyperplane, where 137.8a measured
    ///            the dual to carry information: WORSE than 1 (ITG 1.619x). Also measured, also kept.
    /// WHY IT CANNOT WORK AS IT STANDS, and what a future attempt must fix first: an amount is the
    /// exponential of a potential, and the LP's dual is a median 29.3 RT from the converged one
    /// (WP4) - about 13 decades of amount. That gap is NOT a linearisation error. Re-solving the LP
    /// at the CONVERGED potentials leaves it at 26.5 RT (137.8c), because the LP fixes its objective
    /// and not its dual: on 63 of 76 projects two duals tens of RT apart price the bulk identically
    /// to better than 1e-3. So a future attempt needs a BETTER-DETERMINED dual, not a better formula.
    /// Trailing member: GEMSGUI serialises BASE_PARAM positionally.
    long int LpDualFillout = 0;

    /// pa_FilloutBudget: cap how much the class fill-out may perturb the MASS BALANCE, as a
    /// FRACTION of each element's own bulk amount. Native cold (AIA) path only, right after
    /// DC_RaiseZeroedOff(); HOP/SHP's native leg inherits it. **DEFAULT 0 = OFF.**
    ///
    /// > IT WAS MEANT TO SHIP ON (owner, 2026-09-25, "ship it") AND IT DID NOT, because `ctest` came
    /// back RED at 0.01 and CLAUDE.md s2 says a red gate is not overridden by the ship rule:
    /// **ci.baseline 261 failing checks - 224 ANSWER, 33 COST, 3 STATUS, 1 pH** - plus proposed.aop's
    /// physics check `crossing_T11`. The answer moves are small (Vs ~1e-8 relative, e.g. f_GEOTHERM
    /// 1.594511804e-01 -> 1.594511917e-01) but the scored tolerance is 1e-9, so they count.
    /// THE OPEN QUESTION THAT DECIDES IT (plan v5 138.10): native's stopping test is tolerance-limited
    /// and 20 of 42 projects are count-unstable, so a start change redraws where the loop stops -
    /// these 224 rows may be that known noise rather than a worse answer. Compare each move against
    /// THAT PROJECT'S OWN jitter spread before calling it a regression, as CLAUDE.md s5 requires.
    /// Until that is done the field stays 0.
    /// A CORPUS SWEEP SCORING ONLY G AT 1e-9 FOUND ONE MOVED ROW; ci.baseline scores Vs too and found
    /// 224. The sweep was too narrow - which is why the suite caught this and the sweep did not.
    ///
    /// WHY A JOINT BUDGET AND NOT A PER-SPECIES CAP. DC_RaiseZeroedOff() raises every species the
    /// LP zeroed to a per-class constant, and those constants routinely ask for more of an element
    /// than the system contains. Measured 2026-09-25 over 74 projects: EVERY project over-subscribes
    /// at least one element, by a median of 48.9x and a maximum of 6.68e+09x, and the count of
    /// over-subscribed elements scales with species count (median 60 elements on projects of 200+
    /// species, 6 below that). Capping each species at PhaseInsertionCeiling() - the most the bulk
    /// could supply to THAT species alone - only reaches a median of 3.25x and leaves 52 of 72
    /// projects over, because the constraint being violated is a SUM over species and a per-species
    /// bound cannot enforce a sum (plan v5 138.5).
    ///
    /// THE RULE. The LP solution satisfies A n = b exactly, so every bit of the excess comes from
    /// the raise. With raised_i = sum_j (Y_j - Y_lp_j) * a(j,i), scale each raised species by the
    /// tightest element it consumes, s_j = min_i( f * B_i / raised_i ), capped at 1. Because
    /// raised_i is itself the sum over species, the per-species minimum GUARANTEES
    /// sum_j s_j * raised_ji <= f * B_i - the joint bound holds by construction, not by tuning.
    /// Measured worst ratio comes back at exactly 1 + f to six decimals on every project.
    /// Dimensionless: f is a fraction of the project's own bulk, so it cannot become the fixed
    /// absolute amount that is plan v5 section 76's defect.
    ///
    /// MEASURED AT 0.01, 74 projects, native cold AIA, one draw each (plan v5 138.8):
    /// ITF (mass balance) 571 -> 314, 0.550x, better on 42 projects and worse on 3; ITG (IPM)
    /// 0.984x and a wash per project (24 better, 28 worse, 22 same) - so the win is the
    /// MASS-BALANCE STAGE, not the solve. One G row moves (07PSIna_G_vcomplex_0 @ 80 C) and one
    /// status row moves, f_Solvus_G_test3, which is a coin flip at the SHIPPED settings (4 of 9 OK
    /// under a 1e-15 nudge, Eh spread 1.68 V) and is owner-deferred on that basis (138.6).
    ///
    /// WHY 0.01 - TUNED over six decades, 2026-09-25 (plan v5 138.9), all three corpora, native
    /// cold AIA, one draw per project, against f = 0 as the baseline:
    ///
    ///     f       answers lost   G moved   ITF ratio   ITG ratio   ITF better/worse
    ///     1e-6         3            7        0.657       1.010         52/10
    ///     1e-4         1            6        0.585       0.975         44/4
    ///     1e-3         1            3        0.564       0.948         41/4
    ///     1e-2         1            1        0.550       0.984         42/3     <- default
    ///     1e-1         1            2        0.529       0.991         45/0
    ///     1            1            3        0.538       1.026         42/0
    ///
    /// "TIGHTER IS ALWAYS BETTER" IS FALSE, and the sweep is what shows it: at f = 1e-6 the budget
    /// starves the start and loses THREE answers instead of one, with ITG turning worse. The trend
    /// reverses below about 1e-4. An earlier note here claimed the monotone reading from two points
    /// an order apart; it was wrong and this table replaces it.
    /// 1e-4 ... 1e-1 is a PLATEAU - every value there loses the same single answer (f_Solvus_G_test3,
    /// a coin flip at the shipped settings, 138.6) and cuts ITF by 41-47 %. Within that plateau 0.01
    /// has the FEWEST moved G rows (1, against 3 at 1e-3 and 6 at 1e-4), which is what the ship rule
    /// scores after answers; 1e-3 is marginally cheaper overall (total iterations 0.929 vs 0.963) and
    /// is the value to revisit if cost ever outranks row identity. The default is therefore inside a
    /// measured plateau, not at an untested edge.
    /// GEMS3K_FILLOUT_BUDGET overrides this field for a throwaway arm (negative = unset, so 0
    /// stays a usable arm); that is how the default was chosen and how it should be re-chosen.
    /// Trailing member: GEMSGUI serialises BASE_PARAM positionally.
    double FilloutBudget = 0.;
    /// pa_StabTPD: the tangent-plane (TPD) stability scan of ABSENT multicomponent phases, reported on
    /// the trace's CERT line (FABLE Phase 3 WP6; plan v5 s139.2; owner decision 2026-09-28).
    ///   0  off - nothing computed, the CERT fields read stab_ss=1e300 stab_ph=- stab_n=0 stab_dis=0.
    ///   1  compute and REPORT (default). Runs only where the CERT record runs, i.e. only when
    ///      GEMS3K_NATIVE_TRACE_FILE is set, after the answer has been packed - a production call
    ///      executes none of it, and the scan cannot change what the call returns.
    ///   2  (reserved) gate the certificate on it - NOT implemented; plan s20 owner decision, open.
    /// For each absent non-aqueous multicomponent phase: Michelsen successive substitution
    /// y <- W(y)/sum W(y), W_j = exp(a_j.u - G0_j - fDQF_j - lnGam_j(y)), from every vertex, the ideal
    /// closed form, the current composition, the centroid and a grid (binary) or edge midpoints; stab_ss =
    /// min over phases of the smallest TPD(y) = sum_j y_j (mu_j(y) - a_j.u) seen at ANY composition the
    /// search evaluated (< 0: the phase would lower G - the answer left a phase out). stab_dis counts
    /// phases the solver's own single-point index (pm.Falp <= 0) calls stable while stab_ss < 0 for them.
    /// CORRECTED 2026-09-28: the first version scored only converged stationary points and missed Al2O3-SiO2
    /// rs_ss (TPD -0.016..-0.23 on native's answer, 1100-1900 K) - plan v5 s139.7.
    /// Trailing member: GEMSGUI serialises BASE_PARAM positionally. RAW value: 1 (report-only).
    long int StabTPD = 1;
    /// pa_IpmAugmentedKKT: how the MAIN IPM loop solves its linear system for the dual u
    /// (MakeAndSolveSystemOfLinearEquations(), initAppr = false only; MBR is untouched).
    ///   0  off - normal equations (A_act^T W A_act) u = A_act^T W F, Cholesky then LU (the default until 2026-09-30).
    ///   1  augmented (saddle-point) system, dense LU of size L_act + N:
    ///          [ I          -W A_act ] [ x ]   [ -W F ]
    ///          [ -A_act^T   -D       ] [ u ] = [  0   ]
    ///      Ported from the owner's experimental SolverType == 2 (tmp/ipm_main.cpp, 2026-09-14).
    ///      Eliminating x gives (A^T W A + D) u = A^T W F exactly, so the answer is the normal
    ///      equations' up to D; what changes is that A^T W A is never FORMED, which squares the
    ///      condition number of W^1/2 A. COST: dense O((L_act+N)^3) per iteration - ~1e9 flops per
    ///      iteration on a 1392-species project; use 2 there.
    ///   2  (DEFAULT since 2026-09-30, owner: "on by default if no way to smart detect") the same regularised
    ///      least-squares problem, min |W^1/2 (A u - F)|^2 + u^T D u, by
    ///      Householder QR of [W^1/2 A_act ; D^1/2] - also never forms A^T W A, O(L_act N^2).
    /// D = 0 except on a ZERO ROW. The owner's version put an ABSOLUTE 1e-12 on every row (its size
    /// then depends on pa_DG rescaling and on W), and kept only species with Y > 1e-12 mol where the
    /// assembly keeps Y > min(lowPosNum, DcMinM); the active set here is the assembly's. A uniform
    /// D was measured harmful - see the comment in SolveIpmAugmented() (f_Kaolinite: 1e-12 relative
    /// per row gives 3963 iterations and a G off by 7e-6; D = 0 gives 45 against native's 273).
    /// ZERO-ROW RESCUE: an IC no active species carries has d_i = sum_j W_j a_ji^2 = 0; the normal
    /// equations then fail (E07IPM). Both arms here set D_ii = 1e-12 * max_i d_i instead, which sets that u_i to 0 and
    /// continues - the owner's "the determinant no longer goes to zero". That HIDES a degeneracy rather
    /// than repairing it, so it is recorded: DECIDE "ipmkkt-zerorow" on the first rescue of each
    /// InteriorPointsMethod() call.
    /// FALLBACK: where the augmented solve finds its system singular (QR rank test 1e-14, or LU), that
    /// step is taken by the normal equations instead (DECIDE "ipmkkt-fallback", first per call) - so
    /// the switch can never fail where the default would have gone on. Added after the 2026-09-28c
    /// standard freeze lost three warm SIA answers (07PSIna_G_simple_1/2, CASH+CsSr) to QR declaring
    /// rank deficiency on the first warm step; with it all three are OK again (probe, 2026-09-29).
    /// MEASURED 2026-09-28b, SMOKE ONLY (mode_compare native/AOP/HOP, 27 gems3k projects, one draw, NOT a
    /// freeze): 2 -> native 0 answers lost, iterations x0.994, G moved only on f_/j_CASHNK (1.8e-7 / 8e-9
    /// relative, a known jitter project); 1 -> 0 lost, native x0.960, HOP x1.165 all from j_CASHNK's Optima
    /// leg (4526 -> 6526). Per-project counts move both ways on the lottery projects (f_Kaolinite 276 -> 48,
    /// j_Kaolinite 58 -> 155 under 2), i.e. a redraw. Mode 2's solve was verified to rounding (normal-
    /// equation residual <= 6e-16 on f_Kaolinite). Whether any case needed the zero-row rescue: NOT checked.
    /// DEFAULT FLIP 2026-09-30: gate 2026-09-29-combined STD-OFF -> STD-ON: 0 lost, 1 won, 13 changes <= 2e-6 rel + two
    /// lower-G assemblages, 0.987x; RAW 0 lost, 9 won, 0.851x. A "smart" trigger (QR only where Cholesky of the normal
    /// matrix fails) was probed and rejected: it fired on 1 of 18 projects and kept 1 of the 9 RAW gains - the gains come
    /// from the QR solve's accuracy in general, not from rescuing a failed factorisation.
    /// Trailing member: GEMSGUI serialises BASE_PARAM positionally. RAW value: 0.
    long int IpmAugmentedKKT = 2;
    /// pa_IpmLoopTweaks: bit mask of the owner's three SolverType == 2 main-loop changes, split so each
    /// can be measured alone (default 0 = none). Independent of pa_IpmAugmentedKKT.
    ///   1  step cap - StepSizeEstimate()'s LM clamped to <= 1 before OptimizeStepSize().
    ///   2  activity lag - once pm.PCI < 5e-4, CalculateActivityCoefficients() runs only on every third
    ///      ITG, and the loop may not terminate on an iteration that skipped it (the lnGam in hand
    ///      would not be the one the composition implies). 5e-4 is ABSOLUTE, as in the original.
    ///   4  loose accept - after ITG > 120, accept pm.PCI < 300 * pm.DXM as converged. This LOOSENS
    ///      the stopping test by up to 300x; it reports a state the Dikin test has not certified.
    ///      DECIDE "ipmlooseaccept" when it fires.
    /// MEASURED 2026-09-28b, SMOKE ONLY (same set and harness as pa_IpmAugmentedKKT): 1 -> 0 lost, native
    /// iterations x1.506 (t_Solvus series2 410 -> 1555), HOP x1.220; 2 -> native x1.091, HOP LOSES
    /// j_CASHNK (OK -> FAIL); 4 -> native LOSES 4 (f_/j_TestPNTDB, o_/t_Solvus series2, OK -> FAIL), x0.770
    /// on the rest; 7 with pa_IpmAugmentedKKT = 2 -> native loses 2 (TestPNTDB), x1.662, HOP x2.150. None
    /// is a candidate default; the switch exists so the owner's variant stays reproducible.
    /// Trailing member: GEMSGUI serialises BASE_PARAM positionally. RAW value: 0.
    long int IpmLoopTweaks = 0;
    /// pa_OptimaLineSearch: Optima's merit line search on the UNMASKED error, with this trigger factor
    /// (a step whose error exceeds factor x the previous one is line-searched). 0 = off (Optima's own default,
    /// the behaviour before the field). Reaches Optima::Options::linesearch in both the reduced pre-solve and
    /// the full solve. REQUIRES the local Optima fix in ErrorControl::execute (E updated at the new point
    /// before the comparison) - without it the trigger compares the error with itself and never fires
    /// (plan v5 s139.6). Measured with that fix, AOP cold, 77 projects: factor 1.5 loses 0 answers and
    /// 0 phases, wins T-cement (native's G to 9 digits), median ITG 0.93x, total +8 %. Ships 1.5 (owner,
    /// 2026-09-28). RAW value 0. Trailing member: GEMSGUI serialises BASE_PARAM positionally.
    double OptimaLineSearch = 1.5;
    /// pa_OptimaFDDiagFloor: in the pa_OptimaFDHessian block, which OVERWRITES a basic variable's whole column
    /// (diagonal included) with a finite difference, put the analytic diagonal back where the FD diagonal is
    /// not positive. A basic variable below the amount at which PrimalChemicalPotentials() recomputes F gets
    /// an FD diagonal of exactly 0; with a non-zero coupling that makes the reduced Hessian indefinite by
    /// construction (plan v5 s139.5: resolved negative curvature on 2 of 83 projects, all of it this).
    /// 0 = off (behaviour before the field; ships off, owner 2026-09-28), 1 = on. RAW value 0.
    /// DECIDE fddiagfloor reports the count.
    long int OptimaFDDiagFloor = 0;
    /// pa_OptimaLSStallEscape (2026-09-30, owner's "option 2"): with the line search on, after this many CONSECUTIVE line
    /// searches that leave the error unchanged (relative change <= 1e-8) keep the full step once. 0 = off. Needs an Optima
    /// built with OPTIMA_LINESEARCH_STALL_ESCAPE (optima/install-ls2); against an older install it is inert. Measured on
    /// probes: f_TestPNTDB AOP 1827 -> 493 it at 10 (line search off: 501); T-cement's rescue, 07PSIna_G_edt_2's 656-it speed-up
    /// and j_Solvus unchanged; Cu-Pourbaix's crawl (483) is not a freeze and is not caught - looser tolerances (1e-3, 1e-2)
    /// caught it but lost T-cement. Trailing member: GEMSGUI serialises BASE_PARAM positionally. RAW value: 0.
    /// DEFAULT 10 since 2026-09-30 (owner: "agree"), pending the escape freeze pair; corium diagrams identical to 10 = off
    /// within 0.1 % on all five x {AOP, SHP}.
    long int OptimaLSStallEscape = 10;
    /// pa_OptimaLSWindow (2026-09-30, owner's "option 1", kept opt-in): the line-search trigger compares with the max of the
    /// last N pre-step errors instead of the previous one. 0 = off (default). Helped Cu-Pourbaix (AOP 483 -> 110 it at 10) but
    /// was WORSE than the strict trigger on every corium phase diagram (Al2O3-SiO2 AOP converged 89.9 % -> 38.8 % at 10) and lost
    /// T-cement's rescue at every N - not recommended for phase diagrams. Needs OPTIMA_LINESEARCH_STALL_ESCAPE, as above.
    /// Trailing member: GEMSGUI serialises BASE_PARAM positionally. RAW value: 0.
    long int OptimaLSWindow = 0;
    /// pa_OptimaLSRejectWorse (2026-09-30, gems3k-6f): a line search that ends at or above the pre-step error is discarded
    /// and the full step kept (Optima LineSearchOptions::reject_if_worse). 0 = off (DEFAULT since 2026-10-01, owner: "off by
    /// default if only special difficult projects are affected"), 1 = on; always off in TNode::GEM_run()'s legacy retry.
    /// RC freeze 2026-09-30: clean in STANDARD (0 lost/changed, AOP 1.144x) but MASKED HARM raw (2 fake OKs CSHSnplus mb_rel
    /// 1.767, f_BaSrCarbonate_a0wide AOP/SOP loses native's answer). Recommend per project where the line search crawls. Measured (AOP, line search 1.5, escape 10, objective memo): Cu-Pourbaix fired 341 line searches, 329 ending
    /// ABOVE the pre-step error - a crawl the stall escape cannot see; with the rule 489 -> 149 it, j_Solvus 153 -> 104,
    /// 07PSIna_G_edt_2 687 -> 639, f_TestPNTDB 463 unchanged, same G; solvus sweep 61/61 (6388 -> 6564 it); corium diagrams
    /// within 0.4 % on answers; T-cement OK via the retry. Needs OPTIMA_LINESEARCH_REJECT_WORSE; inert against an older Optima.
    /// Trailing member: GEMSGUI serialises BASE_PARAM positionally. RAW value: 0.
    long int OptimaLSRejectWorse = 0;
    /// pa_OptimaTpdAccept (PROTOTYPE, plan v5 section 140.15 A): when Optima ends not converged / KKT / stability failed
    /// but mass balance holds, accept the state iff every off-bound species is stationary, every at-bound species of a
    /// single-DC/aqueous/gas phase has gradJ >= -kktTol, and every ABSENT non-ideal condensed phase has a direct-TPD search
    /// minimum >= -value (THERMOCHIMICA's phase-level criterion). 0 = off. Measured value 1e-6. RAW value: 0.
#ifndef GEMS3K_DEFAULT_OPTIMA_TPDACCEPT
#define GEMS3K_DEFAULT_OPTIMA_TPDACCEPT 1e-6   // release default 2026-10-01 (RC freeze: clean standard + raw)
#endif
    double OptimaTpdAccept = GEMS3K_DEFAULT_OPTIMA_TPDACCEPT;
    /// pa_OptimaCgSeed (PROTOTYPE, plan v5 section 140.15 B): cold Optima seed by column generation (species Gibbs-LP +
    /// TPD-priced pseudo-compound columns, THERMOCHIMICA's Leveling/PEA); value = TPD tolerance; 0 = off (feasibility-LP
    /// seed). Measured value 1e-6. RAW value: 0.
#ifndef GEMS3K_DEFAULT_OPTIMA_CGSEED
#define GEMS3K_DEFAULT_OPTIMA_CGSEED 1e-6   // release default 2026-10-01 (RC freeze: clean standard + raw)
#endif
    double OptimaCgSeed = GEMS3K_DEFAULT_OPTIMA_CGSEED;
    /// pa_OptimaColdRetry (PROTOTYPE, plan v5 section 140.16 C/C'): a warm Optima call (SOP, SHP) that is not OK is
    /// re-solved cold (AOP) by TNode::GEM_run_optima_cold_retry(); 0 = off, 1 = retry after the full warm budget,
    /// 2 = fail fast (skip the warm call's full-budget re-solves, then retry). Measured value 2. RAW value: 0.
#ifndef GEMS3K_DEFAULT_OPTIMA_COLDRETRY
#define GEMS3K_DEFAULT_OPTIMA_COLDRETRY 2   // release default 2026-10-01 (RC freeze: clean standard + raw)
#endif
    long int OptimaColdRetry = GEMS3K_DEFAULT_OPTIMA_COLDRETRY;
    /// pa_OptimaFinish (PROTOTYPE, plan v5 section 142): when an Optima call ends not converged, a Newton finish on the FIXED
    /// phase set (species amounts and multipliers, equality-constrained, line-searched on G) is run from Optima's last
    /// primal by TMultiBase::PotentialSpaceFinish(); its result is then judged by the same KKT / mass-balance / TPD checks
    /// as Optima's own. 0 = off, 1 = on. RAW value: 0.
#ifndef GEMS3K_DEFAULT_OPTIMA_FINISH
#define GEMS3K_DEFAULT_OPTIMA_FINISH 1   // release default 2026-10-01 (RC freeze: clean standard + raw)
#endif
    long int OptimaFinish = GEMS3K_DEFAULT_OPTIMA_FINISH;
    /// pa_OptimaAcceptRepair (2026-10-01, gems3k-6f; owner: off by default, on per project): when pa_OptimaTpdAccept's
    /// acceptance passes every test except the per-IC relative mass balance (mb_rel > 1), apply MassBalanceReproject() - the
    /// repair the success path already applies (OptimaZeroAbsent = 2, DECIDE optimarepair) - and re-test; restored if it does
    /// not bring mb_rel <= 1 (DECIDE tpdaccept-repair). Needs pa_OptimaTpdAccept > 0. 0 = off (default), 1 = on.
    /// Measured on T-cement (with pa_GAS = 2e-3 in the project): AOP/SOP OK at 11-13 g H2O (were FAIL), G = native's to
    /// 8-10 digits; 10.0 g still fails (no attempt reaches native's answer). Trailing member (GEMSGUI). RAW value: 0.
    long int OptimaAcceptRepair = 0;

    void write(GemDataStream& oss);
    void read(GemDataStream& iss);
};

// ---------------------------------------------------------------------------
// EFFECTIVE values of AUTO-gated settings
// ---------------------------------------------------------------------------
//
// Some pa_* fields are three-valued: a configured 0 does not mean "off", it
// means "decide from the problem". For those, the value in the project file,
// in the calculation trace's SET line and in a benchmark freeze's `# set` line
// is NOT the value the solver ran at, and no settings audit can see the
// difference. That cost a full measurement (plan v5 section 95.4): the whole
// T14 ballast ladder was measured against pa_OptimaPhaseCompaction while every
// rung silently ran 8 dimension-reduction passes, which is what made that field
// look like a no-op - correctly, but for a reason nothing in the record showed.
//
// So the resolution rule lives HERE, in one place, and both the solver and
// native_trace_run_header()'s EFF line call it. They cannot disagree about what
// ran, and a change to the gate shows up in the trace on the next run.
//
// RULE FOR ANY NEW AUTO-GATED FIELD (owner, 2026-09-07c): put its resolution in
// a function like this one, call it from the solver rather than inlining the
// constants, and add it to the EFF line. A knob whose effective value is not in
// the trace is a knob that cannot be measured.

/// Species-count gate above which pa_OptimaDimReduce = 0 (AUTO) turns the
/// dimension-reduction pre-solve ON. Interpolated, not measured - nothing in the
/// corpus lies between j_GEOTHERM (154 species, a 2.7x loss) and 07PSIna_G_mid_1
/// (265 species, a 60x win) - which is why it is a named constant with an
/// explicit off switch rather than a silent hardcode. See BASE_PARAM::OptimaDimReduce.
/// The "total Gibbs energy not computed yet" marker that MultiConstInit() seeds pm.FX with
/// (ipm_simplex.cpp), named here because until 2026-09-22 it was a bare 7777777. under a comment
/// reading "???????" and was being PUBLISHED as an answer.
///
/// The VALUE is deliberately absurd and that is the good part of the design: a converged total
/// Gibbs energy on this corpus is of order -5e3 J, so +7.777777e6 cannot be mistaken for one by
/// eye, and its digits are greppable. Initialising to 0. instead would have looked plausible and
/// hidden the defect indefinitely.
///
/// What was missing is that NOTHING TESTED IT. pm.FX is refreshed by the native descent only, so
/// on every Optima-family mode this marker travelled out through packDataBr() as CNode->Gs - and
/// therefore as GEM_to_MT()'s p_Gs and as the Gs field of every exported -dbr file - while
/// TNode::Get_GibbsEnergy(), a recomputation, returned the right number. A marker whose whole
/// purpose is to say "nobody computed this" is only worth having if something asks; packDataBr()
/// now does. The Optima path also sets pm.FX properly (ipm_optima.cpp), so the guard should be
/// silent - if it ever speaks, a path is returning an answer it never priced.
constexpr double kTotalGibbsEnergyUnset = 7777777.;

constexpr long int kOptimaDimReduceAutoMinDC  = 200;
/// Pass count AUTO selects when the gate opens.
constexpr long int kOptimaDimReduceAutoPasses = 8;

/// Resolves pa_OptimaDimReduce's three-valued setting to the pass count the
/// solver will actually attempt: > 0 is an explicit count, < 0 is explicitly
/// off, 0 is AUTO. `nDC` is the species count (pm.L).
///
/// NOTE this is the pass count only. Whether the pre-solve is REACHED also
/// depends on the call: it runs on a cold Optima leg (pm.pNP == 0), and on a
/// HOP leg only under an EXPLICIT positive setting - AUTO deliberately does not
/// reach HOP's warm verification, which is already O(1) there. See the call site
/// in ipm_optima.cpp.
inline long int optima_dimreduce_passes( long int configured, long int nDC )
{
    if( configured == 0 && nDC >= kOptimaDimReduceAutoMinDC )
        return kOptimaDimReduceAutoPasses;
    return configured;
}

/// Per-pass Optima iteration budget of OptimaReducedPreSolve(): max(2000, pa_IIM),
/// lowered to pa_OptimaPreSolveFirstIters on the FIRST attempt when that is > 0.
/// Never raises the budget. One function for the solver and the EFF trace line.
inline long int optima_presolve_pass_budget( long int configured, long int iim, bool firstAttempt )
{
    const long int full = ( iim > 2000L ) ? iim : 2000L;
    return ( firstAttempt && configured > 0 && configured < full ) ? configured : full;
}


/// Cap that pa_OptimaEarlyStabilityAt = 0 (AUTO) resolves to when the system
/// carries a MULTISITE (sublattice) solid-solution model. 200, and it is not
/// interpolated: it is the value the whole field was gated at on 2026-09-09c
/// (412 rows, 3 corpora, freeze_diff exit 0), and every measurement of this
/// field since section 101 has been taken at it.
constexpr long int kOptimaEarlyStabilityAutoCap = 200;

/// Cap AUTO resolves to on a WARM Optima leg (pm.pNP != 0), where the cold value is
/// pure waste. Measured, plan v5 section 105.4 / 106: on a warm leg the dual gate is
/// INERT - the dual is bit-exactly unchanged at the cap because a warm start inherits
/// a consistent one and nothing has happened yet - so the cap always fires, and a probe
/// that finds nothing costs EXACTLY N. N is therefore the loss, and the warm legs that
/// WIN win at every N tried (f_CASHNK SOP 1001 unarmed -> 201 at N=200 -> 26 at N=25).
constexpr long int kOptimaEarlyStabilityAutoWarmCap = 25;

/// DEFAULT for the early-probe safety net's re-solve strategy, and the single
/// place it is resolved - the solver and the trace both call this, so a knob
/// whose effective value is not in the trace cannot exist here (CLAUDE.md s4,
/// the rule pa_OptimaDimReduce's AUTO gate cost a whole ladder measurement).
///
/// 0  re-solve the net from `initialState` - the behaviour shipped before
///    2026-09-11, so a probe that finds nothing costs exactly N and the run
///    total is `N + full` (plan v5 s105.3).
/// 1  resume from the probe's own end state. LOSES ANSWERS - f_Solvus_G_test3
///    armed at -5, AOP, 9 nudges: 5 of 9 draws against the restart's 8 (s108.4).
/// 2  resume, and if the resumed re-solve does not converge, re-solve once from
///    `initialState` exactly as 0 would have. Bounded worst case, one extra
///    solve on a run that was already failing.
/// 3  as 2, but the fallback is deferred to the END of the retry ladder, where
///    the question "did this call, with every rescue it has, still fail?" is
///    answerable - so the same bounded worst case is paid on strictly fewer
///    rows. On f_CASHNK the phase-extinction tier converges a run whose resumed
///    re-solve stalled, and mode 2 pays 1000 iterations for a guard that buys
///    nothing there; on f_Solvus_G_test3 the tier does NOT rescue it and the
///    extra solve is the answer. Both are `netresolve outcome=stalled`, so the
///    outcome cannot separate them and only the end of the ladder can (s108.6).
///
/// DEFAULT 3 since 2026-09-11, on an owner decision under CLAUDE.md s2's
/// "measured improvements ship ON", bracketed by a freeze either side.
/// GEMS3K_OPTIMA_NET_RESUME overrides it, and 0 restores the old behaviour.
constexpr int kOptimaNetResumeDefault = 3;

/// Effective resume mode for this process. Reads the env override ONCE - the
/// value cannot change within a run, and both callers want the same answer.
inline int optima_net_resume_mode()
{
    static const int mode = []() -> int {
        const char* e = std::getenv( "GEMS3K_OPTIMA_NET_RESUME" );
        return ( e && *e ) ? std::atoi( e ) : kOptimaNetResumeDefault;
    }();
    return mode;
}

/// Is this phase's built-in mixing model a MULTISITE (sublattice) solid solution?
/// Exactly three codes, from m_const_base.h's own comments - Berman/Brown,
/// CALPHAD CEF, and the Modified Bragg-Williams model of Vinograd et al. 2018.
/// Everything else (Van Laar 'V', Guggenheim 'K', Redlich-Kister 'G', Margules,
/// the fluid EoS set, every aqueous model) is single-site or not a solid solution.
inline bool optima_smod_is_multisite( char code )
{
    return code == SM_BERMAN || code == SM_CEF || code == SM_MBW;
}

/// Number of multicomponent phases in the system using a multisite model.
/// `sMod` is MULTI::sMod (char[8] per phase, code in position 0) over FIs.
inline long int optima_multisite_phase_count( char (*sMod)[8], long int FIs )
{
    long int n = 0;
    if( !sMod ) return 0;
    for( long int k = 0; k < FIs; k++ )
        if( optima_smod_is_multisite( sMod[k][0] ) ) n++;
    return n;
}

/// Resolves pa_OptimaEarlyStabilityAt's three-valued setting to the value the
/// solver will actually arm: > 0 is an explicit cap at that iteration, < 0 is
/// the trend form with |value| consecutive falls, and 0 is AUTO - the cap at
/// kOptimaEarlyStabilityAutoCap when the system carries a multisite solid
/// solution, and off otherwise.
///
/// WHY MULTISITE IS THE GATE (owner's proposal, 2026-09-10; plan v5 section 104).
/// The cap's wins are not distributed over the corpus - they are one mechanism.
/// A sublattice model generates end-member sets in which two solution phases can
/// be built from identical stoichiometry at identical G0, i.e. thermodynamically
/// interchangeable, and an interior-point method can only approach the vanishing
/// twin's bound asymptotically (section 102.9(c)). Stopping early is what lets the
/// phase-selection and extinction tiers see that before the budget is spent.
/// MEASURED on the 2026-09-09c freeze, splitting its 19 moved rows by this exact
/// predicate: the 9 rows in multisite projects are net -2198 iterations with a
/// WORST SINGLE ROW of +200 (the probe, which is its bound by construction),
/// while the 10 rows elsewhere are net +1338 and contain every bad row the
/// corpus-wide flip had - a warm restart going 1 -> 570, and a project that fails
/// in both arms moving +9200. The gate keeps the mechanism and drops the noise.
///
/// WHAT IT FORFEITS, stated rather than hidden: 07PSIna_G_simple_0_0_0_150_0 AOP
/// 1032 -> 201, a real win on a project with NO solid solution at all (SIT
/// aqueous + ideal). Its DECIDE record shows a different mechanism - the
/// phase-selection repair loop on a marginal phase, not an interchangeable twin -
/// so it is a second use for the field, not evidence against this gate. That
/// project can still pin the cap explicitly.
///
/// NOTE the evidence base is 5 projects of 68, three of them CASH-family. The
/// predicate is chosen because it names the MECHANISM, not because five points
/// determine a rule.
///
/// WHY AUTO IS ALSO LEG-DEPENDENT (owner decision 2026-09-10, plan v5 section 106).
/// The multisite gate alone makes 4 of its 30 rows cheaper and 5 DEARER by exactly
/// +200, and all five of those are WARM legs paying the probe and finding nothing.
/// The cause is structural: the dual-settled gate is a real test on a cold leg and
/// INERT on a warm one (the dual is bit-exactly unchanged at the cap, because a warm
/// start inherits a consistent dual and nothing has happened yet), so on a warm leg
/// the cap always fires and a losing probe costs exactly N. Making AUTO answer 25
/// there turns four of those +200 into +25 and the fifth (CSHSnplus HOP) into 0,
/// while the warm WINS are untouched because a warm winner wins at every N.
/// This is what AUTO is FOR - "decide from the problem" - and it needs no new
/// BASE_PARAM field and no fourth meaning on this one, which is why it was preferred
/// over both (see DECISIONS 2026-09-10, decision 3).
///
/// There is deliberately no "explicitly off" encoding: < 0 is the trend form and
/// 0 is now AUTO. A caller who wants no early stop on a multisite system sets a
/// cap larger than the iteration budget (pa_IIM), which is self-describing and
/// needs no magic value.
/// `warmLeg` is pm.pNP != 0 - a call whose Optima leg starts from a consistent
/// (primal, dual) pair rather than from a cold LP seed. AUTO returns the short cap
/// there; an EXPLICIT setting is left flat on both legs, because explicit means
/// explicit and a caller who writes 200 is entitled to get 200.
inline long int optima_earlystability_at( long int configured, long int nMultisitePhases,
                                          bool warmLeg )
{
    if( configured == 0 && nMultisitePhases > 0 )
        return warmLeg ? kOptimaEarlyStabilityAutoWarmCap : kOptimaEarlyStabilityAutoCap;
    return configured;
}


typedef struct
{  // MULTI is base structure to Project (local values)
    char
    stkey[EQ_RKLEN+5],   ///< Record key identifying IPM minimization problem
    // NV_[MAXNV], nulch, nulch1, ///< Variant Nr for fixed b,P,T,V; index in a megasystem
    PunE,         ///< Units of energy  { j;  J c C N reserved }
    PunV,         ///< Units of volume  { j;  c L a reserved }
    PunP,         ///< Units of pressure  { b;  B p P A reserved }
    PunT;         ///< Units of temperature  { C; K F reserved }
    long int
    N,        	///< N - number of IC in IPM problem
    NR,       	///< NR - dimensions of R matrix
    L,        	///< L -   number of DC in IPM problem
    Ls,       	///< Ls -   total number of DC in multi-component phases
    LO,       	///< LO -   index of water-solvent in IPM DC list
    PG,       	///< PG -   number of DC in gas phase
    PSOL,     	///< PSOL - number of DC in liquid hydrocarbon phase
    Lads,     	///< Total number of DC in sorption phases included into this system.
    FI,       	///< FI -   number of phases in IPM problem
    FIs,      	///< FIs -   number of multicomponent phases
    FIa,      	///< FIa -   number of sorption phases
    FI1,     ///< FI1 -   number of phases present in eqstate
    FI1s,    ///< FI1s -   number of multicomponent phases present in eqstate
    FI1a,    ///< FI1a -   number of sorption phases present in eqstate
    IT,      ///< It - number of completed IPM iterations
    E,       ///< PE - flag of electroneutrality constraint { 0 1 }
    PD,      ///< PD - mode of calling CalculateActivityCoefficients() { 0 1 2 3 4 }
    PV,      ///< Flag for the volume balance constraint (on Vol IC) - for indifferent equilibria at P_Sat { 0 1 }
    PLIM,    ///< PU - flag of activation of DC/phase restrictions { 0 1 }
    Ec,     ///< CalculateActivityCoefficients() return code: 0 (OK) or 1 (error)
    K2,     ///< Number of IPM loops performed ( >1 up to 6 because of PSSC() )
    PZ,     ///< Indicator of PSSC() status (since r1594): 0 untouched, 1 phase(s) inserted
    ///< 2 insertion done after 5 major IPM loops
    pNP,    ///< Mode of FIA selection: 0-automatic-LP AIA, 1-smart SIA, -1-user's choice
    pESU,   ///< Unpack old eqstate from EQSTAT record?  0-no 1-yes
    pIPN,   ///< State of IPN-arrays:  0-create; 1-available; -1 remake
    pBAL,   ///< State of reloading CSD:  1- BAL only; 0-whole CSD
    tMin,   ///< Type of thermodynamic potential to minimize
    pTPD,   ///< State of reloading thermod data: 0-all  -1-full from database   1-new system 2-no
    pULR,   ///< Start recalc kinetic constraints (0-do not, 1-do )internal
    pKMM, ///< new: State of KinMet arrays: 0-create; 1-available; -1 remake
    ITaia,  ///< Number of IPM iterations completed in AIA mode (renamed from pRR1)
    FIat,   ///< max. number of surface site types
    MK,     ///< IPM return code: 0 - continue;  1 - converged
    W1,     ///< Indicator ofSpeciationCleanup() status (since r1594) 0 untouched, -1 phase(s) removed, 1 some DCs inserted
    is,     ///< is - index of IC for IPN equations ( CalculateActivityCoefficients() )
    js,     ///< js - index of DC for IPN equations ( CalculateActivityCoefficients() )
    next,   ///< for IPN equations (is it really necessary? TW please check!
    sitNcat,    //< Can be re-used
    sitNan,     // Can be re-used
    SolveCallCount ///< Number of calls to MakeAndSolveSystemOfLinearEquations() during this GEM call (reset per call)
    ;
    double
    TC,  	///< Temperature T, min. (0,2000 C)
    TCc, 	///< Temperature T, max. (0,2000 C)
    T,   	///< T, min. K
    Tc,   	///< T, max. K
    P,      ///< Pressure P, min(0,10000 bar)
    Pc,   	///< Pressure P, max.(0,10000 bar)
    VX_,    ///< V(X) - volume of the system, min., cm3
    VXc,    ///< V(X) - volume of the system, max., cm3
    GX_,    ///< Gibbs potential of the system G(X), min. (J)
    GXc,    ///< Gibbs potential of the system G(X), max. (J)
    AX_,    ///< Helmholtz potential of the system F(X)
    AXc,    ///<  reserved
    UX_,  	///< Internal energy of the system U(X)
    UXc,  	///<  reserved
    HX_,    ///< Total enthalpy of the system H(X)
    HXc, 	///<  reserved
    SX_,    ///< Total entropy of the system S(X)
    SXc,   ///<	 reserved
    CpX_,  ///< reserved
    CpXc,  ///< 20 reserved
    CvX_,  ///< reserved
    CvXc,  ///< reserved
    // TKinMet stuff
    kTau,  ///< current time, s (kinetics)
    kdT,   ///< current time step, s (kinetics)

    TMols,      ///< Input total moles in b vector before rescaling
    SMols,      ///< Standart total moles (upscaled) {1000}
    MBX,        ///< Total mass of the system, kg
    FX,    	    ///< Current Gibbs potential of the system in IPM, moles
    IC,         ///< Effective molal ionic strength of aqueous electrolyte
    pH,         ///< pH of aqueous solution
    pe,         ///< pe of aqueous solution
    Eh,         ///< Eh of aqueous solution, V
    DHBM,       ///< balance (relative) precision criterion
    DSM,        ///< min amount of phase DS
    GWAT,       ///< used in ipm_gamma()
    YMET,       ///< reserved
    PCI,        ///< Current value of Dikin criterion of IPM convergence DK>=DX
    CondNum,    ///< Estimated 2-norm condition number of the IPM linear system A (worst case over this GEM call, power/inverse iteration on the Cholesky/LU factors), 0 if not yet computed
    CondNumDiag,///< Cheap proxy: max/min |diagonal| ratio of the IPM linear system A (worst case over this GEM call), 0 if not yet computed
    SolveTimeMs,   ///< Wall-clock time (ms), summed over all calls to MakeAndSolveSystemOfLinearEquations() during this GEM call (matrix build + decomposition + solve + condition-number diagnostics)
    CondNumTimeMs, ///< Wall-clock time (ms), summed over all calls, spent specifically computing the condition-number diagnostics (subset of SolveTimeMs) — isolates the cost of that instrumentation from the rest of the linear solve
    DXM,        ///< IPM convergence criterion threshold DX (1e-5)
    lnP,        ///< log Ptotal
    RT,         ///< RT: 8.31451*T (J/mole/K)
    FRT,        ///< F/RT, F - Faraday constant = 96485.309 C/mol
    Yw,         ///< Current number of moles of solvent in aqueous phase
    ln5551,     ///< ln(55.50837344)
    aqsTail,    ///< v_j asymmetry correction factor for aqueous species
    lowPosNum,  ///< Minimum mole amount considered in GEM calculations (MinPhysAmount = 1.66e-24)
    logXw,      ///< work variable
    logYFk,     ///< work variable
    YFk,        ///< Current number of moles in a multicomponent phase
    FitVar[5];  ///< Internal. FitVar[0] is total mass (g) of solids in the system (sum over the BFC array)
    ///<      FitVar[1], [2] reserved
    ///<       FitVar[4] is the AG smoothing parameter;
    ///<       FitVar[3] is the actual smoothing coefficient
    double
    denW[5],   ///< Density of water, first T, second T, first P, second P derivative for Tc,Pc
    denWg[5],  ///< Density of steam for Tc,Pc
    epsW[5],   ///< Diel. constant of H2O(l)for Tc,Pc
    epsWg[5];  ///< Diel. constant of steam for Tc,Pc

    long int
    *L1,    ///< l_a vector - number of DCs included into each phase [Fi]
    // TSolMod stuff
    *LsMod, ///< Number of interaction parameters. Max parameter order (cols in IPx),
    ///< and number of coefficients per parameter in PMc table [3*FIs]
    *LsMdc, ///<  for multi-site models: [3*FIs] - number of nonid. params per component;
    /// number of sublattices nS; number of moieties nM
    *LsMdc2, ///<  new: [3*FIs] - number of DQF coeffs; reciprocal coeffs per end member;
    /// reserved
    *IPx,   ///< Collected indexation table for interaction parameters of non-ideal solutions
    ///< ->LsMod[k,0] x LsMod[k,1]   over FIs
    *mui,   ///< IC indices in RMULTS IC list [N]
    *muk,   ///< Phase indices in RMULTS phase list [FI]
    *muj,   ///< DC indices in RMULTS DC list [L]

    *LsPhl,  ///< new: Number of phase links; number of link parameters; [Fi][2]
    (*PhLin)[2];  ///< new: indexes of linked phases and link type codes (sum 2*LsPhl[k][0] over Fi)

    /* TSolMod !! arrays and counters to be added (for mixed-solvent electrolyte phase) TW

  ncsolv, /// TW new: number of solvent parameter coefficients (columns in solvc array)
  nsolv,  /// TW new: number of solvent interaction parameters (rows in solvc array)
  *ixsolv, /// new: array of indexes of solvent interaction parameters [nsolv*2]
  *solvc, /// TW new: array of solvent interaction parameters [ncsolv*nsolv]

  ncdiel, /// TW new: number of dielectric constant coefficients (colums in dielc array)
  ndiel,  /// TW new: number of dielectric constant parameters (rows in dielc array)
  *ixdiel /// new: array of indexes of dielectric interaction parameters [ndiel*2]
  *dielc, /// TW new: array of dielectric constant parameters [ncdiel*ndiel]

  ndh,    /// TW new: number of generic DH coefficients (rows in dhc array)
  *dhc,   /// TW new: array of generic DH parameters [ndh]
  */

    // TSorpMod stuff
    long int
    *LsESmo, ///< new: number of EIL model layers; EIL params per layer; CD coefs per DC; reserved  [Fis][4]
    *LsISmo, ///< new: number of surface sites; isotherm coeffs per site; isotherm coeffs per DC; max.denticity of DC [Fis][4]
    *xSMd;   ///< new: denticity of surface species per surface site (site allocation) (-> L1[k]*LsISmo[k][3]] )
    long int  (*SATX)[4]; ///< Setup of surface sites and species (will be applied separately within each sorption phase) [Lads]
    /// link indexes to surface type [XL_ST]; sorbent em [XL_EM]; surf.site [XL-SI] and EDL plane [XL_SP]
    // TKinMet stuff
    long int
    *LsKin,  ///< new: number of parallel reactions nPRk[k]; number of species in activity products nSkr[k];
    /// number of parameter coeffs in parallel reaction term nrpC[k]; number of parameters
    /// per species in activity products naptC[k]; nAscC number of parameter coefficients in As correction;
    /// nFaces[k] number of (separately considered) crystal faces or surface patches ( 1 to 4 ) [Fi][6]
    *LsUpt,  ///< new: number of uptake kinetics model parameters (coefficients) numpC[k];
    /// number of IC element indexes for end members = L1[k]    [Fis][2]
    *xSKrC,  ///< new: Collected array of aq/gas/sorption species indexes used in activity products (-> += LsKin[k][1])
    (*ocPRkC)[2], ///< new: Collected array of operation codes for kinetic parallel reaction terms (-> += LsKin[k][0])
    /// and indexes of faces (surface patches)
    *xICuC;  ///< new: Collected array of IC species indexes used in partition (fractionation) coefficients  ->L1[k]   TBD
    double
    // TSolMod stuff
    *PMc,    ///< Collected interaction parameter coefficients for the (built-in) non-ideal mixing models -> LsMod[k,0] x LsMod[k,2]
    *DMc,    ///< Non-ideality coefficients f(TPX) for DC -> L1[k] x LsMdc[k][0]
    *MoiSN,  ///< End member moiety- site multiplicity number tables ->  L1[k] x LsMdc[k][1] x LsMdc[k][2]
    *SitFr,  ///< Tables of sublattice site fractions for moieties -> LsMdc[k][1] x LsMdc[k][2]

    // Stoichiometry basis
    *A,   ///< DC stoichiometry matrix A composed of a_ji [0:N-1][0:L-1]
    *Awt,    ///< IC atomic (molar) mass, g/mole [0:N-1]

    // Reconsider usage
    *Wb,     ///< Relative Born factors (HKF, reserved) [0:Ls-1]
    *Wabs,   ///< Absolute Born factors (HKF, reserved) [0:Ls-1]
    *Rion,   ///< Ionic or solvation radii, A (reserved) [0:Ls-1]
    *HYM__,    ///< reserved
    *ENT__,    ///< reserved no object

    *H0,     ///< DC pmolar enthalpies, reserved [L]
    *A0,     ///< DC molar Helmholtz energies, reserved [L]
    *U0,     ///< DC molar internal energies, reserved [L]
    *S0,     ///< DC molar entropies, reserved [L]
    *Cp0,    ///< DC molar heat capacity, reserved [L]
    *Cv0__,    ///< DC molar Cv, reserved [L]

    *VL,      ///< ln mole fraction of end members in phases-solutions
    // Old sorption stuff
    *Xcond,   ///< conductivity of phase carrier, sm/m2   [0:FI-1], reserved
    *Xeps,    ///< diel.permeability of phase carrier (solvent) [0:FI-1], reserved
    *Aalp,    ///< Full vector of specific surface areas of phases (m2/g) [0:FI-1]
    *Sigw,    ///< Specific surface free energy for phase-water interface (J/m2)   [0:FI-1]
    *Sigg,  	///< Specific surface free energy for phase-gas interface (J/m2) (not yet used)  [0:FI-1], reserved
    // from here move to --> datach.h
    // TSolMod stuff
    *lPhc,  ///< new: Collected array of phase link parameters (sum(LsPhl[k][1] over Fi)
    *DQFc,  ///< new: Collected array of DQF parameters for DCs in phases -> L1[k] x LsMdc2[k][0]
    //  *rcpc,  ///< new: Collected array of reciprocal parameters for DCs in phases -> L1[k] x LsMdc2[k][1]

    // TSorpMod & TKinMet stuff
    *SorMc, ///< new: Phase-related kinetics and sorption model parameters: [Fis][16]
    ///< in the same order as from Asur until fRes2 in TPhase
    // TSorpMod stuff
    *EImc,  ///< new: Collected EIL model coefficients k -> += LsESmo[k][0]*LsESmo[k][1]
    *mCDc,  ///< new: Collected CD EIL model coefficients per DC k -> += L1[k]*LsESmo[k][2]
    *IsoPc, ///< new: Collected isotherm coefficients per DC k -> += L1[k]*LsISmo[k][2];
    *IsoSc, ///< new: Collected isotherm coeffs per site k -> += LsISmo[k][0]*LsISmo[k][1];
    // TKinMet stuff
    *feSArC, ///< new: Collected array of fractions of surface area related to parallel reactions k-> += LsKin[k][0]
    *rpConC,  ///< new: Collected array of kinetic rate constants k-> += LsKin[k][0]*LsKin[k][2];
    *apConC,  ///< new:!! Collected array of parameters per species involved in activity product terms
    ///  k-> += LsKin[k][0]*LsKin[k][1]*LsKin[k][3];
    *AscpC,   /// new: parameter coefficients of equation for correction of specific surface area k-> += LsKin[k][4]
    *UMpcC  ///< new: Collected array of uptake model coefficients k-> += L1[k]*LsUpt[k][0];
    ;
    // until here move to --> datach.h

    //  Data for old surface comlexation and sorption models (new variant [Kulik,2006])
    double  (*Xr0h0)[2];   ///< mean r & h of particles (- pores), nm  [0:FI-1][2], reserved
    double  (*Nfsp)[MST];  ///< Fractions of the sorbent specific surface area allocated to surface types  [FIs][FIat]
    double  (*MASDT)[MST]; ///< Total maximum site  density per surface type (mkmol/g)  [FIs][FIat]
    double  (*XcapF)[MST]; ///< Capacitance density of Ba EDL layer F/m2 [FIs][FIat]
    double  (*XcapA)[MST]; ///< Capacitance density of 0 EDL layer, F/m2 [FIs][FIat]
    double  (*XcapB)[MST]; ///< Capacitance density of B EDL layer, F/m2 [FIs][FIat]
    double  (*XcapD)[MST]; ///< Eff. cap. density of diffuse layer, F/m2 [FIs][FIat]
    double  (*XdlA)[MST];  ///< Effective thickness of A EDL layer, nm [FIs][FIat], reserved
    double  (*XdlB)[MST];  ///< Effective thickness of B EDL layer, nm [FIs][FIat], reserved
    double  (*XdlD)[MST];  ///< Effective thickness of diffuse layer, nm [FIs][FIat], reserved
    double  (*XlamA)[MST]; ///< Factor of EDL discretness  A < 1 [FIs][FIat], reserved
    double  (*Xetaf)[MST]; ///< Density of permanent surface type charge (mkeq/m2) for each surface type on sorption phases [FIs][FIat]
    double  (*MASDJ)[DFCN];  ///< Parameters of surface species in surface complexation models [Lads][DFCN]
    // Contents defined in the enum below this structure
    // Other data
    double
    *XFs,    ///< Current quantities of phases X_a at IPM iterations [0:FI-1]
    *Falps,  ///< Current Karpov criteria of phase stability  F_a [0:FI-1]
    *Fug,    ///< Demo partial fugacities of gases [0:PG-1]
    *Fug_l,  ///< Demo log partial fugacities of gases [0:PG-1]
    *Ppg_l,  ///< Demo log partial pressures of gases [0:PG-1]

    *DUL,     ///< VG Vector of upper kinetic restrictions to x_j, moles [L]
    *DLL,     ///< NG Vector of lower kinetic restrictions to x_j, moles [L]
    *fDQF,    ///< Increments to molar G0 values of DCs from pure gas fugacities or DQF terms, normalized [L]
    // TKinMet stuff (old DODs, new contents )
    *PUL,  ///< Vector of upper restrictions to multicomponent phases amounts [FIs]
    *PLL,  ///< Vector of lower restrictions to multicomponent phases amounts [FIs]
    *PfFact, /// new: phase surface area - volume shape factor (taken from TKinMet or set from TNode) [FI]
    *PrT,    /// new: Total MWR rate (mol/s) for phases - TKinMet output [FI]
    *PkT,    /// new: Total specific MWR rate (mol/m2/s) for phases - TKinMet output [FI]
    *PvT,    /// new: Total one-dimensional MWR surface propagation velocity (m/s) - TKinMet output [FI]
    //  potentially can be extended to all solution phases?
    *emRd,   /// new: output Rd values (partition coefficients) for end members (in uptake kinetics model) [Ls]
    *emDf,   /// new: output Df values (fractionation coeffs.) for end members (in uptake kinetics model) [Ls]
    //
    *YOF,     ///< Surface free energy parameter for phases (J/g) (to accomodate for variable phase composition) [FI]
    *Vol,     ///< DC molar volumes, cm3/mol [L]
    *MM,      ///< DC molar masses, g/mol [L]
    *Pparc,   ///< Partial pressures or fugacities of pure DC, bar (Pc by default) [0:L-1]
    *Y_m,     ///< Molalities of aqueous species and sorbates [0:Ls-1]
    *Y_la,    ///< log activity of DC in multi-component phases (mju-mji0) [0:L-1]
    *Y_w,     ///< Mass concentrations of DC in multi-component phases,%(ppm)[Ls]
    *Gamma,   ///< DC activity coefficients in molal or other phase-specific scale [0:L-1]
    *lnGmf,   ///< ln of initial DC activity coefficients for correcting G0 [0:L-1]
    *lnGmM,   ///< ln of DC pure gas fugacity (or metastability) coeffs or DDF correction [0:L-1]
    *EZ,      ///< Formula charge of DC in multi-component phases [0:Ls-1]
    *FVOL,    ///< phase volumes, cm3 comment corrected DK 04.08.2009  [0:FI-1]
    *FWGT,    ///< phase (carrier) masses, g                [0:FI-1]
    //
    *G,       ///< Normalized DC energy function c(j), mole/mole [0:L-1]            --> activities.h
    *G0,      ///< Input normalized g0_j(T,P) for DC at unified standard scale[L]   --> activities.h
    *lnGam,   ///< ln of DC activity coefficients in unified (mole-fraction) scale [0:L-1] --> activities.h
    *lnGmo;   ///< Copy of lnGam from previous IPM iteration (reserved)
    double  (*lnSAC)[4]; ///< former lnSAT ln surface activity coeff and Coulomb's term  [Lads][4]

    // TSolMod stuff (detailed output on partial energies of mixing)   --> activities.h
    double *lnDQFt; ///< new: DQF terms adding to overall activity coefficients [Ls_]
    double *lnRcpt; ///< new: reciprocal terms adding to overall activity coefficients [Ls_]
    double *lnExet; ///< new: excess energy terms adding to overall activity coefficients [Ls_]
    double *lnCnft; ///< new: configurational terms adding to overall activity [Ls_]
    // TSorpMod stuff
    double *lnScalT;  ///< new: Surface/volume scaling activity correction terms [Ls_]
    double *lnSACT;   ///< new: ln isotherm-specific SACT for surface species [Ls_]
    double *lnGammF;  ///< new: Frumkin or BET non-electrostatic activity coefficients [Ls_]
    double *CTerms;   ///< new: Coulombic terms (electrostatic activity coefficients) [Ls_]

    double  *B,  ///< Input bulk chem. compos. of the system - b vector, moles of IC[N]
    *U,  ///< IC chemical potentials u_i (mole/mole) - dual IPM solution [N]
    *U_r,  ///< IC chemical potentials u_i (J/mole) [0:N-1]
    *C,    ///< Calculated IC mass-balance deviations (moles) [0:N-1]
    *IC_m, ///< Total IC molalities in aqueous phase (excl.solvent) [0:N-1]
    *IC_lm,	///< log total IC molalities in aqueous phase [0:N-1]
    *IC_wm,	///< Total dissolved IC concentrations in g/kg_soln [0:N-1]
    *BF,    ///< Output bulk compositions of multicomponent phases bf_ai[FIs][N]
    *BFC,   ///< Total output bulk composition of all solid phases [1][N]
    *XF,    ///< Output total number of moles of phases Xa[0:FI-1]
    *YF,    ///< Approximation of X_a in the next IPM iteration [0:FI-1]
    *XFA,   ///< Quantity of carrier in asymmetric phases Xwa, moles [FIs]
    *YFA,   ///< Approximation of XFA in the next IPM iteration [0:FIs-1]
    *Falp;  ///< Karpov phase stability criteria F_a [0:FI-1] or phase stability index (PC==2)

    double (*VPh)[MIXPHPROPS],     ///< Volume properties for mixed phases [FIs]
    (*GPh)[MIXPHPROPS],     ///< Gibbs energy properties for mixed phases [FIs]
    (*HPh)[MIXPHPROPS],     ///< Enthalpy properties for mixed phases [FIs]
    (*SPh)[MIXPHPROPS],     ///< Entropy properties for mixed phases [FIs]
    (*CPh)[MIXPHPROPS],     ///< Heat capacity Cp properties for mixed phases [FIs]
    (*APh)[MIXPHPROPS],     ///< Helmholtz energy properties for mixed phases [FIs]
    (*UPh)[MIXPHPROPS];     ///< Internal energy properties for mixed phases [FIs]

    // old sorption models - EDL models (data for electrostatic activity coefficients)
    double (*XetaA)[MST]; ///< Total EDL charge on A (0) EDL plane, moles [FIs][FIat]
    double (*XetaB)[MST]; ///< Total charge of surface species on B (1) EDL plane, moles[FIs][FIat]
    double (*XetaD)[MST]; ///< Total charge of surface species on D (2) EDL plane, moles[FIs][FIat]
    double (*XpsiA)[MST]; ///< Relative potential at A (0) EDL plane,V [FIs][FIat]
    double (*XpsiB)[MST]; ///< Relative potential at B (1) EDL plane,V [FIs][FIat]
    double (*XpsiD)[MST]; ///< Relative potential at D (2) plane,V [FIs][FIat]
    double (*XFTS)[MST];  ///< Total number of moles of surface DC at surface type [FIs][FIat]
    //
    double *X,  ///< DC quantities at eqstate x_j, moles - primal IPM solution [L]
    *Y,   ///< Copy of x_j from previous IPM iteration [0:L-1]
    *XY,  ///< Copy of x_j from previous loop of Selekt2() [0:L-1]
    *Qp,  ///< Work variables related to non-ideal phases FIs*(QPSIZE=180)
    *Qd,  ///< Work variables related to DC in non-ideal phases FIs*(QDSIZE=60)
    *MU,  ///< mu_j values of differences between dual and primal DC chem.potentials [L]
    *EMU, ///< Exponents of DC increment to F_a criterion for phase [L]
    *NMU, ///< DC increments to F_a criterion for phase [L]
    *W,   ///< Weight multipliers for DC (incl restrictions) in IPM [L]
    *Fx,  ///< Dual DC chemical potentials defined via u_i and a_ji [L]
    *Wx,  ///< Mole fractions Wx of DC in multi-component phases [L]
    *F,   /// <Primal DC chemical potentials defined via g0_j, Wx_j and lnGam_j[L]
    *F0;  ///< Excess Gibbs energies for (metastable) DC, mole/mole [L]
    // Old sorption models
    double (*D)[MST];  ///< Reserved; new work array for calc. surface act.coeff.
    // Name lists
    char (*sMod)[8];   ///< new: Codes for built-in mixing models of multicomponent phases [FIs]
    char (*kMod)[6];  ///< new: Codes for built-in kinetic models [Fi]
    char  (*dcMod)[6];   ///< Codes for PT corrections for dependent component data [L]
    char  (*SB)[MAXICNAME+MAXSYMB]; ///< List of IC names in the system [N]
    char  (*SB1)[MAXICNAME]; ///< List of IC names in the system [N]
    char  (*SM)[MAXDCNAME];  ///< List of DC names in the system [L]
    char  (*SF)[MAXPHNAME+MAXSYMB];  ///< List of phase names in the system [FI]
    char  (*SM2)[MAXDCNAME];  ///< List of multicomp. phase DC names in the system [Ls]
    char  (*SM3)[MAXDCNAME];  ///< List of adsorption DC names in the system [Lads]
    char  *DCC3;   ///< Classifier of DCs involved in sorption phases [Lads]
    char  (*SF2)[MAXPHNAME+MAXSYMB]; ///< List of multicomp. phase names in the syst [FIs]
    char  (*SFs)[MAXPHNAME+MAXSYMB]; ///< List of phases currently present in non-zero quantities [FI]
    char  *pbuf, 	///< Text buffer for table printouts
    // Class codes
    *RLC,   ///< Code of metastability constraints for DCs [L] enum DC_LIMITS
    *RSC,   ///< Units of metastability/kinetic constraints for DCs  [L]
    *RFLC,  ///< Classifier of restriction types for XF_a [FIs]
    *RFSC,  ///< Classifier of restriction scales for XF_a [FIs]
    *ICC,   ///< Classifier of IC { e o h a z v i <int> } [N]
    *DCC,   ///< Classifier of DC { TESKWL GVCHNI JMFD QPR <0-9>  AB  XYZ O } [L]
    *PHC;   ///< Classifier of phases { a g f p m l x d h } [FI]
    char  (*SCM)[MST]; ///< Classifier of built-in electrostatic models applied to surface types in sorption phases [FIs][FIat]
    char  *SATT,  ///< Classifier of applied SACT equations (isotherm corrections) [Lads]
    *DCCW;  ///< internal DC class codes [L]
    // TSorpMod stuff
    char *IsoCt; ///< new: Collected isotherm and SATC codes for surface site types k -> += 2*LsISmo[k][0]

    //  SolutionData *asd; ///< Array of data structures to pass info to TSolMod [FIs]

    long int ITF,       ///< Number of completed IA EFD iterations
    ITG,         ///< Number of completed GEM IPM iterations
    ITau,    /// new: Time iteration for TKinMet class calculations
    IRes1;
    clock_t t_start, t_end;
    double t_elap_sec;  ///< work variables for determining IPM calculation time
    double *Guns;     ///<  mu.L work vector of uncertainty space increments to tp->G + sy->GEX
    double *Vuns;     ///<  mu.L work vector of uncertainty space increments to tp->Vm
    double *tpp_G;    ///< Partial molar(molal) Gibbs energy g(TP) (always), J/mole
    double *tpp_S;    ///< Partial molar(molal) entropy s(TP), J/mole/K
    double *tpp_Vm;   ///< Partial molar(molal) volume Vm(TP) (always), J/bar

    // additional arrays for internal calculation in ipm_main
    double *XU;      ///< dual-thermo calculation of DC amount X(j) from A matrix and u vector [L]
    double (*Uc)[2]; ///< Internal copy of IC chemical potentials u_i (mole/mole) at r-1 and r-2 [N][2]
    double *Uefd;    ///< Internal copy of IC chemical potentials u_i (mole/mole) - EFD function [N]
    char errorCode[100]; ///<  code of error in IPM      (Ec number of error)
    char errorBuf[1024]; ///< description of error in IPM
    double logCDvalues[5]; ///< Collection of lg Dikin crit. values for the new smoothing equation
    double *GamFs;   ///< Copy of activity coefficients Gamma before the first enter in PhaseSelection() [L] new

    double // Iterators for MTP interpolation (do not load/unload for IPM)
    Pai[4],    ///< Pressure P, bar: start, end, increment for MTP array in DataCH , Ptol
    Tai[4],    ///< Temperature T, C: start, end, increment for MTP array in DataCH , Ttol
    Fdev1[2],  ///< Function1 and target deviations for  minimization of thermodynamic potentials
    Fdev2[2];  ///< Function2 and target deviations for  minimization of thermodynamic potentials

    // Experimental: modified cutoff and insertion values (DK 28.04.2010)
    double
    // cutoffs (rescaled to system size)
    XwMinM, ///< Cutoff mole amount for elimination of water-solvent { 1e-13 }
    ScMinM, ///< Cutoff mole amount for elimination of solid sorbent { 1e-13 }
    DcMinM, ///< Cutoff mole amount for elimination of solution- or surface species { 1e-30 }
    PhMinM, ///< Cutoff mole amount for elimination of non-electrolyte condensed phase { 1e-23 }
    ///< insertion values (re-scaled to system size)
    DFYwM,  ///< Insertion mole amount for water-solvent { 1e-6 }
    DFYaqM, ///< Insertion mole amount for aqueous and surface species { 1e-6 }
    DFYidM, ///< Insertion mole amount for ideal solution components { 1e-6 }
    DFYrM,  ///< Insertion mole amount for major solution components (incl. sorbent) { 1e-6 }
    DFYhM,  ///< Insertion mole amount for minor solution components { 1e-6 }
    DFYcM,  ///< Insertion mole amount for single-component phase { 1e-6 }
    DFYsM,  ///< Insertion mole amount used in PhaseSelect() for a condensed phase component  { 1e-7 }
    SizeFactor; ///< factor for re-scaling the cutoffs/insertions to the system size
}
MULTI;

/// Indexation in a row of the pmp->SATX[][] array.
/// [0] - max site density in mkmol/(g sorbent);
/// [1] - species charge allocated to 0 plane;
/// [2] - surface species charge allocated to beta -or third plane;
/// [3] - Frumkin interaction parameter;
/// [4] species denticity or coordination number;
/// [5]  - reserved parameter (e.g. species charge on 3rd EIL plane)
enum IndexationSATX {
    XL_ST = 0, XL_EM = 1, XL_SI = 2, XL_SP = 3
};

/// One phase-assemblage stability violation, as ranked by
/// TMultiBase::WorstPhaseStabilityViolation(). See that method's own doc
/// comment for what a violation is and which thresholds define it.
struct PhStabViolation
{
    long int k = -1;            ///< phase index
    double   viol = 0.;         ///< size past native's threshold, in log10 units
    bool     wasAbsent = false; ///< true: absent but stable; false: present but unstable
    /// The phase's logSI is SATURATED at StabilityIndexes()' own overflow
    /// guard (ipm_chemical.cpp clamps ln_ax_dual to +609 / -608), i.e. it is
    /// a guard value and not a measured driving force. Rank ordering among
    /// clamped entries is therefore arbitrary - see PhStabCensus.
    bool     clamped = false;
};

/// What one phase-assemblage stability scan actually looked at. Exists
/// because WorstPhaseStabilityViolation() otherwise reports only an
/// EXTREMUM, and "the assemblage is consistent" and "every violation found
/// is one the caller may not act on" are the same observable from outside
/// (plan v5 section 91). The clamped count is the load-bearing one: on an
/// UNCONVERGED state most of the corpus's large systems saturate the
/// overflow guard, and a ranking built out of guard values is not a
/// measurement - it orders phases by nothing.
struct PhStabCensus
{
    long int phases = 0;          ///< pm.FI
    long int exempt = 0;          ///< kinetically restricted or twin-exempt, not classified
    long int scanned = 0;         ///< phases actually classified
    long int clamped = 0;         ///< of those, with a saturated logSI
    long int absentStable = 0;    ///< violations of the "should be present" kind
    long int presentUnstable = 0; ///< violations of the "should be absent" kind
};

extern const BASE_PARAM pa_p_;

// Data of MULTI
class TMultiBase
{
    char PAalp_; ///< Flag for using (+) or ignoring (-) specific surface areas of phases
    char PSigm_; ///< Flag for using (+) or ignoring (-) specific surface free energies
    std::shared_ptr<BASE_PARAM> pa_standalone;

    /// Work item 33. ScaleSystemToInternal()/RescaleSystemFromInternal() (ipm_simplex.cpp) guard
    /// the pm.DUL[]/pm.PUL[] scaling with "< 1e6" so the sentinel meaning "no upper limit" is not
    /// itself rescaled. That guard used to be RE-EVALUATED on whatever value was CURRENT at each
    /// call, which is provably wrong once ScFact > 1: a genuine bound below 1e6 that crosses the
    /// sentinel when multiplied, and an untouched value at or above 1e6, then occupy the SAME
    /// post-scale range, and no threshold test on the post-scale value alone can separate them
    /// (the ranges overlap on [1e6, 1e6*ScFact), and the ordinary sentinel 1e6 lies inside it).
    ///
    /// SO THE DECISION IS RECORDED, NOT RE-DERIVED - but recording only the DECISION is not
    /// enough either, because a bound can be REWRITTEN between the two calls: Set_DC_limits(true)
    /// (ipm_main.cpp, the warm pm.pNP path) writes pm.DUL[j] and pm.PUL[k] from INSIDE the scaled
    /// region, for exactly the species that carry metastability restrictions. Replaying a stale
    /// decision on those would be wrong in both directions. So both the pre-scale value and the
    /// value this call LEFT are kept: an entry still holding what we left is restored VERBATIM
    /// (exact - no multiply-then-divide rounding), and an entry that has moved since is a value
    /// written in internal units, which takes the original test.
    ///
    /// Filled fresh by every ScaleSystemToInternal(); consumed and cleared by the matching
    /// RescaleSystemFromInternal(). The two are 1:1 on this instance in every call path
    /// (CalculateEquilibriumState(), CalculateEquilibriumStateOptima()) - never nested, never
    /// re-entered - and an empty or mismatched vector falls back to the original test rather than
    /// reading out of bounds.
    std::vector<double> DUL_preScale_, DUL_postScale_;
    std::vector<double> DLL_preScale_, DLL_postScale_;
    std::vector<double> PUL_preScale_, PUL_postScale_;

    friend class TNode;
protected:
    /// Default logger for ipm chemical
    static std::shared_ptr<spdlog::logger> ipm_logger;

public:
    TNode *node1;

    /// This allocation is used only in standalone GEMS3K
    explicit TMultiBase( TNode* na_ = nullptr );
    virtual ~TMultiBase()
    {
        if(node1) {
           multi_kill();
        }
    }

    virtual void multi_realloc( char PAalp, char PSigm );
    void multi_kill();
    
    virtual BASE_PARAM* base_param() const
    {
       return pa_standalone.get();
    }

    MULTI* GetPM()
    { return &pm; }

    virtual void set_def( int i=0);
    virtual long int testMulti();


    //connection to mass transport
    void to_file( GemDataStream& ff );
    void to_text_file(const std::string& path, bool append=false);
    void solmod_to_text_file(const std::string& path);
    void solmod_to_json_file(const std::string& path);
    void from_file( GemDataStream& ff );
    template<typename TIO>
    void to_text_file_gemipm( TIO& out_format, bool addMui,
                              bool with_comments = true, bool brief_mode = false );
    template<typename TIO>
    void from_text_file_gemipm( TIO& in_format,  DATACH  *dCH );

    /// Writes Multi to a json/key-value string
    /// \param brief_mode - Do not write data items that contain only default values
    /// \param with_comments - Write files with comments for all data entries or as "pretty JSON"
    std::string gemipm_to_string( bool addMui, const std::string& test_set_name, bool with_comments = true, bool brief_mode = false );
    /// Reads Multi structure from a json/key-value string
    bool gemipm_from_string( const std::string& data,  DATACH  *dCH, const std::string& test_set_name );


    ///  Reads the contents of the work instance of the DATABR structure from a stream.
    ///   \param stream    string or file stream.
    ///   \param type_f    defines if the file is in binary format (1), in text format (0) or in json format (2).
    void  read_ipm_format_stream( std::iostream& stream, GEMS3KGenerator::IOModes type_f, DATACH  *dCH, const std::string& test_set_name );

    /// Writes the contents of the work instance of the DATABR structure into a stream.
    ///   \param stream    string or file stream.
    ///   \param type_f    defines if the file is in binary format (1), in text format (0) or in json format (2).
    ///   \param with_comments (text format only): defines the mode of output of comments written before each data tag and  content
    ///                 in the DBR file. If set to true (1), the comments will be written for all data entries (default).
    ///                 If   false (0), comments will not be written;
    ///                         (json format): interpret the flag with_comments=on as "pretty JSON" and
    ///                                   with_comments=off as "condensed JSON"
    ///  \param brief_mode     if true, tells that do not write data items,  that contain only default values in text format
    void  write_ipm_format_stream( std::iostream& stream, GEMS3KGenerator::IOModes type_f,
                                   bool addMui, bool with_comments, bool brief_mode, const std::string& test_set_name );
    virtual void copyMULTI( const TMultiBase& otherMulti );
    /// copyMULTI()'s body. realloc = true is copyMULTI() itself. realloc = false copies VALUES into this
    /// object's EXISTING arrays and allocates nothing: for restoring a snapshot into a live MULTI, whose
    /// TSolMod/TSorpMod/TKinMet objects hold pointers into those arrays (TSolMod::lnGamma = sd->arlnGam),
    /// so reallocating it would leave them dangling. Both objects must have identical dimensions - a
    /// snapshot of this same MULTI taken with copyMULTI(). Used by TNode::GEM_trace_regimes().
    void copyMULTIData( const TMultiBase& otherMulti, bool realloc );
    void read_multi(GemDataStream &ff, DATACH *dCH);
    /// Writing structure MULTI (GEM IPM work structure) to binary file
    void out_multi( GemDataStream& ff  );

    // New functions for TSolMod, TKinMet and TSorpMod parameter arrays
    void getLsModsum( long int& LsModSum, long int& LsIPxSum );
    void getLsMdcsum( long int& LsMdcSum,long int& LsMsnSum,long int& LsSitSum );
    /// Get dimensions from LsPhl array
    void getLsPhlsum( long int& PhLinSum,long int& lPhcSum );
    /// Get dimensions from LsMdc2 array
    void getLsMdc2sum( long int& DQFcSum,long int& rcpcSum );
    /// Get dimensions from LsISmo array
    void getLsISmosum( long int& IsoCtSum,long int& IsoScSum, long int& IsoPcSum,long int& xSMdSum );
    /// Get dimensions from LsESmo array
    void getLsESmosum( long int& EImcSum,long int& mCDcSum );
    /// Get dimensions from LsKin array
    void getLsKinsum( long int& xSKrCSum,long int& ocPRkC_feSArC_Sum,
                      long int& rpConCSum,long int& apConCSum, long int& AscpCSum );
    /// Get dimensions from LsUpot array
    void getLsUptsum(long int& UMpcSum, long int& xICuCSum);

    // EXTERNAL FUNCTIONS
    // MultiCalc
    void Alloc_internal();
    /// Total Gibbs energy G(X) of the converged system, in RT units (moles),
    /// computed identically for every solver path so that results are
    /// comparable across native AIA/SIA and Optima AOP/SOP/ROP.
    ///
    /// Why this exists (plan v5 section 4.1, CLAUDE.md 2026-08-25): total G is
    /// the correct correctness criterion for a Gibbs minimiser, and the
    /// OK/FAIL status is not - `f_/j_Solvus` AOP is a demonstrated case of a
    /// wrong answer reported OK (pH 5.1999 where native and a globalized step
    /// both give 5.1473). `pm.FX` cannot be used for this: it is written only
    /// by native's own IPM loop (ipm_main.cpp), never by
    /// CalculateEquilibriumStateOptima(), so it is stale or meaningless on the
    /// Optima paths.
    ///
    /// Reads the converged primal from pm.Y[] (both solver families leave the
    /// solution there - native via its IPM loop, Optima via the
    /// pm.Y[j]=pm.X[j]=state.x[j] unpack in ipm_optima.cpp) and, as a side
    /// effect, copies Y into X and refreshes pm.XF/pm.XFA. That is a no-op on a
    /// converged state, where X and Y already agree - do not call this mid-solve.
    ///
    /// It deliberately does NOT call GX(0.). GX() reads pm.G[], which is not a
    /// path-independent basis: native's GEM_IPM() resets pm.G[i] = pm.G0[i] on
    /// exit (ipm_main.cpp), dropping the fDQF and F0 activity-coefficient terms,
    /// while CalculateEquilibriumStateOptima() leaves them in. The excess term is
    /// therefore rebuilt here from G0[]+fDQF[]+F0[], which both paths do leave
    /// current and equal - see the implementation in ipm_chemical.cpp for the
    /// measurement that established this.
    double TotalGibbsEnergy();

    // ------------------------------------------------------------------------------------
    // REPORT-ONLY CERTIFICATE INSTRUMENTS (Phase 3 WP1; Docs/PLAN-defaults-and-fallbacks.md
    // s0.1 and its "Two cautions on that gate").
    //
    // Three numbers about the answer a mode returned, computed at exit and printed on the
    // CERT trace record beside the mass-balance fields. They are REPORT-ONLY: none of them
    // enters CERT's mb_pass, and freeze_diff.py's certificate gate does not read them.
    // Caution 1 of that block is why - a field that starts being emitted AND joins the pass
    // rule in one step makes every row the new clause rejects arrive as a `pass -> fail`
    // regression that the scorer cannot tell from a real one. Joining `pass` is a separate,
    // deliberate step that re-baselines inv-cert-standard-zero in the same commit.
    //
    // ALL THREE ARE COMPUTED ONLY WHEN native_trace_file() IS OPEN. With
    // GEMS3K_NATIVE_TRACE_FILE unset - every production call - not one line of this runs, so
    // a caller cannot pay for a diagnostic it never reads. The freeze, which is taken with
    // the trace on, is therefore also the gate on them: WP1's acceptance is a standard freeze
    // row-identical to 2026-09-17-STANDARD-promoted.txt in every scored column.
    //
    // WHY THEY REBUILD Gj RATHER THAN READING pm.G[] OR pm.F[]. Both are path-dependent at
    // exit, and reading either would be the DC_G0() shape - a number that looks like the
    // answer and is one path's private copy (CLAUDE.md s4). native's GEM_IPM() resets
    // pm.G[i] = pm.G0[i] on the way out (ipm_main.cpp, the FORCED_AIA tail), dropping the
    // fDQF and F0 excess terms, while CalculateEquilibriumStateOptima() leaves them in; and
    // pm.F[] is refreshed by the native IPM loop from pm.Y BEFORE the last descent step, so
    // on a native row it belongs to the previous iterate and on an Optima row to Optima's
    // own post-solve refresh. CertPrimalPotentials() therefore rebuilds F from
    // G0[]+fDQF[]+F0[] at the RETURNED pm.X[], exactly as TotalGibbsEnergy() rebuilds G, so
    // one definition serves all six modes.

    /// Primal chemical potentials at the RETURNED amounts pm.X[], into the caller's own
    /// array - a read-only mirror of PrimalChemicalPotentials() that writes no pm.* state
    /// and rebuilds Gj = G0+fDQF+F0 instead of reading pm.G[]. F[j] is left at 0 for a
    /// species below pm.DcMinM or in a phase PrimalChemicalPotentials() would skip, which is
    /// the same set that function leaves at 0. Sized to pm.L.
    void CertPrimalPotentials( std::vector<double>& F ) const;

    /// kkt_max: worst sign-aware reduced-gradient residual in RT over the species that are
    /// present at the answer, with the dual pm.U[] the last linear solve committed.
    /// s_j = F_j - sum_i U_i A(i,j); interior |s_j|, at a lower bound max(-s_j, 0), at an
    /// upper bound max(s_j, 0), and 0 for a species whose box is degenerate (DUL <= DLL - a
    /// kinetically fixed species is an EQUALITY constraint whose multiplier is unrestricted
    /// in sign; treating it as one-sided reported a spurious 2.48 on o_/t_Kaolinite's
    /// Quartz, ipm_optima.cpp's own KKT check). Mirrors that check's sign logic exactly with
    /// one deliberate difference: the log-barrier term -tau/X_j Optima adds for pure-phase
    /// species is NOT subtracted here, because it is Optima's internal objective and not the
    /// thermodynamic one - so on an Optima row this reads up to tau/X_j higher than Optima's
    /// own maxKKTResidual for a pure species near its floor.
    /// \param worstJ index of the species carrying the maximum, -1 if none.
    /// \return the maximum, or -1. if nothing could be scored.
    double CertKktMax( const std::vector<double>& F, long int& worstJ ) const;

    /// dual_free_dirs: how many directions the dual is free along at the answer. The INTERIOR
    /// species (present, and away from both box bounds) are the ones whose s_j = 0 fixes u;
    /// when their stoichiometry columns span rank r < N, the remaining N-r directions leave u
    /// undetermined and every bound-active species' reduced gradient - hence kkt_max above -
    /// depends on where along them the solver happened to stop. Measured on
    /// 07PSIna_G_simple_1 SHP: rank 5 of N = 6, the free direction is the redox one, and
    /// H2(aq) at the floor read s = +0.240 / -0.512 / -2.302 / -2.895 by warm start alone at
    /// identical G (plan v5 s123.6). So dual_free_dirs > 0 is the flag that says kkt_max is
    /// reading a lottery, and a threshold on kkt_max there is a lottery threshold.
    /// Same modified Gram-Schmidt and the same 1e-8 relative rank tolerance as the free-dual
    /// search in ipm_optima.cpp, so the two cannot disagree about the rank.
    /// \return N - rank, with rank and N returned in the out parameters.
    long int CertDualFreeDirs( const std::vector<double>& F, long int& rank, long int& nIC ) const;

    /// The RANK trace record (Phase 3 WP2): what the present species' stoichiometry, and the
    /// two normal-equation matrices each solver stage forms from it, actually look like at the
    /// RETURNED answer pm.X[] - as opposed to dual_free_dirs above, which reads only the
    /// INTERIOR species' columns. rank/sv_ratio/chg_res are properties of the geometry alone,
    /// independent of any solver weight; cond_ipm/cond_mbr fold in the weight each stage's own
    /// linear solve actually applies, so a system near-singular in the plain geometry can still
    /// be well posed once the weight is included (a species riding a tight box collapses its own
    /// column's effective magnitude) - or the reverse.
    struct CertRankReport
    {
        long int rank = 0;            ///< numerical rank of A_present (present-species columns), row-scaled
        long int of = 0;              ///< pm.N - the ambient dimension `rank` is measured against
        long int pres = 0;            ///< number of present species, pm.X[j] > pm.DcMinM
        double sv_ratio = -1.;        ///< sigma_min/sigma_max of the ROW-SCALED A_present; -1 if not computed
        double sv_ratio_raw = -1.;    ///< the same, without row scaling; -1 if not computed
        double chg_res = -1.;         ///< charge row's relative residual off the element rows' span; -1 if E<=0 or N-E<=0
        int chg_span = 0;             ///< 1 when chg_res < 1e-10, i.e. the charge row IS a combination of the element rows
        double cond_ipm = 1e300;      ///< cond(A_p diag(w) A_p^T), w = WeightMultipliers(false)'s shape at pm.X
        double cond_ipm_jac = 1e300;  ///< the same, after symmetric Jacobi scaling
        double cond_mbr = 1e300;      ///< cond(A_p diag(w) A_p^T), w = WeightMultipliers(true)'s shape at pm.X
        double cond_mbr_jac = 1e300;  ///< the same, after symmetric Jacobi scaling
    };

    /// Fills a CertRankReport at the RETURNED pm.X[] - see CertRank() in ipm_main.cpp for the
    /// construction of every field and why each is computed the way it is. REPORT-ONLY: reads
    /// pm.A/pm.X/pm.DLL/pm.DUL/pm.RLC and writes nothing, including pm.W[] (WeightMultipliers()
    /// itself is NOT called - the weight this computes is a local copy of its arithmetic, not a
    /// second call to it, because pm.W[] is live solver scratch a report-only path must not
    /// touch). Leaves `r` at its default (all-1e300/-1) if the guard at the top of the
    /// implementation fails.
    void CertRank( CertRankReport& r ) const;

    /// curv_min: the smallest eigenvalue, over every PRESENT multicomponent non-aqueous
    /// solution phase, of that phase's symmetrised finite-difference curvature block at the
    /// answer. A negative value means the phase converged INSIDE ITS OWN SPINODAL - a point
    /// that satisfies stationarity and mass balance and is a maximum along the unmixing
    /// direction, which no first-order test can see (Michelsen 1982 II). Lead:
    /// j_CASHNK's drifting limb, where native at the shipped pa_DK = 1e-6 stops 1.8e-4 above
    /// the minimum with end-member fractions up to 0.38 away from converged (plan v5 s96.1).
    ///
    /// The block, the present-end-member test and the step h are taken verbatim from the
    /// pa_PhaseHessianFloor site in ipm_optima.cpp so the two measure the same object; the
    /// difference is that this one does NOT floor - SymEigFloorInPlace() returns the repaired
    /// matrix, not its spectrum, so the smallest eigenvalue is computed here by a sweep-only
    /// Jacobi that accumulates no eigenvectors. The aqueous phase is excluded for the same
    /// reason the floor excludes it: FD columns for its many near-floor trace species are
    /// noise.
    ///
    /// MUTATES AND RESTORES. Each column perturbs pm.X[i] by h and re-runs
    /// TotalPhasesAmounts() + CalculateActivityCoefficients(LINK_UX_MODE); pm.X[i] is put
    /// back immediately and XF/XFA/lnGam/Gamma/fDQF/F0 are refreshed at the original X on the
    /// way out, exactly as the Optima site does. pm.G[] is saved and restored verbatim on top
    /// of that, because the refresh would otherwise leave a native row at G0+fDQF+F0 where
    /// GEM_IPM() had reset it to G0 - a difference a later warm call would see. Called from
    /// native_trace_run_result(), i.e. after packDataBr() has already extracted the answer,
    /// so an imperfect restore cannot change what THIS call returns; it could change what a
    /// later warm call starts from, which is what the row-identity freeze gate tests.
    /// \param worstPhase index of the phase carrying the minimum, -1 if none.
    /// \return the minimum, or +1e300 if no phase qualified.
    double CertCurvMin( long int& worstPhase );

    /// pa_StabTPD = 1 (FABLE Phase 3 WP6, plan v5 s139.2): the tangent-plane stability scan of every
    /// ABSENT or TRACE non-aqueous multicomponent phase at the returned answer. Absent = XF <= DSM, or
    /// XF < 1e-6 of the summed phase amounts (trace), or every end-member at or below 1e3 x the certificate's species floor (so Optima's floor-held phases,
    /// kept above DSM by pa_OptimaZeroAbsent = 2, count as absent - the first AOP control built on
    /// XF > DSM alone read 2.5e+02 on floor-held phases). Sorption/ion-exchange/polyelectrolyte phases
    /// are skipped (their activity is not a mole-fraction model). Search: Michelsen successive
    /// substitution from every vertex, the ideal closed form and the current composition; only
    /// every evaluated composition is scored by its TPD directly (a lnTM value from a non-converged start
    /// is not a TPD and read +25.7 RT on a present multi-site phase in the first probe).
    /// \param worstPhase phase carrying the minimum, -1 if none.
    /// \param nScanned   absent phases scanned.
    /// \param nDisagree  phases the single-point index calls stable while the search finds them unstable.
    /// \return min over scanned phases of the smallest TPD(y) evaluated, in RT (< 0: unstable), +1e300 if none.
    /// MUTATES AND RESTORES BY COPY, exactly as CertCurvMin() does, for the same measured reason.
    double CertStabTPD( long int& worstPhase, long int& nScanned, long int& nDisagree );
    /// PROTOTYPE (plan v5 §140.12): composition search for one absent non-ideal phase; see ipm_main.cpp.
    double NativeTpdPhase( long int k, long int p0, std::vector<double>& ybest );

    double CalculateEquilibriumState( /*long int typeMin,*/ long int& NumIterFIA, long int& NumIterIPM );
    void InitalizeGEM_IPM_Data();
    virtual void DC_LoadThermodynamicData( TNode* aNa = nullptr );

#ifdef USE_OPTIMA_SOLVER
    // Equilibrium via the Optima library's general primal-dual interior-
    // point NLP solver, as an alternative to the IPM/MBR loop above - see
    // ipm_optima.cpp. Dispatched from TNode::GEM_run() (node.cpp) for
    // NEED_GEM_AOP/SOP, exactly as CalculateEquilibriumState() is
    // dispatched for NEED_GEM_AIA/SIA. Reads pm.pNP (set by GEM_run()
    // before the call, same flag AIA/SIA already use) to choose a cold
    // (AutoInitialApproximation(), pNP==0, "AOP") or warm (reuse the
    // existing pm.Y[], pNP==1, "SOP") starting point. Any control
    // conditions registered via SetControlCondition_pH()/_Eh() (or a
    // custom EqControlCondition appended directly to
    // optima_control_conditions) are folded into the same joint Newton

    /// Request dn/db sensitivity derivatives from the Optima solver.
    /// Off by default: computing them costs an extra solve of the already-
    /// factorised KKT system per bulk-composition column, and nothing in the
    /// ordinary AOP/SOP/ROP paths needs them.
    ///
    /// This is the plumbing step of the ODML predictor described in
    /// Docs/literature/FINDINGS-Leal2020-ODML-prediction.md - NOT a predictor.
    /// Optima carries the machinery already: Problem::c are sensitivity
    /// parameters and Problem::bec is d(be)/dc, so mapping c to the bulk
    /// composition means bec = I and Sensitivity::xc is exactly dn/db. No
    /// restructuring of the objective or the constraints is involved.
    bool optima_want_sensitivity = false;

    /// dn/db from the last solve when optima_want_sensitivity was set, row-major
    /// [ (L+R) x N ]. Empty if it was not requested or the solve did not reach
    /// the sensitivity step.
    std::vector<double> optima_dndb;
    long int optima_dndb_rows = 0, optima_dndb_cols = 0;

#endif   // USE_OPTIMA_SOLVER
    // ---- NATIVE-PATH MEMBER, deliberately outside the Optima block ----------
    // It was inside it until 2026-09-08, pasted in the middle of
    // CalculateEquilibriumStateOptima()'s own doc comment (which resumes right
    // after this declaration - "// solve; with none registered ..."), which is
    // how it went unnoticed: nothing here reads as belonging to Optima.
    //
    // But PhaseSelectionSpeciationCleanup() uses it on the NATIVE path -
    // ipm_chemical.cpp:1605,1610, assigned in ipm_main.cpp:803 - so with
    // USE_OPTIMA_SOLVER=OFF, which is this project's DEFAULT
    // (GEMS3K/CMakeLists.txt), the library did not compile at all.
    //
    // The block is CLOSED and REOPENED around the declaration rather than the
    // declaration being moved up, so the member keeps its exact position in the
    // class and an Optima build's layout is byte-for-byte unchanged. Moving it
    // would be an ABI change, and gems-benchmark/CLAUDE.md s5 records what those
    // cost: any binary predating one is stale and will not say so.
    /// One flag per phase: has PhaseSelect() already re-inserted this phase at its
    /// COMPOSITION CEILING rather than at the fixed pa_DFYs during this solve?
    /// Sized pm.FI and cleared in GEM_IPM(), so it spans the phase-selection
    /// passes of one solve and nothing more.
    ///
    /// WHY A ONE-SHOT RULE. Inserting min(pa_DFYs, ceiling) instead of skipping is
    /// what stops a trace phase reading ABSENT when it belongs in the assemblage
    /// (plan v5 section 93.3: T14_ball120 loses five, Chromite at logSI +6.39).
    /// But the ceiling consumes the whole element budget, so an unbounded clamp
    /// re-creates the thrash the 2026-09-05c skip was written to stop - measured:
    /// clamping alone left both rescued projects at status 3 (section 77.1), and a
    /// safety FRACTION was swept and came out chaotic and non-monotone, the two
    /// projects disagreeing about which value works (section 77.3). So the phase
    /// gets exactly one budget-sized attempt; if it comes back lost, it is skipped
    /// from then on. That bounds the wasted work at one pass per phase without a
    /// tuned number anywhere.
    std::vector<char> insBudgetTried;
#ifdef USE_OPTIMA_SOLVER
    // solve; with none registered this is a plain equilibrium solve,
    // architecturally equivalent to AIA/SIA but solved via Optima instead
    // of GEMS3K's own IPM/MBR.
    //
    // `reaktoroMode` (default false = AOP/SOP's own established behavior,
    // unchanged): when true, dispatched for NEED_GEM_ROP instead - a
    // faithful port of Reaktoro's OWN equilibrium mechanism onto this same
    // objective/constraint plumbing, not just AOP's seed/options swapped
    // in. Differs from the AOP/SOP path in every one of: (1) initial
    // guess - a uniform tiny seed for every species (Reaktoro's own
    // ChemicalState default), no AutoInitialApproximation() call at all;
    // (2) Hessian - PartiallyExact (Reaktoro's own default:
    // EquilibriumHessian.cpp's approximate() as a base, with columns for
    // Optima-reported *basic* variables, opts.ibasicvars, overwritten by a
    // finite-difference port of Reaktoro's autodiff-exact
    // d(chem.potential)/dn - GEMS3K has no autodiff, so this is an FD
    // port of the same mechanism, not a claim of bit-identical numerics);
    // (3) Optima::Options - left at the library's own untouched defaults
    // (Reaktoro's own EquilibriumOptions.hpp carries a plain Optima::
    // Options with no override), not GEMS3K's pa_p->IIM/OptimaTol/
    // OptimaMaxStepRatio overrides; (4) fallback - Reaktoro's own single
    // retry (EquilibriumSolver.cpp), re-solving ONCE from the ORIGINAL
    // (pre-first-solve) state with backtracksearch.apply_min_max_fix_and_
    // accept toggled - NOT AOP's solvent-dominance reseed retry. See
    // GEMS3K's CLAUDE.md, 2026-08-24, "ROP" for the full scoping. Always a
    // single mode (pm.pNP forced to the cold/AIA-equivalent convention by
    // TNode::GEM_run() for NEED_GEM_ROP) - Reaktoro's own default
    // equilibrate() always starts from the same seed regardless of any
    // prior state, so there is no warm-start ROP variant.
    /// \param runKinetics run the kinetics/metastability time step (RunKineticsStep()). TRUE for a
    ///        normal AOP/SOP call; FALSE for CalculateEquilibriumStateHOP()'s Optima leg, whose native
    ///        leg has already advanced it for this time step - running it twice would double the rate.
    double CalculateEquilibriumStateOptima( long int& NumIterFIA, long int& NumIterIPM, bool reaktoroMode = false,
                                            bool runKinetics = true );

    /// HYBRID: native selects the species (its own IPM/MBR/PSSC pipeline,
    /// cold), then Optima finishes, warm-started from native's converged
    /// primal AND dual - see NODECODECH's own comment in databr.h for why
    /// this is a separate caller-selected mode (dispatched for NEED_GEM_HOP)
    /// rather than something AOP does internally.
    ///
    /// Guarantees the result is never worse than a plain native solve: if
    /// native converges but the Optima leg then fails (throws, for any
    /// reason including pa_DW's hard-error gate), native's own
    /// already-converged state is restored and reported as a soft
    /// BAD_GEM_HOP instead of being lost - before this, the Optima leg's
    /// thrown TError propagated straight to TNode::GEM_run()'s OUTER catch,
    /// which never calls packDataBr(), so a converged native answer (e.g.
    /// 693 iterations on 07PSIna_G_vcomplex @ 80 C) was silently discarded
    /// whenever the warm Optima leg on top of it could not finish. See
    /// GEMS3K/CLAUDE.md and Docs/gems3k-optima-plan-v5.md, section 64.6.
    ///
    /// `warmNative` selects the SHP variant (NEED_GEM_SHP): the native leg
    /// starts warm (pm.pNP = 1) from whatever this node already holds,
    /// instead of cold. It is for sequential work - a sweep or a transport
    /// loop, where the previous point's converged state is sitting right
    /// there and HOP as built throws it away at every call. It is exactly
    /// HOP with native's SIA in place of native's AIA. Its headline result
    /// is not the sweep saving it was built for: on the projects whose
    /// native SIA cannot re-solve their own converged state (plan-v5 29.2
    /// and 60.5) SHP converges anyway, at ~4x fewer iterations than HOP,
    /// because the state it hands native is the OPTIMA leg's rather than
    /// native's own. On a sweep the saving is project-dependent - 2.4x less
    /// wall time than HOP on j_10TH_G_seawater, ~nothing on
    /// j_Solvus_G_series1, a 2.3x LOSS on j_Kaolinite_G_pHtitr, whose
    /// native warm start is itself more expensive than its cold one. Full
    /// tables in NODECODECH's comment (databr.h). Carries a COLD FALLBACK, so a
    /// project whose native SIA refuses its own converged state (ten of
    /// them, plan-v5 section 60.5) degrades to exactly HOP at the cost of
    /// one cheap wasted attempt. See NODECODECH's comment in databr.h.
    double CalculateEquilibriumStateHOP( long int& NumIterFIA, long int& NumIterIPM,
                                         bool warmNative = false );

    /// Registers (or replaces, by name) a pH control condition for the
    /// next CalculateEquilibriumStateOptima() call. Persists across calls
    /// until ClearControlConditions() - not single-shot. `tolerance < 0`
    /// (the default) defers to pa_p->GAS via EqControlCondition::
    /// defaultToleranceFn() instead of a hardcoded constant - see
    /// ipm_optima.h.
    void SetControlCondition_pH( double pH_target, double tolerance = -1. );
    /// Registers (or replaces, by name) an Eh control condition (V) for
    /// the next CalculateEquilibriumStateOptima() call. `tolerance < 0`
    /// (the default) defers to pa_p->GAS, same as SetControlCondition_pH().
    void SetControlCondition_Eh( double Eh_target, double tolerance = -1. );
    /// Removes every registered control condition.
    void ClearControlConditions();

    /// Read-only access to the registered control conditions - after a
    /// CalculateEquilibriumStateOptima() call, each entry's
    /// titrantAmount/achievedValue/targetMet fields hold that call's result.
    const std::vector<EqControlCondition>& GetControlConditions() const
    { return optima_control_conditions; }

    /// Worst-case (infinity-norm) mass-balance residual |A*Y - B| over the
    /// current pm.Y[]/pm.B[] - a standalone diagnostic for callers/tests to
    /// inspect directly. CalculateEquilibriumStateOptima()'s own post-solve
    /// trustworthiness check uses CheckMassBalanceResiduals() instead (the
    /// same per-IC-tolerance check the native solver's own testMulti() path
    /// uses), not this method.
    double OptimaMaxMassBalanceResidual();

    /// Shared solvent-collapse detection/reseed-value computation, used by
    /// both AOP's and ROP's own retry logic inside
    /// CalculateEquilibriumStateOptima() (ipm_optima.cpp) - unified
    /// 2026-08-24 after the two independently-written versions' reseed
    /// formulas had drifted apart (AOP's own included an `otherTotal`
    /// fallback and an `upperBound` cap that an earlier ROP draft lacked).
    /// Detects whether the aqueous solvent (pm.LO) dominates its own phase
    /// by mass in the given trial amounts `x` (length pm.L) - the
    /// confirmed signature of the "aqueous solvent collapses to its
    /// floor" trap (GEMS3K's CLAUDE.md, 2026-08-23/24). Returns false (no
    /// correction needed) if it already dominates or `pm.LO` doesn't
    /// exist for this project; otherwise returns true and sets
    /// `waterSeedOut` to `min(bulk-H/2, bulk-O)` (falling back to the
    /// phase's own other-species total if H/O can't be resolved by name),
    /// capped by `upperBound` if positive.
    ///
    /// Superseded (not removed - see its own declaration comment) by
    /// `DetectPhaseCollapseAndReseed()` below, which generalizes this
    /// aqueous-only check to every multicomponent phase - the underlying
    /// "phase reported absent" cliff in `PrimalChemicalPotentials()`
    /// (`YF[k]<=pm.DSM`) is generic, aqueous just being the first
    /// instance this investigation happened to hit (an extra,
    /// aqueous-specific `Y[pm.LO]<=pm.XwMinM` check sits alongside the
    /// generic one, so aqueous remains a strictly stricter case of the
    /// same trap, not a separate mechanism).
    /// True when the system genuinely contains an aqueous phase.
    /// Tests the phase classifier, NOT pm.LO - see the definition for why.
    bool HasAqueousPhase() const;

    bool DetectSolventCollapseAndReseed( const double* x, double upperBound, double& waterSeedOut );

    /// Generalizes DetectSolventCollapseAndReseed() (see its own
    /// declaration comment for why) to every multicomponent phase
    /// (`0..pm.FIs-1` with more than one species - single-species phases
    /// have their own separate `pm.PhMinM` guard and are not at risk of
    /// the "sparse allocation starves this phase" failure mode this
    /// targets, since there's only ever one species to allocate to).
    /// Scans the given trial amounts `x` (length pm.L) for any phase
    /// whose total sits too close to its own "phase absent" threshold
    /// (`pm.DSM`, or `max(pm.DSM,pm.XwMinM)` for the aqueous phase
    /// specifically) and, for each such phase, computes a physically-
    /// grounded reseed bound - `min` over every IC the phase's own
    /// end-members actually touch of `bIC[i] / (that phase's own largest
    /// stoichiometric coefficient for IC i)` - the same generalization of
    /// the water-specific `min(bulk-H/2, bulk-O)` formula applied to any
    /// phase's own stoichiometry - distributed EVENLY across the phase's
    /// own end-members (a simple, safe default: it doesn't guess which
    /// end-member the true equilibrium favors, it only needs to clear the
    /// phase-absence cliff - Newton's own gradient-driven iteration is
    /// expected to correct the actual mix from there). Appends
    /// `(species index, reseed value)` pairs to `reseedsOut` (only for
    /// entries where the reseed value exceeds the trial amount already
    /// there) and returns true iff at least one was appended.
    /// `excludePhaseIdx` (default -1, none): skip this one phase entirely
    /// - used to avoid double-handling the aqueous phase, which already
    /// has its own dedicated, separately-validated
    /// DetectSolventCollapseAndReseed() call (kept as-is, not replaced,
    /// to avoid any behavior change to the already-validated water-
    /// specific path this generalizes beyond).
    bool DetectPhaseCollapseAndReseed( const double* x,
                                        std::vector<std::pair<long int,double>>& reseedsOut,
                                        long int excludePhaseIdx = -1 );

    /// A genuine LP-feasibility seed for ROP: solves min sum_j(n_j) s.t.
    /// pm.A*n = pm.B, n >= 0 via a small self-contained two-phase dense
    /// simplex (ipm_optima.cpp, anonymous-namespace TwoPhaseSimplexMinSum())
    /// - the same class of computation Reaktoro's own compared-against
    /// iteration counts actually started from (reaktoro_bench.py's
    /// build_seed(), scipy.optimize.linprog), not a bare uniform-tiny
    /// value. Exists to replace ROP's reliance on a uniform seed plus
    /// retries entirely: two independent attempts at reordering those
    /// retries (GEMS3K's CLAUDE.md, 2026-08-24, "Why ROP doesn't just try
    /// the fast (water-first) option first") both broke a real project
    /// (j_CASHNK), because this is a non-convex NLP and different retry
    /// orderings are different Newton starting points that can converge to
    /// different local optima - a single deterministic seed removes that
    /// whole risk class rather than trying a third ordering.
    ///
    /// Returns false (leaving `nOut` unspecified) if the LP itself reports
    /// infeasibility, hits its iteration cap, or - checked internally as a
    /// self-verification before ever trusting the simplex's own output -
    /// the returned point doesn't actually satisfy A*n=b (within tolerance)
    /// and n>=0. The caller is expected to fall back to a simpler seed in
    /// that case, not trust a partial/unverified result.
    ///
    /// Deliberately NOT redistribution-aware on its own: a plain min-sum LP
    /// vertex concentrates mass in as few species as possible (the same
    /// mechanism that made AutoInitialApproximation()'s own LP-simplex seed
    /// starve the aqueous phase in the original solvent-collapse trap,
    /// GEMS3K's CLAUDE.md 2026-08-23/24) - the caller is expected to run
    /// the LP output through DetectSolventCollapseAndReseed()/
    /// DetectPhaseCollapseAndReseed() afterward, same as this method's own
    /// call site in CalculateEquilibriumStateOptima() does.
    bool LPFeasibilitySeed( std::vector<double>& nOut );
    /// PROTOTYPE (plan v5 §140.15): THERMOCHIMICA-style column-generation cold seed; see ipm_optima.cpp.
    bool ColumnGenerationSeed( std::vector<double>& nOut, double tol );
    bool PotentialSpaceFinish( double dcFloor, const std::vector<double>& xlower, const std::vector<double>& xupper );

    /// Dual of the LINEARISED-Gibbs LP (min sum_j G0[j]*n_j s.t. A n = b,
    /// n >= 0), computed with the same simplex as LPFeasibilitySeed() and
    /// therefore just as independent of any solver state. Its value is that
    /// pricing a species against it, s_j = G0[j] - sum_i y_i A[i,j], is a
    /// meaningful measure of how far that species is from being stable -
    /// which the feasibility LP's own dual is not, being an artefact of
    /// "minimise total moles". Used by OptimaReducedPreSolve() to choose a
    /// generous initial active set. Returns false (leaving yOut untouched)
    /// if the LP fails or its dual does not verify against LP optimality.
    bool LPGibbsDual( std::vector<double>& yOut, const double* cost = nullptr );

    /// Phase-assemblage stability scan for the Optima path - refreshes
    /// pm.YF/pm.YFA from pm.Y, calls StabilityIndexes(), and reports the
    /// single worst disagreement between "is this phase in the assemblage"
    /// and "does its stability index say it should be", using native's own
    /// PhaseSelect() thresholds (pa_p->DF / pa_p->DFM, ipm_chemical.cpp).
    ///
    /// Factored out of CalculateEquilibriumStateOptima()'s final
    /// trustworthiness check so that the phase-selection retry tier and that
    /// check share ONE definition of "violation". They must agree exactly:
    /// a retry chasing a violation the final check would not report loops
    /// pointlessly, and a retry blind to one the check does report cannot
    /// fix it. Exemptions (kinetically restricted phases; phases the
    /// caller has already deactivated) live here for the same reason.
    ///
    /// \param presenceThreshold  phase total at/below which a phase counts
    ///        as absent - NOT bare pm.DSM, see the call site.
    /// \param dcFloor  the species lower-bound floor the Optima path built its
    ///        boxes from. A phase every one of whose species sits AT its own
    ///        lower bound was pinned out by the solver and is absent however
    ///        large its total happens to be; without this test the magnitude
    ///        rule scales with the phase's species COUNT while the threshold
    ///        does not, so a wide phase held entirely at the floor is reported
    ///        present. Measured on 07PSIna_G_complex_1: pa_DS = 1e-20 collapses
    ///        presenceThreshold onto dcFloor*10 = 1e-12, and gas_gen's ten
    ///        species at exactly dcFloor = 1e-13 sum to exactly 1e-12 - the
    ///        single reason that project reports BAD rather than OK.
    /// \param exemptSpecies  optional, size >= pm.L when non-null: species
    ///        the caller has deliberately fixed and whose phase must be
    ///        skipped entirely (the interchangeable-twin case).
    /// \param violOut  size of the worst violation, in logSI units past
    ///        the threshold; 0 when none.
    /// \param wasAbsentOut  true if the worst violation is "absent but
    ///        stable", false if "present but unstable".
    /// \param rankedOut  optional: EVERY violation found, ordered by size
    ///        (descending, ties broken by phase index so the order is
    ///        reproducible). A caller that may not act on the worst one can
    ///        walk down this list instead of stopping, which is the whole
    ///        reason it exists - the aqueous phase and any phase in the wrong
    ///        phSelState are unactionable, and before this the loop could not
    ///        tell that apart from "no violation at all".
    ///        NOTE the ordering is only meaningful where the entries are not
    ///        `clamped`: see PhStabCensus.
    /// \param censusOut  optional: what the scan looked at, including how many
    ///        of the scanned phases had a logSI saturated at StabilityIndexes()'
    ///        overflow guard. ANY consumer of this ranking on an unconverged
    ///        state must read that count - measured on 07PSIna_G_complex_1
    ///        @ 80 C, 289 of 314 phases were pegged at the guard, so the
    ///        "worst" violation was picked out of a set that carries no
    ///        ordering information at all (plan v5 section 91).
    /// \return index of the worst violating phase, or -1 if the assemblage
    ///         is self-consistent.
    long int WorstPhaseStabilityViolation( double presenceThreshold, double dcFloor,
                                           const char* exemptSpecies,
                                           double& violOut, bool& wasAbsentOut,
                                           std::vector<PhStabViolation>* rankedOut = nullptr,
                                           PhStabCensus* censusOut = nullptr );

    /// Species-level dimension-reduction pre-solve - see
    /// BASE_PARAM::OptimaDimReduce (above) for what it does and why.
    /// Solves the equilibrium over a reduced set of species, growing that set
    /// by pricing the omitted columns on the reduced dual, and leaves its
    /// answer in pm.Y[] (primal) and pm.U[] (dual) for the ordinary
    /// full-dimension solve to warm-start from and verify.
    /// \param maxPasses  readmission passes allowed (pa_OptimaDimReduce).
    /// \param dcFloor    the same species floor the full path uses.
    /// \param dimTol     the initial-set rule, with pa_OptimaDimReduceTol's own
    ///        sign convention (> 0 an absolute RT threshold, < 0 the rank rule
    ///        -|dimTol| x N, 0 the seed's own support only). Passed rather than
    ///        read from BASE_PARAM so the caller can retry a discarded
    ///        pre-solve under the OTHER rule - see the call site.
    /// \param passBudget  Optima iteration budget of each pass
    ///        (optima_presolve_pass_budget()).
    /// \param iterationsOut  Optima iterations spent here, for pm.ITG.
    /// \param activeOut  size of the final active set, for logging.
    /// \return true if pm.Y[]/pm.U[] now carry a state worth warm-starting
    ///         from; false if the pre-solve was skipped or discarded, in which
    ///         case neither array was modified in a way the caller must undo.
    bool OptimaReducedPreSolve( long int maxPasses, double dcFloor, double dimTol,
                                long int passBudget,
                                long int& iterationsOut, long int& activeOut );
#endif

    // acces for node class
    TSolMod * pTSolMod (int xPH);

    long int CheckMassBalanceResiduals(double *Y );
    double ConvertGj_toUniformStandardState( double g0, long int j, long int k );
    double PhaseSpecificGamma( long int j, long int jb, long int je, long int k, long int DirFlag = 0L ); // Added 26.06.08

    double HelmholtzEnergy( double x );
    double InternalEnergy( double TC, double P );

protected:

#ifdef USE_OPTIMA_SOLVER
    /// Control conditions registered via SetControlCondition_pH()/_Eh()
    /// (or appended directly by a caller building a fully custom
    /// EqControlCondition - see ipm_optima.h), active for the next
    /// CalculateEquilibriumStateOptima() call. `stoich`/`fixedGradientFn`/
    /// `achievedValueFn` are already fully resolved by the time an entry
    /// lands here (index lookups need only CSD, valid immediately after
    /// GEM_init(); the closed-form functions capture resolved indices by
    /// value and read G0[]/T lazily, at actual call time).
    std::vector<EqControlCondition> optima_control_conditions;
#endif

    /// When true, SmoothingFactor() returns exactly 1.0 regardless of
    /// pm.FitVar[3]/[4], disabling the IPM-2 chemical-potential smoothing
    /// (ipm_chemical.cpp, DC_PrimalChemicalPotentialUpdate()'s
    /// "F0 = Fold + dF0 * SmoothingFactor()") for the duration of one
    /// CalculateEquilibriumStateOptima() call.
    ///
    /// Why this exists (see CLAUDE.md, 2026-08-25, plan-v5 Phase A / A.1):
    /// that smoothing blends the current chemical potential with `Fold =
    /// pm.F0[j]` - the value left by the PREVIOUS call - so pm.F0[] is a
    /// persistent accumulator across objective evaluations. Optima's
    /// objective callback invokes CalculateActivityCoefficients(LINK_UX_MODE)
    /// every Newton iteration, whose first statement is SetSmoothingFactor(),
    /// and whose result reaches Optima's gradient via
    /// pm.G[j] = G0[j] + fDQF[j] + pm.F0[j]. With a smoothing factor s < 1
    /// the minimised objective is therefore an exponentially-weighted moving
    /// average over the HISTORY of x evaluated, not a function of x alone -
    /// and no Newton method converges against a moving objective. With s == 1
    /// the blend is an exact algebraic no-op (F0 = Fold + (F0-Fold)*1 = F0).
    ///
    /// Set/cleared ONLY by CalculateEquilibriumStateOptima(). Native
    /// AIA/SIA never touches it, so native behaviour is bit-identical by
    /// construction - deliberately not #ifdef USE_OPTIMA_SOLVER-gated, so
    /// that SmoothingFactor() itself stays free of conditional compilation;
    /// in a non-Optima build this is simply always false.
    ///
    /// Measured scope, before assuming this changes much: across the whole
    /// Resources/gems3k suite only j_CASHNK (s = 0.875..0.123), j_GEOTHERM
    /// (s ~ 0.99998) and tools/Cu-Pourbaix (s = 0.989..0.404) have s != 1 at
    /// all; the other nine projects already evaluate s == 1.0 exactly for
    /// any pm.IT, so for them this flag is provably inert and their results
    /// must stay byte-identical.
    bool optima_disable_smoothing = false;

    /// True only while CalculateEquilibriumStateHOP()'s Optima leg is running
    /// on top of a SUCCESSFUL native solve. Read at exactly one place: the
    /// dimension-reduction gate in CalculateEquilibriumStateOptima(), which is
    /// otherwise cold-start-only (pm.pNP == 0).
    ///
    /// Why the cold-start-only rule does not apply here, which is the whole
    /// content of this flag. That rule reasons that a warm start already
    /// carries the consistent (primal, dual) pair the pre-solve exists to
    /// manufacture, so a reduction in front of it is waste when it discards
    /// and destructive when it settles. True of SOP - an ordinary warm restart
    /// at a nearby composition, where the incoming pair IS an Optima fixed
    /// point. NOT true of HOP: the incoming pair is NATIVE's, and native
    /// decides absence by truncating below pa_DcMin (1e-33) and dropping the
    /// species from its own linear system entirely, while Optima's box floor
    /// is pa_DHB (1e-13, twenty orders higher). So every species native calls
    /// absent arrives sitting AT Optima's floor, carrying a reduced gradient
    /// that need not satisfy Optima's own complementarity test - and on a
    /// large system those absent species are what the max-norm KKT residual
    /// is made of (section 57: one absent gas phase holding 99 % of it). That
    /// is the difference between a warm start Optima can verify in O(1)
    /// iterations and one it cannot verify at all.
    ///
    /// So on the HOP path the pre-solve is not manufacturing a pair - it is
    /// re-expressing native's own assemblage in Optima's terms, at a
    /// dimension where the residual is not dominated by species that are not
    /// there. The initial active set needs no new selection rule for that:
    /// OptimaReducedPreSolve() builds it from whatever pm.Y[] holds on entry,
    /// which on this path is native's converged answer. Section 62.6 measured
    /// that a native-derived assemblage errs only by INCLUDING too much,
    /// which the pre-solve's own pricing repairs at the cost of iterations
    /// rather than correctness. See Docs/gems3k-optima-plan-v5.md section 64.7
    /// / handoff item 8b.
    ///
    /// EXPLICIT OPT-IN ONLY: this flag lets the reduction run on the HOP leg
    /// when pa_OptimaDimReduce > 0, and AUTO (0) does NOT reach it. Measured on
    /// 07PSIna_G_mid_1, which is above the AUTO size gate and does not need the
    /// help: warm verification alone is 2 Optima iterations / 6 ms, and putting
    /// a reduction in front of it costs 35 / 85 ms for the same answer. The
    /// reduction pays here only where the warm verification does not work, and
    /// nothing static predicts which projects those are (section 39.4) - so it
    /// is per-project, like every other knob on this path whose response is not
    /// uniformly favourable. Default HOP behaviour is unchanged by this flag.
    bool optima_hop_leg = false;

    /// Per-leg cost of the last two-leg (HOP/SHP) solve - work item 38.
    ///
    /// CalculateEquilibriumStateHOP() reports NumIterFIA = fiaN + fiaO and
    /// NumIterIPM = ipmN + ipmO, one number each, and TNode::GEM_Iterations()
    /// passes those on. That total is the TRUE cost of the call and stays as it
    /// is - but it makes a per-iteration cost derived from it a BLEND of two
    /// solvers, which is exactly the number the warm standard could not use:
    /// measured over a transport loop, an SHP iteration came out at 223.1 us on
    /// CalcColumn and 511.3 us on LimSeawat1 against 17.5 and 14.0 for a native
    /// one (plan v5 section 130.2a), and those two figures cannot be attributed
    /// to a path while the legs are summed. Native's own per-iteration cost is
    /// nearly constant across the two chemistries (1.25x spread) and the Optima
    /// family's is not (2.3-2.9x), so the blend cannot be calibrated away either.
    ///
    /// Recorded, never acted on: nothing in the solver reads these fields, so
    /// they add no wall-clock-dependent decision - which matters, because every
    /// parallel benchmark freeze depends on the solver having exactly one of
    /// those (pa_OptimaMaxSeconds, off).
    ///
    /// `valid` is false unless the LAST GEM_run() dispatched a two-leg mode;
    /// TNode::GEM_run() clears it before every dispatch so a single-leg call
    /// cannot report a previous HOP call's split. It is set from a DESTRUCTOR, so a
    /// call that THROWS is recorded too, with `failed` true - see the HopLegRecorder
    /// comment in CalculateEquilibriumStateHOP() for what such a record can and
    /// cannot carry.
    struct HopLegSplit
    {
        bool valid = false;        ///< the last solve was HOP or SHP
        /// The call did not return - it left by exception. The iteration counts are still
        /// real work spent (the Optima path's own catch sets them before re-throwing), but
        /// timeOptima is NOT recoverable and is left at 0: this flag is what says that 0 is
        /// not a measurement. A caller summing splits over many solves must read it, or it
        /// is summing only the calls that returned.
        bool failed = false;
        bool warmNative = false;   ///< SHP (warm native leg) rather than HOP
        bool nativeOk = false;     ///< the native leg produced an answer to warm-start from
        long int fiaNative = 0;    ///< MBR iterations, native leg
        long int ipmNative = 0;    ///< IPM descent iterations, native leg
        long int fiaOptima = 0;    ///< MBR-equivalent iterations reported by the Optima leg
        long int ipmOptima = 0;    ///< Optima iterations (including a failed, discarded attempt)
        double timeNative = 0.;    ///< seconds in the native leg, as CalculateEquilibriumState() reports
        double timeOptima = 0.;    ///< seconds in the Optima leg
    };
    HopLegSplit hop_split;

    MULTI pm;
    MULTI *pmp;

    // Internal arrays for the performance optimization  (since version 2.0.0)
    long int sizeN; /*, sizeL, sizeAN;*/
    double *AA;
    double *BB;
    long int *arrL;
    long int *arrAN;

    void Alloc_A_B( long int newN );
    void Free_A_B();
    void Build_compressed_xAN();
    void Free_compressed_xAN();
    void Free_internal();

     // From here move to activities.h or node.h
    long int sizeFIs;     ///< current size of phSolMod
    TSolMod** phSolMod; ///< size current FIs - number of multicomponent phases
    void Alloc_TSolMod( long int newFIs );
    void Free_TSolMod();

    // new - allocation of TsorpMod and TKinMet
    long int sizeFIa;       ///< current size of phSorpMod
    TSorpMod** phSorpMod; ///< size current FIa - number of adsorption phases

    void Alloc_TSorpMod( long int newFIa );
    void Free_TSorpMod();

    long int sizeFI;      ///< current size of phKinMet
    TKinMet** phKinMet; ///< size current FI -   number of phases

    void Alloc_TKinMet( long int newFI );
    void Free_TKinMet();
    // until here move to activities.h or node.h

    // Added for implementation of divergence detection in dual solution 06.05.2011 DK
    long int nNu;  ///< number of ICs in the system
    long int cnr;  ///< current IPM iteration
    long int nCNud; ///< number of IC names for divergent dual chemical potentials
    /// pa_IpmAugmentedKKT: zero-row rescues in the current InteriorPointsMethod() call (reset at its
    /// entry). Only the first one per call emits DECIDE "ipmkkt-zerorow", so a rescue repeated on every
    /// iteration does not flood the trace.
    long int ipmKktRescues = 0;
    /// pa_IpmAugmentedKKT = 1 or 2: the main-loop (initAppr = false) solve for pm.U without forming
    /// A^T W A. \return 0 solved, 2 singular - the caller then takes that step by the normal
    /// equations (DECIDE ipmkkt-fallback). See BASE_PARAM::IpmAugmentedKKT.
    long int SolveIpmAugmented( long int N );
    double *U_mean; ///< Cumulative mean dual solution approximation [nNu]
    double *U_M2;   ///< Cumulative sum of squares [nNu]
    double *U_CVo;  ///< Cumulative Coefficient of Variation for dual solution approximation r-1 [nNu]
    double *U_CV;   ///< Cumulative Coefficient of Variation for r-th dual solution approximation [nNu]
    long int *ICNud; ///< List of IC indexes for divergent dual chemical potentials [nNu]
    void Alloc_uDD( long int newN );
    void Free_uDD();
    void Reset_uDD(long int cr, bool trace = false );
    void Increment_uDD( long int r, bool trace = false );
    long int Check_uDD( long int mode, double DivTol, bool trace = false );

    virtual void get_PAalp_PSigm(char &PAalp, char &PSigm);
    virtual void STEP_POINT( const char* /*str*/);
    virtual void alloc_IPx( long int LsIPxSum );
    virtual void alloc_PMc( long int LsModSum );
    virtual void alloc_DMc( long int LsMdcSum );
    virtual void alloc_MoiSN( long int LsMsnSum );
    virtual void alloc_SitFr( long int LsSitSum );
    virtual void alloc_DQFc( long int DQFcSum );
    virtual void alloc_PhLin( long int PhLinSum );
    virtual void alloc_lPhc( long int lPhcSum );
    virtual void alloc_xSMd( long int xSMdSum );
    virtual void alloc_IsoPc( long int IsoPcSum );
    virtual void alloc_IsoSc( long int IsoScSum );
    virtual void alloc_IsoCt( long int IsoCtSum );
    virtual void alloc_EImc( long int EImcSum );
    virtual void alloc_mCDc( long int mCDcSum );
    virtual void alloc_xSKrC( long int xSKrCSum );
    virtual void alloc_ocPRkC( long int ocPRkC_feSArC_Sum );
    virtual void alloc_feSArC( long int ocPRkC_feSArC_Sum );
    virtual void alloc_rpConC( long int rpConCSum );
    virtual void alloc_apConC( long int apConCSum );
    virtual void alloc_AscpC( long int AscpCSum );
    virtual void alloc_UMpcC( long int UMpcSum );
    virtual void alloc_xICuC( long int xICuCSum );

    virtual void loadData( bool ){}
    virtual bool testTSyst() const;
    virtual bool calculateActivityCoefficients_scripts( long int, long, long, long, long, long, double );
    virtual void initalizeGEM_IPM_Data_GUI();
    virtual void multiConstInit_PN();
    virtual void GEM_IPM_Init_gui1();
    virtual void GEM_IPM_Init_gui2();
   
    void setErrorMessage( long int num, const char *code, const char * msg);
    void addErrorMessage( const char * msg);

    // ipm_chemical.cpp
    void XmaxSAT_IPM2();
    //    void XmaxSAT_IPM2_reset();
    double DC_DualChemicalPotential( double U[], double AL[], long int N, long int j );
    void Set_DC_limits( bool InitState );
    void TotalPhasesAmounts( double X[], double XF[], double XFA[] );
    double DC_PrimalChemicalPotentialUpdate( long int j, long int k );
    double  DC_PrimalChemicalPotential( double G,  double logY,  double logYF,
                                        double asTail,  double logYw,  char DCCW );
    void PrimalChemicalPotentials( double F[], double Y[],
                                   double YF[], double YFA[] );
    double KarpovCriterionDC( double *dNuG, double logYF, double asTail,
                              double logYw, double Wx,  char DCCW );
    void KarpovsPhaseStabilityCriteria();
    void  StabilityIndexes( );
    double DC_GibbsEnergyContribution(   double G,  double x,  double logXF,
                                         double logXw,  char DCCW );
    double GX( double LM  );
    void ConvertDCC();
    long int  getXvolume();

    // ipm_chemical2.cpp
    virtual void GasParcP(){}
    void phase_bcs( long int N, long int M, long int jb, double *A, double X[], double BF[] );
    void phase_bfc( long int k, long int jj );
    double bfc_mass( void );
    void CalculateConcentrationsInPhase( double X[], double XF[], double XFA[],
                                         double Factor, double MMC, double Dsur, long int jb, long int je, long int k );
    void CalculateConcentrations( double X[], double XF[], double XFA[]);
    void IS_EtaCalc();
    long int GouyChapman(  long int jb, long int je, long int k );
    //  Surface activity coefficient terms
    long int SurfaceActivityCoeff( long int jb, long int je, long int jpb, long int jdb, long int k );

    // ipm_chemical3.cpp
    double SmoothingFactor( );
    void SetSmoothingFactor( long int mode ); // new smoothing function (3 variants)
    // Main call for calculation of activity coefficients on IPM iterations
    long int CalculateActivityCoefficients( long int LinkMode );
    // Built-in activity coefficient models
    // Generic solution model calls
    void SolModCreate( long int jb, long int jmb, long int jsb, long int jpb, long int jdb,
                       long int k, long int ipb, char ModCode, char MixCode,
                       /* long int jphl, long int jlphc, */ long int jdqfc/*, long int jrcpc*/ );
    void SolModParPT( long int k, char ModCode );
    void SolModActCoeff( long int k, char ModCode );
    void SolModExcessProp( long int k, char ModCode );
    void SolModIdealProp ( /*long int jb,*/ long int k, char ModCode );
    void SolModStandProp ( /*long int jb,*/ long int k, char ModCode );
    void SolModDarkenProp ( /*long int jb,*/ long int k/*, char ModCode*/ );

    // Specific phase property calculation functions  // obsolete (29.11.10 TW)
    // void IdealGas( long int jb, long int k, double *Zid );
    // void IdealOneSite( long int jb, long int k, double *Zid );
    // void IdealMultiSite( long int jb, long int k, double *Zid );

    // ipm_chemical4.cpp
    // New stuff for TKinMet class implementation
    long int CalculateKinMet( long int LinkMode  );

    /// One kinetics/metastability time step (TKinMet), for EVERY solver path. See the definition in
    /// ipm_chemical4.cpp for why its position relative to ExcludeRedundantDCs() and
    /// ScaleSystemToInternal() is fixed, and why HOP's Optima leg must not run it.
    void RunKineticsStep();
    void KM_Create(long int jb, long int k, long int kc, long int kp, long int kf,
                   long int ka, long int ks, long int kd, long ku, long ki, const char *kmod,
                   long jphl, long jlphc );
    void KM_ParPT( long int k, const char *kMod );
    void KM_InitTime( long int k, const char *kMod );
    void KM_UpdateTime( long int k, const char *kMod );
    void KM_UpdateFSA(long jb, long int k, const char *kMod );
    void KM_ReturnFSA(long int k, const char *kMod );
    void KM_CalcRates( long int k, const char *kMod );
    void KM_InitRates( long int k, const char *kMod );
    void KM_CalcSplit( /*long int jb,*/ long int k, const char *kMod );
    void KM_InitSplit( /*long int jb,*/ long int k, const char *kMod );
    void KM_CalcUptake( /*long int jb,*/ long int k, const char *kMod );
    void KM_InitUptake( /*long int jb,*/ long int k, const char *kMod );
    void KM_SetAMRs( /*long int jb,*/ long int k, const char *kMod );

    // ipm_main.cpp - numerical part of GEM IPM-2
    void GibbsEnergyMinimization();
    void GEM_IPM( long int rLoop );
    long int MassBalanceRefinement( long int WhereCalledFrom );

    /// pa_MbReproject: project pm.X back onto A.X = b with one N x N solve over a
    /// RANK-REVEALING set of the most abundant species (White 1958's Table III note,
    /// with his own "m most abundant" pivot rule corrected - it is exactly singular
    /// on aqueous chemistry; see BASE_PARAM::MbReproject). Native path only, called
    /// once at the end of GibbsEnergyMinimization() and only when the answer has
    /// already failed its own per-IC mass-balance test.
    /// Returns true only if it applied a correction that left every species
    /// non-negative AND strictly reduced the worst relative residual; otherwise it
    /// restores pm.X untouched and returns false.
    /// `amt` is the amount vector to repair - pm.X on the final answer, pm.Y at the
    /// post-PSSC call site (PSSC works on pm.Y). Both are re-synchronised on success.
    bool MassBalanceReproject( double* amt );
    /// Warn when a present phase's AMOUNT is not determined by the minimised energy:
    /// the Gibbs energy is flat enough along a mass-balance-preserving direction that
    /// answers differing in that phase's amount cannot be told apart at the solver's own
    /// energy resolution. Read-only - never changes the answer. See ipm_main.cpp.
    void EnergyDeterminacyCheck();
    /// A species ExcludeRedundantDCs() holds at zero for the current call, with the
    /// metastability settings RestoreRedundantDCs() puts back. (A nested type and two
    /// non-virtual members: no change to the class layout.)
    struct RedundantDCHold { long int j; char rlc; double dll, dul; };
    /// Find REDUNDANT species - identical stoichiometry, class and standard properties at
    /// the current T,P, either twice in one phase or as two single-species phases - warn,
    /// and hold every copy after the first at zero for this call (internal DLL = DUL = 0,
    /// RLC = BOTH_LIM; any starting amount moved onto the kept one). The caller's input
    /// (DATABR dll/dul) is not touched. See ipm_main.cpp.
    std::vector<RedundantDCHold> ExcludeRedundantDCs();
    void RestoreRedundantDCs( const std::vector<RedundantDCHold>& held );
    /// Warn when an element can exist only in ONE multi-component phase and that phase is a trace
    /// amount made up largely of the element (it is held open by it) - a fragile system definition.
    /// Read-only. See ipm_main.cpp.
    void StrandedElementCheck();
    long int InteriorPointsMethod( long int &status/*, long int rLoop*/ );
    void AutoInitialApproximation( );

    // ipm_main.cpp - miscellaneous fuctions of GEM IPM-2
    void MassBalanceResiduals( long int N, long int L, double *A, double *Y,
                               double *B, double *C );
    double OptimizeStepSize( double LM );
    void DC_ZeroOff( long int jStart, long int jEnd, long int k=-1L );
    void DC_RaiseZeroedOff( long int jStart, long int jEnd, long int k=-1L );
    /// pa_LpDualFillout (default 0 = off): size the species the LP zeroed from the LP's own dual
    /// instead of the per-class constants. Native cold path. See the field's doc comment in
    /// BASE_PARAM for every measured number and why it is off; the definition in ipm_main.cpp for
    /// the formula and its three guards (big-M, composition ceiling, class floor).
    void LpDualFillout( const std::vector<double>& yLp );
    /// pa_FilloutBudget (default 0.01): scale the class fill-out so it perturbs each element's
    /// mass balance by at most that fraction of the element's own bulk amount. See the field's doc
    /// comment in BASE_PARAM for the measurement and why a per-species cap cannot do this job.
    void ApplyFilloutBudget( const std::vector<double>& yLp );
    /// The effective pa_FilloutBudget - the field, unless GEMS3K_FILLOUT_BUDGET overrides it.
    /// The call site and the mechanism must both read THIS, never the field directly.
    double FilloutBudgetValue() const;
    /// Effective fill-out mode: the field, unless GEMS3K_LPDUAL_FILLOUT overrides it.
    long int LpFilloutMode() const;
    /// Record of every prediction against the amount the solve converged to, written at the ANSWER
    /// site. Zero cost unless GEMS3K_LPFILL_PROBE is set; this is what keeps plan v5 137.4's and
    /// 137.8's rejecting measurements reproducible.
    void LpFillProbeReport();
    /// Largest amount of a single-species phase the bulk composition can supply,
    /// min_i b_i/a(j,i) over the ordinary IC rows. Used to clamp PSSC's fixed
    /// pure-phase insertion amount (pa_DFYs) so an insertion cannot be infeasible
    /// by construction - see the long comment at the definition in ipm_chemical.cpp.
    double PhaseInsertionCeiling( long int j );
    double RaiseDC_Value( const long int j );
    long int MetastabilityLagrangeMultiplier();
    void WeightMultipliers( bool square );
    long int MakeAndSolveSystemOfLinearEquations( long int N, bool initAppr );
    double DikinsCriterion(  long int N, bool initAppr );
    double StepSizeEstimate( bool initAppr );
    void Restore_Y_YF_Vectors();
    double RescaleToSize( bool standard_size ); // replaced calcSfactor() 30.08.2009 DK
    long int SpeciationCleanup( double AmountThreshold, double ChemPotDiffCutoff ); // added 25.03.10 DK
    long int PhaseSelectionSpeciationCleanup( long int &k_miss, long int &k_unst, long int rLoop );
    long int PhaseSelect( long int &k_miss, long int &k_unst, long int rLoop );
    bool GEM_IPM_InitialApproximation();

    // IPM_SIMPLEX.CPP Simplex method modified with two-sided constraints (Karpov ea 1997)
    void SolveSimplex(long int M, long int N, long int T, double GZ, double EPS,
                      double *UND, double *UP, double *B, double *U,
                      double *AA, long int *STR, long int *NMB );
    void SPOS( double *P, long int STR[],long int NMB[],long int J,long int M,double AA[]);
    void START( long int T,long int *ITER,long int M,long int N,long int NMB[],
                double GZ,double EPS,long int STR[],long int *BASE,
                double B[],double UND[],double UP[],double AA[],double *A,
                double *Q );
    void NEW(long int *OPT,long int N,long int M,double EPS,double *LEVEL,long int *J0,
             long int *Z,long int STR[], long int NMB[], double UP[],
             double AA[], double *A);
    void WORK(double GZ,double EPS,long int *I0, long int *J0,long int *Z,long int *ITER,
              long int M, long int STR[],long int NMB[],double AA[],
              long int BASE[],long int *UNO,double UP[],double *A,double Q[]);
    void FIN(double EPS,long int M,long int N,long int STR[],long int NMB[],
             long int BASE[],double UND[],double UP[],double U[],
             double AA[],double *A,double Q[],long int *ITER);
    double SystemTotalMolesIC( );
    void ScaleSystemToInternal(  double ScFact );
    void RescaleSystemFromInternal(  double ScFact );
    void MultiConstInit(); // from MultiRemake
    void GEM_IPM_Init();

    virtual void load_all_thermodynamic_from_grid(TNode *aNa, double TK, double P);

public:
    /// IC names the caller marked "of interest" (TNode::GEM_set_elements_of_interest()); every other trace IC
    /// is a default seed. Read by StrandedElementCheck() and TNode::GEM_trace_regimes(); never by the solve.
    /// TRAILING member on purpose: TMultiBase is allocated only inside the library (TNode::allocMemory(),
    /// GEM_trace_regimes()), so appending here moves no offset an external binary uses. Copied by copyMULTI().
    std::vector<std::string> elementsOfInterest;

    /// ELEMENT CLASSES for decisions on the Optima path (owner 2026-09-15, DECISIONS.md; HANDOFF-2026-09-15 s4).
    /// An ordinary IC i (not a charge row) is NUMERICALLY TRACE when one numerical-floor amount does not fit
    /// inside its own mass-balance tolerance: B_i*pa_DHB < floor, the floor being pa_OptimaDcFloor if set, else
    /// pa_DHB, in the solve's internal units. With pa_DG > 1e-5 the solve is rescaled to sum(B) = pa_DG, so with
    /// the default floor the test is B_i/sum(B) < 1/pa_DG - unit-free, hence the same on internal and on real
    /// amounts; with pa_DG <= 1e-5 it is B_i < 1 mol. Measured boundary sum(B)/pa_DG: median 0.18 mol over 57
    /// projects; 269 of 1312 ICs major. No new parameter: the boundary moves with each project's own settings.
    /// This is a property of the OPTIMA path's numerics - native's fixed amounts are pa_DcMin-scale and make
    /// every corpus IC major at shipped settings. The dilute-regime check (GEM_trace_regimes) keeps its own
    /// chemical relative rule and does not read this.
    bool ICIsNumericalTrace( long int i ) const;
    /// A DEFAULT SEED is a numerically trace IC the caller did NOT mark of interest, once any IC is marked.
    /// With nothing marked there are no default seeds: every IC is conserved (owner 2026-09-14e). Consumers:
    /// pa_OptimaZeroAbsent's rebalance/undo test (a default seed does not trigger it) and the CERT record
    /// (default seeds reported as mb_seed_rel, not scored in mb_pass).
    bool ICIsDefaultSeed( long int i ) const;
    /// SUB-FLOOR ELEMENTS on the Optima path (owner 2026-09-15: "produce a warning", with information for the user).
    /// Warns when an ordinary IC's bulk amount is smaller than the least its carriers can hold,
    ///   need_i = sum_j a(i,j) * max(DLL_j, floor)   over species with a(i,j) > 0,
    /// i.e. no point inside Optima's box satisfies that IC's balance. The solve still runs; the IC's residual is left
    /// to the post-solve repair, and its dual is free along e_i. Measured on 3Bent-H2O_G_Mont_0_0_1_25_0: Nit 1e-15 mol
    /// against a 1.67e-14 mol floor; survey_dualfree reads the live direction (Nit +1.000) on the AOP answer; under the
    /// unconditional zeroing (2c390db) every Optima row returned mb_rel = 1e13 on Nit - the whole Nit budget missing -
    /// and with the phase-level zeroing + rebalance (plan v5 s123.7) it passes at 1.6e-2. By the formula, not measured:
    /// xGEMS' Material seeds elements at 1e-15 mol, which is sub-floor on this path once total bulk exceeds ~10 mol per
    /// carrier at pa_DHB = 1e-13, pa_DG = 1000. Read-only; called once per Optima call after rescaling, so amounts are
    /// internal and the message converts them to real moles. Native modes are unaffected (their floor is pa_DcMin-scale).
    void SubFloorElementCheck( double dcFloor ) const;
};

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
/// Env-gated event trace of the NATIVE (IPM/MBR/PSSC) solver.
///
/// Returns the open trace stream when GEMS3K_NATIVE_TRACE_FILE=<path> is set in
/// the environment, and nullptr otherwise - so every call site costs one null
/// test on an already-computed pointer. The stream is opened once, in append
/// mode, on first use and deliberately never closed (process lifetime).
///
/// DELIBERATELY NOT #ifndef NDEBUG. GEMS3K's existing staged-snapshot and PSSC
/// logging (added 2026-07-31, ipm_main.cpp) is NDEBUG-gated, and every installed
/// build - gems-benchmark's included - is Release, so that code does not exist
/// in the binary anyone actually runs and has never once been used in this
/// work. Env gating follows GEMS3K_OPTIMA_TRACE_FILE instead, which is usable
/// against a shipped library with no rebuild.
///
/// WHY IT EXISTS. Every statement this branch makes about native is an
/// INFERENCE from outcomes - iteration counts, exact zeros in the answer,
/// residuals - never an observation. The single largest open item is porting a
/// PSSC equivalent to the Optima path, and a decision procedure cannot be
/// ported from its output alone. What this logs is exactly that decision: which
/// phase, which stability index, against which threshold, and what was done.
///
/// Format: one line per EVENT (never per iteration - the psina giants carry
/// 1392 species and a per-iteration dump would be tens of MB), tagged in the
/// first field and otherwise KEY=value so it is greppable and parses on
/// whitespace. Names, not bare indices: the packed arrays pm.SM[]/pm.SF[]/
/// pm.SB[] are reachable from this layer via char_array_to_string(), which is
/// what ipm_optima.cpp already does for its own messages.
FILE* native_trace_file();
/// Suppress (on = true) / re-enable the per-call RUN header and KEY result records; nests.
void native_trace_quiet( bool on );

/// Per-iteration IPM descent record, gated on GEMS3K_IPM_PROBE=<path>. Same
/// zero-cost-when-unset shape as native_trace_file(); see the call site in
/// InteriorPointsMethod() for the measurement it exists to support.
FILE* ipm_probe_file();

/// Per-prediction record for pa_LpDualFillout, gated on GEMS3K_LPFILL_PROBE=<path>. Same
/// zero-cost-when-unset shape as the two above. This is what keeps plan v5 137.4's and 137.8's
/// REJECTING measurements reproducible after the mechanism itself was rejected as a default.
FILE* lpfill_probe_file();

/// Complete run configuration - requested mode, T, P, bulk composition with IC
/// names, and every BASE_PARAM field in force - written into the same file as
/// native_trace_file(), once per GEM_run() call, for EVERY solver mode. Emitted
/// from TNode::GEM_run() rather than from a solver, so the record is uniform
/// across native/AOP/SOP/ROP/HOP/SHP and appears exactly once per call. See the
/// definition in ipm_main.cpp for the format and for why the whole parameter set
/// is dumped rather than the handful the CALL record carries.
void native_trace_run_header( const MULTI& pm, const BASE_PARAM* pa, long int mode );

/// One DECIDE record: a choice the solver MADE during a call, as opposed to what
/// it was configured with (the SET/EFF lines) or what it returned (ANSWER/KEY).
///
/// WHY THIS EXISTS. Two silent wrong-answer channels were found on 2026-09-07c
/// and both had the same shape - the record described something other than what
/// ran (plan v5 sections 95.4 and 96.2). A freeze that carries the configuration
/// and the answer still cannot say WHICH MECHANISMS FIRED, and that is where both
/// bugs were visible: "dimension-reduction pre-solve produced a warm start over
/// 330 of 1287 species in 1117 iterations" was in the log the whole time, and no
/// artefact anyone diffs contained it.
///
/// These went into the trace rather than being scraped from the spdlog stream
/// because the trace is already mode-attributed by position - every DECIDE
/// follows the RUN line of the call that produced it - while log text carries no
/// mode and would have to be re-attributed by a parser that can drift.
///
/// Format: `DECIDE <what> <key>=<value> ...`, key=value so it stays greppable and
/// diffable. NO TIMESTAMPS and no wall time: this is frozen and diffed, and
/// anything that changes between two identical runs makes every diff noise.
/// Zero cost when GEMS3K_NATIVE_TRACE_FILE is unset.
void native_trace_decide( const char* fmt, ... );

/// The OUTCOME KEY of a completed solve, into the same trace file: the present
/// phase assemblage (names, amounts and molar volumes, plus an order-independent
/// hash of the name set), pH, pe and ionic strength. Emitted once per
/// TNode::GEM_run() call, after the dispatch, so a trace carries both the
/// configuration a result was produced at (RUN/BULK/SET) and the regime it
/// reached. Exists for work item 7 / plan v5 section 81.5: a setting cannot be
/// predicted from the INPUT, but may be looked up on the regime the input LEADS
/// TO - which makes the key an output and the lookup memoisation. The fluid-root
/// axis that design also names is deliberately NOT classified here; the per-phase
/// molar volume it would be built from is emitted instead. See the definition
/// (ipm_main.cpp) for why.
void native_trace_run_result( const MULTI& pm, long int mode, long int status, TMultiBase* mb = nullptr );

/// Optima's free-dual search result, stashed for the CERT record (Phase 3 WP1).
/// CertDualFreeDirs() counts the free directions on any path, but only the Optima path
/// SEARCHES them for a dual that satisfies every bound (ipm_optima.cpp, the `dualfree`
/// DECIDE), and only when its plain sign test has already failed. So `dual_resolved` is
/// three-valued: -1 = no search was run on this call, 0 = searched and did not resolve,
/// 1 = searched and resolved. Reset at the run header so a value never carries over from
/// the previous call - the stale-value trap this record exists to expose.
void native_cert_dualfree_reset();
void native_cert_dualfree_set( int resolved );
int  native_cert_dualfree_get();

// ???? syp->PGmax
typedef enum {  // Symbols of thermodynamic potential to minimize
    G_TP    =  'G',   // Gibbs energy minimization G(T,P)
    A_TV    =  'A',   // Helmholts energy minimization A(T,V)
    U_SV    =  'U',   // isochoric-isentropicor internal energy at isochoric conditions U(S,V)
    H_PS    =  'H',   // isobaric-isentropic or enthalpy H(P,S)
    _S_PH   =  '1',   // negative entropy at isobaric conditions and fixed enthalpy -S(P,H)
    _S_UV   =  '2'    // negative entropy at isochoric conditions and fixed internal energy -S(P,H)

} THERM_POTENTIALS;

typedef enum {  // Symbols of thermodynamic potential to minimize
    G_TP_    =  0,   // Gibbs energy minimization G(T,P)
    A_TV_    =  1,   // Helmholts energy minimization A(T,V)
    U_SV_    =  2,   // isochoric-isentropicor internal energy at isochoric conditions U(S,V)
    H_PS_    =  3,   // isobaric-isentropic or enthalpy H(P,S)
    _S_PH_   =  4,   // negative entropy at isobaric conditions and fixed enthalpy -S(P,H)
    _S_UV_   =  5    // negative entropy at isochoric conditions and fixed internal energy -S(P,H)

} NUM_POTENTIALS;

typedef enum {  // Field index into outField structure
    f_pa_PE = 0,  f_PV,  f_PSOL,  f_PAalp,  f_PSigm,
    f_Lads,  f_FIa,  f_FIat
} MULTI_STATIC_FIELDS;

typedef enum {  // Field index into outField structure
    f_sMod = 0,  f_LsMod,  f_LsMdc,  f_B,  f_DCCW,
    f_Pparc,  f_fDQF,  f_lnGmf,  f_RLC,  f_RSC,
    f_DLL,  f_DUL,  f_Aalp,  f_Sigw,  f_Sigg,
    f_YOF,  f_Nfsp,  f_MASDT,  f_C1,  f_C2,
    f_C3,  f_pCh,  f_SATX,  f_MASDJ,  f_SCM,
    f_SACT,  f_DCads,
    // static
    f_pa_DB,  f_pa_DHB,  f_pa_EPS,  f_pa_DK,  f_pa_DF,
    f_pa_DP,  f_pa_IIM,  f_pa_PD,  f_pa_PRD,  f_pa_AG,
    f_pa_DGC,  f_pa_PSM,  f_pa_GAR,  f_pa_GAH,  f_pa_DS,
    f_pa_XwMin,  f_pa_ScMin,  f_pa_DcMin,  f_pa_PhMin,  f_pa_ICmin,
    f_pa_PC,  f_pa_DFM,  f_pa_DFYw,  f_pa_DFYaq,  f_pa_DFYid,
    f_pa_DFYr,  f_pa_DFYh,  f_pa_DFYc,  f_pa_DFYs,  f_pa_DW,
    f_pa_DT,  f_pa_GAS,  f_pa_DG,  f_pa_DNS,  f_pa_IEPS,
    f_pKin,  f_pa_DKIN,  f_mui,  f_muk,  f_muj,
    f_pa_PLLG,  f_tMin,  f_dcMod,
    //new
    f_kMod, f_LsKin, f_LsUpt, f_xICuC, f_PfFact,
    f_LsESmo, f_LsISmo, f_SorMc, f_LsMdc2, f_LsPhl,
    f_pa_PSTALL, f_pa_OptimaTol, f_pa_LogBarrierTau, f_pa_OptimaMaxStepRatio,
    f_pa_PhaseHessianFloor, f_pa_OptimaStallWindow, f_pa_OptimaMaxSeconds,
    f_pa_OptimaFDHessian, f_pa_OptimaMoleFracHessian, f_pa_OptimaPhaseCompaction,
    f_pa_OptimaFDHessianDelay, f_pa_OptimaDcFloor, f_pa_MbClassRule,
    f_pa_MbTrendPhaseDecay, f_pa_OptimaEarlyStabilityAt, f_pa_OptimaDimReduce,
    f_pa_OptimaDimReduceTol, f_pa_MbPivotSplit, f_pa_OptimaZeroAbsent,
    f_pa_OptimaReadmitSeed,
    f_pa_IpmStallWindow, f_pa_MbReproject, f_pa_DeterminacyWarn, f_pa_ColdRetryNudges,
    f_pa_OptimaPreSolveFirstIters, f_pa_LpDualFillout, f_pa_FilloutBudget,
    f_pa_StabTPD, f_pa_IpmAugmentedKKT, f_pa_IpmLoopTweaks,
    f_pa_OptimaLineSearch, f_pa_OptimaFDDiagFloor, f_pa_OptimaLSStallEscape, f_pa_OptimaLSWindow,
    f_pa_OptimaLSRejectWorse,
    f_pa_OptimaTpdAccept, f_pa_OptimaCgSeed, f_pa_OptimaColdRetry, f_pa_OptimaFinish,
    f_pa_OptimaAcceptRepair

} MULTI_DYNAMIC_FIELDS;

enum volume_code {  /* Codes of volume parameter ??? */
    VOL_UNDEF, VOL_CALC, VOL_CONSTR
};

#endif   //_ms_multi_h

