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

#include <cstdio>

class GemDataStream;
class TProfil;
class TNode;

const int  QPSIZE = 180, // earlier 20, 40 SD oct 2005
           QDSIZE = 60;

// Physical constants - see m_param.cpp or ms_param.cpp
extern const double R_CONSTANT, NA_CONSTANT, F_CONSTANT,
    e_CONSTANT,k_CONSTANT, cal_to_J, C_to_K, lg_to_ln, ln_to_lg, H2O_mol_to_kg, Min_phys_amount;

struct BASE_PARAM /// Flags and thresholds for numeric modules
{
   short
           PC,   ///< Mode of PhaseSelect() operation ( 0 1 2 ... ) { 1 }
           PD,   ///< abs(PD): Mode of execution of CalculateActivityCoefficients() functions { 2 }.
                 ///< Modes: 0-invoke, 1-at MBR only, 2-every MBR it, every IPM it. 3-not MBR, every IPM it.
                 ///< if PD < 0 then use test qd_real accuracy mode
           PRD,  ///< Since r1583/r409: Disable (0) or activate (-5 or less) the SpeciationCleanup() procedure { -5 }
           PSM,  ///< Level of diagnostic messages: 0- disabled (no ipmlog file); 1- errors; 2- also warnings 3- uDD trace { 1 }
           DP,   ///< Maximum allowed number of iterations in the MassBalanceRefinement() procedure {  30 }
           DW,   ///< Since r1583: Activate (1) or disable (0) error condition when DP was exceeded { 1 }
           DT,   ///< Since r1583/r409: DHB is relative for all (0) or absolute (-6 or less ) cutoff for major ICs { 0 }
           PLLG, ///< IPM tolerance for detecting divergence in dual solution { 10; range 1 to 1000; 0 disables the detection }
           PE,   ///< Flag for using electroneutrality condition in GEM IPM calculations { 0 1 }
           IIM   ///< Maximum allowed number of iterations in the MainIPM_Descent() procedure up to 9999 { 1000 }
           ;
         double DG,   ///< Standart total moles { 1e5 }
           DHB,  ///< Maximum allowed relative mass balance residual for Independent Components ( 1e-9 to 1e-15 ) { 1e-10 }
           DS,   ///< Cutoff minimum mole amount of stable Phase present in the IPM primal solution { 1e-12 }
           DK,   ///< IPM-2 convergence threshold for the Dikin criterion (may be set in the interval 1e-6 < DK < 1e-4) { 1e-5 }
           DF,   ///< Threshold for the application of the Karpov phase stability criterion: (Fa > DF) for a lost stable phase { 0.01 }
           DFM,  ///< Threshold for Karpov stability criterion f_a for insertion of a phase (Fa < -DFM) for a present unstable phase { 0.1 }
           DFYw, ///< Insertion mole amount for water-solvent { 1e-6 }
           DFYaq,///< Insertion mole amount for aqueous species { 1e-6 }
           DFYid,///< Insertion mole amount for ideal solution components { 1e-6 }
           DFYr, ///< Insertion mole amount for major solution components { 1e-6 }
           DFYh, ///< Insertion mole amount for minor solution components { 1e-6 }
           DFYc, ///< Insertion mole amount for single-component phase { 1e-6 }
           DFYs, ///< Insertion mole amount used in PhaseSelect() for a condensed phase component  { 1e-7 }
           DB,   ///< Minimum amount of Independent Component in the bulk system composition (except charge "Zz") (moles) (1e-17)
           AG,   ///< Smoothing parameter for non-ideal increments to primal chemical potentials between IPM descent iterations { -1 }
           DGC,  ///< Exponent in the sigmoidal smoothing function, or minimal smoothing factor in new functions { -0.99 }
           GAR,  ///< Initial activity coefficient value for major (M) species in a solution phase before LPP approximation { 1 }
           GAH,  ///< Initial activity coefficient value for minor (J) species in a solution phase before LPP approximation { 1000 }
           GAS,  ///< Since r1583/r409: threshold for primal-dual chem.pot.difference (mol/mol) used in SpeciationCleanup() { 1e-3 }.
                 ///< before: Obsolete IPM-2 balance accuracy control ratio DHBM[i]/b[i], for minor ICs { 1e-3 }
           DNS,  ///< Standard surface density (nm-2) for calculating activity of surface species (12.05)
           XwMin,///< Cutoff mole amount for elimination of water-solvent { 1e-9 }
           ScMin,///< Cutoff mole amount for elimination of solid sorbent {1e-7}
           DcMin,///< Cutoff mole amount for elimination of solution- or surface species { 1e-30 }
           PhMin,///< Cutoff mole amount for elimination of  non-electrolyte solution phase with all its components { 1e-10 }
           ICmin,///< Minimal effective ionic strength (molal), below which the activity coefficients for aqueous species are set to 1. { 3e-5 }
           EPS,  ///< Precision criterion of the SolveSimplex() procedure to obtain the AIA ( 1e-6 to 1e-14 ) { 1e-10 }
           IEPS, ///< Convergence parameter of SACT calculation in sorption/surface complexation models { 0.01 to 0.000001, default 0.001 }
           DKIN; ///< Tolerance on the amount of DC with two-side metastability constraints  { 1e-7 }
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
    long int OptimaMoleFracHessian = 0;

    /// PSSC-EQUIVALENT PHASE COMPACTION for the Optima path: number of Newton
    /// iterations to spend on a cheap CLASSIFICATION pass before the real
    /// solve. 0 (default) = off, behaviour unchanged.
    ///
    /// WHY A CLASSIFICATION PASS AT ALL. Native drops absent phases from the
    /// active set on every call (PhaseSelectionSpeciationCleanup(), pa_PC=2);
    /// the Optima path has always carried every absent phase as a live
    /// box-constrained unknown for the whole solve. Measured (plan-v5 section
    /// 20): the ABSOLUTE COUNT of absent phases tracks this solver's cost far
    /// better than problem size - 07PSIna_G_mid_1 has 117 of 122 phases absent
    /// and runs 4184 iterations / 49.5 s against native's 124 / 19.6 ms.
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
    /// This does NOT reproduce native's mechanism, only its reported semantics.
    /// Native's species genuinely leave the problem mid-solve and stop
    /// obstructing it; these leave only the answer. Closing that half is the
    /// dimension reduction (pa_OptimaDimReduce), which omits species from the
    /// vector outright - but its final full-dimension verification pass puts
    /// them all back on the floor before the answer is written, which is the
    /// gap this field closes from the other end.
    ///
    /// DEFAULT 1 (on), set 2026-09-03 on the project owner's decision: matching
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
    long int OptimaZeroAbsent = 1;

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

    void write(GemDataStream& oss);
    void read(GemDataStream& iss);
};


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

extern const BASE_PARAM pa_p_;

// Data of MULTI
class TMultiBase
{
    char PAalp_; ///< Flag for using (+) or ignoring (-) specific surface areas of phases
    char PSigm_; ///< Flag for using (+) or ignoring (-) specific surface free energies
    std::shared_ptr<BASE_PARAM> pa_standalone;

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
    double CalculateEquilibriumStateOptima( long int& NumIterFIA, long int& NumIterIPM, bool reaktoroMode = false );

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

    /// Dual of the LINEARISED-Gibbs LP (min sum_j G0[j]*n_j s.t. A n = b,
    /// n >= 0), computed with the same simplex as LPFeasibilitySeed() and
    /// therefore just as independent of any solver state. Its value is that
    /// pricing a species against it, s_j = G0[j] - sum_i y_i A[i,j], is a
    /// meaningful measure of how far that species is from being stable -
    /// which the feasibility LP's own dual is not, being an artefact of
    /// "minimise total moles". Used by OptimaReducedPreSolve() to choose a
    /// generous initial active set. Returns false (leaving yOut untouched)
    /// if the LP fails or its dual does not verify against LP optimality.
    bool LPGibbsDual( std::vector<double>& yOut );

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
    /// \return index of the worst violating phase, or -1 if the assemblage
    ///         is self-consistent.
    long int WorstPhaseStabilityViolation( double presenceThreshold, double dcFloor,
                                           const char* exemptSpecies,
                                           double& violOut, bool& wasAbsentOut );

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
    /// \param iterationsOut  Optima iterations spent here, for pm.ITG.
    /// \param activeOut  size of the final active set, for logging.
    /// \return true if pm.Y[]/pm.U[] now carry a state worth warm-starting
    ///         from; false if the pre-solve was skipped or discarded, in which
    ///         case neither array was modified in a way the caller must undo.
    bool OptimaReducedPreSolve( long int maxPasses, double dcFloor, double dimTol,
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
    long int InteriorPointsMethod( long int &status/*, long int rLoop*/ );
    void AutoInitialApproximation( );

    // ipm_main.cpp - miscellaneous fuctions of GEM IPM-2
    void MassBalanceResiduals( long int N, long int L, double *A, double *Y,
                               double *B, double *C );
    double OptimizeStepSize( double LM );
    void DC_ZeroOff( long int jStart, long int jEnd, long int k=-1L );
    void DC_RaiseZeroedOff( long int jStart, long int jEnd, long int k=-1L );
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

/// Per-iteration IPM descent record, gated on GEMS3K_IPM_PROBE=<path>. Same
/// zero-cost-when-unset shape as native_trace_file(); see the call site in
/// InteriorPointsMethod() for the measurement it exists to support.
FILE* ipm_probe_file();

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
    f_pa_IpmStallWindow

} MULTI_DYNAMIC_FIELDS;

enum volume_code {  /* Codes of volume parameter ??? */
    VOL_UNDEF, VOL_CALC, VOL_CONSTR
};

#endif   //_ms_multi_h

