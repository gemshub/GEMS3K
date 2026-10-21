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

/// The `{ }` at the end of each field's comment below is the value a standalone GEMS3K run
/// gets when the project file does not set that field (the entry in `pa_p_`,
/// ms_multi_diff.cpp).
struct BASE_PARAM /// Flags and thresholds for numeric modules
{
   short
           PC,   ///< Mode of PhaseSelect() operation ( 0 1 2 ... ) { 2 }
           PD,   ///< abs(PD): Mode of execution of CalculateActivityCoefficients() functions { 2 }.
                 ///< Modes: 0-invoke, 1-at MBR only, 2-every MBR it, every IPM it. 3-not MBR, every IPM it.
                 ///< if PD < 0 then use test qd_real accuracy mode
           PRD,  ///< Since r1583/r409: Disable (0) or activate (-5 or less) the SpeciationCleanup() procedure { -5 }
                 ///< Also sets PSSC's AmountThreshold = 10^-|PRD| (at PRD = 0 the threshold is 1 mol).
                 ///< The cleanup loop runs only when pa_PC == 2 as well.
           PSM,  ///< Level of diagnostic messages: 0- disabled (no ipmlog file); 1- errors; 2- also warnings 3- uDD trace { 1 }
           DP,   ///< Maximum allowed number of iterations in the MassBalanceRefinement() procedure { 130 }
           DW,   ///< Since r1583: Activate (1) or disable (0) error condition when DP was exceeded { 1 }
           DT,   ///< Since r1583/r409: DHB is relative for all (0) or absolute (-6 or less ) cutoff for major ICs { 0 }
           PLLG, ///< IPM tolerance for detecting divergence in dual solution { 30000 }; 1 to 1000 is the working range, 0 disables the detection, |PLLG| >= 30000 also allows complete tracing
           PE,   ///< Flag for using electroneutrality condition in GEM IPM calculations ( 0 or 1 ) { 1 }
           IIM   ///< Maximum allowed number of iterations in the MainIPM_Descent() procedure, at most 32767 { 7000 }
           ;
         double DG,   ///< Standart total moles { 1000. }
           DHB,  ///< Maximum allowed relative mass balance residual for Independent Components ( 1e-9 to 1e-15 ) { 1e-13 }
           DS,   ///< Cutoff minimum mole amount of stable Phase present in the IPM primal solution { 1e-20 }
           DK,   ///< IPM-2 convergence threshold for the Dikin criterion { 1e-6 }.
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
           GAR,  ///< reserved, no effect (not read by GEMS3K; kept for the field order) { 1 }
           GAH,  ///< reserved, no effect (not read by GEMS3K; kept for the field order) { 1000 }
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

    /// pa_PSTALL: 1 (default) = MassBalanceRefinement() stops when it stalls; 0 = it keeps
    /// iterating until pa_DP is used up and then fails.
    /// Fields from here on are GEMS3K additions: they come after every other member and have
    /// default initializers, so the order of the older members is unchanged. Adding one also
    /// needs both counts in ms_multi_format.cpp and the SET line in native_trace_run_header().
    /// In plain words: lets the mass-balance step give up early when it stops improving.
    short PSTALL = 1;

    /// Optima solver settings (ipm_optima.cpp, USE_OPTIMA_SOLVER builds only). Trailing
    /// fields with default member initializers, like PSTALL; keyword-only I/O.
    ///
    /// pa_OptimaTol: Optima's own optimality-error convergence tolerance. Not the same
    /// quantity as pa_DK (the native Dikin criterion). Default 1e-8.
    /// In plain words: how precisely the Optima solver must converge.
    double OptimaTol = 1.0e-8;

    /// pa_LogBarrierTau: logarithmic-barrier penalty weight added to the Optima objective for
    /// pure single-species phases (species whose chemical potential does not depend on
    /// composition). Default 1e-16.
    /// In plain words: a tiny numerical push that keeps pure minerals away from exactly zero.
    double LogBarrierTau = 1.0e-16;

    /// pa_OptimaMaxStepRatio: reserved, no effect. Kept so that the field order of BASE_PARAM
    /// (used by positional project-file serialisation) does not change.
    /// In plain words: an old setting that no longer does anything.
    double OptimaMaxStepRatio = 0.0;

    /// pa_PhaseHessianFloor: lowest allowed curvature of a mixed (non-aqueous) phase in the Optima
    /// solver, as a fraction of its largest curvature. 0 = off. Default 0.01.
    /// In plain words: helps the Optima solver handle mixed phases that tend to split in two.
    double PhaseHessianFloor = 0.01;

    /// pa_OptimaStallWindow: an Optima solve is given up when its best error has not fallen
    /// (relative 1e-8) over this many iterations. 0 = off. Default 200. Values below 100 are
    /// not recommended.
    /// In plain words: stops a solve that has stopped making progress, so a retry can start sooner.
    long int OptimaStallWindow = 200;

    /// pa_OptimaMaxSeconds: time limit in seconds for one Optima solve with its retries.
    /// 0 = off (default). A limit makes results depend on the machine, so set it per project only.
    /// In plain words: a time limit for the Optima solver, off unless you set it.
    double OptimaMaxSeconds = 0.0;

    /// pa_OptimaFDHessian: 1 (default) = the Optima path computes the finite-difference
    /// curvature columns; 0 = skip them (much cheaper, usually the same answer).
    /// In plain words: computes a more exact curvature at extra cost; usually optional.
    long int OptimaFDHessian = 1;

    /// pa_OptimaMoleFracHessian: curvature of the ideal-mixing part of a mixed phase in the Optima
    /// path. 0 (default) = diagonal 1/X[j] only; 1 = the full mole-fraction form.
    /// In plain words: a more exact curvature for mixed phases, which helps some projects and
    /// slows others.
    long int OptimaMoleFracHessian = 0;

    /// pa_OptimaPhaseCompaction: reserved, no effect. Kept so that the field order of BASE_PARAM
    /// (used by positional project-file serialisation) does not change.
    /// In plain words: an old setting that no longer does anything.
    long int OptimaPhaseCompaction = 0;

    /// pa_OptimaFDHessianDelay: number of first iterations that skip the finite-difference
    /// curvature columns; if that attempt does not converge, the solve restarts with them on.
    /// 0 = off (default). Cold starts only.
    /// In plain words: try the fast, approximate method first and switch to the exact one
    /// only if needed.
    long int OptimaFDHessianDelay = 0;

    /// pa_OptimaDcFloor: lowest species amount in the Optima path, in the solver's internal
    /// (rescaled) units. 0 (default) = derive it from pa_DHB; > 0 = use this value.
    /// In plain words: how small an amount the Optima solver can represent for an absent species.
    double OptimaDcFloor = 0.;

    /// pa_MbClassRule: 0 = off (default). Above 0, mass balance in MassBalanceRefinement() is
    /// checked in absolute terms for major elements and in relative terms for trace elements;
    /// the value is the trace/major ratio. Set it per project.
    /// In plain words: checks major elements by absolute and trace elements by relative accuracy.
    double MbClassRule = 0.;

    /// pa_MbTrendPhaseDecay: reserved, no effect. Kept so that the field order of BASE_PARAM
    /// (used by positional project-file serialisation) does not change.
    /// In plain words: an old setting that no longer does anything.
    long int MbTrendPhaseDecay = 0;

    /// pa_OptimaEarlyStabilityAt: ends the first Optima attempt early so the phase-selection
    /// repair can act before the solve has fully converged.
    ///   > 0  stop at iteration N if the dual has settled.
    ///   < 0  stop once a non-solvent mixed phase has fallen for -N evaluations in a row.
    ///   0    AUTO (default): a cap when the system has a multisite solid solution, else off
    ///        (optima_earlystability_at()).
    /// If the cap stopped the attempt and nothing recovered, the call is re-solved in full. To
    /// switch it off, set a cap larger than pa_IIM.
    /// In plain words: lets the solver notice a phase that should disappear early, instead of
    /// after a long, slow dissolution.
    long int OptimaEarlyStabilityAt = 0;

    /// pa_OptimaDimReduce: before the full Optima solve, solve a reduced problem with only the
    /// species likely to be present, add back any that want to form, then use the result as the
    /// start of the full solve. > 0 = limit on add-back passes; 0 = AUTO, on for large systems
    /// (optima_dimreduce_passes()); < 0 = off. Not used with control conditions or in ROP; on the
    /// HOP leg only for an explicit positive value.
    /// In plain words: makes large systems much faster by first solving with only the
    /// species that are likely to be present.
    long int OptimaDimReduce = 0;

    /// pa_OptimaDimReduceTol: which species the reduced problem of pa_OptimaDimReduce starts with.
    ///   > 0  species whose reduced cost in the linear Gibbs problem is below this value (RT).
    ///   < 0  the cheapest species until the set holds |value| x N of them.
    ///   0    only the support of the starting guess.
    /// Default 10. It affects cost, not the answer.
    /// In plain words: how generous the first guess of "species that matter" is.
    double OptimaDimReduceTol = 10.;

    /// pa_MbPivotSplit: reserved, no effect. Kept so that the field order of BASE_PARAM
    /// (used by positional project-file serialisation) does not change.
    /// In plain words: an old setting that no longer does anything.
    long int MbPivotSplit = 0;

    /// pa_OptimaZeroAbsent: how species absent from an Optima answer are reported.
    ///   0  at the floor amount.
    ///   1  absent phases are set to exactly 0 and the element totals are rebalanced; the zeroing
    ///      is undone if the rebalance fails.
    ///   2  (default) floor amounts, plus a repair that projects the answer back onto the mass
    ///      balance if it fails the per-element test; the repair is dropped if it does not help.
    /// Any other non-zero value behaves as 1.
    /// In plain words: decides whether absent species are shown as exactly zero or as a tiny
    /// amount, and fixes the element totals if needed.
    long int OptimaZeroAbsent = 2;

    /// pa_OptimaReadmitSeed: reserved, no effect. Kept so that the field order of BASE_PARAM
    /// (used by positional project-file serialisation) does not change.
    /// In plain words: an old setting that no longer does anything.
    double OptimaReadmitSeed = 0.;

    /// pa_IpmStallWindow: the native IPM loop is accepted as converged when, over this many
    /// iterations, the energy and composition stay flat and pm.PCI only bounces. 0 = off.
    /// Default 30.
    /// In plain words: stops the original solver once the answer has clearly stopped changing.
    short IpmStallWindow = 30;

    /// pa_MbReproject: repairs an answer that fails the per-element mass-balance test by
    /// projecting it back onto A.Y = b; kept only if every species stays non-negative and the
    /// worst residual falls. Native path only. 0 = off, 1 = on (default).
    /// In plain words: a final touch-up that makes the element totals add up exactly.
    short MbReproject = 1;

    /// pa_DeterminacyWarn: warn when the energy fixes a present phase's amount only to worse than
    /// this relative uncertainty (native path, converged calls). 0 or negative = off. Default 1e-2.
    /// In plain words: warns when tiny amounts of a solid in the answer are not reliable.
    double DeterminacyWarn = 1e-2;

    /// pa_ColdRetryNudges: a failed cold native call (NEED_GEM_AIA) is re-solved cold at up to
    /// this many bulk compositions nudged by about 1e-15; the first that works seeds a warm solve
    /// at the exact composition. 0 = off. Default 4.
    /// In plain words: if the original solver fails, retry with a microscopically changed
    /// composition and then return to the exact one.
    long int ColdRetryNudges = 4;

    /// pa_OptimaPreSolveFirstIters: iteration limit for each pass of the reduced pre-solve's first
    /// attempt (never above the normal limit); a pass that hits it is discarded and the full
    /// budget is used instead. 0 = off. Default 6000.
    /// In plain words: gives up early on a preliminary step that is going nowhere.
    long int OptimaPreSolveFirstIters = 6000;

    /// pa_LpDualFillout: reserved, no effect. Kept so that the field order of BASE_PARAM
    /// (used by positional project-file serialisation) does not change.
    /// In plain words: an old setting that no longer does anything.
    long int LpDualFillout = 0;

    /// pa_FilloutBudget: limits how far the class fill-out may change an element's mass balance, as a
    /// fraction of that element's bulk amount. Native cold path only. 0 = off (default). The
    /// environment variable GEMS3K_FILLOUT_BUDGET overrides it.
    /// In plain words: limits how much the starting guess may disturb the element totals.
    double FilloutBudget = 0.;
    /// pa_StabTPD: tangent-plane (TPD) stability scan of absent mixed phases, reported on the
    /// trace's CERT line only. 0 = off, 1 = on (default; runs only when GEMS3K_NATIVE_TRACE_FILE
    /// is set and cannot change the answer), 2 = reserved.
    /// In plain words: a report-only check of whether a mixed phase left out of the answer
    /// should in fact form.
    long int StabTPD = 1;
    /// pa_IpmAugmentedKKT: how the native IPM loop solves its linear system for the dual u.
    ///   0  normal equations.
    ///   1  augmented (saddle-point) system, dense LU; slow on large systems.
    ///   2  (default) regularised least squares by QR.
    /// In modes 1 and 2 a small regularisation is added for an element no active species carries.
    /// If the augmented solve is singular, that step uses the normal equations.
    /// In plain words: a more accurate way to solve each step's equations in the original solver.
    long int IpmAugmentedKKT = 2;
    /// pa_IpmLoopTweaks: reserved, no effect. Kept so that the field order of BASE_PARAM
    /// (used by positional project-file serialisation) does not change.
    /// In plain words: an old setting that no longer does anything.
    long int IpmLoopTweaks = 0;
    /// pa_OptimaLineSearch: Optima's merit line search on the unmasked error, with this
    /// trigger factor (a step whose error exceeds factor x the previous one is line-searched).
    /// 0 = off. Default 1.5. A negative value uses the line search only on a re-run of a failed
    /// call, at |value|. Applies to the reduced pre-solve and the full solve. Requires the
    /// Optima fork's ErrorControl::execute (error updated at the new point before comparing).
    /// In plain words: when a step makes things much worse, try a shorter one.
    double OptimaLineSearch = 1.5;
    /// pa_OptimaFDDiagFloor: reserved, no effect. Kept so that the field order of BASE_PARAM
    /// (used by positional project-file serialisation) does not change.
    /// In plain words: an old setting that no longer does anything.
    long int OptimaFDDiagFloor = 0;
    /// pa_OptimaLSStallEscape: with the line search on, after this many consecutive line
    /// searches that leave the error unchanged (relative change <= 1e-8), keep the full step
    /// once. 0 = off. Default 10. Needs an Optima build with OPTIMA_LINESEARCH_STALL_ESCAPE;
    /// inert otherwise.
    /// In plain words: breaks the Optima solver out when it gets stuck at the same error.
    long int OptimaLSStallEscape = 10;
    /// pa_OptimaLSWindow: reserved, no effect. Kept so that the field order of BASE_PARAM
    /// (used by positional project-file serialisation) does not change.
    /// In plain words: an old setting that no longer does anything.
    long int OptimaLSWindow = 0;
    /// pa_OptimaLSRejectWorse: with the Optima line search on, a line search that ends at or
    /// above the pre-step error is discarded and the full step kept. 0 = off (default),
    /// 1 = on; always off in TNode::GEM_run()'s legacy retry. Needs an Optima build with
    /// OPTIMA_LINESEARCH_REJECT_WORSE; inert otherwise.
    /// In plain words: for difficult projects where the line search keeps making tiny,
    /// useless steps.
    long int OptimaLSRejectWorse = 0;
    /// pa_OptimaTpdAccept: when Optima ends not converged but the mass balance holds, accept
    /// the state if every off-bound species is stationary, every at-bound species of a
    /// single-species, aqueous or gas phase has gradJ >= -kktTol, and every absent non-ideal
    /// condensed phase has a tangent-plane (TPD) search minimum >= -value. 0 = off. Default 1e-6.
    /// In plain words: accepts an answer that is not fully converged when a direct check shows
    /// that no missing phase should form.
#ifndef GEMS3K_DEFAULT_OPTIMA_TPDACCEPT
#define GEMS3K_DEFAULT_OPTIMA_TPDACCEPT 1e-6
#endif
    double OptimaTpdAccept = GEMS3K_DEFAULT_OPTIMA_TPDACCEPT;
    /// pa_OptimaCgSeed: cold Optima seed by column generation - species Gibbs LP plus
    /// pseudo-compound columns priced by a tangent-plane (TPD) search; the value is the TPD
    /// tolerance. 0 = off (feasibility-LP seed). Default 1e-6.
    /// In plain words: a better first guess for the Optima solver, which also considers
    /// mixed (solution) phases.
#ifndef GEMS3K_DEFAULT_OPTIMA_CGSEED
#define GEMS3K_DEFAULT_OPTIMA_CGSEED 1e-6
#endif
    double OptimaCgSeed = GEMS3K_DEFAULT_OPTIMA_CGSEED;
    /// pa_OptimaColdRetry: a warm Optima call (SOP, SHP) that is not OK is re-solved cold
    /// (AOP) by TNode::GEM_run_optima_cold_retry(). 0 = off, 1 = retry after the full warm
    /// budget, 2 = fail fast: skip the warm call's own full-budget re-solves, then retry (default).
    /// In plain words: if starting from the previous answer fails, start again from scratch.
#ifndef GEMS3K_DEFAULT_OPTIMA_COLDRETRY
#define GEMS3K_DEFAULT_OPTIMA_COLDRETRY 2
#endif
    long int OptimaColdRetry = GEMS3K_DEFAULT_OPTIMA_COLDRETRY;
    /// pa_OptimaFinish: when an Optima call ends not converged, run a Newton finish on the
    /// fixed phase set (species amounts and multipliers, equality-constrained, line-searched
    /// on G) from Optima's last point (TMultiBase::PotentialSpaceFinish()); its result is
    /// judged by the same KKT / mass-balance / TPD checks. 0 = off, 1 = on (default).
    /// In plain words: if Optima stops just short, finish the job with the phases it found.
#ifndef GEMS3K_DEFAULT_OPTIMA_FINISH
#define GEMS3K_DEFAULT_OPTIMA_FINISH 1
#endif
    long int OptimaFinish = GEMS3K_DEFAULT_OPTIMA_FINISH;
    /// pa_OptimaAcceptRepair: when pa_OptimaTpdAccept's acceptance passes every test except
    /// the per-IC relative mass balance (mb_rel > 1), apply MassBalanceReproject() and re-test;
    /// the state is restored if mb_rel stays > 1 (DECIDE tpdaccept-repair). Needs
    /// pa_OptimaTpdAccept > 0. 0 = off (default), 1 = on.
    /// In plain words: lets a nearly finished Optima answer be accepted after a small fix of
    /// the element totals.
    long int OptimaAcceptRepair = 0;
    /// pa_NativeOptimaFinish: after a native solve (AIA/SIA) where a trace solid phase still tends to grow,
    /// finish with Optima from native's answer and keep the result only if it is better; also when native
    /// stops at its iteration limit or its linear system breaks down. The value is the most Optima
    /// iterations allowed. 0 = off, default 500.
    long int NativeOptimaFinish = 500;
    /// pa_OptimaCollapseCheck: when checking an Optima answer, judge a phase that has vanished by its pure
    /// species, so a phase that should come back is not taken as correctly absent. 0 = off, 1 = on (default).
    long int OptimaCollapseCheck = 1;
    /// pa_OptimaFastFirst: a cold Optima call (AOP) first tries a quick solve that takes the full step and
    /// clips each amount to its bounds; if that does not end OK, the normal solve runs from the same inputs.
    /// The value is the most Optima iterations the quick try may use. 0 = off, default 60.
    long int OptimaFastFirst = 60;
    /// pa_OptimaFloorCurvature: for a species in a solution phase, use the exact curvature of its ln(amount)
    /// term in the finite-difference Hessian, which the difference step under-states near the floor.
    /// 0 = off, 1 = on, 2 = also a difference step sized to the phase, the ideal-mixing block over all
    /// members of a solid or gas solution, and the barrier curvature of pure phases.
    /// 3 (default) = 1, plus the full mixing block on a nearly vanished gas or fluid phase, or on a vanishing copy
    /// of a twin pair whose other copy is present, so it leaves in a few steps.
    long int OptimaFloorCurvature = 3;
    /// pa_OptimaEmptySolvent: a cold start whose water is at or below its minimum amount gets water added,
    /// as when the solutes outweigh it. 0 = off, 1 = on (default).
    long int OptimaEmptySolvent = 1;
    /// pa_OptimaTpdVeto: the TPD acceptance may not overrule the phase stability check for a solution phase
    /// it cannot test itself, such as a vanished aqueous phase. 0 = off, 1 = on (default).
    long int OptimaTpdVeto = 1;
    /// pa_OptimaEmptyPreSolveSkip: once the reduced pre-solve (pa_OptimaDimReduce) finds no usable species set
    /// in a call, the retries of that same call skip it. 0 = off, 1 = on (default).
    long int OptimaEmptyPreSolveSkip = 1;
    /// pa_OptimaTwinReadmit: a solution phase that the twin-extinction step switched off is tried again after
    /// the solve when the phase stability test says it should form. 0 = off, 1 = on (default).
    long int OptimaTwinReadmit = 1;
    /// pa_OptimaDryRetry: a solve that leaves almost no liquid water, while the bulk could form much more,
    /// is stopped and solved once more from a start that puts hydrogen and oxygen into water first; 0 = off (default), 1 = on
    long int OptimaDryRetry = 0;
    /// pa_OptimaCallCap: most Optima iterations one AOP/SOP call may use before it is retried at nudged
    /// compositions; 0 = no cap, -1 = by system size (optima_call_cap()), N > 0 = N. Default 0.
    long int OptimaCallCap = 0;
    /// pa_OptimaCapNudges: how many 1e-15 nudges of the bulk composition a call that hit its cap is retried at;
    /// if none solves, the call is run once more without the cap. Default 3.
    long int OptimaCapNudges = 3;
    /// pa_OptimaFinishFirst: when a stalled Optima solve has failed, first try the Newton finish on its end point
    /// and keep it if no phase then wants to form; the full re-solves run only otherwise. 0 = off, 1 = on (default).
    long int OptimaFinishFirst = 1;
    /// pa_OptimaFastFrozen: the quick first try of a cold call stops after this many steps in a row that moved
    /// the state and the error only at rounding level. 0 = off, default 3.
    long int OptimaFastFrozen = 3;
    /// pa_OptimaFastDiverge: the quick first try stops 10 steps after its energy jumped above 10 times its start,
    /// if its lowest error since then is still above this many times the lowest before. 0 = off, default 10.
    double OptimaFastDiverge = 10.;
    /// pa_OptimaRedoxFree: before an Optima solve, check whether nothing in the bulk fixes the redox level (every
    /// element in one valence and no spare oxidant or reductant). 0 = off, 1 = log it (default), 2 = also aim the
    /// solve at pe 4 when the caller set no pH or Eh target.
    long int OptimaRedoxFree = 1;
    /// pa_OptimaFastRises: the quick first try of a cold call stops once its energy has gone up by more than 1 %
    /// on this many steps; a try that keeps climbing back rarely finishes. 0 = off, default 10.
    long int OptimaFastRises = 10;
    /// pa_OptimaPreSolveSeed: the first pass of the reduced pre-solve starts from the linear-Gibbs solution and its
    /// element potentials, each amount capped at what its elements allow; a pass that misses 300 steps drops the
    /// pre-solve for that call. 0 = off, 1 = on (default).
    long int OptimaPreSolveSeed = 1;
    /// pa_OptimaColdLevel: where no pre-solve runs, a cold Optima call starts from the linear-Gibbs solution and its
    /// element potentials, and two copies of one solid solution start at different compositions. 0 = off, 1 = on (default).
    long int OptimaColdLevel = 1;
    /// pa_OptimaLevelStart: with only condensed mixed phases, a cold Optima call first levels the system with mixed
    /// phases priced by composition search. If no mixed phase forms, Optima only confirms that answer; otherwise each
    /// mixed phase starts at the composition found. 0 = off, 1 = on (default).
    long int OptimaLevelStart = 1;
    /// pa_OptimaFastCollapse: the quick first try of a cold call stops at once when its energy has grown to more
    /// than this many times its starting value. 0 = off, default 1000.
    long int OptimaFastCollapse = 1000;
    /// pa_OptimaFrozenCut: the full Optima solve stops and hands over to the finish and retries once its energy
    /// and error have not changed at all for this many steps. Needs pa_OptimaStallWindow on; not used when Optima
    /// starts from a native answer (HOP, SHP). 0 = off, default 20.
    long int OptimaFrozenCut = 20;
    /// pa_OptimaDryOut: when the water runs out during an Optima solve, the water phase is taken out and the
    /// solids are solved alone; if that fails, the ordinary solve goes on with the water back. Off when an
    /// element can only be held by the water. Needs pa_OptimaStallWindow on. For systems where water is
    /// used up by design (e.g. cement hydration); it can return a wrong dry answer where the water only
    /// dips during the solve. 0 = off (default), 1 = on.
    long int OptimaDryOut = 0;
    /// pa_OptimaTailStop: an Optima solve whose energy has not changed for this many steps, with all the error
    /// left on trace species, stops and its point is checked for stability and taken as the answer.
    /// Needs pa_OptimaStallWindow on. 0 = off, default 100.
    long int OptimaTailStop = 100;
    /// pa_OptimaTailStopWarm: the same, on a warm call. 0 = off, default 10.
    long int OptimaTailStopWarm = 10;
    /// pa_OptimaTailDryOnly: where the tail stop acts. 0 = in every solve, 1 = only after pa_OptimaDryOut took
    /// the water phase out, 2 (default) = there, and in a warm call whose left-over error sits only on members of
    /// nearly vanished solid solutions.
    long int OptimaTailDryOnly = 2;
    /// pa_OptimaWaterKeep: 1 = a converged answer is kept when the retry with water put back fails;
    /// 2 = also, a dry answer that passed the tail-stop check is not retried with water (default); 0 = off.
    long int OptimaWaterKeep = 2;
    /// pa_OptimaSnapLeaving: in a warm Optima solve, a species that has limited the step this many times in a
    /// row on its way to zero (its target well below its floor) is put on its floor, so it stops creeping down
    /// a little each step. 0 = off, default 2.
    long int OptimaSnapLeaving = 2;
    /// pa_OptimaPartialPivLU: Optima's linear solve uses the faster partial-pivot LU for a matrix of at least this
    /// size, and redoes a step with today's full-pivot LU only when an equation would be dropped. 0 = off, default 100.
    long int OptimaPartialPivLU = 100;
    /// pa_IpmQrSortPivot: with pa_IpmAugmentedKKT = 2, a main-loop step with at least this many active species sorts the rows
    /// of its least-squares system by size and factorises it with column pivoting (faster and more accurate on large systems).
    /// Smaller steps keep the plain loop. 0 = off, default 100.
    long int IpmQrSortPivot = 100;
    /// pa_OptimaLSWarmTrials: on a warm Optima call, the line search tries at most this many step lengths
    /// (cold calls keep Optima's own limit of 20). 0 = off, default 5.
    long int OptimaLSWarmTrials = 5;
    /// pa_TpdStrongBasins: the thorough check of whether a mixed phase can form tries every start composition once
    /// and follows only the corner, ideal and uniform starts and this many best separate ones. 0 = follow every start, default 4.
    long int TpdStrongBasins = 4;
    /// pa_IpmLevelStart: a cold native call on a system without aqueous, gas or sorption mixtures starts from the
    /// leveled solution (the linear solve with mixed phases at their best compositions); 2 = also leave out the mixed
    /// phases it shows cannot form. A call that then fails is solved again from the ordinary start. 0 = off, default 2.
    long int IpmLevelStart = 2;
    /// pa_OptimaFastAcceptCut: a quick first try of a cold Optima call that used up its step cap but ended OK is kept when
    /// no missing pure phase and no missing solid or liquid mixture would form (thorough check); else the full solve runs.
    /// 0 = off, 1 = on (default).
    long int OptimaFastAcceptCut = 1;
    /// pa_OptimaLevelFinish: a cold Optima call whose leveled start holds mixed phases first runs the Newton finish on that
    /// start, and Optima starts from the finished point when it converges. 0 = off, 1 = on (default).
    long int OptimaLevelFinish = 1;
    /// pa_IpmNoProgress: native solver: when the energy has not improved over this many iterations, stop and let Optima take
    /// over from that state (AIA/SIA: the Optima finish, then a cold re-solve without the stop if it fails; HOP/SHP: the
    /// Optima leg, warm, then cold). 0 = off; default 500.
    long int IpmNoProgress = 500;
    /// pa_OptimaFDPureSkip: Optima solver: the Hessian column of a pure-phase species is set to zero without an activity
    /// call (its amount changes no potential). 0 = off, 1 = on (default).
    long int OptimaFDPureSkip = 1;
    /// pa_OptimaFinishStall: the Newton finish gives up when, over this many steps, its largest reduced gradient has not
    /// halved, G has not fallen and its step size has stayed within 1 % (it then reports no convergence, as at its step
    /// limit). 0 = off. Default 15.
    long int OptimaFinishStall = 15;
    /// pa_OptimaTraceRespread: Optima solver: when the species floors of an element add up to more than its amount, its
    /// true amount is spread over its solution species from the converged potentials. 0 = off, 1 = on (default).
    long int OptimaTraceRespread = 1;
    /// pa_FailRescue: a call that ends with an error or a bad result is tried again in steps - warm from the node's previous
    /// answer, nudged compositions, a walk from the previous answer, looser settings, the other solver - each finished at the
    /// exact composition with the normal settings; with 2, if none passes, the best balanced attempt is returned as a bad
    /// result (GEM_rescue_step() = 6). 0 = off, 1 = steps 1-5, 2 = steps 1-6 (default).
    long int FailRescue = 2;
    /// pa_RescueWalkSteps: number of equal steps of the walk in pa_FailRescue; 0 = no walk. Default 8.
    long int RescueWalkSteps = 8;
    /// pa_RescueBudget: after the whole call has used this many iterations of both solvers, pa_FailRescue starts no new
    /// attempt of steps 1-4 (an attempt that passed still gets its finish); step 5 (the other solver) still runs once.
    /// 0 = no limit. Default 100000.
    long int RescueBudget = 100000;
    /// pa_OptimaStallNetFactor: Optima's last re-solve after a stall gets at most this many times the iterations
    /// of the solve that stalled (at least 2000, at most the normal limit). 0 = the normal limit. Default 5.
    long int OptimaStallNetFactor = 5;

    void write(GemDataStream& oss);
    void read(GemDataStream& iss);
};

/// What the quick first try of a cold Optima call has seen so far (pa_OptimaFastFrozen, pa_OptimaFastDiverge,
/// pa_OptimaFastRises, pa_OptimaFastCollapse);
/// reset at the start of each quick try.
struct OptimaQuickTryWatch
{
    int frozenRun = 0;          ///< steps in a row that moved only at rounding level
    double errPrev = -1.;       ///< error at the previous step; < 0 before the first
    long int step = 0;          ///< steps seen by the divergence check
    double f0 = 0.;             ///< energy at the first step
    double preMin = 0.;         ///< lowest error before the energy jump
    double postMin = 0.;        ///< lowest error since the jump
    long int jumpAt = -1;       ///< step of the jump; -1 = none yet
    long int rises = 0;         ///< steps whose energy went up by more than 1 %
    long int riseSteps = 0;     ///< steps seen by the rise count
    double fPrev = 0.;          ///< energy at the previous step
    double cf0 = 0.;            ///< energy at the first step, for the collapse check
    long int cstep = 0;         ///< steps seen by the collapse check
};
extern thread_local OptimaQuickTryWatch g_optimaQuickTry;   // defined in node.cpp

// ---------------------------------------------------------------------------
// Effective values of auto-gated settings
// ---------------------------------------------------------------------------
// Some pa_* fields are three-valued: a configured 0 means "decide from the problem",
// so the value in the project file and the trace's SET line is not the value the solver
// ran at. Each such field is resolved by one function below, called by both the solver
// and native_trace_run_header()'s EFF line, so the two cannot disagree. Any new
// auto-gated field follows the same pattern.

/// Marker that MultiConstInit() seeds pm.FX with (ipm_simplex.cpp): "total Gibbs energy
/// not computed yet". Deliberately absurd (+7.777777e6) so it cannot pass for a result.
/// packDataBr() checks for it, so an answer that was never priced is reported.
constexpr double kTotalGibbsEnergyUnset = 7777777.;

/// Species-count gate at or above which pa_OptimaDimReduce = 0 (AUTO) turns the
/// dimension-reduction pre-solve on.
constexpr long int kOptimaDimReduceAutoMinDC  = 200;
/// Pass count AUTO selects when the gate opens.
constexpr long int kOptimaDimReduceAutoPasses = 8;

/// Resolves pa_OptimaDimReduce's three-valued setting to the pass count the solver
/// attempts: > 0 explicit count, < 0 off, 0 AUTO. `nDC` is the species count (pm.L).
/// Whether the pre-solve is reached also depends on the call: it runs on a cold Optima
/// leg (pm.pNP == 0), and on a HOP leg only under an explicit positive setting.
/// Iteration cap of one Optima call (pa_OptimaCallCap): -1 = max(1000, 20 per species), 0 = none.
/// Used by both the solver and the EFF trace line.
inline long int optima_call_cap( long int configured, long int nDC )
{
    if( configured < 0 )
        return ( 20L * nDC > 1000L ) ? 20L * nDC : 1000L;
    return configured;
}

inline long int optima_dimreduce_passes( long int configured, long int nDC )
{
    if( configured == 0 && nDC >= kOptimaDimReduceAutoMinDC )
        return kOptimaDimReduceAutoPasses;
    return configured;
}

/// Per-pass Optima iteration budget of OptimaReducedPreSolve(): max(2000, pa_IIM),
/// lowered to pa_OptimaPreSolveFirstIters on the first attempt when that is > 0. Never
/// raises the budget. Used by both the solver and the EFF trace line.
inline long int optima_presolve_pass_budget( long int configured, long int iim, bool firstAttempt )
{
    const long int full = ( iim > 2000L ) ? iim : 2000L;
    return ( firstAttempt && configured > 0 && configured < full ) ? configured : full;
}


/// Cap that pa_OptimaEarlyStabilityAt = 0 (AUTO) resolves to on a cold leg when the
/// system carries a multisite (sublattice) solid-solution model.
constexpr long int kOptimaEarlyStabilityAutoCap = 200;

/// Cap that AUTO resolves to on a warm Optima leg (pm.pNP != 0). There the dual is
/// already consistent, so the cap always fires and a probe that finds nothing costs
/// exactly N; a short cap keeps that loss small.
constexpr long int kOptimaEarlyStabilityAutoWarmCap = 25;

/// Default re-solve strategy of the early-probe safety net, resolved in one place for
/// the solver and the trace. GEMS3K_OPTIMA_NET_RESUME overrides it.
/// 0  re-solve from `initialState` (a probe that finds nothing costs exactly N extra).
/// 1  resume from the probe's own end state (can lose answers).
/// 2  resume; if the resumed re-solve does not converge, re-solve once from
///    `initialState` as 0 would.
/// 3  as 2, but the fallback is deferred to the end of the retry ladder, so it is paid
///    only when every other rescue has also failed.
constexpr int kOptimaNetResumeDefault = 3;

/// Effective resume mode for this process. Reads the environment override once.
inline int optima_net_resume_mode()
{
    static const int mode = []() -> int {
        const char* e = std::getenv( "GEMS3K_OPTIMA_NET_RESUME" );
        return ( e && *e ) ? std::atoi( e ) : kOptimaNetResumeDefault;
    }();
    return mode;
}

/// Is this phase's built-in mixing model a multisite (sublattice) solid solution?
/// True for SM_BERMAN, SM_CEF and SM_MBW only.
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

/// Resolves pa_OptimaEarlyStabilityAt's three-valued setting to the value the solver
/// arms: > 0 is an explicit cap at that iteration, < 0 is the trend form with |value|
/// consecutive falls, and 0 is AUTO - a cap at kOptimaEarlyStabilityAutoCap (cold leg)
/// or kOptimaEarlyStabilityAutoWarmCap (warm leg) when the system carries a multisite
/// solid solution, and off otherwise. Multisite models can produce two solution phases
/// with identical stoichiometry and G0 (interchangeable twins), which an interior-point
/// method approaches only asymptotically; stopping early lets the phase-selection and
/// extinction steps resolve that. There is no "off" encoding: to disable it on a
/// multisite system, set a cap larger than pa_IIM.
/// `warmLeg` is pm.pNP != 0. AUTO returns the short cap there; an explicit setting
/// applies to both legs.
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
/// TMultiBase::WorstPhaseStabilityViolation().
struct PhStabViolation
{
    long int k = -1;            ///< phase index
    double   viol = 0.;         ///< size past native's threshold, in log10 units
    bool     wasAbsent = false; ///< true: absent but stable; false: present but unstable
    /// The phase's logSI is saturated at StabilityIndexes()' overflow guard (ln_ax_dual
    /// clamped to +609 / -608), i.e. a guard value, not a measured driving force.
    bool     clamped = false;
};

/// What one phase-assemblage stability scan looked at. The clamped count matters: on an
/// unconverged state many phases saturate the overflow guard, and a ranking of guard
/// values carries no information.
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

    /// ScaleSystemToInternal()/RescaleSystemFromInternal() (ipm_simplex.cpp) must not rescale
    /// the "no upper limit" sentinel (>= 1e6) of pm.DUL[]/pm.DLL[]/pm.PUL[]. The scaling
    /// decision is recorded here: the pre-scale value and the value left after scaling. An
    /// entry still holding the value left is restored verbatim; one rewritten in between (by
    /// Set_DC_limits(true) on the warm path) is in internal units and takes the original test.
    /// Filled by every ScaleSystemToInternal(), consumed and cleared by the matching
    /// RescaleSystemFromInternal(); an empty or mismatched vector falls back to the original test.
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
    /// copyMULTI()'s body. realloc = true is copyMULTI() itself. realloc = false copies values
    /// into this object's existing arrays and allocates nothing, for restoring a snapshot into
    /// a live MULTI whose TSolMod/TSorpMod/TKinMet objects point into those arrays. Both
    /// objects must have identical dimensions. Used by TNode::GEM_trace_regimes().
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
    /// Total Gibbs energy G(X) of the converged system, in RT units (moles), computed the
    /// same way for every solver path so results are comparable across AIA/SIA and
    /// AOP/SOP/ROP. Reads the primal from pm.Y[] and, as a side effect, copies Y into X and
    /// refreshes pm.XF/pm.XFA (a no-op on a converged state - do not call mid-solve).
    /// Does not use GX(0.) or pm.FX, which are path-dependent: the excess term is rebuilt from
    /// G0[]+fDQF[]+F0[].
    /// In plain words: the total energy of the answer, the main measure of whether it is right.
    double TotalGibbsEnergy();

    // ------------------------------------------------------------------------------------
    // Report-only certificate instruments, printed on the CERT trace record. None of them
    // enters CERT's mb_pass. They are computed only when native_trace_file() is open, so
    // production calls do not pay for them. They rebuild Gj = G0+fDQF+F0 at the returned
    // pm.X[] rather than reading pm.G[] or pm.F[], which differ between solver paths at exit.

    /// Primal chemical potentials at the returned amounts pm.X[], into the caller's array. A
    /// read-only mirror of PrimalChemicalPotentials() that writes no pm.* state and rebuilds
    /// Gj = G0+fDQF+F0 instead of reading pm.G[]. F[j] is 0 for a species below pm.DcMinM or
    /// in a phase PrimalChemicalPotentials() would skip. Sized to pm.L.
    void CertPrimalPotentials( std::vector<double>& F ) const;

    /// kkt_max: worst sign-aware reduced-gradient residual in RT over the species present at
    /// the answer, with the dual pm.U[] of the last linear solve. s_j = F_j - sum_i U_i A(i,j);
    /// interior |s_j|, at a lower bound max(-s_j, 0), at an upper bound max(s_j, 0), and 0 for
    /// a species whose box is degenerate (DUL <= DLL, an equality constraint). Same sign logic
    /// as ipm_optima.cpp's KKT check, except that Optima's log-barrier term -tau/X_j is not
    /// subtracted (it is not part of the thermodynamic objective).
    /// \param worstJ index of the species carrying the maximum, -1 if none.
    /// \return the maximum, or -1. if nothing could be scored.
    /// In plain words: how far the answer is from a perfect minimum, as one number.
    double CertKktMax( const std::vector<double>& F, long int& worstJ ) const;

    /// dual_free_dirs: how many directions the dual is free along at the answer. The interior
    /// species (present and away from both bounds) fix u through s_j = 0; when their
    /// stoichiometry columns span rank r < N, the remaining N-r directions leave u
    /// undetermined, and every bound-active species' reduced gradient (hence kkt_max) depends
    /// on where along them the solver stopped. Uses the same modified Gram-Schmidt and 1e-8
    /// relative rank tolerance as the free-dual search in ipm_optima.cpp.
    /// \return N - rank, with rank and N returned in the out parameters.
    /// In plain words: tells whether some element potentials are not pinned down by the answer.
    long int CertDualFreeDirs( const std::vector<double>& F, long int& rank, long int& nIC, double* zzFree = nullptr ) const;

    /// The RANK trace record: the present species' stoichiometry and the two normal-equation
    /// matrices each solver stage forms from it, at the returned answer pm.X[].
    /// rank/sv_ratio/chg_res depend on the geometry alone; cond_ipm/cond_mbr include the
    /// weight each stage's own linear solve applies.
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

    /// Fills a CertRankReport at the returned pm.X[] (see CertRank() in ipm_main.cpp).
    /// Report-only: it changes no solver state. Leaves `r` at its defaults if it cannot run.
    /// In plain words: measures how well-posed the final linear equations are.
    void CertRank( CertRankReport& r ) const;

    /// curv_min: the smallest curvature, over every present mixed non-aqueous phase, of that
    /// phase's finite-difference curvature block at the answer. A negative value means the phase
    /// sits at a maximum along the unmixing direction. Report-only: it restores everything it
    /// touches and runs after packDataBr().
    /// \param worstPhase index of the phase carrying the minimum, -1 if none.
    /// \return the minimum, or +1e300 if no phase qualified.
    /// In plain words: checks that each mixed phase in the answer is truly stable and not
    /// sitting where it would rather split in two.
    double CertCurvMin( long int& worstPhase );

    /// pa_StabTPD = 1: tangent-plane stability scan of every absent or trace non-aqueous mixed
    /// phase at the returned answer. Sorption, ion-exchange and polyelectrolyte phases are
    /// skipped. Report-only: restores what it touches.
    /// \param worstPhase phase carrying the minimum, -1 if none.
    /// \param nScanned   absent phases scanned.
    /// \param nDisagree  phases the single-point index calls stable while the search finds them unstable.
    /// \return min over scanned phases of the smallest TPD(y) in RT (< 0: unstable), +1e300 if none.
    /// In plain words: checks whether a mixed phase left out of the answer should in fact form.
    double CertStabTPD( long int& worstPhase, long int& nScanned, long int& nDisagree );
    /// Prototype: composition search for one absent non-ideal phase; see ipm_main.cpp.
    double NativeTpdPhase( long int k, long int p0, std::vector<double>& ybest );
    /// The leveled solution: species Gibbs LP plus mixed phases at their best compositions,
    /// priced by a tangent-plane search (ipm_level.cpp). Used by the native start and by Optima.
    bool ColumnGenerationSeed( std::vector<double>& nOut, double tol );
    /// The strongest drive of a missing phase to form at the current answer (RT per mole, < 0: it would form): a
    /// missing pure phase by its reduced energy, a missing solid or liquid mixture by the thorough search
    /// (ipm_level.cpp). worstPhase is the phase, -1 if none.
    double AbsentPhaseWorstDrive( long int& worstPhase );

    double CalculateEquilibriumState( /*long int typeMin,*/ long int& NumIterFIA, long int& NumIterIPM );
    void InitalizeGEM_IPM_Data();
    virtual void DC_LoadThermodynamicData( TNode* aNa = nullptr );

#ifdef USE_OPTIMA_SOLVER
    // Optima solver entry points (USE_OPTIMA_SOLVER builds only).

    /// Request dn/db sensitivity derivatives from the Optima solver. Off by default: each
    /// costs an extra solve of the factorised system per bulk-composition column.
    /// Problem::c are mapped to the bulk composition (bec = I), so Sensitivity::xc is dn/db.
    /// In plain words: also computes how the answer would change if the bulk composition
    /// changed slightly.
    bool optima_want_sensitivity = false;

    /// dn/db from the last solve when optima_want_sensitivity was set, row-major
    /// [ (L+R) x N ]. Empty if not requested or the solve did not reach the sensitivity step.
    std::vector<double> optima_dndb;
    long int optima_dndb_rows = 0, optima_dndb_cols = 0;

#endif   // USE_OPTIMA_SOLVER
    /// One flag per phase: has PhaseSelect() already re-inserted this phase at its
    /// composition ceiling rather than at the fixed pa_DFYs during this solve? Sized pm.FI
    /// and cleared in GEM_IPM(). A phase gets one budget-sized insertion; if it is lost
    /// again it is skipped for the rest of the solve. Used on the native path, so it is
    /// declared outside the USE_OPTIMA_SOLVER block (closed and reopened here, which keeps
    /// the member's position in the class).
    std::vector<char> insBudgetTried;
    // Control conditions (pH, Eh targets) exist only with Optima; without it the list stays empty.
    /// Read-only access to the registered control conditions. After a
    /// CalculateEquilibriumStateOptima() call each entry's titrantAmount, achievedValue and
    /// targetMet hold that call's result.
    const std::vector<EqControlCondition>& GetControlConditions() const
    { return optima_control_conditions; }
    /// Puts back control conditions saved with GetControlConditions(), results included (used when a
    /// caller returns an earlier solve's answer).
    void RestoreControlConditions( const std::vector<EqControlCondition>& saved )
    { optima_control_conditions = saved; }
#ifdef USE_OPTIMA_SOLVER
    // Equilibrium via the Optima interior-point solver, as an alternative to the IPM/MBR loop
    // (ipm_optima.cpp). Dispatched from TNode::GEM_run() for NEED_GEM_AOP/SOP: pm.pNP == 0 is a
    // cold start, 1 a warm start from the existing pm.Y[]. Registered control conditions
    // (pH, Eh or custom) are solved in the same Newton system.
    // `referenceMode` (default false) is for NEED_GEM_ROP, a fixed reference setup: uniform
    // tiny seed, PartiallyExact Hessian, Optima's own default options, always cold.
    // In plain words: AOP is the Optima solver from scratch, SOP is the same starting from
    // the previous answer, and ROP is a plain reference setup used only for comparison.
    /// \param runKinetics run the kinetics/metastability step (RunKineticsStep()): true for
    ///        a normal AOP/SOP call, false for the HOP Optima leg, whose native leg already
    ///        ran it for this time step.
    double CalculateEquilibriumStateOptima( long int& NumIterFIA, long int& NumIterIPM, bool referenceMode = false,
                                            bool runKinetics = true );

    /// HOP: native selects the phases (its own IPM/MBR/PSSC, cold), then Optima finishes,
    /// warm-started from native's primal and dual (NEED_GEM_HOP). If native converges but
    /// the Optima leg fails, native's answer is restored and reported as BAD_GEM_HOP, so the
    /// result is never worse than a plain native solve.
    /// `warmNative` selects SHP (NEED_GEM_SHP): the native leg starts warm (pm.pNP = 1) from
    /// the state the node already holds, for sequential work such as sweeps and transport
    /// loops. If that warm native leg fails it falls back to a cold one, i.e. to HOP.
    /// In plain words: runs the original solver first and lets the Optima solver polish
    /// its answer, keeping the original answer if the polish fails.
    double CalculateEquilibriumStateHOP( long int& NumIterFIA, long int& NumIterIPM,
                                         bool warmNative = false );
    /// After a native solve that left a trace solid phase still tending to grow, finishes with Optima
    /// from native's answer (at most pa_NativeOptimaFinish Optima iterations) and keeps the result only
    /// if it is better; otherwise native's answer is put back unchanged. Adds the iterations used.
    double NativeOptimaFinish( long int& NumIterFIA, long int& NumIterIPM );
    /// After a native solve that stopped at its iteration limit (E06IPM): Optima solves from native's last state
    /// (then once from scratch), at most pa_NativeOptimaFinish iterations each. True if it found a clean answer.
    bool NativeOptimaRescue( long int& NumIterFIA, long int& NumIterIPM, double& seconds, bool broken = false );

    /// Registers (or replaces, by name) a pH control condition for the next
    /// CalculateEquilibriumStateOptima() call. Persists until ClearControlConditions().
    /// `tolerance < 0` (the default) uses pa_p->GAS (EqControlCondition::defaultToleranceFn()).
    /// In plain words: asks the Optima solver to reach a given pH by adding a titrant.
    void SetControlCondition_pH( double pH_target, double tolerance = -1. );
    /// Registers (or replaces, by name) an Eh control condition (V) for the next
    /// CalculateEquilibriumStateOptima() call. `tolerance < 0` (the default) uses pa_p->GAS.
    /// In plain words: asks the Optima solver to reach a given redox potential.
    void SetControlCondition_Eh( double Eh_target, double tolerance = -1. );
    /// Removes every registered control condition.
    void ClearControlConditions();


    /// Worst-case (infinity-norm) mass-balance residual |A*Y - B| over pm.Y[]/pm.B[]. A
    /// diagnostic for callers and tests; the solver's own check is CheckMassBalanceResiduals().
    double OptimaMaxMassBalanceResidual();

    /// Detects whether the aqueous solvent (pm.LO) dominates its own phase by mass in the
    /// trial amounts `x` (length pm.L). Returns false if it already dominates or there is no
    /// pm.LO; otherwise returns true and sets `waterSeedOut` to min(bulk-H/2, bulk-O) (or
    /// the phase's other-species total if H/O cannot be resolved by name), capped by
    /// `upperBound` if positive. Used by the retry logic of CalculateEquilibriumStateOptima().
    /// In plain words: detects a starting guess with almost no water and restores it.
    /// True when the system contains an aqueous phase. Tests the phase classifier, not pm.LO.
    bool HasAqueousPhase() const;

    /// Result of CheckRedoxFree(): `free` when, for some choice of one valence per element among those its
    /// species have, the bulk carries no spare electrons, so nothing but traces sets the redox level.
    struct RedoxFreeCheck {
        bool checked = false;     ///< the system has an aqueous phase and a charge element
        bool consistent = false;  ///< every element in the bulk has only one valence
        bool determined = false;  ///< every element's valences are known and the search finished
        bool free = false;
        double excess = 0.;       ///< electron excess of the bulk for the closest choice, mol
        double rel = 0.;          ///< excess over the bulk's total valence amount
    };
    /// Checks whether nothing in the bulk fixes the redox level (pa_OptimaRedoxFree).
    RedoxFreeCheck CheckRedoxFree() const;

    /// Per phase: 1 if another phase has the same end-members (a miscibility gap needs two copies).
    std::vector<char> PhaseHessianTwins() const;
    /// Per-run record of the curvature floor (twin flags, times it fired); see ipm_optima.cpp.
    void PhaseHessianBegin();
    void PhaseHessianEnd();
    /// Warns about copies of a solid solution that ended with the same composition (final answer only).
    void WarnCollapsedTwins();
    std::vector<char> phfTwin;
    std::vector<long int> phfFired;

    bool DetectSolventCollapseAndReseed( const double* x, double upperBound, double& waterSeedOut );

    /// Like DetectSolventCollapseAndReseed(), for every mixed phase (not `excludePhaseIdx`, default
    /// -1 = none): finds phases in the trial amounts `x` that are almost empty and appends
    /// (species index, reseed value) pairs to `reseedsOut`. Returns true if any was appended.
    /// In plain words: rescues a phase that the starting guess has almost emptied.
    /// `stabilityGated` (cold seed only): refill a phase, gas and fluid included, only if the
    /// linearised-Gibbs LP says it should form, split by how favourable each end-member is.
    bool DetectPhaseCollapseAndReseed( const double* x,
                                        std::vector<std::pair<long int,double>>& reseedsOut,
                                        long int excludePhaseIdx = -1,
                                        bool stabilityGated = false );

    /// LP-feasibility seed: a starting point that satisfies pm.A*n = pm.B, n >= 0, found with a
    /// small dense simplex. Returns false (nOut unspecified) if the LP fails or its result does
    /// not verify; the caller then uses a simpler seed.
    /// In plain words: finds a starting point that satisfies the element totals exactly.
    bool LPFeasibilitySeed( std::vector<double>& nOut );
    bool PotentialSpaceFinish( double dcFloor, const std::vector<double>& xlower, const std::vector<double>& xupper );

    /// Dual of the linearised-Gibbs LP (min sum_j G0[j]*n_j s.t. A n = b, n >= 0). Used by
    /// OptimaReducedPreSolve() to choose the initial species set. Returns false (yOut untouched)
    /// if the LP fails or its dual does not verify.
    /// In plain words: a quick rough estimate of which species are likely to matter.
    bool LPGibbsDual( std::vector<double>& yOut, const double* cost = nullptr );

    /// Phase-assemblage stability scan for the Optima path: refreshes pm.YF/pm.YFA from
    /// pm.Y, calls StabilityIndexes(), and reports the worst disagreement between "is this
    /// phase present" and "does its stability index say it should be", using PhaseSelect()'s
    /// thresholds (pa_p->DF / pa_p->DFM). Shared by the phase-selection retry and the final
    /// trustworthiness check, so both use one definition of a violation. Kinetically
    /// restricted and caller-deactivated phases are exempt.
    ///
    /// \param presenceThreshold  phase total at/below which a phase counts as absent.
    /// \param dcFloor  the species lower-bound floor of the Optima boxes. A phase whose
    ///        species all sit at their lower bound counts as absent whatever its total.
    /// \param exemptSpecies  optional, size >= pm.L when non-null: species the caller has
    ///        fixed; their phase is skipped (the interchangeable-twin case).
    /// \param violOut  size of the worst violation, in logSI units past the threshold; 0 when none.
    /// \param wasAbsentOut  true if the worst violation is "absent but stable", false if
    ///        "present but unstable".
    /// \param rankedOut  optional: every violation, ordered by size (descending, ties by
    ///        phase index), so a caller that may not act on the worst one can walk the list.
    ///        The ordering is meaningful only for entries that are not `clamped`.
    /// \param censusOut  optional: what the scan looked at, including how many phases had a
    ///        logSI saturated at StabilityIndexes()' overflow guard.
    /// \return index of the worst violating phase, or -1 if the assemblage is self-consistent.
    /// In plain words: checks whether any phase is missing that should be there, or present
    /// that should not be.
    /// \return index of the worst violating phase, or -1 if the assemblage
    ///         is self-consistent.
    long int WorstPhaseStabilityViolation( double presenceThreshold, double dcFloor,
                                           const char* exemptSpecies,
                                           double& violOut, bool& wasAbsentOut,
                                           std::vector<PhStabViolation>* rankedOut = nullptr,
                                           PhStabCensus* censusOut = nullptr,
                                           const char* probeStage = "check" );

    /// Species-level dimension-reduction pre-solve (pa_OptimaDimReduce). Solves the
    /// equilibrium over a reduced set of species, growing the set by pricing the omitted
    /// columns on the reduced dual, and leaves its answer in pm.Y[] (primal) and pm.U[]
    /// (dual) for the full-dimension solve to warm-start from and verify.
    /// \param maxPasses  readmission passes allowed (pa_OptimaDimReduce).
    /// \param dcFloor    the same species floor the full path uses.
    /// \param dimTol     the initial-set rule, with pa_OptimaDimReduceTol's sign convention
    ///        (> 0 an absolute RT threshold, < 0 the rank rule -|dimTol| x N, 0 the seed's own
    ///        support only). Passed in so the caller can retry a discarded pre-solve under the
    ///        other rule.
    /// \param passBudget  Optima iteration budget of each pass (optima_presolve_pass_budget()).
    /// \param iterationsOut  Optima iterations spent here, for pm.ITG.
    /// \param activeOut  size of the final active set, for logging.
    /// \return true if pm.Y[]/pm.U[] now carry a state worth warm-starting from; false if
    ///         the pre-solve was skipped or discarded (nothing for the caller to undo).
    /// In plain words: first solves a smaller problem with only the likely species, then
    /// hands that answer to the full solve as a starting point.
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

    /// Control conditions registered via SetControlCondition_pH()/_Eh() (or appended
    /// directly as a custom EqControlCondition - see ipm_optima.h), active for the next
    /// CalculateEquilibriumStateOptima() call. Entries are fully resolved when added; the
    /// closed-form functions capture resolved indices by value and read G0[]/T at call time.
    std::vector<EqControlCondition> optima_control_conditions;

    /// When true, SmoothingFactor() returns exactly 1.0 regardless of pm.FitVar[3]/[4],
    /// disabling the IPM-2 chemical-potential smoothing for one
    /// CalculateEquilibriumStateOptima() call. That smoothing blends with pm.F0[] from the
    /// previous evaluation, so with a factor below 1 the objective Optima minimises would
    /// depend on the history of evaluated points; with 1.0 the blend is an exact no-op.
    /// Set and cleared only by CalculateEquilibriumStateOptima(); always false on native paths.
    bool optima_disable_smoothing = false;

    /// True only while CalculateEquilibriumStateHOP()'s Optima leg runs on top of a
    /// successful native solve. Read only by the dimension-reduction gate in
    /// CalculateEquilibriumStateOptima(), which is otherwise cold-start-only (pm.pNP == 0):
    /// with it, the pre-solve may also run on the HOP leg, re-expressing native's assemblage
    /// in Optima's terms (native's absent species arrive at Optima's floor with reduced
    /// gradients Optima's own test may reject). Takes effect only for pa_OptimaDimReduce > 0;
    /// AUTO (0) does not reach it.
    bool optima_hop_leg = false;
    /// Scale factor of this node's previous solve (native or Optima); 0 before the first one.
    double prev_scfact = 0.;

    /// Per-leg cost of the last two-leg (HOP/SHP) solve. NumIterFIA/NumIterIPM report the
    /// sum of both legs; this keeps the legs apart. Recorded only - nothing in the solver
    /// reads it. `valid` is false unless the last GEM_run() dispatched a two-leg mode
    /// (TNode::GEM_run() clears it before every dispatch). Set from a destructor, so a call
    /// that throws is recorded too, with `failed` true.
    /// In plain words: how much of a combined solve's time went to each of its two solvers.
    struct HopLegSplit
    {
        bool valid = false;        ///< the last solve was HOP or SHP
        /// The call left by exception. The iteration counts are still real work spent, but
        /// timeOptima is not recoverable and is left at 0; a caller summing splits must read this.
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
    /// pa_IpmAugmentedKKT: zero-row rescues in the current InteriorPointsMethod() call (reset
    /// at its entry). Only the first one per call emits DECIDE "ipmkkt-zerorow".
    long int ipmKktRescues = 0;
    /// pa_IpmAugmentedKKT = 1 or 2: the main-loop (initAppr = false) solve for pm.U without
    /// forming A^T W A. \return 0 solved, 2 singular - the caller then takes that step by the
    /// normal equations (DECIDE ipmkkt-fallback).
    /// In plain words: an alternative, more robust way to solve the linear equations of each step.
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

    /// One kinetics/metastability time step (TKinMet), for every solver path. Must run
    /// before ExcludeRedundantDCs() and ScaleSystemToInternal(); HOP's Optima leg does not
    /// run it.
    /// In plain words: applies the time-dependent (kinetic) limits before the solve.
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

    /// pa_MbReproject: projects the amounts back onto A.X = b with one N x N solve over a
    /// rank-revealing set of the most abundant species. Native path; called only when the
    /// answer has failed its own per-IC mass-balance test. Returns true only if the
    /// correction kept every species non-negative and strictly reduced the worst relative
    /// residual; otherwise restores the amounts and returns false.
    /// `amt` is the vector to repair - pm.X on the final answer, pm.Y after PSSC. Both are
    /// re-synchronised on success. With keepPartial a repair that improves but does not fully
    /// pass is kept (for amounts that only seed further iterations).
    /// In plain words: a final touch-up that makes the element totals add up exactly.
    bool MassBalanceReproject( double* amt, bool keepPartial = false );
    /// Warns when a present phase's amount is not determined by the minimised energy: G is
    /// flat enough along a mass-balance-preserving direction that answers differing in that
    /// phase's amount cannot be told apart at the solver's energy resolution. Read-only.
    /// In plain words: warns when the answer's split between phases is not well defined.
    void EnergyDeterminacyCheck();
    /// A species held at zero by ExcludeRedundantDCs() for the current call, with the
    /// metastability settings RestoreRedundantDCs() puts back.
    struct RedundantDCHold { long int j; char rlc; double dll, dul; };
    /// Finds redundant species - identical stoichiometry, class and standard properties at
    /// the current T,P, twice in one phase or as two single-species phases - warns, and holds
    /// every copy after the first at zero for this call (internal DLL = DUL = 0,
    /// RLC = BOTH_LIM; any starting amount moved onto the kept one). The caller's DATABR
    /// dll/dul are not touched. RestoreRedundantDCs() undoes it.
    /// In plain words: removes duplicate species for the duration of a calculation.
    std::vector<RedundantDCHold> ExcludeRedundantDCs();
    void RestoreRedundantDCs( const std::vector<RedundantDCHold>& held );
    /// Warns when an element can exist in only one multi-component phase and that phase is
    /// a trace amount made up largely of the element. Read-only.
    /// In plain words: flags a fragile system definition, where a phase exists only to host
    /// one element.
    void StrandedElementCheck();
    /// Warns when a solid phase is present only in a trace amount but still tends to grow:
    /// the answer may stop short of equilibrium. Read-only.
    void TraceGrowthCheck();
    /// Set by TraceGrowthCheck() when a trace phase still tends to grow, with its warning text;
    /// cleared at the start of each native solve.
    bool trace_growth_stuck = false;
    std::string trace_growth_note;
    long int InteriorPointsMethod( long int &status/*, long int rLoop*/ );
    void AutoInitialApproximation( );

    // ipm_main.cpp - miscellaneous fuctions of GEM IPM-2
    void MassBalanceResiduals( long int N, long int L, double *A, double *Y,
                               double *B, double *C );
    double OptimizeStepSize( double LM );
    void DC_ZeroOff( long int jStart, long int jEnd, long int k=-1L );
    void DC_RaiseZeroedOff( long int jStart, long int jEnd, long int k=-1L );
    /// pa_FilloutBudget: scales the class fill-out so it changes each element's mass balance
    /// by at most that fraction of the element's own bulk amount. Native cold path.
    /// In plain words: limits how much the starting guess may disturb the element totals.
    void ApplyFilloutBudget( const std::vector<double>& yLp );
    /// Effective pa_FilloutBudget: the field, unless GEMS3K_FILLOUT_BUDGET overrides it (kept for
    /// the fillout.budget test). Both the call site and the mechanism read this.
    double FilloutBudgetValue() const;
    /// Largest amount of a single-species phase the bulk composition can supply,
    /// min_i b_i/a(j,i) over the ordinary IC rows. Caps PSSC's pure-phase insertion amount
    /// (pa_DFYs) so an insertion is never infeasible by construction.
    /// In plain words: the most of a mineral that the available elements could form.
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
    /// IC names the caller marked "of interest" (TNode::GEM_set_elements_of_interest()); every
    /// other trace IC is a default seed. Read by StrandedElementCheck() and
    /// TNode::GEM_trace_regimes(); never by the solve. Trailing member, so no existing offset
    /// moves. Copied by copyMULTI().
    std::vector<std::string> elementsOfInterest;

    /// True when ordinary IC i (not a charge row) is numerically trace on the Optima path:
    /// one floor amount does not fit inside its own mass-balance tolerance,
    /// B_i*pa_DHB < floor, the floor being pa_OptimaDcFloor if set, else pa_DHB, in internal
    /// units. With pa_DG > 1e-5 this is B_i/sum(B) < 1/pa_DG; otherwise B_i < 1 mol.
    /// GEM_trace_regimes() keeps its own chemical rule and does not use this.
    /// In plain words: tells whether an element is so dilute that the Optima solver treats
    /// its balance as approximate.
    bool ICIsNumericalTrace( long int i ) const;
    /// A default seed is a numerically trace IC that the caller did not mark of interest,
    /// once any IC is marked. With nothing marked there are no default seeds. Used by
    /// pa_OptimaZeroAbsent's rebalance test and by the CERT record (reported as mb_seed_rel,
    /// not scored in mb_pass).
    /// In plain words: a trace element the user did not ask about, whose exact balance is
    /// not checked strictly.
    bool ICIsDefaultSeed( long int i ) const;
    /// Warns, on the Optima path, when an ordinary IC's bulk amount is smaller than the least
    /// its carrier species can hold at the lower bound,
    ///   need_i = sum_j a(i,j) * max(DLL_j, floor)   over species with a(i,j) > 0,
    /// so that no point inside the bounds satisfies that IC's mass balance. The solve still
    /// runs; the residual is left to the post-solve repair. Read-only; called once per Optima
    /// call after rescaling, and reports amounts in real moles.
    /// In plain words: warns when an element is present in such a tiny amount that the
    /// Optima solver cannot represent it exactly.
    void SubFloorElementCheck( double dcFloor ) const;
};

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
/// Event trace of the solvers (IPM/MBR/PSSC decisions, run headers, results).
/// Returns the open trace stream when GEMS3K_NATIVE_TRACE_FILE=<path> is set, nullptr
/// otherwise, so each call site costs one null test. Opened once in append mode and kept
/// open for the process lifetime. Works in Release builds (not NDEBUG-gated).
/// Format: one line per event, tagged in the first field, then KEY=value pairs.
/// In plain words: an optional log file of what the solver decided, switched on by an
/// environment variable.
FILE* native_trace_file();
/// Suppress (on = true) / re-enable the per-call RUN header and KEY result records; nests.
void native_trace_quiet( bool on );

/// Per-iteration IPM descent record, written when GEMS3K_IPM_PROBE=<path> is set;
/// nullptr otherwise.
FILE* ipm_probe_file();

/// Writes the full run configuration - requested mode, T, P, bulk composition with IC
/// names, and every BASE_PARAM field in force - into the trace file, once per
/// TNode::GEM_run() call, for every solver mode.
/// In plain words: records exactly which inputs and settings a run used.
void native_trace_run_header( const MULTI& pm, const BASE_PARAM* pa, long int mode );

/// Writes one DECIDE record: a choice the solver made during the call (as opposed to
/// what it was configured with, SET/EFF, or what it returned, KEY). Format
/// `DECIDE <what> <key>=<value> ...`, no timestamps, so traces of identical runs are
/// identical. Zero cost when GEMS3K_NATIVE_TRACE_FILE is unset.
/// In plain words: notes in the trace which optional mechanisms actually fired.
void native_trace_decide( const char* fmt, ... );

/// Writes the outcome of a completed solve into the trace file: the present phase
/// assemblage (names, amounts, molar volumes and an order-independent hash of the name
/// set), pH, pe and ionic strength. Once per TNode::GEM_run() call, after the dispatch.
/// In plain words: a short fingerprint of the answer, so two runs can be compared at a glance.
void native_trace_run_result( const MULTI& pm, long int mode, long int status, TMultiBase* mb = nullptr );

/// Result of the Optima free-dual search, kept for the CERT trace record. Three-valued:
/// -1 = no search on this call, 0 = searched and not resolved, 1 = searched and resolved.
/// Reset at each run header so a value never carries over from the previous call.
/// In plain words: remembers whether the solver managed to confirm the answer on this call.
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
    f_pa_OptimaAcceptRepair, f_pa_NativeOptimaFinish, f_pa_OptimaCollapseCheck,
    f_pa_OptimaFastFirst, f_pa_OptimaFloorCurvature, f_pa_OptimaEmptySolvent, f_pa_OptimaTpdVeto, f_pa_OptimaEmptyPreSolveSkip, f_pa_OptimaTwinReadmit, f_pa_OptimaDryRetry, f_pa_OptimaCallCap, f_pa_OptimaCapNudges,
    f_pa_OptimaFinishFirst, f_pa_OptimaFastFrozen, f_pa_OptimaFastDiverge, f_pa_OptimaRedoxFree, f_pa_OptimaFastRises,
    f_pa_OptimaPreSolveSeed, f_pa_OptimaColdLevel, f_pa_OptimaLevelStart, f_pa_OptimaFastCollapse, f_pa_OptimaFrozenCut,
    f_pa_OptimaDryOut, f_pa_OptimaTailStop, f_pa_OptimaTailStopWarm, f_pa_OptimaTailDryOnly, f_pa_OptimaWaterKeep,
    f_pa_OptimaSnapLeaving, f_pa_OptimaPartialPivLU, f_pa_IpmQrSortPivot, f_pa_OptimaLSWarmTrials,
    f_pa_TpdStrongBasins, f_pa_IpmLevelStart, f_pa_OptimaFastAcceptCut, f_pa_OptimaLevelFinish, f_pa_IpmNoProgress, f_pa_OptimaFDPureSkip, f_pa_OptimaFinishStall,
    f_pa_OptimaTraceRespread, f_pa_FailRescue, f_pa_RescueWalkSteps, f_pa_RescueBudget, f_pa_OptimaStallNetFactor

} MULTI_DYNAMIC_FIELDS;

enum volume_code {  /* Codes of volume parameter ??? */
    VOL_UNDEF, VOL_CALC, VOL_CONSTR
};

#endif   //_ms_multi_h

