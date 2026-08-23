//-------------------------------------------------------------------
// $Id$
//
/// \file ipm_optima.cpp
/// Chemical equilibrium via the Optima library's general primal-dual
/// interior-point NLP solver, as an alternative to GEMS3K's own
/// IPM/MBR loop. Reuses GEMS3K's own G0[]/activity-coefficient
/// machinery as the objective function fed to Optima; Optima owns the
/// Newton iteration and box/linear-equality handling entirely.
///
/// Two new solver modes are dispatched here from TNode::GEM_run()
/// (node.cpp): AOP (cold/AIA-equivalent start) and SOP (warm/SIA-
/// equivalent start) - see NODECODECH in databr.h. Both call the one
/// CalculateEquilibriumStateOptima() below, which reads pm.pNP (set by
/// GEM_run() before the call, exactly like the native AIA/SIA path)
/// to choose the starting point.
///
/// Optional equilibrium "control conditions" (pH, Eh, and - by design -
/// any future closed-form target such as fixed fugacity or fixed
/// species activity) are folded into the same joint Newton solve as
/// extra unknowns with a fixed objective-gradient, matching Reaktoro's
/// own EquilibriumSpecs mechanism (see ipm_optima.h). With none
/// registered this is a plain equilibrium solve, architecturally
/// equivalent to AIA/SIA but solved via Optima instead of IPM/MBR.
///
/// Scope, deliberately: no phase-stability (PSSC) pass, no volume/T,P-
/// as-unknown control conditions.
//
// Copyright (c) 1992-2026
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

#include "ms_multi.h"

#ifdef USE_OPTIMA_SOLVER

#include "node.h"
#include <Optima/Optima.hpp>
#include <algorithm>
#include <cmath>
#include <ctime>

namespace {
// Molality<->mole-fraction scale correction for the pH formula
// (ipm_chemical2.cpp's lnFmol = log(1000/MMC)), hardcoded for the common
// single/dominant-water-solvent case (MMC = water's molar mass)
const double kLnFmol = std::log(1000.0 / 18.01528);
}

void TMultiBase::SetControlCondition_pH( double pH_target, double tolerance )
{
    if( node1 == nullptr )
        Error( "E93IPM: Optima control condition: ", "SetControlCondition_pH() requires an initialized TNode (call after GEM_init())" );

    long int xH  = node1->IC_name_to_xCH( "H" );
    long int xZz = node1->IC_name_to_xCH( "Zz" );
    long int xHplusDC = node1->DC_name_to_xCH( "H+" );
    if( xH < 0 || xZz < 0 || xHplusDC < 0 )
        Error( "E93IPM: Optima control condition: ", "SetControlCondition_pH(): could not resolve H/Zz IC or H+ DC index in this project" );

    optima_control_conditions.erase(
        std::remove_if( optima_control_conditions.begin(), optima_control_conditions.end(),
                         [](const EqControlCondition& c){ return c.name == "pH"; } ),
        optima_control_conditions.end() );

    EqControlCondition c;
    c.name = "pH";
    c.target = pH_target;
    c.tolerance = tolerance;
    c.stoich = { { xH, 1.0 }, { xZz, 1.0 } };
    // Default tolerance (used only when the caller passes tolerance<0):
    // pa_p->GAS is a chem.pot.-difference threshold (mol/mol); pH's own
    // formula below is linear in Muj with slope -ln_to_lg, so a Muj-space
    // tolerance of GAS maps to a pH-space tolerance of ln_to_lg*GAS.
    c.defaultToleranceFn = [this]() { return std::fabs( ln_to_lg ) * base_param()->GAS; };
    // pH = -ln_to_lg*(Muj - G0[H+] + lnFmol)  =>  Muj_target = G0[H+] - lnFmol - pH_target/ln_to_lg
    // (ipm_chemical2.cpp's own pH formula; Muj = sum_i U[i]*a(H+,i), which is
    // exactly sum_i U[i]*stoich[i] above, since stoich matches H+'s own
    // stoichiometry)
    c.fixedGradientFn = [this, xHplusDC]( double target )
    { return pm.G0[xHplusDC] - kLnFmol - target / ln_to_lg; };
    c.achievedValueFn = [this, xHplusDC]( double achievedMuj )
    { return -ln_to_lg * ( achievedMuj - pm.G0[xHplusDC] + kLnFmol ); };

    optima_control_conditions.push_back( c );
}

void TMultiBase::SetControlCondition_Eh( double Eh_target, double tolerance )
{
    if( node1 == nullptr )
        Error( "E93IPM: Optima control condition: ", "SetControlCondition_Eh() requires an initialized TNode (call after GEM_init())" );

    long int xZz = node1->IC_name_to_xCH( "Zz" );
    if( xZz < 0 )
        Error( "E93IPM: Optima control condition: ", "SetControlCondition_Eh(): could not resolve Zz IC index in this project" );

    optima_control_conditions.erase(
        std::remove_if( optima_control_conditions.begin(), optima_control_conditions.end(),
                         [](const EqControlCondition& c){ return c.name == "Eh"; } ),
        optima_control_conditions.end() );

    EqControlCondition c;
    c.name = "Eh";
    c.target = Eh_target;
    c.tolerance = tolerance;
    // Direct electron/charge (Zz-row) titrant, matching Reaktoro's own
    // qvar.substance="e-" choice. An O2 titrant (mass-balance coupling on
    // the O row only) would have no direct algebraic relationship to Zz's
    // dual/Eh, so it wouldn't admit a closed-form fixed gradient here -
    // the direct charge titrant is what makes the one-shot joint solve
    // possible.
    c.stoich = { { xZz, -1.0 } };
    // Default tolerance (used only when the caller passes tolerance<0):
    // same pa_p->GAS reuse as SetControlCondition_pH(), mapped through Eh's
    // own slope (0.000086*T) instead of pH's. Lazy (reads pm.T at solve
    // time, not here) for the same reason fixedGradientFn below is lazy.
    c.defaultToleranceFn = [this]() { return 0.000086 * pm.T * base_param()->GAS; };
    // Eh = 0.000086*U[Zz]*T  =>  U[Zz]_target = Eh_target/(0.000086*T);
    // achievedMuj = sum_i U[i]*stoich[i] = -U[Zz], so the fixed gradient
    // (pinned to achievedMuj, via the KKT stationarity of this slot) is
    // -U[Zz]_target.
    c.fixedGradientFn = [this]( double target )
    { return -( target / ( 0.000086 * pm.T ) ); };
    c.achievedValueFn = [this]( double achievedMuj )
    { return 0.000086 * ( -achievedMuj ) * pm.T; };

    optima_control_conditions.push_back( c );
}

void TMultiBase::ClearControlConditions()
{
    optima_control_conditions.clear();
}

double TMultiBase::OptimaMaxMassBalanceResidual()
{
    double maxres = 0.;
    for( long int i = 0; i < pm.N; i++ )
    {
        double sum = 0.;
        for( long int j = 0; j < pm.L; j++ )
            sum += pm.A[ i + j*pm.N ] * pm.Y[j];
        maxres = std::max( maxres, std::fabs( sum - pm.B[i] ) );
    }
    return maxres;
}

double TMultiBase::CalculateEquilibriumStateOptima( long int& NumIterFIA, long int& NumIterIPM )
{
    double ScFact = 1.;
    const BASE_PARAM* pa_p = base_param();
    // Reuse GEMS3K's own tolerances rather than adding parallel
    // Optima-specific BASE_PARAM fields: DHB is the numerical DC-amount
    // floor, IIM/DK are the "max iterations"/"convergence tolerance"
    // knobs (below), DW gates hard-error vs. soft-BAD on non-convergence.
    const double dcFloor = std::max( pa_p->DHB, 1e-300 );

    InitalizeGEM_IPM_Data();

    pm.t_start = clock();
    pm.t_end = pm.t_start;
    pm.t_elap_sec = 0.0;
    pm.ITF = pm.ITG = 0;
    pm.Ec = pm.MK = pm.PZ = 0;
    setErrorMessage( 0, "", "" );

    if( pa_p->DG > 1e-5 )
    {
        ScFact = SystemTotalMolesIC();
        ScaleSystemToInternal( ScFact );
    }

    try
    {
        // Allocates/parametrizes each multicomponent phase's TSolMod
        // instance (phSolMod[k]) at the current T,P - a prerequisite for
        // any later LINK_UX_MODE call. The native path gets this from its
        // own GEM_IPM_Init() (ipm_simplex.cpp); this path has no such
        // step, so it must be called explicitly here.
        CalculateActivityCoefficients( LINK_TP_MODE );

        // pm.pNP is set by TNode::GEM_run() before this call, exactly as it
        // is for the native AIA/SIA path: 0 = cold start (AOP, like AIA),
        // 1 = warm start reusing the existing pm.Y[] (SOP, like SIA).
        if( pm.pNP == 0 )
        {
            AutoInitialApproximation();
            DC_RaiseZeroedOff( 0, pm.L );
        }
        TotalPhasesAmounts( pm.Y, pm.YF, pm.YFA );
        for( long int j = 0; j < pm.L; j++ )
            pm.X[j] = pm.Y[j];
        TotalPhasesAmounts( pm.X, pm.XF, pm.XFA );
        CalculateActivityCoefficients( LINK_UX_MODE );

        const long int N = pm.N;
        const long int L = pm.L;
        std::vector<EqControlCondition>& conditions = optima_control_conditions;
        const long int R = (long int)conditions.size();

        // Resolve each active condition's fixed objective-gradient value
        // once, up front (needs this call's G0[]/T, only valid now that
        // InitalizeGEM_IPM_Data() has run) - assign its unknown slot.
        std::vector<double> fixedGrad( R, 0. );
        for( long int k = 0; k < R; k++ )
        {
            conditions[k].slot = L + k;
            fixedGrad[k] = conditions[k].fixedGradientFn( conditions[k].target );
        }

        Optima::Dims dims;
        dims.x  = L + R;
        dims.be = N;

        Optima::Problem problem( dims );

        for( long int i = 0; i < N; i++ )
            for( long int j = 0; j < L; j++ )
                problem.Aex(i, j) = pm.A[ i + j*N ];
        for( long int k = 0; k < R; k++ )
        {
            for( long int i = 0; i < N; i++ )
                problem.Aex(i, L+k) = 0.;
            for( const auto& rc : conditions[k].stoich )
                problem.Aex( rc.first, L+k ) = rc.second;
        }
        for( long int i = 0; i < N; i++ )
            problem.be[i] = pm.B[i];

        // Box constraints: reuse GEMS3K's own DC kinetic-restriction
        // bounds directly, floored by dcFloor (pa_DHB) for well-posedness
        // of the log-based Hessian approximation.
        for( long int j = 0; j < L; j++ )
        {
            problem.xlower[j] = std::max( pm.DLL[j], dcFloor );
            problem.xupper[j] = ( pm.DUL[j] > 0. && pm.DUL[j] < 1e6 )
                                 ? pm.DUL[j] : std::max( pm.SMols, 1.0 ) * 10.;
        }
        // Titrant unknowns are free (can be positive or negative), unlike
        // ordinary species amounts - bounded only generously, as a
        // numerical safety net, not a physical constraint.
        const double titrantBound = std::max( pm.SMols, 1.0 ) * 2.0;
        for( long int k = 0; k < R; k++ )
        {
            problem.xlower[L+k] = -titrantBound;
            problem.xupper[L+k] =  titrantBound;
        }

        // Objective: Gibbs energy, in the same RT-normalized units GEMS3K's
        // own G0[]/G[] already use. fx = chemical potentials, recomputed
        // via GEMS3K's real activity-coefficient machinery every Optima
        // Newton iteration (not frozen/linearized the way MBR's Schur-
        // complement weights are). Virtual control-condition slots get a
        // FIXED (not composition-recomputed) gradient equal to their
        // target's own closed-form chemical potential - the whole
        // mechanism (see ipm_optima.h).
        problem.f = [this, L, R, dcFloor, &fixedGrad]
                    ( Optima::ObjectiveResultRef res, Optima::VectorView x,
                      Optima::VectorView /*p*/, Optima::VectorView /*c*/,
                      Optima::ObjectiveOptions opts )
        {
            for( long int j = 0; j < L; j++ )
                pm.X[j] = x[j];
            TotalPhasesAmounts( pm.X, pm.XF, pm.XFA );
            CalculateActivityCoefficients( LINK_UX_MODE );
            PrimalChemicalPotentials( pm.F, pm.X, pm.XF, pm.XFA );

            double fval = 0.;
            for( long int j = 0; j < L; j++ )
            {
                res.fx[j] = pm.F[j];
                fval += pm.X[j] * pm.F[j];
            }
            for( long int k = 0; k < R; k++ )
            {
                res.fx[L+k] = fixedGrad[k];
                fval += x[L+k] * fixedGrad[k];
            }
            res.f = fval;
            if( opts.eval.fxx )
            {
                res.fxx.setZero();
                for( long int j = 0; j < L; j++ )
                    res.fxx(j,j) = 1.0 / std::max( pm.X[j], dcFloor );
                // Virtual slots get zero direct curvature - a fixed
                // gradient has no second derivative in its own value; all
                // coupling to the rest of the system comes through the
                // shared Aex mass-balance rows.
                res.diagfxx = true;
            }
            res.succeeded = true;
        };

        Optima::State state( dims );
        for( long int j = 0; j < L; j++ )
            state.x[j] = std::max( pm.Y[j], dcFloor );
        for( long int k = 0; k < R; k++ )
            state.x[L+k] = 0.0;

        Optima::Options options;
        // IIM/DK are GEMS3K's own "max iterations"/"convergence tolerance"
        // knobs for the native IPM loop - reused directly here for
        // Optima's equivalent settings rather than adding parallel fields.
        options.maxiters = (unsigned)std::max( (short)1, pa_p->IIM );
        options.convergence.tolerance = pa_p->DK;

        Optima::Solver solver;
        solver.setOptions( options );
        Optima::Result result = solver.solve( problem, state );

        ipm_logger->info( "CalculateEquilibriumStateOptima: pNP={} succeeded={} iterations={} nConditions={}",
                           pm.pNP, result.succeeded, result.iterations, R );

        for( long int j = 0; j < L; j++ )
        {
            pm.Y[j] = state.x[j];
            pm.X[j] = state.x[j];
        }
        // Optima's Lagrange multipliers ye use the opposite sign convention
        // from GEMS3K's own dual U[] - negate here rather than resign every
        // place U[] is used below.
        for( long int i = 0; i < N; i++ )
            pm.U[i] = -state.ye[i];

        TotalPhasesAmounts( pm.X, pm.XF, pm.XFA );
        CalculateActivityCoefficients( LINK_UX_MODE );

        // pH/Eh straight from the dual solution - same closed form as
        // ConCalcDC() (ipm_chemical2.cpp) - always computed, regardless of
        // whether either was an active control condition, as a diagnostic.
        long int protonIdx = -1;
        for( long int j = 0; j < L; j++ )
            if( pm.DCC[j] == DC_AQ_PROTON ) { protonIdx = j; break; }
        if( protonIdx >= 0 )
        {
            double Muj = 0.;
            for( long int i = 0; i < N; i++ )
                Muj += pm.U[i] * pm.A[ i + protonIdx*N ];
            pm.pH = -ln_to_lg * ( Muj - pm.G0[protonIdx] + kLnFmol );
        }
        pm.pe = ln_to_lg * pm.U[N-1];
        pm.Eh = 0.000086 * pm.U[N-1] * pm.T;

        // Commit each active condition's titrant into pm.B[] - without
        // this, CheckMassBalanceResiduals() below (and every downstream
        // consumer of pm.B[], including TNode::packDataBr()'s CNode->bIC[]
        // output) would compare Y against the ORIGINAL, un-titrated bulk
        // composition, not the effective one Optima actually solved
        // against. The caller sees the titrated bulk composition on
        // output.
        // Sign: Optima's linear equality is Aex*x_ext=be, i.e.
        // A*Y + sum_k stoich_k*xi_k = B (the titrant column sits on the
        // SAME side as A*Y, not folded into the RHS) - so the effective
        // bulk composition Y actually balances against is
        // B - sum_k stoich_k*xi_k, hence subtract here, not add.
        for( long int k = 0; k < R; k++ )
        {
            conditions[k].titrantAmount = state.x[L+k];
            for( const auto& rc : conditions[k].stoich )
                pm.B[ rc.first ] -= rc.second * conditions[k].titrantAmount;
        }

        // Verify each active condition actually reached its target - not
        // just that Optima's own KKT residual reports solved. At a
        // bound-active titrant slot, the box constraint's own dual
        // absorbs the KKT residual, so the fixed-gradient condition this
        // mechanism relies on can be unsatisfied even while Optima
        // reports succeeded=true.
        bool allTargetsMet = true;
        std::string targetMissBuf;
        for( long int k = 0; k < R; k++ )
        {
            double achievedMuj = 0.;
            for( const auto& rc : conditions[k].stoich )
                achievedMuj += pm.U[ rc.first ] * rc.second;
            conditions[k].achievedValue = conditions[k].achievedValueFn( achievedMuj );
            double tol = conditions[k].tolerance;
            if( tol < 0. && conditions[k].defaultToleranceFn )
                tol = conditions[k].defaultToleranceFn();
            conditions[k].targetMet = std::fabs( conditions[k].achievedValue - conditions[k].target ) <= tol;
            if( !conditions[k].targetMet )
            {
                allTargetsMet = false;
                targetMissBuf += conditions[k].name + ": target=" + std::to_string(conditions[k].target)
                               + " achieved=" + std::to_string(conditions[k].achievedValue) + "; ";
            }
        }

        // Mass-balance trustworthiness check, reusing GEMS3K's own
        // per-IC tolerance (pm.DHBM), exactly as the native solver's own
        // testMulti()/CheckMassBalanceResiduals() path does - after the
        // pm.B[] commit above, so this compares against the effective
        // (post-titration) bulk composition, not the pre-titration input.
        long int massBalanceBadIC = CheckMassBalanceResiduals( pm.Y );

        if( !result.succeeded && pa_p->DW )
        {
            // DW already gates exactly this decision for the native
            // solver's own "MBR iterations exceeded" case (ipm_main.cpp);
            // reused as-is here rather than adding a parallel flag.
            Error( "E90IPM: Optima solver: ", "Optima::Solver::solve() did not converge (pa_DW forces this to a hard error)" );
        }
        else if( !result.succeeded || !allTargetsMet || massBalanceBadIC >= 0 )
        {
            // Soft failure: a state was produced, but it is not fully
            // trustworthy - mirrors AIA/SIA's own BAD_GEM_* semantics
            // (testMulti() already reads pm.MK below, via TNode::GEM_run()).
            pm.MK = 2;
            std::string buf = "Optima solve ";
            buf += result.succeeded ? "reported success" : "did not converge";
            if( !allTargetsMet )
                buf += std::string("; control condition(s) not met: ") + targetMissBuf;
            if( massBalanceBadIC >= 0 )
                buf += "; mass balance residual exceeds tolerance for IC " + char_array_to_string(pm.SB[massBalanceBadIC],3);
            setErrorMessage( 21, "W21IPM: Optima solver: ", buf.c_str() );
        }
    }
    catch( TError& xcpt )
    {
        if( pa_p->DG > 1e-5 )
            RescaleSystemFromInternal( ScFact );
        NumIterFIA = pm.ITF;
        NumIterIPM = pm.ITG;
        pm.t_end = clock();
        pm.t_elap_sec = double(pm.t_end - pm.t_start)/double(CLOCKS_PER_SEC);
        Error( xcpt.title, xcpt.mess );
    }

    if( pa_p->DG > 1e-5 )
        RescaleSystemFromInternal( ScFact );

    NumIterFIA = pm.ITF;
    NumIterIPM = pm.ITG;
    pm.t_end = clock();
    pm.t_elap_sec = double(pm.t_end - pm.t_start)/double(CLOCKS_PER_SEC);
    return pm.t_elap_sec;
}

#endif // USE_OPTIMA_SOLVER
