//-------------------------------------------------------------------
// $Id$
//
/// \file ipm_optima.h
/// Generic "control condition" record for the Optima-based joint
/// equilibrium solve (TMultiBase::CalculateEquilibriumStateOptima(),
/// ipm_optima.cpp). A control condition (pH, Eh, ...) is an implicit
/// titrant unknown, coupled into the mass balance via its own
/// stoichiometry, whose objective gradient is pinned to a closed-form
/// value derived from the requested target and solved jointly with the
/// ordinary equilibrium unknowns in one Newton system. pH and Eh are
/// the two built-in condition types, constructed by
/// TMultiBase::SetControlCondition_pH()/_Eh(); a further condition
/// (fixed fugacity, fixed activity of a named species, ...) is added by
/// constructing one more EqControlCondition the same way, without
/// touching the Newton-system assembly in ipm_optima.cpp.
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
#ifndef IPM_OPTIMA_H
#define IPM_OPTIMA_H

#ifdef USE_OPTIMA_SOLVER

#include <functional>
#include <string>
#include <utility>
#include <vector>

/// One generic equilibrium control condition - an implicit titrant unknown
/// appended to the primal vector, coupled into the mass-balance system via
/// `stoich`, whose objective-gradient entry is pinned (not recomputed from
/// composition) to `fixedGradientFn(target)`. At the joint KKT-stationary
/// point this forces whatever real species share the same mass-balance
/// row(s) to the chemical potential implied by `target`, solved in one
/// Newton system together with the rest of the equilibrium (see
/// ipm_optima.cpp, TMultiBase::CalculateEquilibriumStateOptima()).
struct EqControlCondition
{
    std::string name;    ///< diagnostic label, e.g. "pH", "Eh" - also used by
                          ///< SetControlCondition_pH()/_Eh() to replace a
                          ///< previously-registered condition of the same name
    double target = 0.;  ///< requested value, in the condition's own natural units (pH units, V, ...)
    double tolerance = -1.; ///< achieved-vs-target tolerance for the post-solve verification check
                             ///< (see CalculateEquilibriumStateOptima()'s targetMet check). A negative
                             ///< value (the default) means "not explicitly set" - the check falls back
                             ///< to defaultToleranceFn() below instead of a hardcoded constant.

    /// Called only when `tolerance < 0`: derives the default from GEMS3K's
    /// own pa_p->GAS ("threshold for primal-dual chem.pot. difference
    /// (mol/mol), also used by PhaseSelectionSpeciationCleanup()"),
    /// converted into this condition's own units via the same linear map
    /// as fixedGradientFn/achievedValueFn. Evaluated at solve time (not at
    /// SetControlCondition_*() time) since the Eh case reads pm.T, which is
    /// only current right before the Newton system is assembled.
    std::function<double()> defaultToleranceFn;

    /// (IC row index, coefficient) pairs: this condition's Aex column, i.e.
    /// how one unit of its titrant unknown contributes to each mass-balance
    /// row. Resolved once, at SetControlCondition_*() time.
    std::vector<std::pair<long int,double>> stoich;

    /// Maps `target` to the fixed objective-gradient value pinned for this
    /// unknown. Evaluated once per solve, right before the Newton system is
    /// assembled (not at SetControlCondition_*() time), since it reads
    /// call-invariant constants (G0[], T) that are only valid once thermo
    /// data has been loaded for this call's (T,P).
    std::function<double(double target)> fixedGradientFn;

    /// Inverse of fixedGradientFn's own linear combination: maps the
    /// achieved `sum_i U[i]*stoich[i]` (the quantity the KKT stationarity
    /// condition pins to fixedGradientFn(target)) back to the condition's
    /// own units, for the post-solve "was the target actually met" check.
    std::function<double(double achievedMuj)> achievedValueFn;

    long int slot = -1;  ///< assigned unknown index (L + k) once active - filled in by CalculateEquilibriumStateOptima(), not by the caller

    // Outputs, filled in by CalculateEquilibriumStateOptima() after solving -
    // meaningless before that (left at their construction-time defaults).
    /// The solved-for titrant unknown, in REAL moles of the recipe (not pa_DG's internal
    /// scale - see the commit loop in CalculateEquilibriumStateOptima()). Sign: the titrant
    /// column sits on the species side, A*Y + stoich*xi = B, so the amount ADDED to the input
    /// bulk composition is -stoich*titrantAmount - for pH (stoich H:+1, Zz:+1) that is
    /// -xi mol of H+; for Eh (stoich Zz:-1) it is +xi on the Zz row, i.e. xi mol of electrons
    /// REMOVED. The returned bIC already carries it: bIC_out = bIC_in - sum_k stoich_k*xi_k.
    double titrantAmount = 0.;
    double achievedValue = 0.; ///< achievedValueFn() evaluated on the final dual solution, in the condition's own units - compare against `target`
    bool targetMet = false;    ///< |achievedValue - target| <= tolerance, evaluated after the solve
};

#endif // USE_OPTIMA_SOLVER

#endif // IPM_OPTIMA_H
