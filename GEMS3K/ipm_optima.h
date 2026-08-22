//-------------------------------------------------------------------
// $Id$
//
/// \file ipm_optima.h
/// Generic "control condition" record for the Optima-based joint
/// equilibrium solve (TMultiBase::CalculateEquilibriumStateOptima(),
/// ipm_optima.cpp). Reaktoro's EquilibriumSpecs implements pH, Eh,
/// fixed fugacity, and fixed species activity all through one
/// mechanism: an implicit titrant unknown, coupled into the mass
/// balance via its own stoichiometry, whose objective gradient is
/// pinned to a closed-form value derived from the requested target,
/// solved jointly with the ordinary equilibrium unknowns in one
/// Newton system (see GEMS3K/CLAUDE.md, "Reaktoro's actual pH/Eh
/// mechanism, read from source"). EqControlCondition is that same
/// generic record for GEMS3K - pH and Eh are just the two built-in
/// condition types constructed by
/// TMultiBase::SetControlCondition_pH()/_Eh(); a future condition
/// (fixed fugacity, fixed activity of a named species, ...) is added
/// by constructing one more EqControlCondition the same way, without
/// touching the Newton-system assembly in ipm_optima.cpp at all.
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
/// row(s) to the chemical potential implied by `target` - ordinary open-
/// system/reservoir thermodynamics, solved in one Newton system together
/// with the rest of the equilibrium (see ipm_optima.cpp,
/// TMultiBase::CalculateEquilibriumStateOptima()).
struct EqControlCondition
{
    std::string name;    ///< diagnostic label, e.g. "pH", "Eh" - also used by
                          ///< SetControlCondition_pH()/_Eh() to replace a
                          ///< previously-registered condition of the same name
    double target = 0.;  ///< requested value, in the condition's own natural units (pH units, V, ...)
    double tolerance = 1e-3; ///< achieved-vs-target tolerance for the post-solve verification check (see CalculateEquilibriumStateOptima()'s targetMet check)

    /// (IC row index, coefficient) pairs: this condition's Aex column,
    /// i.e. how one unit of its titrant unknown contributes to each
    /// mass-balance row. Resolved once, at SetControlCondition_*() time
    /// (index lookups only, valid immediately after GEM_init() - see
    /// ResolveControlConditionIndices() in ipm_optima.cpp).
    std::vector<std::pair<long int,double>> stoich;

    /// Maps `target` to the fixed objective-gradient value pinned for this
    /// unknown. Evaluated once per solve, right before the Newton system is
    /// assembled (not at SetControlCondition_*() time) - it typically reads
    /// call-invariant constants (G0[], T) that are only valid once thermo
    /// data has been loaded for this call's (T,P), which happens inside
    /// CalculateEquilibriumStateOptima() itself, after SetControlCondition_*()
    /// has already returned.
    std::function<double(double target)> fixedGradientFn;

    /// Inverse of fixedGradientFn's own linear combination: maps the
    /// achieved `sum_i U[i]*stoich[i]` (the same quantity the KKT
    /// stationarity condition pins to fixedGradientFn(target)) back to the
    /// condition's own units, for the post-solve "was the target actually
    /// met" verification - see CalculateEquilibriumStateOptima()'s comment
    /// on why Optima::Result::succeeded alone is not sufficient.
    std::function<double(double achievedMuj)> achievedValueFn;

    long int slot = -1;  ///< assigned unknown index (L + k) once active - filled in by CalculateEquilibriumStateOptima(), not by the caller

    // Outputs, filled in by CalculateEquilibriumStateOptima() after solving -
    // meaningless before that (left at their construction-time defaults).
    double titrantAmount = 0.; ///< the solved-for titrant unknown itself, mol (or its condition-specific unit)
    double achievedValue = 0.; ///< achievedValueFn() evaluated on the final dual solution, in the condition's own units - compare against `target`
    bool targetMet = false;    ///< |achievedValue - target| <= tolerance, evaluated after the solve
};

#endif // USE_OPTIMA_SOLVER

#endif // IPM_OPTIMA_H
