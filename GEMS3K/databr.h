//-------------------------------------------------------------------
// $Id$
/// \file databr.h
/// Contains definition of the DATABR structure - data bridge between
/// GEMS3K and another code.
//
/// \struct DATABR databr.h
/// DataBRidge defines the structure of node-dependent data for
/// exchange between the coupled GEM IPM and FMT code parts.
/// DATABR structure is used in TNode and TNodeArray classes.
//
// Copyright (c) 2003-2011 by D.Kulik, S.Dmytriyeva, F.Enzmann, W.Pfingsten
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
#ifndef DataBr_H_
#define DataBr_H_



typedef struct  /// DATABR - template node data bridge structure
{
   long int
     NodeHandle,    ///< Node identification handle, not used in calculaions on TNode level
     NodeTypeHY,    ///< Node type code (hydraulic), not used on TNode level see NODETYPE
     NodeTypeMT,    ///< Node type (mass transport), not used on TNode level see NODETYPE
     NodeStatusFMT, ///< Node status code in FMT part, not used on TNode level see NODECODEFMT
     NodeStatusCH,  ///< Node status code in GEM (input and output) see NODECODECH
     IterDone;      ///< Number of iterations performed by GEM IPM in the last run

/*  these important data array dimensions are provided in the DATACH structure
   long int
    nICb,       ///< Number of Independent Components kept in the DATABR memory structure (<= nIC)
    nDCb,      	///< Number of Dependent Components kept in the DATABR memory structure (<=nDC)
    nPHb,     	///< Number of Phases to be kept in the DATABR structure (<= nPH)
    nPSb,       ///< Number of Phases-solutions (multicomponent phases) to be kept in the DATABR memory structure (<= nPS)
*/
//      Usage of this variable (DB - data bridge)                           MT-DB DB-GEM GEM-DB DB-MT
   double
// \section Chemical scalar variables
    TK,     ///< Node temperature T (Kelvin)                     	         +      +      -     -
    P, 	    ///< Node Pressure P (Pa)                         	             +      +      -     -
    Vs,     ///< Volume V of reactive subsystem  (m3)                       (+)    (+)     +     +
    Vi,     ///< Volume of inert subsystem (m3)          	                 +      -      -     +
    Ms,     ///< Mass of reactive subsystem (kg)         	                 +     (+)     -     -
    Mi,     ///< Mass of inert subsystem (kg)             	                 +      -      -     +

    Gs,     ///< Total Gibbs energy of the reactive subsystem (J/RT) (norm)  -      -      +     +
    Hs, 	///< Total enthalpy of reactive subsystem (J) (reserved)         -      -      +     +
    Hi,     ///< Total enthalpy of inert subsystem (J) (reserved)            +      -      -     +

    IC,     ///< Effective aqueous ionic strength (molal)                    -      -      +     +
    pH,     ///< pH of aqueous solution in the activity scale (-log10 molal) -      -      +     +
    pe,     ///< pe of aqueous solution in the activity scale (-log10 molal) -      -      +     +
    Eh,     ///< Eh of aqueous solution (V)                                  -      -      +     +
    Tm,     ///< Actual total simulation time (s)                            +      +      -     -
    dt      ///< Actual time step (s) - needed for TKinMet, can change!      +      +     (+)   (+)
#ifndef NO_NODEARRAYLEVEL
      ,
// \section  FMT variables (units or dimensionsless) - to be used for storing them
//  at the nodearray level, normally not used in the single-node FMT-GEM coupling
    Dif,    ///< General diffusivity of disolved matter (m2/s)
    Vt,		///< Total volume of the node (m3) (Vs + Vi)
    vp,		///< Advection velocity (in pores)  (m/s)
    eps,	///< Effective (actual) porosity normalized to 1
    Km,		///< Actual permeability (m2)
    Kf,		///< Actual Darcy`s constant (m2/s)
    S,		///< Specific storage coefficient, dimensionless (default 1.0)
            ///<  if not 1.0 then can be used as mass scaling factor relative to Ms for the
            ///<  bulk composition/speciation of reactive sub-system in this node
    Tr,     ///< transmissivity m2/s
    h,		///< Actual hydraulic head (hydraulic potential) (m)
    rho,	///< Actual carrier density for density-driven flow (kg/m3)
    al,		///< Specific longitudinal dispersivity of porous media (m)
    at,		///< Specific transversal dispersivity of porous media (m)
    av,		///< Specific vertical dispersivity of porous media (m)
    hDl,	///< Hydraulic longitudinal dispersivity (m2/s)
    hDt,	///< Hydraulic transversal dispersivity (m2/s)
    hDv,	///< Hydraulic vertical dispersivity (m2/s)
    nto	    ///< Tortuosity factor (dimensionless)
#endif
   ;
// \section Data arrays - dimensions nICb, nDCb, nPHb, nPSb see in the DATACH structure
// exchange of values occurs through lists of indices, e.g. xIC, xDC, xPH from DATACH

//      Usage of this variable (DB = data bridge)                                               MT-DB DB-GEM GEM-DB DB-MT
   double
// IC (stoichiometry units)
    *bIC,  ///< Bulk composition of (reactive part of) the system (moles)[nICb]                      +      +      -     -
    *rMB,  ///< Mass balance residuals (moles) [nICb]                                                -      -      +     +
    *uIC,  ///< Chemical potentials of ICs (dual GEM solution) normalized scale (mol/mol)[nICb]      -      -      +     +
// DC (species) in reactive subsystem
    *xDC,  ///< Speciation - amounts of DCs in equilibrium state - primal GEM solution(moles)[nDCb] (+)    (+)     +     +
    *gam,  ///< Activity coefficients of DCs in their respective phases  [nDCb]                     (+)    (+)     +     +
// Metastability/kinetic controls
    *dul,  ///< Upper additional metastability restrictions AMR on amounts of DCs (moles) [nDCb]     +      +      -     -
    *dll,  ///< Lower AMR on amounts of DCs (moles) [nDCb]                                           +      +      -     -
// Phases in reactive subsystem
    *aPH,  ///< Specific surface areas of phases (m2/kg) [nPHb], can change in TKinMet               +      +     (+)   (+)
    *xPH,  ///< Amounts of phases in equilibrium state (moles) [nPHb]                                -      -      +     +
    *vPS,  ///< Volumes of multicomponent phases (m3)  [nPSb]                                        -      -      +     +
    *mPS,  ///< Masses of multicomponent phases (kg)    [nPSb]                                       -      -      +     +
    *bPS,  ///< Bulk elemental compositions of multicomponent phases (moles) [nPSb][nICb]            -      -      +     +
    *xPA,  ///< Amount of carrier (sorbent or solvent) in multicomponent phases [nPSb]               -      -      +     +
           ///< (+) can be used as input in "smart initial approximation" mode of GEM IPM-2 algorithm
    *bSP,  ///< Output bulk composition of the equilibrium solid part of the system, moles   [nICb]  -      -      +     +
   // Metastability/kinetic controls on phases-solutions (added in devPhase branch)
    *amru, ///< Upper AMRs on amounts of multi-component phases (mol) [nPSb]                         +      +      -     -
    *amrl, ///< Lower AMRs on amounts of multi-component phases (mol) [nPSb]                         +      +      -     -
  *omPH; ///< stability (saturation) indices of phases in log10 scale, can change in GEM [nPHb]     (+)    (+)     +     +
}
DATABR;

typedef DATABR*  DATABRPTR;

/// \enum NODECODECH NodeStatus codes with respect to GEMIPM calculations
/*typedef*/ enum NODECODECH {
 NO_GEM_SOLVER= 0,   ///< No GEM re-calculation needed for this node
 NEED_GEM_AIA = 1,   ///< Need GEM calculation with LPP (automatic) initial approximation (AIA)
 OK_GEM_AIA   = 2,   ///< OK after GEM calculation with LPP AIA
 BAD_GEM_AIA  = 3,   ///< Bad (not fully trustful) result after GEM calculation with LPP AIA
 ERR_GEM_AIA  = 4,   ///< Failure (no result) in GEM calculation with LPP AIA
 NEED_GEM_SIA = 5,   ///< Need GEM calculation with no-LPP (smart) IA, SIA
                     ///<   using the previous speciation (full DATABR lists only)
 OK_GEM_SIA   = 6,   ///< OK after GEM calculation with SIA
 BAD_GEM_SIA  = 7,   ///< Bad (not fully trustful) result after GEM calculation with SIA
 ERR_GEM_SIA  = 8,   ///< Failure (no result) in GEM calculation with SIA
 T_ERROR_GEM  = 9,   ///< Terminal error has occurred in GEMS3K (e.g. memory corruption). Restart is required.
 // "Optima" modes: equilibrium via the Optima library's general primal-
 // dual interior-point NLP solver (TMultiBase::CalculateEquilibriumStateOptima(),
 // ipm_optima.cpp) instead of GEMS3K's own IPM/MBR loop - only meaningful
 // if GEMS3K was built with USE_OPTIMA_SOLVER; otherwise TNode::GEM_run()
 // logs a warning and falls back to the equivalent native
 // AIA/SIA solve (returning OK/BAD/ERR_GEM_AIA/SIA). AOP
 // mirrors AIA (cold/LPP-simplex start), SOP mirrors SIA (warm start
 // reusing the previous speciation)
 NEED_GEM_AOP = 10,  ///< Need GEM calculation via Optima with cold (AIA-equivalent) initial approximation
 OK_GEM_AOP   = 11,  ///< OK after GEM calculation via Optima with cold initial approximation
 BAD_GEM_AOP  = 12,  ///< Bad (not fully trustful) result after GEM calculation via Optima with cold initial approximation
 ERR_GEM_AOP  = 13,  ///< Failure (no result) in GEM calculation via Optima with cold initial approximation
 NEED_GEM_SOP = 14,  ///< Need GEM calculation via Optima with warm (SIA-equivalent) initial approximation
                     ///<   using the previous speciation (full DATABR lists only)
 OK_GEM_SOP   = 15,  ///< OK after GEM calculation via Optima with warm initial approximation
 BAD_GEM_SOP  = 16,  ///< Bad (not fully trustful) result after GEM calculation via Optima with warm initial approximation
 ERR_GEM_SOP  = 17,  ///< Failure (no result) in GEM calculation via Optima with warm initial approximation
 // "ROP" mode: a single, faithful port of Reaktoro's OWN equilibrium
 // mechanism onto GEMS3K's chemistry (same uniform tiny initial guess,
 // same PartiallyExact Hessian strategy, same untouched Optima::Options
 // defaults, same single apply_min_max_fix_and_accept-toggle fallback -
 // see TMultiBase::CalculateEquilibriumStateOptima()'s referenceMode
 // branch, ipm_optima.cpp) - NOT just AOP's own seed/options swapped in.
 // Unlike AOP/SOP there is no cold/warm pair: Reaktoro's own default
 // equilibrate() always starts from the same uniform seed regardless of
 // any previous state, so ROP is a single mode. Only meaningful if built
 // with USE_OPTIMA_SOLVER; otherwise TNode::GEM_run() falls back to AIA,
 // same as AOP/SOP.
 NEED_GEM_ROP = 18,  ///< Need GEM calculation via Optima, using Reaktoro's own mechanism (uniform seed, PartiallyExact Hessian, untouched Optima defaults)
 OK_GEM_ROP   = 19,  ///< OK after GEM calculation via the ROP mechanism
 BAD_GEM_ROP  = 20,  ///< Bad (not fully trustful) result after GEM calculation via the ROP mechanism
 ERR_GEM_ROP  = 21,  ///< Failure (no result) in GEM calculation via the ROP mechanism
 // "HOP" mode: HYBRID - native GEMS3K IPM/MBR first, then Optima warm-started
 // from its result. In series: native does what it is uniquely good at, which
 // is SELECTING THE SPECIES (its line-search objective GX() truncates any
 // amount below pa_DcMin to exactly 0, so it produces a genuine phase
 // assemblage - see the comment at that truncation in ipm_chemical.cpp and
 // plan-v5 section 60), and Optima then finishes from that assemblage, which
 // is what it is good at.
 //
 // WHY THIS IS A SEPARATE, CALLER-SELECTED MODE and not something AOP does
 // internally: an earlier design had AOP call native itself, and that was
 // reverted 2026-08-23 on explicit direction - AOP/SOP must remain a genuinely
 // switchable, standalone alternative to native, not a combination wearing
 // native's name. A caller asking for HOP is asking for both, by name.
 //
 // WHAT IT BUYS, measured: 07PSIna_G_complex_1_0_1_80_0 (1392 species) -
 // the largest project in the corpus and one that NO Optima mode has ever
 // solved, running out its whole budget - converges here in 28 Optima
 // iterations on top of native's 449, at a G agreeing with native's to 7
 // significant figures.
 //
 // If native fails on the system, there is no assemblage to hand over, and
 // this degrades to a plain cold Optima solve (AOP-equivalent) with a logged
 // warning rather than failing outright. Only meaningful if built with
 // USE_OPTIMA_SOLVER; otherwise it falls back to native AIA alone, as AOP/SOP
 // and ROP do.
 NEED_GEM_HOP = 22,  ///< Need GEM calculation via native IPM/MBR first, then Optima warm-started from its result
 OK_GEM_HOP   = 23,  ///< OK after the hybrid native-then-Optima calculation
 BAD_GEM_HOP  = 24,  ///< Bad (not fully trustful) result after the hybrid native-then-Optima calculation
 ERR_GEM_HOP  = 25,  ///< Failure (no result) in the hybrid native-then-Optima calculation

 // SHP - the WARM (SIA-equivalent) counterpart of HOP, and the pair completes
 // the same cold/warm convention every other solver here already has:
 // AIA/SIA, AOP/SOP. "S" is this API's established marker for a smart (warm)
 // initial approximation; HOP starts with no "A" to swap, so the warm member
 // of the pair carries the S in front instead.
 //
 // WHAT IT IS FOR. HOP as built runs its NATIVE leg cold at EVERY call, which
 // is fine for a single equilibrium and plainly wasteful in a sweep or a
 // transport loop, where the previous point's converged state is sitting right
 // there. SHP starts the native leg warm (pm.pNP = 1, native's own SIA path)
 // and hands its result to the same warm Optima leg HOP already uses, so a
 // sequential step costs a warm native solve plus O(1)-O(10) Optima iterations
 // instead of a full cold native solve plus the same.
 //
 // WHAT IT IS NOT: a new heuristic. SHP is exactly HOP with native's SIA in
 // place of native's AIA, so the caller faces here precisely the choice they
 // already face between AIA and SIA, and between AOP and SOP.
 //
 // THE HEADLINE RESULT is not the sweep saving it was built for. Native's own
 // SIA cannot re-solve its own converged state on a documented set of projects
 // - the cold path returns an answer failing native's own mass-balance test,
 // and SIA is the only path that checks it (plan-v5 sections 29.2 and 60.5).
 // SHP converges on those anyway, at ~4x fewer iterations than HOP, because
 // the state it hands native is the OPTIMA leg's, not native's own, and that
 // state satisfies mass balance to machine precision. Measured at zero
 // distance (re-solve the same composition four times on one node,
 // debug-optima-vs-reaktoro/hop_sweep.cpp --bic 4 0.0), warm-step failures and
 // total iterations, G identical to every printed digit in every row:
 //
 //                              native SIA      HOP          SHP
 //   07PSIna_G_iron             3 fails/ 254    0/ 260       0/  80
 //   07PSIna_G_ironsi           3 fails/ 310    0/ 384       0/ 108
 //   10TH_G_00001               3 fails/4141    0/4128       0/1063
 //   Al-species                 3 fails/ 858    0/ 836       0/ 221
 //   FeNaCl_FyGt_Precip_HighpH  3 fails/ 518    0/ 496       0/ 130
 //   FeNaCl_FyGt_TransitionZone 3 fails/2498    0/2476       0/ 625
 //   CSHSnplus                  0 fails/ 104    0/3316       0/ 841
 //
 // and the cold fallback below fires on NONE of them: the warm native leg
 // genuinely succeeds where plain native SIA on native's own answer fails.
 //
 // ON A SWEEP the saving is real but project-dependent, because there SHP can
 // be no better than the native leg it is built on. One TNode stepped through
 // a temperature sweep, iterations and wall time for the whole sweep:
 //
 //   j_10TH_G_seawater      0-80 C,   81 pts   HOP 14490/3500ms  SHP  3903/1467ms
 //   j_Solvus_G_series1   400-700 C, 301 pts   HOP 84125/1113ms  SHP 77068/1083ms
 //   j_Kaolinite_G_pHtitr  25-125 C,  21 pts   HOP  2961/  19ms  SHP  6485/  45ms
 //
 // i.e. 2.4x less wall time where native's warm start is a win (seawater),
 // near-nothing where a second effect eats it (Solvus - see the pm.pNP note in
 // CalculateEquilibriumStateHOP()), and a LOSS where native's own warm start
 // is itself a loss (Kaolinite, whose native_warm costs 5694 iterations on
 // that sweep against native_cold's 2067 - warm is not always cheaper).
 //
 // WHY IT IS SAFE. Native's SIA is documented to REFUSE states its own cold
 // path returns (plan-v5 section 29.2/60.5: on ten projects the cold path
 // returns a state failing its own mass-balance test, and SIA is the only
 // path that checks). So the warm native leg here carries a COLD FALLBACK:
 // if it throws, the leg is retried cold on the same node, and the mode
 // degrades to exactly HOP. Combined with HOP's own floor (if the Optima leg
 // fails, native's answer is restored and reported as BAD), the chain is
 // SHP >= HOP >= native on every project - the worst case is one wasted, and
 // cheap, warm native attempt.
 //
 // On a node that has never solved, there is nothing to warm-start FROM:
 // detected the same way the Optima path detects it (pm.U[] identically zero -
 // it is the dual, so any real solve by any solver leaves it nonzero) and the
 // native leg then runs cold with a logged warning, rather than starting from
 // the .dbr file's stored speciation, which is not a warm start but a cold
 // start from stale data. Same foot-gun, same detector, same reason - see the
 // corresponding block in ipm_optima.cpp.
 NEED_GEM_SHP = 26,  ///< Need the hybrid native-then-Optima calculation with a WARM (SIA) native leg
 OK_GEM_SHP   = 27,  ///< OK after the hybrid calculation with a warm native leg
 BAD_GEM_SHP  = 28,  ///< Bad (not fully trustful) result after the hybrid calculation with a warm native leg
 ERR_GEM_SHP  = 29   ///< Failure (no result) in the hybrid calculation with a warm native leg
} /*NODECODECH*/;


// \typedef NODECODEFMT Node status codes set by the FMT (FluidMassTransport) part
typedef enum {
 No_nodearray  = -1, ///< Indicates that no node transport properties are present in this DATABR and DBR file
 No_transport  = 0,  ///< Chemical calculations only, no transport coupled
 Initial_RUN   = 1,
 OK_Hydraulic  = 2,
 BAD_Hydraulic = 3,  ///< insufficient convergence
 OK_Transport  = 4,
 BAD_Transport = 5,  ///< insufficient convergence
 NEED_RecalcMT = 6,
 OK_MassBal    = 7,
 OK_RecalcPar  = 8,
 Bad_Recalc    = 9,
 Bad_Time_Step = 10  ///< converged but TKinMet finds that time step must be less (returns in dt)

} NODECODEFMT; /// Node status codes set by the FMT (FluidMassTransport) part

typedef enum {  /// Node type codes controlling hydraulic/mass-transport behavior
  normal       = 0, ///< normal node
// boundary condition node
  NBC1source   = 1, ///< Dirichlet source ( constant concentration )
  NBC1sink    = -1, ///< Dirichlet sink
  NBC2source   = 2, ///< Neumann source ( constant gradient )
  NBC2sink    = -2, ///< Neumann sink
  NBC3source   = 3, ///< Cauchy source ( constant flux )
  NBC3sink    = -3, ///< Cauchy sink
  INIT_FUNK    = 4  ///< functional conditions (e.g. input time-depended functions)
} NODETYPE;

typedef enum {  /// Node type codes controlling Indexation Code
   undefi  = 0,    ///<
   nICbi   = -1,   ///< index of Independent Components kept in the DBR file and DATABR memory structure
   nDCbi   = -2,   ///< index of Dependent Components kept in the DBR file and DATABR memory structure
   nPHbi   = -3,   ///< index of Phases to be kept in the DBR file and DATABR structure
   nPSbi   = -4,   ///< index of Phases-solutions (multicomponent phases) to be kept in the DBR file and DATABR memory structure
   nPSbnICbi = -5  ///<  index of multicomponent phases
 } NODEINEX;

typedef enum {  /// Field index into outField structure
f_NodeHandle = 0,f_NodeTypeHY,f_NodeTypeMT,f_NodeStatusFMT,f_NodeStatusCH,
f_IterDone, f_TK, f_P, f_Vs,f_Vi,
f_Ms, f_Mi, f_Hs, f_Hi, f_Gs,
f_IS, f_pH, f_pe, f_Eh,
f_Tm, f_dt,
//#ifndef NO_NODEARRAYLEVEL
f_Dif,f_Vt, f_vp, f_eps,
f_Km, f_Kf, f_S,  f_Tr, f_h,
f_rho,f_al, f_at, f_av, f_hDl,
f_hDt, f_hDv, f_nto,
//#endif
    // dynamic arrays (52-38=14+2new)
f_bIC, f_rMB, f_uIC, f_xDC, f_gam,
f_dll, f_dul, f_aPH, f_xPH, f_vPS,
f_mPS, f_bPS, f_xPA, f_bSP,
f_amru, f_amrl,       f_omph,
// only for VTK format output
f_mPH, f_vPH, f_m_t, f_con, f_mju, f_lga

} DATABR_FIELDS;


#endif

// -----------------------------------------------------------------------------
// end of DataBr_h

