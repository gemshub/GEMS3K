//-------------------------------------------------------------------
// $Id: databr.h 1430 2009-08-27 17:02:11Z gems $
//
// DataBRidge - defines the structure of node-dependent data for
// exchange between the coupled GEM IPM and FMT code parts.
// Requires DATACH.H header and data structure.
// DATABR is used in TNode and TNodeArray classes.
//
//      CH: chemical structure in GEM IPM
//      FMT: fluid mass transport
//
// Copyright (C) 2003,2009 by D.Kulik, S.Dmytriyeva, F.Enzmann, W.Pfingsten
//
// This file is part of GEMIPM2K and GEMS-PSI codes for
// thermodynamic modelling by Gibbs energy minimization
// developed by the Laboratory of Waste Management, Paul Scherrer Institute

// This file may be distributed together with GEMIPM2K source code
// under the licence terms defined in GEMIPM2K.QAL
//
// See also http://gems.web.psi.ch/
// E-mail: gems2.support@psi.ch
//-------------------------------------------------------------------
//
#ifndef _DataBr_H_
#define _DataBr_H_

typedef struct  // DATABR - template node data bridge structure
{
   long int
     NodeHandle,    // Node identification handle, not used in calculaions on TNode level
     NodeTypeHY,    // Node type code (hydraulic), not used on TNode level; see typedef NODETYPE
     NodeTypeMT,    // Node type (mass transport), not used on TNode level; see typedef NODETYPE
     NodeStatusFMT, // Node status code in FMT part, not used on TNode level; see typedef NODECODEFMT
     NodeStatusCH,  // Node status code in GEM (input and output);  see typedef NODECODECH
     IterDone;      // Number of iterations performed by GEM IPM in the last run

/*  these important data array dimensions are provided in the DATACH structure
   long int
    nICb,       // Number of Independent Components kept in the DATABR memory structure (<= nIC)
    nDCb,      	// Number of Dependent Components kept in the DATABR memory structure (<=nDC)
    nPHb,     	// Number of Phases to be kept in the DATABR structure (<= nPH)
    nPSb,       // Number of Phases-solutions (multicomponent phases) to be kept in the DATABR memory structure (<= nPS)
*/
//      Usage of this variable (DB - data bridge)      	               MT-DB DB-GEM GEM-DB DB-MT
   double
// Chemical scalar variables
    TK,     // Node temperature T (Kelvin)                     	          +      +      -     -
    P, 	    // Node Pressure P (Pa)                         	          +      +      -     -
    Vs,     // Volume V of reactive subsystem  (m3)                      (+)    (+)     +     +
    Vi,     // Volume of inert subsystem (m3)          	                  +      -      -     +
    Ms,     // Mass of reactive subsystem (kg)         	                  +     (+)     -     -
    Mi,     // Mass of inert subsystem (kg)             	              +      -      -     +

    Gs,     // Total Gibbs energy of the reactive subsystem (J/RT) (norm) -      -      +     +
    Hs, 	// Total enthalpy of reactive subsystem (J) (reserved)        -      -      +     +
    Hi,     // Total enthalpy of inert subsystem (J) (reserved)           +      -      -     +

    IC,     // Effective aqueous ionic strength (molal)                    -      -      +     +
    pH,     // pH of aqueous solution in the activity scale (-log10 molal) -      -      +     +
    pe,     // pe of aqueous solution in the activity scale (-log10 molal) -      -      +     +
    Eh,     // Eh of aqueous solution (V)                                  -      -      +     +

//  FMT variables (units or dimensionsless) - to be used for storing them
//  at the nodearray level, normally not used in the single-node FMT-GEM coupling
    Tm,     // Actual total simulation time (s)
    dt,     // Actual time step (s)
    Dif,    // General diffusivity of disolved matter (m2/s)
    Vt,		// Total volume of the node (m3) (Vs + Vi)
    vp,		// Advection velocity (in pores)  (m/s)
    eps,	// Effective (actual) porosity normalized to 1
    Km,		// Actual permeability (m2)
    Kf,		// Actual Darcy`s constant (m2/s)
    S,		// Specific storage coefficient, dimensionless (default 1.0)
            //  if not 1.0 then can be used as mass scaling factor relative to Ms for the
            //  bulk composition/speciation of reactive sub-system in this node
    Tr,     // transmissivity m2/s
    h,		// Actual hydraulic head (hydraulic potential) (m)
    rho,	// Actual carrier density for density-driven flow (kg/m3)
    al,		// Specific longitudinal dispersivity of porous media (m)
    at,		// Specific transversal dispersivity of porous media (m)
    av,		// Specific vertical dispersivity of porous media (m)
    hDl,	// Hydraulic longitudinal dispersivity (m2/s)
    hDt,	// Hydraulic transversal dispersivity (m2/s)
    hDv,	// Hydraulic vertical dispersivity (m2/s)
    nto;	// Tortuosity factor (dimensionless)

// Data arrays - dimensions nICb, nDCb, nPHb, nPSb see in the DATACH structure
// exchange of values occurs through lists of indices, e.g. xIC, xDC, xPH from DATACH

//      Usage of this variable (DB = data bridge)                                               MT-DB DB-GEM GEM-DB DB-MT
   double
// IC (stoichiometry units)
    *bIC,  // Bulk composition of (reactive part of) the system (moles)[nICb]                     +      +      -     -
    *rMB,  // Mass balance residuals (moles) [nICb]                                               -      -      +     +
    *uIC,  // Chemical potentials of ICs (dual GEM solution) normalized scale (mol/mol)[nICb]     -      -      +     +
// DC (species) in reactive subsystem
    *xDC,  // Speciation - amounts of DCs in equilibrium state - primal GEM solution(moles)[nDCb] (+)    (+)     +     +
    *gam,  // Activity coefficients of DCs in their respective phases  [nDCb]                     (+)    (+)     +     +
// Metastability/kinetic controls
    *dul,  // Upper metastability constraints on amounts of DCs (moles) [nDCb]                     +      +      -     -
    *dll,  // Lower metastability constraints on amounts of DCs (moles)  [nDCb]                    +      +      -     -
// Phases in reactive subsystem
    *aPH,  // Specific surface areas of phases (m2/kg) [nPHb]                                      +      +      -     -
    *xPH,  // Amounts of phases in equilibrium state (moles) [nPHb]                                -      -      +     +
    *vPS,  // Volumes of multicomponent phases (m3)  [nPSb]                                        -      -      +     +
    *mPS,  // Masses of multicomponent phases (kg)    [nPSb]                                       -      -      +     +
    *bPS,  // Bulk elemental compositions of multicomponent phases (moles) [nPSb][nICb]            -      -      +     +
    *xPA,  // Amount of carrier (sorbent or solvent) in multicomponent phases [nPSb]               -      -      +     +
// (+) can be used as input in "smart initial approximation" mode of GEM IPM-2 algorithm

    // Reserved data array pointer
    *dRes1;
}
DATABR;

typedef DATABR*  DATABRPTR;

typedef enum {  // NodeStatus codes with respect to GEMIPM calculations
 NO_GEM_SOLVER= 0,   // No GEM re-calculation needed for this node
 NEED_GEM_AIA = 1,   // Need GEM calculation with simplex (automatic) initial approximation (IA)
 OK_GEM_AIA   = 2,   // OK after GEM calculation with simplex (automatic) IA
 BAD_GEM_AIA  = 3,   // Bad (not fully trustful) result after GEM calculation with simplex (automatic) IA
 ERR_GEM_AIA  = 4,   // Failure (no result) in GEM calculation with simplex (automatic) IA
 NEED_GEM_SIA = 5,   // Need GEM calculation with no-simplex (smart) IA
                     //   the previous GEM solution (full DATABR lists only)
 OK_GEM_SIA   = 6,   // OK after GEM calculation with no-simplex (smart) IA
 BAD_GEM_SIA  = 7,   // Bad (not fully trustful) result after GEM calculation with no-simplex (smart) IA
 ERR_GEM_SIA  = 8,   // Failure (no result) in GEM calculation with no-simplex (smart) IA
 T_ERROR_GEM  = 9    // Terminal error has occurred in GEMIPM2K (e.g. memory corruption). Restart is required.
} NODECODECH;

typedef enum {  // Node status codes set by the FMT (FluidMassTransport) part
 No_transport  = 0,  // Chemical calculations only, no transport coupled
 Initial_RUN   = 1,
 OK_Hydraulic  = 2,
 BAD_Hydraulic = 3,  // insufficient convergence
 OK_Transport  = 4,
 BAD_Transport = 5,  // insufficient convergence
 NEED_RecalcMT = 6,
 OK_MassBal    = 7,
 OK_RecalcPar  = 8,
 Bad_Recalc    = 9,
 OK_Time       = 10
} NODECODEFMT;

typedef enum {  // Node type codes controlling hydraulic/mass-transport behavior
  normal       = 0, // normal node
// boundary condition node
  NBC1source   = 1, // Dirichlet source ( constant concentration )
  NBC1sink    = -1, // Dirichlet sink
  NBC2source   = 2, // Neumann source ( constant gradient )
  NBC2sink    = -2, // Neumann sink
  NBC3source   = 3, // Cauchy source ( constant flux )
  NBC3sink    = -3, // Cauchy sink
  INIT_FUNK    = 4  // functional conditions (e.g. input time-depended functions)
} NODETYPE;

#endif
// -----------------------------------------------------------------------------
// end of _DataBr_h

