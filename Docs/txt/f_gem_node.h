// ----------------------------------------------------------------------------------
// f_gem_node.h
// C declaration file for external access to TNode class of GEMIPM2K
// from fortran programs
// Part of the GEMS-PSI program package developed at LES PSI
// (C) 2006 S.Dmytriyeva, D.Kulik, W.Pfingsten
//
// This is a documentation file and, in principle, may not be linked to the 
//  compiled program; only the file f_gem_node.cpp should be linked. 
 
#ifndef _fgemnode_
#define _fgemnode_

#include "node.h"

// (1) Initialization of TNode data exchange interface
// This function reads DATACH and MULTI text input files which it takes from 
// the *.lst path file name string_. 
// Returns 0 if success or 1 if error;
//  must be called once before the beginning of the coupled RMT calculation
//
#ifdef __unix
 extern "C" int f_gem_init_( char* string_, unsigned int length_ );
#else
 extern "C" int  __stdcall  F_GEM_INIT(
    char* string_,     // Path to the *.lst file with names of DATACH and MULTI files
    unsigned int length_  // length of the block string_
   );
#endif


// (2) This function gets some elements of the DATACH data structure
// that can be used in the part of the program external to GEMIPM2K
// returns 0 if success or 1 if error
// Parameter list may be extended in future with other DCH elements
//
#ifdef __unix
   extern "C" int  f_gem_get_dch_( int& p_nICb, int& p_nDCb, int& p_nPHb, float* p_A );
#else
  extern "C" int  __stdcall  F_GEM_GET_DCH(  // All parameters are return values
   int& p_nICb,   // Number of Independent Components (ICs) in chemical system
   int& p_nDCb,   // Number of Dependent Components (DCs) in chemical system
   int& p_nPHb,   // Number of Phases (PHs) in chemical system
   float* p_A     // Stoichiometry matrix A with nDCb rows and nICb columns
                  // before calling  F_GEM_GET_DCH(), a memory block of the size
                  // greater than sizeof(float)*nDCb*nICb must be allocated and
                  // a pointer to is passed in p_A
  );
#endif


// (3) Reads DATABR input file (string_) that describes a single node composition,
//   runs the GEM calculation and provides via parameters the GEM input and
//   output data
//
//  returns 0 if success or 1 if error
//  Usually is called once at the beginning of coupled modeling for initializing
//  all nodes
//
#ifdef __unix
  extern "C" int   f_gem_read_node_( char* string_,
#else
  extern "C" int __stdcall   F_GEM_READ_NODE(
  char* string_,        // path (file name) of the DATABR file
  unsigned int length_, // length of the block string_
#endif
   int& p_NodeHandle,    // Node identification handle
   int& p_NodeTypeHY,    // Node type (hydraulic); see typedef NODETYPE
   int& p_NodeTypeMT,    // Node type (mass transport); see typedef NODETYPE
   int& p_NodeStatusFMT, // Node status code FMT; see typedef NODECODEFMT
   int& p_NodeStatusCH,  // Node status code CH;  see typedef NODECODECH
   int& p_IterDone,      // Number of iterations performed by IPM
// Chemical scalar variables
   double &p_T,     // Temperature T, K                        +      +      -     -
   double &p_P,     // Pressure P, bar                         +      +      -     -
   double &p_Vs,    // Volume V of reactive subsystem, cm3     -      -      +     +
   double &p_Vi,    // Volume of inert subsystem  ?            +      -      -     +
   double &p_Ms,    // Mass of reactive subsystem, kg          +      +      -     -
   double &p_Mi,    // Mass of inert part, kg    ?             +      -      -     +
   double &p_Gs,    // Gibbs energy of reactive subsystem (J)  -      -      +     +
   double &p_Hs,    // Enthalpy of reactive subsystem (J)      -      -      +     +
   double &p_Hi,    // Enthalpy of inert subsystem (J) ?       +      -      -     +
   double &p_IC,    // Effective molal aq ionic strength       -      -      +     +
   double &p_pH,    // pH of aqueous solution                  -      -      +     +
   double &p_pe,    // pe of aqueous solution                  -      -      +     +
   double &p_Eh,    // Eh of aqueous solution, V               -      -      +     +
//  FMT variables (units need dimensionsless form)
//   double &p_Tm,    // actual total simulation time
//   double &p_dt,    // actual time step
//   double &p_dt1,   // previous time step
//   double &p_Vt,    // total volume of node (voxel) = dx*dy*dz, m**3
//   double &p_vp,	// advection velocity (in pores) in this node
//   double &p_eps,   // effective (actual) porosity normalized to 1
//   double &p_Km,    // actual permeability, m**2
//   double &p_Kf,    // actual DARCY`s constant, m**2/s
//   double &p_S,	    // specific storage coefficient, dimensionless
//   double &p_Tr,    // transmissivity m**2/s
//   double &p_h,	    // actual hydraulic head (hydraulic potential), m
//   double &p_rho,   // actual carrier density for density driven flow, g/cm**3
//   double &p_al,    // specific longitudinal dispersivity of porous media, m
//   double &p_at,    // specific transversal dispersivity of porous media, m
//   double &p_av,    // specific vertical dispersivity of porous media, m
//   double &p_hDl,   // hydraulic longitudinal dispersivity, m**2/s
//   double &p_hDt,   // hydraulic transversal dispersivity, m**2/s
//   double &p_hDv,   // hydraulic vertical dispersivity, m**2/s
//  double &p_nto,   // tortuosity factor
// Dynamic data - dimensions see in DATACH.H and DATAMT.H structures
// exchange of values occurs through lists of indices, e.g. xDC, xPH
   double  *p_bIC,  // bulk mole amounts of IC[nICb]                +      +      -     -
   double  *p_rMB,  // MB Residuals from GEM IPM [nICb]             -      -      +     +
   double  *p_uIC,  // IC chemical potentials (mol/mol)[nICb]       -      -      +     +
   double  *p_xDC,    // DC mole amounts at equilibrium [nDCb]      -      -      +     +
   double  *p_gam,    // activity coeffs of DC [nDCb]               -      -      +     +
   double  *p_dul,  // upper kinetic restrictions [nDCb]            +      +      -     -
   double  *p_dll,  // lower kinetic restrictions [nDCb]            +      +      -     -
   double  *p_aPH,  // Specific surface areas of phases (m2/g)      +      +      -     -
   double  *p_xPH,  // total mole amounts of phases [nPHb]          -      -      +     +
   double  *p_vPS,  // phase volume, cm3/mol        [nPSb]          -      -      +     +
   double  *p_mPS,  // phase (carrier) mass, g      [nPSb]          -      -      +     +
   double  *p_bPS,  // bulk compositions of phases  [nPSb][nICb]    -      -      +     +
   double  *p_xPA,  // amount of carrier in phases  [nPSb] ??       -      -      +     +
  // What else?
//   double  *p_dRes1
#ifdef __unix
   , unsigned int length_
#endif
 );


// (4) Provides GEM input data to GEMIPM2K kernel (by copying them from
//   parameters to DATABR work structure);
//   runs the GEM-IPM calculation;
//   and provides via parameters
//   the GEM output data back to the calling program
//
//  returns  0 if success or 1 if error
//  Is called on each external iteration for each node
//
#ifdef __unix
   extern "C" int  f_gem_calc_node_(
#else
   extern "C" int  __stdcall   F_GEM_CALC_NODE(
#endif
   int& p_NodeHandle,    // Node identification handle
   int& p_NodeTypeHY,    // Node type (hydraulic); see typedef NODETYPE
   int& p_NodeTypeMT,    // Node type (mass transport); see typedef NODETYPE
   int& p_NodeStatusFMT, // Node status code FMT; see typedef NODECODEFMT
   int& p_NodeStatusCH,  // Node status code CH;  see typedef NODECODECH
   int& p_IterDone,      // Number of iterations performed by IPM
// Chemical scalar variables
   double &p_T,     // Temperature T, K                        +      +      -     -
   double &p_P,     // Pressure P, bar                         +      +      -     -
   double &p_Vs,    // Volume V of reactive subsystem, cm3     -      -      +     +
   double &p_Vi,    // Volume of inert subsystem  ?            +      -      -     +
   double &p_Ms,    // Mass of reactive subsystem, kg          +      +      -     -
   double &p_Mi,    // Mass of inert part, kg    ?             +      -      -     +
   double &p_Gs,    // Gibbs energy of reactive subsystem (J)  -      -      +     +
   double &p_Hs,    // Enthalpy of reactive subsystem (J)      -      -      +     +
   double &p_Hi,    // Enthalpy of inert subsystem (J) ?       +      -      -     +
   double &p_IC,    // Effective molal aq ionic strength       -      -      +     +
   double &p_pH,    // pH of aqueous solution                  -      -      +     +
   double &p_pe,    // pe of aqueous solution                  -      -      +     +
   double &p_Eh,    // Eh of aqueous solution, V               -      -      +     +
//  FMT variables (units need dimensionsless form)
//   double &p_Tm,    // actual total simulation time
//   double &p_dt,    // actual time step
//   double &p_dt1,   // previous time step
//   double &p_Vt,    // total volume of node (voxel) = dx*dy*dz, m**3
//   double &p_vp,	// advection velocity (in pores) in this node
//   double &p_eps,   // effective (actual) porosity normalized to 1
//   double &p_Km,    // actual permeability, m**2
//   double &p_Kf,    // actual DARCY`s constant, m**2/s
//   double &p_S,	    // specific storage coefficient, dimensionless
//   double &p_Tr,    // transmissivity m**2/s
//   double &p_h,	    // actual hydraulic head (hydraulic potential), m
//   double &p_rho,   // actual carrier density for density driven flow, g/cm**3
//   double &p_al,    // specific longitudinal dispersivity of porous media, m
//   double &p_at,    // specific transversal dispersivity of porous media, m
//   double &p_av,    // specific vertical dispersivity of porous media, m
//   double &p_hDl,   // hydraulic longitudinal dispersivity, m**2/s
//   double &p_hDt,   // hydraulic transversal dispersivity, m**2/s
//   double &p_hDv,   // hydraulic vertical dispersivity, m**2/s
//   double &p_nto,   // tortuosity factor
// Dynamic data - dimensions see in DATACH.H and DATAMT.H structures
// exchange of values occurs through lists of indices, e.g. xDC, xPH
   double  *p_bIC,  // bulk mole amounts of IC[nICb]                +      +      -     -
   double  *p_rMB,  // MB Residuals from GEM IPM [nICb]             -      -      +     +
   double  *p_uIC,  // IC chemical potentials (mol/mol)[nICb]       -      -      +     +
   double  *p_xDC,    // DC mole amounts at equilibrium [nDCb]      -      -      +     +
   double  *p_gam,    // activity coeffs of DC [nDCb]               -      -      +     +
   double  *p_dul,  // upper kinetic restrictions [nDCb]            +      +      -     -
   double  *p_dll,  // lower kinetic restrictions [nDCb]            +      +      -     -
   double  *p_aPH,  // Specific surface areas of phases (m2/g)      +      +      -     -
   double  *p_xPH,  // total mole amounts of phases [nPHb]          -      -      +     +
   double  *p_vPS,  // phase volume, cm3/mol        [nPSb]          -      -      +     +
   double  *p_mPS,  // phase (carrier) mass, g      [nPSb]          -      -      +     +
   double  *p_bPS,  // bulk compositions of phases  [nPSb][nICb]    -      -      +     +
   double  *p_xPA  // amount of carrier in phases  [nPSb] ??        -      -      +     +
  // What else?
  // double  *p_dRes1
 );


// (5) Writes a DATABR text file (with file path name provided in string_)
//    from the currently available content of the DATABR work structure.
// can be used for interruption of coupled modeling or for debugging purposes
//
#ifdef __unix
 extern "C" void f_gem_write_node_( char* string_, unsigned int length_ );
#else
 extern "C"  void  __stdcall   F_GEM_WRITE_NODE(
   char* string_,        // path (file name) of the DATABR file
   unsigned int length_, // length of the block fname_
);
#endif

// (6) For detailed examination of GEM work data structure:
// writes GEMIPM internal MULTI data structure into text file
// path name string_ in debugging format (different from MULTI input format).
// This file cannot be read back with F_GEM_INIT()!
//
#ifdef __unix
 extern "C" void f_gem_print_ipm_( char* string_, unsigned int length_ );
#else
 extern "C"  void  __stdcall   F_GEM_PRINT_IPM(
   char* string_,        // path (file name) of the DATABR file
   unsigned int length_, // length of the block fname_
);
#endif


#endif
// end of f_gem_node.h
// -------------------------------------------------------------------------------------
