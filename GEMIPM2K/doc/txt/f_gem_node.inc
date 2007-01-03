c     // GEMIPM2K v. 0.98 (prototype) 12.12.2006
c 
c     // masstransport - GEMIPM2K TNode level coupling
c     // Interface descriptions that can be used in Digital fortran-90
c     // for C++ part, see f_gem_node.h file
c
c     // This is a documentation file as a fortran counterpart to the f_gem_node.h 
c     // C/C++ file 
c
c     // (1) Initialization of TNode data exchange interface
c     // This function reads DATACH and MULTI text input files which it takes from 
c     // the *.lst path file name string_. 
c     // Returns 0 if success or 1 if error;
c     //  must be called once before the beginning of the coupled RMT calculation
c     //
	INTERFACE TO SUBROUTINE F_GEM_INIT( string )
	 CHARACTER* (*) string [REFERENCE]
	END
c
c
c     // (2) This function gets some elements of the DATACH data structure
c     // that can be used in the part of the program external to GEMIPM2K
c     // Parameter list may be extended in future with other DCH elements
c     //
	INTERFACE
	subroutine  F_GEM_GET_DCH( p_nICb, p_nDCb, p_nPHb, p_A )
	 integer p_nICb [REFERENCE]   
         integer p_nDCb [REFERENCE]   
         integer p_nPHb [REFERENCE]
         real*4    p_A( p_nDCb, p_nICb )
	END SUBROUTINE F_GEM_GET_DCH
	END INTERFACE


c     // (3) Reads DATABR input file that describes a single node composition,
c     //   runs the GEM calculation and provides via parameters the GEM input and
c     //   output data
c     //  returns 0 if success or 1 if error
c     //  Usually is called once at the beginning of coupled modeling for initializing
c     //  all nodes
c     //
      INTERFACE 
      SUBROUTINE  F_GEM_READ_NODE( p_file, p_NodeHandle
     *,p_NodeTypeHY,p_NodeTypeMT,p_NodeStatusFMT
     *,p_NodeStatusCH,p_IterDone
     *,p_T, p_P, p_Vs,p_Vi,p_Ms,p_Mi,p_Gs,p_Hs,p_Hi,p_IC,p_pH,p_pe,p_Eh
c     *,p_Tm,p_dt,p_dt1,p_Vt,p_vp,p_eps,p_Km,p_Kf,p_S,p_Tr,p_h,p_rho
c     *,p_al,p_at,p_av,p_hDl,p_hDt,p_hDv,p_nto
     *,p_bIC,p_rMB,p_uIC, p_xDC,p_gam,p_dul,p_dll, p_aPH
     *,p_xPH,p_vPS,p_mPS,p_bPS,p_xPA
     *) 
      CHARACTER* (*) p_file [REFERENCE]  !// File path name 
      integer p_NodeHandle [REFERENCE]   !// Node identification handle
      integer p_NodeTypeHY [REFERENCE]   !// Node type (hydraulic); see typedef NODETYPE
      integer p_NodeTypeMT [REFERENCE]   !// Node type (mass transport); see typedef NODETYPE
      integer p_NodeStatusFMT [REFERENCE]!// Node status code FMT; see typedef NODECODEFMT
      integer p_NodeStatusCH  [REFERENCE]!// Node status code CH;  see typedef NODECODECH
      integer p_IterDone      [REFERENCE]!// Number of iterations performed by IPM (if not need GEM)
c !// Chemical scalar variables
      real*8 p_T    [REFERENCE]  !// Temperature T, K                        +      +      -     -
      real*8 p_P    [REFERENCE]  !// Pressure P, bar                         +      +      -     -
      real*8 p_Vs   [REFERENCE]  !// Volume V of reactive subsystem, cm3     -      -      +     +
      real*8 p_Vi   [REFERENCE]  !// Volume of inert subsystem  ?            +      -      -     +
      real*8 p_Ms   [REFERENCE]  !// Mass of reactive subsystem, kg          +      +      -     -
      real*8 p_Mi   [REFERENCE]  !// Mass of inert part, kg    ?             +      -      -     +
      real*8 p_Gs   [REFERENCE]  !// Gibbs energy of reactive subsystem (J)  -      -      +     +
      real*8 p_Hs   [REFERENCE]  !// Enthalpy of reactive subsystem (J)      -      -      +     +
      real*8 p_Hi   [REFERENCE]  !// Enthalpy of inert subsystem (J) ?       +      -      -     +
      real*8 p_IC   [REFERENCE]  !// Effective molal aq ionic strength           -      -      +     +
      real*8 p_pH   [REFERENCE]  !// pH of aqueous solution                      -      -      +     +
      real*8 p_pe   [REFERENCE]  !// pe of aqueous solution                      -      -      +     +
      real*8 p_Eh   [REFERENCE]  !// Eh of aqueous solution, V                   -      -      +     +
c !//  FMT variables (units need dimensionsless form)
c      real*8 p_Tm  [REFERENCE]  !// actual total simulation time
c      real*8 p_dt  [REFERENCE]  !// actual time step
c      real*8 p_dt1 [REFERENCE]  !// priveous time step
c      real*8 p_Vt  [REFERENCE]  !// total volume of node (voxel) = dx*dy*dz, m**3
c      real*8 p_vp  [REFERENCE]  !// advection velocity (in pores) in this node
c      real*8 p_eps [REFERENCE]  !// effective (actual) porosity normalized to 1
c      real*8 p_Km  [REFERENCE]  !// actual permeability, m**2
c      real*8 p_Kf  [REFERENCE]  !// actual DARCY`s constant, m**2/s
c      real*8 p_S   [REFERENCE]  !// specific storage coefficient, dimensionless
c      real*8 p_Tr  [REFERENCE]  !// transmissivity m**2/s
c      real*8 p_h   [REFERENCE]  !// actual hydraulic head (hydraulic potential), m
c      real*8 p_rho [REFERENCE]  !// actual carrier density for density driven flow, g/cm**3
c      real*8 p_al  [REFERENCE]  !// specific longitudinal dispersivity of porous media, m
c      real*8 p_at  [REFERENCE]  !// specific transversal dispersivity of porous media, m
c      real*8 p_av  [REFERENCE]  !// specific vertical dispersivity of porous media, m
c      real*8 p_hDl [REFERENCE]  !// hydraulic longitudinal dispersivity, m**2/s, diffusities from chemical database
c      real*8 p_hDt [REFERENCE]  !// hydraulic transversal dispersivity, m**2/s
c      real*8 p_hDv [REFERENCE]  !// hydraulic vertical dispersivity, m**2/s
c      real*8 p_nto [REFERENCE]  !// tortuosity factor
c !// Dynamic data - dimensions see in DATACH.H and DATAMT.H structures
c !// exchange of values occurs through lists of indices, e.g. xDC, xPH
      real*8  p_bIC(7) ! (nICb)  !// bulk mole amounts of IC[nICb]                +      +      -     -
      real*8  p_rMB(7) ! (nICb)  !// MB Residuals from GEM IPM [nICb]             -      -      +     +
      real*8  p_uIC(7) ! (nICb)  !// IC chemical potentials (mol/mol)[nICb]       -      -      +     +
      real*8  p_xDC(17) ! (nDCb)  !// DC mole amounts at equilibrium [nDCb]      -      -      +     +
      real*8  p_gam(17) ! (nDCb)  !// activity coeffs of DC [nDCb]               -      -      +     +
      real*8  p_dul(17) ! (nDCb)  !// upper kinetic restrictions [nDCb]           +      +      -     -
      real*8  p_dll(17) ! (nDCb)  !//  lower kinetic restrictions [nDCb]           +      +      -     -
      real*8  p_aPH(20) ! (nPHb)  !// Specific surface areas of phases (m2/g)       +      +      -     -
      real*8  p_xPH(20) ! (nPHb)  !// total mole amounts of phases [nPHb]          -      -      +     +
      real*8  p_vPS(20) ! (nPSb)  !// phase volume, cm3/mol        [nPSb]          -      -      +     +
      real*8  p_mPS(20) ! (nPSb)  !// phase (carrier) mass, g      [nPSb]          -      -      +     +
      real*8  p_bPS(7,20) ! (nICBb,nPSb)  !// bulk compositions of phases  [nPSb][nICb]    -      -      +     +
      real*8  p_xPA(20) ! (nPSb)  !// amount of carrier in phases  [nPSb] ??       -      -      +     +
c !  // What else?
c      real*8  p_dRes1 (7) ! (nICb) ! //should be a reserved array?
      END SUBROUTINE F_GEM_READ_NODE
      END INTERFACE

c     // (4) Provides GEM input data to GEMIPM2K kernel (by copying them from
c     //   parameters to DATABR work structure);
c     //   runs the GEM-IPM calculation;
c     //   and provides via parameters
c     //   the GEM output data back to the calling program
c     //
c     //  returns  0 if success or 1 if error
c     //  Is called on each external iteration for each node
c     //
	INTERFACE
	subroutine  F_GEM_CALC_NODE( p_NodeHandle,p_NodeTypeHY,p_NodeTypeMT
     *,p_NodeStatusFMT,p_NodeStatusCH,p_IterDone
     *,p_T, p_P, p_Vs,p_Vi,p_Ms,p_Mi,p_Gs
     *,p_Hs,p_Hi,p_IC,p_pH,p_pe,p_Eh
c     *,p_Tm,p_dt,p_dt1,p_Vt,p_vp,p_eps,p_Km,p_Kf,p_S,p_Tr,p_h,p_rho
c     *,p_al,p_at,p_av,p_hDl,p_hDt,p_hDv,p_nto
     *,p_bIC,p_rMB,p_uIC, p_xDC,p_gam,p_dul,p_dll, p_aPH
     *,p_xPH,p_vPS,p_mPS,p_bPS,p_xPA 
     *)
      integer p_NodeHandle [REFERENCE]   !// Node identification handle
      integer p_NodeTypeHY [REFERENCE]   !// Node type (hydraulic); see typedef NODETYPE
      integer p_NodeTypeMT [REFERENCE]   !// Node type (mass transport); see typedef NODETYPE
      integer p_NodeStatusFMT [REFERENCE]!// Node status code FMT; see typedef NODECODEFMT
      integer p_NodeStatusCH  [REFERENCE]!// Node status code CH;  see typedef NODECODECH
      integer p_IterDone      [REFERENCE]!// Number of iterations performed by IPM (if not need GEM)
c !// Chemical scalar variables
      real*8 p_T [REFERENCE]    !// Temperature T, K                        +      +      -     -
      real*8 p_P [REFERENCE]    !// Pressure P, bar                         +      +      -     -
      real*8 p_Vs [REFERENCE]   !// Volume V of reactive subsystem, cm3     -      -      +     +
      real*8 p_Vi  [REFERENCE]  !// Volume of inert subsystem  ?            +      -      -     +
      real*8 p_Ms  [REFERENCE]  !// Mass of reactive subsystem, kg          +      +      -     -
      real*8 p_Mi  [REFERENCE]  !// Mass of inert part, kg    ?             +      -      -     +
      real*8 p_Gs  [REFERENCE]  !// Gibbs energy of reactive subsystem (J)  -      -      +     +
      real*8 p_Hs  [REFERENCE]  !// Enthalpy of reactive subsystem (J)      -      -      +     +
      real*8 p_Hi  [REFERENCE]  !// Enthalpy of inert subsystem (J) ?       +      -      -     +
      real*8 p_IC   [REFERENCE] !// Effective molal aq ionic strength       -      -      +     +
      real*8 p_pH   [REFERENCE] !// pH of aqueous solution                  -      -      +     +
      real*8 p_pe  [REFERENCE]  !// pe of aqueous solution                  -      -      +     +
      real*8 p_Eh  [REFERENCE]  !// Eh of aqueous solution, V               -      -      +     +
c !//  FMT variables (units need dimensionsless form)
c      real*8 p_Tm [REFERENCE]   !// actual total simulation time
c      real*8 p_dt  [REFERENCE]  !// actual time step
c      real*8 p_dt1 [REFERENCE]  !// priveous time step
c      real*8 p_Vt  [REFERENCE]  !// total volume of node (voxel) = dx*dy*dz, m**3
c      real*8 p_vp  [REFERENCE]   !// advection velocity (in pores) in this node
c      real*8 p_eps [REFERENCE]  !// effective (actual) porosity normalized to 1
c      real*8 p_Km  [REFERENCE]  !// actual permeability, m**2
c      real*8 p_Kf [REFERENCE]   !// actual DARCY`s constant, m**2/s
c      real*8 p_S	[REFERENCE]  !// specific storage coefficient, dimensionless
c      real*8 p_Tr  [REFERENCE]  !// transmissivity m**2/s
c      real*8 p_h	[REFERENCE]  !  // actual hydraulic head (hydraulic potential), m
c      real*8 p_rho [REFERENCE]  !// actual carrier density for density driven flow, g/cm**3
c      real*8 p_al  [REFERENCE]  !// specific longitudinal dispersivity of porous media, m
c      real*8 p_at [REFERENCE]   !// specific transversal dispersivity of porous media, m
c      real*8 p_av  [REFERENCE]  !// specific vertical dispersivity of porous media, m
c      real*8 p_hDl [REFERENCE]  !// hydraulic longitudinal dispersivity, m**2/s, diffusities from chemical database
c      real*8 p_hDt [REFERENCE]  !// hydraulic transversal dispersivity, m**2/s
c      real*8 p_hDv [REFERENCE]  !// hydraulic vertical dispersivity, m**2/s
c      real*8 p_nto [REFERENCE]  !// tortuosity factor
c !// Dynamic data - dimensions see in DATACH.H and DATAMT.H structures
c !// exchange of values occurs through lists of indices, e.g. xDC, xPH
      real*8  p_bIC(7) ! (nICb)  !// bulk mole amounts of IC[nICb]                +      +      -     -
      real*8  p_rMB(7) ! (nICb)  !// MB Residuals from GEM IPM [nICb]             -      -      +     +
      real*8  p_uIC(7) ! (nICb)  !// IC chemical potentials (mol/mol)[nICb]       -      -      +     +
      real*8  p_xDC(17) ! (nDCb)  !  // DC mole amounts at equilibrium [nDCb]      -      -      +     +
      real*8  p_gam(17) ! (nDCb)  !  // activity coeffs of DC [nDCb]               -      -      +     +
      real*8  p_dul(17) ! (nDCb)  ! // upper kinetic restrictions [nDCb]           +      +      -     -
      real*8  p_dll(17) ! (nDCb)  ! //  lower kinetic restrictions [nDCb]           +      +      -     -
      real*8  p_aPH(20) ! (nPHb)  !// Specific surface areas of phases (m2/g)       +      +      -     -
      real*8  p_xPH(20) ! (nPHb)  !// total mole amounts of phases [nPHb]          -      -      +     +
      real*8  p_vPS(20) ! (nPSb)  !// phase volume, cm3/mol        [nPSb]          -      -      +     +
      real*8  p_mPS(20) ! (nPSb)  !// phase (carrier) mass, g      [nPSb]          -      -      +     +
      real*8  p_bPS(7,20) ! (nICBb,nPSb)  !// bulk compositions of phases  [nPSb][nICb]    -      -      +     +
      real*8  p_xPA(20) ! (nPSb)  !// amount of carrier in phases  [nPSb] ??       -      -      +     +
c !  // What else?
c      real*8  p_dRes1 (7) ! (nICb) ! //should be array Sveta ?
	END SUBROUTINE F_GEM_CALC_NODE
	END INTERFACE

c    // (5) Writes a DATABR text file (with file path name provided in string_)
c    //    from the currently available content of the DATABR work structure.
c    // can be used for interruption of coupled modeling or for debugging purposes
c    //
	INTERFACE TO SUBROUTINE F_GEM_WRITE_NODE( string )
	 CHARACTER* (*) string [REFERENCE]
	END
c
c    // (6) For detailed examination of GEM work data structure:
c    // writes GEMIPM internal MULTI data structure into text file
c    // path name fname in debugging format (different from MULTI input format).
c    // This file cannot be read back with F_GEM_INIT()!
c    //
	INTERFACE TO SUBROUTINE F_GEM_PRINT_IPM( string )
	 CHARACTER* (*) string [REFERENCE]
	END
c
c     // End of f_gem_node.inc file
