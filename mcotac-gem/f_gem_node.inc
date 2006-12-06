c     //masstransport
	INTERFACE TO SUBROUTINE MAIF_START( str_text )
	 CHARACTER* (*) str_text [REFERENCE]
	END

	INTERFACE
	subroutine  MAIF_GET_A( p_nICb, p_nDCb, p_A )
	integer p_nICb [REFERENCE]   
      integer p_nDCb [REFERENCE]   
      real*4    p_A( 17, 7)
	END SUBROUTINE MAIF_GET_A
	END INTERFACE

	INTERFACE
	subroutine  MAIF_CALC( p_NodeHandle,p_NodeTypeHY,p_NodeTypeMT
     *,p_NodeStatusFMT,p_NodeStatusCH,p_IterDone
     *,p_T, p_P, p_Vs,p_Vi,p_Ms,p_Mi,p_Gs
     *,p_Hs,p_Hi,p_IC,p_pH,p_pe,p_Eh
     *,p_Tm,p_dt,p_dt1,p_Vt,p_vp,p_eps,p_Km,p_Kf,p_S,p_Tr,p_h,p_rho
     *,p_al,p_at,p_av,p_hDl,p_hDt,p_hDv,p_nto
     *,p_bIC,p_rMB,p_uIC, p_xDC,p_gam,p_dul,p_dll, p_aPH
     *,p_xPH,p_vPS,p_mPS,p_bPS,p_xPA,p_dRes1
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
      real*8 p_Tm [REFERENCE]   !// actual total simulation time
      real*8 p_dt  [REFERENCE]  !// actual time step
      real*8 p_dt1 [REFERENCE]  !// priveous time step
      real*8 p_Vt  [REFERENCE]  !// total volume of node (voxel) = dx*dy*dz, m**3
      real*8 p_vp  [REFERENCE]   !// advection velocity (in pores) in this node
      real*8 p_eps [REFERENCE]  !// effective (actual) porosity normalized to 1
      real*8 p_Km  [REFERENCE]  !// actual permeability, m**2
      real*8 p_Kf [REFERENCE]   !// actual DARCY`s constant, m**2/s
      real*8 p_S	[REFERENCE]  !// specific storage coefficient, dimensionless
      real*8 p_Tr  [REFERENCE]  !// transmissivity m**2/s
      real*8 p_h	[REFERENCE]  !  // actual hydraulic head (hydraulic potential), m
      real*8 p_rho [REFERENCE]  !// actual carrier density for density driven flow, g/cm**3
      real*8 p_al  [REFERENCE]  !// specific longitudinal dispersivity of porous media, m
      real*8 p_at [REFERENCE]   !// specific transversal dispersivity of porous media, m
      real*8 p_av  [REFERENCE]  !// specific vertical dispersivity of porous media, m
      real*8 p_hDl [REFERENCE]  !// hydraulic longitudinal dispersivity, m**2/s, diffusities from chemical database
      real*8 p_hDt [REFERENCE]  !// hydraulic transversal dispersivity, m**2/s
      real*8 p_hDv [REFERENCE]  !// hydraulic vertical dispersivity, m**2/s
      real*8 p_nto [REFERENCE]  !// tortuosity factor
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
      real*8  p_dRes1 (7) ! (nICb) ! //should be array Sveta ?
	END SUBROUTINE MAIF_CALC
	END INTERFACE


	INTERFACE TO SUBROUTINE  MAIF_READ( p_file, p_NodeHandle
     *,p_NodeTypeHY,p_NodeTypeMT,p_NodeStatusFMT
     *,p_NodeStatusCH,p_IterDone
     *,p_T, p_P, p_Vs,p_Vi,p_Ms,p_Mi,p_Gs,p_Hs,p_Hi,p_IC,p_pH,p_pe,p_Eh
     *,p_Tm,p_dt,p_dt1,p_Vt,p_vp,p_eps,p_Km,p_Kf,p_S,p_Tr,p_h,p_rho
     *,p_al,p_at,p_av,p_hDl,p_hDt,p_hDv,p_nto
     *,p_bIC,p_rMB,p_uIC, p_xDC,p_gam,p_dul,p_dll, p_aPH
     *,p_xPH,p_vPS,p_mPS,p_bPS,p_xPA,p_dRes1) 
      CHARACTER* (*) p_file [REFERENCE]
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
      real*8 p_IC   [REFERENCE] !// Effective molal aq ionic strength           -      -      +     +
      real*8 p_pH   [REFERENCE] !// pH of aqueous solution                      -      -      +     +
      real*8 p_pe  [REFERENCE]  !// pe of aqueous solution                      -      -      +     +
      real*8 p_Eh  [REFERENCE]  !// Eh of aqueous solution, V                   -      -      +     +
c !//  FMT variables (units need dimensionsless form)
      real*8 p_Tm [REFERENCE]   !// actual total simulation time
      real*8 p_dt  [REFERENCE]  !// actual time step
      real*8 p_dt1 [REFERENCE]  !// priveous time step
      real*8 p_Vt  [REFERENCE]  !// total volume of node (voxel) = dx*dy*dz, m**3
      real*8 p_vp  [REFERENCE]   !// advection velocity (in pores) in this node
      real*8 p_eps [REFERENCE]  !// effective (actual) porosity normalized to 1
      real*8 p_Km  [REFERENCE]  !// actual permeability, m**2
      real*8 p_Kf [REFERENCE]   !// actual DARCY`s constant, m**2/s
      real*8 p_S	[REFERENCE]  !// specific storage coefficient, dimensionless
      real*8 p_Tr  [REFERENCE]  !// transmissivity m**2/s
      real*8 p_h	[REFERENCE]  !  // actual hydraulic head (hydraulic potential), m
      real*8 p_rho [REFERENCE]  !// actual carrier density for density driven flow, g/cm**3
      real*8 p_al  [REFERENCE]  !// specific longitudinal dispersivity of porous media, m
      real*8 p_at [REFERENCE]   !// specific transversal dispersivity of porous media, m
      real*8 p_av  [REFERENCE]  !// specific vertical dispersivity of porous media, m
      real*8 p_hDl [REFERENCE]  !// hydraulic longitudinal dispersivity, m**2/s, diffusities from chemical database
      real*8 p_hDt [REFERENCE]  !// hydraulic transversal dispersivity, m**2/s
      real*8 p_hDv [REFERENCE]  !// hydraulic vertical dispersivity, m**2/s
      real*8 p_nto [REFERENCE]  !// tortuosity factor
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
      real*8  p_dRes1 (7) ! (nICb) ! //should be array Sveta ?
      END
