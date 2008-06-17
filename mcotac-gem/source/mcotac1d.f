
      program mcotac1D
      implicit double precision (a-h,o-z)
      include 'gwheader.inc'
c-coeff :commented out ss(nbasis,nsolid) integert stoichiometric arrays, change to flowt

c      include 'f_gem_node.inc'
#ifdef __MPI
      include 'mpif.h'
#endif

      INTERFACE
        subroutine holdat1d (nxmax,fnamec,hb)
        !DEC$ ATTRIBUTES C :: holdat1d
      include 'gwheader.inc'
        double precision hb(NNODEX+2)
      integer nxmax
      character* (*) fnamec

        !DEC$ ATTRIBUTES REFERENCE :: textc, fnamec
        END SUBROUTINE holdat1d
        END INTERFACE

      INTERFACE
        subroutine setpar (npmax,xmin,xmax,partx,nbox)
        !DEC$ ATTRIBUTES C :: setpar
      include 'gwheader.inc'
        double precision xmin,xmax,partx(NUPMAX)
      integer npmax ,nbox
        !DEC$ ATTRIBUTES REFERENCE :: partx
        END SUBROUTINE setpar
        END INTERFACE

      INTERFACE
        subroutine partid (npmax,nbox,xmin,xmax,partib,dx,partx)
        !DEC$ ATTRIBUTES C :: partid
      include 'gwheader.inc'
        double precision xmin,xmax,partx(NUPMAX), dx(NNODEX+2)
      integer npmax ,nbox, partib(NNODEX+2)
        END SUBROUTINE partid
        END INTERFACE

      INTERFACE
        subroutine concver(npmax,nbox,dx,bn,cn,partib,partx,partic
     *,ismooth,m1,m2)
        !DEC$ ATTRIBUTES C :: concver
      include 'gwheader.inc'
        double precision xmin,xmax,partx(NUPMAX), dx(NNODEX+2)
        double precision bn(NBASIS,NNODEX+2),cn(NCOMPL,NNODEX+2)
      double precision partic(NUPMAX,NCOMPL+NBASIS)
      integer npmax ,nbox, partib(NNODEX+2), ismooth, m1,m2
        END SUBROUTINE concver
        END INTERFACE

      INTERFACE
       subroutine concneu (npmax,nbox,nxmax,xminr,xmaxr,dx,bn,cn,partib,
     *                    partx,partic,bo,co,ismooth,m1,m2)
        !DEC$ ATTRIBUTES C :: concneu
      include 'gwheader.inc'
      double precision xminr,xmaxr,bn(NBASIS,NNODEX+2)
      double precision cn(NCOMPL,NNODEX+2)
      double precision partx(NUPMAX), dx(NNODEX+2),bo(NBASIS,nnodex+2)
      double precision co(NCOMPL,nnodex+2),partic(NUPMAX,NCOMPL+NBASIS)
        integer nxmax,npmax ,nbox, partib(nnodex+2), ismooth, m1,m2
        END SUBROUTINE concneu
        END INTERFACE

      INTERFACE
        subroutine hydro1d(nxmax,h0,hb,tx,am,st,por,ir,qw,qbil,text,vx
     *                    ,dx,icyc,texe,time,fname)

        !DEC$ ATTRIBUTES C :: hydro1d
      include 'gwheader.inc'
        double precision h0(NNODEX+2),hb(NNODEX+2),tx(NNODEX+2)
        double precision st(NNODEX+2),por(NNODEX+2),qw(NNODEX+2)
        double precision qbil(NNODEX+2),am(NNODEX+2)
        double precision dx(NNODEX+2), vx(NNODEX+2),texe,time
        integer nxmax, icyc, ir(nnodex+2)
        character*10 text, fname
        END SUBROUTINE hydro1d
        END INTERFACE

      INTERFACE
        subroutine walk2(npmax,nxmax,ncyc,along,aquer,dm,texe,dx,vx
     *,partx,partxo,xmaxr,xminr,partic,bn,cn,partib,ibpstart,x,bo,co,m1
     *,m2)
        !DEC$ ATTRIBUTES C :: walk2
      include 'gwheader.inc'
        double precision xminr,xmaxr,bn(NBASIS,nnodex+2)
        double precision cn(NCOMPL,nnodex+2)
      double precision dx(NNODEX+2),bo(NBASIS,nnodex+2)
      double precision co(NCOMPL,nnodex+2),partic(NUPMAX,NCOMPL+NBASIS)
        double precision along, aquer, dm(nnodex+2),texe
      double precision vx(NNODEX+2),partx(NUPMAX),partxo(NUPMAX)
        double precision x(nnodex+2)
        integer nxmax,npmax ,ncyc, partib(nnodex+2), ismooth, m1,m2,m3,m4
        integer ibpstart
        END SUBROUTINE walk2
        END INTERFACE

      INTERFACE
        subroutine walk2h(npmax,nxmax,ncyc,along,aquer,dm,texe,dx,vx
     *,partx,partxo,xmaxr,xminr,partic,bn,cn,partib,ibpstart,x,bo,co,m1
     *,m2,por)
        !DEC$ ATTRIBUTES C :: walk2
      include 'gwheader.inc'
        double precision xminr,xmaxr,bn(NBASIS,nnodex+2)
       double precision cn(NCOMPL,nnodex+2)
      double precision dx(NNODEX+2),bo(NBASIS,nnodex+2)
      double precision co(NCOMPL,nnodex+2),partic(NUPMAX,NCOMPL+NBASIS)
        double precision along, aquer, dm(nnodex+2),texe
      double precision vx(NNODEX+2),partx(NUPMAX),partxo(NUPMAX)
        double precision x(nnodex+2), por(NNODEX+2)
        integer nxmax,npmax ,ncyc, partib(nnodex+2), ismooth, m1,m2,m3,m4
        integer ibpstart
        END SUBROUTINE walk2h
        END INTERFACE


c<<<<<<<<<<<<<<<FROM GEMS integration<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c MAIN FORTRAN PROGRAM START IS HERE
c variables for gems-buffer
#ifdef __MPI
       DOUBLE PRECISION, ALLOCATABLE :: bn_subdomain(:) !rank 1
       DOUBLE PRECISION, ALLOCATABLE :: cn_subdomain(:) !rank 1
       DOUBLE PRECISION, ALLOCATABLE :: pn_subdomain(:) !rank 1
c and a second buffer 
       DOUBLE PRECISION, ALLOCATABLE :: bn_domain(:) !rank 1
       DOUBLE PRECISION, ALLOCATABLE :: cn_domain(:) !rank 1
       DOUBLE PRECISION, ALLOCATABLE :: pn_domain(:) !rank 1
#endif

      double precision xxyy, pormin, dmin

c time measurements
	double precision time_gemsmpi, time_gemsmpi_start, time_gemsmpi_end 
	INTEGER argc, iinn, i_gems
c	integer CSTR(100)
c	integer c_to_i(100)
c kg44 itergemstime only needed for debug
c	integer itergemstime(250,51)    ! array for output of iterations done per node during every time step
ckg44 we would like to monitor the iterations of gems (smart initial aprox. vs. simplex init.)
	integer itergems, itergemstotal 
c kg44 why is nodeTypes here fixed? should be allocatable
c      integer nodeTypes(52)
c      integer, allocatable :: nodeTypes(:)
ckg44 variable for interactive definition of gems initial aproximation
	integer gems_PIA
c for F_GEM_CALC_NODE      
	integer iNodeF,idum
      integer p_NodeHandle    !// Node identification handle
      integer p_NodeTypeHY    !// Node type (hydraulic); see typedef NODETYPE
      integer p_NodeTypeMT    !// Node type (mass transport); see typedef NODETYPE
      integer p_NodeStatusFMT !// Node status code FMT; see typedef NODECODEFMT
      integer p_NodeStatusCH  !// Node status code CH;  see typedef NODECODECH
      integer p_IterDone      !// Number of iterations performed by IPM (if not need GEM)

	character  argv1,argv2
	character*10 fname10
	character*10 line
      character*1 chch(100)
c	character*100 CSTR_char30
	character*20 dummystring(25)
	character*100 gems_in_ipmf,gems_dbr_f1,gems_dbr_f2
	character*100 gems_dbr_w
	character*20 dummystringb(20)
c in an ideal world this definitions are only used for MPI stuff
#ifdef __MPI
        integer ierr
	integer npes, i_subdomain_length     !// no of procs, number of grid nodes per processor
        integer recvcount, sendcount
#endif
        integer irank, root      !// rank and identifier for root used also outside MPI
	
c12345678901234567890123456789012345678901234567890123456789012345690
      integer itest,ncyc,nxmax,nymax,isteu,inma,ipfile,ntim
      integer npkt,ir(nnodex+2),npmax,nbox,nboxy,partib(nnodex+2)   
      integer ismooth,iortx(5), i_sorb,j_sorb,iche(nnodex+2),j_decay
     *,ialkali,icyc
c-coeff      integer npin, s,ss
      integer npin, s
      integer p_nDCb, p_nICb, p_nPHb, p_nPSb, p_nPH, gsize3
      integer itimestep_tp, dtprstep, kk1,k1
      real t1,t2
      double precision bog(nbasis,nnodex+2),pog(nnodex+2),          ! this are arrays used for conversion gems-mcotac
     &         cog(ncompl,nnodex+2)
      double precision p_T     !// Temperature T, K                        +      +      -     -
      double precision p_P     !// Pressure P, bar                         +      +      -     -
      double precision p_Vs    !// Volume V of reactive subsystem, cm3     -      -      +     +
      double precision p_Vi    !// Volume of inert subsystem  ?            +      -      -     +
      double precision p_Ms    !// Mass of reactive subsystem, kg          +      +      -     -
      double precision p_Mi    !// Mass of inert part, kg    ?             +      -      -     +
      double precision p_Gs    !// Gibbs energy of reactive subsystem (J)  -      -      +     +
      double precision p_Hs    !// Enthalpy of reactive subsystem (J)      -      -      +     +
      double precision p_Hi    !// Enthalpy of inert subsystem (J) ?       +      -      -     +
      double precision p_IC    !// Effective molal aq ionic strength           -      -      +     +
      double precision p_pH    !// pH of aqueous solution                      -      -      +     +
      double precision p_pe    !// pe of aqueous solution                      -      -      +     +
      double precision p_Eh    !// Eh of aqueous solution, V                   -      -      +     +
c !//  FMT variables (units need dimensionsless form)
      double precision p_Tm    !// actual total simulation time
      double precision p_dt    !// actual time step
      double precision p_dt1   !// priveous time step
      double precision p_vp    !// output time control for postprocessing
      double precision p_Vt    !// total volume of node (voxel) = dx*dy*dz, m**3
      double precision p_eps   !// effective (actual) porosity normalized to 1
      double precision p_Km    !// actual permeability, m**2
      double precision p_Kf    !// actual DARCY`s constant, m**2/s
      double precision p_S	  !// specific storage coefficient, dimensionless
      double precision p_Tr    !// transmissivity m**2/s
      double precision p_h	  !  // actual hydraulic head (hydraulic potential), m
      double precision p_rho   !// actual carrier density for density driven flow, g/cm**3
      double precision p_al    !// specific longitudinal dispersivity of porous media, m
      double precision p_at    !// specific transversal dispersivity of porous media, m
      double precision p_av    !// specific vertical dispersivity of porous media, m
      double precision p_hDl   !// hydraulic longitudinal dispersivity, m**2/s, diffusities from chemical database
      double precision p_hDt   !// hydraulic transversal dispersivity, m**2/s
      double precision p_hDv   !// hydraulic vertical dispersivity, m**2/s
      double precision p_nto   !// node Peclet number, dimensionless
      double precision, allocatable ::  p_xDC(:) ! (nDCb)  !  // DC mole amounts at equilibrium [nDCb]      -      -      +     +
      double precision, allocatable ::  p_gam(:) ! (nDCb)  !  // activity coeffs of DC [nDCb]               -      -      +     +
      double precision, allocatable ::  p_aPH(:) ! (nPHb)  !// Specific surface areas of phases (m2/g)       +      +      -     -
      double precision, allocatable ::  p_xPH(:) ! (nPHb)  !// total mole amounts of phases [nPHb]          -      -      +     +
      double precision, allocatable ::  p_vPS(:) ! (nPSb)  !// phase volume, cm3/mol        [nPSb]          -      -      +     +
      double precision, allocatable ::  p_mPS(:) ! (nPSb)  !// phase (carrier) mass, g      [nPSb]          -      -      +     +
      double precision, allocatable ::  p_bPS(:,:) ! (nICBb,nPSb)  !// bulk compositions of phases  [nPSb][nICb]    -      -      +     +
      double precision, allocatable ::  p_xPA(:) ! (nPSb)  !// amount of carrier in phases  [nPSb] ??       -      -      +     +
      double precision, allocatable ::  p_dul(:) ! (nDCb)  ! // upper kinetic restrictions [nDCb]           +      +      -     -
      double precision, allocatable ::  p_dll(:) ! (nDCb)  ! //  lower kinetic restrictions [nDCb]           +      +      -     -
      double precision, allocatable ::  p_bIC(:) ! (nICb)  !// bulk mole amounts of IC[nICb]                +      +      -     -
      double precision, allocatable ::  p_rMB(:) ! (nICb)  !// MB Residuals from GEM IPM [nICb]             -      -      +     +
      double precision, allocatable ::  p_uIC(:) ! (nICb)  !// IC chemical potentials (mol/mol)[nICb]       -      -      +     +

      double precision xminr,xmaxr,de,xnaohmw,xnaohd,xkohmw,xkohd,xmin
      character *79 title
      character *10 dumb,dumc,dump,input,dtdumb(nbasis),dtdumc(ncompl)
      character *10 tdumb(nbasis), dummy, text, fname
      character *10 datei2,datei3,datei4,datei5,datei6
      character *9  datei1
      character *1 ch,ssw
      character*40 dummy_a40
      character*6 dateiv,dateir,dateiw  !kinet
ckinet<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c       kinetics.inc

        common /rkin1/amin,km,eab,eac,kmp,kmd,eabp,eabd,eacp,eacd,
     .         fg,omega,catinh,ratem,a,b,ap,ad,bp,bd,otherj,otheri,
     .         dratemdc,dqdc,dxdc,rng,volmol,pi,por,texe,gespvfi,
     .         gespvfb,phi

        common /rkin2/vout,rout,wout

        common /ikin1/ikin,ifgp,ifgd,nmineq

        double precision amin(nsolid,nnodex+2),km(nsolid)
        double precision eab(nsolid,nbasis),eac(nsolid,ncompl)
        double precision kmp(nsolid)
        double precision eabp(nsolid,nbasis),eacp(nsolid,ncompl)
        double precision kmd(nsolid)
        double precision eabd(nsolid,nbasis),eacd(nsolid,ncompl)
        double precision fg(nsolid),omega(nsolid),catinh(nsolid)
        double precision ratem(nsolid),a(nsolid),b(nsolid)
        double precision ap(nsolid),ad(nsolid),bp(nsolid),bd(nsolid)
        double precision otherj(nbasis),otheri(ncompl)
        double precision dratemdc(nsolid,nbasis)
        double precision dqdc(nsolid,nbasis),dxdc(ncompl,nbasis)
        double precision rng(nsolid),volmol(nsolid),por(nnodex+2)
        double precision gespvfi(nsolid),gespvfb(nsolid)

        double precision vout(nsolid,nnodex+2),rout(nsolid,nnodex+2)
        double precision wout(nsolid,nnodex+2)

        integer ikin(nsolid),ifg(nsolid)
        integer ifgp(nsolid),ifgd(nsolid)
c   kg44 st needs to be defined
        double precision treal, texe
ckinet>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c
        double precision st(nnodex+2),am(nnodex+2)
     * ,dx(nnodex+2),por1(nnodex+2),poro(nnodex+2)
     * ,partx(nupmax),partxo(nupmax),partic(nupmax,nbasis+ncompl)
     * ,dm(nnodex+2)
     * ,vx(nnodex+2), partiv(nnodex+2)
     * ,hb(nnodex+2), por_null(nnodex+2), tx_null(nnodex+2)
     * ,h0(nnodex+2), qw(nnodex+2), qr(nnodex+2)
     * ,qbil(nnodex+2),	tx(nnodex+2)  
     * ,bnflow(nbasis,nnodex+2),cnflow(ncompl,nnodex+2)
     * ,sumbflowt(nbasis,nnodex+2)
     * ,sumcflowt(ncompl,nnodex+2)
     * ,sumqwat(nnodex+2)
     * ,dqwater(nnodex+2)
     * ,sumbcflowt(ncompl,nnodex+2)
	 

      dimension x(nnodex+2),tmp(nnodex+2),bn(nbasis,nnodex+2)
     1,bo(nbasis,nnodex+2),pn(nsolid,nnodex+2),po(nsolid,nnodex+2)
     1,cn(ncompl,nnodex+2),co(ncompl,nnodex+2)
     2,eqconst(ncompl+nsolid,nnodex+2),acb(nbasis,nnodex+2)
     3,acc(ncompl,nnodex+2),q(nsolid,nnodex+2),tmpold(nnodex+2),
     5z(nbasis+nsolid,nbasis+nsolid),re(nbasis+nsolid)
     7,wconst(nbasis),bc2(nbasis),eh(nnodex+2),pHarr(nnodex+2)
     8,pnw(nsolid),pnd(nsolid)
     9,etc(3)
     9,dpn(nsolid,nnodex+2)
     5,index_i1(nbasis),bi_i1(nbasis),gesb_i1(nbasis)
     6,index_i2(nbasis),bi_i2(nbasis),gesb_i2(nbasis)
     7,gesp_i1(nsolid),gesp_i2(nsolid)
     8,gesc_i1(ncompl),gesc_i2(ncompl)
     5,gespvf_i1(nsolid),gespvf_i2(nsolid)

      dimension cs(nnodex+2),tmpk(nnodex+2)
      common /betat/ betac
      common /dismdl/ idismdl,indxca,indxsi,indxcs,csbndry(0:2),
     1                conhcsm(0:3,2,4)
c >>>
      common /in/ itype,idynam,in1,in2,nxmax,kmax,itemp,lnh,li,lb,
     *dtini,dtmax,doa,dob,vo,tprint(25)
     2,num(6),vjb(nbasis),indexi(nbasis),bi(nbasis)
     3,gesbi(nbasis),indexb(nbasis),bc(nbasis),gesbb(nbasis)
     4,con(ncompl+nsolid,4),s(nbasis,ncompl),vjc(ncompl)
     5,ss(nbasis,nsolid),gespi(nsolid),gespb(nsolid),gesci(ncompl)
     6,gescb(ncompl),err,dtmult,dtdiv
     7,iterm,iterj,iterd,ndiv,lne,lnhc,xmax,gc,rw
     8,temp,tempi,tempo,hcw,hcm,tcm
c                                 added 10/08/93
     9,i_bcp_gemx(nbasis+ncompl+nsolid)
c                                 added nov 2004

      common /inc/ dumb(nbasis),dumc(ncompl),dump(nsolid)
c <<<
      common /maxvals/ maxnmax,maxnb,maxnc,maxnp
      common /tempdep/ itmpdep(ncompl+nsolid)
c >>>
      common /titl/ title

c set irank and root to zero , used also outside MPI in seriall version
       root=0
       irank=0
	gems_PIA=1
	write(*,*)"input for initial gems aproximation (AIA:1, PIA5)"
c	read(*,*)gems_PIA
	gems_PIA=5
	if(.not.((gems_PIA.eq.1).or.(gems_PIA.eq.5))) then
	 gems_PIA=1
	endif
	write(*,*)"gems_PIA: ",gems_PIA
c 
ckg44 init several variables in order to make sure they have the correct values
	pormin=1.e+10
        dmin=1.e+10
	icyc=0
	st=0.0
        write(*,*)st
        ir=0
        por=0
        qr=0
        qw=0
        am=0
        h0=0
        itest=0
        ncyc=0
        nxmax=0
        isteu=0
        inma=0
        ipfile=0
        ntim=0
        npin=0
        npkt=0
        ismooth=0
        i_sorb=0
        j_sorb=0
        j_decay=0
        backg=0
        rd=0
        xlambda=0
        aquer=0
        along=0
        vxx=0
        dm0=0
        dtmax=0
        tmult=0
        dx=0
        de=0
c        c_to_i=0
c        c_to_i1=0
c        c_to_i2=0
        iNode=0
        treal=0 
        texe=0
        st=0
         hb=0
c gems init
        p_nPH=0
         p_NodeHandle=0
         p_NodeTypeHY=0
         p_NodeTypeMT=0
         p_NodeStatusFMT=0
         p_NodeStatusCH=0
         p_IterDone=0
         p_T=0 
         p_P=0
        p_nICb=0 
        p_nDCb=0 
        p_nPHb=0
        p_nPSb=0

        p_Vs=0
        p_Vi=0
        p_Ms=0
        p_Mi=0
         p_Gs=0
        p_Hs=0
        p_Hi=0
        p_IC=0
         p_pH=0
         p_pe=0
         p_Eh=0
          p_Tm=0
        p_dt=0
        p_dt1 =0
         p_ot =0
         p_Vt =0
         p_eps =0
          p_Km =0
          p_Kf =0
          p_S =0
          p_Tr =0
           p_h =0
           p_rho =0
           p_al =0
           p_at=0
          p_av =0
          p_hDl =0
           p_hDt =0
           p_hDv =0
           p_nPe =0
            p_dRes1 =0
            p_dRes2=0
           

c
#ifdef __MPI
C mpi init
      call mpi_init(ierr)
      if (ierr.eq.0) then
         write(*,*)'mpi_init successfull'
      else
         write(*,*)'mpi_init failed'
         stop
      endif
      call mpi_comm_size(MPI_COMM_WORLD,npes, ierr)
      if (ierr.eq.0) then
         write(*,*)'mpi_comm_size successfull'
      else
         write(*,*)'mpi_comm_size failed'
         stop
      endif
      call mpi_comm_rank(MPI_COMM_WORLD,irank,ierr)
      if (ierr.eq.0) then
         write(*,*)'mpi_comm_rank successfull'
      else
         write(*,*)'mpi_comm_rank failed'
         stop
      endif
      print*,npes,irank,' = number of procs and process rank'
c
#endif

c init for gems iterations
	itergems=0
	itergemstotal=0

c<<<<<<  system time initialisation for CPU consumption purposes
      time_initial=secnds(0.)
      if (irank.eq.root) write(*,*)"time initial: ",time_initial
c      pause "initial time"
ccc init itimestep_tp
      itimestep_tp=0
c    init k1
      k1=1
c     dtprstep needs to be initialized
      dtprstep = 0

      i_gems=1	                 ! =0  no GEMS calculation   =1 GEMS calculation
      igems_rw=0               ! =1  read/write gems IO old performed
      i_output=1

      line="a1b2c3d4e5"

c   From GEMS INTEGRATION
c   these files are necessary for GEMS command line execution to start running
      gems_in_ipmf ="MCOTAC-GEM/MySystem-dat.lst" ! list file including nodes' chemical systems
      gems_dbr_f1="MCOTAC-GEM/MySystem-dbr-0-0000.dat"    ! chemical system at the boundary
      gems_dbr_f2="MCOTAC-GEM/MySystem-dbr-0-0001.dat"    

	k1=1                    !time cycle counter

c04<<<<<<<<<<<<<<<<<<<<
      ikin=0
      i_puls=0
      j_chaindecay=0                          !3 member chain
       xlambda1=dlog(2.d0)/2.449e+5/3.15576e+7  ! mother       =m1-2
       xlambda2=dlog(2.d0)/7.702e+4/3.15576e+7  ! daughter     =m1-1
       xlambda3=dlog(2.d0)/1.6e+5/3.15576e+7    ! grant daughte=m1
c       xlambda1=dlog(2.)/10./3.15576e+7  ! 576e7  ! mother       =m1-2
c       xlambda2=dlog(2.)/50./3.15576e+7  ! 576e7 ! daughter     =m1-1
c       xlambda3=dlog(2.)/1./3.15576e+7   ! 576e7  ! grant daughte=m1
      pi=dacos(-1.0d0)
c04>>>>>>>>>>>>>>>>>>>>
c      ialkali=1    
      ialkali=0  
c      ipor=1
      ipor=1
      imodbound=0	 
c nov 2002>>>>>>>>
      ipuls=0
      yzarea=0.0102           !cross sectional area
c      ihydro=1
      ihydro=0
c nov 2002<<<<<<<           
      itest=0
c  --------------
c  set max values
c  --------------
      maxnmax = nnodex
      maxnb = nbasis
      maxnc = ncompl
      maxnp = nsolid 
      npmax = nupmax
c      chemgen=.2

c  --------------------------------------------
c  set constant value used in subroutine eqcont
c    betac = 1/(R*ln10), R is the gas constant
c  N.B. R is in units of kJ/mole K
c  --------------------------------------------
      betac = 1./(8.3143*dlog(10.d0))

c  **************************
c  read input parameters
c  **************************
c04      call datin

c04>>>>>>>
       call datin(itype,idynam,in1,in2,nxmax,kmax,itemp,lnh,li
     +,lb,li_i1,li_i2
     1,dtini,dtmax,doa,dob,vo,tprint
     2,num,vjb,indexi,bi,gesbi,indexb
     3,bc,gesbb,index_i1,bi_i1,gesb_i1,index_i2,bi_i2,gesb_i2
     3,con,s,vjc
     4,ss,gespi,gespb,gesp_i1,gesp_i2,gesc_i1,gesc_i2
     5,gesci,gescb,err,dtmult,dtdiv
     5,gespvf_i1,gespvf_i2
     6,iterm,iterj,iterd,ndiv,lne,lnhc,xmax,gc,rw
     7,temp,tempi,tempo,hcw,hcm,tcm
     8,dumb,dumc,dump
     9,maxnmax,maxnb,maxnc,maxnp
     1,itmpdep
     1,title
c     2,betac
     3,idismdl,indxca,indxsi,indxcs,csbndry,
     1                conhcsm
     9,i_bcp_gemx)
c                                 added nov 2004
      if (irank.eq.root) write (*,*)gespi

      m1=num(1)
      m2=num(2)
      m3=num(3)
      m4=m1+m3
      num(4)=m4
c <<<
      m5=num(5)
c >>>
      m6=m2+m3

      num(6)=m6

c  **************************
c  read iche(nnodex) array used to have heterogeneous chemistry
c  **************************
      fname= "iche01.dat"
	if (irank.eq.root) write(*,*) 'fname = ',fname,nxmax
c      call holdat1d(nxmax,"iche01.dat",hb,text)
#ifdef __GNU
      call holdat1d(%val(nxmax),fname // char(0),hb)
#else
      call holdat1d(nxmax,fname // char(0),hb)
#endif
c allocate memory for nodeTypes
c      if (.not.ALLOCATED(nodeTypes)) then
c	ALLOCATE(nodeTypes(nxmax))
c      endif
c assign values
c      do 1331 ih=1,nxmax
c      iche(ih)=int(hb(ih))
c      nodeTypes(ih) = iche(ih) ! 1
c1331  continue
	if(m3.gt.0)then 
      call solid(m3,pnw,pnd,etc,xnaohmw,xnaohd,xkohmw,xkohd)
      endif
c      write(*,*)iche
c      call wegdat1d(nxmax,"iccccc.dat",iche,"iicc")
c       pause
c  **************************************
c  read input locations for output of breakthrough locations
c  **************************************
      open(20,file ='input_1.dat')
      read(20,*)
      read(20,'(40x,5i5)')(iortx(iort),iort=1,5) 
      close(20)
	if (irank.eq.root) then 
      write(*,'(a22,5(1x,i4))')'breakthrough at nodes '
     *         ,(iortx(iort),iort=1,5)
        endif
      do 3 iort=1,5
      if(iortx(iort).gt.nxmax) then
        write(*,*)'one node out of x range'
        write(*,*)iort,iortx(iort),nxmax
        stop
      endif
  3   continue 


      do 4 m=1,m1
        l=0
        do 7 j=1,10
         if((dumb(m)(j:j)).ne.' ')then
          l=l+1
          dummy(l:l)=(dumb(m)(j:j))
         endif
        if((dumb(m)(j:j)).eq.' '.and.j.gt.2)then
         l=l+1
         dummy(l:l)=' '
        endif
  7     continue
        write(dtdumb(m),'(a4,a6)')'dtb_',dummy
         write(tdumb(m),'(a4,a6)'),'stb_',dummy
  4   continue
      do 5 m=1,m2
        l=0
        do 8 j=1,10
         if((dumc(m)(j:j)).ne.' ')then
          l=l+1
          dummy(l:l)=(dumc(m)(j:j))
         endif
        if((dumc(m)(j:j)).eq.' '.and.j.gt.2)then
         l=l+1
         dummy(l:l)=' '
        endif
  8     continue
        write(dtdumc(m),'(a4,a6)')'dtc_',dummy
  5   continue
	if (irank.eq.root) then 
        write(*,*)(dtdumb(m), m=1,m1)
        write(*,*)(tdumb(m), m=1,m1)
        write(*,*)(dtdumc(m), m=1,m2)
        endif
cpause        pause

c
c  ***************************************
c  assign initial temperature and porosity to each node
c  ***************************************
   43 do 42 n=1,nxmax+2
      por(n)=0.32                  !for initial porcalc 
      poro(n)=0.32       
c      por(n)=phi              ! kinet
      sumqwat(n)=0.
      dqwater(n)=0.
   42 continue
      do 142, n=1,nxmax
  142 tmp(n)=temp
      if (itemp.eq.0) tb=temp

c
c  **************************
c  initialize chemical arrays
c  **************************
   44 do 145 i=1,m1
      bn(i,1)=gesbb(i)
  145 continue
      if (m2.eq.0) go to 47
      do 146 i=1,m2
      cn(i,1)=gescb(i)
  146 continue
   47 if (m3.eq.0) go to 49
      do 48 i=1,m3
      pn(i,1)=gespb(i)
   48 continue
   49 do 54 n=2,nxmax+2
      do 150 i=1,m1
      bn(i,n)=gesbi(i)
      acb(i,n)=1.
  150 continue
      if (m3.eq.0) go to 60
      do 151 i=1,m3
      j=m1+i
      pn(i,n)=gespi(i)
  151 continue
   60 if (m2.eq.0) go to 54
      do 52 i=1,m2
      cn(i,n)=gesci(i)
      acc(i,n)=1.
   52 continue
   54 continue


c   input for dynamic calculations
      call inparf(itest,ncyc,nxmax,isteu,inma,
     *ipfile,ntim,npin,npkt,ismooth,i_sorb,j_sorb,j_decay
     *,backg,rd
     *,xlambda,aquer,along,vxx,dm0,dtmax,tmult,dx,de)
	if (irank.eq.root) then 
      Write(*,*)itest,ncyc,nxmax,isteu,inma,ipfile
     *,ntim,npin,npkt,ismooth,i_sorb,j_sorb,j_decay
     *,backg,rd
     *,xlambda,aquer,along,vxx,dm0,dtmax,tmult,dx,de
       endif
cpause	pause 
c
c
c  *******************************************************************
c  determine the equilibrium concentrations for the initial conditions
c  *******************************************************************
c      write(*,*)'itest isorb=',itest,nxmax,i_sorb 
c04      call initial(itemp,li,lb,in1,in2,lnh,err,nxmax,tb,temp,vo
c04     1,x,itype,num,con,eqconst,tmp,bi,bc,indexi,indexb
c04     2,q,acb,acc,bn,pn,cn,vjb,vjc,s,ss,bc2,lne,eh,idismdl,cs,tmpk
c04     3,iche,i_sorb,ialkali)
c >>>                                  ^       ^  ^ added 9/87

c-coeff      call initial(itemp,li,lb,li_i1,li_i2,in1,in2,lnh,err,nxmax,tb
c-coeff      +,temp,vo
c-coeff      1,x,itype,num,con,eqconst,tmp,bi,bc,indexi,indexb
c-coeff      +,index_i1,bi_i1,gesb_i1,index_i2,bi_i2,gesb_i2
c-coeff      +,gesp_i1,gesp_i2,gesc_i1,gesc_i2
c-coeff      +,gespvf_i1,gespvf_i2
c-coeff      2,q,acb,acc,bn,pn,cn,vjb,vjc,s,ss,lne,eh,idismdl,cs,tmpk
c-coeff      3,iche,i_sorb,ialkali,dumb,dumc,dump,itmpdep)

c  *******************************************************************
c  initial porosity calculation from solids concentration
c  *******************************************************************
       do 1315, nspezx=1,nxmax
       if (iche(nspezx).eq.2)then
            if(m3.gt.0.and.ipor.gt.0)then                                    !  .and.ipor.gt.0
            call porcalc(nspezx,m3,pn,pnw,pnd,etc,por
     *         ,cn,i_sorb,xnaohmw,xnaohd,xkohmw,xkohd)
	  if (por(nspezx).le.1.e-6) por(nspezx)=1.e-6   ! make sure porosity does not get zero
c now change diffusion coefficient
	    dm(nspezx)=dm0*por(nspezx)
            endif
       endif
 1315  continue
	if (irank.eq.root) then 
       write(*,*)'porosities'
       write(*,'(85(f5.3,1x))')(por(nspezx),nspezx=1,nxmax)
	write(*,*)'end porosities'
cpause       pause
      write(*,*)(pn(1,nx),nx=1,nxmax)
        endif
cpause	pause
c      t2=secnds(t1)
c      write (6,2200) t2

c>>>>>02-2003 modified boundary on the right side
c      if(imodbound.gt.0)then  
c       do 1316, nspezx=1,nxmax+1
c       if (nspezx.ge.45)then
c         bn(4,nspezx)=0.
c         bn(5,nspezx)=0.
c       endif
c 1316  continue
c      endif
c>>>>>02-2003 modified boundary on the right side


c03      if (ihydro.eq.1)then                             !if ihydro = 1


c<<<<<<<140895      START HYDROLOGY
c  input ir - array h0 - array  tt - array  s - array
#ifdef __GNU
      call holdat1d(%val(nxmax+2),"ir0001.dat" //char(0),hb)
#else
      call holdat1d(nxmax+2,"ir0001.dat" //char(0),hb)
#endif


      do 3330 ih=1,nxmax+1
3330  ir(ih)=int(hb(ih))
      	if (irank.eq.root) write(*,*)ir
cpause	pause
#ifdef __GNU
      call holdat1d(%val(nxmax+2),"ss0001.dat"//char(0),st)
#else
      call holdat1d(nxmax+2,"ss0001.dat"//char(0),st)
#endif
c	if (irank.eq.root) write(*,*)st
#ifdef __GNU
	call holdat1d(%val(nxmax+2),"por001.dat"//char(0),por)
#else
	call holdat1d(nxmax+2,"por001.dat"//char(0),por)
#endif
	
	if (irank.eq.root) write(*,*)por(1:nxmax)
cpause	pause
	do 3331 ih=1,nxmax+2
c2003      tx_null(ih)= 1.28E-10*(1.-por(ih))**2/por(ih)**3.       !exp 4 specific
      tx_null(ih)= 1.28E-10*(1.-por(ih))**2/por(ih)**3.       !exp 4 specific
      tx(ih)= tx_null(ih)*por(ih)**3/(1.-por(ih))**2
 3331 continue 
#ifdef __GNU
      call holdat1d(%val(nxmax+2),"qr0001.dat"//char(0),qr)
#else
      call holdat1d(nxmax+2,"qr0001.dat"//char(0),qr)
#endif
      	if (irank.eq.root) write(*,*)qr
cpause	pause
#ifdef __GNU
      call holdat1d(%val(nxmax+2),"qn0001.dat"//char(0),qw)
#else
      call holdat1d(nxmax+2,"qn0001.dat"//char(0),qw)
#endif
	if (irank.eq.root) write(*,*)qw
cpause	pause
#ifdef __GNU
      call holdat1d(%val(nxmax+2),"am0001.dat"//char(0),am)
#else
      call holdat1d(nxmax+2,"am0001.dat"//char(0),am)
#endif
	if (irank.eq.root) write(*,*)am
cpause	pause
#ifdef __GNU
      call holdat1d(%val(nxmax+2),"h00001.dat"//char(0),h0)
#else
      call holdat1d(nxmax+2,"h00001.dat"//char(0),h0)
#endif
	if (irank.eq.root) write(*,*)h0
cpause	pause
      do 3332 ih=1,nxmax+2
         hb(ih)=h0(ih)
         qw(ih)=qw(ih)+qr(ih)
c      write(*,*)'main', ih,tx(ih)
 3332  continue

c03      endif                                          !endif ihydro = 1



cpause      pause
c04<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

c<<<<<<<<<<<<<<<FROM GEMS integration<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
cccc -t "cal-dol-boun1-dch.dat", "cal-dol-boun1-dbr-1.dat"  are also in ipmfiles-dat.lst normally
c     open data-ch file for initialisation and names (only once)

	if(i_gems.eq.1) then

      p_NodeHandle=1
c     open data bridge file initially for initialising the spatial distribution of chemical systems
c  first read is for boundary conditons node 1


      if (F_GEM_INIT( gems_in_ipmf ).eq.1) then
		write(*,*) "GEMS init failed"
		stop
       else
	if (irank.eq.root) write(*,*) "GEMS init ok"
       endif

c      do 50 i=1,20
c	  chch(i)=char(CSTR(i))
c  50	continue
c      write(CSTR_char30,'(30a1)')(CSTR(i),i=1,30) 




c   read 2. for p_ variables        
c        if(igems_rw.eq.1)
c      call F_GEM_GET_DCH( p_nICb, p_nDCb, p_nPHb, p_A )
      call F_GEM_GET_DCH( p_nICb, p_nDCb, p_nPHb, p_nPSb, p_nPH)
      	if (irank.eq.root) write(*,*) 'gemsread_dch.dat done'
        write(*,*)p_nICb, p_nDCb, p_nPHb, p_nPSb, p_nPH
        
c        p_nPSb=43
c       gsize3=43
c      double precision  p_xDC(gsize1) ! (nDCb)  !  // DC mole amounts at equilibrium [nDCb]      -      -      +     +
	if (.NOT.allocated(p_xDC)) then
         allocate(p_xDC(p_nDCb))
        endif
c      double precision  p_gam(gsize1) ! (nDCb)  !  // activity coeffs of DC [nDCb]               -      -      +     +
	if (.NOT.allocated(p_gam)) then
         allocate(p_gam(p_nDCb))
        endif
c
c      double precision  p_aPH(gsize3) ! (nPHb)  !// Specific surface areas of phases (m2/g)       +      +      -     -
	if (.NOT.allocated(p_aPH)) then
         allocate(p_aPH(p_nPHb))
        endif
cc      double precision  p_xPH(gsize3) ! (nPHb)  !// total mole amounts of phases [nPHb]          -      -      +     +
	if (.NOT.allocated(p_xPH)) then
         allocate(p_xPH(p_nPHb))
        endif
c      double precision  p_dul(gsize1) ! (nDCb)  ! // upper kinetic restrictions [nDCb]           +      +      -     -
	if (.NOT.allocated(p_dul)) then
         allocate(p_dul(p_nDCb))
        endif
c      double precision  p_dll(gsize1) ! (nDCb)  ! //  lower kinetic restrictions [nDCb]           +      +      -     -
	if (.NOT.allocated(p_dll)) then
         allocate(p_dll(p_nDCb))
        endif
c      double precision  p_bIC(gsize2) ! (nICb)  !// bulk mole amounts of IC[nICb]                +      +      -     -
	if (.NOT.allocated(p_bIC)) then
         allocate(p_bIC(p_nICb))
        endif
c      double precision  p_rMB(gsize2) ! (nICb)  !// MB Residuals from GEM IPM [nICb]             -      -      +     +
	if (.NOT.allocated(p_rMB)) then
         allocate(p_rMB(p_nICb))
        endif
c      double precision  p_uIC(gsize2) ! (nICb)  !// IC chemical potentials (mol/mol)[nICb]       -      -      +     +
	if (.NOT.allocated(p_uIC)) then
         allocate(p_uIC(p_nICb))
        endif
c      double precision  p_vPS(gsize3) ! (nPSb)  !// phase volume, cm3/mol        [nPSb]          -      -      +     +
	if (.NOT.allocated(p_vPS)) then
         allocate(p_vPS(p_nPSb))
        endif
c      double precision  p_mPS(gsize3) ! (nPSb)  !// phase (carrier) mass, g      [nPSb]          -      -      +     +
	if (.NOT.allocated(p_mPS)) then
         allocate(p_mPS(p_nPSb))
        endif
c      double precision  p_bPS(gsize2,gsize3) ! (nICBb,nPSb)  !// bulk compositions of phases  [nPSb][nICb]    -      -      +     +
	if (.NOT.allocated(p_bPS)) then
         allocate(p_bPS(p_nICb,p_nPSb))          ! because of different indexing in fortran and c for arrays!!
        endif
c      double precision  p_xPA(gsize3) ! (nPSb)  !// amount of carrier in phases  [nPSb] ??       -      -      +     +
	if (.NOT.allocated(p_xPA)) then
         allocate(p_xPA(p_nPSb))
        endif

	p_xDC=0.0
      iNode=0
      p_NodeHandle=1
      p_NodeStatusFMT = 1
      p_NodeStatusCH=1

	call F_GEM_READ_NODE( gems_dbr_f1, p_NodeHandle,p_NodeTypeHY
     *,p_NodeTypeMT,p_NodeStatusFMT,p_NodeStatusCH,p_IterDone
     *,p_T, p_P,p_Vs,p_Vi,p_Ms,p_Mi,p_Gs
     *,p_Hs,p_Hi,p_IC,p_pH,p_pe,p_Eh
     *,p_bIC,p_rMB,p_uIC,p_xDC,p_gam
     *,p_dul, p_dll, p_aPH,p_xPH, p_vPS,p_mPS,p_bPS,p_xPA
     *)
	if (irank.eq.root) then 
      write(*,*)'iNode,p_NodeHandle,p_NodeTypeHY,p_NodeTypeMT
     *,p_NodeStatusFMT,p_NodeStatusCH,p_IterDone'
      write(*,*)iNode,p_NodeHandle,p_NodeTypeHY,p_NodeTypeMT
     *,p_NodeStatusFMT,p_NodeStatusCH,p_IterDone,p_T, p_P
	write(*,*)"p_XDC 0", p_xDC
	endif
c      pause
c  here loop to read in additional geochmical systems to define other nodes....
c      aa-initial-dbr-1.dat to aa-initial-dbr-50.dat and aa-boundary-dbr-0.dat are available)
c  <<<<<<<   tranfer of GEMS nomenclature to MCOTAC naming
c     bn=
c     cn=
c     pn=
ckg44
c      if(i_output.eq.1)open(35, file='mco_out.out')
c kg44 
c normalize everything to grid size! dxx x dxx x dxx 
c and for GEMS it is better to work with bigger numbers ;-) so multipy with 10e3
c      gridvol=dxx * dxx * dxx *1000.0
	gridvol=1.0
      write(*,*)"scaling 1:",gridvol,gridvol/p_Vs
      do 1690 n=1,nxmax/2
	 gridvol=dx(n)   ! normalized !
	do 1691 ib=1,m1-1    !charge is last parameter in the list of bn  24.01.2005 but not transported
	bn(ib,n)=p_xDc(i_bcp_gemx(ib))/p_Vs*gridvol    !2)
         bog(ib,n)=bn(ib,n)              ! initial copy of masses        
 1691    continue
	do 1692 ic=1,m2
	   cn(ic,n)=p_xDc(i_bcp_gemx(m1+ic))/p_Vs*gridvol  !10)
           cog(ic,n)=cn(ic,n)
 1692    continue
	do 1693 ip=1,m3
 	   pn(ip,n)=p_xDc(i_bcp_gemx(m1+m2+ip))/p_Vs*gridvol     ! 12)/1.
 1693     continue
c transform the b and c vector to concentrations j_sorb+1 is water!
	do ib = 1, m1-1
 	   bn(ib,n)=bn(ib,n)/bn(j_sorb+1,n)
	enddo
	do ic = 1, m2
 	   cn(ic,n)=cn(ic,n)/bn(j_sorb+1,n)
              ! initial copy of masses        
	enddo
 1690  continue


c     open data bridge file initially for initialising the spatial distribution of chemical systems
c  second read is for initial conditons nodes 2 to nxmax
      
      iNode=0
      p_NodeHandle=1
      p_NodeStatusFMT = 1
      p_NodeStatusCH=1
	call F_GEM_READ_NODE( gems_dbr_f2, p_NodeHandle,p_NodeTypeHY
     *,p_NodeTypeMT,p_NodeStatusFMT,p_NodeStatusCH,p_IterDone
     *,p_T, p_P,p_Vs,p_Vi,p_Ms,p_Mi,p_Gs
     *,p_Hs,p_Hi,p_IC,p_pH,p_pe,p_Eh
     *,p_bIC,p_rMB,p_uIC,p_xDC,p_gam
     *,p_dul, p_dll, p_aPH,p_xPH, p_vPS,p_mPS,p_bPS,p_xPA
     *)

	gridvol=1.0
      write(*,*)"scaling 2:",gridvol,gridvol/p_Vs

      do 1695 n=nxmax/2+1,nxmax+2
	 gridvol=dx(n)   ! normalized !
	do 1696 ib=1,m1-1
	bn(ib,n)=p_xDc(i_bcp_gemx(ib))/p_Vs*gridvol   !   2)   ! i_bcp_gemx(1)=2
         bog(ib,n)=bn(ib,n)         
 1696 continue
	do 1697 ic=1,m2
	cn(ic,n)=p_xDc(i_bcp_gemx(m1+ic))/p_Vs*gridvol  ! i_bcp_gemx(7)=10
 	   cog(ic,n)=cn(ic,n)
 1697 continue
	do 1698 ip=1,m3
	pn(ip,n)=p_xDc(i_bcp_gemx(m1+m2+ip))/p_Vs*gridvol    !12)/1.    ! i_bcp_gemx(12)=12
 1698 continue
c transform the b and c vector to concentrations j_sorb+1 is water!
	do ib = 1, m1-1
 	   bn(ib,n)=bn(ib,n)/bn(j_sorb+1,n)
	enddo
	do ic = 1, m2
 	   cn(ic,n)=cn(ic,n)/bn(j_sorb+1,n)
	enddo
 1695 continue
c here we update porosities from GEMS molar volumes!
c         f_gem_get_molar_volume(int& i, double& Tc, double& P)
	Tc_dummy=25.0
        P_dummy = 1.0
	do n=1,nxmax+2
          poro(n)=por(n)
	  por(n)=0.0
	  do  ip=1,m3
      dum1=f_gem_get_molar_volume(i_bcp_gemx(m1+m2+ip),Tc_dummy,P_dummy)
           dum2= pn(ip,n)
            por(n)=por(n) + dum1*dum2
          enddo
c	 por(n)=1-por(n)/abs((dx(n+1)-dx(n-1))*0.5)   ! normalized !
           gridvol=dx(n)   ! normalized !
	   por(n)=1-por(n)*gridvol*0.1  ! normalized  ...factor 0.1 from definition of molar volume in GEMS!!!
	  if (por(n).le.1.e-6) por(n)=1.e-6   ! make sure porosity does not get zero
c now change diffusion coefficient
        dm(n)=dm0*por(n)
        por_null(n)=por(n)
        poro(n)=por(n)
c2003      tx_null(ih)= 1.28E-10*(1.-por(ih))**2/por(ih)**3.       !exp 4 specific
      tx_null(n)= 1.28E-10*(1.-por(n))**2/por(n)**3.       !exp 4 specific
      tx(n)= tx_null(n)*por(n)**3/(1.-por(n))**2
c change amount of water in the system ....water should be at bn(m1-1)
        enddo
	write(*,*) "porosity update:", por(1:nxmax+2)

	if (irank.eq.root) then 
         write(*,'(13(e8.2,1x))')(bn(ib,2),ib=1,m1),(cn(ic,2),ic=1,m2)
     *,(pn(ip,2),ip=1,m3)
ccc      write(*,'(13(e8.2,1x))')(gemsxDc(ib),ib=1,gemsnDCb)
         write(*,*)' 2 p_xDc(ib) '
         write(*,'(13(e8.2,1x))')(p_xDc(ib),ib=1,p_nDCb)
         endif

#ifdef __MPI
       call MPI_BARRIER (MPI_COMM_WORLD,ierr)
c we define the size of subintervalls for buffer variables
	i_subdomain_length = abs((nxmax)/npes)
	if (irank.eq.root) then 
         write(*,*) 'nodes/procs: ', (nxmax)/npes
         write(*,*) 'nxmax, procs',nxmax, npes
         write(*,*) 'i_subdomain_length', i_subdomain_length
        endif
c  make sure grid size and no of processors fit together
        if (mod((nxmax-2),npes).ne.0) then
         write(*,*) 'parallelization error: check numer of processors'
            write(*,*)'i_subdomain_length is not an integer'
            write(*,*)'mod((nxmax),npes)',mod((nxmax),npes)
            call mpi_finalize(ierr)
            if (ierr.eq.0) then
              write(*,*)'mpi_finalize successfull'
            else
              write(*,*)'mpi_finalize failed'
              stop
            endif
         stop
        endif

c now we allocate memory for the buffers..
c we assume that if bn_subdomain is already allocated the others are also existing
	if (.NOT.allocated(bn_subdomain)) then
         allocate(bn_subdomain((m1-1)*i_subdomain_length))
         allocate(cn_subdomain(m2*i_subdomain_length))
         allocate(pn_subdomain(m3*i_subdomain_length))
c and a second buffer because W. refuses to work with allocate
         allocate(bn_domain((m1-1)*(nxmax)))
         allocate(cn_domain(m2*(nxmax)))
         allocate(pn_domain(m3*(nxmax)))
c	 write(*,*) 'allocated buffers: irank,nxmax,i_subdom',
c     &   irank,nxmax,i_subdomain_length 
c	 write(*,*) m1, m2, m3
         root = 0
	 bn_domain=0.0
	 cn_domain=0.0
         pn_domain=0.0
	 bn_subdomain=0.0
	 cn_subdomain=0.0
         pn_subdomain=0.0
	endif
#endif


	endif                ! end if for i_gems


c>>>>>>>>>>>>>>>>>>>>FROM GEMS integration>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c**************************************************************************


      if (idynam.eq.0) go to 500



c i_sorb for sorption as a complexation
      if(i_sorb.gt.0)then
	if (irank.eq.root) write(*,*)'Sorption for complexes greater 
     *   than No.',i_sorb,'assumed' 
cpause         pause
      endif
      if(j_sorb.gt.0)then
	if (irank.eq.root)write(*,*)'Sorption for basis species greater 
     *   than No.',j_sorb,'assumed' 
cpause         pause
      endif

      mges=m1+m2+m3

c   set grid and subgrid
      x(1)=0.-dx(1)
c init ix and iy
      ix=0
      iy=0
c      uabs= abs(vx(1,1))
      do 1320 ix=2,nxmax+2
       x(ix)=x(ix-1)+dx(ix)
 1320 continue
	if (irank.eq.root) write(*,*)'dx',dx(1:nxmax)  
      do 1322 ix=1,nxmax+2
       vx(ix) = vxx
       dm(ix) = dm0
	if (irank.eq.root) write(*,*)'vx',ix,iy,vx(ix)  
 1322 continue

      xmin=x(1)
      xmax=x(nxmax)
      xminr=x(2)
      xmaxr=x(nxmax-1)
c        write(*,*)'xminr xmaxr',xminr, xmaxr
c	pause
      nbox=nxmax-1
      npmax=npin
      if(ismooth.eq.1)then
        ibpstart=1
      else 
        ibpstart=npmax/nbox
      endif     
	if (irank.eq.root) write(*,*)'ibpstart =',ibpstart, npmax, nbox
      de=de*tmult
c      write(*,*)'npbox xmin xmax',nbox, xmin,xmax,xminr,xmaxr
c       pause
c   set particles in the grid
#ifdef __GNU
	if (irank.eq.root) write(*,*)'before setpar',npmax,xmin,xmax,nbox
       call setpar(%val(npmax),%val(xmin),%val(xmax),partx,
     *             %val(nbox))
#else
      call setpar(npmax,xmin,xmax,partx,nbox)
#endif
c      do 1328 ip=1,npmax
c      write(*,*)'i partx ',ip,partx(ip)
c 1328 continue
c      do 1529 ib =1,nxmax
c      do 1519 id =1,maxnb
c 1519 write(*,'(a6,1x,3i8,1x,e10.4)')'c_b x y j',ib,ic,id,bn(id,ib)
c      do 1518 idc =1,maxnc
c 1518 write(*,'(a6,1x,3i8,1x,e10.4)')'c_c x y j',ib,ic,idc,cn(idc,ib)
c 1529 continue
c>>>>>>>>>>>>>>>>>>>>> nov 2002

c      write(*,*)'main  vor hydro'
      if (ihydro.eq.1)then
#ifdef __GNU
      call hydro1d(%val(nxmax),h0,hb,tx,am,st,
     *   por,ir,qw,qbil,text,vx,dx,
     *   %val(icyc),%val(texe),%val(time),fname)
#else
      call hydro1d(nxmax,h0,hb,tx,am,st,
     *   por,ir,qw,qbil,text,vx,dx,icyc,texe,time,fname)
#endif
	if (irank.eq.root) then 
           write(*,*)'hb',hb
cpause	pause
           write(*,*)'vx',vx
cpause	pause
	endif
	if (irank.eq.root) then
      open(28,file='arrays.dat',form='formatted')
       write(28,'(a4,85(i10,1x))')
     *       'ir  ',(ir(nspezx),nspezx=1,nxmax+2)
       write(28,'(a4,85(i10,1x))')
     *       'iche',(iche(nspezx),nspezx=1,nxmax+2)
       write(28,'(a4,85(e10.3,1x))')
     *       'por ',(por(nspezx),nspezx=1,nxmax+2)
       write(28,'(a4,85(e10.3,1x))')
     *       'tx  ',(tx(nspezx),nspezx=1,nxmax+2)
       write(28,'(a4,85(e10.3,1x))')
     *       'h0  ',(hb(nspezx),nspezx=1,nxmax+2)
       write(28,'(a4,85(e10.3,1x))')
     *       'am  ',(am(nspezx),nspezx=1,nxmax+2)
       write(28,'(a4,85(e10.3,1x))')
     *       'st  ',(st(nspezx),nspezx=1,nxmax+2)
       write(28,'(a4,85(e10.3,1x))')
     *       'qw  ',(qw(nspezx),nspezx=1,nxmax+2)
       write(28,'(a4,85(e10.3,1x))')
     *       'dx  ',(dx(nspezx),nspezx=1,nxmax+2)
       write(28,'(a4,85(e10.3,1x))')
     *       'hb  ',(hb(nspezx),nspezx=1,nxmax+2)
       write(28,'(a4,85(e10.3,1x))')
     *       'vx  ',(vx(nspezx),nspezx=1,nxmax+2)
      close(28)
	endif       ! endif irank.eq.root
      endif                                       ! endif ihydro=1
c>>>>>>>>140895    END HYDROLOGY,  initially 
c       stop

c <<<<< nov 2002



c***  input of concentration arrays at a certain time 
c      write(*,*)'do you want to start calculation at a certain time?'
c      write(*,*)'(Y/N)'
c        pause      
Ckg44 removed this for batch jobs and testing
c       read(*,*)ssw
        ssw='n'
      if(ssw.eq.'y'.or.ssw.eq.'Y')then
        write(*,*)'give: filename by number of tprint(k1)'
        read(*,*)k1
        kk1=k1+96
        time=tprint(k1)
        write(*,*)'kk1',kk1,'time',time
        ch=char(kk1)

c****   read in 2D concentration arrays
c****   solids and porosity 
        write(datei2,'(a5,a1,a4)')'c_2da',ch,'.dat' 
        open(16,file=datei2,form='formatted',status='unknown')
        do 1416 j1=1,m3
        read(16,1476)(pn(j1,nix),nix=1,nxmax)
 1415   continue
        read(16,*)
 1416   continue
 1476   format(21(1x,e10.4))
        read(16,1476)(por(nix),nix=1,nxmax)
 1417   continue
        read(16,'(7x,e20.10)')time
        close(16)
        write(*,*)'data read from',datei2
c****   solutes (basis species and complexes)
        write(datei2,'(a5,a1,a4)')'c_2ds',ch,'.dat' 
        open(16,file=datei2,form='formatted',status='unknown')
        do 1422 j1=1,m1
        read(16,1476)(bn(j1,nix),nix=1,nxmax)
 1423   continue
        read(16,*)dumb(m1)
 1422   continue
        do 1420 j1=1,m2
        read(16,1476)(cn(j1,nix),nix=1,nxmax)
 1421   continue
        read(16,*)dumc(m2)
 1420   continue
        read(16,'(7x,e20.10)')time
        close(16)
        write(*,*)'data read from',datei2
      else
c
c  *********************************************************************
c                 start dynamic mode calculations
c  *********************************************************************
c
        time=0.
        k1=1
        itprint=1
      endif
c
c  **************************************
c  write header of f(t) files
c  **************************************
	if (irank.eq.root) then
      OPEN(25,file='conc1t.dat',form='formatted',status='unknown')
      OPEN(11,file='conc2t.dat',form='formatted',status='unknown')
      OPEN(12,file='conc3t.dat',form='formatted',status='unknown')
      OPEN(13,file='conc4t.dat',form='formatted',status='unknown')
      OPEN(14,file='conc5t.dat',form='formatted',status='unknown')
      write(25,1111)' time     ',(dumb(iw1),iw1=1,m1),
     *  (dumc(iw2),iw2=1,m2),(dump(iw3),iw3=1,m3),'  c/s    ',
     *  ' porosity ','  pH      ','tot._Na_aq','tot._K_aq '
     * ,'  Ca_total','  vx      ','  tx      ' 
     * ,(dtdumb(iw1),iw1=1,m1),(dtdumc(iw2),iw2=1,m2)            ! bnflow , cnflow
     * ,(tdumb(iw1),iw1=1,m1), 'sumqwater ',' delta-t  '         ! sumbnflow, sumwaterflow
     * ,' delta-q  ',(dtdumb(iw1),iw1=1,m1)

 1111   format(2x,80((a10),1x)) 
      write(11,1111)' time    ',(dumb(iw1),iw1=1,m1),
     *  (dumc(iw2),iw2=1,m2),(dump(iw3),iw3=1,m3),'  c/s    ',
     *  ' porosity ','  pH      ','tot._Na_aq','tot._K_aq '
     * ,'  Ca_total','  vx      ' ,'  tx      '  
     * ,(dtdumb(iw1),iw1=1,m1),(dtdumc(iw2),iw2=1,m2)            ! bnflow , cnflow
     * ,(tdumb(iw1),iw1=1,m1), ' sumqwater',' delta-t  '         ! sumbnflow, sumwaterflow
     * ,' delta-q  ',(dtdumb(iw1),iw1=1,m1)
 
      write(12,1111)' time    ',(dumb(iw1),iw1=1,m1),
     *  (dumc(iw2),iw2=1,m2),(dump(iw3),iw3=1,m3),'  c/s    ',
     *  ' porosity ','  pH      ','tot._Na_aq','tot._K_aq '
     * ,'  Ca_total','  vx      ','  tx      '   
     * ,(dtdumb(iw1),iw1=1,m1),(dtdumc(iw2),iw2=1,m2)            ! bnflow , cnflow
     * ,(tdumb(iw1),iw1=1,m1), ' sumqwater',' delta-t  '         ! sumbnflow, sumwaterflow
     * ,' delta-q  ',(dtdumb(iw1),iw1=1,m1)
 
      write(13,1111)' time    ',(dumb(iw1),iw1=1,m1),
     *  (dumc(iw2),iw2=1,m2),(dump(iw3),iw3=1,m3),'  c/s    ',
     *  ' porosity ','  pH      ','tot._Na_aq','tot._K_aq '
     * ,'  Ca_total' ,'  vx      ' ,'  tx      '  
     * ,(dtdumb(iw1),iw1=1,m1),(dtdumc(iw2),iw2=1,m2)            ! bnflow , cnflow
     * ,(tdumb(iw1),iw1=1,m1), ' sumqwater',' delta-t  '         ! sumbnflow, sumwaterflow
     * ,' delta-q  ',(dtdumb(iw1),iw1=1,m1)

      write(14,1111)' time  ',(dumb(iw1),iw1=1,m1),
     *  (dumc(iw2),iw2=1,m2),(dump(iw3),iw3=1,m3),'  c/s    ',
     *  ' porosity ','  pH      ','tot._Na_aq','tot._K_aq '
     * ,'  Ca_total','  vx      ','  tx      '  
     * ,(dtdumb(iw1),iw1=1,m1),(dtdumc(iw2),iw2=1,m2)            ! bnflow , cnflow
     * ,(tdumb(iw1),iw1=1,m1), ' sumqwater' ,' delta-t  '        ! sumbnflow, sumwaterflow
     * ,' delta-q  ',(dtdumb(iw1),iw1=1,m1) 
c
	endif
      texe=dtmax

 235  continue                                              !  next time step interval
      itimestep_tp=itimestep_tp+1
	if (irank.eq.root) then 
      write(*,'(a4,1x,i3,1x,6(e8.2,1x))')'235 ',itimestep_tp,
     *bn(1,1),pn(1,1),pn(2,1),bn(1,2),pn(1,2),pn(2,2)
	endif
      texe=dtmax
c      write(*,*)'DIM F90', nnodex,nbasis,ncompl,nsolid,nupmax,npmax
      told=time
      time=time+texe
c      time=time+dtmax

	if (irank.eq.root) write(*,*)'time',time,told,texe,tprint(k1),k1
c	pause
      deltn=texe
      if (time.le.tprint(k1)) go to 245
      time=tprint(k1)
      texe=time-told
  245 continue

cccc      restzeit = delt
c 2002 c 2002
c      restzeit=tprint(k1)-tprint(k1-1)
c      if (restzeit.eq.0) goto 3333
 335  continue 

      if (ihydro.eq.1)then
c  NEW HYDRAULIC HEAD AND FLOW FIELD 
c  calculate new conductivity depending on  porosity ?? dpor >  xxx 
       do 3432 ih=2,nxmax+1
c dec2002        tx(ih)=tx(ih)*(1.-(por(ih)-poro(ih)))
      tx(ih)=tx_null(ih)*por(ih)**3/(1.-por(ih))**2

c2002      write(*,*)'tx ',ih,tx(ih),por(ih),poro(ih)
 3432  continue
c       pause
       icyc=icyc+1
#ifdef __GNU
       call hydro1d(%val(nxmax),h0,hb,tx,am,st,
     *   por,ir,qw,qbil,text,vx,dx,
     *   %val(icyc),%val(texe),%val(time),fname)
#else
       call hydro1d(nxmax,h0,hb,tx,am,st,
     *   por,ir,qw,qbil,text,vx,dx,icyc,texe,time,fname)
#endif
c       write(*,*)'nach hydro',icyc
      endif                                     ! endif ihydro=1

c  **********************************
c  set old values equal to new values
c  **********************************  
c>>>>>02-2003 modified boundary on the right side
c      if(imodbound.gt.0) then
c       do 1317, nspezx=1,nxmax+1
c       if (nspezx.ge.41)then
c         bn(4,nspezx)=0.
c         bn(5,nspezx)=0.
c       endif
c 1317  continue
c      endif
c>>>>>02-2003 modified boundary on the right side
      do 240 n=1,nxmax-1
        if (j_decay.gt.0)then
          if(m3.gt.0)then
             dissolvef=1
             dp1=pn(1,n)-po(1,n)
             if(dp1.gt.0)bn(j_decay,n)=bn(j_decay,n)+dp1*dissolvef
             dp2=pn(2,n)-po(2,n)
             if(dp2.gt.0)bn(j_decay,n)=bn(j_decay,n)+dp2*dissolvef
             dp3=pn(3,n)-po(3,n)
             if(dp3.gt.0)bn(j_decay,n)=bn(j_decay,n)+dp3*dissolvef
          endif
          bn(j_decay,n)=bn(j_decay,n)*exp((-xlambda*texe))
        endif
 
        do 236 j=1,m1
c               write(*,*)' n j bo bn ',n,j,bo(j,n),bn(j,n)
        if (j_sorb.le.j) bo(j,n)=bn(j,n)         
  236   continue
        if (m2.eq.0) go to 241
        do 237 j=1,m2
        co(j,n)=cn(j,n)
  237   continue
c2003        if(i_sorb.gt.0)cn(i_sorb+1,n)=
c2003     * cn(i_sorb+1,n)*poro(n)/por(n)
c2003        if(i_sorb.gt.0)cn(i_sorb+2,n)=
c2003     *     cn(i_sorb+2,n)*poro(n)/por(n)
  241   if (m3.eq.0) go to 240
        do 238 j=1,m3
cc2004b        pn(j,n)=pn(j,n)*poro(n)/por(n)
c2003        pn(j,n)=pn(j,n)*por(n)/poro(n)
        po(j,n)=pn(j,n)
  238   continue
        tmpold(n)=tmp(n)
        poro(n)=por(n)
c kg44
c	  if (i_output.eq.1.and.n.eq.2)then
c      	    write(35,*) 'node',n,'vor transport' 
c	     write(35,*) 'DCb'
c	    write(35,'(20(e12.6,1x))')(p_xDc(ib),ib=1,p_nDCb)
c        endif
  240 continue
c2002<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c         if(ipor.gt.0)then
c          if (k1.gt.3)then
c           dtnminmin=dtmax
c           dtdxmin=xmax
c           ndtmin=nxmax
c           vmin=vx(1)
c           do 221 n1=2,nxmax-1
c            if(vx(n1).ne.0.)dtnmin1=dx(n1)/vx(n1)/2.
c            if(dm(n1).ne.0.)dtnmin2=dx(n1)**2/2/dm(n1)
c             dtnmin=dtnmin1+dtnmin2
c            if(dtnminmin.le.dtnmin)then
c               dtnminmin=dtnmin
c               dtdxmin=dx(n1)
c               vmin=vx(n1) 
c               ndtmin=n1
c            endif
c  221      continue
c           dtmax=dtnminmin/2.
c           texe=dtmax
c            dtmax=dx(2)/2./vx(2)
c            texe=dx(2)/2./vx(2)
c          endif
c          write(*,*)'k3-ddtmax=',k1, texe,dtdxmin,vmin,ndtmin
c         pause
c         endif
c2002>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c      tmax= abs(dx(1)/( (2.0 *(along*uabs+dm)/dx(1))+uabs))
c      tmax=  dtmax
c      ixxxx=mod(irest+1,2)
c      write(*,*)'ixxxx irest restzeit',ixxxx,irest,restzeit

c ** tasaechlich ausgef. zeitschritt
ccccc      texe = min(restzeit,dtmax)
ccc       texe = min(texe,dtmax)
c       write(*,*)'texe ',texe,restzeit,dtmax,time,told,texe


c*** realtime im walk-loop
      treal=treal+texe

c*** zeit, die noch abzuarbeiten ist 
      restzeit = restzeit -texe
c      write(*,*)'rest ..texe dtmax tp',irest,restzeit,texe,dtmax
c     *,tprint(k1)
 
c      write(*,'(a21,1x,i2,1x,i2,1x,3(e10.4,1x))')
c     *  'k1 it texe time treal',k1,it,texe,time,treal      
c      it=it+1
c       pause 'time'


c   assign concentration at time t to particles
c   1. nuber of particles in each grid cell
c      write (*,*) npmax, nbox,xmin,xmax
#ifdef __GNU
       call partid(%val(npmax),%val(nbox),%val(xmin),%val(xmax),
     *            partib,dx,partx)
#else
      call partid(npmax, nbox,xmin,xmax,partib,dx,partx)
#endif
c      write(*,*)'partid hinter'

c   2. "concentrations of species" assigned to particles 

#ifdef __GNU
      call concver(%val(npmax),%val(nbox),dx,bn,cn,partib,partx,partic,
     * %val(ismooth) ,%val(m1),%val(m2))
#else
      call concver(npmax,nbox,dx,bn,cn,partib,partx,partic,ismooth
     * ,m1,m2)
#endif

c       write(*,*)'conver'

      do 1329 ib=1,npmax
      partxo(ib)=partx(ib)

c      do 1329 ic =1,35
c      write(*,'(a6,1x,2i8,3(1x,e10.4))')'conver',ib,ic,partic(ib,ic),
c     * partx(ib),party(ib)
 1329 continue
c
c  ************************************************************
c  calculate new values of conc. and temp. as functions of time
c  ************************************************************
c
c   move particles during dt
	if (irank.eq.root) then 
      write (*,*) 'walk time:  treal texe', treal, texe
         endif
#ifdef __GNU
      call walk2h(%val(npmax),%val(nxmax),%val(ncyc),%val(along),
     * %val(aquer),dm,%val(texe),dx,vx
     *  ,partx,partxo, %val(xmaxr),%val(xminr),partic,bn,cn,partib,
     *  %val(ibpstart),x,bo,co,%val(m1),%val(m2),por)
#else
      call walk2h(npmax,nxmax,ncyc,along,aquer,dm,texe,dx,vx
     *  ,partx,partxo, xmaxr,xminr,partic,bn,cn,partib,
     *  ibpstart,x,bo,co,m1,m2,por)
#endif
c kg44 walker with variable porosity
c      call walk2h(npmax,nxmax,ncyc,along,aquer,dm,texe,dx,vx
c     *  ,partx,partxo, xmaxr,xminr,partic,bn,cn,partib,
c     *  ibpstart,x,bo,co,m1,m2,por)

c**assign concentrations at t+dt to grid  (including boundary conditions)
c**particle in which nbox

c      write (*,*)'xminr xmaxr',xminr,xmaxr
c      write(*,'(6(e10.4,1x))')(bn(i,1),i=1,m1)
c      write(*,'(6(e10.4,1x))')(bn(i,2),i=1,m1)
c      write(*,'(6(e10.4,1x))')(bn(i,3),i=1,m1)
c      write(*,'(6(e10.4,1x))')(bn(i,11),i=1,m1)
c       pause
c**new concentration in each box
c	if (irank.eq.root) then 
c	do n=1,nxmax
c	 write(*,*)" vectors before concneut, n=",n
c         write(*,'(13(e8.2,1x))')(bn(ib,n),ib=1,m1),(cn(ic,n),ic=1,m2)
c     *,(pn(ip,2),ip=1,m3)
c         enddo
c         endif

#ifdef __GNU
      call concneu(%val(npmax),%val(nbox),%val(nxmax),%val(xminr),
     *  %val(xmaxr),dx,bn,cn,partib,partx,
     *  partic,bo,co,%val(ismooth),%val(m1),%val(m2))
#else
      call concneu(npmax,nbox,nxmax,xminr,
     *  xmaxr,dx,bn,cn,partib,partx,
     *  partic,bo,co,ismooth,m1,m2)
#endif
c	if (irank.eq.root) then 
c	do n=1,nxmax
c	 write(*,*)" vectors after concneut, n=",n
c        write(*,'(13(e8.2,1x))')(bn(ib,n),ib=1,m1),(cn(ic,n),ic=1,m2)
c     *,(pn(ip,2),ip=1,m3)
c         enddo
c         endif

c kg44 changed output format
c	if (irank.eq.root) then 
c      write(*,*)'cneu',itimestep_tp,
c     *bn(1,1),pn(2,1),pn(3,1)
c	endif
c      write(*,'(a4,1x,i3,1x,6(e12.6,1x))')'cneu',itimestep_tp,
c     *bn(1,1),pn(2,1),pn(3,1)
c      write(*,'(a4,1x,i3,1x,6(e12.6,1x))')'cneu',itimestep_tp,
c     *bn(1,2),pn(2,2),pn(3,2)




c>>>>>02-2003 modified boundary on the right side
c      if(imodbound.gt.0)then  
c       do 1318, nspezx=1,nxmax+1
c       if (nspezx.ge.45)then
c         bn(4,nspezx)=0.
c         bn(5,nspezx)=0.
c       endif
c 1318  continue
c      endif 
c>>>>>02-2003 modified boundary on the right side


c<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c   here GEMS calculations for each noce inclusive data transfer
c     write dbr files
c<<<<<< Gems databridge files have to be generated in gems format for each node     
c      bn, cn,pn ->  independent element masses for each node new except for boundary nodes
c  calculation bn, cn pn, to xDC (transported/not transportted)
c  then new xDC trasnfereed to 'bIC' total idepandent masses      

1558  CONTINUE ! NO TRANSPORT DONE BEFORE 

c  <<<<<<<   transfer of GEMS nomenclature to MCOTAC naming
c	do 1699 ib=1,m1
          time_gemsmpi_start=secnds(0.)
	if (i_gems. eq. 1) then     

c ********************************************************
c   neue "total-massen" in den zellen
c ********************************************************
c   if sorption is included for a basis species , basis species with index 
c   greater than j_sorb are not moved 
c   if sorption is included as a complexation reaction the sorbed complexes
c   with index greater than i_sorb are not moved : reset after transport
c   for this species
c****************************
c transform the b and c vector to back to absolute values j_sorb+1 is water!

	do n=1,nxmax+2
c set water!
              bn(j_sorb+1,n)=por(n)*bog(j_sorb+1,n)/por_null(n)
          do ii=1,j_sorb
   	      bn(ii,n)=bn(ii,n)*bn(j_sorb+1,n)
           enddo
            do jj=1,m2
        	     cn(jj,n)=cn(jj,n)*bn(j_sorb+1,n)
             enddo
	enddo


c reset itergems
	itergems=0

#ifdef __MPI
c f irst scatter the inital dataset among the processors
c root obtains the full dataset

        
        if (irank.eq.root) then        
	  do n=1,nxmax
            do ib=1,m1-1
	  	bn_domain(ib+(n-1)*(m1-1))=bn(ib,n)
            enddo
            do ic=1,m2
	  	cn_domain(ic+(n-1)*m2)=cn(ic,n)
            enddo
            do ip=1,m3
	  	pn_domain(ip+(n-1)*m3)=pn(ip,n)  !//we reuse pn_domain for Po_domain
            enddo
          enddo
	endif


c        write(*,*)' MPI Scatter:', MPI_COMM_WORLD
c	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
c         write(*,*) 'barrier ok'
c	if (irank.eq.root) write(*,*)bn_domain
        sendcount = i_subdomain_length*(m1-1)
        recvcount = i_subdomain_length*(m1-1)
        call MPI_SCATTER(bn_domain, sendcount, 
     &      MPI_DOUBLE_PRECISION,
     &	    bn_subdomain, recvcount, MPI_DOUBLE_PRECISION,
     &      root, MPI_COMM_WORLD, ierr);
c	write(*,*)'irank',irank,'values',bn_subdomain
        sendcount = i_subdomain_length*m2
        recvcount = i_subdomain_length*m2
        call MPI_SCATTER(cn_domain, sendcount, 
     &      MPI_DOUBLE_PRECISION,
     &	    cn_subdomain, recvcount, MPI_DOUBLE_PRECISION,
     &      root, MPI_COMM_WORLD, ierr);
        sendcount = i_subdomain_length*m3
        recvcount = i_subdomain_length*m3
        call MPI_SCATTER(pn_domain, sendcount, 
     &      MPI_DOUBLE_PRECISION,
     &	    pn_subdomain, recvcount, MPI_DOUBLE_PRECISION,
     &      root, MPI_COMM_WORLD, ierr);
c

	do 1555 n=1,  i_subdomain_length                 !node loop for GEMS after Transport step
c      goto 1556  ! only node 2 with old gems  values
	do 1596 ib=1,m1-1   ! last value is charge zzz and is mapped to gems 
	p_xDc(i_bcp_gemx(ib))=bn_subdomain(ib+(n-1)*(m1-1))   ! 2)  index_gems(1,...,m1,m+1,m2,   m3) index_mcotac(m1+m2+m3)
 1596 continue
	do 1597 ic=1,m2
	p_xDc(i_bcp_gemx(m1+ic))=cn_subdomain(ic+(n-1)*m2)
 1597 continue
	do 1598 ip=1,m3
	p_xDc(i_bcp_gemx(m1+m2+ip))= pn_subdomain(ip+(n-1)*m3)  !  gemsxDc(12)    ! pn(1,n)   solids not used for transport
c      if(n.eq.2.and.ip.eq.2)p_xDc(i_bcp_gemx(m1+m2+ip))= po(ip,n)/10.
 1598 continue

c	write(*,*)"p_bic",n, p_bIC
c	write(*,*)"node",n,"p_xDc",(p_xDc(ib),ib=1,p_nDCb)


      iNode=  n
      p_NodeHandle=  n
      p_NodeStatusCH= gems_PIA    ! 1 : with simplex PIA; 5 smart PIA
      p_NodeStatusFMT = 1
c<<<<<<  system time initialisation for CPU consumption purposes
c      time_gemsstart=RTC()
      time_gemsstart=secnds(0.)

	idum = F_GEM_CALC_NODE( p_NodeHandle,p_NodeTypeHY,p_NodeTypeMT
     *,p_NodeStatusFMT,p_NodeStatusCH,p_IterDone,p_T, p_P
     *,p_Vs,p_Vi,p_Ms,p_Mi,p_Gs,p_Hs,p_Hi,p_IC,p_pH,p_pe,p_Eh
     *,p_bIC,p_rMB,p_uIC,p_xDC,p_gam, p_dul, p_dll, p_aPH
     *,p_xPH,p_vPS,p_mPS,p_bPS,p_xPA
     *)
c	if (idum.ne.1) then 
c	   write(*,*)"GEMS problem ", idum
c	   stop
c	endif

c  monitor gems iterations
        itergems=itergems+p_IterDone
        itergemstotal=itergemstotal+p_IterDone
	
c	write(*,*)"itimestep_tp,n: ",itimestep_tp,n,p_NodeHandle,
c     &              p_NodeStatusCH,p_IterDone

c      time_gemsend=RTC()
      time_gemsend=secnds(0.)
      time_gemstotal=time_gemstotal+(time_gemsend-time_gemsstart)
c      time_gemstotal=time_gemstotal+ secnds(time_gemsstart)

c kg44
c	if (i_output.eq.1.and.n.eq.2)then
c            write(35,*) 'node',n,'nach GEMS' 
c	    write(35,*) 'DCb', '#######   ',time_gemstotal
c	    write(35,'(20(e12.6,1x))')(p_xDc(ib),ib=1,p_nDCb)
c	    write(35,*) 'ICb'
c	    write(35,'(10(e18.12,1x))')(p_bIC(ib),ib=1,p_nICb)
c	write(35,*) 'b_bPS'
c	    write(35,'(10(e8.2,1x))')(p_bPS(ib),ib=1,p_nDCb)
c	write(35,*) 'xPH'
c	    write(35,'(10(e8.2,1x))')(p_xPH(ib),ib=1,p_nICb)
c	endif
c  <<<<<<<   tranfer of GEMS nomenclature to MCOTAC naming
c     bn=
c     cn=
c     pn=
c      do 1695 n=2,nxmax
	do 1796 ib=1,m1-1
	bn_subdomain(ib+(n-1)*(m1-1))=p_xDc(i_bcp_gemx(ib))
 1796 continue
	do 1797 ic=1,m2
	cn_subdomain(ic+(n-1)*m2)=p_xDc(i_bcp_gemx(m1+ic))
 1797 continue
	do 1798 ip=1,m3
	pn_subdomain(ip+(n-1)*m3)=p_xDc(i_bcp_gemx(m1+m2+ip))
 1798 continue
c kg44
c	if (i_output.eq.1.and.n.eq.2)then
c      	write(35,*) 'node',n,'weit nach GEMS bei 1798' 
c	    write(35,*) 'DCb'
c	    write(35,'(20(e12.6,1x))')(p_xDc(ib),ib=1,p_nDCb)
c      endif

c kg44 only needed for debug
c      itergemstime(itimestep_tp,n)=p_IterDone

c      write(*,*)itimestep_tp,n,itergemstime(itimestep_tp,n),p_IterDone
c	pause
	p_IterDone=0

 1555 continue                 ! end node loop for GEMS after Transport step             


c now do MPI_GATHER
        sendcount = i_subdomain_length*(m1-1)
        recvcount = i_subdomain_length*(m1-1)
      call MPI_AllGather(bn_subdomain, sendcount, 
     &     MPI_DOUBLE_PRECISION,
     &	    bn_domain, recvcount, MPI_DOUBLE_PRECISION,
     &	    MPI_COMM_WORLD,ierr)

        sendcount = i_subdomain_length*m2
        recvcount = i_subdomain_length*m2
      call MPI_AllGather(cn_subdomain, sendcount, 
     &      MPI_DOUBLE_PRECISION,
     &	    cn_domain, recvcount, MPI_DOUBLE_PRECISION,
     &	    MPI_COMM_WORLD,ierr)

        sendcount = i_subdomain_length*m3
        recvcount = i_subdomain_length*m3
      call MPI_AllGather(pn_subdomain, sendcount, 
     &      MPI_DOUBLE_PRECISION,
     &	    pn_domain, recvcount, MPI_DOUBLE_PRECISION,
     &	    MPI_COMM_WORLD,ierr)

	  do n=1,nxmax
            do ib=1,m1-1
	  	bn(ib,n)=bn_domain(ib+(n-1)*(m1-1))
            enddo
            do ic=1,m2
	  	cn(ic,n)=cn_domain(ic+(n-1)*m2)
            enddo
            do ip=1,m3
	  	pn(ip,n)=pn_domain(ip+(n-1)*m3)  
            enddo
          enddo
c       call MPI_BARRIER (MPI_COMM_WORLD,ierr)
c       sendcount = nbasis*nnodex
c       call MPI_BCAST(bn,sendcount,MPI_DOUBLE_PRECISION,
c     &                    root,MPI_COMM_WORLD,ierr)
c       sendcount = ncompl*nnodex      
c       call MPI_BCAST(cn,sendcount,MPI_DOUBLE_PRECISION,
c     &                    root,MPI_COMM_WORLD,ierr)
c       sendcount = nsolid*nnodex      
c       call MPI_BCAST(pn,sendcount,MPI_DOUBLE_PRECISION,
c     &                     root,MPI_COMM_WORLD,ierr)
C

#else       
c	pause "node loop start"
c
	do 1555 n=1,  nxmax                  !node loop for GEMS after Transport step
c      goto 1556  ! only node 2 with old gems  values
	do 1596 ib=1,m1-1

	p_xDc(i_bcp_gemx(ib))=bn(ib,n)   ! 2)  index_gems(1,...,m1,m+1,m2,   m3) index_mcotac(m1+m2+m3)
 1596 continue
	do 1597 ic=1,m2
	p_xDc(i_bcp_gemx(m1+ic))=cn(ic,n)
 1597 continue
	do 1598 ip=1,m3
	p_xDc(i_bcp_gemx(m1+m2+ip))= pn(ip,n)  !  gemsxDc(12)    ! pn(1,n)   solids not used for transport
c      if(n.eq.2.and.ip.eq.2)p_xDc(i_bcp_gemx(m1+m2+ip))= po(ip,n)/10.
 1598 continue

	iNode=  n
      p_NodeHandle=  n
      p_NodeStatusCH= gems_PIA    ! 1 : with simplex PIA; 5 smart PIA
      p_NodeStatusFMT = 1
c<<<<<<  system time initialisation for CPU consumption purposes
c      time_gemsstart=RTC()
      time_gemsstart=secnds(0.)

c	write(*,*)"p_bic",n, p_bIC
c	    write(*,*)"node",n,"p_xDc",(p_xDc(ib),ib=1,p_nDCb)

      idum= F_GEM_CALC_NODE( p_NodeHandle,p_NodeTypeHY,p_NodeTypeMT
     *,p_NodeStatusFMT,p_NodeStatusCH,p_IterDone,p_T, p_P
     *,p_Vs,p_Vi,p_Ms,p_Mi,p_Gs,p_Hs,p_Hi,p_IC,p_pH,p_pe,p_Eh
     *,p_bIC,p_rMB,p_uIC,p_xDC,p_gam, p_dul, p_dll, p_aPH
     *,p_xPH,p_vPS,p_mPS,p_bPS,p_xPA
     *)
c	if (idum.ne.1) then 
c	   write(*,*)"GEMS problem ", idum
c	   stop
c	endif

c  monitor gems iterations
        itergems=itergems+p_IterDone
        itergemstotal=itergemstotal+p_IterDone

c	write(*,*)"itimestep_tp,n: ",itimestep_tp,n,p_NodeHandle,
c     &              p_NodeStatusCH,p_IterDone

c      time_gemsend=RTC()
      time_gemsend=secnds(0.)
      time_gemstotal=time_gemstotal+(time_gemsend-time_gemsstart)
c      time_gemstotal=time_gemstotal+ secnds(time_gemsstart)
ckg44
c	if (i_output.eq.1.and.n.eq.2)then
c      	write(35,*) 'node',n,'nach GEMS' 
c	    write(35,*) 'DCb', '#######   ',time_gemstotal
c	    write(35,'(20(e12.6,1x))')(p_xDc(ib),ib=1,p_nDCb)
c	    write(35,*) 'ICb'
c	    write(35,'(10(e18.12,1x))')(p_bIC(ib),ib=1,p_nICb)
c	write(35,*) 'b_bPS'
c	    write(35,'(10(e8.2,1x))')(p_bPS(ib),ib=1,p_nDCb)
c	write(35,*) 'xPH'
c	    write(35,'(10(e8.2,1x))')(p_xPH(ib),ib=1,p_nICb)
c	endif
c  <<<<<<<   tranfer of GEMS nomenclature to MCOTAC naming
c     bn=
c     cn=
c     pn=
c      do 1695 n=2,nxmax
	do 1796 ib=1,m1-1
	bn(ib,n)=p_xDc(i_bcp_gemx(ib))
 1796 continue
	do 1797 ic=1,m2
	cn(ic,n)=p_xDc(i_bcp_gemx(m1+ic))
 1797 continue
	do 1798 ip=1,m3
	pn(ip,n)=p_xDc(i_bcp_gemx(m1+m2+ip))
 1798 continue


c	if (irank.eq.root) then 
c	 write(*,*)" vectors after transport after chemistry"
c         write(*,'(13(e8.2,1x))')(bn(ib,n),ib=1,m1),(cn(ic,n),ic=1,m2)
c     *,(pn(ip,2),ip=1,m3)
ccc      write(*,'(13(e8.2,1x))')(gemsxDc(ib),ib=1,gemsnDCb)
c         write(*,*)' 2 p_xDc(ib) '
c         write(*,'(13(e8.2,1x))')(p_xDc(ib),ib=1,p_nDCb)
c         endif
c         pause
ckg44
c	if (i_output.eq.1.and.n.eq.2)then
c      	write(35,*) 'node',n,'weit nach GEMS bei 1798' 
c	    write(35,*) 'DCb'
c	    write(35,'(20(e12.6,1x))')(p_xDc(ib),ib=1,p_nDCb)
c      endif

c kg44 only needed for debug
c      itergemstime(itimestep_tp,n)=p_IterDone

c      write(*,*)itimestep_tp,n,itergemstime(itimestep_tp,n),p_IterDone
c	pause
	p_IterDone=0

 1555 continue                 ! end node loop for GEMS after Transport step             
#endif           

c  
ckg44    print out itergems
      write(*,*)'proc, t-step,gems iterations, total iterations'
      write(*,*)irank, itimestep_tp,itergems, 
     &          itergemstotal," CPU time:", time_gemstotal
c here we update porosities from GEMS molar volumes!
c         f_gem_get_molar_volume(int& i, double& Tc, double& P)
	Tc_dummy=25.0
        P_dummy = 1.0
	do n=1,nxmax+2
	  por(n)=0.0
	  do  ip=1,m3
      dum1=f_gem_get_molar_volume(i_bcp_gemx(m1+m2+ip),Tc_dummy,P_dummy)
           dum2= pn(ip,n)
            por(n)=por(n) + dum1*dum2
          enddo
c	 por(n)=1-por(n)/abs((dx(n+1)-dx(n-1))*0.5)   ! normalized !
         gridvol=dx(n)   ! normalized !
         por(n)=1-por(n)*gridvol*0.1  ! normalized  ...factor 0.1 from definition of molar volume in GEMS!!!
         if (por(n).le.1.e-6) por(n)=1.e-6   ! make sure porosity does not get zero	 
c
c now change diffusion coefficient
	 dm(n)=dm0*por(n)
c transform the b and c vector to concentrations j_sorb+1 is water!
	 pormin=min(por(n),pormin)
         dmin=min(dm(n),dmin)


c set water!
              bn(j_sorb+1,n)=por(n)*bog(j_sorb+1,n)/por_null(n)
c
          do ii=1,j_sorb
 	      bn(ii,n)=bn(ii,n)/bn(j_sorb+1,n)
           enddo
           do jj=1,m2
  	          cn(jj,n)=cn(jj,n)/bn(j_sorb+1,n)
           enddo

c         end loop over nodes
        enddo
c	if (irank.eq.root) write(*,*) "porosity update:", por(1:nxmax)
	if (irank.eq.root) write(*,*) "min porosity:",pormin
     &                             ," min diffusion: ",dmin

c

      endif        ! i_gems eq.1  
      time_gemsmpi_end=secnds(0.)
      time_gemsmpi=time_gemsmpi+(time_gemsmpi_end-time_gemsmpi_start)


c<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c          here MCOTAC-chem calculations at each node
      if(i_gems.ne.1)then 
c*** speciation at each node with diffrent total concentrations
      itest=itest+1
c      write(*,*) 'dddd',itest,time,treal,texe
c      pause      
      do 3322 nspez=2,nxmax-1
      if(iche(nspez).eq.2)then

c ********************************************************
c   neue "total-massen" in den zellen
c ********************************************************
c   if sorption is included for a basis species , basis species with index 
c   greater than j_sorb are not moved 
c   if sorption is included as a complexation reaction the sorbed complexes
c   with index greater than i_sorb are not moved : reset after transport
c   for this species
c****************************
         do 677 ii=1,m1
           if(ii.gt.j_sorb.and.j_sorb.gt.0) then
            bbx=bo(ii,nspez)
           else  
            bbx=bn(ii,nspez)
           endif
            do 678 jj=1,m2
                    if(jj.gt.i_sorb.and.i_sorb.gt.0) then
                     bbx=bbx+s(ii,jj)*co(jj,nspez)
c2003              bbx=bbx+s(ii,jj)*co(jj,nspez)/poro(nspez)
c       write(*,*)'i_sorb',dumc(jj),jj,cn(jj,nspez),co(jj,nspez)
                    else
                     bbx=bbx+s(ii,jj)*cn(jj,nspez)
c       write(*,*)'nicht sorb',dumc(jj),jj,cn(jj,nspez),co(jj,nspez)
                    endif
  678       continue
            do 679 kk=1,m3
c2003              bbx=bbx+ss(ii,kk)*po(kk,nspez)  
              bbx=bbx+ss(ii,kk)*po(kk,nspez)                  ! /por(nspez)  2004
  679       continue
              bi(ii)=bbx
c             write(*,*)'nspez bbx m1 ',nspez,bbx,ii
 677     continue
c            if(nspez.eq.3)write (*,*)'init2 v',nspez,pn(1,3)


c*******************************************************  
c speciation
c******************************************************* 
c        if (isio2.gt.0.and. nspez.eq.2)then
c             cccc=pn(1,nspez)
c        endif

cxz      if(i_gems.ne.1)then
	time_initreadstart=secnds(0.)

c-coeff        call init2(itemp,li,lb,in1,in2,lnh,err,nspez,tb,
c-coeff     1     temp,vo,x,itype,num,con,eqconst,tmp,bi,bc,
c-coeff     2          indexi,indexb,q,acb,acc,bn,pn,cn,vjb,vjc,
c-coeff     3          s,ss,bc2,lne,eh,idismdl,cs,tmpk,bo,co,po
c-coeff     4          ,i_sorb,ialkali,dumb,dumc,dump,itmpdep)
c******************************************************* 
c            write(*,*)'init2 n ',nspez
c            if(nspez.eq.3)then 
c	itt=itt+1
c		  write (*,*)'init2 n',treal,time,restzeit,texe,itt 
c            pause
c
c	      endif 
	time_initreadend=secnds(0.)
	time_initreadtotal= time_initreadtotal+
     *                    (time_initreadend-time_initreadstart)
cxz      endif
         do 687 ii=1,m1
                bi(ii)=0.
 687     continue

         if (m3.eq.0) go to 402

c  *******************************************
c  check if a solid is starting to precipitate
c  or has dissolved at any node
c  *******************************************
	if (irank.eq.root) then

        do 390 i=1,m3
        if (pn(i,nspez).le.(zero)) go to 330
cc2004        if (pn(3,nspez).gt.0.) then
cc2004             isio2=isio2+1
cc2004       endif
        if (po(i,nspez).le.(zero)) write(6,2000)nspez,
     1  dump(i),treal,tmp(nspez),eqconst(i,nspez)
 2000   format(/,1x,'at node',i4,1x,a10,
     1  ' started precipitating at time', 1pe12.4,
     1  /,4x,'temp =',0pf8.3,', sol prod =',1pe12.4)
        go to 390
  330   if (po(i,nspez).gt.(zero)) write(6,2050) nspez,
     1  dump(i),treal,tmp(nspez),eqconst(i,nspez)
 2050   format(/,1x,'at node',i4,1x,' solid',a10,' has dissolved at 
     1  time',1pe12.4,/,4x,'temp =',0pf8.3,', sol prod =',1pe12.4)
  390   continue

       endif                ! endif printout on root only


c************************************************
c    calculation of the apparent porosity resulting 
c    from the solid concentrations
c    molvolume[mol/l]/solid_density[kg/l]
c      *solid[mol/fluidvolume in cell nspez]
c      /1000[kg->g]/fluidvolume in cell nspez[l]
c*************************************************
cc2004a      if (iche(nspez).eq.2)then
c<<<<<<<2003
c        do 1860 j=1,m2
c          if(jj.gt.i_sorb.and.i_sorb.gt.0) then
c                cn(jj,nspez)=cn(jj,nspez)*por(nspez)
c          endif
c 1860   continue
cc204        if(m3.gt.0)then
cc204         do 1861 i=1,m3
cc204          pn(i,nspez)= pn(i,nspez)*por(nspez)
cc204 1861    continue

c>>>>2003
cc2004a        call porcalc(nspez,m3,pn,pnw,pnd,etc,por
cc2004a     *   ,cn,i_sorb,xnaohmw,xnaohd,xkohmw,xkohd)
c2004      endif
cc2004a      endif
      
  402   continue
      endif                                     ! if (iche(nspez).eq.2)
c2000>>>>>>>>>>>>>>>> puls

       if(i_puls.eq.1) then
        if(time.gt.1000*3600*24*365.25)then
         do 1424 i1=1,m1
          cn(i1,1)=1.e-20
          bn(i1,1)=1.e-20
c          cn(i1,2)=1.e-20
c          bn(i1,2)=1.e-20
 1424    continue
        else
         do 1425, i1=1, m1
          cn(i1,1)=co(i1,1)
          bn(i1,1)=bo(i1,1)
c          cn(i1,2)=co(i1,2)
c          bn(i1,2)=bo(i1,2)
c          cn(2,2)=co(i1,1)
c          bn(2,1)=1.

 1425    continue
        endif
       endif                                       !i_puls =1

 3322 continue                                     !  nspez=2,nxmax-1
      endif                                     ! i_gems.ne.1 
c 4444 continue
c dec2002<<<<<<<<<
      do 1880 nx=1,nxmax
      do 1881 j1=1,m1
       if(j_sorb.eq.0.or.(j_sorb.ne.0.and.j1.le.j_sorb))then
         bnflow(j1,nx)=1000.*por(nx)*((-1.)*(along*vx(nx)
     *   +dm(nx))*(bn(j1,nx+1)-bn(j1,nx))/dx(nx)+vx(nx)*bn(j1,nx))
         sumbflowt(j1,nx)= sumbflowt(j1,nx)+bnflow(j1,nx)
         sumbcflowt(j1,nx)=sumbflowt(j1,nx)
c         if (nx.eq.13.or.nx.eq.12.or.nx.eq.2)then
c           write(*,*)nx,j1, bnflow(j1,nx) ,sumbflowt(j1,nx),
c     *bn(j1,nx+1),bn(j1,nx),am(nx),por(nx),dm(nx),dx(nx)
c           pause
c         endif
       endif
      do 1882 i1=1,m2
       if(i_sorb.eq.0.or.(i_sorb.ne.0.and.i1.le.i_sorb))then 
         cnflow(i1,nx)=1000.*por(nx)*((-1.*along*vx(nx)+dm(nx))
     *   *(cn(i1,nx+1)-cn(i1,nx))/dx(nx)+vx(nx)*cn(i1,nx))
         sumcflowt(i1,nx)= sumcflowt(i1,nx)+cnflow(i1,nx)
       sumbcflowt(j1,nx)=sumbcflowt(j1,nx)+s(j1,i1)*cnflow(i1,nx)
       endif
 1882 continue
 1881 continue
      dqwater(nx)=am(nx)*por(nx)*vx(nx)*texe*100.               ! 100  : exp_4 specific q in g
      sumqwat(nx)=sumqwat(nx)+dqwater(nx)                       ! vx in m/s - am is in cm^2 cross section
c      if (nx.eq.2)then
c      write(*,'(2(i2,1x),1x,7(e8.2,1x))')k1,
c     *nx,am(nx),por(nx),vx(nx),texe,treal,qwater,sumqwat(nx)
c      pause
c      endif

 1880 continue
c dec20002>>>>>>>>>>>>
c      write(*,*)restzeit,pn
c	pause
c      if( restzeit .gt.0.) goto 335
c 3333 continue
c<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c   here GEMS 
c  <<<<<<<   tranfer of GEMS nomenclature to MCOTAC naming to check solids
c     open data bridge file initially for initialising the spatial distribution of chemical systems
c  loop read for nodes 2 to nxmax
c      gems_dbr_f="aa-initial-dbr-1.dat"


      if(i_gems.eq.1) then
c      do 1799 n=2,nxmax


cc      itergemstime(itimestep_tp,n)=gemsIterDone
cc	gemsIterDone=0

c 1799 continue
c      write(*,*)(bn(ib,2),ib=1,m1)
c      write(*,*)(cn(ic,2),ic=1,m2)
c      write(*,*)(pn(ip,2),ip=1,m3)
c	pause "basis complex solids at n=2 to nxmax"            


c  *******************************************gems
c  check if a solid is starting to precipitate
c  or has dissolved at any node
c  *******************************************gems
       if (irank.eq.root) then
        do 1390 i=1,m3
        if (pn(i,n).le.(zero)) go to 1330
        if (po(i,n).le.(zero)) write(6,2000)n,
     1  dump(i),treal,tmp(n),eqconst(i,n)
        go to 1390
 1330   if (po(i,n).gt.(zero)) write(6,2050) n,
     1  dump(i),treal,tmp(n),eqconst(i,n)
 1390   continue
        endif   ! end printout root only

	endif                                       ! if i_gems=1 

	if (irank.eq.root) then 
      write(*,*)'n_ge',itimestep_tp,
     *bn(1,1),pn(1,1),pn(2,1),bn(1,2),pn(1,2),pn(2,2)
         endif
cgems	 pause "next time step"
c      if(itimestep_tp.gt.10)stop

c  end GEMS reading and re-naming after re-equilibration


c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


c**** output fuer t(k1)
c2002      if (time.eq.tprint(k1)) then
      if (treal.eq.tprint(k1)) then
        kk1=k1+96
        ttt=treal/31557600
	if (irank.eq.root) write(*,*)'kk1',kk1,'time',treal
        ch=char(kk1)

        if (irank.eq.root) then
        write(datei1,'(a4,a1,a4)')'conc',ch,'.dat' 
        open(15,file=datei1,form='formatted',status='unknown')
c       Volume fractions, rates, and omegas       -kinet-
        write(dateiv,'(a1,a1,a4)')'v',ch,'.dat' 
        open(17,file=dateiv,form='formatted',status='unknown')
        write(dateir,'(a1,a1,a4)')'r',ch,'.dat' 
        open(18,file=dateir,form='formatted',status='unknown')
        write(dateiw,'(a1,a1,a4)')'w',ch,'.dat' 
        open(19,file=dateiw,form='formatted',status='unknown')

        write(15,1577)' x       ',(dumb(iw1),iw1=1,m1),
     *  (dumc(iw2),iw2=1,m2),(dump(iw3),iw3=1,m3),' c/s     ',
     *  ' porosity ','  pH      ','   vx     ','  heads   '
 1577   format(2x,50((a10),1x))
c       Volume fractions, rates, and omegas       -kinet-
        write(17,1577)' x       ',(dump(iw3),iw3=1,m3)
        write(18,1577)' x       ',(dump(iw3),iw3=nmineq+1,m3)
        write(19,1577)' x       ',(dump(iw3),iw3=nmineq+1,m3)
        endif    ! end printout root only

        if (irank.eq.root) then
        do 1513 nn1=1,nxmax
        if(nn1.eq.1)then
                xspez=x(nn1)+dx(1)
       else
               xspez=x(nn1)+dx(1)      !   feb 2003    /2
       endif
      if(lnhc.gt.0)then                  !2008   and.acc(lnhc,nn1).gt.0.and.cn(lnhc,nn1).gt0.
          pHarr(nn1)=-dlog10(acc(lnhc,nn1)*cn(lnhc,nn1))
         else
c15112004          pHarr(nn1)=-dlog10(acb(lnh,nn1)*bn(lnh,nn1))
         endif

c      Vol. fractions, rates, and log omega   -kinet-
        if (nn1.eq.1) goto 1513
        write(17,1579)xspez,(vout(i70,nn1),i70=1,m3)   !Vol. fraction
        write(18,1579)xspez,(rout(i71,nn1),i71=nmineq+1,m3)  !rate
        write(19,1579)xspez,(wout(i72,nn1),i72=nmineq+1,m3) !log omega
 1513   write(15,1579)xspez,(bn(i5,nn1),i5=1,m1),
     *                      (cn(i6,nn1),i6=1,m2)
     *                     ,(pn(i7,nn1),i7=1,m3)
     *                     ,cs(nn1),por(nn1)
     *               ,pHarr(nn1), vx(nn1),hb(nn1)
 1579   format(1x,70(e10.3,1x))
        write(15,*)'time = ',treal ,'(s)',ttt , '(J)'
        close(15)
        close(17)   ! kinet
        close(18)   ! kinet
        close(19)   ! kinet
         endif     ! printout root only
        if (irank.eq.root) then
c****  write out particle positions
        write(datei2,'(a5,a1,a4)')'p_ibx',ch,'.dat' 
        open(16,file=datei2,form='formatted',status='unknown')
 1514   write(16,1578)(partib(nix),nix=1,nxmax)
 1578   format(51(1x,i4))
        write(16,*)'time = ',treal
        close(16)
c****  write out 2D concentration arrays
c****  solids and porosity 
        write(datei2,'(a5,a1,a4)')'c_2da',ch,'.dat' 
        open(16,file=datei2,form='formatted',status='unknown')
        do 1516 j1=1,m3
        write(16,1576)(pn(j1,nix),nix=1,nxmax)
 1515   continue
        write(16,*)'  '
 1516   continue
 1576   format(51(1x,e10.4))
        write(16,1576)(por(nix),nix=1,nxmax)
 1517   continue
        write(16,*)'time = ',treal
        close(16)

c****   solutes (basis species and complexes)
        write(datei2,'(a5,a1,a4)')'c_2ds',ch,'.dat' 
        open(16,file=datei2,form='formatted',status='unknown')
        do 1522 j1=1,m1
        write(16,1576)(bn(j1,nix),nix=1,nxmax)
 1523   continue
        write(16,*)dumb(m1),'above'
 1522   continue
        do 1520 j1=1,m2
        write(16,1576)(cn(j1,nix),nix=1,nxmax)
 1521   continue
        write(16,*)dumc(m2),'above'
 1520   continue
        write(16,*)'time = ',treal
        close(16)

       endif       ! end output if root only

c        delt=deltn
        k1=k1+1
        itprint=1
      endif
      if(kmax.gt.k1.and.k1.eq.1)dtprstep=tprint(k1)/10.                     !windows compiler tprint(0) - array limit
      if(kmax.gt.k1.and.k1.gt.1)dtprstep=(tprint(k1)-tprint(k1-1))/10.

c      if((mod(int(time/365.25/24/3600),10)).eq.0
c      if(kmax.gt.k1.and.k1.gt.1)
c     * write(*,*)'k1= ',k1,'dt pr-step= ',dtprstep,texe
c     *,tprint(k1),tprint(k1-1)
c      write(*,*)'k1= ',k1,'dt pr-step= ',dtprstep,texe

c kg44 ----------------this are two very long if ....
       if(k1.gt.1) then
        if (time.gt.(tprint(k1-1)+itprint*dtprstep)
     *.or.time.eq.tprint(k1-1))then
c      if(treal.gt.0.)then
           itprint=itprint+1

           ttt=time/1.
           tty=treal/31557600.
	if (irank.eq.root) then 
         write(*,'(a6,2(e10.2,a6,2x), 2x,a7,i10)')
     *   'time= ',treal,' [sec]',tty,'  [yr]', 'iccyle= ',icyc
         write(*,*)
c     *   'solids iort(1)',(pn(ii,iortx(1)),ii=1,m3),cs(iortx(1)),
c     *    (eqconst(ii,iortx(1)),ii=1,m3)
         endif
c	if (tty.ge.1.e4) pause
c------------------------------------------------------------- 19/09/96
C pH, tot. K and tot Na in aqueous phase 
      if(i_sorb.gt.0.and.ialkali.gt.0)then
      j_K=m1
      j_Na=m1-1
      tot_Na = 0.
      tot_K = 0.
      pH = 0.
c      write(*,*)lnhc,iortx(1),acc(lnhc,iortx(1)),cn(lnhc,iortx(1))
c04      if(lnhc.gt.0)ph = - dlog10(acc(lnhc,iortx(1))*cn(lnhc,iortx(1)))
      if(lnhc.gt.0.or.lnh.gt.0)then                 ! if no pH then no calculation
       if(lnhc.gt.0)then
        ph = - dlog10(acc(lnhc,iortx(1))*cn(lnhc,iortx(1)))
       else
        ph = - dlog10(acb(lnh,iortx(1))*bn(lnh,iortx(1)))
       endif   
      endif

      tot_K = bn(j_k,iortx(1))
      tot_Na = bn(j_Na,iortx(1))
      do 775 i=1,m2
       if(i.ne.i_sorb+2)tot_K=tot_K+s(j_k,i)*cn(i,iortx(1))
c      write(*,*)i,j_k,i_sorb+2,s(j_k,i),cn(i,iortx(1)),tot_K
       if(i.ne.i_sorb+1)tot_Na=tot_Na+s(j_Na,i)*cn(i,iortx(1))
  775 continue
       endif
c---------------------------------------------------------- 
       tot_Ca=bn(1,iortx(1))+s(1,1)*cn(1,iortx(1))
         if (irank.eq.root) write(25,2300)treal
     *               ,(bn(ii,iortx(1)),ii=1,m1)
     *               ,(cn(ii,iortx(1)),ii=1,m2)
     *               ,(pn(ii,iortx(1)),ii=1,m3)
     *               ,cs(iortx(1)),por(iortx(1))
     *               ,pH,tot_Na,tot_K,tot_Ca
     *               ,vx(iortx(1)),tx(iortx(1))
     *               ,(bnflow(ii,iortx(1)),ii=1,m1)
     *               ,(cnflow(ii,iortx(1)),ii=1,m2)
     *               ,(sumbflowt(ii,iortx(1)),ii=1,m1)
     *               ,sumqwat(iortx(1)), texe
     *               ,dqwater(iortx(1))
     *               ,(sumbcflowt(ii,iortx(1)),ii=1,m1)

      if(i_sorb.gt.0.and.ialkali.gt.0)then
      tot_Na = 0.
      tot_K = 0.
      pH = 0.
c04      if(lnhc.gt.0)ph = - dlog10(acc(lnhc,iortx(2))*cn(lnhc,iortx(2)))
      if(lnhc.gt.0.or.lnh.gt.0)then                 ! if no pH then no calculation
       if(lnhc.gt.0)then
        ph = - dlog10(acc(lnhc,iortx(2))*cn(lnhc,iortx(2)))
       else
        ph = - dlog10(acb(lnh,iortx(2))*bn(lnh,iortx(2)))
       endif   
      endif

      tot_K = bn(j_k,iortx(2))
      tot_Na = bn(j_Na,iortx(2))
      do 776 i=1,m2
       if(i.ne.i_sorb+2)tot_K=tot_K+s(j_k,i)*cn(i,iortx(2))
       if(i.ne.i_sorb+1)tot_Na=tot_Na+s(j_Na,i)*cn(i,iortx(2))
  776 continue
      endif

       tot_Ca=bn(1,iortx(2))+s(1,1)*cn(1,iortx(2))
        if (irank.eq.root) write(11,2300)treal
     *               ,(bn(ii,iortx(2)),ii=1,m1)
     *               ,(cn(ii,iortx(2)),ii=1,m2)
     *               ,(pn(ii,iortx(2)),ii=1,m3)
     *               ,cs(iortx(2)),por(iortx(2))
     *               ,pH,tot_Na,tot_K,tot_Ca
     *               ,vx(iortx(2)),tx(iortx(2))
     *               ,(bnflow(ii,iortx(2)),ii=1,m1)
     *               ,(cnflow(ii,iortx(2)),ii=1,m2)
     *               ,(sumbflowt(ii,iortx(2)),ii=1,m1)
     *               ,sumqwat(iortx(2)), texe
     *               ,dqwater(iortx(2))
     *               ,(sumbcflowt(ii,iortx(2)),ii=1,m1)

      if(i_sorb.gt.0.and.ialkali.gt.0)then
      tot_Na = 0.
      tot_K = 0.
      pH = 0.
c0-4      if(lnhc.gt.0)ph = - dlog10(acc(lnhc,iortx(3))*cn(lnhc,iortx(3)))
      if(lnhc.gt.0.or.lnh.gt.0)then                 ! if no pH then no calculation
       if(lnhc.gt.0)then
        ph = - dlog10(acc(lnhc,iortx(3))*cn(lnhc,iortx(3)))
       else
        ph = - dlog10(acb(lnh,iortx(3))*bn(lnh,iortx(3)))
       endif   
      endif

      tot_K = bn(j_k,iortx(3))
      tot_Na = bn(j_Na,iortx(3))
      do 777 i=1,m2
       if(i.ne.i_sorb+2)tot_K=tot_K+s(j_k,i)*cn(i,iortx(3))
       if(i.ne.i_sorb+1)tot_Na=tot_Na+s(j_Na,i)*cn(i,iortx(3))
  777 continue
      endif

       tot_Ca=bn(1,iortx(3))+s(1,1)*cn(1,iortx(3))
        if (irank.eq.root) write(12,2300)treal
     *               ,(bn(ii,iortx(3)),ii=1,m1)
     *               ,(cn(ii,iortx(3)),ii=1,m2)
     *               ,(pn(ii,iortx(3)),ii=1,m3)
     *               ,cs(iortx(3)),por(iortx(3))
     *               ,pH,tot_Na,tot_K,tot_Ca
     *               ,vx(iortx(3)),tx(iortx(3))
     *               ,(bnflow(ii,iortx(3)),ii=1,m1)
     *               ,(cnflow(ii,iortx(3)),ii=1,m2)
     *               ,(sumbflowt(ii,iortx(3)),ii=1,m1)
     *               ,sumqwat(iortx(3)), texe
     *               ,dqwater(iortx(3))
     *               ,(sumbcflowt(ii,iortx(3)),ii=1,m1)

      if(i_sorb.gt.0.and.ialkali.gt.0)then
      tot_Na = 0.
      tot_K = 0.
      pH = 0.
c04      if(lnhc.gt.0)ph = - dlog10(acc(lnhc,iortx(4))*cn(lnhc,iortx(4)))
      if(lnhc.gt.0.or.lnh.gt.0)then                 ! if no pH then no calculation
       if(lnhc.gt.0)then
        ph = - dlog10(acc(lnhc,iortx(4))*cn(lnhc,iortx(4)))
       else
        ph = - dlog10(acb(lnh,iortx(4))*bn(lnh,iortx(4)))
       endif   
      endif

      tot_K = bn(j_k,iortx(4))
      tot_Na = bn(j_Na,iortx(4))
      do 778 i=1,m2
       if(i.ne.i_sorb+2)tot_K=tot_K+s(j_k,i)*cn(i,iortx(4))
       if(i.ne.i_sorb+1)tot_Na=tot_Na+s(j_Na,i)*cn(i,iortx(4))
  778 continue
      endif

       tot_Ca=bn(1,iortx(4))+s(1,1)*cn(1,iortx(4))
        if (irank.eq.root) write(13,2300)treal
     *               ,(bn(ii,iortx(4)),ii=1,m1)
     *               ,(cn(ii,iortx(4)),ii=1,m2)
     *               ,(pn(ii,iortx(4)),ii=1,m3)
     *               ,cs(iortx(4)),por(iortx(4))
     *               ,pH,tot_Na,tot_K,tot_Ca
     *               ,vx(iortx(4)),tx(iortx(4))
     *               ,(bnflow(ii,iortx(4)),ii=1,m1)
     *               ,(cnflow(ii,iortx(4)),ii=1,m2)
     *               ,(sumbflowt(ii,iortx(4)),ii=1,m1)
     *               ,sumqwat(iortx(4)), texe
     *               ,dqwater(iortx(4))
     *               ,(sumbcflowt(ii,iortx(4)),ii=1,m1)

      if(i_sorb.gt.0.and.ialkali.gt.0)then
      tot_Na = 0.
      tot_K = 0.
      pH = 0.
c04      if(lnhc.gt.0)ph = - dlog10(acc(lnhc,iortx(5))*cn(lnhc,iortx(5)))
      if(lnhc.gt.0.or.lnh.gt.0)then                 ! if no pH then no calculation
       if(lnhc.gt.0)then
        ph = - dlog10(acc(lnhc,iortx(5))*cn(lnhc,iortx(5)))
       else
        ph = - dlog10(acb(lnh,iortx(5))*bn(lnh,iortx(5)))
       endif   
      endif

      tot_K = bn(j_k,iortx(5))
      tot_Na = bn(j_Na,iortx(5))
      do 779 i=1,m2
       if(i.ne.i_sorb+2)tot_K=tot_K+s(j_k,i)*cn(i,iortx(5))
       if(i.ne.i_sorb+1)tot_Na=tot_Na+s(j_Na,i)*cn(i,iortx(5))
  779 continue
      endif

       tot_Ca=bn(1,iortx(5))+s(1,1)*cn(1,iortx(5))
        if (irank.eq.root) write(14,2300)treal
     *               ,(bn(ii,iortx(5)),ii=1,m1)
     *               ,(cn(ii,iortx(5)),ii=1,m2)
     *               ,(pn(ii,iortx(5)),ii=1,m3)
     *               ,cs(iortx(5)),por(iortx(5))
     *               ,pH,tot_Na,tot_K,tot_Ca
     *               ,vx(iortx(5)),tx(iortx(5))
     *               ,(bnflow(ii,iortx(5)),ii=1,m1)
     *               ,(cnflow(ii,iortx(5)),ii=1,m2)
     *               ,(sumbflowt(ii,iortx(5)),ii=1,m1)
     *               ,sumqwat(iortx(5)), texe
     *               ,dqwater(iortx(5))
     *               ,(sumbcflowt(ii,iortx(5)),ii=1,m1)
 2300    format(70(1x,e10.3))
c kg44  here the two very long if are finished... 
      endif
      endif  

c  *********************************************************
c  write out the species concentrations
c  *********************************************************
      if (lne.gt.0) call ehcalc(1,nxmax,lne,tmp,bn,eh)

      if (k1.gt.kmax) go to 500
c      write(*,*)(pn(1,nx),nx=1,nxmax)
c	pause "next time step"
      goto 235
  500 continue
          if (irank.eq.root) then
C end of file for breakthrough curves 
         xloc=x(1)+iortx(1)*dx(1)
         write(25,*)'at location ',xloc,' [m]'
         close (25)
         xloc=x(1)+iortx(2)*dx(1)
         write(11,*)'at location ',xloc,' [m]'
         close (11)
         xloc=x(1)+iortx(3)*dx(1)
         write(12,*)'at location ',xloc,' [m]'
         close (12)
         xloc=x(1)+iortx(4)*dx(1)
         write(13,*)'at location ',xloc,' [m]'
         close (13)
         xloc=x(1)+iortx(5)*dx(1)
         write(14,*)'at location ',xloc,' [m]'
         close (14)
       endif     ! end print root only
      time_end=secnds(1.)
      time_interval=time_end-time_initial
      write(*,*) 'Proc No: ',irank,
     &       'time_end  time_interval',time_end, time_interval
      write(*,*) 'Proc No: ',irank,
     &       'time steps calculated timestep_tp = ',itimestep_tp
      write(*,*) 'Proc No: ',irank,
     &        'time used during GEMS calling = ',time_gemstotal
      write(*,*) 'Proc No: ',irank,
     &'time used for Gems-Loop (incl. communication) = ',time_gemsmpi

c      write(*,*) 'time used during GEMS read = ',time_gemsreadtotal
c      write(*,*) 'time used during GEMS write = ',time_gemswritetotal
c      write(*,*) 'time used during MCOTAC-chem-calc=',time_initreadtotal
c kg44 this is not needed..only for debug
c      open (25, file = "iter_array.grd")
c	write(25,'(A4)')"DSAA"
c	write(25,'(2(i2,1x))')80,  51 ! grid dimension
c	write(25,'(2(i2,1x))')1,80    ! xmin xmax
c	write(25,'(2(i2,1x))')1,51   ! ymin ymax
c	write(25,'(2(i3,1x))')0,200   ! zin zmax
c	do 2408 n=51,1,-1
c2408	write(25,'(80(i3,1x))')(itergemstime(it,n),it=1,80)
c      close (25)

c kg44
c      if(i_output.eq.1) close(35)


#ifdef __MPI
C mpi finalize
      call mpi_finalize(ierr)
      if (ierr.eq.0) then
         write(*,*)'mpi_finalize successfull'
      else
         write(*,*)'mpi_finalize failed'
	deallocate (bn_domain) 
        deallocate (cn_domain) 
        deallocate (pn_domain) 
        deallocate (bn_subdomain)
        deallocate (cn_subdomain)
        deallocate (pn_subdomain)
         stop
      endif
c get rid of the buffer for gems
	deallocate (bn_domain) 
        deallocate (cn_domain) 
        deallocate (pn_domain) 
        deallocate (bn_subdomain)
        deallocate (cn_subdomain)
        deallocate (pn_subdomain)
c
#endif

    
cpause	pause "ende"
      stop
c  *********************************************************************

      end



c  *********************************************************************
c     read solid properties

      subroutine solid(m3,pnw,pnd,etc,xnaohmw,xnaohd,xkohmw,xkohd)
      implicit double precision (a-h,o-z)
      include 'gwheader.inc'

      integer m3
      character*10 dumchar
      double precision pnw(nsolid),pnd(nsolid),etc(3)

      open (30,file='solids.dat')

      read (30,*)
      do 10, i=1,m3
      read (30,'(a10,2(e10.5))')dumchar,pnw(i),pnd(i)
c	write(*,'(a10,2(e10.5))')dumchar,pnw(i),pnd(i)
  10  continue
      read (30,'(a10,3(e10.5))')dumchar,etc(1),etc(2),etc(3)
c      write(*,'(a10,3(e10.5))')dumchar,etc(1),etc(2),etc(3)
      read (30,'(a10,2(e10.5))')dumchar,xnaohmw,xnaohd
c      write(*,'(a10,2(e10.5))')dumchar,xnaohmw,xnaohd
      read (30,'(a10,2(e10.5))')dumchar,xkohmw,xkohd
c      write(*,'(a10,2(e10.5))')dumchar,xkohmw,xkohd
c      pause
      close(30)
      return

      end


c  *********************************************************************
c     calculate porosity from solids amount density and molweight

      subroutine porcalc(nspezx,m3,pn,pnw,pnd,etc,por
     *   ,cn,i_sorb,xnaohmw,xnaohd,xkohmw,xkohd)
      implicit double precision (a-h,o-z)

      include 'gwheader.inc'
      integer m3
      dimension pn(nsolid,nnodex+2),pnw(nsolid),pnd(nsolid),etc(3)
      dimension por(nnodex+2),cn(ncompl,nnodex+2)

      xxnaoh=0.
      xxkoh=0.
      porsum=0.
      do 10, i=1,m3
      porsum=porsum+pn(i,nspezx)*pnw(i)/pnd(i)/1000.
c      write(*,'(i3,1x,4(e10.4,1x))')i,pn(i,nspezx),pnw(i),pnd(i),porsum
   10 continue
      por(nspezx)=1.-(etc(3)*etc(1)/etc(2)/1000.+ porsum
     *   ) *por(nspezx)

      if (i_sorb.gt.0)then
         xxnaoh=cn(i_sorb+1,nspezx)*xnaohmw/xnaohd/1000.
         xxkoh=cn(i_sorb+2,nspezx)*xkohmw/xkohd/1000.
           por(nspezx)=por(nspezx)-xxnaoh  -xxkoh
      endif
c       write(*,*)nspezx, por(nspezx)
c      pause
      return

      end
      


c*==========================================================================*
c
c                            Unterprogramm: holdat1df
c               Lesen eines Feldes aus Datei Fname.dat
c
c
c    nx, ny   : Anzahl der Knoten in X- bzw. in Y-Richtung
c    fname    : File-Name
c    hb[i][j] : einzulesendes Feld
c    text     : Bezeichnung des Feldes
c
c   return =  Fehlernummer ierr
c
c    07.12.01
c===========================================================================*/

      subroutine holdat1df(nxmax,cname,hb,text)
       include 'gwheader.inc'

        double precision hb(NNODEx+2)
        integer nxmax,nymax,ihb(NNODEx+2)
        character*10 text,  cname

        open(31, file=cname)
        read (31, *)nxmaxx,faktor
        write(*,*)'hol nx  faktor', nxmax,faktor, cname
         read(31, *)(ihb(i),i=1,nxmax)
c         read(31, 1010)(ihb(i),i=1,nxmax)
   10   continue
 1010   format (520i4)
        do 20, i=1,nxmax
        hb(i)= faktor*ihb(i)
 20     continue
        close (31)
        return
        end
