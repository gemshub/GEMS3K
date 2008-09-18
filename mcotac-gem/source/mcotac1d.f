
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

      double precision, allocatable ::  ph_domain(:) 
      double precision, allocatable ::  ph_subdomain(:) 
      double precision, allocatable ::  eh_domain(:) 
      double precision, allocatable ::  eh_subdomain(:) 

      double precision, allocatable ::  gems_iterations_domain(:) 
      double precision, allocatable ::  gems_iterations_subdomain(:) 


#endif

      double precision xxyy, pormin, dmin
c debug for failing GEMS iterations
	double precision maxmola,maxmolb
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
	double precision, allocatable ::  gems_iterations(:) ! number of Iterations for each node
        integer, allocatable ::  node_distribution(:) ! list with ordered nodes (ascending numbers of GEMS iterations)
        integer i_reorder
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
      double precision bog(nbasis,nnodex+2),pog(nsolid,nnodex+2),          ! this are arrays used for conversion gems-mcotac
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
        double precision texe
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
c 
ckg44 init several variables in order to make sure they have the correct values
        irestart=0
	t_outint=1e+10
	idum=0
	pormin=1.e+10
        dmin=1.e+10
	icyc=0
	st=0.0
	pHarr=0.0
	eh=0.0
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
      igems_rw=0               ! =1  read/write gems IO old performedvtk
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

	if(m3.gt.0)then 
      call solid(m3,pnw,pnd,etc,xnaohmw,xnaohd,xkohmw,xkohd)
      endif

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
   43 do 42 n=1,nxmax
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
   49 do 54 n=2,nxmax
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

c pressure and temperature values for gems
	Tc_dummy=25.0
        P_dummy = 1.0
c   input for dynamic calculations
      call inparf(itest,ncyc,nxmax,isteu,inma,
     *ipfile,ntim,npin,npkt,ismooth,i_sorb,j_sorb,j_decay
     *,backg,rd
     *,xlambda,aquer,along,vxx,dm0,dtmax,tmult,dx,de,gems_PIA
     *,Tc_dummy,P_dummy,irestart,it_out)

	if (irank.eq.root) then 
      Write(*,*)itest,ncyc,nxmax,isteu,inma,ipfile
     *,ntim,npin,npkt,ismooth,i_sorb,j_sorb,j_decay
     *,backg,rd
     *,xlambda,aquer,along,vxx,dm0,dtmax,tmult,dx,de,gems_PIA
     *,Tc_dummy,P_dummy, irestart, it_out

       endif

c	gems_PIA=1
c	write(*,*)"input for initial gems aproximation (AIA:1, PIA5)"
c	read(*,*)gems_PIA
c	gems_PIA=5
	if(.not.((gems_PIA.eq.1).or.(gems_PIA.eq.5))) then
	 gems_PIA=1
	endif
	write(*,*)"gems_PIA: ",gems_PIA





c03      if (ihydro.eq.1)then                             !if ihydro = 1


c<<<<<<<140895      START HYDROLOGY
c  input ir - array h0 - array  tt - array  s - array
#ifdef __GNU
      call holdat1d(%val(nxmax),"ir0001.dat" //char(0),hb)
#else
      call holdat1d(nxmax,"ir0001.dat" //char(0),hb)
#endif


      do 3330 ih=1,nxmax
3330  ir(ih)=int(hb(ih))
      	if (irank.eq.root) write(*,*)ir
cpause	pause
#ifdef __GNU
      call holdat1d(%val(nxmax),"ss0001.dat"//char(0),st)
#else
      call holdat1d(nxmax,"ss0001.dat"//char(0),st)
#endif
c	if (irank.eq.root) write(*,*)st
#ifdef __GNU
	call holdat1d(%val(nxmax),"por001.dat"//char(0),por)
#else
	call holdat1d(nxmax,"por001.dat"//char(0),por)
#endif
	
	if (irank.eq.root) write(*,*)por(1:nxmax)
cpause	pause
	do 3331 ih=1,nxmax
c2003      tx_null(ih)= 1.28E-10*(1.-por(ih))**2/por(ih)**3.       !exp 4 specific
      tx_null(ih)= 1.28E-10*(1.-por(ih))**2/por(ih)**3.       !exp 4 specific
      tx(ih)= tx_null(ih)*por(ih)**3/(1.-por(ih))**2
 3331 continue 
#ifdef __GNU
      call holdat1d(%val(nxmax),"qr0001.dat"//char(0),qr)
#else
      call holdat1d(nxmax,"qr0001.dat"//char(0),qr)
#endif
      	if (irank.eq.root) write(*,*)qr
cpause	pause
#ifdef __GNU
      call holdat1d(%val(nxmax),"qn0001.dat"//char(0),qw)
#else
      call holdat1d(nxmax,"qn0001.dat"//char(0),qw)
#endif
	if (irank.eq.root) write(*,*)qw
cpause	pause
#ifdef __GNU
      call holdat1d(%val(nxmax),"am0001.dat"//char(0),am)
#else
      call holdat1d(nxmax,"am0001.dat"//char(0),am)
#endif
	if (irank.eq.root) write(*,*)am
cpause	pause
#ifdef __GNU
      call holdat1d(%val(nxmax),"h00001.dat"//char(0),h0)
#else
      call holdat1d(nxmax,"h00001.dat"//char(0),h0)
#endif
	if (irank.eq.root) write(*,*)h0
cpause	pause
      do 3332 ih=1,nxmax
         hb(ih)=h0(ih)
         qw(ih)=qw(ih)+qr(ih)
c      write(*,*)'main', ih,tx(ih)
 3332  continue

c03      endif                                          !endif ihydro = 1



cpause      pause
c04<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

c<<<<<<<<<<<<<<<FROM GEMS integration<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



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
c make sure we have an equilibrated file
	iNode=  n
      p_NodeHandle=  n
      p_NodeStatusCH= gems_PIA    ! 1 : with simplex PIA; 5 smart PIA
      p_NodeStatusFMT = 1
c<<<<<<  system time initialisation for CPU consumption purposes
c      time_gemsstart=RTC()
      time_gemsstart=secnds(0.)

c	write(*,*)"p_bic",n, p_bIC
c	    write(*,*)"node",n,"p_xDc",(p_xDc(ib),ib=1,p_nDCb)

      call F_GEM_CALC_NODE( p_NodeHandle,p_NodeTypeHY,p_NodeTypeMT
     *,p_NodeStatusFMT,p_NodeStatusCH,p_IterDone,p_T, p_P
     *,p_Vs,p_Vi,p_Ms,p_Mi,p_Gs,p_Hs,p_Hi,p_IC,p_pH,p_pe,p_Eh
     *,p_bIC,p_rMB,p_uIC,p_xDC,p_gam, p_dul, p_dll, p_aPH
     *,p_xPH,p_vPS,p_mPS,p_bPS,p_xPA,idum,idebug
     *)

	if (idum.ne.1) then 
	   write(*,*)"GEMS problem for first dbr file", idum
	  write(*,*) "P_IterDone: ",p_IterDone
	  write(*,*)" Please look into ipmlog.txt"
c           pause
c	   stop
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
	ref_vol=p_Vs
     
      do 1690 n=1,nxmax/2
	 gridvol=dx(1)/dx(n)*ref_vol/p_Vs   ! normalized !
         write(*,*)"scaling 1:",gridvol,p_Vs,dx(1),dx(2)
	do 1691 ib=1,m1-1    !charge is last parameter in the list of bn  24.01.2005 but not transported
	bn(ib,n)=p_xDc(i_bcp_gemx(ib))*gridvol    !2)
         bog(ib,n)=bn(ib,n)              ! initial copy of masses        
 1691    continue
	do 1692 ic=1,m2
	   cn(ic,n)=p_xDc(i_bcp_gemx(m1+ic))*gridvol  !10)
           cog(ic,n)=cn(ic,n)
 1692    continue
	do 1693 ip=1,m3
 	   pn(ip,n)=p_xDc(i_bcp_gemx(m1+m2+ip))*gridvol     ! 12)/1.
 	   pog(ip,n)=pn(ip,n)
 1693     continue



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

c make sure we have an equilibrated file
	iNode=  n
      p_NodeHandle=  n
      p_NodeStatusCH= gems_PIA    ! 1 : with simplex PIA; 5 smart PIA
      p_NodeStatusFMT = 1
c<<<<<<  system time initialisation for CPU consumption purposes
c      time_gemsstart=RTC()
      time_gemsstart=secnds(0.)

c	write(*,*)"p_bic",n, p_bIC
c	    write(*,*)"node",n,"p_xDc",(p_xDc(ib),ib=1,p_nDCb)

      call F_GEM_CALC_NODE( p_NodeHandle,p_NodeTypeHY,p_NodeTypeMT
     *,p_NodeStatusFMT,p_NodeStatusCH,p_IterDone,p_T, p_P
     *,p_Vs,p_Vi,p_Ms,p_Mi,p_Gs,p_Hs,p_Hi,p_IC,p_pH,p_pe,p_Eh
     *,p_bIC,p_rMB,p_uIC,p_xDC,p_gam, p_dul, p_dll, p_aPH
     *,p_xPH,p_vPS,p_mPS,p_bPS,p_xPA,idum,idebug
     *)
	if (idum.ne.1) then 
	   write(*,*)"GEMS problem for second dbr file", idum
	  write(*,*) "P_IterDone: ",p_IterDone
	  write(*,*)" Please look into ipmlog.txt"
c           pause
c	   stop
	endif

	write(*,*)"p_XDC 1", p_xDC
	gridvol=1.0

      do 1695 n=nxmax/2+1,nxmax
	 gridvol=dx(1)/dx(n)* ref_vol/p_Vs   ! normalized !
       write(*,*)"scaling 2:",gridvol,p_Vs

	do 1696 ib=1,m1-1
	bn(ib,n)=p_xDc(i_bcp_gemx(ib))*gridvol   !   2)   ! i_bcp_gemx(1)=2
         bog(ib,n)=bn(ib,n)         
 1696 continue
	do 1697 ic=1,m2
	cn(ic,n)=p_xDc(i_bcp_gemx(m1+ic))*gridvol  ! i_bcp_gemx(7)=10
 	   cog(ic,n)=cn(ic,n)
 1697 continue
	do 1698 ip=1,m3
	pn(ip,n)=p_xDc(i_bcp_gemx(m1+m2+ip))*gridvol    !12)/1.    ! i_bcp_gemx(12)=12
 	   pog(ip,n)=pn(ip,n)
 1698 continue


 1695 continue
c here we update porosities from GEMS molar volumes!
c         f_gem_get_molar_volume(int& i, double& Tc, double& P)

	do n=1,nxmax
	  por(n)=0.0
          poro(n)=por(n)
	  do  ip=1,m3
      dum1=f_gem_get_molar_volume(i_bcp_gemx(m1+m2+ip),Tc_dummy,P_dummy)
           dum2= pn(ip,n)
            por(n)=por(n) + dum1*dum2
          enddo
c	 por(n)=1-por(n)/abs((dx(n+1)-dx(n-1))*0.5)   ! normalized !
           gridvol=dx(n)/dx(1)*1.e-4/ref_vol *0.1  ! normalized ... convert from cm^3 to m^3   factor 0.1 from definition of molar volume in GEMS!!!

	   por(n)=1-gridvol*por(n)  ! normalized  ...f
	  if (por(n).le.1.e-6) por(n)=1.e-6   ! make sure porosity does not get zero
c now change diffusion coefficient
        dm(n)=dm0
        por_null(n)=por(n)
        poro(n)=por(n)
c2003      tx_null(ih)= 1.28E-10*(1.-por(ih))**2/por(ih)**3.       !exp 4 specific
      tx_null(n)= 1.28E-10*(1.-por(n))**2/por(n)**3.       !exp 4 specific
      tx(n)= tx_null(n)*por(n)**3/(1.-por(n))**2
c convert aqueous species to concentrations
c
          do ii=1,j_sorb
 	      bn(ii,n)=bn(ii,n)/por(n)
           enddo
           do jj=1,m2
  	       cn(jj,n)=cn(jj,n)/por(n)
           enddo
 
       enddo
	if (irank.eq.root) then 
	write(*,*) "porosity update:", por(1:nxmax)
         write(*,'(13(e8.2,1x))')(bn(ib,2),ib=1,m1),(cn(ic,2),ic=1,m2)
     *,(pn(ip,2),ip=1,m3)
ccc      write(*,'(13(e8.2,1x))')(gemsxDc(ib),ib=1,gemsnDCb)
         write(*,*)' 2 p_xDc(ib) '
         write(*,'(13(e8.2,1x))')(p_xDc(ib),ib=1,p_nDCb)
         endif


#ifdef __GNU
	idum=vtkout(%val(itimestep_tp),%val(time),%val(nxmax),%val(m1),
     &     %val(m2),%val(m3),dx,bn,cn,pn,por,eh,pHarr,dumb, dumc, dump)
#else
	idum=vtkout(itimestep_tp,time,nxmax,m1,
     &     m2,m3,dx,bn,cn,pn,por,eh,pHarr,dumb,dumc,dump)
#endif
c	pause



#ifdef __MPI
       call MPI_BARRIER (MPI_COMM_WORLD,ierr)
c we define the size of subintervalls for buffer variables
	i_subdomain_length = abs((nxmax)/npes)
	if (irank.eq.root) then 
         write(*,*) 'nodes/procs: ', (nxmax)/npes
         write(*,*) 'nxmax, procs',nxmax, npes
         write(*,*) 'i_subdomain_length', i_subdomain_length
        endif
c allocate memory for arrays (gems_iterations,node_distribution) 
	if (.NOT.allocated(gems_iterations)) then
         allocate(gems_iterations(nxmax))
        endif
	if (.NOT.allocated(node_distribution)) then
         allocate(node_distribution(nxmax))
        endif
	do i=1,nxmax
	  gems_iterations=1.0
	  node_distribution(i)=i
	enddo
	
        write(*,*) node_distribution
c	pause
	i_reorder=100   ! set counter for load balancing --counting is done backward ;-)
c  make sure grid size and no of processors fit together
        if (mod((nxmax),npes).ne.0) then
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
c   ph_domain
	allocate(ph_domain(nxmax))
        allocate(ph_subdomain(i_subdomain_length))
c   pe_domain
	allocate(eh_domain(nxmax))
        allocate(eh_subdomain(i_subdomain_length))
c gems iterations is already allocated..now we take the subdomain
        allocate(gems_iterations_domain(nxmax))
        allocate(gems_iterations_subdomain(i_subdomain_length))

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
         ph_domain=0.0
         ph_subdomain=0.0
         pe_domain=0.0
         pe_subdomain=0.0
	endif
#endif

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
      do 1320 ix=2,nxmax
       x(ix)=x(ix-1)+dx(ix)
 1320 continue
	if (irank.eq.root) write(*,*)'dx',dx(1:nxmax)  
      do 1322 ix=1,nxmax
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
     *       'ir  ',(ir(nspezx),nspezx=1,nxmax)
       write(28,'(a4,85(i10,1x))')
     *       'iche',(iche(nspezx),nspezx=1,nxmax)
       write(28,'(a4,85(e10.3,1x))')
     *       'por ',(por(nspezx),nspezx=1,nxmax)
       write(28,'(a4,85(e10.3,1x))')
     *       'tx  ',(tx(nspezx),nspezx=1,nxmax)
       write(28,'(a4,85(e10.3,1x))')
     *       'h0  ',(hb(nspezx),nspezx=1,nxmax)
       write(28,'(a4,85(e10.3,1x))')
     *       'am  ',(am(nspezx),nspezx=1,nxmax)
       write(28,'(a4,85(e10.3,1x))')
     *       'st  ',(st(nspezx),nspezx=1,nxmax)
       write(28,'(a4,85(e10.3,1x))')
     *       'qw  ',(qw(nspezx),nspezx=1,nxmax)
       write(28,'(a4,85(e10.3,1x))')
     *       'dx  ',(dx(nspezx),nspezx=1,nxmax)
       write(28,'(a4,85(e10.3,1x))')
     *       'hb  ',(hb(nspezx),nspezx=1,nxmax)
       write(28,'(a4,85(e10.3,1x))')
     *       'vx  ',(vx(nspezx),nspezx=1,nxmax)
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
c        ssw='n'
      if(irestart.eq.1)then               ! 1 read backup file - 0 start from t=0
c      if (irank.eq.root)then
       call readdump(bn,cn,pn,por,dm,
     * time,itimestep_tp,tprint,k1,m1,m2,m3,nxmax)
        kk1=k1+96
        write(*,*)'kk1',k1,'time',time, "dtmax",dtmax
        ch=char(kk1)
        if(k1.ge.kmax)stop 'end of defined calculation time reached'
c	endif
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
     *  ' porosity ','  pH      '

 1111   format(2x,80((a10),1x)) 
      write(11,1111)' time    ',(dumb(iw1),iw1=1,m1),
     *  (dumc(iw2),iw2=1,m2),(dump(iw3),iw3=1,m3),'  c/s    ',
     *  ' porosity ','  pH      '
 
      write(12,1111)' time    ',(dumb(iw1),iw1=1,m1),
     *  (dumc(iw2),iw2=1,m2),(dump(iw3),iw3=1,m3),'  c/s    ',
     *  ' porosity ','  pH      '
 
      write(13,1111)' time    ',(dumb(iw1),iw1=1,m1),
     *  (dumc(iw2),iw2=1,m2),(dump(iw3),iw3=1,m3),'  c/s    ',
     *  ' porosity ','  pH      '

      write(14,1111)' time  ',(dumb(iw1),iw1=1,m1),
     *  (dumc(iw2),iw2=1,m2),(dump(iw3),iw3=1,m3),'  c/s    ',
     *  ' porosity ','  pH      '
c
	endif
      texe=dtmax

 235  continue                                              !  next time step interval
      itimestep_tp=itimestep_tp+1
c	if (irank.eq.root) then 
c      write(*,'(a4,1x,i3,1x,6(e8.2,1x))')'235 ',itimestep_tp,
c     *bn(1,1),pn(1,1),pn(2,1),bn(1,2),pn(1,2),pn(2,2)
c	endif
      texe=dtmax
c      write(*,*)'DIM F90', nnodex,nbasis,ncompl,nsolid,nupmax,npmax
      told=time
      time=time+texe
c      time=time+dtmax

	if (irank.eq.root) 
     &   write(*,*)'time',itimestep_tp,time,told,texe,tprint(k1),k1
c	pause
      deltn=texe
      if (time.le.tprint(k1)) go to 245
      time=tprint(k1)
      texe=time-told
  245 continue

 335  continue 

      if (ihydro.eq.1)then
c  NEW HYDRAULIC HEAD AND FLOW FIELD 
c  calculate new conductivity depending on  porosity ?? dpor >  xxx 
       do 3432 ih=1,nxmax
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
      do 240 n=1,nxmax
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



c      write(*,*)'rest ..texe dtmax tp',irest,texe,dtmax
c     *,tprint(k1)
 


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
      write (*,*) 'walk time:  time texe', time, texe
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


#ifdef __GNU
      call concneu(%val(npmax),%val(nbox),%val(nxmax),%val(xminr),
     *  %val(xmaxr),dx,bn,cn,partib,partx,
     *  partic,bo,co,%val(ismooth),%val(m1),%val(m2))
#else
      call concneu(npmax,nbox,nxmax,xminr,
     *  xmaxr,dx,bn,cn,partib,partx,
     *  partic,bo,co,ismooth,m1,m2)
#endif





c>>>>>0modified boundary on both sides
c      if(imodbound.gt.0)then  
            do ib=1,j_sorb
	  	bn(ib,1)=bog(ib,1)/por(1)
	  	bn(ib,2)=bog(ib,2)/por(2)
	  	bn(ib,nxmax)=bog(ib,nxmax)/por(nxmax)
	  	bn(ib,nxmax-1)=bog(ib,nxmax-1)/por(nxmax-1)
            enddo
            do ic=1,m2
	  	cn(ic,1)=cog(ic,1)/por(1)
	  	cn(ic,2)=cog(ic,2)/por(2)
	  	cn(ic,nxmax)=cog(ic,nxmax)/por(nxmax)
	  	cn(ic,nxmax-1)=cog(ic,nxmax-1)/por(nxmax-1)
            enddo
            do ip=1,m3
	  	pn(ip,1)=pog(ip,1)  
	  	pn(ip,2)=pog(ip,2)  
	  	pn(ip,nxmax)=pog(ip,nxmax)  
	  	pn(ip,nxmax-1)=pog(ip,nxmax-1)  
            enddo


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

	do n=1,nxmax
c set water!
              bn(j_sorb+1,n)=por(n)*bog(j_sorb+1,n)/por_null(n)
          do ii=1,j_sorb
   	        bn(ii,n)=bn(ii,n)*por(n)
           enddo
            do jj=1,m2
        	 cn(jj,n)=cn(jj,n)*por(n)
             enddo
	enddo


c reset itergems
	itergems=0

#ifdef __MPI
c f irst scatter the inital dataset among the processors
c root obtains the full dataset
c for better load balancing, it is possible to manipulate
c the node ordering with node_distribution
        
        if (irank.eq.root) then        
	  do n=1,nxmax,npes
          do i=1,npes
            do ib=1,m1-1
      bn_domain(ib+(i-1+n-1)*(m1-1))=bn(ib,node_distribution(i-1+n))
            enddo
            do ic=1,m2
	  cn_domain(ic+(i-1+n-1)*m2)=cn(ic,node_distribution(n+i-1))
            enddo
            do ip=1,m3
	  pn_domain(ip+(i-1+n-1)*m3)=pn(ip,node_distribution(n+i-1))  
            enddo
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
	maxmola=0.0 ! for debugging of GEMS 
	do 1555 n=1,  i_subdomain_length                 !node loop for GEMS after Transport step
c      goto 1556  ! only node 2 with old gems  values
	do 1596 ib=1,m1-1   ! last value is charge zzz and is mapped to gems 
	p_xDc(i_bcp_gemx(ib))=bn_subdomain(ib+(n-1)*(m1-1))   ! 2)  index_gems(1,...,m1,m+1,m2,   m3) index_mcotac(m1+m2+m3)
        maxmola=max(maxmola,p_xDc(i_bcp_gemx(ib)))
 1596 continue
	do 1597 ic=1,m2
	p_xDc(i_bcp_gemx(m1+ic))=cn_subdomain(ic+(n-1)*m2)
        maxmola=max(maxmola,p_xDc(i_bcp_gemx(m1+ic)))
 1597 continue
	do 1598 ip=1,m3
	p_xDc(i_bcp_gemx(m1+m2+ip))= pn_subdomain(ip+(n-1)*m3)  !  gemsxDc(12)    ! pn(1,n)   solids not used for transport
        maxmola=max(maxmola,p_xDc(i_bcp_gemx(m1+m2+ip)))
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

c	call F_GEM_CALC_NODE( p_NodeHandle,p_NodeTypeHY,p_NodeTypeMT
c     *,p_NodeStatusFMT,p_NodeStatusCH,p_IterDone,p_T, p_P
c     *,p_Vs,p_Vi,p_Ms,p_Mi,p_Gs,p_Hs,p_Hi,p_IC,p_pH,p_pe,p_Eh
c     *,p_bIC,p_rMB,p_uIC,p_xDC,p_gam, p_dul, p_dll, p_aPH
c     *,p_xPH,p_vPS,p_mPS,p_bPS,p_xPA,idum,idebug
c     *)
	idum=1
	if (idum.ne.1) then 
	   write(*,*)"GEMS problem ", idum
	  write(*,*) "P_IterDone: ",p_IterDone
	  write(*,*)" Please look into ipmlog.txt"
c            pause
c	   stop
	endif

c	if (idum.ne.1) then 
c	   write(*,*)"GEMS problem ", idum
c	   stop
c	endif

	pH_subdomain(n)=p_pH
	eh_subdomain(n)=p_Eh

c  monitor gems iterations
        itergems=itergems+p_IterDone
        itergemstotal=itergemstotal+p_IterDone

	gems_iterations_subdomain(n)=
     &   gems_iterations_subdomain(n)+p_IterDone    ! conversion from integer to double

c	write(*,*)"itimestep_tp,n: ",itimestep_tp,n,p_NodeHandle,
c     &              p_NodeStatusCH,p_IterDone

c      time_gemsend=RTC()
      time_gemsend=secnds(0.)
      time_gemstotal=time_gemstotal+(time_gemsend-time_gemsstart)
	maxmolb=0.0
	do 1796 ib=1,m1-1
	bn_subdomain(ib+(n-1)*(m1-1))=p_xDc(i_bcp_gemx(ib))
	maxmolb=max(bn_subdomain(ib+(n-1)*(m1-1)),maxmolb)
 1796 continue
	do 1797 ic=1,m2
	cn_subdomain(ic+(n-1)*m2)=p_xDc(i_bcp_gemx(m1+ic))
	maxmolb=max(cn_subdomain(ic+(n-1)*m2),maxmolb)
 1797 continue
	do 1798 ip=1,m3
	pn_subdomain(ip+(n-1)*m3)=p_xDc(i_bcp_gemx(m1+m2+ip))
	maxmolb=max(pn_subdomain(ip+(n-1)*m3),maxmolb)
 1798 continue

	if((maxmola.gt.100.0).or.(maxmola.gt.100.0).or.
     &     (maxmolb/maxmola.gt.2.0).or.(idum.ne.1)) then
	   write(*,*) "Problem with maxmol or idum!", maxmola,maxmolb
	   STOP
	endif

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

        sendcount = i_subdomain_length
        recvcount = i_subdomain_length
      call MPI_AllGather(ph_subdomain, sendcount, 
     &     MPI_DOUBLE_PRECISION,
     &	    ph_domain, recvcount, MPI_DOUBLE_PRECISION,
     &	    MPI_COMM_WORLD,ierr)
      call MPI_AllGather(eh_subdomain, sendcount, 
     &     MPI_DOUBLE_PRECISION,
     &	    eh_domain, recvcount, MPI_DOUBLE_PRECISION,
     &	    MPI_COMM_WORLD,ierr)


	  do n=1,nxmax,npes
	  do i=1,npes
	    pHarr(node_distribution(n+i-1))=ph_domain(n+i-1)            
	    eh(node_distribution(n+i-1))=eh_domain(n+i-1)            
            do ib=1,m1-1
	  bn(ib,node_distribution(n+i-1))=bn_domain(ib+(n-1+i-1)*(m1-1))
            enddo
            do ic=1,m2
	  cn(ic,node_distribution(n+i-1))=cn_domain(ic+(n-1+i-1)*m2)
            enddo
            do ip=1,m3
	  pn(ip,node_distribution(n+i-1))=pn_domain(ip+(n-1+i-1)*m3)  
            enddo
          enddo
          enddo
c	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
c         write(*,*) 'barrier ok'

	    if (irank.eq.root) write(*,*) "i_reorder =",i_reorder

	if (i_reorder.le.0) then    ! we reorder only every 10th step
        sendcount = i_subdomain_length
        recvcount = i_subdomain_length
      call MPI_AllGather(gems_iterations_subdomain, sendcount, 
     &     MPI_DOUBLE_PRECISION,
     &	    gems_iterations_domain, recvcount, MPI_DOUBLE_PRECISION,
     &	    MPI_COMM_WORLD,ierr)

	if (irank.eq.root) then
	  do n=1,nxmax,npes
	  do i=1,npes
                gems_iterations(node_distribution(n+i-1))=
     &          gems_iterations_domain(n+1-1)
          enddo
          enddo

c now calculate new node_distribution
c	pause
c	write(*,*)"gems_iterations ", gems_iterations
c	write(*,*)"node_distribution before", node_distribution
	write(*,*)"node_distribution before", node_distribution
	write(*,*)"gems_iterations ", gems_iterations
	call sortrx(nxmax,1.0/gems_iterations, node_distribution)   ! we would like to get the result in ascending order, therefore use 1.0/iteration
	write(*,*)"node_distribution after"
	write(*,*)"gems_iterations ", 
     &        (gems_iterations(node_distribution(i)),i=1,nxmax)
	endif !endif for  root

            i_reorder=100   ! reset i_reorder
           gems_iterations_subdomain=0
c distribute mpi data
	call MPI_BCAST(node_distribution,nxmax, MPI_INTEGER,root,
     &                 MPI_COMM_world,ierr) 

	else
	    i_reorder=i_reorder-1
        endif ! endif for i_reorder process  
#else       
c	pause "node loop start"
c
	do 1555 n=1,  nxmax                  !node loop for GEMS after Transport step
c      goto 1556  ! only node 2 with old gems  values
	do 1596 ib=1,m1-1
	p_xDc(i_bcp_gemx(ib))=bn(ib,n)   ! 2)  index_gems(1,...,m1,m+1,m2,   m3) index_mcotac(m1+m2+m3)
        maxmola=max(maxmola,p_xDc(i_bcp_gemx(ib)))
 1596 continue
	do 1597 ic=1,m2
	p_xDc(i_bcp_gemx(m1+ic))=cn(ic,n)
        maxmola=max(maxmola,p_xDc(i_bcp_gemx(ic+m1)))
 1597 continue
	do 1598 ip=1,m3
	p_xDc(i_bcp_gemx(m1+m2+ip))= pn(ip,n)  !  gemsxDc(12)    ! pn(1,n)   solids not used for transport
        maxmola=max(maxmola,p_xDc(i_bcp_gemx(ip+m1+m2)))
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

c      call F_GEM_CALC_NODE( p_NodeHandle,p_NodeTypeHY,p_NodeTypeMT
c     *,p_NodeStatusFMT,p_NodeStatusCH,p_IterDone,p_T, p_P
c     *,p_Vs,p_Vi,p_Ms,p_Mi,p_Gs,p_Hs,p_Hi,p_IC,p_pH,p_pe,p_Eh
c     *,p_bIC,p_rMB,p_uIC,p_xDC,p_gam, p_dul, p_dll, p_aPH
c     *,p_xPH,p_vPS,p_mPS,p_bPS,p_xPA,idum,idebug
c     *)
	idum=1
	if (idum.ne.1) then 
	   write(*,*)"GEMS problem ", idum
	  write(*,*) "P_IterDone: ",p_IterDone
	  write(*,*)" Please look into ipmlog.txt"
c           pause
c	   stop
	endif
	pHarr(n)=p_pH
        eh(n)=p_Eh
c  monitor gems iterations
        itergems=itergems+p_IterDone
        itergemstotal=itergemstotal+p_IterDone


      time_gemsend=secnds(0.)
      time_gemstotal=time_gemstotal+(time_gemsend-time_gemsstart)
c
	maxmolb=0.0
	do 1796 ib=1,m1-1
	bn(ib,n)=p_xDc(i_bcp_gemx(ib))
	maxmolb=max(bn(ib,n),maxmolb)
 1796 continue
	do 1797 ic=1,m2
	cn(ic,n)=p_xDc(i_bcp_gemx(m1+ic))
	maxmolb=max(cn(ic,n),maxmolb)
 1797 continue
	do 1798 ip=1,m3
	pn(ip,n)=p_xDc(i_bcp_gemx(m1+m2+ip))
	maxmolb=max(pn(ip,n),maxmolb)
 1798 continue

	if((maxmola.gt.100.0).or.(maxmola.gt.100.0).or.
     &     (maxmolb/maxmola.gt.2.0).or.(idum.ne.1)) then
	   write(*,*) "Problem with maxmol or idum!", maxmola,maxmolb
	   STOP
	endif


	p_IterDone=0

 1555 continue                 ! end node loop for GEMS after Transport step             
#endif           

c  
ckg44    print out itergems
      write(*,*)'proc, t-step,gems iterations, total iterations'
      write(*,*)irank, itimestep_tp,itergems, 
     &          itergemstotal," CPU time:", time_gemstotal


	do n=1,nxmax
	  por(n)=0.0
	  do  ip=1,m3
      dum1=f_gem_get_molar_volume(i_bcp_gemx(m1+m2+ip),Tc_dummy,P_dummy)
           dum2= pn(ip,n)
            por(n)=por(n) + dum1*dum2
          enddo
c	 por(n)=1-por(n)/abs((dx(n+1)-dx(n-1))*0.5)   ! normalized !
           gridvol=dx(n)/dx(1)*1.e-4/ref_vol *0.1  ! normalized ... convert from cm^3 to m^3   factor 0.1 from definition of molar volume in GEMS!!!
	   por(n)=1-gridvol*por(n)  ! normalized  ...f
         if (por(n).le.1.e-6) por(n)=1.e-6   ! make sure porosity does not get zero	 
c
c now change diffusion coefficient
c	 dm(n)=dm0*por(n)
c transform the b and c vector to concentrations j_sorb+1 is water!
	 pormin=min(por(n),pormin)
         dmin=min(dm(n),dmin)


c
          do ii=1,j_sorb
 	       bn(ii,n)=bn(ii,n)/por(n)
           enddo
           do jj=1,m2
  	        cn(jj,n)=cn(jj,n)/por(n)
           enddo

c         end loop over nodes
        enddo
c	if (irank.eq.root) write(*,*) "porosity update:", por(1:nxmax)
	if (irank.eq.root) write(*,*) "min porosity:",pormin
     &                             ," min diffusion: ",dmin

c

      time_gemsmpi_end=secnds(0.)
      time_gemsmpi=time_gemsmpi+(time_gemsmpi_end-time_gemsmpi_start)

c reset boundary value

c<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c          here MCOTAC-chem calculations at each node
c kg44 deleted

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
c     *nx,am(nx),por(nx),vx(nx),texe,time,qwater,sumqwat(nx)
c      pause
c      endif

 1880 continue


c**** output fuer t(k1)
c2002      if (time.eq.tprint(k1)) then
      if (time.eq.tprint(k1)) then
        kk1=k1+96
        ttt=time/31557600
c	if (irank.eq.root) write(*,*)'kk1',kk1,'time',time
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
        write(15,*)'time = ',time ,'(s)',ttt , '(J)'
        close(15)
        close(17)   ! kinet
        close(18)   ! kinet
        close(19)   ! kinet
         endif     ! printout root only
        if (irank.eq.root) then
c****  write out particle positions
        write(datei2,'(a5,a1,a4)')'p_ibx',ch,'.dat' 


       endif       ! end output if root only

c        delt=deltn
        k1=k1+1
        itprint=1
      endif
      if(kmax.gt.k1.and.k1.eq.1)dtprstep=tprint(k1)/10.                     !windows compiler tprint(0) - array limit
      if(kmax.gt.k1.and.k1.gt.1)dtprstep=(tprint(k1)-tprint(k1-1))/10.



c kg44 ----------------this are two very long if ....
c I do not have any idea how to change the settings for output in the files
c I just use every timestep 
c       if(k1.gt.1) then
c        if (time.gt.(tprint(k1-1)+itprint*dtprstep)
c     *.or.time.eq.tprint(k1-1))then
c      if(time.gt.0.)then
           itprint=itprint+1

           ttt=time/1.
           tty=time/31557600.
	if (irank.eq.root) then 
         write(*,*)
     *   'time= ',time,' [sec]',tty,'  [yr]', 'iccyle= ',icyc
         write(*,*)
c     *   'solids iort(1)',(pn(ii,iortx(1)),ii=1,m3),cs(iortx(1)),
c     *    (eqconst(ii,iortx(1)),ii=1,m3)
         endif
c	if (tty.ge.1.e4) pause
c------------------------------------------------------------- 19/09/96
C pH, tot. K and tot Na in aqueous phase 
         if (irank.eq.root) write(25,2300)time
     *               ,(bn(ii,iortx(1)),ii=1,m1)
     *               ,(cn(ii,iortx(1)),ii=1,m2)
     *               ,(pn(ii,iortx(1)),ii=1,m3)
     *               ,cs(iortx(1)),por(iortx(1))
     *               ,pHarr(iortx(1))

        if (irank.eq.root) write(11,2300)time
     *               ,(bn(ii,iortx(2)),ii=1,m1)
     *               ,(cn(ii,iortx(2)),ii=1,m2)
     *               ,(pn(ii,iortx(2)),ii=1,m3)
     *               ,cs(iortx(2)),por(iortx(2))
     *               ,pHarr(iortx(2))

        if (irank.eq.root) write(12,2300)time
     *               ,(bn(ii,iortx(3)),ii=1,m1)
     *               ,(cn(ii,iortx(3)),ii=1,m2)
     *               ,(pn(ii,iortx(3)),ii=1,m3)
     *               ,cs(iortx(3)),por(iortx(3))
     *               ,pHarr(iortx(3))

        if (irank.eq.root) write(13,2300)time
     *               ,(bn(ii,iortx(4)),ii=1,m1)
     *               ,(cn(ii,iortx(4)),ii=1,m2)
     *               ,(pn(ii,iortx(4)),ii=1,m3)
     *               ,cs(iortx(4)),por(iortx(4))
     *               ,pHarr(iortx(4))

        if (irank.eq.root) write(14,2300)time
     *               ,(bn(ii,iortx(5)),ii=1,m1)
     *               ,(cn(ii,iortx(5)),ii=1,m2)
     *               ,(pn(ii,iortx(5)),ii=1,m3)
     *               ,cs(iortx(5)),por(iortx(5))
     *               ,pHarr(iortx(5))
 2300    format(70(1x,e10.3))
c kg44  here the two very long if are finished... 
c      endif
c      endif  

c  *********************************************************
c  write out the species concentrations
c  *********************************************************
c**** output fuer t als backup
      if (irank.eq.root)then
         if (mod(itimestep_tp,it_out).eq.0) then
      call writedump(bn,cn,pn,por,dm,
     *     time,itimestep_tp,tprint,k1,m1,m2,m3,nxmax)
#ifdef __GNU
	idum=vtkout(%val(itimestep_tp),%val(time),%val(nxmax),%val(m1),
     &     %val(m2),%val(m3),dx,bn,cn,pn,por,eh,pHarr,dumb, dumc, dump)
#else
	idum=vtkout(itimestep_tp,time,nxmax,m1,
     &     m2,m3,dx,bn,cn,pn,por,eh,pHarr,dumb,dumc,dump)
#endif

         endif                          ! write out backso
      endif 
c**** end output fuer t als backup
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


c  *********************************************************
c  subroutine writes out the species concentrations dm por as backup
c  to file backup.dat
c  *********************************************************

       subroutine writedump(bn,cn,pn,por,dm,
     *     time,itimestep_tp,tprint,k1,m1,m2,m3,nxmax)
      include 'gwheader.inc'
        character *10 dumb,dumc,dump

        double precision bn(nbasis,nnodex+2),cn(ncompl,nnodex+2)
     *, pn(nsolid,nnodex+2),por(nnodex+2),dm(nnodex+2)
     * ,time,tprint(25),dtmax
        integer itimestep_tp,k1,m1,m2,m3,nxmax
        common /inc/ dumb(nbasis),dumc(ncompl),dump(nsolid)

c****  write out concentration arrays
c****   solutes (basis species and complexes)
        open(17,file="backup.dat",form='formatted',status='unknown')
        do 1722 j1=1,m1
        write(17,*)(bn(j1,nix),nix=1,nxmax)
        write(17,*)dumb(j1),'above'
 1722   continue
        do 1720 j1=1,m2
        write(17,*)(cn(j1,nix),nix=1,nxmax)
        write(17,*)dumc(j1),'above'
 1720   continue
c****  solids and porosity
        do 1716 j1=1,m3
        write(17,*)(pn(j1,nix),nix=1,nxmax)
        write(17,*)dump(j1),'above  '
 1716   continue
 1776   format(85(1x,e10.4))
        write(17,*)(por(nix),nix=1,nxmax+2)
        write(17,*)'porosity above  '
        write(17,*)(dm(nix),nix=1,nxmax+2)
 1775   format(85(1x,e10.4))
        write(17,*)'diffusion coefficient above  '
        write(17,*)time,itimestep_tp,tprint(k1),k1
        close(17)

        end

c  *********************************************************
c  subroutine reads  the species concentrations dm por from backup
c  file backup.dat
c  *********************************************************

       subroutine readdump(bn,cn,pn,por,dm,
     *         time,itimestep_tp,tprint,k1,m1,m2,m3,nxmax)
      include 'gwheader.inc'
        character *10 dumb,dumc,dump
        double precision bn(nbasis,nnodex+2),cn(ncompl,nnodex+2)
     *,pn(nsolid,nnodex+2),por(nnodex+2),dm(nnodex+2)
     * ,time,tprint(25),dtmax
        integer itimestep_tp,k1,m1,m2,m3,nxmax
        common /inc/ dumb(nbasis),dumc(ncompl),dump(nsolid)

c****  write out concentration arrays
c****   solutes (basis species and complexes)
        open(17,file="backup.dat",form='formatted',status='unknown')
        do 1722 j1=1,m1
        read(17,*)(bn(j1,nix),nix=1,nxmax)
        read(17,*)
 1722   continue
        do 1720 j1=1,m2
        read(17,*)(cn(j1,nix),nix=1,nxmax)
        read(17,*)
 1720   continue
c****  solids and porosity
        do 1716 j1=1,m3
        read(17,*)(pn(j1,nix),nix=1,nxmax)
        read(17,*)
 1716   continue
 1776   format(85(1x,e10.4))
        read(17,*)(por(nix),nix=1,nxmax+2)
        read(17,*)
        read(17,*)(dm(nix),nix=1,nxmax+2)
 1775   format(85(1x,e10.4))
        read(17,*)
        read(17,*)time,itimestep_tp,tprint(k1),k1
        close(17)

        end

