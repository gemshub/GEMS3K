

      program mcotac1D
      implicit double precision (a-h,o-z)
      include 'gwheader.inc'
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
      integer npmax ,nbox, partib(NNODEX)
        END SUBROUTINE partid
        END INTERFACE

      INTERFACE
        subroutine concver(npmax,nbox,dx,bn,cn,partib,partx,partic
     *,ismooth,m1,m2)
        !DEC$ ATTRIBUTES C :: concver
      include 'gwheader.inc'
        double precision xmin,xmax,partx(NUPMAX), dx(NNODEX+2)
        double precision bn(NBASIS,NNODEX),cn(NCOMPL,NNODEX)
      double precision partic(NCOMPL+NBASIS,NUPMAX)
      integer npmax ,nbox, partib(NNODEX), ismooth, m1,m2
        END SUBROUTINE concver
        END INTERFACE

      INTERFACE
       subroutine concneu (npmax,nbox,nxmax,xminr,xmaxr,dx,bn,cn,partib,
     *                    partx,partic,bo,co,ismooth,m1,m2)
        !DEC$ ATTRIBUTES C :: concneu
      include 'gwheader.inc'
        double precision xminr,xmaxr,bn(NBASIS,NNODEX),cn(NCOMPL,NNODEX)
      double precision partx(NUPMAX), dx(NNODEX+2),bo(NBASIS,NNODEX)
      double precision co(NCOMPL,NNODEX),partic(NCOMPL+NBASIS,NUPMAX)
        integer nxmax,npmax ,nbox, partib(NNODEX), ismooth, m1,m2
        END SUBROUTINE concneu
        END INTERFACE

      INTERFACE
        subroutine hydro1d(nxmax,h0,hb,tx,am,st,por,ir,qw,qbil,text,vx
     *                    ,dx,icyc,texe,time,fname)

        !DEC$ ATTRIBUTES C :: hydro1d
      include 'gwheader.inc'
        double precision h0(NNODEX),hb(NNODEX),tx(NNODEX),am(NNODEX)
        double precision st(NNODEX),por(NNODEX),qw(NNODEX),qbil(NNODEX)
        double precision dx(NNODEX+2), vx(NNODEX+2),texe,time
        integer nxmax, icyc, ir(NNODEX)
        character*10 text, fname
        END SUBROUTINE hydro1d
        END INTERFACE

      INTERFACE
        subroutine walk2(npmax,nxmax,ncyc,along,aquer,dm,texe,dx,vx
     *,partx,partxo,xmaxr,xminr,partic,bn,cn,partib,ibpstart,x,bo,co,m1
     *,m2)

        !DEC$ ATTRIBUTES C :: walk2
      include 'gwheader.inc'
        double precision xminr,xmaxr,bn(NBASIS,NNODEX),cn(NCOMPL,NNODEX)
      double precision dx(NNODEX+2),bo(NBASIS,NNODEX)
      double precision co(NCOMPL,NNODEX),partic(NCOMPL+NBASIS,NUPMAX)
        double precision along, aquer, dm(nnodex+2),texe
      double precision vx(NNODEX+2),partx(NUPMAX),partxo(NUPMAX)
        double precision x(NNODEX)
        integer nxmax,npmax ,ncyc, partib(NNODEX), ismooth, m1,m2,m3,m4
        integer ibpstart
        END SUBROUTINE walk2
        END INTERFACE


c<<<<<<<<<<<<<<<FROM GEMS integration<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c MAIN FORTRAN PROGRAM START IS HERE

      double precision xxyy , xarray(10)

	INTEGER argc, iinn, i_gems, nNodes
	integer CSTR(100)
	integer c_to_i(100)
c kg44 itergemstime only needed for debug
c	integer itergemstime(250,51)    ! array for output of iterations done per node during every time step
      integer nodeTypes(51)
c for F_GEM_CALC_NODE      
	integer iNodeF
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
	character*100 CSTR_char30
	character*20 dummystring(25)
	character*100 gems_in_ipmf,gems_dbr_f1,gems_dbr_f2
	character*100 gems_dbr_w
	character*20 dummystringb(20)
c in an ideal world this definitions are only used for MPI stuff
#ifdef __MPI
        integer ierr
#endif
	
c12345678901234567890123456789012345678901234567890123456789012345690
      integer itest,ncyc,nxmax,nymax,isteu,inma,ipfile,ntim
      integer npkt,ir(nnodex+2),npmax,nbox,nboxy,partib(nnodex)   
      integer ismooth,iortx(5), i_sorb,j_sorb,iche(nnodex+2),j_decay
     *,ialkali,icyc
      integer npin, s,ss
	integer p_nICb, p_nPHb  
      integer p_nDCb  
      integer itimestep_tp, dtprstep, k1
      real t1,t2
      real p_A_trans(17,7)          ! GEMS invert coeff.-matrix
      real p_A(7,17)          ! GEMS invert coeff.-matrix
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
c !// Dynamic data - dimensions see in DATACH.H and DATAMT.H structures
c !// exchange of values occurs through lists of indices, e.g. xDC, xPH
      double precision  p_xDC(17) ! (nDCb)  !  // DC mole amounts at equilibrium [nDCb]      -      -      +     +
      double precision  p_gam(17) ! (nDCb)  !  // activity coeffs of DC [nDCb]               -      -      +     +
      double precision  p_aPH(20) ! (nPHb)  !// Specific surface areas of phases (m2/g)       +      +      -     -
      double precision  p_xPH(20) ! (nPHb)  !// total mole amounts of phases [nPHb]          -      -      +     +
      double precision  p_vPS(20) ! (nPSb)  !// phase volume, cm3/mol        [nPSb]          -      -      +     +
      double precision  p_mPS(20) ! (nPSb)  !// phase (carrier) mass, g      [nPSb]          -      -      +     +
      double precision  p_bPS(7,20) ! (nICBb,nPSb)  !// bulk compositions of phases  [nPSb][nICb]    -      -      +     +
      double precision  p_xPA(20) ! (nPSb)  !// amount of carrier in phases  [nPSb] ??       -      -      +     +
      double precision  p_dul(17) ! (nDCb)  ! // upper kinetic restrictions [nDCb]           +      +      -     -
      double precision  p_dll(17) ! (nDCb)  ! //  lower kinetic restrictions [nDCb]           +      +      -     -
      double precision  p_bIC(7) ! (nICb)  !// bulk mole amounts of IC[nICb]                +      +      -     -
      double precision  p_rMB(7) ! (nICb)  !// MB Residuals from GEM IPM [nICb]             -      -      +     +
      double precision  p_uIC(7) ! (nICb)  !// IC chemical potentials (mol/mol)[nICb]       -      -      +     +


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

        double precision amin(nsolid,nnodex),km(nsolid)
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

        double precision vout(nsolid,nnodex),rout(nsolid,nnodex)
        double precision wout(nsolid,nnodex)

        integer ikin(nsolid),ifg(nsolid)
        integer ifgp(nsolid),ifgd(nsolid)
ckinet>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c
      dimension dx(nnodex+2),am(nnodex+2),por1(nnodex+2),poro(nnodex+2)
     * ,partx(nupmax),partxo(nupmax),partic(nupmax,nbasis+ncompl)
     * ,dm(nnodex+2)
     * ,vx(nnodex+2), partiv(nnodex)
     * ,hb(nnodex+2), por_null(nnodex+2), tx_null(nnodex+2)
     * ,h0(nnodex+2), st(nnodex+2), qw(nnodex+2), qr(nnodex+2)
     * ,qbil(nnodex+2),	tx(nnodex+2)  
     * ,bnflow(nbasis,nnodex+2),cnflow(ncompl,nnodex+2)
     * ,sumbflowt(nbasis,nnodex+2)
     * ,sumcflowt(ncompl,nnodex+2)
     * ,sumqwat(nnodex+2)
     * ,dqwater(nnodex+2)
     * ,sumbcflowt(ncompl,nnodex+2)
	 

      dimension x(nnodex),tmp(nnodex),bn(nbasis,nnodex)
     1,bo(nbasis,nnodex),pn(nsolid,nnodex),po(nsolid,nnodex)
     1,cn(ncompl,nnodex),co(ncompl,nnodex)
     2,eqconst(ncompl+nsolid,nnodex),acb(nbasis,nnodex)
     3,acc(ncompl,nnodex),q(nsolid,nnodex),tmpold(nnodex),
     5z(nbasis+nsolid,nbasis+nsolid),re(nbasis+nsolid)
     7,wconst(nbasis),bc2(nbasis),eh(nnodex),pHarr(nnodex)
     8,pnw(nsolid),pnd(nsolid)
     9,etc(3)
     9,dpn(nsolid,nnodex)
     5,index_i1(nbasis),bi_i1(nbasis),gesb_i1(nbasis)
     6,index_i2(nbasis),bi_i2(nbasis),gesb_i2(nbasis)
     7,gesp_i1(nsolid),gesp_i2(nsolid)
     8,gesc_i1(ncompl),gesc_i2(ncompl)
     5,gespvf_i1(nsolid),gespvf_i2(nsolid)

      dimension cs(nnodex),tmpk(nnodex)
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

c<<<<<<  system time initialisation for CPU consumption purposes
      time_initial=secnds(0.)
      write(*,*)time_initial
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
	nNodes=51
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
      ipor=0
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
      write (*,*)gespi

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
	write(*,*) 'fname = ',fname,nxmax
c      call holdat1d(nxmax,"iche01.dat",hb,text)
#ifdef __GNU
      call holdat1d(%val(nxmax),fname // char(0),hb)
#else
      call holdat1d(nxmax,fname // char(0),hb)
#endif
      do 1331 ih=1,nxmax
      iche(ih)=int(hb(ih))
      nodeTypes(ih) = iche(ih) ! 1
1331  continue
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
      read(20,'(40x,5i3)')(iortx(iort),iort=1,5) 
      close(20)
      write(*,'(a22,5(1x,i4))')'breakthrough at nodes '
     *         ,(iortx(iort),iort=1,5)
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
        write(*,*)(dtdumb(m), m=1,m1)
        write(*,*)(tdumb(m), m=1,m1)
        write(*,*)(dtdumc(m), m=1,m2)
cpause        pause

c
c  ***************************************
c  assign initial temperature and porosity to each node
c  ***************************************
   43 do 42 n=1,nxmax+2
      por(n)=0.32                  !for initial porcalc 
      poro(n)=0.32       
c      por(n)=phi              ! kinet
      por(n)=1              ! kinet
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

c   input for dynamic calculations
      call inparf(itest,ncyc,nxmax,isteu,inma,
     *ipfile,ntim,npin,npkt,ismooth,i_sorb,j_sorb,j_decay
     *,backg,rd
     *,xlambda,aquer,along,vxx,dm0,dtmax,tmult,dx,de)
      Write(*,*)itest,ncyc,nxmax,isteu,inma,ipfile
     *,ntim,npin,npkt,ismooth,i_sorb,j_sorb,j_decay
     *,backg,rd
     *,xlambda,aquer,along,vxx,dm0,dtmax,tmult,dx,de
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

      call initial(itemp,li,lb,li_i1,li_i2,in1,in2,lnh,err,nxmax,tb
     +,temp,vo
     1,x,itype,num,con,eqconst,tmp,bi,bc,indexi,indexb
     +,index_i1,bi_i1,gesb_i1,index_i2,bi_i2,gesb_i2
     +,gesp_i1,gesp_i2,gesc_i1,gesc_i2
     +,gespvf_i1,gespvf_i2
     2,q,acb,acc,bn,pn,cn,vjb,vjc,s,ss,lne,eh,idismdl,cs,tmpk
     3,iche,i_sorb,ialkali,dumb,dumc,dump,itmpdep)

c  *******************************************************************
c  initial porosity calculation from solids concentration
c  *******************************************************************
       do 1315, nspezx=1,nxmax
       if (iche(nspezx).eq.2)then
            if(m3.gt.0.and.ipor.gt.0)then                                    !  .and.ipor.gt.0
            call porcalc(nspezx,m3,pn,pnw,pnd,etc,por
     *         ,cn,i_sorb,xnaohmw,xnaohd,xkohmw,xkohd)
            endif
       endif
 1315  continue
       write(*,'(85(f5.3,1x))')(por(nspezx),nspezx=1,nxmax)
cpause       pause
      write(*,*)(pn(1,nx),nx=1,nxmax)
cpause	pause
c      t2=secnds(t1)
c      write (6,2200) t2

c>>>>>02-2003 modified boundary on the right side
      if(imodbound.gt.0)then  
       do 1316, nspezx=1,nxmax+1
       if (nspezx.ge.45)then
         bn(4,nspezx)=0.
         bn(5,nspezx)=0.
       endif
 1316  continue
      endif
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
      write(*,*)ir
cpause	pause
#ifdef __GNU
      call holdat1d(%val(nxmax+2),"ss0001.dat"//char(0),st)
#else
      call holdat1d(nxmax+2,"ss0001.dat"//char(0),st)
#endif
      write(*,*)st
#ifdef __GNU
	call holdat1d(%val(nxmax+2),"por001.dat"//char(0),por)
#else
	call holdat1d(nxmax+2,"por001.dat"//char(0),por)
#endif
      write(*,*)por
cpause	pause
	do 3331 ih=1,nxmax+2
      por_null(ih)=por(ih)
c2003      tx_null(ih)= 1.28E-10*(1.-por(ih))**2/por(ih)**3.       !exp 4 specific
      tx_null(ih)= 1.28E-10*(1.-por(ih))**2/por(ih)**3.       !exp 4 specific
      tx(ih)= tx_null(ih)*por(ih)**3/(1.-por(ih))**2
 3331 continue 
#ifdef __GNU
      call holdat1d(%val(nxmax+2),"qr0001.dat"//char(0),qr)
#else
      call holdat1d(nxmax+2,"qr0001.dat"//char(0),qr)
#endif
      write(*,*)qr
cpause	pause
#ifdef __GNU
      call holdat1d(%val(nxmax+2),"qn0001.dat"//char(0),qw)
#else
      call holdat1d(nxmax+2,"qn0001.dat"//char(0),qw)
#endif
       write(*,*)qw
cpause	pause
#ifdef __GNU
      call holdat1d(%val(nxmax+2),"am0001.dat"//char(0),am)
#else
      call holdat1d(nxmax+2,"am0001.dat"//char(0),am)
#endif
      write(*,*)am
cpause	pause
#ifdef __GNU
      call holdat1d(%val(nxmax+2),"h00001.dat"//char(0),h0)
#else
      call holdat1d(nxmax+2,"h00001.dat"//char(0),h0)
#endif
      write(*,*)h0
cpause	pause
      do 3332 ih=1,nxmax+2
         hb(ih)=h0(ih)
         qw(ih)=qw(ih)+qr(ih)
c      write(*,*)'main', ih,tx(ih)
 3332  continue

c03      endif                                          !endif ihydro = 1



cpause      pause
c04<<<<<<<<<<<
      p_NodeHandle=1
c     open data bridge file initially for initialising the spatial distribution of chemical systems
c  first read is for boundary conditons node 1
      do 45 i=1,10
	CSTR(i)=96+i
  45	xarray(i)=i**i
	argc=5
	iinn=3
      xxyy=1.23456789
	FNAME10="abcdevwxyz"
	CSTR_char30="xbtdefghijklmnopqrst1234567890"
      read(CSTR_char30,'(100(a1))')(chch(i),i=1,100)

      write(*,*)gems_in_ipmf
      write(*,*)gems_dbr_f1
      write(*,*)gems_dbr_f2

      write (*,*)'FORTRAN defined in C++ argc', argc
      write (*,*)'FORTRAN integer        iinn', iinn
      write (*,*)'FORTRAN double         xxyy', xxyy
	write (*,*)'FORTRAN char*10   FNAME10  ', fname10
      write (*,*)'FORTRAN double array (10)  ', xarray
	write (*,*)'FORTRAN int array CSTR(10) ', (CSTR(L),L=1,10)
	write (*,*)'FORTRAN char*10     line   ', line
	write (*,*)'FORTRAN char*10 c_to_i ',  c_to_i
      

      write(*,*)'nnode',nNodes
	write(*,*)'cto',c_to_i1
	write(*,*)'cto',c_to_i2
	write(*,*)'nodetype',nodeTypes

c	pause "F_GEM_INIT"

       nNodes= nxmax-1    !   1 
      call F_GEM_INIT( gems_in_ipmf )
      write(*,*)'nNodes =', nNodes

c	pause "F_GEM_INIT called"

      do 50 i=1,20
	  chch(i)=char(CSTR(i))
  50	continue
      write(CSTR_char30,'(30a1)')(CSTR(i),i=1,30) 


c<<<<<<<<<<<<<<<FROM GEMS integration<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
cccc -t "cal-dol-boun1-dch.dat", "cal-dol-boun1-dbr-1.dat"  are also in ipmfiles-dat.lst normally
c     open data-ch file for initialisation and names (only once)

c   read 2. for p_ variables        
c        if(igems_rw.eq.1)
      call F_GEM_GET_DCH( p_nICb, p_nDCb, p_nPHb, p_A )

      write(*,*) 'gemsA', p_A
c	pause "F_GEM_GET_DCH nach A "

      write(*,*) 'gemsread_dch.dat done'

	p_A_trans = transpose (p_A)               !invert stoichiometric coeff. matrix
      write(*,*) 'gemsA', p_A_trans
c	pause "F_GEM_GET_DCH nach A_trans "

      
cc2005      endif  ! igems_rw=0


      iNode=0
      p_NodeHandle=1
      p_NodeStatusFMT = 1
      p_NodeStatusCH=1
      if(i_gems.eq.1) then
	call F_GEM_READ_NODE( gems_dbr_f1, p_NodeHandle,p_NodeTypeHY
     *,p_NodeTypeMT,p_NodeStatusFMT,p_NodeStatusCH,p_IterDone
     *,p_T, p_P,p_Vs,p_Vi,p_Ms,p_Mi,p_Gs
     *,p_Hs,p_Hi,p_IC,p_pH,p_pe,p_Eh
c     *,p_Tm,p_dt,p_dt1,p_Vt,p_vp,p_eps,p_Km,p_Kf,p_S,p_Tr,p_h,p_rho,p_al
c     *,p_at
c     *,p_av,p_hDl,p_hDt,p_hDv,p_nto 
     *,p_bIC,p_rMB,p_uIC,p_xDC,p_gam
     *,p_dul, p_dll, p_aPH,p_xPH, p_vPS,p_mPS,p_bPS,p_xPA
     *)
      write(*,*)'iNode,p_NodeHandle,p_NodeTypeHY,p_NodeTypeMT
     *,p_NodeStatusFMT,p_NodeStatusCH,p_IterDone'
      write(*,*)iNode,p_NodeHandle,p_NodeTypeHY,p_NodeTypeMT
     *,p_NodeStatusFMT,p_NodeStatusCH,p_IterDone,p_T, p_P
c      pause

      write(*,*)'p_T, p_P,p_Vs,p_Vi,p_Ms,p_Mi,p_Gs'
      write(*,*)p_T, p_P,p_Vs,p_Vi,p_Ms,p_Mi,p_Gs
c	pause
     *
      write(*,*)p_Hs,p_Hi,p_IC,p_pH,p_pe,p_Eh,
     * p_Tm
     *,p_dt,p_dt1,p_ot,p_Vt,p_eps,p_Km,p_Kf,p_S,p_Tr,p_h,p_rho,p_al,p_at
     *,p_av,p_hDl,p_hDt,p_hDv,p_nPe,p_xDC,p_gam,p_xPH,p_vPS,p_mPS,p_bPS
     *,p_xPA,p_bIC,p_rMB,p_uIC,p_dRes1,p_dRes2
c      pause

c  here loop to read in additional geochmical systems to define other nodes....
c      aa-initial-dbr-1.dat to aa-initial-dbr-50.dat and aa-boundary-dbr-0.dat are available)
c  <<<<<<<   tranfer of GEMS nomenclature to MCOTAC naming
c     bn=
c     cn=
c     pn=
      if(i_output.eq.1)open(35, file='mco_out.out')

      do 1690 n=1,1
	do 1691 ib=1,m1-1    !charge is last parameter in the lict of bn  24.01.2005 but not transported
ccc	bn(ib,n)=gemsxDc(i_bcp_gemx(ib))         !2)
        write(*,*) ' debug', i_bcp_gemx
        write(*,*) ' debug', ib, i_bcp_gemx(ib)
	bn(ib,n)=p_xDc(i_bcp_gemx(ib))         !2)
 1691 continue
	do 1692 ic=1,m2
ccc	cn(ic,n)=gemsxDc(i_bcp_gemx(m1+ic))       !10)
	cn(ic,n)=p_xDc(i_bcp_gemx(m1+ic))       !10)
 1692 continue
	do 1693 ip=1,m3
ccc	pn(ip,n)=gemsxDc(i_bcp_gemx(m1+m2+ip))     ! 12)/1.
	pn(ip,n)=p_xDc(i_bcp_gemx(m1+m2+ip))     ! 12)/1.
c	pn(2,n)=gemsxDc(13)/1.
 1693 continue
 1690 continue
      write(*,'(13(e8.2,1x))')(bn(ib,1),ib=1,m1),(cn(ic,1),ic=1,m2)
     *,(pn(ip,1),ip=1,m3)
ccc      write(*,'(13(e8.2,1x))')(gemsxDc(ib),ib=1,gemsnDCb)
      write(*,*)' 1 p_xDc(ib) '
	write(*,'(13(e8.2,1x))')(p_xDc(ib),ib=1,p_nDCb)

c      write(*,*)(bn(ib,1),ib=1,m1)
c      write(*,*)(cn(ic,1),ic=1,m2)
c      write(*,*)(pn(ip,1),ip=1,m3)
c	pause "basis complex solids at n=1"
      endif

c     open data bridge file initially for initialising the spatial distribution of chemical systems
c  second read is for initial conditons nodes 2 to nxmax
      
      if(i_gems.eq.1) then
      iNode=0
      p_NodeStatusFMT = 1
      p_NodeStatusCH=1
      if(i_gems.eq.1) then
	call F_GEM_READ_NODE( gems_dbr_f2, p_NodeHandle,p_NodeTypeHY
     *,p_NodeTypeMT,p_NodeStatusFMT,p_NodeStatusCH,p_IterDone
     *,p_T, p_P,p_Vs,p_Vi,p_Ms,p_Mi,p_Gs,p_Hs,p_Hi,p_IC,p_pH,p_pe,p_Eh
c     *,p_Tm,p_dt,p_dt1,p_Vt,p_vp
c     *,p_eps,p_Km,p_Kf,p_S,p_Tr,p_h,p_rho,p_al,p_at
c     *,p_av,p_hDl,p_hDt,p_hDv,p_nto
     *,p_bIC,p_rMB,p_uIC, p_xDC,p_gam
     *,p_dul, p_dll, p_aPH,p_xPH, p_vPS,p_mPS,p_bPS,p_xPA
     *)
	 write(*,*) 'inode= ' ,inode
      if (i_output.eq.1)then
	   write(35,*) 'node','2 and rest' 
	   write(35,'(20(e18.12,1x))')(p_xDc(ib),ib=1,p_nDCb)
	   write(35,'(10(e18.12,1x))')(p_bIC(ib),ib=1,p_nICb)
         write(35,*) 'end node 2 and rest'
      endif
	 endif

      do 1695 n=2,nxmax-1
	do 1696 ib=1,m1-1
	bn(ib,n)=p_xDc(i_bcp_gemx(ib))    !   2)   ! i_bcp_gemx(1)=2
 1696 continue
	do 1697 ic=1,m2
	cn(ic,n)=p_xDc(i_bcp_gemx(m1+ic))   ! i_bcp_gemx(7)=10
 1697 continue
	do 1698 ip=1,m3
	pn(ip,n)=p_xDc(i_bcp_gemx(m1+m2+ip))     !12)/1.    ! i_bcp_gemx(12)=12
 1698 continue
 1695 continue
      if(i_output.eq.1)then
       write(35,*)' bn n=2  '
       write(35,'(13(e8.2,1x))')(bn(ib,2),ib=1,m1)
       write(35,*)' cn n=2  '
       write(35,'(13(e8.2,1x))')(cn(ic,2),ic=1,m2)
       write(35,*)' pn n=2  '
       write(35,'(13(e8.2,1x))')(pn(ip,2),ip=1,m3)

       write(35,*)' 2 p_xDc(ib) '
       write(35,'(13(e8.2,1x))')(p_xDc(ib),ib=1,p_nDCb)
	 write (35,*)' end of initial calc ' 
      endif
      write(*,'(13(e8.2,1x))')(bn(ib,2),ib=1,m1),(cn(ic,2),ic=1,m2)
     *,(pn(ip,2),ip=1,m3)
ccc      write(*,'(13(e8.2,1x))')(gemsxDc(ib),ib=1,gemsnDCb)
      write(*,*)' 2 p_xDc(ib) '
      write(*,'(13(e8.2,1x))')(p_xDc(ib),ib=1,p_nDCb)

c      write(*,*)(cn(ic,2),ic=1,m2)
c      write(*,*)(pn(ip,2),ip=1,m3)
c	pause "basis complex solids at n=2 to nxmax"

cgems	stop



      do 51 i=1,20
	  chch(i)=char(c_to_i(i))
  51	continue
      endif 
cc      write(CSTR_char30,'(30a1)')(chch(i),i=1,30) 
cc	write (*,*)'FORTRAN char*10 c_to_i ',c_to_i,CSTR_char30
cpause	pause "pause"

c>>>>>>>>>>>>>>>>>>>>FROM GEMS integration>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>



      if (idynam.eq.0) go to 500



c i_sorb for sorption as a complexation
      if(i_sorb.gt.0)then
         write(*,*)' Sorption for complexes greater 
     *   than No.',i_sorb,'assumed' 
cpause         pause
      endif
      if(j_sorb.gt.0)then
         write(*,*)' Sorption for basis species greater 
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
       write(*,*)'dx',dx  
      do 1322 ix=1,nxmax+2
       vx(ix) = vxx
       dm(ix) = dm0
       write(*,*)'vx',ix,iy,vx(ix)  
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
      write(*,*)'ibpstart =',ibpstart, npmax, nbox
      de=de*tmult
c      write(*,*)'npbox xmin xmax',nbox, xmin,xmax,xminr,xmaxr
c       pause
c   set particles in the grid
#ifdef __GNU
      write(*,*)'before setpar',npmax,xmin,xmax,nbox
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
      write(*,*)'hb',hb
cpause	pause
      write(*,*)'vx',vx
cpause	pause
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
      endif                                       ! endif ihydro=1
c>>>>>>>>140895    END HYDROLOGY,  initially 
c       stop

c <<<<< nov 2002



c***  input of concentration arrays at a certain time 
      write(*,*)'do you want to start calculation at a certain time?'
      write(*,*)'(Y/N)'
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
      OPEN(10,file='conc1t.dat',form='formatted',status='unknown')
      OPEN(11,file='conc2t.dat',form='formatted',status='unknown')
      OPEN(12,file='conc3t.dat',form='formatted',status='unknown')
      OPEN(13,file='conc4t.dat',form='formatted',status='unknown')
      OPEN(14,file='conc5t.dat',form='formatted',status='unknown')
      write(10,1111)' time     ',(dumb(iw1),iw1=1,m1),
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

      texe=dtmax

 235  continue                                              !  next time step interval
      itimestep_tp=itimestep_tp+1
      write(*,'(a4,1x,i3,1x,6(e8.2,1x))')'235 ',itimestep_tp,
     *bn(1,1),pn(1,1),pn(2,1),bn(1,2),pn(1,2),pn(2,2)
      texe=dtmax
c      write(*,*)'DIM F90', nnodex,nbasis,ncompl,nsolid,nupmax,npmax
      told=time
      time=time+texe
c      time=time+dtmax

      write(*,*)'time',time,told,texe,tprint(k1),k1
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
      if(imodbound.gt.0) then
       do 1317, nspezx=1,nxmax+1
       if (nspezx.ge.41)then
         bn(4,nspezx)=0.
         bn(5,nspezx)=0.
       endif
 1317  continue
      endif
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
        bo(j,n)=bn(j,n)         
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
	  if (i_output.eq.1.and.n.eq.2)then
      	write(35,*) 'node',n,'vor transport' 
	     write(35,*) 'DCb'
	    write(35,'(20(e12.6,1x))')(p_xDc(ib),ib=1,p_nDCb)
        endif
  240 continue
c2002<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c2003         if(ipor.gt.0)then
c2003          if (k1.gt.3)then
c2003           dtnminmin=dtmax
c2003           dtdxmin=xmax
c2003           ndtmin=nxmax
c2003           vmin=vx(1)
c2003           do 221 n1=2,nxmax-1
c2003            if(vx(n1).ne.0.)dtnmin1=dx(n1)/vx(n1)/2.
c2003            if(dm(n1).ne.0.)dtnmin2=dx(n1)**2/2/dm(n1)
c2003             dtnmin=dtnmin1+dtnmin2
c2003            if(dtnminmin.le.dtnmin)then
c2003               dtnminmin=dtnmin
c2003               dtdxmin=dx(n1)
c2003               vmin=vx(n1) 
c2003               ndtmin=n1
c2003            endif
c2003  221      continue
c2003           dtmax=dtnminmin/2.
c          texe=dtnminmin
c            dtmax=dx(2)/2./vx(2)
c            texe=dx(2)/2./vx(2)
c2003          endif
c          write(*,*)'k3-ddtmax=',k1, texe,dtdxmin,vmin,ndtmin
c         pause
c2003         endif
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
      write (*,*) 'walk time:  treal texe', treal, texe
#ifdef __GNU
      call walk2(%val(npmax),%val(nxmax),%val(ncyc),%val(along),
     * %val(aquer),dm,%val(texe),dx,vx
     *  ,partx,partxo, %val(xmaxr),%val(xminr),partic,bn,cn,partib,
     *  %val(ibpstart),x,bo,co,%val(m1),%val(m2))
#else
      call walk2(npmax,nxmax,ncyc,along,aquer,dm,texe,dx,vx
     *  ,partx,partxo, xmaxr,xminr,partic,bn,cn,partib,
     *  ibpstart,x,bo,co,m1,m2)
#endif
c**assign concentrations at t+dt to grid  (including boundary conditions)
c**particle in which nbox

c      write (*,*)'xminr xmaxr',xminr,xmaxr
c      write(*,'(6(e10.4,1x))')(bn(i,1),i=1,m1)
c      write(*,'(6(e10.4,1x))')(bn(i,2),i=1,m1)
c      write(*,'(6(e10.4,1x))')(bn(i,3),i=1,m1)
c      write(*,'(6(e10.4,1x))')(bn(i,11),i=1,m1)
c       pause
c**new concentration in each box
#ifdef __GNU
      call concneu(%val(npmax),%val(nbox),%val(nxmax),%val(xminr),
     *  %val(xmaxr),dx,bn,cn,partib,partx,
     *  partic,bo,co,%val(ismooth),%val(m1),%val(m2))
#else
      call concneu(npmax,nbox,nxmax,xminr,
     *  xmaxr,dx,bn,cn,partib,partx,
     *  partic,bo,co,ismooth,m1,m2)
#endif
c kg44 changed output format
      write(*,*)'cneu',itimestep_tp,
     *bn(1,1),pn(2,1),pn(3,1)
      write(*,*)'cneu',itimestep_tp,
     *bn(1,2),pn(2,2),pn(3,2)

c      write(*,'(a4,1x,i3,1x,6(e12.6,1x))')'cneu',itimestep_tp,
c     *bn(1,1),pn(2,1),pn(3,1)
c      write(*,'(a4,1x,i3,1x,6(e12.6,1x))')'cneu',itimestep_tp,
c     *bn(1,2),pn(2,2),pn(3,2)




c>>>>>02-2003 modified boundary on the right side
      if(imodbound.gt.0)then  
       do 1318, nspezx=1,nxmax+1
       if (nspezx.ge.45)then
         bn(4,nspezx)=0.
         bn(5,nspezx)=0.
       endif
 1318  continue
      endif 
c>>>>>02-2003 modified boundary on the right side

c********  decay and chain decay
c********  decay
       do 1444 n=1, nxmax
       if(j_decay.gt.0) then
          bn(j_decay,n)=bn(j_decay,n)*dexp((-xlambda*texe))
       endif

c********  chain decay
       if(j_chaindecay.eq.1) then
c       goto 1234       

       co_m2_2=co(m2-2,n)
       co_m2_1=co(m2-1,n)
       co_m2_0=co(m2,n)
       bn_m2_2=bo(m2-2,n)
       bn_m2_1=bo(m2-1,n)
       bn_m2_0=bo(m2,n)

       dbndec1=bn_m2_2*(1.-dexp(-xlambda1*texe))
       co(m2-2,n)=co(m2-2,n)-dbndec1
       if(bn(m1-2,n).le.zero) bn(m1-2,n)=0.          !zero was 1.e-40

       dcndec1=co(m2-2,n)*(1.-dexp(-xlambda1*texe))
       co(m2-2,n)=co(m2-2,n)-dcndec1
       if(co(m2-2,n).le.zero) co(m2-2,n)=0.          !zero was 1.e-40

       dbndec2=bn_m2_1-
     *        ((xlambda1*bn_m2_2*dexp(-texe*(xlambda1
     *         -xlambda2)) / (-xlambda1+xlambda2)
     *        -(xlambda1*bn_m2_2+bn_m2_1*xlambda1
     *        -bn_m2_1*xlambda2) / (-xlambda1+xlambda2)))
     *                   *dexp(-xlambda2*texe)
       bn(m1-1,n)=bn_m2_1-dbndec2                    !+dcndec1  !+dbndec1   !+dbndec1
       if(bn(m1-1,n).le.zero) bn(m1-1,n)=0.          !zero was 1.e-40

       dcndec2=co(m2-1,n)-
     *        ((xlambda1*co_m2_2*dexp(-texe*(xlambda1
     *         -xlambda2)) / (-xlambda1+xlambda2)
     *        -(xlambda1*co_m2_2+co(m2-1,n)*xlambda1
     *        -co(m2-1,n)*xlambda2) / (-xlambda1+xlambda2)))
     *                   *dexp(-xlambda2*texe)
       co(m2-1,n)=co(m2-1,n) - dcndec2                    !+dcndec1  
       if(co(m2-1,n).le.zero) co(m2-1,n)=0.          !zero was 1.e-40

       dbndec3=bn_m2_0-
     *  (xlambda2*(xlambda1*bn_m2_2*dexp(texe*
     *  (xlambda3-xlambda1))/(xlambda3-xlambda1) 
     *  -xlambda1*bn_m2_2*dexp(texe*(-xlambda2+xlambda3))
     *  /(-xlambda2+xlambda3)-bn_m2_1*xlambda1
     *  *dexp(texe*(-xlambda2+xlambda3))/(-xlambda2+xlambda3)
     *  +bn_m2_1*xlambda2*dexp(texe*(-xlambda2+xlambda3))
     *  /(-xlambda2+xlambda3))/(-xlambda1+xlambda2)+
     *  (xlambda1*bn_m2_2*xlambda2+xlambda2*bn_m2_1
     *  *xlambda1-xlambda2*bn_m2_1*xlambda3-
     *  bn_m2_0*xlambda2*xlambda3+bn_m2_0*xlambda2*xlambda1
     *  +bn(m1,n)*xlambda3**2-bn_m2_0*xlambda1*xlambda3)/
     *  ((xlambda3-xlambda1)*(-xlambda2+xlambda3)))
     *  *dexp(-xlambda3*texe)
       bn(m1,n)=bn_m2_0-dbndec3                !+dcndec2        !dbndec2   !+dbndec2
       if(bn(m1,n).le.zero) bn(m1,n)=0.          !zero was 1.e-40

       dcndec3=co_m2_0-
     *  (xlambda2*(xlambda1*co_m2_2*dexp(texe*
     *  (xlambda3-xlambda1))/(xlambda3-xlambda1) 
     *  -xlambda1*co_m2_2*dexp(texe*(-xlambda2+xlambda3))
     *  /(-xlambda2+xlambda3)-co_m2_1*xlambda1
     *  *dexp(texe*(-xlambda2+xlambda3))/(-xlambda2+xlambda3)
     *  +co_m2_1*xlambda2*dexp(texe*(-xlambda2+xlambda3))
     *  /(-xlambda2+xlambda3))/(-xlambda1+xlambda2)+
     *  (xlambda1*co_m2_2*xlambda2+xlambda2*co_m2_1
     *  *xlambda1-xlambda2*co_m2_1*xlambda3-
     *  co_m2_0*xlambda2*xlambda3+co_m2_0*xlambda2*xlambda1
     *  +co_m2_0*xlambda3**2-co_m2_0*xlambda1*xlambda3)/
     *  ((xlambda3-xlambda1)*(-xlambda2+xlambda3)))
     *  *dexp(-xlambda3*texe)
       co(m2,n)=co_m2_0-dcndec3                !+dcndec2 
       if(co(m2,n).le.zero) co(m2,n)=0.          !zero was 1.e-40

 1234  continue

       goto 1235
       dbndec1=bn(m1-2,n)*xlambda1*texe
       bn(m1-2,n)=bn(m1-2,n)-dbndec1                ! evtl  bn(m1-2,n)=bn(m1-1,n)-dbndec1  ???2004
       if(bn(m1-2,n).le.zero) bn(m1-2,n)=0.          !zero was 1.e-40
       dcndec1=cn(m2-2,n)*xlambda1*texe
       cn(m2-2,n)=cn(m2-2,n)-dcndec1            !error ??  m2-2 ??        =cn(m2-2,n)2004
       if(cn(m2-2,n).le.zero) cn(m2-2,n)=0.          !zero was 1.e-40

       dbndec2=bn(m1-1,n)*xlambda2*texe
       bn(m1-1,n)=bn(m1-1,n)-dbndec2+dbndec1+dcndec1
       if(bn(m1-1,n).le.zero) bn(m1-1,n)=0.          !zero was 1.e-40
       dcndec2=cn(m2-1,n)*xlambda2*texe
       cn(m2-1,n)=cn(m2-1,n)-dcndec2+dcndec1            !error ??  m2-2 ??
       if(cn(m2-1,n).le.zero) cn(m2-1,n)=0.          !zero was 1.e-40

       dbndec3=bn(m1,n)*xlambda3*texe
       bn(m1,n)=bn(m1,n)-dbndec3+dbndec2+dcndec2
       if(bn(m1,n).le.zero) bn(m1,n)=0.

c       goto 1235
       dcndec3=cn(m2,n)*xlambda3*texe
       cn(m21,n)=cn(m2,n)-dcndec3+dcndec2
       if(cn(m2,n).le.zero) cn(m2,n)=0.
       
 1235  continue
       endif
 1444  continue
c>>>>>>>>>>>>>>>>>>>>>>>  chain decay
c      goto 4444


c<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c   here GEMS calculations for each noce inclusive data transfer
c     write dbr files
c<<<<<< Gems databridge files have to be generated in gems format for each node     
c      bn, cn,pn ->  independent element masses for each node new except for boundary nodes
c  calculation bn, cn pn, to xDC (transported/not transportted)
c  then new xDC tarsnfereed to 'bIC' total idepandent masses      

1558  CONTINUE ! NO TRANSPORT DONE BEFORE 

c  <<<<<<<   transfer of GEMS nomenclature to MCOTAC naming
c	do 1699 ib=1,m1
	if (i_gems. eq. 1) then     ! goto 1557
	do 1555 n=2,  nxmax-1                  !node loop for GEMS after Transport step
c      goto 1556  ! only node 2 with old gems  values
	do 1596 ib=1,m1-1

	p_xDc(i_bcp_gemx(ib))=bn(ib,n)   ! 2)  index_gems(1,...,m1,m+1,m2,   m3) index_mcotac(m1+m2+m3)
 1596 continue
	do 1597 ic=1,m2
	p_xDc(i_bcp_gemx(m1+ic))=cn(ic,n)
 1597 continue
	do 1598 ip=1,m3
	p_xDc(i_bcp_gemx(m1+m2+ip))= po(ip,n)  !  gemsxDc(12)    ! pn(1,n)   solids not used for transport
c      if(n.eq.2.and.ip.eq.2)p_xDc(i_bcp_gemx(m1+m2+ip))= po(ip,n)/10.
 1598 continue

c      write(*,*)n,(p_xDc(ib),ib=1,p_nDCb)
c      pause "xDC after transport at n"
      do 1680 ii=1, p_nICb-1    !gems   -1 because last is charge and should be zero
      sum=0.
      do 1681 jj=1, p_nDCb    !gems           
       sum=sum+ p_A(ii,jj) *p_xDc(jj)
 1681 continue

      p_bIC(ii)=sum

c      do 31 i=1, gemsnPSb  
cc      write (*,*) (gemsbPS(i,j), j=1, gemsnICb)   !ICb) 
c  31  continue
cxx	gemsbPS(ii)=sum
c       write(*,*)ii,p_bIC(ii) 
 1680 continue
c      goto 1556  ! only node 2 with old gems  values

c      do 1684 ii=1, p_nICb   !gems
c        sum=0.
c        do 1685 jj=1 , p_nDCb  !  -2  14-02-2005     !gems         

c         sum=sum+ p_A_trans(ii,jj) *p_xDc(jj)
c 1685   continue
c   	  p_bPS(1,ii)=sum            ! only one phse here       ! ( nPSb,nICb) ???
c 1684 continue
 
      do 1682 ii=1, p_nPHb  !gems
        sum=0.
c        do 1683 jj=11,13                              !  gemsnDCb-gemsnPHb+1, gemsnDCb               
cc   	  gemsxPH(ii)=gemsxDc(10+ii)            ! only one phse here       ! ( nPSb,nICb) ???
cfalsch   	  p_xPH(ii)=p_xDc(13+ii)            ! only one phse here       ! ( nPSb,nICb) ???
   	  p_xPH(ii)=p_xDc(p_nDCb-p_nPHb+ii)          !13+ii)            ! only one phse here       ! ( nPSb,nICb) ???
c        write(*,*)ii, gemsxPH(ii)
 1682 continue
cgems      pause "XXXX"
	if (i_output.eq.1.and.n.eq.2)then
      	write(35,*) 'node',n,'vor GEMS' 
	     write(35,*) 'DCb'
	    write(35,'(20(e12.6,1x))')(p_xDc(ib),ib=1,p_nDCb)
	write(35,*) 'ICb'
	    write(35,'(10(e18.12,1x))')(p_bIC(ib),ib=1,p_nICb)
c	     write(35,*) 'b_bPS'
c	    write(35,'(10(e8.2,1x))')(p_bPS(ib),ib=1,p_nDCb)
c	     write(35,*) 'xPH'
c	    write(35,'(10(e8.2,1x))')(p_xPH(ib),ib=1,p_nICb)
	endif


cc      gemsNodeStatusCH =1    ! need GEMS AIA  ??
      p_NodeStatusCH =1    ! need GEMS AIA  
c      p_NodeStatusCH = 5    ! uses PIA (Smart Initial Approximation) to accelerate calculatios
c   array boundaries  * gemsnDCb,gemsnPHb,gemsnPSb,gemsnICb
c                       gemsA(MaxDCN,MaxICN) 
c      write(*,*)n,(gemsbIC(ib),ib=1,gemsnICb)

c	stop


cfalsch 1555 continue   
c      pause "write dbr"
c<<<<<<<<<<GEMS caclculations for all nodes have to be prepared here

ccx      if(i_gems.eq.1)then

cc	do 1411 nspez=2,nxmax-1
 1556 continue

	iNode=  n
      p_NodeHandle=  n
      p_NodeStatusCH= 1    ! 1 : with simplex PIA; 5 smart PIA
      p_NodeStatusFMT = 1
c<<<<<<  system time initialisation for CPU consumption purposes
c      time_gemsstart=RTC()
      time_gemsstart=secnds(0.)

	call F_GEM_CALC_NODE( p_NodeHandle,p_NodeTypeHY,p_NodeTypeMT
     *,p_NodeStatusFMT,p_NodeStatusCH,p_IterDone,p_T, p_P
     *,p_Vs,p_Vi,p_Ms,p_Mi,p_Gs,p_Hs,p_Hi,p_IC,p_pH,p_pe,p_Eh
c     *,p_Tm,p_dt,p_dt1
c     *,p_Vt,p_vp, p_eps,p_Km,p_Kf,p_S,p_Tr,p_h,p_rho,p_al,p_at
c     *,p_av,p_hDl,p_hDt,p_hDv,p_nto
     *,p_bIC,p_rMB,p_uIC,p_xDC,p_gam, p_dul, p_dll, p_aPH
     *,p_xPH,p_vPS,p_mPS,p_bPS,p_xPA
c     *,p_dRes1
     *)

c      time_gemsend=RTC()
      time_gemsend=secnds(0.)
      time_gemstotal=time_gemstotal+(time_gemsend-time_gemsstart)
c      time_gemstotal=time_gemstotal+ secnds(time_gemsstart)

	if (i_output.eq.1.and.n.eq.2)then
      	write(35,*) 'node',n,'nach GEMS' 
	    write(35,*) 'DCb', '#######   ',time_gemstotal
	    write(35,'(20(e12.6,1x))')(p_xDc(ib),ib=1,p_nDCb)
c	    write(35,*) 'ICb'
	    write(35,'(10(e18.12,1x))')(p_bIC(ib),ib=1,p_nICb)
c	write(35,*) 'b_bPS'
c	    write(35,'(10(e8.2,1x))')(p_bPS(ib),ib=1,p_nDCb)
c	write(35,*) 'xPH'
c	    write(35,'(10(e8.2,1x))')(p_xPH(ib),ib=1,p_nICb)
	endif
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
	if (i_output.eq.1.and.n.eq.2)then
      	write(35,*) 'node',n,'weit nach GEMS bei 1798' 
	    write(35,*) 'DCb'
	    write(35,'(20(e12.6,1x))')(p_xDc(ib),ib=1,p_nDCb)
      endif

c kg44 only needed for debug
c      itergemstime(itimestep_tp,n)=p_IterDone

c      write(*,*)itimestep_tp,n,itergemstime(itimestep_tp,n),p_IterDone
c	pause
	p_IterDone=0

 1555 continue                 ! end node loop for GEMS after Transport step             
 1557 continue
      endif        ! i_gems eq.1  
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

        call init2(itemp,li,lb,in1,in2,lnh,err,nspez,tb,
     1     temp,vo,x,itype,num,con,eqconst,tmp,bi,bc,
     2          indexi,indexb,q,acb,acc,bn,pn,cn,vjb,vjc,
     3          s,ss,bc2,lne,eh,idismdl,cs,tmpk,bo,co,po
     4          ,i_sorb,ialkali,dumb,dumc,dump,itmpdep)
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

        do 1390 i=1,m3
        if (pn(i,n).le.(zero)) go to 1330
        if (po(i,n).le.(zero)) write(6,2000)n,
     1  dump(i),treal,tmp(n),eqconst(i,n)
        go to 1390
 1330   if (po(i,n).gt.(zero)) write(6,2050) n,
     1  dump(i),treal,tmp(n),eqconst(i,n)
 1390   continue


	endif                                       ! if i_gems=1 
      write(*,'(a4,1x,i3,1x,6(e8.2,1x))')'n_ge',itimestep_tp,
     *bn(1,1),pn(1,1),pn(2,1),bn(1,2),pn(1,2),pn(2,2)
cgems	 pause "next time step"
c      if(itimestep_tp.gt.10)stop

c  end GEMS reading and re-naming after re-equilibration


c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


c**** output fuer t(k1)
c2002      if (time.eq.tprint(k1)) then
      if (treal.eq.tprint(k1)) then
        kk1=k1+96
        ttt=treal/31557600
        write(*,*)'kk1',kk1,'time',treal
        ch=char(kk1)
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


        do 1513 nn1=1,nxmax
        if(nn1.eq.1)then
                xspez=x(nn1)+dx(1)
       else
               xspez=x(nn1)+dx(1)      !   feb 2003    /2
       endif
      if(lnhc.gt.0.and.acc(lnhc,nn1).gt.0.and.cn(lnhc,nn1).gt.0.)then                  !2003   and.acc(lnhc,nn1).gt.0.and.cn(lnhc,nn1).gt0.
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
         write(*,'(a6,2(e10.2,a6,2x), 2x,a7,i10)')
     *   'time= ',treal,' [sec]',tty,'  [yr]', 'iccyle= ',icyc
         write(*,'(a14,1x,10(e10.4,1x))')
     *   'solids iort(1)',(pn(ii,iortx(1)),ii=1,m3),cs(iortx(1)),
     *    (eqconst(ii,iortx(1)),ii=1,m3)
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
         write(10,2300)treal
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
        write(11,2300)treal
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
        write(12,2300)treal
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
        write(13,2300)treal
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
        write(14,2300)treal
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
C end of file for breakthrough curves 
         xloc=x(1)+iortx(1)*dx(1)
         write(10,*)'at location ',xloc,' [m]'
         close (10)
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
      time_end=secnds(1.)
      time_interval=time_end-time_initial
      write(*,*) 'time_end  time_interval',time_end, time_interval
      write(*,*) 'time steps calculated timestep_tp = ',itimestep_tp
      write(*,*) 'time used during GEMS calling = ',time_gemstotal
      write(*,*) 'time used during GEMS read = ',time_gemsreadtotal
      write(*,*) 'time used during GEMS write = ',time_gemswritetotal
      write(*,*) 'time used during MCOTAC-chem-calc=',time_initreadtotal

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

      if(i_output.eq.1) close(35)


#ifdef __MPI
C mpi finalize
      call mpi_finalize(ierr)
      if (ierr.eq.0) then
         write(*,*)'mpi_finalize successfull'
      else
         write(*,*)'mpi_finalize failed'
         stop
      endif
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
      dimension pnw(nsolid),pnd(nsolid),etc(3)
      open (30,file='solids.dat')

      read (30,*)
      do 10, i=1,m3
      read (30,'(a10,2(e10.5))')dumchar,pnw(i),pnd(i)
      write(*,'(a10,2(e10.5))')dumchar,pnw(i),pnd(i)
  10  continue
      read (30,'(a10,3(e10.5))')dumchar,etc(1),etc(2),etc(3)
      write(*,'(a10,3(e10.5))')dumchar,etc(1),etc(2),etc(3)
      read (30,'(a10,2(e10.5))')dumchar,xnaohmw,xnaohd
      write(*,'(a10,2(e10.5))')dumchar,xnaohmw,xnaohd
      read (30,'(a10,2(e10.5))')dumchar,xkohmw,xkohd
      write(*,'(a10,2(e10.5))')dumchar,xkohmw,xkohd
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
      dimension pn(nsolid,nnodex),pnw(nsolid),pnd(nsolid),etc(3)
      dimension por(nnodex+2),cn(ncompl,nnodex)

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

        real*8 hb(NNODEx+2)
        integer nxmax,nymax,ihb(NNODEx+2)
        character*10 text,  cname

        open(31, file=cname)
        read (31, *)nxmaxx,faktor
        write(*,*)'hol nx  faktor', nxmax,faktor, cname
         read(31, 1010)(ihb(i),i=1,nxmax)
   10   continue
 1010   format (83i4)
        do 20, i=1,nxmax
        hb(i)= faktor*ihb(i)
 20     continue
        close (31)
        return
        end
