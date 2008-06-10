c***********************************************************************
c***********************************************************************
c
c04      subroutine datin
      subroutine datin(itype,idynam,in1,in2,nxmax,kmax,itemp,lnh,li
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

      implicit double precision (a-h,o-z)
c
c  *********************************************************************
c  this subroutine reads the input data and prints it out
c  *********************************************************************
c

      include 'gwheader.inc'
c      include 'kinetics.inc'    ! kinet

c-coeff      integer s,ss,itest,ncyc,nxmax,ny,isteu,inma,ipfile,ntim,npin
      integer s,itest,ncyc,nxmax,ny,isteu,inma,ipfile,ntim,npin
c04      real*4 t1,t2
      character *79 title
      character *10 dumb(nbasis),dumc(ncompl),dump(nsolid),input
      character *9 datei1
      character*10 datei2
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
ckinet>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c <<<
c04      common /betat/ betac
c04      common /dismdl/ idismdl,indxca,indxsi,indxcs,csbndry(0:2),
c04     1                conhcsm(0:3,2,4)
      dimension csbndry(0:2),conhcsm(0:3,2,4)


c >>>
c04      common /in/ itype,idynam,in1,in2,nxmax,kmax,itemp,lnh,li,lb,
c04     *dtini,dtmax,doa,dob,vo,tprint(25)
c04     2,num(6),vjb(nbasis),indexi(nbasis),bi(nbasis)
c04     3,gesbi(nbasis),indexb(nbasis),bc(nbasis),gesbb(nbasis)
c04     3,con(ncompl+nsolid,4),s(nbasis,ncompl),vjc(ncompl),
c04     4ss(nbasis,nsolid),gespi(nsolid),gespb(nsolid)
c04     5,gesci(ncompl),gescb(ncompl),err,dtmult,dtdiv
c04     6,iterm,iterj,iterd,ndiv,lne,lnhc,xmax,gc,rw
c04     7,temp,tempi,tempo,hcw,hcm,tcm
c04      common /inc/ dumb(nbasis),dumc(ncompl),dump(nsolid)
c04      common /maxvals/ maxnmax,maxnb,maxnc,maxnp
c04      common /tempdep/ itmpdep(ncompl+nsolid)
c04      common /titl/ title

      dimension tprint(25),num(6),vjb(nbasis),indexi(nbasis)
     1,bi(nbasis),gesbi(nbasis),indexb(nbasis),bc(nbasis)
     2,gesbb(nbasis),con(ncompl+nsolid,4),s(nbasis,ncompl)
     3,vjc(ncompl),ss(nbasis,nsolid),gespi(nsolid),gespb(nsolid)
     4,gesci(ncompl),gescb(ncompl),itmpdep(ncompl+nsolid)
     5,index_i1(nbasis),bi_i1(nbasis),gesb_i1(nbasis)
     6,index_i2(nbasis),bi_i2(nbasis),gesb_i2(nbasis)
     7,gesp_i1(nsolid),gesp_i2(nsolid)
     8,gesc_i1(ncompl),gesc_i2(ncompl)
     5,gespvf_i1(nsolid),gespvf_i2(nsolid)
     9,i_bcp_gemx(nbasis+ncompl+nsolid)

ckg44 idatin is not inialized! 
c      idatin=idatin+1
c      write(*,*)'subroutine datin aufgerufen',idatin

      open(7,file='initial.out')
c

      OPEN(4, FILE='TEST.INI',STATUS='UNKNOWN')


c
c  **************************
c  read physical parameters
c  **************************
      read (4,'(A50)') title
      write (*,*) title
      write (7,*) title
      write (6,1011)
      write (7,1011)
      read (4,910) itype,idynam,idismdl
      write(6,*) 'itype,idynam,idismdl',itype,idynam,idismdl
      write(7,*) 'itype,idynam,idismdl',itype,idynam,idismdl
      read (4,910) in1,in2
      write(6,*) 'in1,in2',in1,in2
      write(7,*) 'in1,in2',in1,in2
      read (4,915) err,dtmult,dtdiv,iterm,iterj,iterd,ndiv
      write(6,*)'err,dtmult,dtdiv,iterm,iterj,iterd,ndiv',      
     +err,dtmult,dtdiv,iterm,iterj,iterd,ndiv
      read (4,920) nxmax,xmax,gc,rw
      write (6,*)'nxmax,xmax,gc,rw',nxmax,xmax,gc,rw
      write (7,*)'nxmax,xmax,gc,rw',nxmax,xmax,gc,rw
      read (4,920) kmax,dtini,dtmax
      read (4,900) (tprint(k),k=1,kmax)
      write(*,*)kmax,(tprint(k),k=1,kmax)
      write(7,*)kmax,(tprint(k),k=1,kmax)
      read (4,900) doa,dob,vo

	read (4,910) 
c                 ^(ibnd(k),k=1,2)
      read (4,910) itemp
      if (itemp.eq.0) read (4,900) temp
      if (itemp.eq.1) read (4,1002) temp,tempi,hcw,hcm,tcm
      if (itemp.eq.2) read (4,901) tempi,tempo
      read (4,910) (num(i),i=1,3),nmineq,lnh,lne,lnhc     !kinet
      write(7,910) (num(i),i=1,3),nmineq,lnh,lne,lnhc     !kinet
c04      read (4,910) (num(i),i=1,3),lnh,lne,lnhc
c                                         ^added 10/8/93      
c <<<
c      if (num(1).gt.maxnb) call error (2, num(1))
c      if (num(2).gt.maxnc) call error (3, num(2))
c     -----------------------------------------------------------------
c     num(5) is the total number of solids in the dissolution model.
c     one of the built-in assumptions of the dissolution model is that
c     there are a total of three solids in the model, only two of which
c     are present in each region of C/S.
c     -----------------------------------------------------------------
      if (idismdl.gt.0) then
         num(5) = 3
         nmineq = 3         !idismdl + kinetic possible 2004
      else
         num(5) = 0
      end if
c     -----------------------------------------------------------------
c     num(3) is input as the number of solids not in the dissolution
c     model, but it's value is changed below to include the number of
c     solids in the dissolution model.
c     -----------------------------------------------------------------
      num(3) = num(3) + num(5)
c      if (num(3).gt.maxnp) call error (4, num(3))
c >>>
      m1=num(1)
      m2=num(2)
      m3=num(3) 
c <<<
      m5=num(5)
c >>>
c  ************************************
c  write out some of the constants used
c  ************************************
c      if (ibnd(1).eq.1) write (6,2061)
c     if (ibnd(1).eq.0) write (6,2062)
c      write (6,2071)
c      if (itemp.eq.0) write (6,2065) temp
      if (itemp.eq.0) write (7,2065) temp
c      if (itemp.eq.1) write (6,2067) temp,hcw,hcm,tcm
c      if ((itemp.eq.1).and.(ibnd(2).eq.0)) write (6,2068) tempi
c      if ((itemp.eq.1).and.(ibnd(2).eq.1)) write (6,2069) tempi
c      if (itemp.eq.2) write (6,2070) tempi,tempo
c      if (itemp.eq.1) write (6,2073)
c      if (itype.eq.0) write (6,600) doa,dob,vo,nxmax
c      if (itype.eq.1) write (6,601) doa,dob,vo,nxmax
c      write (6,610) kmax,(tprint(i),i=1,kmax)
c      write (6,620) err,dtmult,dtdiv,iterm,iterj,iterd,ndiv
      write (7,610) kmax,(tprint(i),i=1,kmax)
      write (7,620) err,dtmult,dtdiv,iterm,iterj,iterd,ndiv
      li=0
      lb=0
c
c  *****************************************************
c  read in information on basis species
c  *****************************************************
      nspecm=m1

      do 17 j=1,nspecm
      read (4,950) dumb(j),vjb(j),indexi(j),bi(j),gesbi(j),indexb(j),
     1bc(j),gesbb(j),index_i1(j),bi_i1(j),gesb_i1(j)
     2,index_i2(j),bi_i2(j),gesb_i2(j),i_bcp_gemx(j)
      write(*,950) dumb(j),vjb(j),indexi(j),bi(j),gesbi(j),indexb(j),
     1bc(j),gesbb(j),index_i1(j),bi_i1(j),gesb_i1(j)
     2,index_i2(j),bi_i2(j),gesb_i2(j),i_bcp_gemx(j)
      if (indexi(j).eq.4) li=j
      if (indexb(j).eq.4) lb=j
      if (index_i1(j).eq.4) li_i1=j
      if (index_i2(j).eq.4) li_i2=j
  17  continue
c
c  ***********************************************
c  write out the conditions on the basis species
c  ***********************************************
      write (6,800) nspecm
      write (7,800) nspecm
      do 20 j=1,nspecm
      write (6,805) dumb(j),vjb(j),indexi(j),bi(j),indexb(j),bc(j)
     1,index_i1(j),bi_i1(j),index_i2(j),bi_i2(j)
      write (7,805) dumb(j),vjb(j),indexi(j),bi(j),indexb(j),bc(j)
     1,index_i1(j),bi_i1(j),index_i2(j),bi_i2(j)
   20 continue
      ncmplxp=m2
c
c  ********************************************************
c  if dissociation of water is included, read and write out
c  information about the oh ion
c  ---- deleted
c  ********************************************************

   21 if (m2.eq.0) go to 35
c
c  *****************************************************
c  read in information on complexes 
c  *****************************************************
      do 27 i=1,m2
      kk=i+m3
      read (4,2005) dumc(i),itmpdep(kk),(con(kk,j),j=1,4)
     *             ,i_bcp_gemx(m1+i)
c >>>                       ^ added 9/87
      write (*,2005) dumc(i),itmpdep(kk),(con(kk,j),j=1,4)
     *             ,i_bcp_gemx(m1+i)
c >>>                       ^ added 9/87
 2005 format (a10,i5,4e12.5,i10)

      read (4,1004) (s(j,i),j=1,m1)
      vj=0.
      do 22 j=1,m1
   22 vj=vj+s(j,i)*vjb(j)
      vjc(i)=vj
c      write(*,*)'resultant charge for complex i',i,vjc(i)
      read (4,960) gesci(i),gescb(i)
   27 continue
c
c  *****************************************
c  write out information about the complexes
c  *****************************************
c      write (6,802) ncmplxp
      write (7,802) ncmplxp

      do 30 i=1,m2
      write (6,820) dumc(i),vjc(i),(s(j,i),j=1,m1)
      write (7,820) dumc(i),vjc(i),(s(j,i),j=1,m1)
   30 continue
      do 31 i=1,m2
ckg44 gesc_il and gesc_i2 are not initialized (never read?)
c   31 write(*,960) gesci(i),gescb(i),gesc_i1(i),gesc_i2(i)
   31 write(*,960) gesci(i),gescb(i)
c
   35 continue
c
c  *****************************
c  read in information on solids
c  *****************************
      phi=1.                              !phi  = porosity     kinet
c >>> go to below changed from 41 to 60
      if (m3.eq.0) go to 60
      write (6,840) m3
      write (7,840) m3
      do 40 i=1,m3
      read (4,3005) dump(i),itmpdep(i),(con(i,j),j=1,4)
     *             ,i_bcp_gemx(m1+m2+i)
c >>>                       ^ added 9/87
 3005 format (a10,i5,4e12.5,i50)
      write(*,3005)dump(i),itmpdep(i),(con(i,j),j=1,4)
     *,i_bcp_gemx(m1+m2+i)
      write(7,3005)dump(i),itmpdep(i),(con(i,j),j=1,4)
     *,i_bcp_gemx(m1+m2+i)
c04      read (4,1004) (ss(j,i),j=1,m1)
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>    kinet
      write(7,1005)dump(i),itmpdep(i),(con(i,j),j=1,4)
c-coeff      read (4,1004) (ss(j,i),j=1,m1)
      read (4,*) (ss(j,i),j=1,m1)
c-coeff       write(*,1004) (ss(j,i),j=1,m1)
c-coeff       write(7,1004) (ss(j,i),j=1,m1)
      write(*,*) (ss(j,i),j=1,m1)
cc2004      pause
c
c  **************************************
c  read/write out information about the solids
c  **************************************
c  Init.&bound. cond. (vol. fraction), #grains/m3rock, vol. mol (cm3/mol)
      read (4,*) gespvfi(i),gespvfb(i),gespvf_i1(i),gespvf_i2(i)
     1   ,rng(i),volmol(i) 
      phi=phi-gespvfi(i)   !Porosity    -pep-

c04>>>>>>>>>>>>>>>                          to ensure imdex=3 compatibility with kinetic version
      gespi(i)=gespvfi(i)*1000/volmol(i) 
      gespb(i)=gespvfb(i)*1000/volmol(i) 
      gesp_i1(i)=gespvf_i1(i)*1000/volmol(i) 
      gesp_i2(i)=gespvf_i2(i)*1000/volmol(i) 
      write(*,*)i,gespi(i),gespb(i),gesp_i1(i),gesp_i2(i)
c04<<<<<<<<<<<<<<<<<

      write(*,2510) gespvfi(i),gespvfb(i),gespvf_i1(i),gespvf_i2(i)
     1   ,rng(i),volmol(i),phi
      write(7,2510) gespvfi(i),gespvfb(i),gespvf_i1(i),gespvf_i2(i)
     1   ,rng(i),volmol(i),phi
 2510 format(7(1x,e8.2)) 

c     Mineral kinetics parameters    -pep-
      read(4,*) ikin(i)        !0-equilibrium, 1-kinetics
      write(*,*) 'ikin = ', ikin(i) 
      write(7,*) 'ikin = ', ikin(i) 
      if (ikin(i).eq.1) then
c      read k(mol/m2/sec),flag for f(dG), and exponents for f(dG)
         read(4,*) kmp(i),ifgp(i),ap(i),bp(i) 
         write(*,2511) kmp(i),ifgp(i),ap(i),bp(i)
         write(7,2511) kmp(i),ifgp(i),ap(i),bp(i)
 2511    format(2x,e8.2,i4,2(1x,e8.2))
         read(4,*) (eabp(i,j),j=1,m1)  !exponents for cat./inh. (basis)
         write(*,2512) (eabp(i,j),j=1,m1)  !exponents for cat./inh. (basis)
         write(7,2512) (eabp(i,j),j=1,m1)  !exponents for cat./inh. (basis)
 2512    format(2x,20(f3.1,1x))
         read(4,*) (eacp(i,j),j=1,m2)  !exponents for cat./inh. (complexes)
         write(*,2512) (eacp(i,j),j=1,m2)  !exponents for cat./inh. (complexes)
         write(7,2512) (eacp(i,j),j=1,m2)  !exponents for cat./inh. (complexes)
         read(4,*) kmd(i),ifgd(i),ad(i),bd(i)
         write(*,2511) kmd(i),ifgd(i),ad(i),bd(i)
         write(7,2511) kmd(i),ifgd(i),ad(i),bd(i)
         read(4,*) (eabd(i,j),j=1,m1)
         write(*,2512) (eabd(i,j),j=1,m1)
         write(7,2512) (eabd(i,j),j=1,m1)
         read(4,*) (eacd(i,j),j=1,m2)
         write(*,2512) (eacd(i,j),j=1,m2)
         write(7,2512) (eacd(i,j),j=1,m2)
      else
         read(4,*)
         read(4,*)
         read(4,*)
         read(4,*)
         read(4,*)
         read(4,*)
      endif
c<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<    kinet
c      read (4,960) gespi(i),gespb(i)

c04      read (4,960) gespi(i),gespb(i)
c
c  **************************************
c  write out information about the solids
c  **************************************
c      write (6,971) dump(i),(ss(j,i),j=1,m1)
c-coeff      write (7,971) dump(i),(ss(j,i),j=1,m1)
      write (7,*) dump(i),(ss(j,i),j=1,m1)
   40 continue
c
c  -----------------------------------------------------------------------------
c  read and write information about dissolution model.
c  note that the following has been ASSUMED:
c    the range of C/S values in the solid has been divided into three regions
c      (four if you count C/S = 0 as a region),
c    the number of components in each region is two (except for region 0, where
c      there is only one component),
c    the identity of the components (Ca(OH)2 = C, SiO2 = S, CaH2SiO4 = CS),
c    there are no more than four constants in the variable log^Kso expressions.
c  ----------------------------------------------------------------------------
c <<<
      if (idismdl.eq.1) then
         read (4,1100) csbndry(0),csbndry(1),csbndry(2)
         read (4,1105) indxca,indxsi,indxcs
c ------ check index values read above
         mc = m3 - m5
c         if (indxca.le.mc) call error (5,indxca)
c         if (indxsi.le.mc) call error (6,indxsi)
c         if (indxcs.le.mc) call error (7,indxcs)
         read (4,1110) (conhcsm(0,1,j),j=1,4)
         do 50 i = 1,m5
            read (4,1110) (conhcsm(i,1,j),j=1,4)
            read (4,1110) (conhcsm(i,2,j),j=1,4)
   50    continue

         write (6,1200)
         write (6,1205) csbndry(0),csbndry(1),csbndry(2)
         write (6,1210) indxca,indxsi,indxcs
         write (6,1215)
         write (6,1220) (conhcsm(0,1,j),j=1,4)
         write (7,1200)
         write (7,1205) csbndry(0),csbndry(1),csbndry(2)
         write (7,1210) indxca,indxsi,indxcs
         write (7,1215)
         write (7,1220) (conhcsm(0,1,j),j=1,4)
         do 55 i = 1,m5
            write (6,1221) i,(conhcsm(i,1,j),j=1,4)
            write (6,1222) i,(conhcsm(i,2,j),j=1,4)
            write (7,1221) i,(conhcsm(i,1,j),j=1,4)
            write (7,1222) i,(conhcsm(i,2,j),j=1,4)
   55    continue
      end if

      CLOSE(4)

   60 return
c
c  *****************
c  format statements
c  *****************
  600 format(/,1x,'dispersivity =',1pe10.2,/,1x,'diffusion coefficient =
     1',1pe10.2,/,1x,'fluid velocity =',1pe10.2,//,1x,'number of nodes =
     2',i4)
  601 format(/,1x,'dispersivity =',1pe10.2,/,1x,'diffusion coefficient =
     1',1pe10.2,/,1x,'fluid velocity =',1pe10.2,'/ r',//,1x,'number of n
     2odes =',i4)
  610 format(/,1x,'number of print outs =',i4,4x,'given at the following
     1 times:',/,4x,8(1pe9.2),(/,4x,8(1pe9.2)))
  620 format(/,1x,'err=',1pe10.2,5x,'dtmult=',0pf5.2,5x,'dtdiv=',f5.2,/,
     11x,'iterm=',i4,5x,'iterj=',i4,5x,'iterd=',i4,5x,'ndiv=',i4)
  800 format(/,1x,i4,1x,'basis species',10x,'initial conditions',7x,'bou
     1ndary conditions',/,16x,'charge',7x,'type',7x,'value',9x,'type',7x
     2,'value')
  802 format(/,1x,i4,1x,'complexes',/,16x,'charge',12x,'stoichiometry')
c  805 format(5x,a10,f5.0,6x,i5,4x,1pe12.4,4x,i5,4x,1pe12.4)
  805 format(5x,a10,f5.0,4(6x,i5,4x,1pe12.4))
  820 format(5x,a10,f5.0,2x,20i5)
  840 format(/,1x,i4,1x,'solids',/,21x,'stoichiometry')
  900 format (10x,7e10.3)
  901 format (10x,2f10.3,i5)
  910 format (10x,7i10)
  915 format (10x,3e10.3,4i5)
  920 format (10x,i10,6e10.3)
  950 format (a10,f5.0,i5,2f10.3,i5,2f10.3,i5,2f10.3,i5,2f10.3,i10)        ! 4 columns
  960 format (10x,4f10.3)
  971 format (5x,a10,1x,10i5)
 1000 format (1x,a79)
 1002 format (10x,f10.3,4e10.2)
c <<<
 1003 format (a10,f5.0,i5,4e12.5)
c >>>
 1004 format (10x,20i5)
c >>> 1005 format (a10,4e12.5) - changed 9/87
 1005 format (a10,i5,4e12.5)
 1011 format(/,1x,'unit of mass is kilogram',10x,'unit of length is mete
     1r',/,1x,'unit of time is second',12x,'unit of electrical potential
     2 is volt',/,1x,'unit of temperature is celsius',4x,'unit of heat i
     3s joule',/,1x,'exception: units for standard enthaply of ',
     4'reaction (input by user) are kJ/mol')
c <<<
 1100 format (10x,3f10.4)
 1105 format (10x,3i5)
 1110 format (10x,4e12.5)
 1200 format (1x,/,1x,'dissolution of hydrated calcium-silicates ',
     1'included in this simulation',/,1x,'there are four regions in ',
     2'the model',/,1x,'in region 0, the solid is comprised of one ',
     3'component',/,1x,'in regions I, II and III, the solid is ',
     4'comprised of two components ')
 1205 format (1x,/,1x,'the boundary between regions 0 and I is ',
     1'located at C/S = ',f4.2,/,1x,'the boundary between regions I ',
     2'and II is located at C/S = ',f4.2,/,1x,'the boundary between ',
     2'regions II and III is located at C/S = ',f4.2)
 1210 format (1x,/,1x,'solid ',i2,' has been designated Ca(OH)2',
     1/,1x,'solid ',i2,' has been designated SiO2',
     1/,1x,'solid ',i2,' has been designated CaH2SiO4')
 1215 format (1x,/,1x,'the constants in the expressions used to ',
     1'calculate log^Kso as C/S varies are',/,1x,'given below')
 1220 format (1x,/,6x,'region 0, component 1: ',4(1pe12.4))
 1221 format (1x,/,6x,'region ',i1,', component 1: ',4(1pe12.4))
 1222 format (1x,/,6x,'region ',i1,', component 2: ',4(1pe12.4))
c >>>
 2061 format(/,1x,'inner boundary conditions are constant solute fluxes'
     1)
 2062 format(/,1x,'inner boundary conditions are constant solute concent
     1rations')
 2065 format(1x,'this simulation is isothermal, temperature =',f5.0)
 2067 format(1x,'simulation of mixing with heat transport, with initial 
     1temperature =',f5.0,/,4x,'heat capacities of fluid and solid are',
     21p,g10.2,2x,'and',g10.2,/,4x'thermal conductivity of solid is',g10
     3.2)
 2068 format(1x,'inner boundary condition is constant temperature =',f5.
     10)
 2069 format(1x,'inner boundary condition is constant flux of fluid with
     1 temperature =',f5.0)
 2070 format(1x,'fixed temperature gradient:'/,4x,f10.3,2x,'at inner bou
     1ndary',/,4x,f10.3,2x,'at outer boundary')
 2071 format(1x,'basis species concentrations at outer boundary are held
     1 at initial values')
 2073 format(1x,'temperature is held at initial value at outer boundary'
     1)
      end
c

c***********************************************************************
c***********************************************************************
c
      subroutine datout(time,itemp,nxmax,lnh,num,x,tmp,bn,cn,pn,acb,acc
     1,eqconst,lne,eh,cs,lnhc)
c >>>                     ^ added 9/87

      implicit double precision (a-h,o-z)

      include 'gwheader.inc'

c
c  *********************************************************************
c  this subroutine prints results of the calculations.
c  *********************************************************************
c
c-coeff      integer s,ss,itest,ncyc,nxmax,ny,isteu,inma,ipfile,ntim,npin
      integer s,itest,ncyc,nxmax,ny,isteu,inma,ipfile,ntim,npin
c      integer ir(21,1),npmax,nbox,partib(21),ismooth ,npkt   
c      real*4 t1,t2
      character *79 title
      character *10 dumb,dumc,dump,input
      character *9 datei1
      character*10 datei2
      character *1 ch
 
       dimension x(nnodex+2),tmp(nnodex+2),bn(nbasis,nnodex+2)
     1,cn(ncompl,nnodex+2),pn(nsolid,nnodex+2),num(6)
     2,eqconst(ncompl+nsolid,nnodex+2),acb(nbasis,nnodex+2)
     3,acc(ncompl,nnodex+2),cs(nnodex+2),tmpk(nnodex+2),eh(nnodex+2)
c <<<
      common /betat/ betac
      common /dismdl/ idismdl,indxca,indxsi,indxcs,csbndry(0:2),
     1                conhcsm(0:3,2,4)

      common /inc/ dumb(nbasis),dumc(ncompl),dump(nsolid)
      common /maxvals/ maxnmax,maxnb,maxnc,maxnp
      common /tempdep/ itmpdep(ncompl+nsolid)
      common /titl/ title


      idatout=idatout+1
      write(*,*)'subroutine datout aufgerufen',idatout
c
      m1=num(1)
      m2=num(2)
      m3=num(3)
      m6=m2+m3
c
c      write (6,530) time
c
      npri=m1/5
      i1=1
      i4=5
      if (m1.lt.5) i4=m1
      k2=0
c
c  basis species
c
      ny =1 
  405 k1=k2+1
      k2=k1+i4-1
      if (k2.gt.5.or.lne.eq.0) go to 411
c      write (6,620) (dumb(j),j=k1,k2)
      do 410 n=1,nxmax
c      write (6,628) x(n),eh(n),(bn(j,n),j=k1,k2)
  410 continue
      go to 413
  411 write (6,619) (dumb(j),j=k1,k2)
      do 412 n=1,nxmax
c      write (6,625) x(n),(bn(j,n),j=k1,k2)
  412 continue
  413 i1=i1+1
      if (i1.le.npri) go to 405
      i4=m1-k2
      if (i4.gt.0) go to 405

c new 11/8/93
c  ---------------------------------------
c  pH if dissociation of water in included 
c  ---------------------------------------
      if (lnh.gt.0) then
         if (lnhc.eq.0) then  ! H+ is a basis species, OH- is a complex
c            write (6,633)
  633 format (/,5x,'concentrations'/,7x,'x',8x,'pH ')
            do 416 n = 1,nxmax
               ph = - dlog10(acb(lnh,n)*bn(lnh,n))
c               write (6,625) x(n),ph
  416       continue
         else  ! H+ is listed as a complex, OH- is a basis species
c            write (6,633)
            do 418 n = 1,nxmax
               ph = - dlog10(acc(lnhc,n)*cn(lnhc,n))
               write (6,625) x(n),ph
  418       continue
         end if
      end if
c new 11/8/93      

c      
c
c  aqueous phase complexes
c
      if (m2.eq.0) go to 445
      npri=m2/5
      i1=1
      i4=5
      if (m2.lt.5) i4=m2
      k2=0
  420 k1=k2+1
      k2=k1+i4-1
      write (6,619) (dumc(j),j=k1,k2)
      do 421 n=1,nxmax
      write (6,625) x(n),(cn(j,n),j=k1,k2)
  421 continue
      i1=i1+1
      if (i1.le.npri) go to 420
      i4=m2-k2
      if (i4.gt.0) go to 420
c
c  solids
c
  445 if (m3.eq.0) go to 452
      npri=m3/5
      i1=1
      i4=5
      if (m3.lt.5) i4=m3
      k2=0
  450 k1=k2+1
      k2=k1+i4-1
      write (6,619) (dump(j),j=k1,k2)
      do 451 n=1,nxmax
      write (6,625) x(n),(pn(j,n),j=k1,k2)
  451 continue
      i1=i1+1
      if (i1.le.npri) go to 450
      i4=m3-k2
      if (i4.gt.0) go to 450
c <<<
      if (idismdl.eq.1) then
         write (6,640)
         do 453 n = 1,nxmax
            rlogc = dlog10(eqconst(indxca,n))
            rlogs = dlog10(eqconst(indxsi,n))
            rlogcs = dlog10(eqconst(indxcs,n))
            if (cs(n).lt.0.) then
               write (6,642) x(n),rlogc,rlogs,rlogcs
            else
               write (6,641) x(n),cs(n),rlogc,rlogs,rlogcs
            end if
  453    continue
      end if
c >>>
c
c  temperatures and equilibrium constants
c
  452 if (itemp.eq.0) go to 490
c      if (lnh.gt.0.or.m6.gt.0) write (6,623)
      if (lnh.gt.0) go to 460
c      write (6,624)
      do 455 n=1,nxmax
c      write (6,626) x(n),tmp(n)
  455 continue
      go to 470
  460 write (6,627)
      do 465 n=1,nxmax
      write (6,629) x(n),tmp(n),bn(lnh,n)
  465 continue
  470 if (m2.eq.0) go to 480
      npri=m2/5
      i1=1
      i4=5
      if (m2.lt.5) i4=m2
      k2=m3
  475 k1=k2+1
      k2=k1+i4-1
      k1d=k1-m3
      k2d=k2-m3
      write (6,621) (dumc(j),j=k1d,k2d)
      do 476 n=1,nxmax
      write (6,625) x(n),(eqconst(j,n),j=k1,k2)
  476 continue
      i1=i1+1
      if (i1.le.npri) go to 475
      i4=m2-k2d
      if (i4.gt.0) go to 475
  480 if (m3.eq.0) go to 500
      npri=m3/5
      i1=1
      i4=5
      if (m3.lt.5) i4=m3
      k2=0
  485 k1=k2+1
      k2=k1+i4-1
      write (6,621) (dump(j),j=k1,k2)
      do 486 n=1,nxmax
      write (6,625) x(n),(eqconst(j,n),j=k1,k2)
  486 continue
      i1=i1+1
      if (i1.le.npri) go to 485
      i4=m3-k2
      if (i4.gt.0) go to 485
      go to 500
  490 write (6,630) tmp(1)
      if (lnh.gt.0) write (6,631) bn(lnh,1)
      if (m2.eq.0) go to 495
      do 492 i=1,m2
      j=i+m3
      write (6,632) dumc(i),eqconst(j,1)
  492 continue
c <<<
c     print the equilibrium constants only for the solids not included
c     in the dissolution model.  those for the solids in the dissolution
c     model are printed in loop 453 above.
  495 m3 = m3 - num(5)
c 495 if (m3.eq.0) go to 500 
c >>>
      if (m3.eq.0) go to 500
      do 497 i=1,m3
      write (6,632) dump(i),eqconst(i,1)
  497 continue

  500 return
c
c  *****************
c  format statements
c  *****************
  530 format (//,1x,'time=',1pg12.4,10x,'concentrations of species in mo
     1les/liter solution')
  619 format (/,5x,'concentrations'/,7x,'x',4x,5(1x,a10))
  620 format (/,5x,'concentrations'/,7x,'x',8x,'eh ',5(1x,a10))
  621 format (/,5x,'equilibrium constants',/,7x,'x',4x,5(1x,a10))
  623 format (/,8x,'temperatures and equilibrium constants')
  624 format (/,7x,'x',6x,'temp')
  625 format (1x,1p,6e11.3)
  626 format (1x,1pe11.3,0pf7.2)
  627 format (/,7x,'x',6x,'temp',5x,'bn h2ok')
  628 format (1x,1pe11.3,0pf7.3,1p,5e11.3)
  629 format (1x,1pe11.3,0pf7.2,1pe11.3)
  630 format (/,14x,'temp =',f7.2)
  631 format (14x,'bn h2ok =',1pe11.3)
  632 format (1x,a10,'eqconst ='1pe11.3)
c <<<
  640 format (1x,/,5x,'ratio of calcium to silicon in solid phase',/,
     15x,'calculated for dissolution model',/,5x,'C=Ca(OH)2, S=SiO2, ',
     2'CS=CaH2SiO4',//,7x,'x',9x,'C/S',6x,'log^K_C',4x,'log^K_S',3x,
     3'log^K_CS')
  641 format (1x,5(1pe11.3))
  642 format (1x,(1pe11.3),2x,'undefined',3(1pe11.3))
c >>>
      end


