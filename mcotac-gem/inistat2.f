c***********************************************************************
c***********************************************************************
c
c04      subroutine initial(itemp,li,lb,in1,in2,lnh,err,nxmax,tb,
c04     1temp,vo,x,itype,num,con,eqconst,tmp,bi,bc,indexi,indexb
c04     2,q,acb,acc,bn,pn,cn,vjb,vjc,s,ss,bc2,lne,eh,idismdl,cs,tmpk
c04     3,iche,i_sorb,ialkali)
c >>>                                          ^       ^  ^ added 9/87
      subroutine initial(itemp,li,lb,li_i1,li_i2,in1,in2,lnh,err,nxmax  
     +,tb,temp,vo,x,itype,num,con,eqconst,tmp,bi,bc,indexi,indexb
     +,index_i1,bi_i1,gesb_i1,index_i2,bi_i2,gesb_i2
     +,gesp_i1,gesp_i2,gesc_i1,gesc_i2
     +,gespvf_i1,gespvf_i2
     2,q,acb,acc,bn,pn,cn,vjb,vjc,s,ss,lne,eh,idismdl,cs,tmpk
     3,iche,i_sorb,ialkalic,dumb,dumc,dump,itmpdep)



      implicit double precision (a-h,o-z)
      include 'gwheader.inc'
c      include 'kinetics.inc'    ! kinet
c
c  *********************************************************************
c  this subroutine calculates the equilibrium concentrations of all
c  species at each node for the initial and inner boundary or influx
c  conditions.
c  *********************************************************************
c
      integer s,ss,iche(nnodex),i_sorb,ialkali
      character *10 dumb,dumc,dump

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

        real*8 amin(nsolid,nnodex),km(nsolid)
        real*8 eab(nsolid,nbasis),eac(nsolid,ncompl)
        real*8 kmp(nsolid)
        real*8 eabp(nsolid,nbasis),eacp(nsolid,ncompl)
        real*8 kmd(nsolid)
        real*8 eabd(nsolid,nbasis),eacd(nsolid,ncompl)
        real*8 fg(nsolid),omega(nsolid),catinh(nsolid)
        real*8 ratem(nsolid),a(nsolid),b(nsolid)
        real*8 ap(nsolid),ad(nsolid),bp(nsolid),bd(nsolid)
        real*8 otherj(nbasis),otheri(ncompl)
        real*8 dratemdc(nsolid,nbasis)
        real*8 dqdc(nsolid,nbasis),dxdc(ncompl,nbasis)
        real*8 rng(nsolid),volmol(nsolid),por(nnodex+2)
        real*8 gespvfi(nsolid),gespvfb(nsolid)

        real*8 vout(nsolid,nnodex),rout(nsolid,nnodex)
        real*8 wout(nsolid,nnodex)

        integer ikin(nsolid),ifg(nsolid)
        integer ifgp(nsolid),ifgd(nsolid)
ckinet>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      dimension num(6),con(ncompl+nsolid,4)
     1,eqconst(ncompl+nsolid,nnodex),tmp(nnodex),bi(nbasis)
     1,bc(nbasis),indexi(nbasis),indexb(nbasis)
     5,index_i1(nbasis),bi_i1(nbasis),gesb_i1(nbasis)
     6,index_i2(nbasis),bi_i2(nbasis),gesb_i2(nbasis)
     7,gesp_i1(nsolid),gesp_i2(nsolid)
     8,gesc_i1(ncompl),gesc_i2(ncompl)
     9,index(nbasis),wconst(nbasis)
     1,q(nsolid,nnodex),acb(nbasis,nnodex),acc(ncompl,nnodex)
     3,bn(nbasis,nnodex),pn(nsolid,nnodex),cn(ncompl,nnodex)
     3,vjb(nbasis),vjc(ncompl),s(nbasis,ncompl)
     4,ss(nbasis,nsolid),bc2(nbasis),x(nnodex),eh(nnodex)
      dimension cs(nnodex),tmpk(nnodex)
     1,dumb(nbasis),dumc(ncompl),dump(nsolid)
     2,itmpdep(ncompl+nsolid)
     +,nchem_index(4)
     +,gespvf_i1(nsolid),gespvf_i2(nsolid)

c      include 'kinetics.inc'    ! kinet

c04      common /inc/ dumb(nbasis),dumc(ncompl),dump(nsolid)
c04      common /tempdep/ itmpdep(ncompl+nsolid)

      open(7,file='initial.out')
c      iinitial=iinitial+1
c      write(*,*)'subroutine initial aufgerufen',iinitial, nxmax
c
      do 9 i=1,4
   9  nchem_index(i)=0

c04      iche_index=0
      m1=num(1)
      m2=num(2)
      m3=num(3)
      m5=num(5)
      m6=num(6)

c** two different starting compositiones : boundrary and within the model
c04      n_chem = 1
      n_chem_max=maxval(iche)
      write(*,*) 'n_chem_max',n_chem_max

      do 68 ni=1,nxmax
      if (iche(ni).eq.1.and.nchem_index(1).eq.0) nchem_index(1)=ni       
      if (iche(ni).eq.2.and.nchem_index(2).eq.0) nchem_index(2)=ni       
      if (iche(ni).eq.3.and.nchem_index(3).eq.0) nchem_index(3)=ni        
      if (iche(ni).eq.4.and.nchem_index(4).eq.0) nchem_index(4)=ni       
  68  continue


      do 69 nchem=1,n_chem_max
c      write(*,*)'anf 69',nchem,nchem_index
c     pause
c       nchem_index(1)=1       
c
c  *******************************
c  calculate equilibrium constants
c  *******************************
      if(idismdl.eq.1)call eqconcs(nchem_index(nchem)
     +,num,cs,eqconst,pn)
 
      call eqcont(nchem_index(nchem),con,num,lnh,tmp,eqconst,
     1             itmpdep,tmpk,idismdl)

c******* chemical conditions for initial distribution
c****** 1. boundary
      if (iche(nchem_index(nchem)).eq.1) then
        do 2 i=1,m1
                wconst(i)=bc(i)
                index(i)=indexb(i)
    2   continue
        lchrg=lb
      endif
      
c***** 2. initial 1
      if(iche(nchem_index(nchem)).eq.2) then
        do 3 i=1,m1
                wconst(i)=bi(i)
                index(i)=indexi(i)
    3   continue
        lchrg=li
      endif
c***** 3. initial 2
      if(iche(nchem_index(nchem)).eq.3) then
        do 4 i=1,m1
                wconst(i)=bi_i1(i)
                index(i)=index_i1(i)
    4   continue
        lchrg=li_i1
      endif
c***** 4. initial 3
      if(iche(nchem_index(nchem)).eq.4) then
        do 5 i=1,m1
                wconst(i)=bi_i2(i)
                index(i)=index_i2(i)
    5   continue
        lchrg=li_i2
      endif

c
c  ******************************
c  equilibrate nodes with starting chemistry
c  ******************************
c      write(*,*)'vorstatic',nchem,nchem_index(nchem)
c      pause
c04test      nnchem=nchem_index(nchem)
      call static(nchem_index(nchem),wconst,index,lchrg,in1,in2
     1,num,eqconst,q,acb,acc,bn,pn,cn,vjb,vjc,s,ss,lnh,err
     2,lne,con,i_sorb,ialkali,dumb,dumc,dump)
      if (lne.gt.0) call ehcalc(nchem_index(nchem),nchem_index(nchem)
     +                         ,lne,tmp,bn,eh)
c
c  *******************************************
c  write out the boundary or influx conditions
c  *******************************************
   10 tbound=tb
      if (lne.eq.0) write (6,790) tbound
      if (lne.gt.0) write (6,791) tbound,eh(nchem_index(nchem))
      write (6,780)
      if (lne.eq.0) write (7,790) tbound
      if (lne.gt.0) write (7,791) tbound,eh(nchem_index(nchem))
      write (7,780)
c
c  calculate the total aqueous phase concentrations of the species
c  
      do 180 j=1,m1
      cmtotaq=bn(j,nchem_index(nchem))
      if (m2.eq.0) go to 176
      do 175 i=1,m2
       cmtotaq=cmtotaq+s(j,i)*cn(i,nchem_index(nchem))
  175 continue
c  
c  set boundary conditions
c  
  176  continue
      write (6,700) dumb(j),cmtotaq,bn(j,nchem_index(nchem))
     +,acb(j,nchem_index(nchem))
      write (7,700) dumb(j),cmtotaq,bn(j,nchem_index(nchem))
     +,acb(j,nchem_index(nchem))
  180 continue
      if (m2.eq.0) go to 191
c  
c  write out information about the complexes
c  
      write (6,740)
      do 190 i=1,m2
      write (7,750) dumc(i),vjc(i),cn(i,nchem_index(nchem))
  190 write (6,750) dumc(i),vjc(i),cn(i,nchem_index(nchem))
c
  191 if (m3.eq.0) go to 210
c  
c  write out information about the solids
c  
      write (6,841)
      write (7,841)
      do 205 i=1,m3
      write (6,845) dump(i),pn(i,nchem_index(nchem))
      write (7,845) dump(i),pn(i,nchem_index(nchem))
  205 continue
c
  210 continue
c      pause
c  ******************************
c  set boundary condition to specified nodes
c  ******************************
      do 50 ni=1,nxmax
      if (iche(ni).eq.iche(nchem_index(nchem)))then
c        write(*,*)'iche()=nchem_index()',ni,iche(ni),nchem_index(nchem)
        if (lne.gt.0) eh(ni)=eh(nchem_index(nchem))
	  if(idismdl.eq.1)cs(ni)=cs(nchem_index(nchem))
        do 51 j=1,m1
           bn(j,ni)=bn(j,nchem_index(nchem))
           acb(j,ni)=acb(j,nchem_index(nchem))
   51      continue
        if (m6.eq.0) go to 50
        do 52 j=1,m6
           eqconst(j,ni)=eqconst(j,nchem_index(nchem))
   52   continue
        if (m2.eq.0) go to 54
        do 53 j=1,m2
          cn(j,ni)=cn(j,nchem_index(nchem))
          acc(j,ni)=acc(j,nchem_index(nchem))
   53   continue
   54   if (m3.eq.0) go to 56
        do 55 j=1,m3
           pn(j,ni)=pn(j,nchem_index(nchem))
           select case (iche(ni))
             case (1)
             por(ni)=por(ni)-gespvfb(i)
             case (2)
             por(ni)=por(ni)-gespvfi(i)
             case (3)
             por(ni)=por(ni)-gespvf_i1(i)
             case (4)
             por(ni)=por(ni)-gespvf_i2(i)
           end select

   55   continue
   56   continue
        write(*,*)'node ', ni,'set to chem. conditions index ',nchem
c      else
c first node with iche(ni) not equal 1 is used for n_chemindex : iche(ni)= 2
c      if (iche_index.eq.0)iche_index=1
c      if (iche(ni).eq.2.and.nchem_index(2).eq.0) nchem_index(2)=ni       
c      if (iche(ni).eq.3.and.nchem_index(3).eq.0) nchem_index(3)=ni        
c      if (iche(ni).eq.4.and.nchem_index(4).eq.0) nchem_index(4)=ni       
      endif      
  50  continue
c      n_chem = iche_index
c      n_chem = iche_index
c      write(*,*)'ende 69',nchem,n_chem ,nchem_index
  69  continue
cpause      pause
c  ******************************
c  equilibrate nodes with initial chemistry  with  iche(ni)= 2
c  ******************************
c       n_chem = iche_index
c      if(idismdl.eq.1)call eqconcs(n_chem,num,cs,eqconst,pn)
c      call eqcont (n_chem,con,num,lnh,tmp,eqconst,
c     1             itmpdep,tmpk,idismdl)

c      call static(n_chem,bi,bc,indexi,indexb,li,lb,in1,in2
c     1,num,eqconst,q,acb,acc,bn,pn,cn,vjb,vjc,s,ss,lnh,err
c     2,lne,con,i_sorb,ialkali,dumb,dumc,dump)
c      if (lne.gt.0) call ehcalc(n_chem,n_chem,lne,tmp,bn,eh)
c
c  ************************************************
c  set all nodes to the initial chemical conditions
c  ************************************************
      goto 123

  123  continue
       time=0.
     
c      call datout(time,itemp,nxmax,lnh,num,x,tmp,bn,cn,pn,acb,acc,
c     1 eqconst,lne,eh,cs,lnhc)
c      call datout(time,itemp,nxmax,lnh,num,x,tmp,bn,cn,pn,acb,acc
c     1,eqconst,lne,eh,cs,lnhc,dumb,dumc,dump,itmpdep)

c >>>                ^ added 9/87

 300  continue
      close(7)
c      write(*,*)(pn(1,nn),nn=1,nxmax)
c      write(*,*)(pn(2,nn),nn=1,nxmax)

c      pause
      return
c
c  ****************************************
c  format statements for subroutine initial
c  ****************************************
  700 format(5x,a10,1p,e13.4,e15.4,6x,g12.4)
  740 format(5x,'complexes',5x,'charge',7x,'aqueous conc.')
  750 format(5x,a10,f9.0,6x,1pe13.4)
  760 format(/,22x,'initial conditions',/,24x,'temp =',f7.2)
  761 format(/,22x,'initial conditions',/,24x,'temp =',f7.2,/,26x,'eh ='
     1,f7.3)
  780 format(21x,'total',11x,'basis',9x,'activity',/,5x,'basis',7x,'aque
     1ous conc.',2x,'concentration',4x,'coefficient')
  790 format(/,16x,'boundary or influx conditions',/,24x,'temp =',f7.2)
  791 format(/,16x,'boundary or influx conditions',/,24x,'temp =',f7.2,/
     1,26x,'eh =',f7.3)
  841 format(5x,'solids',5x,'equivalent moles/liter solution')
  845 format(5x,a10,1pe13.4)
      end


c***********************************************************************
c***********************************************************************
c
c04      subroutine static(n1,n2,bi,bc,indexi,indexb,li,lb,in1,in2,num
c04     1,eqconst,q,acb,acc,bn,pn,cn,vjb,vjc,s,ss,lnh,err,lne,con,i_sorb
c04     2,ialkali)

      subroutine static(n1,wconst,index,lchrg,in1,in2,num
     1,eqconst,q,acb,acc,bn,pn,cn,vjb,vjc,s,ss,lnh,err,lne,con,i_sorb
     2,ialkali,dumb,dumc,dump)
 
      implicit double precision (a-h,o-z)

      include 'gwheader.inc'
c      include 'kinetics.inc'    ! kinet


c
c  *********************************************************************
c  this subroutine calculates the equilibrium distribution of species
c  for the initial conditions and inner boundary conditions.
c  unknowns are concentrations of basis species and solids.
c  after solution, concentrations of complexes are calculated from mass
c  action relations.
c  this subroutine initializes the array q, the activity products of
c  basis species in equilibrium with solids.
c  *********************************************************************
c
      character*1 ssw
      integer s,ss,i_sorb,ialkali
      character *10 dumb,dumc,dump
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

        real*8 amin(nsolid,nnodex),km(nsolid)
        real*8 eab(nsolid,nbasis),eac(nsolid,ncompl)
        real*8 kmp(nsolid)
        real*8 eabp(nsolid,nbasis),eacp(nsolid,ncompl)
        real*8 kmd(nsolid)
        real*8 eabd(nsolid,nbasis),eacd(nsolid,ncompl)
        real*8 fg(nsolid),omega(nsolid),catinh(nsolid)
        real*8 ratem(nsolid),a(nsolid),b(nsolid)
        real*8 ap(nsolid),ad(nsolid),bp(nsolid),bd(nsolid)
        real*8 otherj(nbasis),otheri(ncompl)
        real*8 dratemdc(nsolid,nbasis)
        real*8 dqdc(nsolid,nbasis),dxdc(ncompl,nbasis)
        real*8 rng(nsolid),volmol(nsolid),por(nnodex)
        real*8 gespvfi(nsolid),gespvfb(nsolid)

        real*8 vout(nsolid,nnodex),rout(nsolid,nnodex)
        real*8 wout(nsolid,nnodex)

        integer ikin(nsolid),ifg(nsolid)
        integer ifgp(nsolid),ifgd(nsolid)
ckinet>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      dimension wconst(nbasis),index(nbasis),num(6),
     1eqconst(ncompl+nsolid,nnodex),ss(nbasis,nsolid),bi(nbasis)
     2,bc(nbasis),acb(nbasis,nnodex),acc(ncompl,nnodex)
     3,re(nbasis+nsolid),vjb(nbasis),vjc(ncompl),s(nbasis,ncompl)
     4,bn(nbasis,nnodex),pn(nsolid,nnodex),cn(ncompl,nnodex)
     5,z(nbasis+nsolid,nbasis+nsolid),q(nsolid,nnodex)
     6,indexi(nbasis),indexb(nbasis)
     7,con(ncompl+nsolid,4)
     8,dumb(nbasis),dumc(ncompl),dump(nsolid)
c
c04      common /inc/ dumb(nbasis),dumc(ncompl),dump(nsolid)

c      istatic=istatic+1
      write(*,*)'subroutine static aufgerufen',istaticnxmax, i_sorb,
     *ialkali

      m1=num(1)
      m2=num(2)
      m3=num(3)
      m4=num(4)
c
       n=n1

c      if (n.eq.1) then
c        do 2 i=1,m1
c                wconst(i)=bc(i)
c                index(i)=indexb(i)
c                write(*,*)'static n=1 index bc ',i,bc(i),index(i)
c    2   continue
c        lchrg=lb
c        in=in2
c      endif
c      if(n.gt.1) then
c        do 4 i=1,m1
c                wconst(i)=bi(i)
c                index(i)=indexi(i)
c                write(*,*)'indexb n>1',i,index(i),n,bi(i)
c    4   continue
c        lchrg=li
c        in=in1
c      endif
c
c  initialize residue array to zero
  5   do 6 i=1,m4
              re(i)=0.
    6 continue
c
c      do 1122 i1=1,m1
c 1122   write(*,*)'i  index(i)  test',i1,index(i1),n

      iter=0
   10 iter=iter+1
c      write(*,*)'iter =',iter
c
c  if convergence is not obtained within 100 iterations, the
c  calculations are stopped
      if (iter.gt.in1) go to 69
c
c  **************************
c  residues for basis species
c  **************************
c
c  calculate activity coefficients
c
      call actco(n,num,bn,cn,vjb,vjc,acb,acc,lne)
c2004      i_sorb=6
      if(i_sorb.gt.0.and.ialkali.gt.0)then
c  *******************************
c   Modification for alkalihydroxide dissolution
       bnnn= acb(m1-1,n)*
c      acc(i_sorb+1,n)*
     * acb(lnh,n)*bn(lnh,n)
       eqconst(m3+i_sorb+1,n) =(10**(-0.18))
     *   *con(m3+i_sorb+1,1)+con(m3+i_sorb+1,1)/bnnn

       bnnn= acb(m1,n)*
     * acb(lnh,n)*bn(lnh,n)
       eqconst(m3+i_sorb+2,n) =(10**(-0.46))
     *   *con(m3+i_sorb+2,1)+con(m3+i_sorb+2,1)/bnnn

c      write(*,*)'eqNaOHs= ',eqconst(m3+i_sorb+1,n)
c     *,con(m3+i_sorb+1,1),acb(lnh,n),bn(lnh,n)
c      write(*,*)'eqKOHs= ',eqconst(m3+i_sorb+2,n)
c     *,con(m3+i_sorb+2,1),bnnn,iter,n
c	pause
      endif
c
c  set conditions to those specified for each basis species
c
      do 20 j=1,m1
c       write(*,*)' loop 20 j=  zz=', j,zz
        re(j)=wconst(j)-bn(j,n)
        if (m2.eq.0) go to 12
c
c  if index = 1 the concentration of the basis species only is given 
c
        if (index(j).eq.1) go to 12
        do 11 i=1,m2
   11   re(j)=re(j)-s(j,i)*cn(i,n)
   12   if (m3.eq.0) go to 20
c ppppppppppppppppppppppppppppppppppppppppppppppppppppppp
c        goto 20
c ppppppppppppppppppppppppppppppppppppppppppppppppppppppp
c
c  index must = 3 for the concentration of a solid to be included in the
c  specification for the species.
c
        if ((index(j).eq.1).or.(index(j).eq.2)) go to 20
        do 15 i=1,m3
c       write(*,*) 'stati pn ij',i, pn(i,n), j, re(j)
c   15  re(j)=re(j)-ss(j,i)*pn(i,n)
        re(j)=re(j)-ss(j,i)*pn(i,n)
c        write(*,*)'stati pn i',n,i,j, re(j),index(j)
  15   continue
c
c ppppppppppppppppppppppppppppppppppppppppppppppppppppppp
   20 continue
c
c  provision for use of a charge balance to specify one of the 
c  concentrations
c
      if (lchrg.eq.0) go to 30
      re(lchrg)=0.
      do 25 j=1,m1
c      write(*,*)'charge'
        if (j.eq.lne) go to 25
        re(lchrg)=re(lchrg)-bn(j,n)*vjb(j)
   25 continue
      if (m2.eq.0) go to 30
      do 26 i=1,m2
   26   re(lchrg)=re(lchrg)-cn(i,n)*vjc(i)
c
   30 continue
      if (m3.eq.0) go to 41
c ppppppppppppppppppppppppppppppppp
c      goto 41
c ppppppppppppppppppppppppppppppppp
c
c  *******************
c  residues for solids
c  *******************
      do 40 k=1,m3
       i=m1+k
c
c  calculate q, activity products of basis species 
c  in equilibrium with solids
c
       prod=1.
       do 36 l=1,m1
            if (ss(l,k).eq.0)go to 36 
c 
c                if (l.eq.lne) GOTO 36
c ^test 080394       
            prod=prod*(acb(l,n)*bn(l,n))**ss(l,k)
   36  continue
       q(k,n)=prod
c
       re(i)=prod-eqconst(k,n)
       if((pn(k,n).le.zero).and.(re(i).le.0.))re(i)=pn(k,n)

c              write(*,*) ' i reeee (i) k pn(k)',i,re(i),k,pn(k,n)
 
   40 continue
c
   41 continue
c
c  *************************************
c  calculate elements of jacobian matrix
c  *************************************
c
c  initialize elements to zero
c
      do 43 i=1,m4
      do 43 j=1,m4
c                   ^^ m4 ersetzt durch 30 ,m4 hat mit dissolution mod. wert 4 !!
      z(i,j)=0.
   43 continue
c
      do 52 i=1,m1
c
c  basis rows, basis columns
c

      do 49 j=1,m1
      zz=0.
      if (i.eq.j) zz=-1.
c     write(*,*)'i,j,zz sum    1' ,i,j,zz,sum
      if (i.eq.lchrg) zz=-vjb(j)
c     write(*,*)'i,j,zz sum    2' ,i,j,zz,sum
      if ((m2.eq.0).or.(index(i).eq.1)) go to 47
      sum=0.
      do 46 k=1,m2
      coef=s(i,k)
      if (i.eq.lchrg) coef=vjc(k)
      if ((s(j,k).eq.0).or.(coef.eq.0.)) go to 46
      mk=m3+k
      prod=eqconst(mk,n)*coef*s(j,k)*acb(j,n)/acc(k,n)
      do 44 l=1,m1
      if ((l.eq.j).or.(s(l,k).eq.0).or.(bn(l,n).eq.0.)) go to 44
c      write(*,*)'prod k l',k,l,acb(l,n),bn(l,n)
      prod=prod*(acb(l,n)*bn(l,n))**s(l,k)
   44 continue
      if ((s(j,k).eq.0).or.(bn(j,n).eq.0.)) go to 45
      ism1=s(j,k)-1
      prod=prod*(acb(j,n)*bn(j,n))**ism1
   45 sum=sum+prod
   46 continue
      zz=zz-sum
   47 z(i,j)=zz
   49 continue
c
      if ((m3.eq.0).or.(index(i).lt.3)) go to 52
c ppppppppppppppppppppppppppppppppppppp
c      goto 52
c ppppppppppppppppppppppppppppppppppppp
c
c  basis rows, solids columns
c
      do 50 k=1,m3
      j=m1+k
      z(i,j)=-ss(i,k)
c       write(*,*)'m3ij,z',i,j,z(i,j) 
   50 continue

c ppppppppppppppppppppppppppppppppppppp
   52 continue
c       write(*,*)'mmmmmmmm33333',m3
      if (m3.eq.0) go to 60
c ppppppppppppppppppppppppppppppppppppp
c      goto 60
c ppppppppppppppppppppppppppppppppppppp
c
      do 59 k=1,m3
      i=m1+k
      if ((q(k,n).le.eqconst(k,n)).and.(pn(k,n).le.
     *zero)) go to 57
c
c  solids rows, basis columns
c
      do 56 j=1,m1
c >>> if (ss(j,k).gt.0) go to 53 <<< correct 9/25/87


cccccc
c      if (ss(j,k).ne.0.or.ss(j,k).gt.0) go to 53

      if (ss(j,k).ne.0) go to 53

      prod=0.
      go to 55
   53 prod=ss(j,k)*acb(j,n)
      do 54 l=1,m1
      if ((l.eq.j).or.(ss(l,k).eq.0)) go to 54
      prod=prod*(acb(l,n)*bn(l,n))**ss(l,k)
   54 continue
c >>> if ((s(j,k).eq.0).or.(bn(j,n).eq.0.)) go to 55 <<< corrected 9/25/87
      if ((ss(j,k).eq.0).or.(bn(j,n).eq.0.)) go to 55
      ism1=ss(j,k)-1
      prod=prod*(acb(j,n)*bn(j,n))**ism1
   55 z(i,j)=prod
   56 continue
      go to 59
c
c  solids rows, solids columns
c
   57 z(i,i)=1.
c
   59 continue
c ppppppppppppppppppppppppppppppppppppp
c
c  *******************************************************
c  calculate corrections and new values of concentrations;
c  check convergence
c  *******************************************************
   60 continue
   
      call simq(m4,re,z,n,isimflag)
c
      if (iter.gt.in1-1) write (6,101) n,iter
      nfail=0
      do 61 j=1,m1
      cc1=bn(j,n)
      bn(j,n)=cc1-re(j)
c090394      write(*,*)'simq bn',bn(6,2),re(6)
c  110193 
      if (bn(j,n).lt.zero) bn(j,n)=zero
c      IF(bn(j,n).gt.1)then 
c       WRITE(*,*)'NONS bn(J,N)>1 at',J,N,iter,bn(j,n)
c       bn(j,n)=bn(j,n)/10.
c      endif 
c  oder doch ?
c      if (iter.gt.in1) write (6,102) j,cc1,bn(j,n),re(j)
cccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (iter.gt.in1-1) then
        write (6,103)j,cc1,bn(j,n),re(j),bc(j)
 103  format(i2,4(1x,e10.4))
      endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (cc1.eq.0.) go to 61
      aaa=dabs(re(j))
      bbb=cc1*err
      if (aaa.gt.bbb)THEN
       nfail=nfail+1
c090394       WRITE(*,*)'FEHLER  !!!! RE = ',RE(J), 'J= ',J,err,m1,m2
c090394     *  ,m3,m4,m5,m6
      ENDIF 
c090394      aaaa=dabs(wconst(6)-bn(6,2))
c090394      if (j.eq.6)write(*,*)'aaaa',j,wconst(6),bn(6,n),aaaa
c090394      if(aaaa.gt.err)then
c090394      nfail=nfail+1 
c090394      write(*,*)'nfail'
c090394      endif
   61 continue
      if (m3.eq.0) go to 63
c ppppppppppppppppppppppppppppppppppppp
      goto 63
c ppppppppppppppppppppppppppppppppppppp
      do 62 j=1,m3
      k=m1+j
      cc1=pn(j,n)
      pn(j,n)=cc1-re(k)

cc11  newinput

c      if (n.eq.2)write(*,*)'cc1 n=2',cc1,pn(j,n),re(k),k

c      if(pn(j,n).gt.1000.)write(*,*)'groesser1000 jnk alt',j,n,k,
c     *cc1,re(k)

      if (pn(j,n).lt.zero.or.eqconst(j,n).eq.rlarge)pn(j,n)=0.
      if (iter.gt.in1-1) write (6,102) k,cc1,pn(j,n),re(k)
      if (cc1.eq.0.) go to 62
      aaa=dabs(re(k))
      bbb=cc1*err
      if (aaa.gt.bbb)THEN
       nfail=nfail+1
c       WRITE(*,*)'FEHLER  !!!! RE = ',RE(K), 'K= ',K
      ENDIF     
      if (iter.gt.in1-1) then
         write (6,104)k,cc1,pn(j,n),re(k)
  104           format(i2,3(1x,e10.4))
      endif



   62 continue

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     User input into iteration
      if (iter.gt.in1-1) then
        write(*,*)'try a new guess for a basic species j  (y/n)?'
        read(*,*)ssw
        if(ssw.eq.'y'.or.ssw.eq.'Y')then
                write(*,*)'give:  j,bn(j,n)'
                read(*,*)j,bn(j,n)
                iter=0
        else
                write(*,*)'try a new guess for a solid k= ?'
                read(*,*)ssw
                if(ssw.eq.'y'.or.ssw.eq.'Y')then
                        write(*,*)'give:  k,pn(k-m1,n)'
                        read(*,*)k,pn(k-m1,n)
                        iter=0
                endif
        endif
      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c ppppppppppppppppppppppppppppppppppppp
   63 if (m2.eq.0) go to 66
      do 65 j=1,m2
      i=m3+j
      prod=1.
      do 64 l=1,m1
      if (s(l,j).eq.0.or.bn(l,n).eq.0) go to 64
c                     .or.bn(l,n).eq.0 
c                     ^new 08 02 95 2004  
      prod=prod*(acb(l,n)*bn(l,n))**s(l,j)
   64 continue
      cn(j,n)=eqconst(i,n)*prod/acc(j,n)
      if (cn(j,n).lt.zero) cn(j,n)=zero
c      write(*,*)'cn',j,cn(j,n)
   65 continue
c
   66 if (nfail.gt.0) go to 10
c
   90 continue

      return
   69 iterm=iter-1
      write (6,100) n,iterm
      stop
  100 format(1x,'no convergence in static at node',i4,' after iter ='i3)
  101 format(4x,'node='i4,5x,'iter=',i5,5x,'values for this iteration',/
     1,/,4x,'species no.',10x,'old value',10x,'new value',10x,'error',/)
  102 format(/,4x,i10,10x,1pe12.4,9x,e12.4,9x,e12.4)

      end

c***********************************************************************
c***********************************************************************
c
      subroutine init2(itemp,li,lb,in1,in2,lnh,err,nspezx,tb,
     1temp,vo,x,itype,num,con,eqconst,tmp,bi,bc,indexi,
     2indexb,q,acb,acc,bn,pn,cn,vjb,vjc,s,ss,bc2,lne,eh,idismdl,
     3cs,tmpk,bo,co,po,i_sorb,ialkali,dumb,dumc,dump,itmpdep)
c >>>                                          ^       ^  ^ added 9/87


      implicit double precision (a-h,o-z)
      include 'gwheader.inc'

c
c  *********************************************************************
c  this subroutine calculates the equilibrium concentrations of all species 
c  *********************************************************************
c
      integer s,ss,iche(nnodex),i_sorb,ialkali
      character *10 dumb,dumc,dump
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

        real*8 amin(nsolid,nnodex),km(nsolid)
        real*8 eab(nsolid,nbasis),eac(nsolid,ncompl)
        real*8 kmp(nsolid)
        real*8 eabp(nsolid,nbasis),eacp(nsolid,ncompl)
        real*8 kmd(nsolid)
        real*8 eabd(nsolid,nbasis),eacd(nsolid,ncompl)
        real*8 fg(nsolid),omega(nsolid),catinh(nsolid)
        real*8 ratem(nsolid),a(nsolid),b(nsolid)
        real*8 ap(nsolid),ad(nsolid),bp(nsolid),bd(nsolid)
        real*8 otherj(nbasis),otheri(ncompl)
        real*8 dratemdc(nsolid,nbasis)
        real*8 dqdc(nsolid,nbasis),dxdc(ncompl,nbasis)
        real*8 rng(nsolid),volmol(nsolid),por(nnodex)
        real*8 gespvfi(nsolid),gespvfb(nsolid)

        real*8 vout(nsolid,nnodex),rout(nsolid,nnodex)
        real*8 wout(nsolid,nnodex)

        integer ikin(nsolid),ifg(nsolid)
        integer ifgp(nsolid),ifgd(nsolid)
ckinet>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      dimension num(6),con(ncompl+nsolid,4)
     1,eqconst(ncompl+nsolid,nnodex),tmp(nnodex),bi(nbasis)
     1,bc(nbasis),indexi(nbasis),indexb(nbasis)
     1,q(nsolid,nnodex),acb(nbasis,nnodex),acc(ncompl,nnodex)
     3,bn(nbasis,nnodex),pn(nsolid,nnodex),cn(ncompl,nnodex)
     3,vjb(nbasis),vjc(ncompl),s(nbasis,ncompl)
     4,ss(nbasis,nsolid),bc2(nbasis),x(nnodex),eh(nnodex)
     5,cs(nnodex),tmpk(nnodex)
     6,bo(nbasis,nnodex),po(nsolid,nnodex),co(ncompl,nnodex)
     7,dumb(nbasis),dumc(ncompl),dump(nsolid)
     8,itmpdep(ncompl+nsolid)

c04      common /inc/ dumb(nbasis),dumc(ncompl),dump(nsolid)
c      common /tempdep/ itmpdep(ncompl+nsolid)
c
      common /vels/vx,vxx     !kinet
      real*8 vx(nnodex)       !kinet

c      iinit2=iinit2+1
c      write(*,*)'subroutine initi2 aufgerufen',iinit2
c

      m1=num(1)
      m2=num(2)
      m3=num(3)
      m5=num(5)
      m6=num(6)
      if ((m6.eq.0).and.(lnh.eq.0)) go to 5
c
c  *******************************
c  calculate equilibrium constants
c  *******************************
c <<<
c04      if(idismdl.eq.1)call eqconcs (nspezx,num,cs,eqconst,pn)

c      call eqcont (nspezx,con,num,lnh,tmp,eqconst
c     1             ,itmpdep,tmpk,idismdl)
c  for last two rows at t=25, no further coefficients

c  ******************************
c  equilibrate nodes 1 through n2
c  ******************************

C       WRITE(*,*)'INIT2',BN(1,1), NSPEZ
    5 call stati2(nspezx,bi,bc,indexi,indexb,li,lb,in1,in2,
     1num,eqconst,q,acb,acc,bn,pn,cn,vjb,vjc,s,ss,lnh,err,lne,
     2         istat2flag,isimflag,bo,po,con,i_sorb,idismdl,cs
     3,ialkali,dumb,dumc,dump)

c 210793
      if (lne.gt.0)call ehcalc(nspezx,nspezx,lne,tmp,bn,eh)

c  110193
      if(istat2flag.eq.1.or.isimflag.eq.1) then
c  following only if no convergence up to id = 4444



      write(*,*)'istat2flag =',istat2flag,'isimflag =',isimflag

c  *******************************************
c  write out the boundary or influx conditions
c  *******************************************
   10 tbound=tb
      if (itemp.eq.2) tbound=tmp(nspezx)
      if (lne.eq.0) write (6,790) tbound
      if (lne.gt.0) write (6,791) tbound,eh(nspezx)
      write (6,780)
c
c  calculate the old aqueous phase concentrations of the species
c  
      do 180 j=1,m1
      cmtotaq=bo(j,nspezx)
      if (m2.eq.0) go to 176
      do 175 i=1,m2
  175 cmtotaq=cmtotaq+s(j,i)*co(i,nspezx)
c  
c  set boundary conditions
c  
  176 continue
      write (6,700) dumb(j),cmtotaq,bo(j,nspezx)
     *            ,acb(j,nspezx)
  180 continue
      if (m2.eq.0) go to 191
c  
c  write out information about the complexes
c  
      write (6,740)
      do 190 i=1,m2
  190 write (6,750) dumc(i),vjc(i),co(i,nspezx)
c
  191 if (m3.eq.0) go to 210
c  
c  write out information about the solids
c  
      write (6,841)
      do 205 i=1,m3
      write (6,845) dump(i),po(i,nspezx)
  205 continue
c
  210 continue
c      if (itemp.eq.2) go to 220

      if (lne.eq.0) write (6,760) tmp(nspezx)
      if (lne.gt.0) write (6,761) tmp(nspezx),eh(nspezx)
      write (6,780)
c
c  calculate the total amount of the species in the aqueous phase
c  
      do 80 j=1,m1
      cmtotaq=bn(j,nspezx)
      if (m2.eq.0) go to 70
      do 60 i=1,m2
   60 cmtotaq=cmtotaq+s(j,i)*cn(i,nspezx)
c  
c  write out this information
c  
   70 write (6,700) dumb(j),cmtotaq,bn(j,nspezx),acb(j,nspezx)
   80 continue
      if (m2.eq.0) go to 91
c  
c  write out information about the complexes
c  
      write (6,740)
      do 90 i=1,m2
      write (6,750) dumc(i),vjc(i),cn(i,nspezx)
   90 continue
   91 if (m3.eq.0) go to 300
c
c  write out information about the solids
c
      write (6,841)
      do 150 i=1,m3
      write (6,845) dump(i),pn(i,nspezx)
  150 continue
c
c  ****************************************
c  format statements for subroutine initial
c  ****************************************
  700 format(5x,a10,1p,e13.4,e15.4,6x,g12.4)
  740 format(5x,'complexes',5x,'charge',7x,'aqueous conc.')
  750 format(5x,a10,f9.0,6x,1pe13.4)
  760 format(/,22x,'no conv. conditions',/,24x,'temp =',f7.2)
  761 format(/,22x,'no conv. conditions',/,24x,'temp =',f7.2,/,26x,'eh =
     1',f7.3)
  780 format(21x,'total',11x,'basis',9x,'activity',/,5x,'basis',7x,'aque
     1ous conc.',2x,'concentration',4x,'coefficient')
  790 format(/,16x,'     old   values  ',/,24x,'temp =',f7.2)
  791 format(/,16x,'  old conditions',/,24x,'temp =',f7.2,/
     1,26x,'eh =',f7.3)
  841 format(5x,'solids',5x,'equivalent moles/liter solution')
  845 format(5x,a10,1pe13.4)

c 220  continue


c      call datout(time,itemp,nxmax,lnh,num,x,tmp,bn,cn,pn,acb,acc,
c     1eqconst,lne,eh,cs,lnhc)
      stop
      endif
c  before, only if no convergence up to id = 4444
c  id = 4444

ccccccccccccccccccccc
c
c      if (lne.gt.0) call ehcalc(nspez,lne,tmp,bn,eh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  300 return
      end


c***********************************************************************
c***********************************************************************
c
      subroutine stati2(nspezx,bi,bc,indexi,indexb,li,lb,in1,
     1in2,num,eqconst,q,acb,acc,bn,pn,cn,vjb,vjc,s,ss,lnh,err,
     2lne,istat2flag,isimflag,bo,po,con,i_sorb,idismdl,cs
     3,ialkali,dumb,dumc,dump)

      implicit double precision (a-h,o-z)
      include 'gwheader.inc'
c      include 'kinetics.inc'    ! kinet
c
c  *********************************************************************
c  this subroutine calculates the equilibrium distribution of species
c  for the initial conditions and inner boundary conditions.
c  unknowns are concentrations of basis species and solids.
c  after solution, concentrations of complexes are calculated from mass
c  action relations.
c  this subroutine initializes the array q, the activity products of
c  basis species in equilibrium with solids.
c  *********************************************************************
c
      character*1 ssw
      character*10 dumb,dumc,dump
      integer s,ss,ialkali
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

        real*8 amin(nsolid,nnodex),km(nsolid)
        real*8 eab(nsolid,nbasis),eac(nsolid,ncompl)
        real*8 kmp(nsolid)
        real*8 eabp(nsolid,nbasis),eacp(nsolid,ncompl)
        real*8 kmd(nsolid)
        real*8 eabd(nsolid,nbasis),eacd(nsolid,ncompl)
        real*8 fg(nsolid),omega(nsolid),catinh(nsolid)
        real*8 ratem(nsolid),a(nsolid),b(nsolid)
        real*8 ap(nsolid),ad(nsolid),bp(nsolid),bd(nsolid)
        real*8 otherj(nbasis),otheri(ncompl)
        real*8 dratemdc(nsolid,nbasis)
        real*8 dqdc(nsolid,nbasis),dxdc(ncompl,nbasis)
        real*8 rng(nsolid),volmol(nsolid),por(nnodex)
        real*8 gespvfi(nsolid),gespvfb(nsolid)

        real*8 vout(nsolid,nnodex),rout(nsolid,nnodex)
        real*8 wout(nsolid,nnodex)

        integer ikin(nsolid),ifg(nsolid)
        integer ifgp(nsolid),ifgd(nsolid)
ckinet>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      dimension wconst(nbasis),index(nbasis),num(6)
     1,eqconst(ncompl+nsolid,nnodex),ss(nbasis,nsolid),bi(nbasis)
     2,bc(nbasis),acb(nbasis,nnodex),acc(ncompl,nnodex)
     3,re(nbasis+nsolid),vjb(nbasis),vjc(ncompl),s(nbasis,ncompl)
     4,bn(nbasis,nnodex),pn(nsolid,nnodex),cn(ncompl,nnodex)
     5,z(nbasis+nsolid,nbasis+nsolid),q(nsolid,nnodex)
     6,indexi(nbasis),indexb(nbasis),bo(nbasis,nnodex)
     7,po(nsolid,nnodex)
     8,con(ncompl+nsolid,4),cs(nnodex)
     9,dumb(nbasis),dumc(ncompl),dump(nsolid)
c      common /inc/ dumb(nbasis),dumc(ncompl),dump(nsolid)
      common /vels/vx,vxx     !kinet
      real*8 vx(nnodex)       !kinet


c      istati2=istati2+1
c      write(*,*)'subroutine stati2 aufgerufen',istati2

c
      nfail2=0
      m1=num(1)
      m2=num(2)
      m3=num(3)
      m4=num(4)
      iterneu=0
        do 4 i=1,m1
                wconst(i)=bi(i)
    4   continue
        lchrg=li
c        in=in1
  
c
c  initialize residue array to zero
  5   do 6 i=1,m4
              re(i)=0.
    6 continue

      iter=0
   10 iter=iter+1
c
c  if convergence is not obtained within "in2" iterations, the
c  calculations are stopped
      if (iter.gt.in2) go to 69
c
c  **************************
c  residues for basis species
c  **************************
c
c  calculate activity coefficients
c


c      if(iter.eq.1.or.iter.eq.5.or.iter.eq.10)then
      call actco(nspezx,num,bn,cn,vjb,vjc,acb,acc,lne)
c      Write(*,*)'actco aufgerufen iter nx ny',iter,nspezx
c      endif
      if(i_sorb.gt.0.and.ialkali.gt.0)then
c  *******************************
c   Modification for alkalihydroxide dissolution
       bnnn= acb(m1-1,nspezx)*
c      acc(i_sorb+1,nspezx)*
     * acb(lnh,nspezx)*bn(lnh,nspezx)
       eqconst(m3+i_sorb+1,nspezx) =(10**(-0.18))
     *   *con(m3+i_sorb+1,1)+con(m3+i_sorb+1,1)/bnnn

       bnnn= acb(m1,nspezx)*
     * acb(lnh,nspezx)*bn(lnh,nspezx)
       eqconst(m3+i_sorb+2,nspezx) =(10**(-0.46))
     *   *con(m3+i_sorb+2,1)+con(m3+i_sorb+2,1)/bnnn

c      write(*,*)'eqNaOHs= ',eqconst(m3+i_sorb+1,nspezx)
c     *,con(m3+i_sorb+1,1),acb(lnh,nspezx),bn(lnh,nspezx)
c      write(*,*)'eqKOHs= ',eqconst(m3+i_sorb+2,nspezx)
c     *,con(m3+i_sorb+2,1),bnnn,iter,nspezx
      endif
      if(idismdl.eq.1.and.(iter.eq.1.or.iter.eq.5.or.iter.eq.10))
     *call eqconcs (nspezx,num,cs,eqconst,pn)
c
c  set conditions to those specified for each basis species
c
      do 20 j=1,m1
c       write(*,*)' loop 20 j=  zz=', j,zz
        re(j)=wconst(j)-bn(j,nspezx)
          if (nspezx.eq.2.and.dabs(re(j)).gt.1.e+32)then
       write(*,'(a4,1x,i2,1x,3(e10.4,1x),i2)')'re_b',
     * j,re(j),bn(j,nspezx),wconst(j)
           pause
          endif
 
        if (m2.eq.0) go to 12
c
c  if index = 1 the concentration of the basis species only is given 
c
c        if (index(j).eq.1) go to 12
        do 11 i=1,m2

      re(j)=re(j)-s(j,i)*cn(i,nspezx)
   11  continue
   12   if (m3.eq.0) go to 20
c
c  index must = 3 for the concentration of a solid to be included in the
c  specification for the species.
c
        do 15 i=1,m3
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  kinet
c this from below
        prod=1.
        do 36 l=1,m1
         if (ss(l,i).eq.0.) go to 36                      ! .or.bn(l,nspezx).eq.0.
         prod=prod*(acb(l,nspezx)*bn(l,nspezx))**
     *        ss(l,i)
c          if (nspezx.eq.2.and.dabs(prod).gt.1.e+50)then
c       write(*,'(a4,1x,i2,1x,3(e10.4,1x),i2)')'prod',
c     * l,prod,bn(l,nspezx),acb(l,nspezx),ss(l,i)
c           pause
c          endif
   36   continue
        q(i,nspezx)=prod


c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  kinet
           if (ikin(i).eq.0) then    
  
              re(j)=re(j)-ss(j,i)*pn(i,nspezx)

           else

              omega(i)=q(i,nspezx)/eqconst(i,nspezx)
              if (omega(i).ne.0.) wout(i,nspezx)=dlog10(omega(i))        ! 2004
              if (omega(i).le.1.0 .and. po(i,nspezx).le.0.0) then
                 po(i,nspezx)=0.0
                 ratem(i)=0.0
                 rout(i,nspezx)=0.0
c                 re(j)=re(j)-ss(j,i)*po(i,nspezx)
                 goto 15
              endif
              if (omega(i).gt.1.0 .and. po(i,nspezx).le.0.0) then
                 po(i,nspezx)=0.0
                 radius=1.0d-8  !r=0.01 micrometers
                 amin(i,nspezx)=4.0*pi*radius*radius*rng(i)
              endif

              if (omega(i).ge.1.0) then
                 km(i)=kmp(i)
                 ifg(i)=ifgp(i)
                 a(i)=ap(i)
                 b(i)=bp(i)
                 do 16 nps=1,m1
                    eab(i,nps)=eabp(i,nps)
 16              continue   
                 do 17 nss=1,m2
                    eac(i,nss)=eacp(i,nss)
 17              continue   
              else
                 km(i)=kmd(i)
                 ifg(i)=ifgd(i)
                 a(i)=ad(i)
                 b(i)=bd(i)
                 do 18 nps=1,m1
                    eab(i,nps)=eabd(i,nps)
 18              continue   
                 do 19 nss=1,m2
                    eac(i,nss)=eacd(i,nss)
 19              continue   
              endif
              
           if (ifg(i).eq.1) then        !f(dG)=|omega**a-1|**b
              fg(i)=dabs(omega(i)**a(i)-1.0)**b(i)
           else if (ifg(i).eq.2) then   !f(dG)=|ln omega|**a
              fg(i)=dabs(log(omega(i)))**a(i)
           endif

           catinh(i)=1.0
           do 420 l=1,m1
              if (eab(i,l).eq.0.0) then
                 catinht=1.0
              else
                 catinht=(bn(l,nspezx)*acb(l,nspezx))**eab(i,l)
              endif
              catinh(i)=catinh(i)*catinht
 420       continue
           do 430 l=1,m2
              if (eac(i,l).eq.0.0) then
                 catinht=1.0
              else
                 catinht=(cn(l,nspezx)*acc(l,nspezx))**eac(i,l)
              endif
              catinh(i)=catinh(i)*catinht
 430       continue

c          Rate in (mol/dm**3rock/sec)
           ratem(i)=amin(i,nspezx)*km(i)*catinh(i)*fg(i)*1.0d-3
           rout(i,nspezx)=ratem(i)
              
           re(j)=re(j)-ss(j,i)*(po(i,nspezx)+ratem(i)*texe/por(nspezx))

           endif       
c<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< kinet
c        re(j)=re(j)-ss(j,i)*pn(i,nspezx)   !kinet
  15   continue
c

c
   20 continue
c      do 2223 ij=1,m1+m2
c2223      write(*,*)ij,re(ij)
c      pause
   22 continue
c
c  provision for use of a charge balance to specify one of the 
c  concentrations
c
      if (lchrg.eq.0) go to 30
      re(lchrg)=0.
      do 25 j=1,m1
        if (j.eq.lne) go to 25
        re(lchrg)=re(lchrg)-bn(j,nspezx)*vjb(j)
   25 continue
      if (m2.eq.0) go to 30
      do 26 i=1,m2
   26   re(lchrg)=re(lchrg)-cn(i,nspezx)*vjc(i)
c
   30 continue
      if (m3.eq.0) go to 41
c
c  *******************
c  residues for solids
c  *******************
      do 40 k=1,m3
        i=m1+k
c
c  calculate q, activity products of basis species 
c  in equilibrium with solids
c
c>>>>>>>>> kinet to below under 25
c04        prod=1.
c04        do 36 l=1,m1
c04         if (ss(l,k).eq.0.or.bn(l,nspezx).eq.0.) go to 36       ! 2003    .or.bn(l,nspez).eq.0.
c04         prod=prod*(acb(l,nspezx)*bn(l,nspezx))**
c04     *        ss(l,k)
c          if (nspezx.eq.2.and.dabs(prod).gt.1.e+50)then
c       write(*,'(a4,1x,i2,1x,3(e10.4,1x),i2)')'prod',
c     * l,prod,bn(l,nspezx),acb(l,nspezx),ss(l,k)
c           pause
c          endif
c04   36   continue
c04        q(k,nspezx)=prod
c<<<<<<<<<<<<<<<<<<< kinet to below under 25
c
        if (ikin(k).ne.0) goto 40  !kinetics  
c04        re(i)=prod-eqconst(k,nspezx)
        re(i)=q(k,nspezx)-eqconst(k,nspezx)  !IMPORTANT, now q, not prod.
        if ((pn(k,nspezx).le.0.).and.(re(i).le.0.))then
           re(i)=pn(k,nspezx)
c          if (nspezx.eq.2)write(*,*)k,re(i),prod,pn(k,nspezx)
        endif
c       if (nspezx.eq.2) then
c        write(*,'(a20,1x,i2,1x,4(e10.4,1x))') 'resolids m3 pn prod re'
c     *   ,k,pn(k,nspezx),prod,re(i),eqconst(k,nspezx)
c        endif

  40  continue
c
   41 continue
c
c  *************************************
c  calculate elements of jacobian matrix
c  *************************************

c not for every iteration a new jacobian
c      if(mod(nfail2,2).ne.0)goto 4444
c      nfail2=nfail2+1
c
c  initialize elements to zero
c
      do 43 i=1,m4
      do 43 j=1,m4
c                   ^^ m4 ersetzt durch 30 ,m4 hat mit dissolution mod. wert 4 !!
      z(i,j)=0.
   43 continue
c
      do 52 i=1,m1
c
c  basis rows, basis columns
c
      do 49 j=1,m1
      zz=0.
      if (i.eq.j) zz=-1.
      if (i.eq.lchrg) zz=-vjb(j)
      if (m2.eq.0) go to 47
      sum=0.
      do 46 k=1,m2
      coef=s(i,k)
      if (i.eq.lchrg) coef=vjc(k)
      if ((s(j,k).eq.0).or.(coef.eq.0.)) go to 46
      mk=m3+k
      prod=eqconst(mk,nspezx)*coef*s(j,k)*
     *acb(j,nspezx)/acc(k,nspezx)

      do 44 l=1,m1
      if ((l.eq.j).or.(s(l,k).eq.0.).or.
     * (bn(l,nspezx).eq.0..and.s(l,k).lt.0.)) go to 44   ! 2003 .or.bn(l,nspezx).eq.0.)
      prod=prod*(acb(l,nspezx)*bn(l,nspezx))
     ***s(l,k)
   44 continue
      if((s(j,k).eq.0).or.(bn(j,nspezx).eq.0.))goto 45
      ism1=s(j,k)-1
      prod=prod*(acb(j,nspezx)*bn(j,nspezx))
     ***ism1
   45 sum=sum+prod

      dxdc(k,l)=prod/coef   !dXk/dCl    kinet

   46 continue
      zz=zz-sum
   47 z(i,j)=zz
c       if(dabs(z(i,j)).gt.1000)write(*,*)i,j,zz

   49 continue
c
      if (m3.eq.0) go to 52
c
c  basis rows, solids columns
c
      do 50 k=1,m3
      if (ikin(k).ne.0) goto 50     !kinetics     
      j=m1+k
      z(i,j)=-ss(i,k)

   50 continue

   52 continue

      if (m3.eq.0) go to 60
c
      do 59 k=1,m3
      i=m1+k
c      if(nspezx.eq.2.and.dabs(pn(1,2)).gt.0.)write(*,*)
c     * k,pn(k,nspezx),q(k,nspezx),eqconst(k,nspezx)
      if ((q(k,nspezx).le.eqconst(k,nspezx)).and.
     *(pn(k,nspezx).le.0.)) go to 57
c
c  solids rows, basis columns
c
      do 56 j=1,m1
c >>> if (ss(j,k).gt.0) go to 53 <<< correct 9/25/87
      if (ss(j,k).ne.0) go to 53
      prod=0.
      go to 55
   53 prod=ss(j,k)*acb(j,nspezx)
      do 54 l=1,m1
      if ((l.eq.j).or.(ss(l,k).eq.0)) go to 54
      prod=prod*(acb(l,nspezx)*bn(l,nspezx))
     *      **ss(l,k)
   54 continue
      if((ss(j,k).eq.0).or.(bn(j,nspezx).eq.0.))goto 55
      ism1=ss(j,k)-1
      prod=prod*(acb(j,nspezx)*bn(j,nspezx))
     *     **ism1
c04  55 z(i,j)=prod  kinet
   55 if (ikin(k).eq.0) z(i,j)=prod  !Only if equilibrium  kinet

      dqdc(k,j)=prod    !dqmin/dcj    kinet

   56 continue
      go to 59
c
c  solids rows, solids columns
c
   57 z(i,i)=1.
c
   59 continue

c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> kinet
c     Rate terms in the jacobian matrix (basis rows-basis columns)
c     dRatem/dCj          -pep-

      do 450 k=1,m3
         if (ikin(k).eq.0) goto 450
         if (omega(k).le.1.0 .and. po(k,nspezx).le.0.0) then
            do 452 jjj=1,m1
               dratemdc(k,jjj)=0.0
 452        continue
            goto 450
         endif
         sdcatinhdc=0.0
         sdfgdc=0.0

c        Calculate dP(a**n)/dc terms

         do 500 jjn=1,m1  !First primary (basis) species
            if (eab(k,jjn).eq.0.0) goto 500
            product=1.0
            do 510 ii=1,m2  !product of the secondary species terms
               if (eac(k,ii).eq.0.0) goto 510
               product=product*(cn(ii,nspezx)*acc(ii,nspezx))
     .              **eac(k,ii)
 510        continue
            do 520 jj=1,m1
               if (jj.eq.jjn) goto 520
               if (eab(k,jj).eq.0.0) goto 520
               product=product*(bn(jj,nspezx)*acb(jj,nspezx))
     .              **eab(k,jj)
 520        continue
            otherj(jjn)=product*acb(jjn,nspezx)**eab(k,jjn)
     .           *eab(k,jjn)*bn(jjn,nspezx)**(eab(k,jjn)-1.0)
 500     continue

         do 530 iin=1,m2  !Now secondary (complex) species
            if (eac(k,iin).eq.0.0) goto 530
            product=1.0d0
            do 540 jj=1,m1 !product of the primary species terms
               if (eab(k,jj).eq.0.0) goto 540
               product=product*(bn(jj,nspezx)*acb(jj,nspezx))
     .              **eab(k,jj)
 540        continue
            do 550 ii=1,m2
               if (ii.eq.iin) goto 550
               if (eac(k,ii).eq.0.0) goto 550
               product=product*(cn(ii,nspezx)*acc(ii,nspezx))
     .              **eac(k,ii)
 550        continue
            otheri(iin)=product*acc(iin,nspezx)**eac(k,iin)
     .           *eac(k,iin)*cn(iin,nspezx)**(eac(k,iin)-1.0)

 530     continue


c  **** Now calculate f(dG) and cat/inh terms ******* -pep-

         do 460 j=1,m1

c          Delta G terms
            if (omega(k).lt.1.0d-15) then   !if (dg < -20 kcal/mol)
               dfgdc=0.0
            else if (ifg(k).eq.1) then       !f(dG)=|omega**a-1|**b
               oam1=omega(k)**a(k)-1.0
               dfgdc=b(k)*dabs(oam1)**(b(k)-1.0)*
     .          a(k)*omega(k)**(a(k)-1.0)
     .          /eqconst(k,nspezx)*dqdc(k,j)
               if (oam1.lt.0.0) dfgdc=-dfgdc
            else if (ifg(k).eq.2) then  !f(dG)=|ln omega|**a
               dfgdc=a(k)*dabs(dlog(omega(k)))**(a(k)-1.0)
     .          /omega(k)/eqconst(k,nspezx)*dqdc(k,j)
               if (dlog(omega(k)).lt.0) dfgdc=-dfgdc
            endif


c          cat/inh terms
            term=0.0
            do 470 jjn=1,m1
c              First the primary species terms
               if (eab(k,jjn).eq.0.0) goto 470
               term=term+otherj(jjn)  !dP/dcj
               do 480 iin=1,m2
c                 Now secondary species terms
                  if (s(jjn,iin).eq.0.0) goto 480
                  if (eac(k,iin).eq.0.0) goto 480
                  term=term+otheri(iin)*dxdc(iin,jjn) !dP/dXi * dXi/dCj
 480           continue
 470        continue
            
            sdcatinhdc=sdcatinhdc+term
            sdfgdc=sdfgdc+dfgdc

 460     continue

         dratemdc(k,j)=amin(k,nspezx)*km(k)*1.0d-3*
     .    (sdcatinhdc*fg(k) + catinh(k)*sdfgdc)

 450  continue


c     Add terms to basis rows-basis columns elements of jacob. matrix
      do 492 i=1,m1
         do 494 j=1,m1
            do 496 k=1,m3
               if (ikin(k).eq.0) goto 496

c              Kinetics
               z(i,j)=z(i,j)-ss(i,k)/por(nspezx)*texe*dratemdc(k,j)

c              Old mineral content.  This is a constant !!!
c               z(i,j)=z(i,j)-ss(i,k)*po(k,nspezx)

 496        continue
 494     continue
 492  continue

c<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< kinet
c
c  *******************************************************
c  calculate corrections and new values of concentrations;
c  check convergence
c  *******************************************************

 4444 continue
c        do 4434 ik=1,m1+m3
c 4434    write(*,'(10(1x,e10.4))')(z(ik,ij),ij=1,10)
c        write(*,*)'nfail2',nfail2,nspezx, iter
c        pause

c       nfail2=nfail2+1
   60 call simq(m4,re,z,n,isimflag)
c      if(isimflag.eq.1) return
c
c
c        do 4711 ikk=1,m1+m3
c          if (nspezx.eq.2.and.dabs(re(ikk)).gt.1.e+50)then
c       write(*,'(a4,1x,i2,1x,1(e10.4,1x))')'simq',
c     * ikk,re(ikk)
c       write(*,*)(bn(ikm,2),ikm=1,6)
c       write(*,*)(pn(ikm,2),ikm=1,3)
c       write(*,*)
c       do 4712 ikn=1,m1+m3
c       write(*,'(9(e10.2,1x))')(z(ikn,ikm),ikm=1,m1+m3)
c 4712  continue
c           pause
c          endif
c 4711 continue
      if (iter.gt.in2-1.and.iterneu.ge.10) write (6,101) nspezx,iter
 101  format(4x,'node=',i4,5x,'iter=',i5,5x,'values for this iteration'
     1,/,4x,'no species',2x,'old value',2x,'new value',2x,'error',2x,
     2'last old',2x,'total',/)
      nfail=0
      do 61 j=1,m1
      cc1=bn(j,nspezx)
      bn(j,nspezx)=(cc1-re(j))
c 110193       
c      if (bn(j,nspezx).lt.zero) bn(j,n)=0.
      if (bn(j,nspezx).lt.zero) bn(j,nspezx)=zero                  !2003  bn(j,nspezx)

c04      IF(dabs(bn(j,nspezx)).gt.1)then
c      WRITE(*,*)'NONSE BN(J,N)>1 at',bn(j,nspezx),J,N,iter
c04       bn(j,nspezx)=bo(j,nspezx)/(2.*2.)
c04      endif
      IF(bn(j,nspezx).le.0)then
c        new automatic guess
         bn(j,nspezx)=bo(j,nspezx)/1.
c       ^ 250495
      endif

      if (iter.gt.in2-1.and.iterneu.ge.10) then
c                   ^ 170594      
        write (6,102)dumb(j),j,cc1,bn(j,nspezx),re(j),
     * bo(j,nspezx),bi(j)
 102    format(a10,1x,i2,8(2x,e8.2))
      endif
      if (cc1.eq.0.) go to 61
      aaa=dabs(re(j))
      bbb=cc1*err
      if (iter.gt.in2-1.and.iterneu.ge.10) then
      write(*,'(a27,3(2x,e8.2))')'err aaa=abs(re) bbb=err*bn',
     *          err,aaa,bbb
      endif
      if (aaa.gt.bbb) then
        nfail=nfail+1
        nfailba_max=j
c max(nfailba_max,j)
c        write(*,*)'nfailba', nfailba_max
      endif
   61 continue

      if (m3.eq.0) go to 63
      do 62 j=1,m3
      if (ikin(j).eq.0) then  !equilibrium   kinet
      k=m1+j
      cc1=pn(j,nspezx)
      pn(j,nspezx)=cc1-re(k)
      if(pn(j,nspezx).le.0.or.eqconst(j,nspezx).eq.rlarge)
     * pn(j,nspezx)=0.
c                          ^new 10_96 
      if (iter.gt.in2-1.and.iterneu.ge.10) then
c                   ^ 170594
       xxxx=dlog10(eqconst(j,nspezx))
       write (6,104)dump(j),k,cc1,pn(j,nspezx),re(k),
     *       po(j,nspezx),cs(nspezx),xxxx
  104           format(a10,1x,i2,6(1x,e10.4))
      endif

      if (cc1.eq.0.) go to 62
      aaa=dabs(re(k))
      bbb=cc1*err
      if (iter.gt.in2-1.and.iterneu.ge.10) then
      write(*,'(a27,3(1x,e10.4))')'err aaa=abs(re) bbb=err*pn',
     *          err,aaa,bbb
      endif 
      if (aaa.gt.bbb)then
         nfail=nfail+1
        nfailso_max=j
c max(nfailso_max,j)
      endif
      else         ! kinet
         if (iter.gt.in2-1.and.iterneu.ge.10) then
c                   ^ 170594
            write (6,104)dump(j),k,cc1,pn(j,nspezx),re(k),
     *       po(j,nspezx)
        endif
      endif       !kinet

   62 continue

c**************************************************
c     User input into iteration
      if (iter.gt.in2-1) then
c<<<<<
       if(iterneu.lt.10)then
c        bn(2,nspezx)=1.
c        do 1323 ierr=1,m1
          if(mod(iterneu,2).eq.0)then
             write (*,*)'changen-x0',nfailba_max,
     *            bn(nfailba_max,nspezx-1)
c2004xx             bn(nfailba_max,nspezx)=bn(nfailba_max,nspezx-1)/2.
           do 1818 ib=1,m1
 1818		 bn(ib,nspezx)=bn(ib,1)  
           do 1817 ib=1,m2
 1817		 cn(ib,nspezx)=cn(ib,1)  
           do 1816 ib=1,m3
 1816		 pn(ib,nspezx)=pn(ib,1)  
          endif
          if(mod(iterneu,2).eq.1)then
c2004xx             bn(nfailba_max,nspezx)=bn(nfailba_max,nspezx+1)/2.
           do 1819 ib=1,m1
 1819		 bn(ib,nspezx)=bn(ib,nnodex)  
           do 1820 ib=1,m2
 1820		 cn(ib,nspezx)=cn(ib,nnodex)  
           do 1821 ib=1,m3
 1821		 pn(ib,nspezx)=pn(ib,nnodex)  
             write (*,*)'changen+nnodex',nfailba_max,
     *                      bn(nfailba_max,nspezx+1)
          endif
 1323    continue
c  ^  240495
        write(*,'(1x,a9,i4,2(1x,a8,i3))')'iterneu= ',
     *        iterneu,'nspezx= ',nspezx
        iterneu=iterneu+1
        iter=0
       else  
c>>>>> 170594
         write(*,*)'try a new guess for a basic species j  (y/n)?'
     *     ,'fail ',nfailba_max,nfailso_max
         write(*,*)'nespex= ',nspezx
         read(*,*)ssw
         if(ssw.eq.'y'.or.ssw.eq.'Y')then
                write(*,*)'give:  j,bn(j,nspezx)'
                read(*,*)j,bn(j,nspezx)
                iter=0
         else
                write(*,*)'try a new guess for a solid k= ?'
                read(*,*)ssw
                if(ssw.eq.'y'.or.ssw.eq.'Y')then
                   write(*,*)'give:  k,pn(k-m1,nspezx)'
                   read(*,*)k,pn(k-m1,nspezx)
                   iter=0
                endif
                 
         endif
         
c<<<<<  
       write(*,*)'Change error tolerance ?'
       read(*,*)ssw
       if(ssw.eq.'y'.or.ssw.eq.'Y')then
         write(*,*)'give:  error tol'
         read(*,*)err
         errflag=1
       endif
                    
       iterneu=0
       Write(*,*)'iterneu= 0 gesetzt'
c       iter=0
       endif
c>>>>>> 170594
      endif


   63 if (m2.eq.0) go to 66
      do 65 j=1,m2
      i=m3+j
      prod=1.
      do 64 l=1,m1
c      if (s(l,j).eq.0.) goto 64      
      if (s(l,j).eq.0.or.bn(l,nspezx).eq.0) go to 64
c                     .or.bn(l,nspezx).eq.0
c                     ^new 04 11 93
c                     ^new 08 02 95
      prod=prod*(acb(l,nspezx)*bn(l,nspezx))**s(l,j)
   64 continue
      cn(j,nspezx)=eqconst(i,nspezx)*prod/
     *acc(j,nspezx)
c      if (cn(j,nspezx).lt.zero) cn(j,nspezx)=zero
      if (cn(j,nspezx).lt.zero) cn(j,nspezx)=zero
   65 continue


c
   66  if (nfail.gt.0) go to 10

c<<<<<<
       if(errflag.eq.1)then
c      ^191094
         write(*,*)'Change error tolerance back again? Y/N'
         read(*,*)ssw
         if(ssw.eq.'y'.or.ssw.eq.'Y')then
           write(*,*)'give:  error tol'
           read(*,*)err
           errflag=0
         endif
       endif
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  kinet
c     Update of min. content for kinetics 
c2004xx      pornew=1.0
      do 699 j=1,m3
       if (ikin(j).eq.1) then  !Kinetics   -pep-
c       Use rate to update the mineral content    -pep-
        pn(j,nspezx)=po(j,nspezx)+ratem(j)/por(nspezx)*texe
c       dpdp=pn(j,nspezx)-po(j,nspezx)
     
c       if(nspezx.eq.2)write(*,'(1x,i2,3(e10.4))')
c     *                       j,dpdp,por(nspezx),poro(nspezx)
       endif
       if (pn(j,nspezx).lt.0.0) pn(j,nspezx)=0.0
       volfr=pn(j,nspezx)*volmol(j)*por(nspezx)*1.0d-3
c       volfr_o=po(j,nspezx)*volmol(j)*por(nspezx)*1.0d-3
       vout(j,nspezx)=volfr
       amin(j,nspezx)=(4.0*pi*rng(j))**(1.0/3.0)*(3.0*volfr)**(2.0/3.0)
c       amin_o=(4.0*pi*rng(j))**(1.0/3.0)*(3.0*volfr_o)**(2.0/3.0)

c2004xx       pornew=pornew-vout(j,nspezx)  !Update in porosity  -pep-
c        write(*,'(a1,1x,i2,1x,i2,5(1xe10.4))')'1',nspezx,j,
cc     *  pn(j,nspezx),po(j,nspezx),ratem(j),por(nspezx),texe
c      write(*,*)'2',volfr,pn(j,nspezx),volmol(j)
c      write(*,*)'3 amin ',(amin(j,nspezx),j=1,6),amin_o
c        pause
c      endif
 699  continue
c      write(*,*)'3 amin ',(amin(j,nspezx),j=1,6)
c        pause
c      if (pornew.gt.0.2)then
c        write(*,*)nspez,por(nspezx),poro(nspezx),pronew
c        pause
c      endif
c      por(nspezx)=pornew
c2004xx      vx(nspezx)=vxx/por(nspezx)  !Update in pore velocity   -pep-
c      write(*,*)'n 4vx',nspezx,vx(nspezx)
c      pause
c<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<      kinet

      return
   69 iterm=iter-1
      write (6,100) nspezx,iterm
 100  format(1x,'no convergence in stati2 at node',i4,
     *' after iter ='i6)
      write (6,*) nspezx,iterm
      istat2flag=1
      return
c      call datout
c      stop

      end
c***********************************************************************
c***********************************************************************
c

