
c ***********************************************************************
c ***********************************************************************
c
      subroutine actco(nspezx,num,bn,cn,vjb,vjc,acb
     *,acc,lne)
      implicit double precision (a-h,o-z)
      include 'gwheader.inc'
c
c  **********************************************************
c  calculates activities based on the davies equation.
c  the activities are stored in the matrices acb and acc.
c  **********************************************************
c
      dimension num(6),bn(nbasis,nnodex),cn(ncompl,nnodex)
     1,vjb(nbasis),vjc(ncompl),
     1acb(nbasis,nnodex),acc(ncompl,nnodex)

c      iactco=iactco+1
c      write(*,*)'subroutine actco aufgerufen',iactco,nspezx,nspezy
c     *,nbegin,nend
c
      m1=num(1)
      m2=num(2)
c      do 50 n=nbegin,nend
c      do 50 ny=1,1
c
c  calculation of ionic strength
c
      ctot=0.
      si=0.
      do 10 j=1,m1
c       write(*,*)'actcom1',j,bn(j,n),vjb(j)
        if (j.eq.lne) go to 10
          
        cc=(bn(j,nspezx))
c        if(cc.gt.2.)write(*,*)'j= ',j,'cc =',cc

        vjj=vjb(j)
c       write(*,*)'j= ',j,'cc =',cc,'ctot= ',ctot
        ctot=ctot+vjj*cc
        si=si+cc*vjj*vjj
   10 continue
      if (m2.eq.0) go to 16
      do 15 i=1,m2
c       write(*,*)'actco m2',i,cn(i,n),vjc(i)

        vjj=vjc(i)
        cc=cn(i,nspezx)
c        if(cc.gt.2.)write(*,*)'i= ',i,'cc =',cc
        ctot=ctot+vjj*cc
c2004xx       ctot=0.
c       ^Fiona MINEQL compatible ???
c        write(*,*)'i= ',i,'cc =',cc,'ctot= ',ctot
        
        si=si+cc*vjj*vjj
   15 continue
c
c  although the net charge of the species being considered in the
c  calculations may not be equal to zero, for an accurate calculation
c  of the ionic strength, all the chemical species must be considered.
c  therefore, it is assumed here that the balance of the charge is made
c  up by a chemical species of a +1 or -1 charge and a concentration
c  of ctot.  if this is not reasonable, then the actual calculation
c  should be done with a net charge of 0.
c
   16 si=si+dabs(ctot)
      cionst=0.5*si
      
c      write(*,*)'n ionic strength',n,cionst

      si=sqrt(cionst)
      fi=0.5*(si/(1.+si)-0.3*cionst)
c                         ^ 0.2 damit 06.04.94 MINEQL (fiona compatible)      
c
c  calculation of activity coefficients
c
      do 20 j=1,m1
        vjj=vjb(j)
        if (j.eq.lne) vjj=0.
c         if(g.lt.0)write(*,*)'actco fi<0 ',fi,cionst,vjj
        g=fi*vjj*vjj
        ag=dabs(g)
        if (ag.gt.10.) then
           g=0.
c           write(*,*)'g for i= ',j,' >10'
        endif
        if (ag.lt.0.) then
           g=0.
         write(*,*)'ag lower 0 for i= ',j,(bn(iii,nspezx),
     * iii=1,m1)
        pause
        endif
        acb(j,nspezx)=10.**(-g)

   20 continue
      if (m2.eq.0) go to 50
      do 25 j=1,m2
        vjj=vjc(j)
        gcx=fi*vjj*vjj
c        if(gcx.lt.0)write(*,*)'actco fi<0 ',fi,cionst,vjj
        agcx=dabs(gcx)
        if (agcx.gt.10.)then
           gcx=0.
c           write(*,*)'agcx for i= ',j,' >10'
        endif   
        if (agcx.lt.0.) then
           gcx=0.
           write(*,*)'agcx lower 0 for i= ',j
        endif
        acc(j,nspezx)=10.**(-gcx)
 
   25 continue
   50 continue
      return
      end
c

c***********************************************************************
c***********************************************************************
c
      subroutine ehcalc(n1,n2,lne,tmp,bn,eh)
      implicit double precision (a-h,o-z)
cc  equilibrium constants are calculated from the expression
c
c           log k = c1 + c2/tk + c3*(log tk) + c4*tk
c
c  where k is the equilibrium constant at absolute temperature tk and the
c  constants c1,...,c4 are input items.  For isothermal simulations the
c  appropriate log k value can be input for c1, and the value zero input
c  for c2,...,c4.
c
c  variable oxidation potential is treated by defining a hypothetical
c  electron activity to be a basis species.  The basis species for a
c  multivalent element is taken to be the chemical form of the highest
c  oxidation state of the element.  Chemical forms in lower oxidation
c  states are realized by combining the oxidized basis species with
c  electrons in a manner analogous to the formation of aqueous complexes.
c  Input concentrations of electrons are the activities related to
c  oxidation potential (eh) by
c
c                eh(volts) = -1.9841412e-04*tk*log(e-)
c
c  where tk is the absolute temperature and (e-) is the hypothetical
c  electron activity.  Variable oxidation potentials in thcc can be
c  turned off by specifying input parameter lne = 0 and not identifying
c  (e-) as a basis species.
c
c  *********************************************************************
c  subroutine ehcalc calculates the oxidation potential (eh) of the 
c  fluid at temperature tmp and current electron activity.
c  eh is computed from
c
c                 eh = -1.9841412e-04*tk*log(e-)
c
c  where tk = tmp + 273.15 and (e-) is the electron activity.
c  *********************************************************************
c
      include 'gwheader.inc'

      dimension tmp(nnodex),bn(nbasis,nnodex)
     1,acb(nbasis,nnodex),eh(nnodex)
      data a1,tcel /1.9841412e-04,273.15/
c
      iehcalc=iehcalc+1
      write(*,*)'subroutine ehcalc aufgerufen',iehcalc, n1,n2
c


      do 10 node=n1,n2
        tk=tmp(node)+tcel
        bnl=bn(lne,node)
        if (bnl.le.0.) write(6,20) bnl,node
        eh(node)=-a1*tk*dlog10(bn(lne,node))
c2003        eh(node)=-a1*tk*log10(bn(lne,node))
   10 continue
      return
   20 format(/,5x,'nonphysical electron activity =',1pe12.4,2x,'at node'
     1,i4)
      end

c ---------------------------------------------------------------------- 
      subroutine eqconcs (nspezx,num,cs,eqconst,pn)

      implicit double precision (a-h,o-z)
c
c ---------------------------------------------------------------------- 
csubroutine calculates log Kso as a f (C/S), U. Berner, see TM-45-87-10
c
c in Berner's model, incongruent dissolution of hydrated calcium-
c silicates (cement) is modelled by assuming that the cement is 
c comprised of two solids, each of which dissolves congruently.  what 
c the two solids are depends on the ratio of calcium to silicon (C/S)
c in the solid.  
c
c for each solid that dissolved congruently, 
c
c        log Kso = log a^l - log a^s
c
c where a^l is the activity of the component in the fluid phase
c (actually written as a sum of the activities of the ions that form
c it) and a^s is the activity of the solid phase, which is not 
c necessarily unity.  Berner defines 
c        log Kso = log Kso + log a^s = log a^l 
c the assumptions made in writing this subroutine are described in the
c program header.
c N.B. in the code below, the ^Kso values of components that are not
c      present in a given region are set to large values in order to 
c      prevent the solid from forming.
c ---------------------------------------------------------------------- 

      include 'gwheader.inc'

      dimension num(6),cs(nnodex),eqconst(ncompl+nsolid,nnodex)
     1,pn(nsolid,nnodex)
      common /dismdl/ idismdl,indxca,indxsi,indxcs,
     1                csbndry(0:2),conhcsm(0:3,2,4)
c >>>                              ^ named shortened from conhcsm in main
c
      ieqconcs=ieqconcs+1
c      write(*,*)'subr. eqconcs aufgerufen',ieqconcs,nspezx
c
c      rlarge = 99.
c      zero = 1.0e-37
c
c ------------------------------------
c calculate C/S at nodes n1 through n2
c ------------------------------------
         pncan = pn(indxca,nspezx)
         pnsin = pn(indxsi,nspezx)
         pncsn = pn(indxcs,nspezx)
cccccccccc
         if((pncan.eq.0.).and.(pnsin.eq.0.)
     *        .and.(pncsn.eq.0.)) then
c           if ((pncan.eq.0.).and.(pnsin.eq.0.)) then
c ----------if the last condition is used there is no switching to SiO2 !!
c --------- all solid concentrations are zero
c --------- set cs(n) to a negative value to use as a flag below
            cs(nspezx) = -999.
            icscase=1
         else
           if (pnsin.le.0..and.pncsn.ne.0.) then
c ------------ regions II and III
               cs(nspezx)=(pncsn+pncan)/pncsn
               icscase=2
           endif
c ------------ region I
           if(pncsn.ne.0..and.pnsin.ne.0.) then
               cs(nspezx)=pncsn/(pncsn+pnsin)
               icscase=3
           endif
           if(pncan.le.0.and.pncsn.gt.0.and.pnsin.le.0.) then
               cs(nspezx)=1.
               icscase=4
           end if
           if(pncan.le.0.and.pncsn.le.0..and.pnsin.gt.0.) then
               cs(nspezx)=0.
               icscase=5
            end if
         end if
cccccccccc
c         end if
c         write(*,*)nspezx,'icscase    ',icscase
c -------------------
c calculate log^K_so
c -------------------
         if (cs(nspezx).lt.0.) then
c --------- all solid concentrations are zero
c --------- set all equilibrium constants to constant values
            term = conhcsm(0,1,1)
            eqconst(indxsi,nspezx) = rlarge
             term = conhcsm(3,1,1)
            eqconst(indxca,nspezx) = 10.**term
            term = conhcsm(3,2,1)
            eqconst(indxcs,nspezx) = 10.**term
         else
            if (cs(nspezx).gt.csbndry(2)) then
c ------------ region III: Ca(OH)2 (comp 1) and CaH2Si04 (comp 2) present
               term = conhcsm(3,1,1)
               eqconst(indxca,nspezx) = 10.**term
               term = conhcsm(3,2,1)
               eqconst(indxcs,nspezx) = 10.**term
               eqconst(indxsi,nspezx) = rlarge
            else
             if (cs(nspezx).gt.csbndry(1)) then
c --------------- region II: Ca(OH)2 (comp 1) and CaH2Si04 (comp 2) present
               term=conhcsm(2,1,1)+conhcsm(2,1,2)/(cs(nspezx)+
     1              conhcsm(2,1,3))
               eqconst(indxca,nspezx) = 10.**term
               term = conhcsm(2,2,1) 
               eqconst(indxcs,nspezx) = 10.**term
               eqconst(indxsi,nspezx) = rlarge
              else
               if (cs(nspezx).gt.csbndry(0)) then
c ------------------ region I: SiO2 (comp 1) and CaH2Si04 (comp 2) present
                  csn = cs(nspezx)
                  term = conhcsm(1,1,1) + conhcsm(1,1,2)/
     1                  (csn + conhcsm(1,1,3))
                  eqconst(indxsi,nspezx) = 10.**term
                  term = conhcsm(1,2,1) 
     1                   - ((1. - csn)/csn)
     2                   *(conhcsm(1,2,2) + conhcsm(1,2,3)
     3                   /(csn + conhcsm(1,2,4)))
                  eqconst(indxcs,nspezx) = 10.**term
                  eqconst(indxca,nspezx) = rlarge
                else
c ------------------ C/S = 0: only SiO2 present  
                  term = conhcsm(0,1,1) 
                  eqconst(indxsi,nspezx) = 10.**term
                  eqconst(indxca,nspezx) = rlarge
                  eqconst(indxcs,nspezx) = rlarge
                end if
             end if
           end if
         end if
        if (cs(nspezx).eq.0.)then
c ------------------ C/S = 0: only SiO2 present  
                  term = conhcsm(0,1,1) 
                  eqconst(indxsi,nspezx) = 10.**term
                  eqconst(indxca,nspezx) = rlarge
                  eqconst(indxcs,nspezx) = rlarge
        endif
   10 continue
      return
      end

c ----------------------------------------------------------------------
c >>> subroutine eqcon rewritten 9/87 and renamed eqcont.  the argument 
c >>> list was changed too.  m6 was replaced by num and the last four
c >>> arguments are new.
c ----------------------------------------------------------------------
c
      subroutine eqcont (nspezx,con,num,lnh,tmp,
     1                   eqconst,itmpdep,tmpk,idismdl)
      implicit double precision (a-h,o-z)
c ----------------------------------------------------------------------
c subroutine eqcont calculates the dissociation constant of water, the
c formation constants of aqueous complexes, and the solubility products
c of solids at the temperature tmp.  there are two possible forms for
c calculating the constants.  which is used depends on the control
c variable itmpdep.
c
c if itmpdep equals 1 the constants are calculated from expressions of 
c the form
c
c            log K = c1 + c2/Tk + c3*(log Tk) + c4*Tk
c
c where c1, c2, c3 and c4 are calculated from data.
c
c if itmpdep equals 2, the constants are calculated from expressions of
c the form
c
c            log K(Tk) = log K(To) + beta(Tk)*dHo(To)
c
c where To is the reference temperature (K), dHo is the standard enthalpy
c of the reaction (kJ/mol) and 
c
c            beta(Tk) = (1/To - 1/Tk) / (R*ln10)
c
c in both expressions above Tk is the temperature on the Kelvin scale,
c
c            Tk = tmp + 273.15.
c
c in the original version of this subroutine (THCC), only log K
c expressions of the first form were included.
c ----------------------------------------------------------------------
c

      include 'gwheader.inc'

      dimension itmpdep(ncompl+nsolid),num(6)
     1,con(ncompl+nsolid,4),eqconst(ncompl+nsolid,nnodex)
     2,tmp(nnodex),tmpk(nnodex)

c --- betac = 1/(R*ln10)
      common /betat/ betac
c      betac = 1./(8.3143*dlog(10.d0))

c      ieqcont=ieqcont+1
c      write(*,*)'eqcont aufgerufen',ieqcont,nspezx,nspezy


c
      m2 = num(2)
      m3 = num(3)
      m5 = num(5)
c
c ------------------------
c convert from deg. C to K
c ------------------------
c      do 10 node = n1,n2
c      do 10 ny =1,1
c      write(*,*) node,con(j,1),con(j,2),tk,con(j,3)
         tmpk(nspezx) = tmp(nspezx) + 273.15
   10 continue
c
c ---------------------------------------------------
c calculate temperature-dependent formation constants 
c for the complexes.  
c ---------------------------------------------------
      if (m2.gt.0) then
         j1 = m3 + 1
         j2 = m2 + m3
         do 50 j = j1,j2
c      write(*,*)'step1'
c            if (itmpdep(j).eq.1) then
c                  write(*,*)'step4'
c               do 30 node = n1,n2
c               do 30 ny=1,1
                tk = tmpk(nspezx)
c                write(*,*) j,con(j,1),con(j,2),tk,con(j,3)
                term = con(j,1)+con(j,2)/tk+con(j,3)*dlog10(tk)
c2003                term = con(j,1)+con(j,2)/tk+con(j,3)*log10(tk)
     1                   +con(j,4)*tk
c                write(*,*)'step2',term
                  eqconst(j,nspezx) = 10.**term
c                write(*,*)'compl', j,nspezx,
c     *eqconst(j,nspezx),tk
   30          continue
c            else
c                  write(*,*)'step5'
c             do 40 node = n1,n2
c             do 40 ny = 1,2 
c              tk = tmpk(node)
cc               write(*,*)'kkk', j,con(j,1),con(j,2),tk,con(j,3)
c              term = con(j,1)+(1./con(j,2)- 
c     *        1./tk)*betac*con(j,3)
c              eqconst(j,node) = 10.**term
cc              write(*,*)'step3',term
c   40        continue
c            end if
   50    continue
      end if
c
c ---------------------------------------------------
c calculate temperature-dependent solubility products
c for solids not included in the dissolution model.
c ---------------------------------------------------
      m = m3 - m5
      
c      if(idismdl.gt.0)m=4
c      write(*,*) 'mmm',m,m3,m5
C     ^ damit solids und inkongruente dissolution  10 05 94

      if (m.gt.0) then
         do 80 j = (1+m5),(m3+m5)
c            if (itmpdep(j).eq.1) then
c               do 60 node = n1,n2
c               do 60 ny=1,1
c                write(*,*)'eqcont node, ny',node,ny
                tk = tmpk(nspezx)
c      write(*,*) j,con(j,1),con(j,2),tk,con(j,3)
                term = con(j,1)+con(j,2)/tk+con(j,3)*dlog10(tk)
c2003                term = con(j,1)+con(j,2)/tk+con(j,3)*log10(tk)
     1                +con(j,4)*tk
                  eqconst(j,nspezx) = 10.**term
c                write(*,*)'sol', j,nspezx,
c     *eqconst(j,nspezx),tk
   60          continue
c            else
c             do 70 node = n1,n2
c             do 70 ny = 1,2
c              tk = tmpk(node)
cc      write(*,*) j,con(j,1),con(j,2),tk,con(j,3)
c              term = con(j,1) + (1./con(j,2) - 1./tk)*betac*con(j,3)
c              eqconst(j,node) = 10.**term
c   70        continue
c            end if
   80    continue
      end if
c
      return
      end

c***********************************************************************
c***********************************************************************
c
      subroutine simq(n,y,z,nspez,isimflag)
      implicit double precision (a-h,o-z)
c
c  ****************************************************************
c  this subroutine solves the equation zx=y for x.
c  it is used in calculating an equilibrium distribution of species 
c  such as in the initial conditions and the boundary conditions.
c  the jacobian is stored in z and the residues are in y.  
c  ****************************************************************
c

      include 'gwheader.inc'

      dimension y(nbasis+nsolid),z(nbasis+nsolid,nbasis+nsolid)


c      isimq=isimq+1
c      write(*,*)'subroutine simq aufgerufen',isimq
c


c
c     provision for n=1
c
c      if(n.ne.1) go to 50
c      y(1)=y(1)/z(1,1)
c      return
c   50 continue
c
c     element of elimination
c
      NN=N
      n1=n-1
      do 10 m=1,n1
      zmax=0.
      imax=0
c
c     find max of column
c
      do 20 i=m,n

c       write(*,'(A15,2(1x,i2),2e10.2)')'i,m,z(i,m),zmax',i,m,
c     *z(i,m),zmax
c
        if(dabs(z(i,m)).lt.zmax) go to 20
        imax=i
        zmax=dabs(z(i,m))
c<<<   
c      if(imax.eq.0) then
c      imax=i
c      zmax=1.e-200
c      endif
c>>> ^ 11 05 1994
   20 continue
   
      if(imax.ne.0) go to 30
c
c     error return
c
c2003      write (6,1000)
c2003      write (6,900) m
c2003     write(*,*)'zmax= i n ',zmax,i,n
      DO 66 M1=1,NN
      write (6,*) (z(m1,I),i=1,n)
  66  CONTINUE
       write(*,*)'nspez =', nspez
       isimflag=1
       return
c
   30 continue
c
c     row interchange
c
      if(imax.eq.m) go to 35
      v=y(m)
      y(m)=y(imax)
      y(imax)=v
      do 40 j=m,n
        v=z(m,j)
        z(m,j)=z(imax,j)
        z(imax,j)=v
   40 continue
   35 continue
c
c     diagonalize
c
      m1=m+1

      do 70 i=m1,n
      if (z(i,m).eq.0.)then
         v=0.
      else
       if (z(m,m).eq.0.) then
c2003      write(*,*)'simq_dia z mm =0',m1,n,y(n),z(m,m)
       v=0.
       else 
       v=z(i,m)/z(m,m)
       endif
      endif
c     ^ if - endif noetig, da division nicht wie bei i860 mit software erfolgt  
        y(i)=y(i)-v*y(m)
        do 70 j=m,n
                z(i,j)=z(i,j)-v*z(m,j)
   70   continue

   10 continue
c
c     back substitute
c
c       write(*,*)'simq2 z mm =0',m1,m,y(n),z(n,n)
      if (z(n,n).eq.0.)then
c2003       write(*,*)'simq2 z mm =0',m1,m,y(n),z(n,n)
      else 
       y(n)=y(n)/z(n,n)
      endif
c     ^ if - endif noetig, da division nicht wie bei i860 mit software erfolgt  
      n1=n-1
      do 100 k=1,n1
        i=n-k
        i1=i+1
        do 90 j=i1,n
   90      y(i)=y(i)-y(j)*z(i,j)
        if(z(i,i).eq.0)then
          y(i)=0.
c2003         write(*,*)'simq4  y(i)=0 gesetzt'
        else
          y(i)=y(i)/z(i,i)
        endif
  100 continue
c      write(*,*)'ende simq'
      return
  900 format (1x,'m=',i10,/,4x,'(z(i,m),i=1,n) values are')
  990 format (1x,40(1X,1pe9.1))
 1000 format(1x,'singular jacobian matrix in simq')
      end

