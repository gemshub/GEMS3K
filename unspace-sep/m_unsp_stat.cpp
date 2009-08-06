#include <math.h>
#include <stdio.h>

#ifndef IPMGEMPLUGIN

#include "m_unspace.h"
#include "m_syseq.h"
#include "visor.h"

#else

#include "ms_unspace.h"
#include "nodearray.h"

#endif


// !!! internal using syp->GEX, syp->B, syp->Dcl, mup->Laq, mup->Pg, mup->Ll
// if( syp->PGEX != S_OFF )
//    if( fabs( syp->GEX[pmp->muj[j]] - pmp->GEX[j]*pmp->RT ) >= 0.001 )
//        break;
// tpp->G, tpp->S, tpp->Vm must be internal LagranInterp() build
// syp->Guns, syp->Vuns internal work arrays, set up 0 before calc

//=========================================================================
// Generation arrays part
//=========================================================================

#ifdef IPMGEMPLUGIN
  extern  bool load;
#endif

// make EqStat key  && calculate records
void TUnSpace::unsp_eqkey()
{

	double calculation_time;
    long int NumPrecLoops = 0, NumIterFIA = 0, NumIterIPM = 0;

#ifndef IPMGEMPLUGIN

    // calculate EqStat record (Thermodynamic&Equlibria)
    pmu->pTPD = 0;
    //     pmu->pNP = 0;
    vstr buf(40);

    sprintf(buf, "%.4d", usp->q);
    memset(usp->timep, 0, 5 );
    strncpy(usp->timep, buf, 4 );
    Gcvt( usp->Tc, 6,  usp->TCp );
    Gcvt( usp->Pc, 6,  usp->Pp );
    Gcvt( usp->Vc, 6,  usp->Bnamep );

 //  rt[RT_SYSEQ].MakeKey( RT_UNSPACE,  usp->stkey, RT_UNSPACE, 0, RT_UNSPACE,1,
 //       RT_UNSPACE, 2, K_IMM, usp->timep, K_IMM, usp->Bnamep,
 //       K_IMM, usp->Pp, K_IMM, usp->TCp, RT_UNSPACE, 7, K_END );
   rt[RT_SYSEQ].MakeKey( RT_UNSPACE,  usp->stkey, RT_UNSPACE, 0, RT_UNSPACE,1,
        RT_UNSPACE, 2, RT_UNSPACE, 3, K_IMM, usp->Bnamep,
        K_IMM, usp->Pp, K_IMM, usp->TCp, K_IMM, usp->timep, K_END );
   rt[RT_SYSEQ].Find(usp->stkey);

// calc current SyStat 16/02/2007
     pmu->TCc = usp->Tc;
     pmu->Pc = usp->Pc;
         
     calculation_time = TProfil::pm->calcMulti( NumPrecLoops, NumIterFIA, NumIterIPM );
    //TProfil::pm->CalcEqstat( false ); // 16/02/2007
// Later: to implement account for calculation time and numbers of iterations/loops

// Added 05.06.2008 by DK to save g0 variations done by UnSpace
// by incrementing them to GEX in SysEq records (needed for diagnostic SyEq records) 
  int k, jb, je, j, jj; 
  double Gg, Ge, Go, pGo;   
  for( k = 0, je=0; k < pmu->FI; k++ )
  {
   	jb = je;
   	je += pmu->L1[k];
   	for( j=jb; j<je; j++ )
    {    		
       Gg = 0.0;    //   This part had to be changed after integrating Ge into pmp->G0 
       jj = pmu->muj[j]; //        DK    07.03.2008,  16.05.2008
       Go = TProfil::pm->tpp->G[jj]; //  G0(T,P) value taken from MTPARM
       if( syu->Guns )  // This is used mainly in UnSpace calculations
           Gg = syu->Guns[jj];    // Increment to G0 set by UnSpace
// Here we determine what was the initial syp->GEX[jj] increment
       pGo = TProfil::pm->pmulti->Cj_init_calc( Go+Gg, j, k );
       Ge = ( pmu->G0[j] - pGo )* pmu->RT; 
       if( TProfil::pm->syp->GEX ) // Setting Gg + Ge to syp->GEX[jj]
    	   TProfil::pm->syp->GEX[jj] = Gg + Ge;   // Ge is part of pmp->G0 since 07.03.2008 (DK)  
//       if( syu->GEX && syu->PGEX != S_OFF )
//           syu->GEX[jj] = Gg + Ge;       
    } // j loop
  } // k loop	
    if( usp->PsSY != S_OFF )
       TSysEq::pm->CmSave();           // save SysEq record with results to DB
    if( usp->stl )
       memcpy( usp->stl+usp->q, usp->stkey, EQ_RKLEN );
 
    for( k = 0, je=0; k < pmu->FI; k++ )
    {
     	jb = je;
     	je += pmu->L1[k];
     	for( j=jb; j<je; j++ )
      {    		
         Gg = 0.0;    //   This part had to be changed after integrating Ge into pmp->G0 
         jj = pmu->muj[j]; //        DK    07.03.2008,  16.05.2008
         Go = TProfil::pm->tpp->G[jj]; //  G0(T,P) value taken from MTPARM
         if( syu->Guns )  // This is used mainly in UnSpace calculations
             Gg = syu->Guns[jj];    // Increment to G0 set by UnSpace
  // Here we determine what was the initial syp->GEX[jj] increment
         pGo = TProfil::pm->pmulti->Cj_init_calc( Go+Gg, j, k );
         Ge = ( pmu->G0[j] - pGo ) * pmu->RT; 
         if( TProfil::pm->syp->GEX ) // Setting Gg + Ge to syp->GEX[jj]
      	   TProfil::pm->syp->GEX[jj] = Ge;   
  //       if( syu->GEX && syu->PGEX != S_OFF )
  //           syu->GEX[jj] = Ge;       
      } // j loop
    } // k loop
    
#else
    // calculate EqStat record (Thermodynamic&Equlibria)
     load = false;
    
    TNode::na->pCNode()->NodeStatusCH = NEED_GEM_AIA; // activating GEM IPM for automatic initial
                                      // approximation
 // re-calculating equilibria by calling GEMIPM
    TNode::na->GEM_run( 1., false );
#endif    

}


// building arrays for make pay off matrix
void TUnSpace::buildTestedArrays()
{
 int i;
 short Ip;

#ifndef IPMGEMPLUGIN
 showMss = 1L;
#endif
 
 for( Ip=0; Ip<usp->Q; Ip++)
 {
    usp->q = Ip;

#ifndef IPMGEMPLUGIN

   pVisor->Message( window(), GetName(),
             "Generation of EqStat records\n"
                 "Please, wait...", usp->q, usp->Q);
#endif
   
 // setup uncertainty point data
     NexT( Ip );
     unsp_eqkey();
if(!Ip) 
	TProfil::pm->outMultiTxt( "ipm_1.txt"  );


 // set zeros for tested arrays
 if(usp->PsGen[0]== S_ON  )
 {   for( i=0; i<usp->L; i++)
    {
       usp->vY[Ip*usp->L+i]  = 0.;
       usp->vGam[Ip*usp->L+i]= 0.;
       usp->vG[Ip*usp->L+i]= 0.;
    }
    for( i=0; i<usp->Fi; i++)
         usp->vYF[Ip*usp->Fi+i] = 0.;
    for( i=0; i<usp->N; i++)
    {
       usp->vU[Ip*usp->N+i] = 0.;
       usp->vMol[Ip*usp->N+i] = 0.;
    }
    if(usp->Ls/*TProfil::pm->mup->Pg*/)
       for( i=0; i<usp->Ls; i++)
            usp->vFug[Ip*usp->Ls+i] = 0.;

 // copy data from multy
    for( i=0; i<pmu->L; i++)
    {
      usp->vY[Ip*usp->L+(pmu->muj[i])]  =
                                         pmu->X[i];
      usp->vGam[Ip*usp->L+(pmu->muj[i])]=
                                         pmu->lnGam[i];

//      usp->vG[Ip*usp->L+(pmu->muj[i])]=
//                                     pmu->G[i];
     }

     for( i=0; i<pmu->L; i++)
     { int ii = (pmu->muj[i]);
#ifndef IPMGEMPLUGIN
       double xx = (double)(syu->Guns[ii]);
       xx += (double)(TProfil::pm->syp->GEX[ii]);  // Added by DK 05.06.2008
//            xx += (float)(pmu->GEX[i])/pmu->RT; // syu->GEX[ii]
       xx += TProfil::pm->tpp->G[ii];
#else
       double xx = (pmu->Guns[ii]);
//            xx += (float)(pmu->GEX[i])/pmu->RT; // syu->GEX[ii]
            xx += pmu->tpp_G[ii];
#endif
      usp->vG[Ip*usp->L+ii]= xx;
     }

    for( i=0; i<pmu->FI; i++)
       usp->vYF[Ip*usp->Fi+(pmu->muk[i])] =
                                           pmu->XF[i];

    for( i=0; i<pmu->N; i++)
    {
       usp->vU[Ip*usp->N+(pmu->mui[i])] = pmu->U[i];
       usp->vMol[Ip*usp->N+(pmu->mui[i])] = pmu->IC_m[i];
    }
    usp->vpH[Ip][0] = pmu->pH;
    usp->vpH[Ip][1] = pmu->Eh;
    usp->vpH[Ip][2] = pmu->IC;
    if(usp->Ls/*TProfil::pm->mup->Pg*/)
       for( i=0; i<pmu->Ls; i++)
          usp->vFug[Ip*usp->Ls+(pmu->muj[i])] =
                pmu->Y_la[i];

    usp->vT[Ip]= pmu->TCc;
    usp->vP[Ip]= pmu->Pc;
    usp->vV[Ip]= pmu->VXc;
  }

 // added for copy of input data
  if( usp->PsGen[1]== S_ON )
    for( i=0; i<usp->L; i++)
#ifndef IPMGEMPLUGIN
        if( TProfil::pm->tpp->S )
            usp->vS[i] = TProfil::pm->tpp->S[i];
#else
            if( tpp_S )   // not initalazed
                usp->vS[i] = tpp_S[i];// pmu->tpp_S
#endif

  if( usp->PsGen[5]== S_ON )
    for( i=0; i<usp->L; i++)
    {
#ifndef IPMGEMPLUGIN
        double xx = syu->Vuns[i];
                xx += TProfil::pm->tpp->Vm[i];
#else
       double xx = pmu->Vuns[i];
               xx += pmu->tpp_Vm[i];
#endif
      usp->vmV[Ip*usp->L+i]= xx;
    }

  if( usp->PsGen[2]== S_ON ) // 16/02/2007
  {
    for( i=0; i<usp->N; i++)
      usp->vB[Ip*usp->N+i] = 0.;
    for( i=0; i<pmu->N; i++)
    {    double xx = pmu->B[i];
         usp->vB[Ip*usp->N+(pmu->mui[i])] = xx;
    }
  }

//  if(usp->PsGen[6]== S_ON && usp->Ls )   // new by DK
//    for( i=0; i<usp->Ls; i++)
//    {    usp->vNidP[i] = ?????;
//    }

 
 }
}

//===============================================
// From Kostya  (not for changed)

/*  select initial values x0 */
void TUnSpace::UNIFORM0( int reg )
// reg: 0 - uniform, 1 - normal
// usp->OVR - ravnomernye sluchainye chusla
// usp->OVN - normal sluchainye chusla
{
  int i,j,k;
  double R;

  if(!reg)
    k=usp->nGR;
  else
    k=usp->nGN;

  for(i=0;i<k;i++)
  {   j=rand();
      R=ceil(24359738368.*j/RAND_MAX + 10000000000.);
      if(!fmod(R,2))
         R=R+1.;
       if(!reg)
         usp->OVR[i+1]=R;
       else
         usp->OVN[i+1]=R;
    }
}

void TUnSpace::BELOV(int QT, int NG, float *OVB)
// QT - kol-vo points
// NG - kol-vo groups
// OVB - chisla po zakonu Belova
{
   int QT2,i,j,k,D,B,W,MC,m2,*LI,*LI1;

   QT2=(QT-1)/2;

   LI = new int[QT2+1];
   LI1 = new int[QT2+1];

   W=floor(float(NG)/float(QT2));
   if(W)
   {  k=1;
      for(i=1;i<=QT2;i++)
         for(j=0;j<=W-1;j++)
            OVB[i+j*QT2]=i;
   }
   else
      k=0;
   W*=QT2;
   if(W>=NG)
    {
        delete[] LI;
        delete[] LI1;
        return;
    }

  OVB[W+1]=k*QT2+1;
  for(i=1;i<=QT2;i++)
   LI[i]=WL(i,OVB[W+1],QT);
  for(D=W+2;D<=NG;D++)
  { MC=0;
    for(j=1;j<=QT;j++)
    { for(i=1;i<=QT2;i++)
         LI1[i]=WL(i,j,QT)+LI[i];
      B=LI1[1];
      for(i=2;i<=QT2;i++)
        if(LI1[i]<B)
           B=LI1[i];
      m2=0;
      if(B==MC)
       for(i=1;i<=D-1;i++)
         if(OVB[i]==OVB[D])
         { OVB[D]=j;
            m2=1;
            break;
         }
       if(!m2&&B>MC)
         { MC=B;
           OVB[D]=j;
         }
      }
      for(i=1;i<=QT2;i++)
        LI[i]+=WL(i,OVB[D],QT);
   }
  delete[] LI;
  delete[] LI1;
}

int TUnSpace::WL(int i,int j,int QT)
{ int W;
  W=i*j-floor(float(i*j)/float(QT))*QT;  //  ??????????
  if( QT-W < W )
      W=QT-W;
  return(W);
}

// Belov point
float TUnSpace::PNTBEL(int J,int QT,int OV)
// J - number point
// QT - numbers points
// OV -  number of Belov for group
{ float TB;
  TB=fabs(float(J*OV)-floor(float(J*OV)/float(QT))*float(QT));
  TB+=(QT-1)/2;
  while(1)
   if(TB>=QT) TB-=QT;
    else break;
  return(TB);
}

// uniform point
double TUnSpace::ravrand(double *x)
{ double m35=34359738368., m36=68719476736., m37=137438953472.;
  float a=0.,b=1.;
  *x=*x*5.;
  if(*x>=m37) *x=*x-m37;
  if(*x>=m36) *x=*x-m36;
  if(*x>=m35) *x=*x-m35;
 return(*x/m35*(b-a)+a);
}

// normal point
double TUnSpace::norrand(double *x)
{ double R1=0.;
  int j;
  for(j=0;j<101;j++)
    R1+=ravrand(x);
      R1=(R1-101./2.)/pow(101./12.,0.5);
           R1=1./6.*(R1-(-3.0));
           if(R1<0.) R1=0.;
           if(R1>1.) R1=1.;
           return(R1);
/*return(1./9.*(R1-(-4.5)));*/
}

// Next point space uncertainty
void  TUnSpace::NexT(int J )
// J - point of sample
{
  int i,j,k1,k2,k3;
  float R;
  double x, xx;

#ifdef IPMGEMPLUGIN
 DATABR* dBR = TNode::na->pCNode();
 //TC = TNode::na->cTC();
 //P = TNode::na->cP()/bar_to_Pa;
#endif
 
  k1 = k2 = k3 = 1;
  for(i=1; i<=usp->nG; i++)
  { if(!usp->PbD[i-1])
    {  R = PNTBEL( J, usp->Q, usp->OVB[k1]);
      k1++;
      if(usp->Q > 1)
         R/=(usp->Q-1);
    }
   if(usp->PbD[i-1]==1)
     {  x = usp->OVR[k2];
        R = ravrand(&x);
        usp->OVR[k2] = x;
        k2++;
     }
   if(usp->PbD[i-1]==2)
     {  x = usp->OVN[k3];
        R = norrand(&x);
        usp->OVN[k3] = x;
        k3++;
     }
//#ifndef IPMGEMPLUGIN // test for compare
   usp->ncp[J*usp->nG+(i-1)]=R;
//#else
//   R = usp->ncp[J*usp->nG+(i-1)];
//#endif
   if( usp->PsGen[0]== S_ON || usp->PsGen[1]== S_ON || usp->PsGen[5]== S_ON )
    for( j=0; j<usp->L; j++)
     if( ( usp->PsGen[0]== S_ON && usp->NgLg[j]==i) ||
          ( usp->PsGen[1]== S_ON && usp->NgLs[j]==i) ||
        ( usp->PsGen[5]== S_ON && usp->NgLv[j]==i ) )
     {

      if(usp->PsGen[0]== S_ON && usp->NgLg[j]==i)
      {
        xx = 2*usp->IntLg[j][0]*R-usp->IntLg[j][0];
#ifndef IPMGEMPLUGIN
        syu->Guns[j] = xx;
#else
        pmu->Guns[j] = xx;
#endif
      }
/*      if(usp->PsGen[1]== S_ON&& usp->NgLs[j]==i)
      {  xx = usp->Ss[j][0] - usp->IntLs[j][0]+2*usp->IntLs[j][0]*R;
//         usp->R1[ii*17+3]= xx;
      }
*/
      if(usp->PsGen[5]== S_ON&& usp->NgLv[j]==i)
      {
        xx = 2*usp->IntLv[j][0]*R-usp->IntLv[j][0];
#ifndef IPMGEMPLUGIN
        syu->Vuns[j] = xx;
#else
        pmu->Vuns[j] = xx;
#endif
      }
    }
    if( usp->PsGen[2]== S_ON )
     for(j=0; j<usp->N; j++) // 16/02/2007
     {  if( usp->NgNb[j] == i)
        { xx = usp->Bs[j][0] - usp->IntNb[j][0] + 2*usp->IntNb[j][0]*R;
          if( xx < 0 )
            xx = 0.;
#ifndef IPMGEMPLUGIN
          syu->B[j]= xx;
#else
          dBR->bIC[TNode::na->IC_xCH_to_xDB( j )] = xx;
          syu_B[j]= xx;
#endif
        }
     }
     /*for( j=0; j<pmu->N; j++)
     { int jj= pmu->mui[j];
       if( usp->NgNb[jj] == i)
        { xx = usp->Bs[jj][0] - usp->IntNb[jj][0] + 2*usp->IntNb[jj][0]*R;
          if( xx < 0 )
            xx = 0.;
          pmu->B[j]= xx;
        }
    } */

    if( usp->PsGen[3]== S_ON && usp->NgT==i )
     { xx = usp->T[0] - usp->IntT[0] + 2* usp->IntT[0]*R;
       usp->Tc = xx;
       if(usp->Tc<0.)
          usp->Tc=0.;
#ifdef IPMGEMPLUGIN
         dBR->TC = usp->Tc; 
#endif          
     }
   if( usp->PsGen[4]== S_ON && usp->NgP==i )
     {
       xx = usp->P[0]- usp->IntP[0] + 2*usp->IntP[0]*R;
       usp->Pc = xx;
       if(usp->Pc<0.) usp->Pc = 1e-5;
#ifdef IPMGEMPLUGIN
         dBR->P = usp->Pc; 
#endif          
     }
   }
}

//=========================================================================
// Analyze part
//=========================================================================

// calculated pay-off matrix  && analise data
void TUnSpace::analyseArrays( )
{
     int k, Ngr;

     usp->ob=0;

     setPhaseAssemb();   // set up phase assemblages data
     Ngr = sel_parg( usp->sv );

     for( k=0; k<usp->Q; k++ )
        if(  usp->sv[k] == Ngr /*&& !filters( k ) */ )
                // sv[i] < 0  if !filters
                // Ngr > 0  aqlways
            usp->ob++;
     switch( usp->Pa_OF )
     {
       case UNSP_OF_A:
       case UNSP_OF_B:
       default:   Un_criteria( );     // biuld payoff matrix and analyse
                  out_QT( Ngr );  // put data for result arrays
                  break;

     }
}


short TUnSpace::kol_in_sol( int j )
//  return: kol-vo solutions s naborom phases kak v j-solution
{
return usp->PhNum[ abs(usp->sv[j])-1];
}

int TUnSpace::kolgrup()
// return: kol-vo phase groups
{
return usp->nPhA;
}


int TUnSpace::kvant_parag( int type /*short *mas*/ )
{   // sol[i] - kol-vo solutions c usp->quanCx[i][type] naborom phases
    // return:  kol-vo solutions v quantile

  int i,j,I=0,S=0,sum,kvant;
  short maxkol,*sol;

  kvant = (int)(usp->quan_lev*usp->Q);
  if(kvant<1)
       kvant=1;
  sol= new  short[kvant];

  for( i=0; i<kvant; i++)
   sol[i]= kol_in_sol( usp->quanCx[i][type] );
  for( i=0; i<kvant; i++)
  { sum=0;
    for( j=0; j<kvant ; j++)
     if( sol[j] == sol[i] )
       sum++;
    if( sum > S)
      { S=sum; I=i; }
  }
maxkol = sol[I];
delete[] sol;
return(I); //return(maxkol);
}

// indexes into array <mas> po quantile
void TUnSpace::kvant_index(float per,int N,double *mas, short type /**ind*/)
{
  double A0=-1,A1=fabs(mas[0]);
  short i,j,kvant;

  kvant=(int)(per*N);
  if(kvant<1)
      kvant=1;
  for(j=0; j<kvant; j++)
  { for(i=0 ;i<N; i++)
     if(fabs(mas[i])>A0 && fabs(mas[i])<=A1)
     { A1=fabs(mas[i]);
       usp->quanCx[j][type]=i;   //ind[j]=i;
       usp->quanCv[j][type]=mas[i];
     }
    A0=A1;
    for(i=0;i<N;i++)
     if(fabs(mas[i])>A1)
      { A1=fabs(mas[i]); break; }
   }
}


int TUnSpace::sel_parg( short *sv )
{
  int i, j,I,J, maxgr=1;

 if( usp->Pa_Crit == UNSP_CRIT_PA ) // select statistic
 {
  J = 0;
  for( j=0; j<usp->Q; j++)
     if(abs(sv[j])==1)
          J++;
  for(i=1; i<kolgrup(); i++ )
  { I=0;
    for( j=0; j<usp->Q; j++ )
     if(abs(sv[j])==i+1)
        I++;
    if(I>J)
    { maxgr = i+1 ; J=I; }
  }
 }
 else
  if(usp->Pa_Crit <= UNSP_CRIT_HOME_QAN)
  { maxgr = kvant_parag( usp->Pa_Crit - UNSP_CRIT_LAPL_QAN );
    maxgr++;
  }
return(maxgr);

}

int TUnSpace::filters( int k )
// return: 1-not, 0-yes (propuskaet)
{
  int Filtr=0,i;
  if( usp->Pa_f_pH == S_ON &&  ( usp->vpH[k][0] < usp->pH_lo ||
                         usp->vpH[k][0] > usp->pH_up ))
    Filtr=1;
  if( usp->Pa_f_Eh == S_ON &&  ( usp->vpH[k][1] < usp->Eh_lo ||
                         usp->vpH[k][1] > usp->Eh_up ))
    Filtr=1;
  if( usp->Pa_f_IC == S_ON &&  ( usp->vpH[k][2] < usp->IC_lo ||
                         usp->vpH[k][2] > usp->IC_up ))
    Filtr=1;
  if( usp->Pa_f_fug == S_ON )
     for(i=0; i<usp->Ls; i++)
      if(  usp->fug_up[i] &&
          ( usp->vFug[k*usp->Ls+i] < usp->fug_lo[i] ||
            usp->vFug[k*usp->Ls+i] > usp->fug_up[i] ))
       Filtr=1;
  if( usp->Pa_f_mol == S_ON )
    for(i=0; i<usp->N;i++)
       if( usp->m_t_up[i] &&
        ( usp->vMol[k*usp->N+i] < usp->m_t_lo[i] ||
          usp->vMol[k*usp->N+i] > usp->m_t_up[i]))
       Filtr=1;
  // phases
  if( usp->Pa_f_pha == S_ON )
     for(i=0; i<usp->N; i++)
      if(  usp->PhAndx[(abs(usp->sv[k])-1)*usp->N+i] != usp->f_PhA[i] )
         Filtr=1;

return(Filtr);
}

//=========================================================================
// Calculation part for the payoffs matrix

// Returns an element e(t,q) of the pay-off matrix using function A
//    Eq. 8 on p.18,19 of TM-44-04-01 (with v_qt_j if PvSi == ON or a 
//    simplified function with v_q_j if PvSi == OFF ) 
// i=t, j=q
double TUnSpace::ePO( int t, int q )
{
  double PM,/*rab,*/RG, RGT, ln5551 = 4.0165339;
  int i,k,ii,z,GF=-1,WF=-1,i1,j1;
  short Laq_;
#ifndef IPMGEMPLUGIN
  Laq_ = TProfil::pm->mup->Laq;
#else
  Laq_ = mup_Laq; 
#endif
  
  if(usp->PvSi == S_ON)
    { j1=q; i1=t; i=t; }
  else
    { j1=q; i1=t; i=q; }

  RG = R_CONSTANT;
  RGT = RG*(usp->vT[j1]+273.15);
    
  if( Laq_ )
     WF=0;

#ifndef IPMGEMPLUGIN
  if( TProfil::pm->mup->Pg )
#else
  if( mup_Pg )
#endif
  { 
     if( Laq_ )
         GF=0;
    else GF=1;
   }
   PM=0.;
   ii=0;
   for( z=0; z<usp->Fi; z++)
   {
#ifndef IPMGEMPLUGIN
    for( k=ii; k<ii+(TProfil::pm->mup->Ll[z]); k++)
       if (syu->Dcl[k]!=S_OFF )
#else
     for( k=ii; k<ii+(mup_Ll[z]); k++)
//         if (syu_Dcl[k]!=S_OFF )
#endif
       { if( z==WF )
         { if( k < ( Laq_-1) &&
               usp->vYF[i*usp->Fi] >1e-19 &&
               usp->vY[i*usp->L+( Laq_ )-1] >1e-19 &&
               usp->vY[i*usp->L+ k ] >1e-19 )

             PM += (usp->vG[j1*usp->L+k]/RGT + ln5551 - usp->vY[i*usp->L+(Laq_)-1]/
                    usp->vYF[i*usp->Fi] + log(usp->vY[i*usp->L+k]) - log(usp->vY[i*usp->L+(Laq_)-1]) +
                    1. + usp->vGam[i*usp->L+k]) * usp->vY[i1*usp->L+k];

           if( k == (Laq_-1) &&
               usp->vY[i*usp->L+(Laq_-1)] > 1e-19 &&
               usp->vYF[i*usp->Fi+z] > 1e-19 )

              PM += (usp->vG[j1*usp->L+k]/RGT + log( usp->vY[i*usp->L+k] ) - log( usp->vYF[i*usp->Fi+z] ) -
                  ( usp->vY[i*usp->L+k] / usp->vYF[i*usp->Fi+z]) - ( usp->vYF[i*usp->Fi+z] / usp->vY[i*usp->L+k]) +
                  2.+ usp->vGam[i*usp->L+k]) * usp->vY[i1*usp->L+k];

          continue;
         }
        if ( z==GF)
        { if( usp->vY[i*usp->L+k] > 1e-19 &&  usp->vYF[i*usp->Fi+z] >1e-19 )
            PM += ( usp->vG[j1*usp->L+k] / RGT  +
                    log(usp->vY[i*usp->L+k]) - log(usp->vYF[i*usp->Fi+z]) +
                    usp->vGam[i*usp->L+k] + log(usp->vP[i]) ) * (usp->vY[i1*usp->L+k]);
         continue;
        }
        if( z>WF && z>GF &&
            usp->vY[i*usp->L+k] > 1e-19 &&
            usp->vYF[i*usp->Fi+z] > 1e-19 )
           PM += ( usp->vG[j1*usp->L+k] / RGT  + log(usp->vY[i*usp->L+k])
           - log(usp->vYF[i*usp->Fi+z])  + usp->vGam[i*usp->L+k] ) * (usp->vY[i1*usp->L+k] );
       }
#ifndef IPMGEMPLUGIN
    ii += TProfil::pm->mup->Ll[z];
#else
    ii += mup_Ll[z];
#endif
     }

//     for( z=0; z<usp->N; z++) // 16/02/2007
//       PM -= syu->B[z] * (usp->vU[i1*usp->N+z]);
       for( z=0; z<pmu->N; z++) //  B[z] is zero behind MULTI
         PM -= pmu->B[z] * (usp->vU[i1*usp->N+(pmu->mui[z])]); 
         // To check that vB[] is properly used here except pmu->B[] !!!! 

   return(PM);
}

// calculate the element e(t,q) of the pay-off matrix B
// Eq 9 on p. 20 of TM 44-04-01 (with v_qt,j if PvSi == ON and v_q_j otherwise)
// i=t, j=q
double TUnSpace::ePO1( int t, int q )
{
	  double /*rab,*/RG, RGT, ln5551 = 4.0165339;
	  int i,k,ii,z,zi,GF=-1,WF=-1,i1,j1, L1F;
  short Laq_;
#ifndef IPMGEMPLUGIN
  Laq_ = TProfil::pm->mup->Laq;
#else
  Laq_ = mup_Laq; 
#endif

  if(usp->PvSi == S_ON)
    { j1=q; i1=t; i=t; }
  else
    { j1=q; i1=t; i=q; }

  RG = R_CONSTANT;
  RGT = RG*(usp->vT[j1]+273.15);
    
  if( Laq_ )
     WF=0;

#ifndef IPMGEMPLUGIN
  if( TProfil::pm->mup->Pg )
#else
  if( mup_Pg )
#endif
  { if(!Laq_)
         GF=0;
    else GF=1;
   }

  double PM=0.;   // Sum for primal potentials for species
  double DM=0.,    // Sum for dual potential for species 
         DDM=0.,  // Sum for differences between PM and DM for a single species 
         Dif, Ps, Ds, Ddm1;   // Difference, primal, dual potential for a species 
  ii=0;
  for( z=0; z<usp->Fi; z++)
  {
#ifndef IPMGEMPLUGIN
L1F = TProfil::pm->mup->Ll[z];
  for( k=ii; k<ii+(TProfil::pm->mup->Ll[z]); k++)
  {
	  if (syu->Dcl[k]==S_OFF )
		 continue; 
#else
L1F = mup_Ll[z];
  for( k=ii; k<ii+(mup_Ll[z]); k++)
//    if (syu_Dcl[k]!=S_OFF )
  {	  
#endif
      Ps = Ds = 0.;
      if( z==WF )
      {
         if( k < (Laq_-1) &&  usp->vYF[i*usp->Fi] >1e-19 &&
               usp->vY[i*usp->L+(Laq_)-1] >1e-19 && usp->vY[i*usp->L+ k ] >1e-19  )
         {
        	 Ps = usp->vG[j1*usp->L+k]/RGT + ln5551 - usp->vY[i*usp->L+(Laq_-1)]/
                    usp->vYF[i*usp->Fi] + log(usp->vY[i*usp->L+k]) - log(usp->vY[i*usp->L+(Laq_-1)]) +
                   1. + usp->vGam[i*usp->L+k];
             for( zi=0; zi<usp->N; zi++)
                Ds += usp->A[k*usp->N+zi]* (usp->vU[i1*usp->N+zi]);
         }
         if( k == (Laq_-1) && usp->vY[i*usp->L+(Laq_-1)] > 1e-19 &&
               usp->vYF[i*usp->Fi+z] > 1e-19 )
         {
           	 Ps = usp->vG[j1*usp->L+k]/RGT + log( usp->vY[i*usp->L+k] ) - log( usp->vYF[i*usp->Fi+z] ) -
                 ( usp->vY[i*usp->L+k] / usp->vYF[i*usp->Fi+z]) - ( usp->vYF[i*usp->Fi+z] / usp->vY[i*usp->L+k]) +
                 2.+ usp->vGam[i*usp->L+k];
           	for( zi=0; zi<usp->N; zi++)
           	     Ds += usp->A[k*usp->N+zi]* (usp->vU[i1*usp->N+zi]);
         }
         goto SPECIES;
        }
        if ( z==GF)
        { 
        	if( (usp->vY[i*usp->L+k] > 1e-19) && usp->vYF[i*usp->Fi+z] >1e-19 )
        	{	
              Ps = ( usp->vG[j1*usp->L+k] / RGT + log(usp->vY[i*usp->L+k]) - log(usp->vYF[i*usp->Fi+z]) +
                    usp->vGam[i*usp->L+k] + log(usp->vP[i]) );
              for( zi=0; zi<usp->N; zi++)
                 Ds += usp->A[k*usp->N+zi]* (usp->vU[i1*usp->N+zi]);
        	}
        	goto SPECIES;
        }
        if( z>WF && z>GF && L1F > 1 )
        {	 // Multicomponent condensed solution phase
        	if (usp->vY[i*usp->L+k] > 1e-19 &&  usp->vYF[i*usp->Fi+z] > 1e-19 )     
            {  
        	   Ps = ( usp->vG[j1*usp->L+k] / RGT + log(usp->vY[i*usp->L+k])
                   - log(usp->vYF[i*usp->Fi+z]) + usp->vGam[i*usp->L+k] );
        	   for( zi=0; zi<usp->N; zi++)
        	                   Ds += usp->A[k*usp->N+zi]* (usp->vU[i1*usp->N+zi]);
            }
            goto SPECIES; 
        }
        if( z>WF && z>GF && L1F == 1 &&  usp->vYF[i*usp->Fi+z] > 1e-19 )  // Single-component condensed phase
        {	
        	Ps = usp->vG[j1*usp->L+k] / RGT;        	
        	for( zi=0; zi<usp->N; zi++)
               Ds += usp->A[k*usp->N+zi]* (usp->vU[i1*usp->N+zi]);
        }        
SPECIES: Dif = Ps - Ds; 
//         if( Dif < 10. )
        	 DDM += Dif; 
//         DDM += Dif;
         PM += Ps; DM += Ds;         
      }  // k loop
#ifndef IPMGEMPLUGIN
    ii += TProfil::pm->mup->Ll[z];
#else
    ii += mup_Ll[z];
#endif
    } // z loop

    Ddm1 = PM - DM;
    return(DDM);
}

  // calculate the element e(t,q) of the pay-off matrix B
  // Eq 10 on p. 20 of TM 44-04-01 (with v_qt,j if PvSi == ON and v_q_j otherwise)
  // i=t, j=q
  double TUnSpace::ePO2( int t, int q )
  {
  	  double RG, RGT, ln5551 = 4.0165339;
  	  int i,k,ii,z,GF=-1,WF=-1,i1,j1, L1F;
    short Laq_;
  #ifndef IPMGEMPLUGIN
    Laq_ = TProfil::pm->mup->Laq;
  #else
    Laq_ = mup_Laq; 
  #endif

    if(usp->PvSi == S_ON)
      { j1=q; i1=t; i=t; }
    else
      { j1=q; i1=t; i=q; }

    RG = R_CONSTANT;
    RGT = RG*(usp->vT[j1]+273.15);
      
    if( Laq_ )
       WF=0;

  #ifndef IPMGEMPLUGIN
    if( TProfil::pm->mup->Pg )
  #else
    if( mup_Pg )
  #endif
    { if(!Laq_)
           GF=0;
      else GF=1;
     }

    double PM=0.;   // Sum for primal potentials for species
    double DM=0.,    // Sum for dual potential for species 
           DDM=0.,  // Sum for differences between PM and DM for a single species 
           Dif, Ps, Ds, Ddm1;   // Difference, primal, dual potential for a species 
    ii=0;
    for( z=0; z<usp->Fi; z++)
    {
  #ifndef IPMGEMPLUGIN
  L1F = TProfil::pm->mup->Ll[z];
    for( k=ii; k<ii+(TProfil::pm->mup->Ll[z]); k++)
    {
  	  if (syu->Dcl[k]==S_OFF )
  		 continue; 
  #else
  L1F = mup_Ll[z];
    for( k=ii; k<ii+(mup_Ll[z]); k++)
  //    if (syu_Dcl[k]!=S_OFF )
    {	  
  #endif
        Ps = Ds = 0.;
        if( z==WF )
        {
           if( k < (Laq_-1) &&  usp->vYF[i*usp->Fi] >1e-19 && usp->vY[i1*usp->L+ k ] >1e-19 &&
                 usp->vY[i*usp->L+(Laq_)-1] >1e-19 && usp->vY[i*usp->L+ k ] >1e-19  )
           {
          	 Ps = usp->vG[j1*usp->L+k]/RGT + ln5551 - usp->vY[i*usp->L+(Laq_-1)]/
                      usp->vYF[i*usp->Fi] + log(usp->vY[i*usp->L+k]) - log(usp->vY[i*usp->L+(Laq_-1)]) +
                     1. + usp->vGam[i*usp->L+k];
         	 Ds = usp->vG[i1*usp->L+k]/RGT + ln5551 - usp->vY[i1*usp->L+(Laq_-1)]/
                      usp->vYF[i1*usp->Fi] + log(usp->vY[i1*usp->L+k]) - log(usp->vY[i1*usp->L+(Laq_-1)]) +
                     1. + usp->vGam[i1*usp->L+k];
           }
           if( k == (Laq_-1) && usp->vY[i*usp->L+(Laq_-1)] > 1e-19 && usp->vY[i1*usp->L+(Laq_-1)] > 1e-19 &&
                 usp->vYF[i*usp->Fi+z] > 1e-19 )
           {
             	 Ps = usp->vG[j1*usp->L+k]/RGT + log( usp->vY[i*usp->L+k] ) - log( usp->vYF[i*usp->Fi+z] ) -
                   ( usp->vY[i*usp->L+k] / usp->vYF[i*usp->Fi+z]) - ( usp->vYF[i*usp->Fi+z] / usp->vY[i*usp->L+k]) +
                   2.+ usp->vGam[i*usp->L+k];
            	 Ds = usp->vG[i1*usp->L+k]/RGT + log( usp->vY[i1*usp->L+k] ) - log( usp->vYF[i1*usp->Fi+z] ) -
                   ( usp->vY[i1*usp->L+k] / usp->vYF[i1*usp->Fi+z]) - ( usp->vYF[i1*usp->Fi+z] / usp->vY[i1*usp->L+k]) +
                   2.+ usp->vGam[i1*usp->L+k];
           }
           goto SPECIES;
          }
          if ( z==GF)
          { 
          	if( (usp->vY[i*usp->L+k] > 1e-19) && usp->vYF[i*usp->Fi+z] >1e-19 && 
          			(usp->vY[i1*usp->L+k] > 1e-19) && usp->vYF[i1*usp->Fi+z] >1e-19 )
          	{	
                Ps = ( usp->vG[j1*usp->L+k] / RGT + log(usp->vY[i*usp->L+k]) - log(usp->vYF[i*usp->Fi+z]) +
                      usp->vGam[i*usp->L+k] + log(usp->vP[i]) );
                Ds = ( usp->vG[i1*usp->L+k] / RGT + log(usp->vY[i1*usp->L+k]) - log(usp->vYF[i1*usp->Fi+z]) +
                      usp->vGam[i1*usp->L+k] + log(usp->vP[i1]) );
          	}
          	goto SPECIES;
          }
          if( z>WF && z>GF && L1F > 1 )
          {	 // Multicomponent condensed solution phase
          	if (usp->vY[i*usp->L+k] > 1e-19 &&  usp->vYF[i*usp->Fi+z] > 1e-19 &&
          			usp->vY[i1*usp->L+k] > 1e-19 &&  usp->vYF[i1*usp->Fi+z] > 1e-19 )     
            {  
          	   Ps = ( usp->vG[j1*usp->L+k] / RGT + log(usp->vY[i*usp->L+k])
                     - log(usp->vYF[i*usp->Fi+z]) + usp->vGam[i*usp->L+k] );
          	   Ds = ( usp->vG[i1*usp->L+k] / RGT + log(usp->vY[i1*usp->L+k])
          	                      - log(usp->vYF[i1*usp->Fi+z]) + usp->vGam[i1*usp->L+k] );
            }
            goto SPECIES; 
          }
          if( z>WF && z>GF && L1F == 1 &&  usp->vYF[i*usp->Fi+z] > 1e-19 
        		  &&  usp->vYF[i1*usp->Fi+z] > 1e-19 )  // Single-component condensed phase
          {	
          	Ps = usp->vG[j1*usp->L+k] / RGT;  
          	Ds = usp->vG[i1*usp->L+k] / RGT;
          }        
  SPECIES: Dif = Ps - Ds; 
  //         if( Dif < 10. )
          	 DDM += Dif; 
  //         DDM += Dif;
           PM += Ps; DM += Ds;         
        }  // k loop
  #ifndef IPMGEMPLUGIN
      ii += TProfil::pm->mup->Ll[z];
  #else
      ii += mup_Ll[z];
  #endif
      } // z loop
      Ddm1 = PM - DM;
      return(DDM);
}  
 
// calculate the element e(t,q) of the pay-off matrix  using
// Eq 11 on p. 21 of TM 44-04-01 
// i=t, j=q    
double TUnSpace::ePO3( int i,int j )
{  
   double DM=0.,PM=0.;
   int k,z;
   for(k=0;k<usp->L;k++)
   {
	  if(usp->PvSi == S_ON )   // experimental 11.06.08 DK  
	     if(!( usp->vY[i*usp->L+k] > 1e-19 && usp->vY[j*usp->L+k] > 1e-19 ) )
		    continue;
	  DM = 0.;
	  for(z=0;z<usp->N;z++)
         DM += usp->A[k*usp->N+z]*(*(usp->vU+j*usp->N+z))
             - usp->A[k*usp->N+z]*(*(usp->vU+i*usp->N+z));
	  PM += DM; 
   }
/*    
    for(k=0;k<usp->L;k++)
     for(z=0;z<usp->N;z++)
       R+=usp->A[k*usp->N+z]*(*(usp->vU+j*usp->N+z));
    PM+=R;
    R=0;
    for(k=0;k<usp->L;k++)
     for(z=0;z<usp->N;z++)
       R+=usp->A[k*usp->N+z]*(*(usp->vU+i*usp->N+z));
     PM-=R;
*/
 return(PM);
}

// calculate the element e(t,q) of the pay-off matrix  using
// Eq 12 on p. 21 of TM 44-04-01 
// i=t, j=q 
double TUnSpace::ePO4( int i,int j )
{  double PM=0.;
   int z;
   for( z=0; z<usp->N; z++)
        PM+= (usp->vU[j*usp->N+z]) - (usp->vU[i*usp->N+z]);
 return(PM);
}

// calculate the element e(t,q) of the pay-off matrix  using
// Eq 13 on p. 21 of TM 44-04-01 
// i=t, j=q 
double TUnSpace::ePO5( int i,int j )
{  double PM=0.;
   int z;
//   for( z=0; z<usp->N; z++) //
//     PM += syu->B[z]*(usp->vU[j*usp->N+z]) -
//           syu->B[z]*(usp->vU[i*usp->N+z]);
   for( z=0; z<pmu->N; z++) //  B[z] is zero behind MULTI
      PM += pmu->B[z]*(usp->vU[j*usp->N+(pmu->mui[z])]) -
            pmu->B[z]*(usp->vU[i*usp->N+(pmu->mui[z])]);
 return(PM);
}

double TUnSpace::g_uw( int j, double *U, float *A,
       double lgGAM,double x,double xjw, double Xw)
// dual solution water
{ double R,gk;
  int i;
  R=0.;
  for( i=0; i<usp->N; i++)
    R+= A[j*usp->N+i]*( U[i] * R_CONSTANT * 298.15);
  gk = R + R_CONSTANT * 298.15 *
       (-lgGAM - log(55.51) - log(x/Xw) + log(xjw/Xw) + xjw/Xw-1 );
return(gk);
}



double TUnSpace::g_u( int j, double *U,float *A,
                      double lgGAM,double x,double Xw)
// dual solution others
{ double R,gk;
  int i;
  R=0.;
  for(i=0;i<usp->N;i++)
   R+= A[j*usp->N+i] * (U[i] * R_CONSTANT * 298.15);
  gk= R + R_CONSTANT * 298.15 * (-lgGAM + log(x/Xw));
return(gk);
}


// calculate pay-off matrix with different functions
void TUnSpace::Un_criteria()
{
   short  t, q, jj;
   double R;
   double Kr=0.;

   for( t=0; t<usp->Q; t++ )
   {
#ifndef IPMGEMPLUGIN
   pVisor->Message( window(), GetName(),
             "Generation of pay off matrix\n"
                 "Please, wait...", t, usp->Q);
#endif
  


    for( q=0; q<usp->Q; q++ )
    {
      usp->t = t;
      usp->q = q;
      switch( usp->Pa_OF )
     {
       case UNSP_OF_A:  R = ePO( t,q ); break;
       case UNSP_OF_B:  R = ePO1( t,q ); break;
       case UNSP_OF_C:  R = ePO2( t,q ); break;
       case UNSP_OF_D:  R = ePO3( t,q ); break;
       case UNSP_OF_E:  R = ePO4( t,q ); break;
       case UNSP_OF_F:  R = ePO5( t,q ); break;
     }
        if(!q)
         {
            // save  pay-off matrix
            if( usp->PvPOM == S_ON )
            {
             usp->pmr = usp->POM + t*usp->Q ;
#ifndef IPMGEMPLUGIN
             aObj[ o_unpmr].SetPtr( usp->pmr );
#endif
             usp->POM[ t*usp->Q+q ] = R;
            }
            if( usp->PvPOR == S_ON )
             usp->POR[ q ] = R;
           usp->Zmin[t] = usp->Zmax[t]= R;
         if( usp->Pa_Zcp == S_OFF )
           usp->Zcp[t] = R;
         else 
           usp->Zcp[t] = fabs(R);
         usp->ZmaxAbs[t] = fabs(R);
           
         }
        else
        {
            // save  pay-off matrix
            if( usp->PvPOM == S_ON )
             usp->POM[ t*usp->Q+q ] = R;
            if( usp->PvPOR == S_ON )
             usp->POR[ q ] = R;

             if( R >usp->Zmax[t] )
                 usp->Zmax[t] = R;
             if( R < usp->Zmin[t] )
                 usp->Zmin[t] = R;
             if( fabs(R) > usp->ZmaxAbs[t] )
                 usp->ZmaxAbs[t] = fabs(R);
         if( usp->Pa_Zcp == S_OFF )
             usp->Zcp[t] += R;
          else
             usp->Zcp[t] += fabs(R);
         }
        Kr += fabs(R);
        if(!t )
         usp->Prob[q] = fabs(R);
        else
         usp->Prob[q] += fabs(R);
    }

    if( usp->Q)        // computing mean
      usp->Zcp[t] /= (double)usp->Q;
   }

// posibility for Homenuk

   for( q=0; q<usp->Q;q++)
       usp->Prob[q] /= Kr;

//  Wald criteria
   R = usp->ZmaxAbs[0];
   jj=0;
   for( t=1; t<usp->Q; t++ )
    if( usp->ZmaxAbs[t] < R )
    { jj=t;
      R= usp->ZmaxAbs[t];
    }
   usp->Wald = jj;
   usp->CrW = R;
   usp->nw = kol_in_sol(jj );
   kvant_index( usp->quan_lev, usp->Q, usp->ZmaxAbs, 2 );


//  Hurwitz criteria (a=0.5)
    R = fabs( usp->Zmin[0] *.5 + usp->Zmax[0]* .5 );
    jj=0;
    for( t=1; t<usp->Q; t++ )
      if( fabs( usp->Zmin[t]*.5 + usp->Zmax[t]*.5 ) < R )
        { jj = t;
          R = fabs( usp->Zmin[t]*.5 + usp->Zmax[t]*.5 );
        }
   usp->Hurw = jj;
   usp->CrH = R;
   usp->nh = kol_in_sol( jj );
   for( t=0; t<usp->Q; t++)
    usp->Hom[t]= fabs( usp->Zmin[t]*.5 + usp->Zmax[t]*.5 ); // usp->Hom - working array
   kvant_index( usp->quan_lev, usp->Q, usp->Hom, 1 );

//  Laplace criteria
if( usp->nGB > 0)
  R = 1e6;
else
  R= fabs( usp->Zcp[0]);
   jj=0;
  for( t=0; t<usp->Q; t++ )
     if(fabs( usp->Zcp[t] ) < R )
      { jj = t;
        R = fabs(usp->Zcp[t]);
      }
   usp->Lapl = jj;
   usp->CrL = R;
   usp->nl = kol_in_sol( jj );

//   Temporarily commented out DK 10.06.2008
//   for( t=0; t<usp->Q; t++)
//    usp->Zcp[t]= fabs( usp->Zcp[t]);
   kvant_index( usp->quan_lev, usp->Q, usp->Zcp, 0/*usp->quanLap*/ );

// 2 prohod for Homenuk

   for( t=0; t<usp->Q; t++ )
   {
    for( q=0; q<usp->Q; q++ )
    {
       switch( usp->Pa_OF )
       {
         case UNSP_OF_A:  R = ePO( t,q ); break;
         case UNSP_OF_B:  R = ePO1( t,q ); break;
         case UNSP_OF_C:  R = ePO2( t,q ); break;
         case UNSP_OF_D:  R = ePO3( t,q ); break;
         case UNSP_OF_E:  R = ePO4( t,q ); break;
         case UNSP_OF_F:  R = ePO5( t,q ); break;
       }

        if(!q)
          usp->Hom[t] = R*usp->Prob[q];
        else
          usp->Hom[t] += R*usp->Prob[q];
    }
     usp->Hom[t] = fabs( usp->Hom[t] );
   }

//  Homenuk criteria
   R = usp->Hom[0];
   jj=0;
   for( t=1; t<usp->Q; t++ )
     if( usp->Hom[t] < R )
      { jj = t;
        R = usp->Hom[t];
      }
   usp->Homen = jj;
   usp->CrHom = R;
   usp->nHom = kol_in_sol( jj );
   kvant_index( usp->quan_lev, usp->Q, usp->Hom, 3/*usp->quanHom*/ );

}



//=========================================================================

// output table with set phases: numbers and %
void TUnSpace::setPhaseAssemb( )
{ int i, j, fl, np, k, numPH;

// set array usp->sv Indices of PhA ids in sampled variants
// and nPhA Number of different phase assemblages in all sampled variants
  for(i=0; i<usp->Q; i++ )
    usp->sv[i] = -1;

  for( i=0, usp->nPhA=0; i<usp->Q; i++)
   {
     if( usp->sv[i] < 0 )
     { usp->sv[i] = usp->nPhA;
       for( j=i+1; j<usp->Q; j++)
       { fl=0;
         for( k=0; k<usp->Fi; k++ )
          if( (usp->vYF[i*usp->Fi+k] && !usp->vYF[j*usp->Fi+k])  ||
             (!usp->vYF[i*usp->Fi+k] && usp->vYF[j*usp->Fi+k])   )
          { fl=1; break; }
         if(!fl)
           usp->sv[j] = usp->nPhA;
       }
       usp->nPhA++;
     }
   }

// alloc memory for nPhA size
  phase_lists_new();

for(np=0; np<usp->nPhA; np++ )
    for( i=0; i<usp->N; i++ )
       usp->PhAndx[np*usp->N+i] = -1;

// set data to arrays
for(np=0; np<usp->nPhA; np++ )
{
   usp->PhNum[np] = 0;
   for( i=0; i<usp->Q; i++)
     if( np == usp->sv[i])
     {
       if( !usp->PhNum[np] )
       {
         gstring lst = "";
         numPH=0;

         for( k=0; k<usp->Fi; k++ )
          if( usp->vYF[i*usp->Fi+k] > 1e-19 )
          {
            usp->PhAndx[np*usp->N+numPH ] = (short)k;
            numPH++;
#ifndef IPMGEMPLUGIN
            lst += gstring(
               TProfil::pm->mup->SF[k]+MAXSYMB+MAXPHSYMB, 0, MAXPHNAME);
#else
            lst += gstring(
               mup_SF[k]+MAXSYMB, 0, MAXPHNAME);
#endif
            lst.strip();
            lst += ";";
          }
         sprintf( usp->PhAID[np], "GR%2d", np );
         strncpy( usp->PhAlst[np], lst.c_str(), 80 );
       }
       usp->PhNum[np]++;
     }
   usp->PhAfreq[np] = (float)usp->PhNum[np]/usp->Q*100;
  }

// not from 0
   for( k=0; k<usp->Q; k++ )
   {   usp->sv[k]++;
       if( filters( k ) )
          usp->sv[k] *= -1;
   }

}
//====================================================================


