//-------------------------------------------------------------------
// $Id: ipm_simplex.cpp 705 2006-04-28 19:39:01Z gems $
//
// Copyright (C) 1992-2007 K.Chudnenko, I.Karpov, D.Kulik, S.Dmitrieva
//
// Implementation of parts of the Interior Points Method (IPM) module
// for convex programming Gibbs energy minimization, described in:
// (Karpov, Chudnenko, Kulik (1997): American Journal of Science
//  v.297 p. 767-806)
//
// This file is part of a GEM-Selektor (GEMS) v.2.x.x program
// environment for thermodynamic modeling in geochemistry and of the
// standalone GEMIPM2K code (define IPMGEMPLUGIN).
//
// This file may be distributed under the terms of the GEMS-PSI
// QA Licence (GEMSPSI.QAL)
//
// See http://gems.web.psi.ch/ for more information
// E-mail: gems2.support@psi.ch
//-------------------------------------------------------------------
//

#include "m_param.h"

// Calculation of LPP-based automatic initial approximation of the primal vector x
// using the modified simplex method with two-side constraints on x
//
void TMulti::AutoInitialApproximation( )
{
    long int T,Q,*STR=0,*NMB=0;
    long int i,j,k;
    double GZ,EPS,*DN=0,*DU=0,*AA=0,*B1=0;

    try
    {  // Allocation of work arrays

        for( i=0; i<pmp->N; i++) // added SD 15/07/2009
             pmp->U[i] = 0.;

    	pmp->Ec=0;
        Q=pmp->L;
        DN= new double[Q];
        DU= new double[Q+pmp->N];
        B1= new double[pmp->N];
        ErrorIf( !DN || !DU || !B1, "AutoInitialApproximation()", "Memory alloc error." );
        for( i=0; i<pmp->N; i++)
             DU[i+Q] = 0.;
        EPS = TProfil::pm->pa.p.EPS; //  13.10.00  KC  DK
        GZ = 1./EPS;    

        T=0; // Calcuation of all non-zero values in A and G arrays
        for(i=0;i<pmp->L;i++)
            if(fabs(pmp->G[i])>1E-19)
                T++;
        for(j=0;j<pmp->N;j++)
            for(i=0;i<pmp->L;i++)
                if(fabs(*(pmp->A+i*pmp->N+j))>1E-19)
                    T++;
        if( pmp->PLIM ) // Setting constraints on x elements
            Set_DC_limits(  DC_LIM_INIT );

        for(i=0;i<Q;i++)
        {
            DN[i]=pmp->DLL[i];
            DU[i]=/*1e+6;//*/pmp->DUL[i];
        }

        // for(i=Q;i<Q+pmp->N;i++) DU[i]=0.;
        // Allocation of arrays on T
        AA= new double[T];
        STR= new long int[T];
        NMB= new long int[Q+1];
        ErrorIf( !AA || !STR || !NMB, "AutoInitialApproximation()",
            "Memory allocation error #2");
        for( k=0; k<T; k++)
         STR[k] = 0;

        // if( wn[W_EQCALC].status )
        //   aSubMod[MD_EQCALC]->ModUpdate(" SIMPLEX Approximation ");

        // Copying vector b
        for(j=0;j<pmp->N;j++)
            B1[j]=pmp->B[j];
        k=-1;
        for(i=0;i<pmp->L;i++)
        {  // Loading non-zero values
            if(fabs(pmp->G[i])>1E-19)
            {
                k++;
                AA[k]=-pmp->G[i];
                NMB[i]=k+1;
            }
            else NMB[i]=k+2;
            for(j=0;j<pmp->N;j++)
                if(fabs(*(pmp->A+i*pmp->N+j))>1E-19)
                {
                    k++;
                    AA[k]=*(pmp->A+i*pmp->N+j);
                    STR[k]=j+1;
                }
        }
        NMB[Q]=T+1;
        // Calling generic simplex solver
        SolveSimplex(pmp->N,Q,T,GZ,EPS,DN,DU,B1,pmp->U,AA,STR,NMB);

        // unloading simplex solution into a copy of x vector
        for(i=0;i<pmp->L;i++)
            pmp->Y[i]=DU[i];

        // calculating initial quantities of phases
        TotalPhasesAmounts( pmp->Y, pmp->YF, pmp->YFA );

#ifdef Use_qd_real
        if( TProfil::pm->pa.p.PD > 0)
        {
#endif
            pmp->FX = GX( 0.0 ); // calculation of initial G(X) value
            MassBalanceResiduals( pmp->N, pmp->L, pmp->A, pmp->Y, pmp->B, pmp->C );
#ifdef Use_qd_real
        }
        else
        {
            unsigned int old_cw;
            fpu_fix_start(&old_cw);
            pmp->FX = to_double(qdGX( 0.0 )); // calculation of initial G(X) value
            fpu_fix_end(&old_cw);
            qdMassBalanceResiduals( pmp->N, pmp->L, pmp->A, pmp->Y, pmp->B, pmp->C );
        }
#endif

//        	for(long int i=0; i<pmp->N; i++)
//             cout << i << " C " << pmp->C[i] << " B " << pmp->B[i] << endl;


        // Deleting work arrays
        if( DN) delete[]DN;
        if( DU) delete[]DU;
        if( AA) delete[]AA;
        if( B1) delete[]B1;
        if( STR) delete[]STR;
        if( NMB) delete[]NMB;
    }
    catch( TError& xcpt )
    {
        if( DN) delete[]DN;
        if( DU) delete[]DU;
        if( AA) delete[]AA;
        if( B1) delete[]B1;
        if( STR) delete[]STR;
        if( NMB) delete[]NMB;
        Error( xcpt.title.c_str(), xcpt.mess.c_str());
    }

}

// Generic simplex method with two sided constraints (c) K.Chudnenko 1992
//  SPOS function
//
void TMulti::SPOS( double *P, long int STR[],long int NMB[],long int J,long int M,double AA[])
{
    long int I,K;
    K=0;
    for(I=0; I<=M; I++)
    {
        if( I==STR[NMB[J]+K-1])
        {
            *(P+I)=AA[NMB[J]+K-1];
            if( NMB[J]+K+1!=NMB[J+1])
                K++;
        }
        else *(P+I)=0.;
    }
}

// Generic simplex method with two sided constraints (c) K.Chudnenko 1992
//  START function
//
void TMulti::START( long int T,long int *ITER,long int M,long int N,long int NMB[],
           double GZ,double EPS,long int STR[],long int *BASE, double B[],
           double UND[],double UP[],double AA[],double *A, double *Q )
{
    long int I,J;

    for( I=0;I<M;I++)
        UP[N+I]=0.;
    for( J=0;J<N;J++)
    {
        if(fabs(UP[J])<EPS)
            UP[J]=0.;
        else
        {
            UP[J]-=UND[J];
            if( fabs(UP[J])<EPS)
                UP[J]=EPS;
            else if( UP[J]<0.)
                Error("E00IPM: SolveSimplex()", "Inconsistent LP problem (negative UP[J] value(s) in START()) ");
        }
        SPOS(Q, STR, NMB, J, M, AA);
        for( I=0;I<M;I++)
            B[I]-=Q[I+1]*UND[J];
    }
    for( I=0;I<M;I++)
    {
        if( B[I]<0.)
        {
            B[I]=fabs(B[I]);
            for( J=0;J<T;J++)
                if(STR[J]==I)
                    AA[J]=-AA[J];
        }
    }
    *A=0.;
    *ITER=0;
    for( I=0;I<M;I++)
    {
        *A-=GZ*B[I];
        *(A+(I+1)*(M+1))=B[I];
        for( J=0;J<M;J++)
            *(A+(I+1)*(M+1)+J+1)=0.;
        *(A+(I+1)*(M+1)+I+1)=1.;
        BASE[I]=N+I;
        *(A+I+1)=-GZ;
    }
}

// Generic simplex method with two sided constraints (c) K.Chudnenko 1992
//  NEW function
//
void TMulti::NEW(long int *OPT,long int N,long int M,double EPS,double *LEVEL,long int *J0,
                  long int *Z,long int STR[], long int NMB[], double UP[],
                  double AA[], double *A)
{
    long int I,J,J1;
    double MAX,A1;
    double *P;
    P= new double[M+1];
    ErrorIf( !P, "SolveSimplex()", "At NEW: memory allocation error ");
    J1=*J0;
    MAX=0.;
    for( J=J1+1;J<=N;J++)
    {
        SPOS( P, STR, NMB, J-1, M, AA);
        A1=-P[0];
        for( I=1;I<=M;I++)
            A1+=P[I]*(*(A+I));
        if(fabs(A1)>MAX)
        {
            if(UP[J-1]>=-EPS && A1<-EPS)
            {
                *Z=1;
                goto MK3;
            }
            else if(UP[J-1]<-EPS && A1>EPS)
            {
                *Z=0;
                goto MK3;
            }
            else continue;
MK3:
            MAX=fabs(A1);
            *J0=J;
            if( MAX>=*LEVEL)
                goto MK4;
        }
    }

    for( J=1;J<J1;J++)
    {
        SPOS(P, STR, NMB, J-1, M, AA);
        A1=-P[0];
        for( I=1;I<=M;I++)
            A1+=P[I]*(*(A+I));
        if(fabs(A1)>MAX)
        {
            if(UP[J-1]>=-EPS && A1<-EPS)
            {
                *Z=1;
                goto MK3A;
            }
            else if(UP[J-1]<-EPS && A1>EPS)
            {
                *Z=0;
                goto MK3A;
            }
            else continue;
MK3A:
            MAX=fabs(A1);
            *J0=J;
            if( MAX>=*LEVEL)
                goto MK4;
        }
    }
    *LEVEL=MAX/2;
    if( *LEVEL<EPS)
        *LEVEL=EPS;
MK4:
    if( MAX<EPS)
        *OPT=1;
    else *OPT=0;
    delete[] P;
}


// Generic simplex method with two sided constraints (c) K.Chudnenko 1992
//  WORK function
//
void TMulti::WORK(double GZ,double EPS,long int *I0, long int *J0,long int *Z,long int *ITER,
                   long int M, long int STR[],long int NMB[],double AA[],
                   long int BASE[],long int *UNO,double UP[],double *A,double Q[])
{
    double MIM,A1;
    long int UN,J,I;
    double *P;
    P=  new double[M+1];
    ErrorIf( !P, "SolveSimplex()", "At WORK: memory allocation error. ");
    *UNO=0;
    *ITER=*ITER+1;
    J=*J0-1;
    SPOS(P, STR, NMB, J, M, AA);
    for( I=0;I<=M;I++)
    {
        Q[I]=0.;
        for( J=1;J<=M;J++)
            Q[I]+=*(A+I*(M+1)+J)*P[J];
    }
    Q[0]-=P[0];
    UN=1;
    MIM=0.;
    for( I=1;I<=M;I++)
    {
        if(*Z==1)
            A1=GZ;
        else A1=-GZ;
        if( (*Z==1 && Q[I]>EPS)||(*Z==0 && Q[I]<-EPS))
        {
            UN=0;
            A1=*(A+I*(M+1))/Q[I];
        }
        else if((*Z==1 && Q[I]<-EPS)||(*Z==0 &&Q[I]>EPS))
        {
            J=BASE[I-1];
            if(fabs(UP[J])>EPS)
            {
                UN=0;
                A1=((*(A+I*(M+1)))-fabs(UP[J]))/Q[I];
            }
        }
        if(I==1||((*Z==1 && A1<MIM)||(*Z==0 && A1>MIM)))
        {
            MIM=A1;
            *I0=I;
        }
    }
    if( UN==1 && fabs(UP[*J0-1])<EPS)
    {
        *UNO=1;
        delete[] P;
        return;
    }
    if((fabs(MIM)<fabs(UP[*J0-1]))||(fabs(UP[*J0-1])<EPS))
    {
        J=BASE[*I0-1];
        if((*Z==1&&Q[*I0]>0.)||(*Z==0&&Q[*I0]<0.))
            UP[J]=fabs(UP[J]);
        else UP[J]=-fabs(UP[J]);
        if(*Z==1)
            *(A+(*I0)*(M+1))=MIM;
        else *(A+(*I0)*(M+1))=MIM+fabs(UP[*J0-1]);

        for( I=0;I<=M;I++)
            if( I!=*I0)
                *(A+I*(M+1))-=Q[I]*MIM;

        BASE[*I0-1]=*J0-1;

        A1=1.E0/Q[*I0];

        for( J=1;J<=M;J++)
            *(A+(*I0)*(M+1)+J)*=A1;
        for( I=0;I<*I0;I++)
            for( J=1;J<=M;J++)
                *(A+I*(M+1)+J)-= Q[I]*(*(A+(*I0)*(M+1)+J));
        for( I=*I0+1;I<=M;I++)
            for( J=1;J<=M;J++)
                *(A+I*(M+1)+J)-=Q[I]*(*(A+(*I0)*(M+1)+J));
    }
    else
    {
        for( I=0;I<=M;I++)
        {
            *(A+I*(M+1))-=UP[*J0-1]*Q[I];
        }
        UP[*J0-1] = -UP[*J0-1]; // Fixed error SD 23/01/2009
    }
    delete[] P;
}

// Generic simplex method with two sided constraints (c) K.Chudnenko 1992
//  FIN function
//
void TMulti::FIN(double EPS,long int M,long int N,long int STR[],long int NMB[],
                  long int BASE[],double UND[],double UP[],double U[],
                  double AA[],double *A,double Q[],long int * /*ITER*/)
{
    long int /* K,*/I,J;
    double *P;
    P=  new double[M+1];
    ErrorIf( !P, "SolveSimplex()", "At FIN: memory allocation error. ");

    for( J=0;J<N;J++)
    {
        if( UP[J]>-EPS )
            UP[J]=0.;
        else UP[J]=fabs(UP[J]);
    }

    for( I=1;I<=M;I++)
    {
        UP[BASE[I-1]]=*(A+I*(M+1));
        Q[I]=0.;
    }
    Q[0]=0.;
    for( J=1;J<=N;J++)
    {
        SPOS( P, STR, NMB, J-1, M, AA);
        UP[J-1]+=UND[J-1];
        for( I=0;I<=M;I++)
            Q[I]+=UP[J-1]*P[I];
    }
    for( I=1;I<=M;I++)
        U[I-1] -= *(A+I);  // was =- *(A+I)

    delete[] P;
}

// Generic simplex method with two sided constraints (c) K.Chudnenko 1992
//  Main function
//
//  M  - number of independent components
//  N  - number of unknowns
//  T  - dimension of a work vector AA[] containing all non-zero
//        values of vector GT[] and A[][] matrix (over lines)
//  GZ - Limiting value of the unknown
//  EPS - precision (convergence) criterion (default 1e-9)
//  UND - vector of lower constraints to unknowns
//  UP - input vector of upper constraints to unknowns;
//        output vector of unknowns (simplex solution) (N+M)
//  B -  M input values of independent components (bulk composition)
//  U -  output vector of the dual solution (M)
//  AA - work array (T)
//  STR - markup vector of values in AA array (T)
//  NMB - indices of values in AA
// returns 0 if OK;
//         1 if inconsistent input constraints;
//        -1 if memory allocation error;
//
void TMulti::SolveSimplex(long int M, long int N, long int T, double GZ, double EPS,
                      double *UND, double *UP, double *B, double *U,
                      double *AA, long int *STR, long int *NMB )
{
    long int IT=200,I0=0,J0=0,Z,UNO,OPT=0,ITER, i;
    double LEVEL;
    long int *BASE=0;
    double *A=0,*Q=0;
    try
    {
        A=  new double[(M+1)*(M+1)];
        Q=  new double[M+1];
        BASE=  new long int[M];
        ErrorIf( !A || !Q || !BASE, "SolveSimplex()", "Memory allocation error ");

        fillValue(A, 0., (M+1)*(M+1) );
        fillValue(Q, 0., (M+1) );
        fillValue(BASE, 0L, (M) );

        LEVEL=GZ;
        START( T, &ITER, M, N, NMB, GZ, EPS, STR, BASE, B,  UND, UP, AA, A, Q );

        for( i=0; i<IT; i++ )   // while(1) fixed  03.11.00
        {
            NEW( &OPT, N, M,EPS, &LEVEL, &J0, &Z, STR, NMB, UP, AA, A);
            if( OPT)
                goto FINISH;  // Converged
            WORK( GZ, EPS, &I0, &J0, &Z, &ITER, M, STR, NMB, AA, BASE, &UNO, UP, A, Q);
            if( UNO)
                goto FINISH; // Solution at boundary of the constraints polyhedron
        }
        if( EPS > 1.0e-6 )
        {
         Error( "E01IPM: SolveSimplex()",
             "LP solution cannot be obtained with sufficient precision" );
        }
FINISH: FIN( EPS, M, N, STR, NMB, BASE, UND, UP, U, AA, A, Q, &ITER);
        delete[] A;
        delete[] Q;
        delete[] BASE;
    }
    catch( TError& xcpt )
    {
        if( A) delete[]A;
        if( Q) delete[]Q;
        if( BASE) delete[]BASE;
        Error( xcpt.title.c_str(), xcpt.mess.c_str());
    }

    // Done
}


//-----------------------------------------------------------------------
// Main call to GEM IPM calculation of equilibrium state in MULTI
// (with internal re-scaling of the system)
//
double TMulti::CalculateEquilibriumState( long int typeMin, long int& NumIterFIA, long int& NumIterIPM )
{
 // const char *key;
  double ScFact=1.;

//#ifndef IPMGEMPLUGIN
//  key = rt[RT_SYSEQ].UnpackKey();
//#else
//  key = "GEMIPM2K";
//#endif

  InitalizeGEM_IPM_Data();

  pmp->t_start = clock();
  pmp->t_end = pmp->t_start;
  pmp->t_elap_sec = 0.0;
  pmp->ITF = pmp->ITG = 0;

//  to_text_file( "MultiDump1.txt" );   // Debugging

if( TProfil::pm->pa.p.DG > 1e-5 )
{
   ScFact = SystemTotalMolesIC();
   ScaleSystemToInternal( ScFact );
}

try{
       switch( pmp->tMin )
       {
       case  A_TV_:
       case  U_SV_:
       case  H_PS_:
       case  _S_PH_:
       case  _S_UV_:
                break;
       case  G_TP_:
       default: GibbsEnergyMinimization();
            break;
        }
  }
  catch( TError& xcpt )
  {

      if( TProfil::pm->pa.p.DG > 1e-5 )
         RescaleSystemFromInternal( ScFact );
//      to_text_file( "MultiDump2.txt" );   // Debugging

      NumIterFIA = pmp->ITF;
      NumIterIPM = pmp->ITG;
      pmp->t_end = clock();
      pmp->t_elap_sec = double(pmp->t_end - pmp->t_start)/double(CLOCKS_PER_SEC);

     Error( xcpt.title, xcpt.mess);
  }

  if( TProfil::pm->pa.p.DG > 1e-5 )
       RescaleSystemFromInternal(  ScFact );

//  to_text_file( "MultiDump3.txt" );   // Debugging

  NumIterFIA = pmp->ITF;
  NumIterIPM = pmp->ITG;
  pmp->t_end = clock();
  pmp->t_elap_sec = double(pmp->t_end - pmp->t_start)/double(CLOCKS_PER_SEC);

  return pmp->t_elap_sec;
}


// Calculate total IC mole amounts in b vector and
// return the scaling factor
double TMulti::SystemTotalMolesIC( )
{
  double ScFact, mass_temp = 0.0;

  for( int i=0; i<pmp->N - pmp->E; i++ ) //?????
         mass_temp +=pmp->B[i];

  pmp->TMols = mass_temp;

  pmp->SMols = TProfil::pm->pa.p.DG;
  ScFact = pmp->SMols/pmp->TMols;

  return ScFact;
}

// Resizes MULTI (GEM IPM work structure) data into internally scaled constant mass
void TMulti::ScaleSystemToInternal(  double ScFact )
{
 long int i, j, k;

  pmp->SizeFactor = ScFact;
  pmp->MBX *= ScFact;
  pmp->FitVar[0] *= ScFact;  // added: bugfix by DK 31.05.2010
  pmp->VXc *= ScFact;
  pmp->HXc *= ScFact;
  pmp->HX_ *= ScFact;
  pmp->FX  *= ScFact;
  pmp->Yw  *= ScFact;  // added 08.06.10 DK

  for( j=0; j<pmp->L; j++ )
  {
    if(	pmp->DUL[j] < 1e6  )
       pmp->DUL[j] *= ScFact;

    // if( pmp->DLL[j] > 0.0  )
       pmp->DLL[j] *= ScFact;

        pmp->Y[j] *= ScFact;
        pmp->X[j] *= ScFact;
        pmp->XY[j] *= ScFact;
        pmp->XU[j] *= ScFact;
  }

  for( i=0; i<pmp->N; i++ )
  {
        pmp->B[i] *= ScFact;
  //      pmp->C[i] *= ScFact;
        pmp->BFC[i] *= ScFact;
  }

  for( k=0; k<pmp->FI; k++ )
  {
    pmp->XFs[k] *= ScFact;
    pmp->XF[k] *= ScFact;
    pmp->YF[k] *= ScFact;
    pmp->FVOL[k] *= ScFact;
    pmp->FWGT[k] *= ScFact;
  }

  for( k=0; k<pmp->FIs; k++ )
  {
      pmp->XFA[k] *= ScFact;
      pmp->YFA[k] *= ScFact;

      if( pmp->PUL )
        if( pmp->PUL[k] < 1e6  )
         pmp->PUL[k] *= ScFact;

      if( pmp->PLL )
      // if( pmp->PLL[k] > 0.0  )
         pmp->PLL[k] *= ScFact;

      for( i=0; i<pmp->N; i++ )
           pmp->BF[k*pmp->N+i] *= ScFact;
  }
 if( pmp->FIat > 0 /*&& pm.Lads > 0*/ && pmp->FIs > 0 )
  {
  for( k=0; k<pmp->FIs; k++ )
     for( j=0; j<MST; j++ )
          {
              pmp->XetaA[k][j] *= ScFact;
              pmp->XetaB[k][j] *= ScFact;
              pmp->XetaD[k][j] *= ScFact;
              pmp->XFTS[k][j] *= ScFact;
          }
}

 pmp->SizeFactor = ScFact;
}

// Re-scaling the internal constant-mass MULTI system definition
// back to real size
void TMulti::RescaleSystemFromInternal(  double ScFact )
{
 long int i, j, k;

  pmp->MBX /= ScFact;
  pmp->FitVar[0] /= ScFact;  // added: bugfix by DK 31.05.2010
  pmp->VXc /= ScFact;
  pmp->HXc /= ScFact;
  pmp->HX_ /= ScFact;
  pmp->FX  /= ScFact;
  // added SV
  //pmp->YFk /= ScFact;
  pmp->Yw /= ScFact;  // added 08.06.10 DK

  for( j=0; j<pmp->L; j++ )
  {
    if(	pmp->DUL[j] < 1e6  )
       pmp->DUL[j] /= ScFact;

    // if( pmp->DLL[j] > 0.0  )
       pmp->DLL[j] /= ScFact;

        pmp->Y[j] /= ScFact;
        pmp->X[j] /= ScFact;
        pmp->XY[j] /= ScFact;
        pmp->XU[j] /= ScFact;
    //    pmp->MU[j] /= ScFact;
    //    pmp->W[j] /= ScFact;
  }

  for( i=0; i<pmp->N; i++ )
  {
        pmp->B[i] /= ScFact;
        pmp->C[i] /= ScFact;
        pmp->BFC[i] /= ScFact;
  }

  for( k=0; k<pmp->FI; k++ )
  {
    pmp->XFs[k] /= ScFact;
    pmp->XF[k] /= ScFact;
    pmp->YF[k] /= ScFact;
    pmp->FVOL[k] /= ScFact;
    pmp->FWGT[k] /= ScFact;
  }

  for( k=0; k<pmp->FIs; k++ )
  {
      pmp->XFA[k] /= ScFact;
      pmp->YFA[k] /= ScFact;

      if( pmp->PUL )
        if( pmp->PUL[k] < 1e6  )
         pmp->PUL[k] /= ScFact;

      if( pmp->PLL )
      // if( pmp->PLL[k] > 0.0  )
         pmp->PLL[k] /= ScFact;

      for( i=0; i<pmp->N; i++ )
           pmp->BF[k*pmp->N+i] /= ScFact;
  }

  if( pmp->FIat > 0 /*&& pm.Lads > 0*/ && pmp->FIs > 0 )
   for( k=0; k<pmp->FIs; k++ )
          for( j=0; j<MST; j++ )
          {
              pmp->XetaA[k][j] /= ScFact;
              pmp->XetaB[k][j] /= ScFact;
              pmp->XetaD[k][j] /= ScFact;
              pmp->XFTS[k][j]  /= ScFact;
          }

  pmp->SizeFactor = 1.;   // using in TNode class
}

//========================================================================================
// Multi initialization part 10/05/2010

// Before Calculations
//Calculation by IPM (preparing for calculation, unpacking data)
// parameter "key" contains SysEq record key
// In IPM
void TMulti::InitalizeGEM_IPM_Data( ) // Reset internal data formerly MultiInit()
{

   MultiConstInit();

#ifndef IPMGEMPLUGIN
   // for GEMIPM unpackDataBr( bool uPrimalSol );
   // to define quantities

   //   MultiKeyInit( key ); //into PMtest

   if( pmp->pBAL < 2  )
   {
     // Allocating list of phases currently present in non-zero quantities
     MultiSystemInit( );
   }

   // Allocating list of phases currently present in non-zero quantities
     if( !pmp->SFs )
        pmp->SFs = (char (*)[MAXPHNAME+MAXSYMB])aObj[ o_wd_sfs].Alloc(
                    pmp->FI, 1, MAXPHNAME+MAXSYMB );

   // no old solution => must be simplex
      if( pmp->pESU == 0 )
           pmp->pNP = 0;

      TProfil::pm->CheckMtparam(); //load tpp structure

//      cout << "Init pmp->pTPD " << pmp->pTPD << endl;
//      cout << "Init pmp->P " << pmp->P << endl;
//      cout << "Init pmp->T " << pmp->T << endl;

   if( pmp->pTPD < 2 )
   {
      DC_LoadThermodynamicData(); // Loading thermodynamic data into MULTI structure
     // In future build/rebuild internal lookup arrays
   }

#else

      DC_LoadThermodynamicData(); // Loading thermodynamic data into MULTI structure

#endif

   Alloc_internal();

  // calculate mass of the system
   pmp->MBX = 0.0;
  for(int i=0; i<pmp->N; i++ )
   pmp->MBX += pmp->B[i] * pmp->Awt[i];
   pmp->MBX /= 1000.;

   RescaleToSize( true );  // Added to set default cutoffs/inserts 30.08.2009 DK

   if(  pmp->pNP )
    {  //  Smart Initial Approximation mode
       long int j,k;

#ifndef IPMGEMPLUGIN
       loadData( false );  // unpack SysEq record into MULTI
#endif

     bool AllPhasesPure = true;   // Added by DK on 09.03.2010
     // checking if all phases are pure
     for( k=0; k < pmp->FI; k++ )
     if( pmp->L1[k] > 1 )
        AllPhasesPure = false;

     for( j=0; j< pmp->L; j++ )
         pmp->X[j] = pmp->Y[j];
//       pmp->IC = 0.;  //  Problematic statement!  blocked 13.03.2008 DK
      TotalPhasesAmounts( pmp->X, pmp->XF, pmp->XFA );
      CalculateConcentrations( pmp->X, pmp->XF, pmp->XFA);  // 13.03.2008  DK
       // test multicomponent phases and load data for mixing models
       if( pmp->FIs && AllPhasesPure == false )
       {
           // Load activity coeffs for phases-solutions
         int k, jb, je=0;
         for( k=0; k<pmp->FIs; k++ )
         { // loop on solution phases
            jb = je;
            je += pmp->L1[k];
            // Load activity coeffs for phases-solutions
            for( j=jb; j< je; j++ )
            {
               pmp->lnGmo[j] = pmp->lnGam[j];
               if( fabs( pmp->lnGam[j] ) <= 84. )
      //                pmp->Gamma[j] = exp( pmp->lnGam[j] );
                      pmp->Gamma[j] = PhaseSpecificGamma( j, jb, je, k, 0 );
               else pmp->Gamma[j] = 1.0;
             } // j
          }  // k
       }
    }
}


// Setup/copy flags and thresholds for numeric modules to TMulti structure
// Do it before calculations
void TMulti::MultiConstInit() // from MultiRemake
{
  SPP_SETTING *pa = &TProfil::pm->pa;

  pmp->FI1 = 0;
  pmp->FI1s = 0;
  pmp->FI1a = 0;
  pmp->ITF = 0; pmp->ITG = 0;
  pmp->PD = abs(pa->p.PD);
  pmp->Ec = pmp->K2 = pmp->MK = 0;
  pmp->W1 = 0;
  pmp->is = 0;
  pmp->js = 0;
  pmp->next  = 0;
  pmp->ln5551 = log( H2O_mol_to_kg );             // constant corrected 30.08.2008
  pmp->lowPosNum = Min_phys_amount;               // = 1.66e-24 mol
  pmp->logXw = -16.;
  pmp->logYFk = -9.;
  pmp->DXM = pa->p.DK;

  //  ???????
  pmp->FX = 7777777.;
  pmp->pH = pmp->Eh = pmp->pe = 0.0;
  pmp->YMET = 0;                      // always 0.0 ????
  pmp->PCI = 1.0;
  pmp->FitVar[4] = 1.0;

#ifndef IPMGEMPLUGIN
  pmp->PZ = 0; // IPM default
//  pmp->FitVar[0] = pa->aqPar[0]; // setting T,P dependent b_gamma parameters
  pmp->pH = pmp->Eh = pmp->pe = 0.0;
#else
  pmp->PZ = pa->p.DW;  // in IPM
//  pmp->FitVar[0] = 0.0640000030398369;
#endif

}

//Calculation by IPM (internal step initialization)
void TMulti::GEM_IPM_Init()
{
   int i,j,k;

   for( i=0; i<pmp->N; i++ )
     pmp->Uefd[i] = 0.;

   bool AllPhasesPure = true;   // Added by DK on 09.03.2010
  // checking if all phases are pure
  for( k=0; k < pmp->FI; k++ )
    if( pmp->L1[k] > 1 )
        AllPhasesPure = false;

   if(!pmp->pNP) // Simplex initial approximation to be done
    {
        for( j=0; j<pmp->L; j++ )
        {                           // cleaning work vectors
                pmp->X[j] = pmp->Y[j] = pmp->lnGam[j] = pmp->lnGmo[j] = 0.0;
                pmp->Gamma[j] = 1.0;
                pmp->MU[j] = 0.;
                pmp->XU[j] = 0.;
                pmp->EMU[j] = 0.;
                pmp->NMU[j] = 0.;
                pmp->W[j] = 0.;
                pmp->F[j] = 0.;
                pmp->F0[j] = 0.;
           }
    }

    // recalculating kinetic restrictions for DC amounts
//     if( pmp->pULR && pmp->PLIM )
//          Set_DC_limits(  DC_LIM_INIT );

    if( pmp->FIs && AllPhasesPure == false )   /// line must be tested !pmp->FIs
    {
#ifndef IPMGEMPLUGIN
//????? problematic point
       if( pmp->pIPN <=0 )  // mixing models finalized in any case (AIA or SIA)
       {
             // not done if these models are already present in MULTI !
           pmp->PD = abs(TProfil::pm->pa.p.PD);
           SolModLoad();   // Call point to loading scripts for mixing models
       }
       pmp->pIPN = 1;
#endif
        // Calc Eh, pe, pH,and other stuff
       if( pmp->E && pmp->LO && pmp->pNP )
       {
            CalculateConcentrations( pmp->X, pmp->XF, pmp->XFA);
            IS_EtaCalc();
            if( pmp->Lads )  // Calling this only when sorption models are present
            {
               int k, jb, je=0;
               for( k=0; k<pmp->FIs; k++ )
               { // loop on solution phases
                  jb = je;
                  je += pmp->L1[k];
                  if( pmp->PHC[k] == PH_POLYEL || pmp->PHC[k] == PH_SORPTION )
                  {
                       if( pmp->PHC[0] == PH_AQUEL && pmp->XF[k] > pmp->DSM
                           && (pmp->XFA[0] > pmp->ScMinM && pmp->XF[0] > pmp->XwMinM )) // fixed 30.08.2009 DK
                           GouyChapman( jb, je, k );  // getting PSIs - electrical potentials on surface planes
                  }
               }
           }
       }
       CalculateActivityCoefficients( LINK_TP_MODE);
       // Computing DQF, FugPure and G wherever necessary; Activity coeffs are restored from lnGmo
    }
    else {  // no multi-component phases
        pmp->PD = 0;
        pmp->pNP = 0;
        pmp->pIPN = 1;
    }

    // recalculating kinetic restrictions for DC amounts
     if( pmp->pULR && pmp->PLIM )
          Set_DC_limits(  DC_LIM_INIT );

#ifndef IPMGEMPLUGIN
    // dynamic work arrays - loading initial data
    for( k=0; k<pmp->FI; k++ )
    {
        pmp->XFs[k] = pmp->YF[k];
        pmp->Falps[k] = pmp->Falp[k];
        memcpy( pmp->SFs[k], pmp->SF[k], MAXPHNAME+MAXSYMB );
    }
#endif
}

//===========================================================================================
// Calls to minimization of other system potentials (A, ...)





//--------------------- End of ipm_simplex.cpp ---------------------------
