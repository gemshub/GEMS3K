//-------------------------------------------------------------------
// $Id: ipm_simplex.cpp 705 2006-04-28 19:39:01Z gems $
//
// Copyright (C) 1992-2000 K.Chudnenko, I.Karpov, D.Kulik, S.Dmitrieva
//
// Implementation of parts of the Interior Points Method (IPM) module
// for convex programming Gibbs energy minimization, described in:
// (Karpov, Chudnenko, Kulik (1997): American Journal of Science
//  v.297 p. 767-806)
//
// This file is part of a GEM-Selektor (GEMS) v.2.x.x program
// environment for thermodynamic modeling in geochemistry
// Uses: GEM-Vizor GUI DBMS library, gems/lib/gemvizor.lib
//
// This file may be distributed under the terms of the GEMS-PSI
// QA Licence (GEMSPSI.QAL)
//
// See http://les.web.psi.ch/Software/GEMS-PSI for more information
// E-mail: gems2.support@psi.ch; chud@igc.irk.ru
//-------------------------------------------------------------------
//

#include "m_param.h"

/* Calculation of Simplex initial approximation of prime vector x
*  using modified simplex method with two-side constraints on x
*/
void TMulti::SimplexInitialApproximation( )
{
    int T,Q,*STR=0,*NMB=0;
    register int i,j,k;
    double GZ,EPS,*DN=0,*DU=0,*AA=0,*B1=0;
    /*  char DCC, ICC, PHC; */
    try
    {  // Allocation of work arrays
        pmp->Ec=0;
        Q=pmp->L;
        DN= new double[Q];
        DU= new double[Q+pmp->N];
        B1= new double[pmp->N];
        ErrorIf( !DN || !DU || !B1, "SimplexInitialApproximation()", "Memory alloc erorr." );
        memset(DN, 0, sizeof(double)*Q );
        memset(DU, 0, sizeof(double)*(Q+pmp->N) );
        memset(B1, 0, sizeof(double)*pmp->N );

        // GZ = 1e+9; /* 1e6 */
        // EPS = 1e-9;    pmp->EPS; 1e-6;  Precision ?
        EPS = TProfil::pm->pa.p.EPS; //  added 13.10.00  KVC  DAK
        GZ = 1./EPS;    //  added 13.10.00

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
            DU[i]=pmp->DUL[i];
        }

        /* for(i=Q;i<Q+pmp->N;i++) DU[i]=0.; */
        // Allocation of arrays on T
        AA= new double[T];
        STR= new int[T];
        NMB= new int[Q+1];
        ErrorIf( !AA || !STR || !NMB, "SimplexInitialApproximation", "Memory alloc error #2");
        memset(AA, 0, sizeof(double)*T );
        memset(STR, 0, sizeof(int)*T );
        memset(NMB, 0, sizeof(int)*(Q+1) );

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
        Simplex(pmp->N,Q,T,GZ,EPS,DN,DU,B1,pmp->U,AA,STR,NMB);

        // unloading simplex solution into a copy of x vector
        for(i=0;i<pmp->L;i++)
            pmp->Y[i]=DU[i];

        // calculating initial quantities of phases
        TotalPhases( pmp->Y, pmp->YF, pmp->YFA );

        pmp->FX = GX( 0.0 ); // calculation of initial G(X) value
        MassBalanceDeviations( pmp->N, pmp->L, pmp->A, pmp->Y, pmp->B, pmp->C );
        // Freeing work arrays
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

// Generic simplex method with two sided constraints
// (c) K.Chudnenko 1992
//  SPOS function
// P -
// STR -
// NMB -
// J -
// M -
// AA -
//
void TMulti::SPOS( double *P, int STR[],int NMB[],int J,int M,double AA[])
{
    int I,K;
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

// Generic simplex method with two sided constraints
// (c) K.Chudnenko 1992
//  START function
// T -
// ITER -
// M
// N
// STR -
// NMB -
// J -
// M -
// AA -
// ......
//
void TMulti::START( int T,int *ITER,int M,int N,int NMB[],
                     double GZ,double EPS,int STR[],int *BASE,
                     double B[],double UND[],double UP[],double AA[],double *A,
                     double *Q )
{
    int I,J;

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
                Error("E00IPM: Simplex", "Negative UP[J] value(s) in START() ");
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

// Generic simplex method with two sided constraints
// (c) K.Chudnenko 1992
//  NEW function   parameters
//  OPT -
// .....
//
void TMulti::NEW(int *OPT,int N,int M,double EPS,double *LEVEL,int *J0,
                  int *Z,int STR[], int NMB[], double UP[],
                  double AA[], double *A)
{
    int I,J,J1;
    double MAX,A1;
    double *P;
    P= new double[M+1];
    ErrorIf( !P, "Simplex", "At NEW memory alloc error. ");
    memset(P, 0, sizeof(double)*(M+1) );
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


// Generic simplex method with two sided constraints
// (c) K.Chudnenko 1992
//  WORK function
// GZ -
//.....
//
void TMulti::WORK(double GZ,double EPS,int *I0, int *J0,int *Z,int *ITER,
                   int M, int STR[],int NMB[],double AA[],
                   int BASE[],int *UNO,double UP[],double *A,double Q[])
{
    double MIM,A1;
    int UN,J,I;
    double *P;
    P=  new double[M+1];
    ErrorIf( !P, "Simplex", "At WORK memory alloc error. ");
    memset(P, 0, sizeof(double)*(M+1) );
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
    if(fabs(MIM)<fabs(UP[*J0-1])||fabs(UP[*J0-1])<EPS)
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
            UP[*J0-1]=-UP[*J0-1];
        }
    }
    delete[] P;
}

// Generic simplex method with two sided constraints
// (c) K.Chudnenko 1992
//  FIN function
//
//  EPS -
// ....
//
void TMulti::FIN(double EPS,int M,int N,int STR[],int NMB[],
                  int BASE[],double UND[],double UP[],double U[],
                  double AA[],double *A,double Q[],int * /*ITER*/)
{
    int /* K,*/I,J;
    double *P;
    P=  new double[M+1];
    ErrorIf( !P, "Simplex", "At FIN memory alloc error. ");
    memset(P, 0, sizeof(double)*(M+1) );
    for( J=0;J<N;J++)
    {
        if( UP[J]>-EPS)
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
        U[I-1]=-*(A+I);
    delete[] P;
}

// Generic simplex method with two sided constraints
// (c) K.Chudnenko 1992
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
void TMulti::Simplex(int M, int N, int T, double GZ, double EPS,
                      double *UND, double *UP, double *B, double *U,
                      double *AA, int *STR, int *NMB )
{
    int IT=200,I0=0,J0=0,Z,UNO,OPT=0,ITER, i;
    double LEVEL;
    int *BASE=0;
    double *A=0,*Q=0;

    try
    {
        A=  new double[(M+1)*(M+1)];   /* A[pm->N][pm->N],Q[pm->N]; */
        Q=  new double[M+1];
        BASE=  new int[M];          /* BASE[pm->N-1]  */
        ErrorIf( !A || !Q || !BASE, "Simplex", "Memory alloc error ");
        memset(A, 0, sizeof(double)*(M+1)*(M+1) );
        memset(Q, 0, sizeof(double)*(M+1) );
        memset(BASE, 0, sizeof(int)*(M) );

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
         Error( "E01IPM: Simplex",
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

//--------------------- End of ipm_simplex.cpp ---------------------------
