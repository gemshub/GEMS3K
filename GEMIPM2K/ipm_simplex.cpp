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

// Calculation of Simplex initial approximation of the primal vector x
// using modified simplex method with two-side constraints on x
//
void TMulti::SimplexInitialApproximation( )
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
        ErrorIf( !DN || !DU || !B1, "SimplexInitialApproximation()", "Memory alloc error." );
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
        ErrorIf( !AA || !STR || !NMB, "SimplexInitialApproximation()",
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
        Simplex(pmp->N,Q,T,GZ,EPS,DN,DU,B1,pmp->U,AA,STR,NMB);

        // unloading simplex solution into a copy of x vector
        for(i=0;i<pmp->L;i++)
            pmp->Y[i]=DU[i];

        // calculating initial quantities of phases
        TotalPhases( pmp->Y, pmp->YF, pmp->YFA );

        if( TProfil::pm->pa.p.PD > 0)
        {    pmp->FX = GX( 0.0 ); // calculation of initial G(X) value
             MassBalanceResiduals( pmp->N, pmp->L, pmp->A, pmp->Y, pmp->B, pmp->C );
        }
        else
        {
#ifdef Use_qd_real
            unsigned int old_cw;
            fpu_fix_start(&old_cw);
#endif
            pmp->FX = to_double(qdGX( 0.0 )); // calculation of initial G(X) value
#ifdef Use_qd_real
            fpu_fix_end(&old_cw);
#endif
                  qdMassBalanceResiduals( pmp->N, pmp->L, pmp->A, pmp->Y, pmp->B, pmp->C );

         }
        	
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
                Error("E00IPM: Simplex()", "Inconsistent LP problem (negative UP[J] value(s) in START()) ");
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
    ErrorIf( !P, "Simplex()", "At NEW: memory allocation error ");
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
    ErrorIf( !P, "Simplex()", "At WORK: memory allocation error. ");
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
    ErrorIf( !P, "Simplex()", "At FIN: memory allocation error. ");

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
        U[I-1]=-*(A+I);

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
void TMulti::Simplex(long int M, long int N, long int T, double GZ, double EPS,
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
        ErrorIf( !A || !Q || !BASE, "Simplex()", "Memory allocation error ");

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
         Error( "E01IPM: Simplex()",
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
