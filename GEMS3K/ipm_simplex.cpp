//-------------------------------------------------------------------
// $Id$
//
/// \file ipm_simplex.cpp
/// Implementation of parts of the Interior Points Method (IPM) algorithm
/// for convex programming Gibbs energy minimization,including modified
/// simplex LP method with two-side constraints
//
// Copyright (c) 1992-2012 K.Chudnenko, D.Kulik, S.Dmitrieva
// <GEMS Development Team, mailto:gems2.support@psi.ch>
//
// This file is part of the GEMS3K code for thermodynamic modelling
// by Gibbs energy minimization <http://gems.web.psi.ch/GEMS3K/>
//
// GEMS3K is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation, either version 3 of
// the License, or (at your option) any later version.

// GEMS3K is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with GEMS3K code. If not, see <http://www.gnu.org/licenses/>.
//-------------------------------------------------------------------
//

#include "ms_multi.h"
#include "activities.h"
#include "kinetics.h"
#include "v_service.h"

/// Calculation of LPP-based automatic initial approximation of the primal vector x.
/// Use the modified simplex method with two-side constraints on x
//
void TMultiBase::AutoInitialApproximation( )
{
    long int T,Q,*STR=0,*NMB=0;
    long int i,j,k;
    double GZ,EPS,*DN=0,*DU=0,*AA1=0,*B1=0;

    try
    {  // Allocation of work arrays

        for( i=0; i<pm.N; i++) // added SD 15/07/2009
             pm.U[i] = 0.;

        pm.Ec=0;
        Q=pm.L;
        DN= new double[Q];
        DU= new double[Q+pm.N];
        B1= new double[pm.N];
        ErrorIf( !DN || !DU || !B1, "AutoInitialApproximation()", "Memory alloc error." );
        for( i=0; i<pm.N; i++)
             DU[i+Q] = 0.;
        EPS = base_param()->EPS; //  13.10.00  KC  DK
        GZ = 1./EPS;

        T=0; // Calcuation of all non-zero values in A and G arrays
        for(i=0;i<pm.L;i++)
            if(fabs(pm.G[i])>1E-19)
                T++;
        for(j=0;j<pm.N;j++)
            for(i=0;i<pm.L;i++)
                if(fabs(*(pm.A+i*pm.N+j))>1E-19)
                    T++;
//        if( pm.PLIM ) // Setting constraints on x elements
//            Set_DC_limits(  DC_LIM_INIT );
//            Set_DC_limits( false );    // test 29.07.15

        for(i=0;i<Q;i++)
        {
            DN[i]=pm.DLL[i];
            DU[i]=/*1e+6;//*/pm.DUL[i];
        }

        // for(i=Q;i<Q+pm.N;i++) DU[i]=0.;
        // Allocation of arrays on T
        AA1= new double[T];
        STR= new long int[T];
        NMB= new long int[Q+1];
        ErrorIf( !AA1 || !STR || !NMB, "AutoInitialApproximation()",
            "Memory allocation error #2");
        for( k=0; k<T; k++)
         STR[k] = 0;

        // if( wn[W_EQCALC].status )
        //   aSubMod[MD_EQCALC]->ModUpdate(" SIMPLEX Approximation ");

        // Copying vector b
        for(j=0;j<pm.N;j++)
            B1[j]=pm.B[j];
        k=-1;
        for(i=0;i<pm.L;i++)
        {  // Loading non-zero values
            if(fabs(pm.G[i])>1E-19)
            {
                k++;
                AA1[k]=-pm.G[i];
                NMB[i]=k+1;
            }
            else NMB[i]=k+2;
            for(j=0;j<pm.N;j++)
                if(fabs(*(pm.A+i*pm.N+j))>1E-19)
                {
                    k++;
                    AA1[k]=*(pm.A+i*pm.N+j);
                    STR[k]=j+1;
                }
        }
        NMB[Q]=T+1;
        // Calling generic simplex solver
        SolveSimplex(pm.N,Q,T,GZ,EPS,DN,DU,B1,pm.U,AA1,STR,NMB);

        // unloading simplex solution into a copy of x vector
        for(i=0;i<pm.L;i++)
            pm.Y[i]=DU[i];

        // calculating initial quantities of phases
        TotalPhasesAmounts( pm.Y, pm.YF, pm.YFA );
        pm.FX = GX( 0.0 ); // calculation of initial G(X) value
        MassBalanceResiduals( pm.N, pm.L, pm.A, pm.Y, pm.B, pm.C );

        if(ipm_logger->should_log(spdlog::level::trace)) {
            ipm_logger->trace("AutoInitialApproximation C = {}", to_string(pm.C, pm.N));
            ipm_logger->trace("AutoInitialApproximation B = {}", to_string(pm.B, pm.N));
        }
        // Deleting work arrays
        if( DN) delete[]DN;
        if( DU) delete[]DU;
        if( AA1) delete[]AA1;
        if( B1) delete[]B1;
        if( STR) delete[]STR;
        if( NMB) delete[]NMB;
    }
    catch( TError& xcpt )
    {
        if( DN) delete[]DN;
        if( DU) delete[]DU;
        if( AA1) delete[]AA1;
        if( B1) delete[]B1;
        if( STR) delete[]STR;
        if( NMB) delete[]NMB;
        Error( xcpt.title, xcpt.mess );
    }
}

/// Generic simplex method with two sided constraints (c) K.Chudnenko 1992
///  SPOS function
//
void TMultiBase::SPOS( double *P, long int STR[],long int NMB[],long int J,long int M,double AA1[])
{
    long int I,K;
    K=0;
    for(I=0; I<=M; I++)
    {
        if( I==STR[NMB[J]+K-1])
        {
            *(P+I)=AA1[NMB[J]+K-1];
            if( NMB[J]+K+1!=NMB[J+1])
                K++;
        }
        else *(P+I)=0.;
    }
}

/// Generic simplex method with two sided constraints (c) K.Chudnenko 1992
///  START function
//
void TMultiBase::START( long int T,long int *ITER,long int M,long int N,long int NMB[],
           double GZ,double EPS,long int STR[],long int *BASE, double B[],
           double UND[],double UP[],double AA1[],double *A, double *Q )
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
        SPOS(Q, STR, NMB, J, M, AA1);
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
                    AA1[J]=-AA1[J];
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

/// Generic simplex method with two sided constraints (c) K.Chudnenko 1992
///  NEW function
//
void TMultiBase::NEW(long int *OPT,long int N,long int M,double EPS,double *LEVEL,long int *J0,
                  long int *Z,long int STR[], long int NMB[], double UP[],
                  double AA1[], double *A)
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
        SPOS( P, STR, NMB, J-1, M, AA1);
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
        SPOS(P, STR, NMB, J-1, M, AA1);
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


/// Generic simplex method with two sided constraints (c) K.Chudnenko 1992
///  WORK function
//
void TMultiBase::WORK(double GZ,double EPS,long int *I0, long int *J0,long int *Z,long int *ITER,
                   long int M, long int STR[],long int NMB[],double AA1[],
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
    SPOS(P, STR, NMB, J, M, AA1);
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

/// Generic simplex method with two sided constraints (c) K.Chudnenko 1992
///  FIN function
//
void TMultiBase::FIN(double EPS,long int M,long int N,long int STR[],long int NMB[],
                  long int BASE[],double UND[],double UP[],double U[],
                  double AA1[],double *A,double Q[],long int * /*ITER*/)
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
        SPOS( P, STR, NMB, J-1, M, AA1);
        UP[J-1]+=UND[J-1];
        for( I=0;I<=M;I++)
            Q[I]+=UP[J-1]*P[I];
    }
    for( I=1;I<=M;I++)
        U[I-1] -= *(A+I);  // was =- *(A+I)

    delete[] P;
}

/// Generic simplex method with two sided constraints (c) K.Chudnenko 1992.
///  Main function
//
///  \param M  - number of independent components
///  \param N  - number of unknowns
///  \param T  - dimension of a work vector AA[] containing all non-zero
///        values of vector GT[] and A[][] matrix (over lines)
///  \param GZ - Limiting value of the unknown
///  \param EPS - precision (convergence) criterion (default 1e-9)
///  \param UND - vector of lower constraints to unknowns
///  \param UP - input vector of upper constraints to unknowns;
///        output vector of unknowns (simplex solution) (N+M)
///  \param B -  M input values of independent components (bulk composition)
///  \param U -  output vector of the dual solution (M)
///  \param AA - work array (T)
///  \param STR - markup vector of values in AA array (T)
///  \param NMB - indices of values in AA
/// \return 0 if OK;
///         1 if inconsistent input constraints;
///        -1 if memory allocation error;
//
void TMultiBase::SolveSimplex(long int M, long int N, long int T, double GZ, double EPS,
                      double *UND, double *UP, double *B, double *U,
                      double *AA1, long int *STR, long int *NMB )
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
        START( T, &ITER, M, N, NMB, GZ, EPS, STR, BASE, B,  UND, UP, AA1, A, Q );

        for( i=0; i<IT; i++ )   // while(1) fixed  03.11.00
        {
            NEW( &OPT, N, M,EPS, &LEVEL, &J0, &Z, STR, NMB, UP, AA1, A);
            if( OPT)
                goto FINISH;  // Converged
            WORK( GZ, EPS, &I0, &J0, &Z, &ITER, M, STR, NMB, AA1, BASE, &UNO, UP, A, Q);
            if( UNO)
                goto FINISH; // Solution at boundary of the constraints polyhedron
        }
        if( EPS > 1.0e-6 )
        {
         Error( "E01IPM: SolveSimplex()",
             "LP solution cannot be obtained with sufficient precision" );
        }
FINISH: FIN( EPS, M, N, STR, NMB, BASE, UND, UP, U, AA1, A, Q, &ITER);
        delete[] A;
        delete[] Q;
        delete[] BASE;
    }
    catch( TError& xcpt )
    {
        if( A) delete[]A;
        if( Q) delete[]Q;
        if( BASE) delete[]BASE;
        Error( xcpt.title, xcpt.mess);
    }

    // Done
}


//-----------------------------------------------------------------------
/// Main call to GEM IPM calculation of equilibrium state in MULTI
/// (with internal re-scaling of the system).
//
double TMultiBase::CalculateEquilibriumState(  long int& NumIterFIA, long int& NumIterIPM )
{
  double ScFact=1.;

  InitalizeGEM_IPM_Data();

  pm.t_start = clock();
  pm.t_end = pm.t_start;
  pm.t_elap_sec = 0.0;
  pm.ITF = pm.ITG = 0;

 // New: Run of TKinMet class library
  ipm_logger->trace("kMM: {}  ITau: {}   kTau: {}   kdT: {}", pm.pKMM, pm.ITau, pm.kTau, pm.kdT);
  if( pm.pKMM < 2 )
  {
    if( pm.ITau < 0 || pm.pKMM != 1 )
    {
      /*  KMretCode = */ CalculateKinMet( LINK_TP_MODE ); // Re-create TKinMet class instances
        pm.ITau = 0; pm.pKMM = 1;
      /*  KMretCode = */ CalculateKinMet( LINK_IN_MODE ); // Initial state calculation of rates
    }
//    if( pm.ITau == 0 )
//    {
//        KMretCode = CalculateKinMet( LINK_IN_MODE ); // Initial state calculation of rates
//    }
    else if( pm.ITau >= 0 ) {
      /*  KMretCode = */ CalculateKinMet( LINK_PP_MODE ); // Calculation of rates and metast.constraints at time step
    }
//  switch(KMretCode)
//  {
//        case 0L:
//
//  }
//  to_text_file( "MultiDump1.txt" );   // Debugging
  }

    if( base_param()->DG > 1e-5 )
    {
        ScFact = SystemTotalMolesIC();
        ScaleSystemToInternal( ScFact );
    }

try{
       switch( pm.tMin )
       {
       case  A_TV_:
           {   
	     // kg44 this is not yet properly implemented ....removed it
	     /*
               GoldenSection gsData( pm.Pai[0], pm.Pai[1], pm.Pai[3], pm.Fdev1[1], A_P);
               pm.P = gsData.getMinimum();
               A_P(pm.P);
	     */
           }
           break;
       case  U_SV_:
           {      
	     // kg44 this is not yet properly implemented...removed it 
	     /*
                  GoldenSectionTwo gsData( pm.Tai[0], pm.Tai[1], pm.Tai[3],
                                             pm.Pai[0], pm.Pai[1], pm.Pai[3],
                                             pm.Fdev1[1], U_TP);
                  gsData.getMinimum();
                  pm.TC = gsData.getMinX();
                  pm.P = gsData.getMinY();
                  U_TP(pm.TC, pm.P);
	     */
              }
            break;
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

      if( base_param()->DG > 1e-5 )
         RescaleSystemFromInternal( ScFact );
//      to_text_file( "MultiDump2.txt" );   // Debugging

      NumIterFIA = pm.ITF;
      NumIterIPM = pm.ITG;
      pm.t_end = clock();
      pm.t_elap_sec = double(pm.t_end - pm.t_start)/double(CLOCKS_PER_SEC);

     Error( xcpt.title, xcpt.mess);
  }

  if( base_param()->DG > 1e-5 )
       RescaleSystemFromInternal(  ScFact );
//  to_text_file( "MultiDump3.txt" );   // Debugging

  NumIterFIA = pm.ITF;
  NumIterIPM = pm.ITG;
  pm.t_end = clock();
  pm.t_elap_sec = double(pm.t_end - pm.t_start)/double(CLOCKS_PER_SEC);

  return pm.t_elap_sec;
}


/// Calculate total IC mole amounts in b vector and
/// return the scaling factor
double TMultiBase::SystemTotalMolesIC( )
{
  double ScFact, mass_temp = 0.0;

  for( int i=0; i<pm.N - pm.E; i++ ) //?????
         mass_temp +=pm.B[i];

  pm.TMols = mass_temp;

  pm.SMols = base_param()->DG;
  ScFact = pm.SMols/pm.TMols;

  return ScFact;
}

/// Resizes MULTI (GEM IPM work structure) data into internally scaled constant mass
void TMultiBase::ScaleSystemToInternal(  double ScFact )
{
 long int i, j, k;

  pm.SizeFactor = ScFact;
  pm.MBX *= ScFact;
  pm.FitVar[0] *= ScFact;  // added: bugfix by DK 31.05.2010
  pm.VXc *= ScFact;
  pm.HXc *= ScFact;
  pm.HX_ *= ScFact;
  pm.FX  *= ScFact;
  pm.Yw  *= ScFact;  // added 08.06.10 DK

  for( j=0; j<pm.L; j++ )
  {
    if(	pm.DUL[j] < 1e6  )
       pm.DUL[j] *= ScFact;

    // if( pm.DLL[j] > 0.0  )
       pm.DLL[j] *= ScFact;

        pm.Y[j] *= ScFact;
        pm.X[j] *= ScFact;
        pm.XY[j] *= ScFact;
        pm.XU[j] *= ScFact;
  }

  for( i=0; i<pm.N; i++ )
  {
        pm.B[i] *= ScFact;
  //      pm.C[i] *= ScFact;
        pm.BFC[i] *= ScFact;
  }

  for( k=0; k<pm.FI; k++ )
  {
    pm.XFs[k] *= ScFact;
    pm.XF[k] *= ScFact;
    pm.YF[k] *= ScFact;
    pm.FVOL[k] *= ScFact;
    pm.FWGT[k] *= ScFact;
  }

  for( k=0; k<pm.FIs; k++ )
  {
      pm.XFA[k] *= ScFact;
      pm.YFA[k] *= ScFact;

      if( pm.PUL )
        if( pm.PUL[k] < 1e6  )
         pm.PUL[k] *= ScFact;

      if( pm.PLL )
      // if( pm.PLL[k] > 0.0  )
         pm.PLL[k] *= ScFact;

      for( i=0; i<pm.N; i++ )
           pm.BF[k*pm.N+i] *= ScFact;
  }
 if( pm.FIat > 0 /*&& pm.Lads > 0*/ && pm.FIs > 0 )
  {
  for( k=0; k<pm.FIs; k++ )
     for( j=0; j<MST; j++ )
          {
              pm.XetaA[k][j] *= ScFact;
              pm.XetaB[k][j] *= ScFact;
              pm.XetaD[k][j] *= ScFact;
              pm.XFTS[k][j] *= ScFact;
          }
}

 pm.SizeFactor = ScFact;
}

/// Re-scaling the internal constant-mass MULTI system definition
/// back to real size
void TMultiBase::RescaleSystemFromInternal(  double ScFact )
{
 long int i, j, k;

  pm.MBX /= ScFact;
  pm.FitVar[0] /= ScFact;  // added: bugfix by DK 31.05.2010
  pm.VXc /= ScFact;
  pm.HXc /= ScFact;
  pm.HX_ /= ScFact;
  pm.FX  /= ScFact;
  // added SV
  //pm.YFk /= ScFact;
  pm.Yw /= ScFact;  // added 08.06.10 DK

  for( j=0; j<pm.L; j++ )
  {
    if(	pm.DUL[j] < 1e6  ) {
       pm.DUL[j] /= ScFact;
    }
    // if( pm.DLL[j] > 0.0  )
       pm.DLL[j] /= ScFact;

        pm.Y[j] /= ScFact;
        pm.X[j] /= ScFact;
        pm.XY[j] /= ScFact;
        pm.XU[j] /= ScFact;
    //    pm.MU[j] /= ScFact;
    //    pm.W[j] /= ScFact;
  }

  for( i=0; i<pm.N; i++ )
  {
        pm.B[i] /= ScFact;
        pm.C[i] /= ScFact;
        pm.BFC[i] /= ScFact;
  }

  for( k=0; k<pm.FI; k++ )
  {
    pm.XFs[k] /= ScFact;
    pm.XF[k] /= ScFact;
    pm.YF[k] /= ScFact;
    pm.FVOL[k] /= ScFact;
    pm.FWGT[k] /= ScFact;
  }

  for( k=0; k<pm.FIs; k++ )
  {
      pm.XFA[k] /= ScFact;
      pm.YFA[k] /= ScFact;

      if( pm.PUL )
        if( pm.PUL[k] < 1e6  )
         pm.PUL[k] /= ScFact;

      if( pm.PLL )
      // if( pm.PLL[k] > 0.0  )
         pm.PLL[k] /= ScFact;

      for( i=0; i<pm.N; i++ )
           pm.BF[k*pm.N+i] /= ScFact;
  }

  if( pm.FIat > 0 /*&& pm.Lads > 0*/ && pm.FIs > 0 )
   for( k=0; k<pm.FIs; k++ )
          for( j=0; j<MST; j++ )
          {
              pm.XetaA[k][j] /= ScFact;
              pm.XetaB[k][j] /= ScFact;
              pm.XetaD[k][j] /= ScFact;
              pm.XFTS[k][j]  /= ScFact;
          }

  pm.SizeFactor = 1.;   // using in TNode class
}

//========================================================================================
// Multi initialization part 10/05/2010

// Before Calculations
/// Calculation by IPM (preparing for calculation, unpacking data)
/// In IPM
void TMultiBase::InitalizeGEM_IPM_Data() // Reset internal data formerly MultiInit()
{

   MultiConstInit();

   initalizeGEM_IPM_Data_GUI();

   if( pm.pTPD < 2 )
    {
        //  pm.pTPD == 0 => changed T or P;   pm.pTPD == 1 => changed system;
      DC_LoadThermodynamicData(); // Loading thermodynamic data into MULTI structure
     }
    Alloc_internal();
    Alloc_uDD( pm.N );      // Added 06.05.2011 DK

  // calculate mass of the system
  pm.MBX = 0.0;
  for(int i=0; i<pm.N; i++ )
     pm.MBX += pm.B[i] * pm.Awt[i];
  pm.MBX /= 1000.;

   RescaleToSize( true );  // Added to set default cutoffs/inserts 30.08.2009 DK

   if(  pm.pNP )
    {  //  Smart Initial Approximation mode
       long int j,k;

       loadData( false );  // unpack SysEq record into MULTI

     bool AllPhasesPure = true;   // Added by DK on 09.03.2010
     // checking if all phases are pure
     for( k=0; k < pm.FI; k++ )
     if( pm.L1[k] > 1 )
        AllPhasesPure = false;

     for( j=0; j< pm.L; j++ )
         pm.X[j] = pm.Y[j];
//       pm.IC = 0.;  //  Problematic statement!  blocked 13.03.2008 DK
      TotalPhasesAmounts( pm.X, pm.XF, pm.XFA );
      CalculateConcentrations( pm.X, pm.XF, pm.XFA);  // 13.03.2008  DK
       // test multicomponent phases and load data for mixing models
       if( pm.FIs && AllPhasesPure == false )
       {
           // Load activity coeffs for phases-solutions
         int jb, je=0;
         for( int kk=0; kk<pm.FIs; kk++ )
         { // loop on solution phases
            jb = je;
            je += pm.L1[kk];
            // Load activity coeffs for phases-solutions
            for( j=jb; j< je; j++ )
            {
               pm.lnGmo[j] = pm.lnGam[j];
               if( fabs( pm.lnGam[j] ) <= 84. )
      //                pm.Gamma[j] = exp( pm.lnGam[j] );
                      pm.Gamma[j] = PhaseSpecificGamma( j, jb, je, kk, 0 );
               else pm.Gamma[j] = 1.0;
             } // j
          }  // kk
       }
    }
}


/// Setup/copy flags and thresholds for numeric modules to TMulti structure.
/// Do it before calculations
void TMultiBase::MultiConstInit() // from MultiRemake
{
  pm.FI1 = 0;
  pm.FI1s = 0;
  pm.FI1a = 0;
  pm.ITF = 0; pm.ITG = 0;
  pm.PD = base_param()->PD;
  pm.Ec = pm.K2 = pm.MK = 0;
  pm.W1 = 0;
  pm.is = 0;
  pm.js = 0;
  pm.next  = 0;
  pm.ln5551 = log( H2O_mol_to_kg );             // constant corrected 30.08.2008
  pm.lowPosNum = Min_phys_amount;               // = 1.66e-24 mol
  pm.logXw = -16.;
  pm.logYFk = -9.;
  pm.DXM = base_param()->DK;

  //  ???????
  pm.FX = 7777777.;
  if( pm.pH < -15. || pm.pH > 16.  )   // Check for trash in pH - bugfix 19.06.2013
      pm.pH = pm.Eh = pm.pe = 0.0;
  pm.YMET = 0;                      // always 0.0 ????
  pm.PCI = 1.0;
  pm.FitVar[4] = 1.0;

    multiConstInit_PN();

}

/// Calculation by IPM (internal step initialization)
void TMultiBase::GEM_IPM_Init()
{
   int i,j,k;

   for( i=0; i<pm.N; i++ )
     pm.Uefd[i] = 0.;

   bool AllPhasesPure = true;   // Added by DK on 09.03.2010
  // checking if all phases are pure
  for( k=0; k < pm.FI; k++ )
    if( pm.L1[k] > 1 )
        AllPhasesPure = false;

   if(!pm.pNP) // Simplex initial approximation to be done
    {
        for( j=0; j<pm.L; j++ )
        {                           // cleaning work vectors
                pm.X[j] = pm.Y[j] = pm.lnGam[j] = pm.lnGmo[j] = 0.0;
                pm.Gamma[j] = 1.0;
                pm.MU[j] = 0.;
                pm.XU[j] = 0.;
                pm.EMU[j] = 0.;
                pm.NMU[j] = 0.;
                pm.W[j] = 0.;
                pm.F[j] = 0.;
                pm.F0[j] = 0.;
           }
    }

    // recalculating kinetic restrictions for DC amounts
//     if( pm.pULR && pm.PLIM )
//          Set_DC_limits(  DC_LIM_INIT );

    if( pm.FIs && AllPhasesPure == false )   /// line must be tested !pm.FIs
    {
        GEM_IPM_Init_gui1();
      
        // Calc Eh, pe, pH,and other stuff
       if( pm.E && pm.LO && pm.pNP )
       {
            CalculateConcentrations( pm.X, pm.XF, pm.XFA);
            IS_EtaCalc();
            if( pm.Lads )  // Calling this only when sorption models are present
            {
               int jb, je=0;
               for(int kk=0; kk<pm.FIs; kk++ )
               { // loop on solution phases
                  jb = je;
                  je += pm.L1[kk];
                  if( pm.PHC[kk] == PH_POLYEL || pm.PHC[kk] == PH_SORPTION )
                  {
                       if( pm.PHC[0] == PH_AQUEL && pm.XF[kk] > pm.DSM
                           && (pm.XFA[0] > pm.ScMinM && pm.XF[0] > pm.XwMinM )) // fixed 30.08.2009 DK
                           GouyChapman( jb, je, kk );  // getting PSIs - electrical potentials on surface planes
                  }
               }
           }
       }
       CalculateActivityCoefficients( LINK_TP_MODE);
       // Computing DQF, FugPure and G wherever necessary; Activity coeffs are restored from lnGmo
    }
    else {  // no multi-component phases
        pm.PD = 0;
        pm.pNP = 0;
        pm.pIPN = 1;
    }

    // recalculating kinetic restrictions for DC amounts
//     if( pm.pULR && pm.PLIM )
//          Set_DC_limits(  DC_LIM_INIT );
//        Set_DC_limits( true );  // test 29.07.15

  GEM_IPM_Init_gui2();
}


TSolMod * TMultiBase::pTSolMod (int xPH){
    return this->phSolMod[xPH];
}

//===========================================================================================
// Calls to minimization of other system potentials (A, ...)

/// Calls to minimization of other system potentials - HelmholtzEnergy.
/// Calc function for Method of golden section (only in GEMS ).
double TMultiBase::HelmholtzEnergy( double x )
{
    pm.P = x;
    DC_LoadThermodynamicData(); // from lookup arrays
    GibbsEnergyMinimization();
    return( pm.VXc - pm.VX_ ); // VXc current value, VX_ value from key
}
//
// kg44: not properly implemented...removed it
/*
double A_P( double x, double )
{
#ifndef IPMGEMPLUGIN
  return TProfil::pm->HelmholtzEnergy(x);
#endif
}
*/

/// Calls to minimization of other system potentials - InternalEnergy.
/// Calc function for Method of golden section (only in GEMS ).
double TMultiBase::InternalEnergy( double TC, double P )
{
    pm.P = P;
    pm.TC = TC;
    DC_LoadThermodynamicData(); // from lookup arrays
    GibbsEnergyMinimization();
    return( (pm.VXc - pm.VX_)+(pm.SXc - pm.SX_) ); // VXc current value, VX_ value from key
}

// kg44: not properly implemented...removed it
/*
double U_TP( double TC, double P)
{
#ifndef IPMGEMPLUGIN
  return TProfil::pm->InternalEnergy(  TC,  P );
#endif
}
*/


//--------------------- End of ipm_simplex.cpp ---------------------------
