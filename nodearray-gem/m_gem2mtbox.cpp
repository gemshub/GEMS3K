//-------------------------------------------------------------------
// $Id: m_gem2mtbox.cpp 968 2007-12-13 13:23:32Z gems $
//
// Implementation of TInteg/TGEM2MT classes, calculation functions
//
// Rewritten from C to C++ by S.Dmytriyeva  
// Copyright (C) 1995, 2008  S. Dmytriyeva, D.Kulik 
//
// This file is part of a GEM-Selektor library for thermodynamic
// modelling by Gibbs energy minimization
// Uses: GEM-Vizor GUI DBMS library, gems/lib/gemvizor.lib
//
// This file may be distributed under the terms of the GEMS-PSI
// QA Licence (GEMSPSI.QAL)
//
// See http://gems.web.psi.ch/ for more information
// E-mail: gems2.support@psi.ch
//-------------------------------------------------------------------
//
#include <math.h>
#include <stdio.h>
#include <iomanip>
#include "io_arrays.h"


#ifndef IPMGEMPLUGIN

#include "m_gem2mt.h"
#include "nodearray.h"
#include "service.h"
#include "visor.h"

#else

#include <time.h>
#include "ms_gem2mt.h"
#include "nodearray.h"

#endif

#define dMB( q, i) ( dm[ (q)*mtp->Nf + (i)] ) 

#define MB( q, i)  ( m[ (q)*mtp->Nf + (i)] )

#define g(q,f,i)   ( mtp->gfc[ (q)*(mtp->nPG*mtp->Nf)+ (f)*mtp->Nf + (i)] )
#define y(q,f,i)   ( mtp->yfb[ (q)*(mtp->nPG*mtp->Nf)+ (f)*mtp->Nf + (i)] )

#define ord(kk)    ( ROUND(mtp->FDLf[kk][0]) )
#define v(kk)      ( (mtp->FDLf[kk][1]) )

#define H(f, i)    ( mtp->BSF[ (f) * mtp->Nf + ( i )] )

// Returns MGP index from MGP identifier MGPid
//   or -1 if the identifier was not found in the MGP id list
int TGEM2MT::LookUpXMGP( const char* MGPid )
{   
	int found = -1;
	// Check if the first character is 0 1 2 3 4 5 6 7 8 9 
	// If so, this is index of elemental flux with composition from BSF table
	for( int f=0; f < mtp->nPG; f++ )
	{
		if( strncmp( mtp->MPGid[f], MGPid, MAXSYMB ) )
			continue; 
	    found = f;
	    break;
	}
	return found;
}

// calculate 1-step from system of equation 
void TGEM2MT::Solut( double *m, double *dm, double t )
{ 
  int kk, q, p, f, fe=-1, i;
  char FLXid[MAXSYMB+1], MGPid[MAXSYMB+1];
  int FLXorder, FLXtype;
  double fRate;
  bool sinkOut; 
  
  // Zeroing element mass derivatives off
  for( q=0; q <mtp->nC; q++ )
	for(i =0; i< mtp->Nf; i++ )
      dMB(q,i) = 0.;

  for(kk=0; kk<mtp->nFD; kk++ )  // Looking through the list of fluxes 
  {
    fe = -1;
	q = mtp->FDLi[kk][0];    // index of source box
	p = mtp->FDLi[kk][1];    // index of receiving box 
//	FLXorder = mtp->FDLf[kk][0];  // flux order code
FLXorder = ord(kk);
//	FLXtype = mtp->FDLop[kk][1];   // flux type code
	strncpy( MGPid, mtp->FDLmp[kk], MAXSYMB );  
	MGPid[MAXSYMB] = 0;                           // MGP identifier
	f = LookUpXMGP( MGPid );                      // Getting MGP index 
    if( f<0 ) // Elemental production flux? 	
    {
    	sscanf( MGPid, "%d", &fe );  // Reading BSF row index 
        if( fe < 0 || fe >= mtp->nSFD )
           Error( "BOXFLUX:", "Wrong MPG identifier or BSF row index!" );	
    }
    strncpy( FLXid, mtp->FDLid[kk], MAXSYMB );  
	FLXid[MAXSYMB] = 0; 						  // Flux identifier
		
	if( q >= 0 && f >= 0 )
	{            // Normal MGP flux from box q to box p
// NB: Negative v(f) means "production" in q box and "consumption" in p box 
		if( p < 0 ) 
			sinkOut = true; // This is a sinkout flux from box q to nowhere (if v > 0)
		                    //  or production of MGP in box q (if v < 0) 
        switch( FLXorder )
        {
          case 0:  // Zero-order flux 
        	 for(i=0; i<mtp->Nf; i++ )
        	 {	  
        	    fRate = v(kk) * g(q,f,i);
        		dMB(q,i) -=  fRate;
        	    if( !sinkOut)
        	       dMB(p,i) +=  fRate;
        	 }
        	 break; 
          case 1:  // First-order to source flux 
             for(i=0; i<mtp->Nf; i++ )
	         {	  
                 fRate = v(kk) * g(q,f,i) * MB(q,i); 
            	dMB(q,i) -=  fRate;
            	if( !sinkOut)
            	   dMB(p,i) +=  fRate;
	         }
             break; 
          case 2:  // Second-order to source flux 
             for(i=0; i<mtp->Nf; i++ )
	         {	  
                 fRate = v(kk) * g(q,f,i) * MB(q,i)*MB(q,i);
            	dMB(q,i) -=  fRate;
            	if( !sinkOut)
            	   dMB(p,i) +=  fRate;
	         }
             break;
          case 3:  // Second-order to source and sink flux 
        	         // (proportional to product of masses in both boxes)  
             if( sinkOut )
            	 break; 
        	 for(i=0; i<mtp->Nf; i++ )
	         {	  
                 fRate = v(kk) * g(q,f,i) * MB(q,i)*MB(p,i);
            	dMB(q,i) -=  fRate;
            	dMB(p,i) +=  fRate;
	         }
             break;
           case -1:  // First-order to receiver flux only 
             if( sinkOut )
              	 break;        
        	 for(i=0; i<mtp->Nf; i++ )
          	 {	  
                 fRate = v(kk) * g(q,f,i) * MB(p,i);
               	dMB(q,i) -=  fRate;
               	dMB(p,i) +=  fRate;
          	 }
             break;
          case -2:  // Second-order to receiver flux 
        	  if( sinkOut )
        	     break;                	          	 
        	 for(i=0; i<mtp->Nf; i++ )
          	 {	  
                 fRate = v(kk) * g(q,f,i) * MB(p,i)*MB(p,i);
                dMB(q,i) -=  fRate;
           	    dMB(p,i) +=  fRate;
          	 }
             break;                      
          default:  // Other orders or flux types - to be done! 
        	 break; 
        }
        ;  
    } 
	else { 
		if( q >= 0 && f < 0 )
		 {    // This is an elemental flux ( fe row in BSF table )         
	// NB: Negative v(f) means "production" in q box and "consumption" in p box 
			if( p < 0 ) 
				sinkOut = true; // This is an elemental sinkout flux from box q to nowhere (if v > 0)
			                    //  or elemental production in box q (if v < 0) 
	        switch( FLXorder )
	        {
	          case 0:  // Zero-order flux 
	        	 for(i=0; i<mtp->Nf; i++ )
	        	 {	  
	        	    fRate = v(kk) * H(fe,i);
	        		dMB(q,i) -=  fRate;
	        	    if( !sinkOut)
	        	       dMB(p,i) +=  fRate;
	        	 }
	        	 break; 
	          case 1:  // First-order to source flux 
	             for(i=0; i<mtp->Nf; i++ )
		         {	  
	            	 fRate = v(kk) * H(fe,i) * MB(q,i); 
	            	dMB(q,i) -=  fRate;
	            	if( !sinkOut)
	            	   dMB(p,i) +=  fRate;
		         }
	             break; 
	          case 2:  // Second-order to source flux 
	             for(i=0; i<mtp->Nf; i++ )
		         {	  
	            	 fRate = v(kk) * H(fe,i) * MB(q,i)*MB(q,i);
	            	dMB(q,i) -=  fRate;
	            	if( !sinkOut)
	            	   dMB(p,i) +=  fRate;
		         }
	             break;
	          case 3:  // Second-order to source and sink flux 
	        	         // (proportional to product of masses in both boxes)  
	             if( sinkOut )
	            	 break; 
	        	 for(i=0; i<mtp->Nf; i++ )
		         {	  
	        		 fRate = v(kk) * H(fe,i) * MB(q,i)*MB(p,i);
	            	dMB(q,i) -=  fRate;
	            	dMB(p,i) +=  fRate;
		         }
	             break;
	           case -1:  // First-order to receiver flux only 
	             if( sinkOut )
	              	 break;        
	        	 for(i=0; i<mtp->Nf; i++ )
	          	 {	  
	        		 fRate = v(kk) * H(fe,i) * MB(p,i);
	               	dMB(q,i) -=  fRate;
	               	dMB(p,i) +=  fRate;
	          	 }
	             break;
	          case -2:  // Second-order to receiver flux 
	        	  if( sinkOut )
	        	     break;                	          	 
	        	 for(i=0; i<mtp->Nf; i++ )
	          	 {	  
	        		 fRate = v(kk) * H(fe,i) * MB(p,i)*MB(p,i);
	                dMB(q,i) -=  fRate;
	           	    dMB(p,i) +=  fRate;
	          	 }
	             break;                      
	          default:  // Other orders or flux types - nothing to be done yet!
	        	 break; 
	        }
	        ;
		}
	 }
  }   // end kk loop
}

#undef dMB
#undef MB

#define Mb( q, i)  ( mtp->MB[ (q)*mtp->Nf + (i)])

#define dMb( q, i)  (mtp->dMB[(q)*mtp->Nf + (i)])

// Calculate new equilibrium states in the boxes for tcur = t
//  Ni 
//  pr
//  tcur
//  step
//
bool
TGEM2MT::CalcNewStates(  int Ni, int pr, double tcur, double step)
{
  int q, i, f, k;
  double MGPfactor = 0.;
  bool iRet = true;
  FILE* diffile = NULL;

#ifndef IPMGEMPLUGIN
       iRet = pVisor->Message( window(), GetName(),
           "Calculating Reactive Mass Transport Box-Flux Simulation \n"
           "Please, wait (may take long)...", nstep, mtp->ntM );
       if( iRet )
        Error("GEM2MT Flux-box model", "Cancel by user");
#endif
  
  mtp->dTau = step; 
  mtp->cTau = tcur;
  
/*  if( mtp->PvMO != S_OFF )
 {
   // Preparations: opening output files for monitoring 1D profiles
    diffile = fopen( "ICdif-log.dat", "w+" );   //  Element amount diffs for t and t-1
    ErrorIf( !diffile, "GEM2MT Flux-box model",
    "Error writing monitoring file (ICdif-log.dat)");
 }
*/
 clock_t t_start, t_end, t_out, t_out2;
 clock_t outp_time = (clock_t)0;
 t_start = clock();

 if( Ni >= 0)
 { // Set up new reservoir states at cTau
   for( q=0; q <mtp->nC; q++ )
	 for(i =0; i< mtp->Nf; i++ )
	 {
		node1_bIC(q, i) += dMb( q, i) / nodeCH_ICmm( i ) * mtp->dTau;
		if( i == mtp->Nf-1 )
			node1_bIC(q, i) = 0.0;  // Provisorial - zeroing-off charge
		else if( node1_bIC(q, i) < 1e-12 )   
   			     node1_bIC(q, i) = 1e-12;    // preventing negative amounts of ICs
	 }
 }  

 
// Calculate new reservoir states at tcur	
// Calculation of chemical equilibria in all nodes at the beginning
// with the Simplex initial approximation
 if( mtp->PvSIA != S_ON )  
     CalcIPM( NEED_GEM_AIA, 0, mtp->nC, diffile );
 else 
	 CalcIPM( NEED_GEM_PIA, 0, mtp->nC, diffile );
 
 if( Ni >= 0 )
 { // Here one has to compare old and new equilibrium phase assemblage
   // and pH/pe in all nodes and decide if the time step was Ok or it
   // should be decreased. If so then the nodes from C0 should be
   // copied to C1 (to be implemented)
 

  // Output of the results if step accepted
   if( mtp->PvMO != S_OFF && diffile )
   {
    t_out = clock();
    na->logDiffsIC( diffile, mtp->ct, mtp->cTau/(365*86400), mtp->nC, 10 );
        // logging differences after the MT iteration loop
    t_out2 = clock();
    outp_time += ( t_out2 -  t_out);
  }
 }
#ifndef IPMGEMPLUGIN
   // time step accepted - Copying nodes from C1 to C0 row
      pVisor->Update();
      CalcGraph();
#endif
   
  // copy node array for T0 into node array for T1
  mtp->oTau = mtp->cTau;
  copyNodeArrays();

  // Calculation of current box reactive IC masses in g   
     for( q=0; q <mtp->nC; q++ )
		 for(i=0; i<mtp->Nf; i++ )
			 Mb( q, i) = node1_bIC( q, i ) * nodeCH_ICmm( i );
  
  // Calculation of MPG bulk compoisitions
     for( q=0; q <mtp->nC; q++ )
	   for(f=0; f<mtp->nPG; f++ )
		 for(i=0; i<mtp->Nf; i++ )
   		   y(q,f,i) = 0.0;
     for( q=0; q <mtp->nC; q++ )
	   for(f=0; f<mtp->nPG; f++ )
	     for( k=0; k<mtp->FIf; k++) 
	     {  
	    	MGPfactor = mtp->PGT[k*mtp->nPG+f];
	        if( fabs(MGPfactor) > 0.0 )   // For now, moles only!!!
	        {  	  
	    	  if( k < na->pCSD()->nPSb )
	    	     for(i=0; i<mtp->Nf; i++ )
			        y(q,f,i) += node1_bPS( q, k, i ) * MGPfactor;
	    	  else
	    		 for(i=0; i<mtp->Nf; i++ )
	 			    y(q,f,i) += node1_bPH( q, k, i ) * MGPfactor;
	    	  y(q,f,mtp->Nf-1) = 0.0;    // provisorial - zeroing-off charge
	        } 
	     }
  //  Calculation of MPG IC distribution coefficients   
     for( q=0; q <mtp->nC; q++ )
	   for(f=0; f<mtp->nPG; f++ )
	   { for(i=0; i<mtp->Nf; i++ )
		    	 g(q,f,i) = y(q,f,i)/node1_bIC( q, i );
			  g(q,f,mtp->Nf-1) = 0.0;    // Provisorial for charge
	   } 
     
   return iRet;  
}

#undef dMB
#undef MB

//Calculate record
bool TGEM2MT::CalcBoxModel( char mode )
{
  try{
	 
    bool iRet = false;

    // Init part ?????  
    // mtp->dx = mtp->cLen/mtp->nC;
    mtp->dTau = mtp->Tau[STEP_];;
    mtp->oTau =  mtp->cTau = mtp->Tau[START_];
    // mtp->cTau = 0;
    mtp->ct = 0;
#ifndef IPMGEMPLUGIN
    mtp->gfc = (double *)aObj[ o_mtgc].Alloc(  mtp->nC*mtp->nPG, mtp->Nf, D_);
#else
    if( mtp->gfc )
	    delete[] mtp->gfc; 
    mtp->gfc = new double[mtp->nC * mtp->nPG * mtp->Nf];
#endif
    if( mtp->yfb )
	    delete[] mtp->yfb; 
    mtp->yfb = new double[mtp->nC * mtp->nPG * mtp->Nf];
    if( tt )
    	delete[] tt; 
    tt = new double[mtp->nC * mtp->Nf][9];
    
#ifndef IPMGEMPLUGIN
      if(  mtp->PvMSg != S_OFF && vfQuestion(window(),
             GetName(), "Use graphic monitoring?") )
        {
            RecordPlot( 0 );
        }
#endif

  // calculate inital states
  // Calculation of chemical equilibria in all nodes at the beginning
  // with the Simplex initial approximation
  CalcNewStates( -1,  0, mtp->cTau, mtp->dTau );  
     
  // calc part  
   nfcn = nstep = naccept = nrejct = 0;
   
#ifndef IPMGEMPLUGIN
       iRet = pVisor->Message( window(), GetName(),
           "Calculating Reactive Mass Transport (RMT)\n"
           "Please, wait (may take long)...", nstep, mtp->ntM );
       if( iRet )
        Error("GEM2MT Flux-box model", "Cancel by user");
#endif
  INTEG( 1e-3, /*mtp->cdv,*/ mtp->dTau, mtp->Tau[START_], mtp->Tau[STOP_] );

#ifndef IPMGEMPLUGIN
    pVisor->CloseMessage();
#endif
    
  }
  catch( TError& xcpt )
   {
#ifndef IPMGEMPLUGIN
       vfMessage(window(), xcpt.title, xcpt.mess);
#else
       cerr << xcpt.title.c_str() << "  " <<  xcpt.mess.c_str() << endl;
#endif
       return 1;
   }
  catch(...)
  {
#ifndef IPMGEMPLUGIN
       vfMessage(window(), "CalcBoxModel", "Unknown exception");
#else
       cerr << "CalcBoxModel Unknown exception" << endl;
#endif
      return -1;
  }
  return 0;
}

//--------------------------------------------------------------------
// Integration process

const int NMAX = 800;
const int KM = 8;
const double UROUND = 1.73e-18;
const double FAC1 = 2.e-2;
const double FAC2 = 4.0;
const double FAC3 = .9;
const double FAC4 = .8;
const double SAFE1 = .65;
const double SAFE2 = .94;
const double MAXSTEP = 1.7;
//const int MAXINTEGEXPR = 200;

//double *x;
//double *dx;
//double *tv;

//static double tt[9][ MAXINTEGEXPR ]; change 14/12/2007 [ MAXINTEGEXPR ][9]
static double hh[9];
static double w[9];
static double err;
static double epsd4;
static int    nj[9]={2,4,6,8,10,12,14,16,18};
static double a1[9]={3.0,7.0,13.0,21.0,31.0,43.0,57.0,73.0,91.0};


// internal point j calculation
void TGEM2MT::MIDEX( int j, double t, double h )
{
    double *z1=0, *z2=0, *dz=0;
    double hj, scal, fac, facmin, expo, ys, v1, v2;
    int i,m,mm,l,n;
    try
    {
        n = mtp->nC*mtp->Nf;
        z1 = new double[n];
        z2 = new double[n];
        dz = new double[n];
        memset( z1, 0, sizeof(double)*n );
        memset( z2, 0, sizeof(double)*n );
        memset( dz, 0, sizeof(double)*n );

        hj = h / (double)nj[ j ];
        // 
        for( i=0; i < n; i++ )
        {
            z1[ i ] = x[ i ];
            z2[ i ] = x[ i ] + hj * dx[ i ];
        }
        // 
        m = nj[ j ] - 1;
        for( mm=0; mm < m; mm++ )
        {
            Solut(  z2, dz, t+hj*(double)(mm+1) );
            for( i=0; i < n; i++ )
            {
                ys = z1[ i ];
                z1[ i ] = z2[ i ];
                z2[ i ] = ys + hj * dz[ i ] * 2.;
            }
        }
        // 
        Solut(  z2, dz, t+h );
        for( i=0; i < n; i++ )
            tt[ i ][ j ] = ( z1[ i ] + z2[ i ] + hj * dz[ i ] ) / 2.;
        nfcn += nj[ j ];
        //
        if(j == 0)
            goto VEL;
        for( l=j; l >= 1; l-- )
        {
            fac = pow( (double)nj[ j ] / (double)nj[ l-1 ], 2.) - 1.;
            for( i=0; i < n; i++ )
                tt[ i ][ l-1 ] = tt[ i ][ l ] +
                                 ( tt[ i ][ l ] - tt[ i ][ l-1 ] ) / fac;
        }
        err = 0e0;
        for( i=0; i < n; i++ )
        { // 
            v2 = fabs( tt[ i ][ 0 ] );
            v1 = fabs( x[ i ] );
            v1 = max( v1, v2 );
            v2 = max( 1.e-6, UROUND / epsd4 );
            scal = max( v1, v2 );
            err += pow( (tt[ i ][ 0 ] - tt[ i ][ 1 ] ) / scal, 2. );
        }
        err = pow( err / (double)n, .5 );
        //
        expo = 1.e0 / (double)( 2 * ( j + 1 )-1 );
        facmin = pow( FAC1, expo );
        v1 = max( facmin, pow( err/epsd4, expo ) / SAFE2 );
        fac = min( FAC2/facmin, v1 );
        fac= 1.e0 / fac;
        hh[ j ] = min( h * fac, MAXSTEP );
        w[ j ] = a1[ j ] / hh[ j ];

VEL:
        delete[] z1;
        delete[] z2;
        delete[] dz;
    }
    catch( ... )
    {
        if( z1 ) delete[] z1;
        if( z2 ) delete[] z2;
        if( dz ) delete[] dz;
        Error( "INTEG ", "Error in MIDEX!");
    }
}


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//  Subroutine INTEG - numerical integration of the system of n ordinary 
//    differential equations of the form 
//		 dxi/dt = f( x1,x2, ... ,xn,t )  i=1,...,n
//               xi( t0 )=xi0
//  over the time interval [t_begin, t_end] with minimal time step length 
//     step and precision eps 
//  Uses implicit method of rational extrapolation (??)
// 
void
TGEM2MT::INTEG( double eps, double& step, double t_begin, double t_end )
{
    double  t, h1, h, v1;
    int     j,i,reject,last,kc,kopt,k;

    x = mtp->MB;
    dx = mtp->dMB;
    tv = &mtp->cTau;

    h = step;
    epsd4 = eps * SAFE1;
    v1 = -log10(eps)*.6+1.5;
    v1 = min( 8.0, v1 );
    k = (int)max( 3.0, v1 ) - 1;
    t = t_begin;
    h1 = t_end-t;
    v1 = min( MAXSTEP, h1/2. );
    h = min( h, v1 );
    //CalcNewStates( 0, k, t, h ); // 14/12/2007 ????? may be done before in calc
    err = w[ 0 ] = 0.0;
    reject = last = 0;   // false

    //
    while( fabs( h1 ) >= UROUND )
    {
        v1 = min( h1, MAXSTEP);
        h = min( h, v1 );
        if( h >= ( h1 - UROUND ) )  last = 1;      // true
        Solut(  x, dx, t );
        nfcn++;
        if (( nstep == 0 )||( last ))     // 1
        {
            nstep++;
            for( j=0; j <= k; j++ )
            {
                kc=j;
                MIDEX( j, t, h );
                if( ( j > 0 ) && ( err <= eps ) ) goto l60;
            }
            goto l55;
        }
        //
l30:
        nstep++;
        ErrorIf( nstep >= mtp->ntM /*MaxIter*/, "INTEG", "No convergence - too many iterations" ); // 14/12/2007 !!!!
        kc = k-1;
        for( j=0; j <= kc; j++ )
            MIDEX( j, t, h);
        // 
        if( !( k == 1 || reject ) )
        {
            if( err <= eps ) goto l60;
            if( ( err/eps ) >
                    ( pow( (double)(nj[k+1]*nj[k])/4., 2. ) ) )  goto l100;
        }
        MIDEX( k, t, h );
        kc = k;
        if( err <= eps ) goto l60;
        // 
l55:
        if( err/eps > pow( (double)(nj[k])/2.0, 2.0 ) ) goto l100;
        kc = k + 1;
        MIDEX( kc, t, h );
        if( err > eps ) goto l100;
        // 
l60:
        t += h;
        step = h;
//        mtp->cTau = t;
        for( i=0; i < mtp->nC*mtp->Nf; i++ )
            x[i] = tt[i][0];
        Solut(  x, dx, t );
        naccept++;
        CalcNewStates( naccept, kc, t, h );
        //
        if( kc == 1 )
        {
            kopt = 2;
            if( reject ) kopt = 1;
        }
        else if( kc <= k )
        {
            kopt = kc;
            if( w[kc-1] < w[kc]*FAC3 ) kopt = kc - 1;
            if( w[kc] < w[kc-1]*FAC3 )
                kopt = min( (kc+1) , (KM-1) );
        }
        else
        {
            kopt = kc-1;
            if( kc > 2 && w[kc-2] < w[kc-1]*FAC3 )
                kopt = kc - 2;
            if( w[kc] < w[kopt]*FAC3 )
                kopt = min( kc, (KM-1) );
        }
        // 
        if( reject )
        {
            k = min( kopt, kc );
            h = min( h, hh[ k ] );
            reject = 0;           // false
        }
        else
        {  // 
            if( kopt <= kc )
                h = hh[ kopt ];
            else if( kc < k && ( w[kc] < ( w[kc-1] * FAC4 ) ))
                h = hh[kc] * a1[kopt+1] / a1[kc];
            else h = hh[kc] * a1[kopt] / a1[kc];
            k = kopt;
        }
        h1=t_end-t;
    }     // end while
    *tv = t;
    step = h;
    return;
l100:
    k = min( k, kc );
    if( k > 1 && ( w[k-1] < w[k] * FAC3 ) )
        k--;
    nrejct++;
    h = hh[k];
    reject = 1;
    goto l30;
}
//====================================================================
extern bool _comment;

void TGEM2MT::to_text_file( fstream& ff, bool with_comments )
{
  _comment = with_comments;
  
  TPrintArrays  prar(ff);

   if( _comment )
   {  ff << "# GEMIPM2K v. 2.2.3" << endl;
      ff << "# Prototype 31.03.2008" << endl;
      ff << "# Comments can be marked with # $ ;" << endl << endl;
      ff << "# Template for the Gem2mt data file" << endl;
      ff << "# (should be read before the DATACH, the IPM-DAT and DATABR files)" << endl << endl;
      ff << "#Section (scalar): Controls and dimensionalities of the Gem2mt operation" << endl;
   }
   if( _comment )
      ff << "# Code of GEM2MT mode of operation { S F A D T V W }" << endl;
   ff << left << setw(17) << "<Mode> " << "\'"<<  mtp->PsMode << "\'"<< endl;
   if( _comment )
        ff << "# number of nodes nC (1D mass-transport or box-flux models only)" << endl;
   ff << left << setw(7) << "<Size> " <<   mtp->nC /* << " 1" << " 1" */ << endl;
   if( _comment )
       ff << "# Maximum allowed number of time iteration steps (default 1000)" << endl;
   ff << left << setw(7) << "<MaxSteps> " <<   mtp->ntM << endl;
   if( _comment )
        ff << "# Physical time iterator";
   prar.writeArray(  "Tau", mtp->Tau, 3 );
   ff << endl;
   
   if( mtp->PsMode == GMT_MODE_F )
   {  if( _comment )
       ff << "# Use phase groups definitions, flux definition list & source fluxes and elemental stoichiometries";
      prar.writeArray(  "PvFDL", &mtp->PvPGD, 3, 1 );
      if( _comment )
        ff << "\n#  number of MPG flux definitions (0 or >1)" <<  endl;
      ff << left << setw(7) << "<nFD> " <<   mtp->nFD << endl;
      if( _comment )
        ff << "# number of mobile phase groups (0 or >1)" << endl;
      ff << left << setw(7) << "<nPG> " << mtp->nPG <<  endl;
      if( _comment )
        ff << "# number of source flux definitions (0 or < nFD )" << endl;
      ff << left << setw(7) << "<nSFD> " << mtp->nSFD <<  endl;
      if( _comment )
        ff << "# nICb number of ICs in  (DATABR) for setting box-fluxes" << endl;
      ff << left << setw(7) << "<Nf> " << mtp->Nf <<  endl;
      if( _comment )
        ff << "# nPHb number of phases in (DATABR) for setting box-fluxes" << endl;
      ff << left << setw(7) << "<FIf> " << mtp->FIf <<  endl;
   }     

   if( mtp->PsMode == GMT_MODE_W || mtp->PsMode == GMT_MODE_V )
   {  if( _comment )
       ff << "# Use array of grid point locations? (+ -)" << endl;
       ff << left << setw(17) << "<PvGrid> " << "\'"<<  mtp->PvGrid << "\'"<< endl;
       if( _comment )
        ff << "# res Number of allocated particle types (< 20 ? )" << endl;
       ff << left << setw(7) << "<nPTypes> " << mtp->nSFD <<  endl;
       if( _comment )
         ff << "# res Number of particle statistic properties (for monitoring) >= anPTypes" << endl;
       ff << left << setw(7) << "<nProps> " << mtp->nProps <<  endl;
       if( _comment )
         ff << "# spatial dimensions of the medium defines topology of nodes";
      prar.writeArray(  "LSize", mtp->sizeLc, 3 );
      ff << endl;
   }     

   if( _comment )
    ff << "# Advection/diffusion mass transport: initial fluid advection velocity (m/sec)" << endl;
   ff << left << setw(7) << "<fVel> " << mtp->fVel <<  endl;
   if( _comment )
    ff << "# column length (m)" << endl;
   ff << left << setw(7) << "<cLen> " << mtp->cLen <<  endl;
   if( _comment )
    ff << "# time step reduction factor" << endl;
   ff << left << setw(7) << "<tf> " << mtp->tf <<  endl;
   if( _comment )
    ff << "# cutoff factor for differences (1e-9)" << endl;
   ff << left << setw(7) << "<cdv> " << mtp->cdv <<  endl;
   if( _comment )
    ff << "# cutoff factor for minimal amounts of IC in node bulk compositions (1e-12)" << endl;
   ff << left << setw(7) << "<cez> " << mtp->cez <<  endl;
   if( _comment )
    ff << "# initial value of longitudinal dispersivity (m), usually 1e-3" << endl;
   ff << left << setw(7) << "<al_in> " << mtp->al_in <<  endl;
   if( _comment )
    ff << "# initial general diffusivity (m2/sec), usually 1e-9" << endl;
   ff << left << setw(7) << "<Dif_in> " << mtp->Dif_in <<  endl;
        
   ff<< "\n<END_DIM>\n";

 // dynamic arrays - must follow static data
   if( _comment )
   {   ff << "\n## Task configuration section ";
       ff << "\n#  Array of indexes of initial system variants for distributing to nodes";
   }
   prar.writeArray(  "DiCp", mtp->DiCp[0], mtp->nC*2, 2);
   if( _comment )
    ff << "\n# Hydraulic parameters for nodes in mass transport model";
   prar.writeArray(  "HydP", mtp->HydP[0], mtp->nC*SIZE_HYDP, SIZE_HYDP);
       
   if( mtp->PsMode == GMT_MODE_W || mtp->PsMode == GMT_MODE_V )
   {  if( _comment )
       ff << "\n# Array of initial mean particle type numbers per node";
       prar.writeArray(  "NPmean", mtp->NPmean, mtp->nPTypes );
       if( _comment )
        ff << "\n# Minimum average total number of particles of each type per one node";
       prar.writeArray(  "nPmin", mtp->nPmin, mtp->nPTypes );
        if( _comment )
         ff << "\n# Maximum average total number of particles of each type per one node";
       prar.writeArray(  "nPmax", mtp->nPmax, mtp->nPTypes );
       if( _comment )
         ff << "\n# Array of particle type definitions at t0 or after interruption";
       prar.writeArray(  "ParTD", mtp->ParTD[0], mtp->nPTypes*6, 6 );
      if( mtp->PvGrid != S_OFF )
      {
          if( _comment )
            ff << "\n# Array of grid point locations";
          prar.writeArray(  "grid", mtp->grid[0], mtp->nC*3, 3 );
      }
   }     
   if( mtp->PsMode == GMT_MODE_F )
   {  
    if( mtp->PvFDL != S_OFF )
    {
 	   if( _comment )
        ff << "\n# Indexes of nodes where this flux begins and ends";
        prar.writeArray(  "FDLi", mtp->FDLi[0], mtp->nFD*2,2 );
 	   if( _comment )
        ff << "\n# Part of the flux defnition list (flux order, flux rate, MPG quantities)";
        prar.writeArray(  "FDLf", mtp->FDLf[0], mtp->nFD*4, 4 );
 	   if( _comment )
        ff << "\n# ID of fluxes";
        prar.writeArray(  "FDLid", mtp->FDLid[0], mtp->nFD, MAXSYMB );
 	   if( _comment )
        ff << "\n# Operation codes (letters) flux type codes";
        prar.writeArray(  "FDLop", mtp->FDLop[0],  mtp->nFD, MAXSYMB  );
 	   if( _comment )
        ff << "\n# ID of MPG to move in this flux";
        prar.writeArray(  "FDLmp", mtp->FDLmp[0], mtp->nFD, MAXSYMB  );
    }
	if( mtp->PvPGD != S_OFF )
	{
	  if( _comment )
        ff << "\n# Units for setting phase quantities in MPG";
      prar.writeArray(  "UMPG", mtp->UMPG, mtp->FIf, 1 );
	  if( _comment )
         ff << "\n# Quantities of phases in MPG ";
	  prar.writeArray(  "PGT", mtp->PGT, mtp->FIf*mtp->nPG, mtp->nPG );
	  if( _comment )
	     ff << "\n# ID list of mobile phase groups";
	  prar.writeArray(  "MPGid", mtp->MPGid[0], mtp->nPG, MAXSYMB );
	}
   if( mtp->PvSFL != S_OFF )
   {
	  if( _comment )
	     ff << "\n# Table of bulk compositions of source fluxes";
      prar.writeArray(  "BSF", mtp->BSF, mtp->nSFD*mtp->Nf, mtp->Nf );
   }
 }     
   if( _comment )
     ff << "\n\n# End of file"<< endl;
 }

// --------------------- end of m_gem2mtbox.cpp ---------------------------


