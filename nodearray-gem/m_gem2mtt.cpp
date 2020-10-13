//-------------------------------------------------------------------
// $Id: m_gem2mtt.cpp 663 2005-12-13 16:27:14Z gems $
//
// Implementation of TGEM2MT class, mass transport functions
//
// Copyright (C) 2005,2007  S.Dmytriyeva, D.Kulik
//
// This file is part of a GEM-Selektor library for thermodynamic
// modelling by Gibbs energy minimization
// Uses: GEM-Selektor GUI GUI DBMS library, gems/lib/gemvizor.lib
//
// This file may be distributed under the GPL v.3 license

//
// See http://gems.web.psi.ch/ for more information
// E-mail: gems2.support@psi.ch
//-------------------------------------------------------------------
//

#include <cmath>
#include <cstdio>
#ifdef useOMP
#include <omp.h>
#endif

#ifndef IPMGEMPLUGIN

#include "m_gem2mt.h"
#include "visor.h"
#include "stepwise.h"

#else

#include <ctime>
#include "m_gem2mt.h"
#include "nodearray.h"

#endif

// ===========================================================

// Copying nodes from C1 to C0 row
void  TGEM2MT::copyNodeArrays()
{
  //  Getting direct access to TNodeArray class data
  DATACH* CH = na->pCSD();       // DataCH structure
  DATABRPTR* C0 = na->pNodT0();  // nodes at current time point
  DATABRPTR* C1 = na->pNodT1();  // nodes at previous time point
  double dc;

  for (long int ii=1; ii<mtp->nC; ii++)    // node iteration
   {
     bool NeedCopy = false;
     for(long int ic=0; ic < CH->nICb-1; ic++) // do not check charge
     {
        dc = C0[ii]->bIC[ic] - C1[ii]->bIC[ic];
        if( fabs( dc ) > std::min( mtp->cdv, (C1[ii]->bIC[ic] * 1e-3)))
                  NeedCopy = true;
     }
     if( NeedCopy )
       na->CopyNodeFromTo( ii, mtp->nC, C1, C0 );
   }  // ii    end of node iteration loop
}


// put HydP
void  TGEM2MT::putHydP( DATABRPTR* C0 )
{
    for( long int jj=0; jj<mtp->nC; jj ++)
    {
        C0[jj]->NodeTypeHY = mtp->DiCp[jj][1];
#ifdef NODEARRAYLEVEL
        if( mtp->HydP )
        { C0[jj]->Vt = mtp->HydP[jj][0];
            C0[jj]->vp = mtp->HydP[jj][1];
            C0[jj]->eps = mtp->HydP[jj][2];
            C0[jj]->Km = mtp->HydP[jj][3];
            C0[jj]->al = mtp->HydP[jj][4];
            C0[jj]->Dif = mtp->HydP[jj][5];
            C0[jj]->hDl = C0[jj]->al*C0[jj]->vp+C0[jj]->Dif;
            C0[jj]->nto = mtp->HydP[jj][6];
        }
#endif
    }
}

#ifndef IPMGEMPLUGIN


//-------------------------------------------------------------------
// NewNodeArray()  makes work DATACH structure
// Allocates a new RMT node array of DATABR structures and
// reads in the data from MULTI (using nC, Bn and DiCp, DDc, HydP)
//
void  TGEM2MT::NewNodeArray()
{
    long int nit;

    // generate Tval&Pval arrays
    if( mtp->PsTPai != S_OFF )
        gen_TPval();

    na->InitCalcNodeStructures( mtp->nICb, mtp->nDCb,  mtp->nPHb,
                                mtp->xIC, mtp->xDC, mtp->xPH, mtp->PsTPpath != S_OFF,
                                mtp->Tval, mtp->Pval,
                                mtp->nTai,  mtp->nPai, mtp->Tai[3], mtp->Pai[3]  );
    DATACH* data_CH = na->pCSD();

    // put DDc
    if( data_CH->DD && mtp->DDc )
        for( long int jj=0; jj<data_CH->nDCs; jj ++)
            data_CH->DD[jj*data_CH->nPp*data_CH->nTp] = mtp->DDc[jj];

    // Distribute data from initial systems into nodes as prescribed in DiCp[*][0]

    for( mtp->kv = 0, nit=0; mtp->kv < mtp->nIV; mtp->kv++ )
    {
        // Make new SysEq record
        gen_task( true );
        // calculate EqStat (Thermodynamic&Equlibria)
        //  TProfil::pm->pmp->pTPD = 0;
        calc_eqstat( true );

        for( long int q=0; q<mtp->nC; q++)
        {
            mtp->qc = q;

            pVisor->Message( window(), GetName(),
                             "Initial calculation of equilibria in nodes. "
                             "Please, wait...", nit, mtp->nC);

            if( mtp->DiCp[q][0] == mtp->kv  )
            {
                nit++;
                gen_task( false );
                // set T, P for multi and calculate current EqStat
                calc_eqstat( false );

                if( mtp->PsMode == RMT_MODE_S )
                {
                    // empty current node
                    DATABR* data_BR = na->reallocDBR( q, mtp->nC,  na->pNodT0());
                    data_BR->TK = TMulti::sm->GetPM()->TCc+C_to_K; //25
                    data_BR->P = TMulti::sm->GetPM()->Pc*bar_to_Pa; //1
                    for(long int i1=0; i1<mtp->nICb; i1++ )
                        data_BR->bIC[i1] = TMulti::sm->GetPM()->B[ mtp->xIC[i1] ];
                }
                else // Save databr
                {
                    na->SaveToNode( q, mtp->nC,  na->pNodT0());
                }
                na->pNodT0()[q]->NodeHandle = q;
            }
        } // q
        mt_next();      // Generate work values for the next EqStat rkey
    } // mtp->kv

    pVisor->CloseMessage();

    DATABRPTR* C0 = na->pNodT0();  // nodes at current time point
    for( long int ii=0; ii<mtp->nC; ii++)
        if(  C0[ii] == nullptr )
            Error( "NewNodeArray() error:" ," Undefined boundary condition!" );

    // put HydP
    putHydP( C0 );

    for (long int ii=0; ii<mtp->nC; ii++)    // node iteration
    {
        na->CopyNodeFromTo( ii, mtp->nC, C0, na->pNodT1() );
    }  // ii    end of node iteration loop

    /* test output generated structure
    ProcessProgressFunction messageF = [](const std::string& , long ){
        //std::cout << "TProcess GEM3k output" <<  message.c_str() << point << std::endl;
        return false;
    };

    auto dbr_list =  na->genGEMS3KInputFiles(  "Te_start/Test-dat.lst", messageF, mtp->nC, 0, false, false, false, false );
    */
}

#endif

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// This function (for 1D-case only) analyses the node array
// to decide where to use AIA and where PIA for GEM calculations
// in each node.
// if IAmode == NEED_GEM_AIA then all nodes are set to AIA
// If IAmode == NEED_GEM_PIA then the adjacent nodes with different
// reactive phase assemblage (very different pe, pH), as well as their
// neighbors from each side, are set to AIA. Other nodes are set to PIA.
// returns: Number of nodes set to AIA
//
long int TGEM2MT::CheckPIAinNodes1D( char IAmode, long int start_node, long int end_node )
{
       long int nSetAIA = 0;
       // Getting direct access to data
       DATACH* CH = na->pCSD();       // DataCH structure
       DATABRPTR* C0 = na->pNodT0();  // nodes at previous time point
//       DATABRPTR* C1 = na->pNodT1();  // nodes at current time point
       bool* iaN = na->piaNode();      // indicators for IA in the nodes

       start_node =std::max( start_node, 0L );
       end_node = std::min( end_node, mtp->nC-1 );

       // Initializing iaNode vector
       if( IAmode == NEED_GEM_SIA && CH->nDCb == CH->nDC )
       {
          for( long int ii = start_node; ii<= end_node; ii++ )
    		 iaN[ii] = false;    // potentially need PIA
       }
       else { // setting all nodes to AIA GEM calculations
          for( long int ii = 0; ii< mtp->nC; ii++ )
     	  {
     		  iaN[ii] = true;
     		  nSetAIA++;
     	  }
          return nSetAIA;
       }

       // This is done only if PIA mode is requested!
       // First variant: 1-neighbour algorithm
       for( long int ii = start_node; ii<= end_node; ii++) // node iteration
	   {
         // mtp->qc = ii;
	      switch( C0[ii]->NodeTypeMT )
	     {
	        case normal: // normal node
	             // Checking phase assemblage
	        	break;
	                     // boundary condition node
	        case NBC1source:  //1: Dirichlet source ( constant concentration )
	        	    break;
	        case NBC1sink:    // -1: Dirichlet sink
	        case NBC2source:  // 2: Neumann source ( constant gradient )
	        case NBC2sink:    // -2: Neumann sink
	        case NBC3source:  // 3: Cauchy source ( constant flux )
	        case NBC3sink:    // -3: Cauchy sink
	        case INIT_FUNK:   // 4: functional conditions (e.g. input time-depended functions)
	        default:
	        			iaN[ii] = true;
	        			nSetAIA++;
	        			continue;
	     }
         if( (ii == 0 || ii == mtp->nC-1 ) ) // &&
//	    	 (C0[ii]->NodeTypeMT == normal || C0[ii]->NodeTypeMT == NBC1source) )
//	     {
//	     	iaN[ii] = true;
//	     	nSetAIA++;
	     	continue;
//	     }
	     // checking pair of adjacent nodes for phase assemblage differences
         for( long int kk=0; kk<CH->nPHb; kk++ )
	     {
            if( approximatelyZero( C0[ii]->xPH[kk] ) && approximatelyZero( C0[ii-1]->xPH[kk] ) )
	    		continue;    // Phase kk is absent in both node systems
	    	if( C0[ii]->xPH[kk] > 0.0 && C0[ii-1]->xPH[kk] > 0.0 )
	    		continue;    // Phase kk is present in both node systems
	    	// At least one phase is absent in one and present in another node system
	    	goto DIFFERENT;
	     }
	     // Should we check difference in pe or pH?

	     continue;
DIFFERENT:
         if( iaN[ii] == false )
         {
	         iaN[ii] = true;
	         nSetAIA++;
         }
         if( iaN[ii-1] == false )
         {
	         iaN[ii-1] = true;
	         nSetAIA++;
         }
// Activating first neighbours - experimental algorithm
         if( iaN[ii+1] == false )
         {
	         iaN[ii+1] = true;
	         nSetAIA++;
         }
         if( ii > start_node+1 && iaN[ii-2] == false )
         {
	         iaN[ii-2] = true;
	         nSetAIA++;
         }
	   } // loop ii

      // End of checking
       return nSetAIA;
}

//   Here we call a loop on GEM calculations over nodes
//   parallelization should affect this loop only
//   mode can be NEED_GEM_PIA (smart algorithm) or NEED_GEM_AIA
//   (forced AIA, usually used at the beginning of coupled run)
//   return code   true   Ok
//                 false  Error in GEMipm calculation part
//
bool TGEM2MT::CalcIPM( char mode, long int start_node, long int end_node, FILE* updiffile )
{
    bool iRet = true;

    FILE* diffile = nullptr;
    if( mtp->PsMO != S_OFF )
      diffile = updiffile;

    start_node = std::max( start_node, 0L );
    end_node = std::min( end_node, mtp->nC-1 );

    for( long int ii = start_node; ii<= end_node; ii++) // node iteration
    {
      node1_Tm( ii ) = mtp->cTau;     // set time and time step (for TKinMet)
      node1_dt( ii ) = mtp->dTau;
   }

#ifdef useOMP
   double  t0 = omp_get_wtime();
#else
   clock_t t_start, t_end;
   t_start = clock();
#endif


   na->CalcIPM_List( TestModeGEMParam(mode, mtp->PsSIA, mtp->ct, mtp->cdv, mtp->cez ), start_node, end_node, diffile );

   /* test output generated structure
   ProcessProgressFunction messageF = [](const std::string& , long ){
       //std::cout << "TProcess GEM3k output" <<  message.c_str() << point << std::endl;
       return false;
   };

   auto dbr_list =  na->genGEMS3KInputFiles(  "Te_point/Test-dat.lst", messageF, mtp->nC, 0, false, false, true, false );
   */

#ifdef useOMP
   double  t1 = omp_get_wtime();
   mtp->TimeGEM += t1-t0;
#else
   double clc_sec = CLOCKS_PER_SEC;
   t_end = clock();
   mtp->TimeGEM = ( t_end- t_start )/clc_sec;
#endif

   // Here dt can be analyzed over nodes - if changed anywhere by TKinMet
   //   then the overall time step should be reduced and the MT loop started all over

   mtp->qc = end_node;
   return iRet;
}

// The mass transport iteration initial time step
void TGEM2MT::MassTransAdvecStart()
{
    if( mtp->tf <= 0 ) // foolproof
        mtp->tf = 1.;
    mtp->dx = mtp->sizeLc[0]/mtp->nC;  // was mtp->cLen/mtp->nC; changed 15.12.2011 DK
    mtp->dTau = 0.5*(mtp->dx/mtp->fVel)*1./mtp->tf;
    mtp->oTau = 0.;
    mtp->cTau = mtp->Tau[START_];
    // mtp->cTau = 0;
    mtp->ct = 0;
}

// The mass transport iteration initial time step (Crank-Nicolson scheme)
void TGEM2MT::MassTransCraNicStart()
{
    if( mtp->tf <= 0 ) // foolproof
        mtp->tf = 1.;
    mtp->dx = mtp->sizeLc[0]/mtp->nC;  // was mtp->cLen/mtp->nC; changed 15.12.2011 DK
    mtp->dTau = (mtp->dx/mtp->fVel)/mtp->tf;
    mtp->oTau = 0.;
    mtp->cTau = mtp->Tau[START_];
    // mtp->cTau = 0;
    mtp->ct = 0;
}


// The mass transport iteration initial time step (random-walk scheme)
void TGEM2MT::MassTransParticleStart()
{
    double dt_adv=0.0, dt_dif=0.0;
    if( mtp->tf <= 0 )
        mtp->tf = 1.;
    mtp->dx = mtp->sizeLc[0]/mtp->nC;  // was mtp->cLen/mtp->nC; changed 15.12.2011 DK
    if(fabs(mtp->fVel) > 1e-19)
        dt_adv = (mtp->dx / mtp->fVel);
    // Courant - Neumann criteria
    dt_dif = (mtp->dx*mtp->dx)/2./(mtp->al_in*mtp->fVel + mtp->eps_in*mtp->Dif_in/mtp->nto_in); // bugfix 10.10.19 DK
    if(fabs(dt_adv) > 1e-19)
        mtp->dTau = std::min( dt_adv, dt_dif ) / mtp->tf; // advection and diffusion
    else
        mtp->dTau = dt_dif / mtp->tf;    // only diffusion
//    mtp->dTau = 0.5*(mtp->dx/mtp->fVel)*1./mtp->tf;  // Courant criterion
    mtp->oTau = 0.;
    mtp->cTau = mtp->Tau[START_];
    // mtp->cTau = 0;
    mtp->ct = 0;
    pa_mt->setUpCounters();

}


// The mass transport iteration time step
void TGEM2MT::MassTransParticleStep( bool CompMode )
{
   mtp->ct += 1;
   mtp->oTau = mtp->cTau;
   mtp->cTau += mtp->dTau;

   pa_mt->GEMPARTRACK( mtp->PsMode, CompMode, mtp->oTau, mtp->cTau );
}


// The mass transport iteration time step
// CompMode: true: Transport via dependent components in Aq phase (except H2O)
//           false: Transport via dissolved indepedent components (with H2O)
//
void TGEM2MT::MassTransAdvecStep( bool CompMode )
{
 double c0, c1, cm1,  cm2, /*cmax,*/ c12, cm12,
        /*charge, c0new,*/ dc, cr, aji, fmolal;      // some help variables
 long int ii, ic, jc;

 //  Getting direct access to TNodeArray class data
 DATACH* CH = na->pCSD();       // DataCH structure
 // DATABRPTR* C0 = na->pNodT0();  // nodes at current time point
 // DATABRPTR* C1 = na->pNodT1();  // nodes at previous time point

   mtp->ct += 1;
   mtp->cTau += mtp->dTau;
   cr = mtp->fVel * mtp->dTau/mtp->dx;

   for( ii = 2; ii< mtp->nC-1; ii++) // node iteration, -1 the right boundary is open ....
   {
     mtp->qc = ii;
     if( CompMode == true )
     {  // Advective mass transport over DC gradients in Aq phase
         fmolal = 55.5084/node1_xDC( ii, CH->nDCinPH[0]-1 );
    	 for( jc=0; jc < CH->nDCinPH[0]-1; jc++ )
     	 {   // splitting for dependent components except H2O@
     		 // It has to be checked on minimal allowed c0 value
    		 c0  = node1_xDC( ii, jc )* fmolal;
     		 c1  = node1_xDC( ii+1, jc )* fmolal;
     		 cm1 = node1_xDC( ii-1, jc )* fmolal;
     		 cm2 = node1_xDC( ii-2, jc )* fmolal;
       		 if( c0 < 1e-20 && c1 < 1e-20 && cm1 < 1e-20 && cm2 < 1e-20 )
      			continue;
     		 c12=((c1+c0)/2)-(cr*(c1-c0)/2)-((1-cr*cr)*(c1-2*c0+cm1)/6);
    	     cm12=((c0+cm1)/2)-(cr*(c0-cm1)/2)-((1-cr*cr)*(c0-2*cm1+cm2)/6);
    	     dc = cr*(c12-cm12);
             if( fabs( dc ) < mtp->cdv )   // *c0   Insignificant fractional difference?
            	 continue;
    	     dc /= fmolal;
             // Checking the new DC amount
    	     if( (c0/fmolal - dc) > mtp->cez )
    	     {  // the amount of DC remains positive 
    	    	 node1_xDC( ii, jc ) -= dc;			 
            	 for( ic=0; ic<CH->nICb; ic++)  // incrementing independent components
            	 {
                     aji = na->DCaJI( jc, ic );
                     if( noZero( aji ) )
            		     node1_bIC(ii, ic) -= aji * dc;
            	 }
    	     }
    	     else {  // setting the DC amount to minimum
    	    	 node1_xDC( ii, jc ) = mtp->cez;
    	    	 dc = c0/fmolal - mtp->cez;
    	    	 for( ic=0; ic<CH->nICb; ic++)  // incrementing independent components
            	 {
                     aji = na->DCaJI( jc, ic );
                     if( noZero( aji ) )
            		     node1_bIC(ii, ic) -= aji * dc;
            	 }
    	     }	 
    	 } // loop over DC
     }
     else { // Advective mass transport over IC gradients in Aq phase
    	 double niw; // IC amount in H2O
    	 fmolal = 55.5084/node1_xDC( ii, CH->nDCinPH[0]-1 ); // molality factor
    	 for( ic=0; ic < CH->nICb; ic++)  // splitting for independent components
    	 {                        
                 niw = node1_xDC( ii, CH->nDCinPH[0]-1 )* na->DCaJI( CH->nDCinPH[0]-1, ic );
    		    // IC amount in H2O
    		 // Chemical compositions may become inconsistent with time
    		 // It has to be checked on minimal allowed c0 value
    		 c0  = (node1_bPS( ii, 0, ic ) - niw )*fmolal;    //C1[ii]->bPS[0*CH->nICb + ic];
             c1  = (node1_bPS( ii+1, 0, ic) - niw )*fmolal;   //C1[ii+1]->bPS[0*CH->nICb + ic];
    		 cm1 = (node1_bPS( ii-1, 0, ic ) - niw )*fmolal;  //C1[ii-1]->bPS[0*CH->nICb + ic];
             cm2 = (node1_bPS( ii-2, 0, ic ) - niw )*fmolal;  //C1[ii-2]->bPS[0*CH->nICb + ic];

    		 // Finite-difference calculation (suggested by FE )
             c12=((c1+c0)/2)-(cr*(c1-c0)/2)-((1-cr*cr)*(c1-2*c0+cm1)/6);
             cm12=((c0+cm1)/2)-(cr*(c0-cm1)/2)-((1-cr*cr)*(c0-2*cm1+cm2)/6);
             dc = cr*(c12-cm12);
             if( fabs( dc ) < mtp->cdv )  // *c0
                    continue;
    		 dc /= fmolal; 
             // Checking the new IC amount
             if( (node1_bPS(ii, 0, ic) - dc) > mtp->cez )
    		 {  // New IC amount is positive
    			 node1_bPS( ii, 0, ic ) -= dc;
    			 node1_bIC(ii, ic) -= dc;
    		 }
    		 else { // Setting the new IC amount to threshold  
    			 node1_bPS( ii, 0, ic ) = mtp->cez;
    			 node1_bIC(ii, ic) = mtp->cez;
    		 }
    		 /*if( dc >= C1[i]->bIC[ic] )
 			 {
    			fprintf( diffile, "\nError in Mass Transport calculation part :" );
    			fprintf( diffile, " Node= %-8d  Step= %-8d  IC= %s ", i, t, CH->ICNL[ic] );
    			fprintf(diffile, "\n Attempt to set new amount to %lg (old: %lg  Diff: = %lg ) ",
         		C1[i]->bIC[ic]-dc, C1[i]->bIC[ic], dc);
    			BC_error = true;
 			 } */
    	 } // loop over IC
     }
	 // checking charge balance
     //charge = node1_bIC(ii, CH->nICb-1 );
     node1_bIC(ii, CH->nICb-1 ) = 0.0;		// debugging
   } // end of loop over nodes
}

// The mass transport iteration time step (Crank-Nicolson scheme, cf Alina Yapparova 2015)
// CompMode: true: Transport via dependent components in Aq phase (except H2O)
//           false: Transport via dissolved indepedent components (with H2O)
//
void TGEM2MT::MassTransCraNicStep( bool CompMode )
{
 double c0, c1, cm1,  cm2, /*cmax,*/ c12, cm12,
        /*charge, c0new,*/ dc, cr, aji, fmolal;      // some help variables
 long int ii, ic, jc;

 //  Getting direct access to TNodeArray class data
 DATACH* CH = na->pCSD();       // DataCH structure
 // DATABRPTR* C0 = na->pNodT0();  // nodes at current time point
 // DATABRPTR* C1 = na->pNodT1();  // nodes at previous time point

   mtp->ct += 1;
   mtp->cTau += mtp->dTau;
   cr = mtp->fVel * mtp->dTau/mtp->dx;

   for( ii = 2; ii< mtp->nC-1; ii++) // node iteration, -1 the right boundary is open ....
   {
     mtp->qc = ii;
     if( CompMode == true )
     {  // Advective mass transport over DC gradients in Aq phase
         fmolal = 55.5084/node1_xDC( ii, CH->nDCinPH[0]-1 );
         for( jc=0; jc < CH->nDCinPH[0]-1; jc++ )
         {   // splitting for dependent components except H2O@
             // It has to be checked on minimal allowed c0 value
             c0  = node1_xDC( ii, jc )* fmolal;
             c1  = node1_xDC( ii+1, jc )* fmolal;
             cm1 = node1_xDC( ii-1, jc )* fmolal;
             cm2 = node1_xDC( ii-2, jc )* fmolal;
             if( c0 < 1e-20 && c1 < 1e-20 && cm1 < 1e-20 && cm2 < 1e-20 )
                continue;
             c12=((c1+c0)/2)-(cr*(c1-c0)/2)-((1-cr*cr)*(c1-2*c0+cm1)/6);
             cm12=((c0+cm1)/2)-(cr*(c0-cm1)/2)-((1-cr*cr)*(c0-2*cm1+cm2)/6);
             dc = cr*(c12-cm12);
             if( fabs( dc ) < mtp->cdv )   // *c0   Insignificant fractional difference?
                 continue;
             dc /= fmolal;
             // Checking the new DC amount
             if( (c0/fmolal - dc) > mtp->cez )
             {  // the amount of DC remains positive
                 node1_xDC( ii, jc ) -= dc;
                 for( ic=0; ic<CH->nICb; ic++)  // incrementing independent components
                 {
                     aji = na->DCaJI( jc, ic );
                     if( noZero( aji ) )
                         node1_bIC(ii, ic) -= aji * dc;
                 }
             }
             else {  // setting the DC amount to minimum
                 node1_xDC( ii, jc ) = mtp->cez;
                 dc = c0/fmolal - mtp->cez;
                 for( ic=0; ic<CH->nICb; ic++)  // incrementing independent components
                 {
                     aji = na->DCaJI( jc, ic );
                     if( noZero( aji ) )
                         node1_bIC(ii, ic) -= aji * dc;
                 }
             }
         } // loop over DC
     }
     else { // Advective mass transport over IC gradients in Aq phase
         double niw; // IC amount in H2O
         fmolal = 55.5084/node1_xDC( ii, CH->nDCinPH[0]-1 ); // molality factor
         for( ic=0; ic < CH->nICb; ic++)  // splitting for independent components
         {
                 niw = node1_xDC( ii, CH->nDCinPH[0]-1 )* na->DCaJI( CH->nDCinPH[0]-1, ic );
                // IC amount in H2O
             // Chemical compositions may become inconsistent with time
             // It has to be checked on minimal allowed c0 value
             c0  = (node1_bPS( ii, 0, ic ) - niw )*fmolal;    //C1[ii]->bPS[0*CH->nICb + ic];
             c1  = (node1_bPS( ii+1, 0, ic) - niw )*fmolal;   //C1[ii+1]->bPS[0*CH->nICb + ic];
             cm1 = (node1_bPS( ii-1, 0, ic ) - niw )*fmolal;  //C1[ii-1]->bPS[0*CH->nICb + ic];
             cm2 = (node1_bPS( ii-2, 0, ic ) - niw )*fmolal;  //C1[ii-2]->bPS[0*CH->nICb + ic];

             // Finite-difference calculation (suggested by FE )
             c12=((c1+c0)/2)-(cr*(c1-c0)/2)-((1-cr*cr)*(c1-2*c0+cm1)/6);
             cm12=((c0+cm1)/2)-(cr*(c0-cm1)/2)-((1-cr*cr)*(c0-2*cm1+cm2)/6);
             dc = cr*(c12-cm12);
             if( fabs( dc ) < mtp->cdv )  // *c0
                    continue;
             dc /= fmolal;
             // Checking the new IC amount
             if( (node1_bPS(ii, 0, ic) - dc) > mtp->cez )
             {  // New IC amount is positive
                 node1_bPS( ii, 0, ic ) -= dc;
                 node1_bIC(ii, ic) -= dc;
             }
             else { // Setting the new IC amount to threshold
                 node1_bPS( ii, 0, ic ) = mtp->cez;
                 node1_bIC(ii, ic) = mtp->cez;
             }
             /*if( dc >= C1[i]->bIC[ic] )
             {
                fprintf( diffile, "\nError in Mass Transport calculation part :" );
                fprintf( diffile, " Node= %-8d  Step= %-8d  IC= %s ", i, t, CH->ICNL[ic] );
                fprintf(diffile, "\n Attempt to set new amount to %lg (old: %lg  Diff: = %lg ) ",
                C1[i]->bIC[ic]-dc, C1[i]->bIC[ic], dc);
                BC_error = true;
             } */
         } // loop over IC
     }
     // checking charge balance
     //charge = node1_bIC(ii, CH->nICb-1 );
     node1_bIC(ii, CH->nICb-1 ) = 0.0;		// debugging
   } // end of loop over nodes
}




/*
// A very simple example of finite difference transport algorithm
// Contributed on 22.08.2014 by Alina Yapparova, Chair of Reservoir Engineering,
// Montanuniveritaet Leoben, Austria
//
double TMyTransport::OneTimeStepRun( long int *ICndx, long int nICndx )
{
    double column_length  = 0.5; // [m]
    double dx = column_length/(nNodes-1);

    //constant velocity field
    double v = 1.e-8; // velocity [m/s]
    // stability requirement: dt<=dx/velocity, so we can choose any coefficient k<1
    double k = 0.5; // k = dt/dx*v
    // calculate dt
    double dt = k*dx/v;
    // and print dt into the output file
    // cout<<"\tdt = " << dt << "[s]" << endl;

    //Finite difference approximation for the equation dc/dt + v*dc/dx = 0
    // explicit time, left spatial derivative
    // (c^{n+1}_{i} - c^{n}_{i})/dt + v*(c^{n}_{i} - c^{n}_{i-1})/dx = 0
    // c^{n+1}_{i} = (1 - k)*c^{n}_{i} + k*c^{n}_{i-1}
    long int in;
    for(  in=1; in< nNodes; in++ )
    {
        for (int i=0; i<nICndx; i++)
        {
            abIC[in][ICndx[i]] = (1 - k)*abPS[in][ICndx[i]] + k*abPS[in-1][ICndx[i]] + abSP[in][ICndx[i]];
            // where abPS is the total amount of an independent component ICndx[i] in the aqueous phase (after GEMS computation)
            // abSP is the total amount of an independent component ICndx[i] in ALL solid phases
            // abIC is the total amount of ICndx[i] that will serve as an input constraint for GEMS computation at the next time level
            // NB: more precisely, one should write abPS[in][n*nIC + ICndx[i]], where 0=<n<=nPS, and transport each phase-solution separately
            //but as far as an aqueous phase is the first in the list (n=0), this simplified indexing will work for transport of aq phase
        }
    }
    return dt;
}

// Finite difference transport algorithm (Crank-Nicolson scheme)
// Contributed on 27.07.2015 by Alina Yapparova, Chair of Reservoir Engineering,
// Montanuniveritaet Leoben, Austria
//
double TMyTransport::OneTimeStepRun_CN( long int *ICndx, long int nICndx )
{
    double column_length  = 0.5; // [m]
    double dx = column_length/(nNodes-1);

    //constant velocity field
    double v = 1.e-8; // velocity [m/s]
    // stability requirement: unconditionally stable
    double dt = 1000000.; //[s]
    double k = dt*v/(2*dx);

    //Finite difference approximation for the equation dc/dt + v*dc/dx = 0
    // semi-implicit time, left spatial derivative
    // (1 + k)*c^{n+1}_{i} - k*c^{n+1}_{i-1} =  (1 - k)*c^{n}_{i} + k*c^{n}_{i-1}
    // c^{n+1}_{i} = 1/(1+k)*(k*c^{n+1}_{i-1} + (1 - k)*c^{n}_{i} + k*c^{n}_{i-1})
    long int in;
    for(  in=1; in< nNodes; in++ )
    {
        for (int i=0; i<nICndx; i++)
        {
            abIC[in][ICndx[i]] = 1/(1+k)*( k*(abIC[in-1][ICndx[i]] - abSP[in-1][ICndx[i]]) + (1 - k)*abPS[in][ICndx[i]] + k*abPS[in-1][ICndx[i]] ) + abSP[in][ICndx[i]];
            // where abPS is the total amount of an independent component ICndx[i] in the aqueous phase (after GEMS computation)
            // abSP is the total amount of an independent component ICndx[i] in ALL solid phases
            // abIC is the total amount of ICndx[i] that will serve as an input constraint for GEMS computation at the next time level
            // NB: more precisely, one should write abPS[in][n*nIC + ICndx[i]], where 0=<n<=nPS, and transport each phase-solution separately
            //but as far as an aqueous phase is the first in the list (n=0), this simplified indexing will work for transport of aq phase
        }
    }
    return dt;
}

*/






// Main call for the mass transport iterations in 1D case
//
// mode is NEED_GEM_AIA or NEED_GEM_PIA (see DATABR.H) mode of GEM initial approximation
// If NEED_GEM_PIA then the program will try smart initial approximation for the nodes
// (PIA when possible, AIA in the vicinity of fronts, pH and pe barriers).
// returns true, if cancelled/interrupted by the user; false if finished Ok
bool TGEM2MT::Trans1D( char mode )
{
  bool iRet = false;
  bool CompMode = false;   // Component transport mode: true: DC; false: IC
  long int nStart = 0, nEnd = mtp->nC;
  // long int NodesSetToAIA;
// std::string Vmessage;

FILE* logfile = nullptr;
FILE* ph_file = nullptr;
FILE* diffile = nullptr;

#ifndef IPMGEMPLUGIN
   showMss = 0L;
#endif

if( mtp->PvDDc == S_ON && mtp->PvDIc == S_OFF )  // Set of DC transport using record switches
	CompMode = true; 
if( mtp->PvDDc == S_OFF && mtp->PvDIc == S_ON )  // Set of IC transport using record switches
	CompMode = false; 

if( mtp->PsMO != S_OFF )
{
    std::string fname;
#ifndef IPMGEMPLUGIN
    fname = pVisor->userGEMDir();
#endif
    // Preparations: opening output files for monitoring 1D profiles
logfile = fopen( ( fname + "ICaq-log.dat").c_str(), "w+" );    // Total dissolved element molarities
if( !logfile)
  return iRet;
ph_file = fopen( ( fname + "Ph-log.dat" ).c_str(), "w+" );   // Mole amounts of phases
if( !ph_file)
  return iRet;
diffile = fopen( ( fname + "ICdif-log.dat").c_str(), "w+" );   //  Element amount diffs for t and t-1
if( !diffile)
  return iRet;
}

// time scales testing
#ifdef useOMP
double  t0 = omp_get_wtime();
#else
clock_t t_start, t_end;
t_start = clock();
#endif
double otime = 0.;
mtp->TimeGEM = 0.0;

   if( mtp->iStat!= AS_RUN  )
   {  switch( mtp->PsMode )
     {
       case RMT_MODE_A:   // A: 1D advection (numerical) coupled FMT finite-differences model
                         MassTransAdvecStart();
                           nStart = 1; nEnd = mtp->nC-1;
                 break;
       case RMT_MODE_C:   // C: 1D advection (numerical) coupled FMT Crank-Nicolson scheme model
                         MassTransCraNicStart();
                           nStart = 1; nEnd = mtp->nC-1;
                 break;
       case RMT_MODE_W:   // W: random-walk advection-diffusion coupled FMT scoping model
                         MassTransParticleStart();
                 break;
       case RMT_MODE_F:   // F Flow-throught model ( with reactors MGP fluxes )
                         BoxFluxTransportStart();
                break;
      default: // more mass transport models here
               break;
     }
//
       // Calculation of chemical equilibria in all nodes at the beginning
       // with the LPP AIA
       if( mtp->PsSIA != S_ON )
           CalcIPM( NEED_GEM_AIA, nStart, nEnd, diffile );
       else
           CalcIPM( NEED_GEM_SIA,nStart, nEnd, diffile );
       ///       CalcIPM( NEED_GEM_AIA, nStart, nEnd, diffile );
   }
   mtp->iStat = AS_READY;

   if( mtp->PsMode == RMT_MODE_F )
      CalcMGPdata();

if( mtp->PsMO != S_OFF )
  otime += PrintPoint( 2, diffile, logfile, ph_file );

if( mtp->PsVTK != S_OFF )
   otime += PrintPoint( 0 );


#ifndef IPMGEMPLUGIN

bool UseGraphMonitoring = false;
char buf[300];

   if( mtp->PsSmode == S_OFF )
    if(  mtp->PvMSg != S_OFF && vfQuestion(window(),
             GetName(), "Use graphic monitoring?") )
        {
            RecordPlot( nullptr );
            UseGraphMonitoring = true;
        }
#endif

//  This loop contains the mass transport iteration time step
     do {   // time iteration step

#ifndef IPMGEMPLUGIN

       if( mtp->ct > 0)
            CalcStartScript();

       sprintf(buf, "   time %lg; step %ld ", mtp->cTau, mtp->ct );
       Vmessage = "Simulating Reactive Transport: ";
       Vmessage += buf;
       Vmessage += ". Please, wait (may take time)...";


       if( mtp->PsSmode != S_OFF  )
       {
         STEP_POINT2();
       }
       else
           iRet = pVisor->Message( window(), GetName(),Vmessage.c_str(),
                           mtp->ct, mtp->ntM, UseGraphMonitoring );

#endif
      if( iRet )
           break;

      //  the mass transport iteration time step
      switch( mtp->PsMode )
      {
          case RMT_MODE_A: MassTransAdvecStep( CompMode );
                           break;
          case RMT_MODE_C: MassTransCraNicStep( CompMode );
                       break;
          case RMT_MODE_W: MassTransParticleStep( CompMode );
                           break;
          case RMT_MODE_F: FlowThroughBoxFluxStep();
                           break;
//          case RMT_MODE_D:
          default: // more mass transport models here
                  break;
      }

        // The analysis of GEM IA modes in nodes - optional
//        NodesSetToAIA = CheckPIAinNodes1D( mode, nStart, nEnd );
         
        //   Here we call a loop on GEM calculations over nodes
        //   parallelization should affect this loop only
        CalcIPM( mode, nStart, nEnd, diffile );

if( mtp->PsMO != S_OFF )
  otime += PrintPoint( 3, diffile, logfile, ph_file );

          // Here one has to compare old and new equilibrium phase assemblage
          // and pH/pe in all nodes and decide if the time step was Ok or it
          // should be decreased. If so then the nodes from C0 should be
          // copied to C1 (to be implemented)

#ifndef IPMGEMPLUGIN
          // time step accepted - Copying nodes from C1 to C0 row
          pVisor->Update();
          CalcGraph();
#endif
          // copy node array for T0 into node array for T1
          copyNodeArrays();
          // copy particle array ?
          if( mtp->PsMode == RMT_MODE_W )
            pa_mt->CopyfromT1toT0();

          if( mtp->PsMode == RMT_MODE_F )
             CalcMGPdata();

if( mtp->PsMO != S_OFF )
   otime += PrintPoint( 4, diffile, logfile, ph_file );

if( mtp->PsVTK != S_OFF )
   otime += PrintPoint( 0 );

 } while ( mtp->cTau < mtp->Tau[STOP_] && mtp->ct < mtp->ntM );


#ifdef useOMP
double  t1 = omp_get_wtime();
double dtime = ( t1- t0 );
#else
double clc_sec = CLOCKS_PER_SEC;
t_end = clock();
double dtime = ( t_end- t_start )/clc_sec;
#endif


if( mtp->PsMO != S_OFF )
{
fprintf( diffile,
  "\nTotal time of calculation %lg s;  Time of output %lg s;  Whole run time %lg s;  Pure GEM run time %lg s\n",
    (dtime-otime),  otime, dtime, mtp->TimeGEM );
fclose( logfile );
fclose( ph_file );
fclose( diffile );
}

#ifndef IPMGEMPLUGIN

   pVisor->CloseMessage();


#endif

  return iRet;
}


// plotting the record -------------------------------------------------
//Added one point to graph
double TGEM2MT::PrintPoint( long int nPoint, FILE* diffile, FILE* logfile, FILE* ph_file )
{
    long int evrt =10;

#ifdef useOMP
    double  t0 = omp_get_wtime();
#else
    clock_t t_out, t_out2;
    t_out = clock();
#endif


    // from BoxEqStatesUpdate (BoxFlux ) not tested and used
    if( nPoint == 1 )
    {
       if( diffile )
       {
           na->logDiffsIC( diffile, mtp->ct, mtp->cTau, mtp->nC, 10 );
           // logging differences after the MT iteration loop
       }
   }


   // from  Trans1D
   if( nPoint == 2 )
   {
     na->logDiffsIC( diffile, mtp->ct, mtp->cTau, mtp->nC, 1 );
     na->logProfileAqIC( logfile, mtp->ct, mtp->cTau, mtp->nC, 1 );
     na->logProfilePhMol( ph_file, pa_mt, mtp->ct, mtp->cTau, mtp->nC, 1 );
   }

   if( nPoint == 3 )
   {
       na->logDiffsIC( diffile, mtp->ct, mtp->cTau, mtp->nC, evrt );
   }

   if( nPoint == 4 )
   {
       na->logProfileAqIC( logfile, mtp->ct, mtp->cTau, mtp->nC, evrt );
       na->logProfilePhMol( ph_file, pa_mt, mtp->ct, mtp->cTau, mtp->nC, evrt );
   }

   // write to VTK
   if( nPoint == 0  && mtp->PsVTK != S_OFF )
   {
       char buf[200];

       sprintf( buf, "%05ld", mtp->ct);
       std::string name = buf;

       strip(name);
       name += ".vtk";

       name = pathVTK + prefixVTK + name;

       std::fstream out_br(name.c_str(), std::ios::out );
       ErrorIf( !out_br.good() , name, "VTK text make error");
       na->databr_to_vtk(out_br, nameVTK.c_str(), mtp->cTau, mtp->ct, mtp->nVTKfld, mtp->xVTKfld );
   }

#ifdef useOMP
   double  t1 = omp_get_wtime();
   return ( t1- t0 );
#else
   double clc_sec = CLOCKS_PER_SEC;
   t_out2 = clock();
   return ( t_out2 -  t_out)/clc_sec;
#endif

}

// --------------------- end of m_gem2mtt.cpp ---------------------------

