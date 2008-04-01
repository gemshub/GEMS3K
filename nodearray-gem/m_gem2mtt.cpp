//-------------------------------------------------------------------
// $Id: m_gem2mtt.cpp 663 2005-12-13 16:27:14Z gems $
//
// Implementation of TGEM2MT class, mass transport functions
//
// Copyright (C) 2005,2007  S.Dmytriyeva, D.Kulik
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

// ===========================================================

// Copying nodes from C1 to C0 row
void  TGEM2MT::copyNodeArrays()
{
  //  Getting direct access to TNodeArray class data
  DATACH* CH = na->pCSD();       // DataCH structure
  DATABRPTR* C0 = na->pNodT0();  // nodes at current time point
  DATABRPTR* C1 = na->pNodT1();  // nodes at previous time point
  double dc;

  for (int ii=1; ii<mtp->nC; ii++)    // node iteration
   {
     bool NeedCopy = false;
     for(int ic=0; ic < CH->nICb-1; ic++) // do not check charge
     {
        dc = C0[ii]->bIC[ic] - C1[ii]->bIC[ic];
        if( fabs( dc ) > min( mtp->cdv, (C1[ii]->bIC[ic] * 1e-3)))
                  NeedCopy = true;
     }
     if( NeedCopy )
       na->CopyNodeFromTo( ii, mtp->nC, C1, C0 );
   }  // ii    end of node iteration loop
}

#ifndef IPMGEMPLUGIN

//-------------------------------------------------------------------
// NewNodeArray()  make worked DATACH structure
// Allocates a new FMT node array of DATABR structures and
// reads in the data from MULTI (using nC, Bn and DiCp, DDc, HydP)
//-------------------------------------------------------------------
void  TGEM2MT::NewNodeArray()
{
 // generate Tval&Pval arrays
 if( mtp->PsTPai != S_OFF )
    gen_TPval();

 na->MakeNodeStructures( mtp->nICb, mtp->nDCb,  mtp->nPHb,
      mtp->xIC, mtp->xDC, mtp->xPH,
      mtp->Tval, mtp->Pval,
      mtp->nTai,  mtp->nPai, mtp->Tai[3], mtp->Pai[3]  );
 DATACH* data_CH = na->pCSD();

 // put DDc
 if( data_CH->DD && mtp->DDc )
  for( int jj=0; jj<data_CH->nDCs; jj ++)
      data_CH->DD[jj*data_CH->nPp*data_CH->nTp] = mtp->DDc[jj];

 for( mtp->kv = 0; mtp->kv < mtp->nIV; mtp->kv++ )
 {
   pVisor->Message( window(), GetName(),
      "Generation of EqStat records\n"
           "Please, wait...", mtp->kv, mtp->nIV);

  // Make new Systat record & calculate it
     gen_task();
  // Save databr
     na->packDataBr();
  //
   for( int jj=0; jj<mtp->nC; jj ++)
    if( mtp->DiCp[jj][0] == mtp->kv )
     {    na->setNodeHandle( jj );
          na->MoveWorkNodeToArray( jj, mtp->nC,  na->pNodT0());
          na->CopyWorkNodeFromArray( jj, mtp->nC,  na->pNodT0() );
     }
    mt_next();      // Generate work values for the next EqStat rkey

 } // mtp->kv

 pVisor->CloseMessage();

 DATABRPTR* C0 = na->pNodT0();  // nodes at current time point
 for( int ii=0; ii<mtp->nC; ii++)
      if(  C0[ii] == 0 )
      Error( "NewNodeArray() error:" ," Undefined boundary condition!" );

 // put HydP
  for( int jj=0; jj<mtp->nC; jj ++)
   {
        C0[jj]->NodeTypeHY = mtp->DiCp[jj][1];
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
   }

  for (int ii=0; ii<mtp->nC; ii++)    // node iteration
   {
       na->CopyNodeFromTo( ii, mtp->nC, C0, na->pNodT1() );
   }  // ii    end of node iteration loop

}

#endif


// This function (for 1D-case only) analyses the node array
// to decide where to use AIA and where PIA for GEM calculations
// in each node.
// if IAmode == NEED_GEM_AIA then all nodes are set to AIA
// If IAmode == NEED_GEM_PIA then the adjacent nodes with different
// reactive phase assemblage (very different pe, pH), as well as their
// neighbors from each side, are set to AIA. Other nodes are set to PIA.
// returns: Number of nodes set to AIA
//
int TGEM2MT::CheckPIAinNodes1D( char IAmode, int start_node, int end_node )
{
       int nSetAIA = 0;
       // Getting direct access to data
       DATACH* CH = na->pCSD();       // DataCH structure
       DATABRPTR* C0 = na->pNodT0();  // nodes at previous time point
//       DATABRPTR* C1 = na->pNodT1();  // nodes at current time point
       bool* iaN = na->piaNode();      // indicators for IA in the nodes

	   start_node = max( start_node, 0 );
       end_node = min( end_node, (int)mtp->nC-1 );

       // Initializing iaNode vector
       if( IAmode == NEED_GEM_PIA && CH->nDCb == CH->nDC )
       {
    	  for( int ii = start_node; ii<= end_node; ii++ )
    		 iaN[ii] = false;    // potentially need PIA
       }
       else { // setting all nodes to AIA GEM calculations
     	  for( int ii = 0; ii< (int)mtp->nC; ii++ )
     	  {
     		  iaN[ii] = true;
     		  nSetAIA++;
     	  }
          return nSetAIA;
       }

       // This is done only if PIA mode is requested!
       // First variant: 1-neighbour algorithm
       for( int ii = start_node; ii<= end_node; ii++) // node iteration
	   {
    	 // mtp->qc = (short)ii;
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
	     if( (ii == 0 || ii == (int)mtp->nC-1 ) ) // &&
//	    	 (C0[ii]->NodeTypeMT == normal || C0[ii]->NodeTypeMT == NBC1source) )
//	     {
//	     	iaN[ii] = true;
//	     	nSetAIA++;
	     	continue;
//	     }
	     // checking pair of adjacent nodes for phase assemblage differences
	     for( int kk=0; kk<CH->nPHb; kk++ )
	     {
	    	if( C0[ii]->xPH[kk] == 0.0 && C0[ii-1]->xPH[kk] == 0.0 )
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
bool TGEM2MT::CalcIPM( char mode, int start_node, int end_node, FILE* diffile )
{
   int Mode, ic, RetCode=OK_GEM_AIA;
   bool NeedGEM = false;
   bool iRet = true;
   double dc; // difference (decrement) to concentration/amount
   //  Getting direct access to TNodeArray class data

   DATACH* CH = na->pCSD();       // DataCH structure
   DATABRPTR* C0 = na->pNodT0();  // nodes at previous time point
   DATABRPTR* C1 = na->pNodT1();  // nodes at current time point
   bool* iaN = na->piaNode();     // indicators for IA in the nodes

   if(mtp->PvSIA == S_OFF )
   {  for(int ii=0; ii<(int)mtp->nC; ii++ )  // AIA mode forced 
	    iaN[ii] = true;
   }
   else 
	 for(int ii=0; ii<(int)mtp->nC; ii++ )  // PIA mode allowed (if mode==NEED_GEM_PIA) 
	   iaN[ii] = false;

   start_node = max( start_node, 0 );
   end_node = min( end_node, (int)mtp->nC-1 );

   for( int ii = start_node; ii<= end_node; ii++) // node iteration
   {
     mtp->qc = (short)ii;
     if(mtp->PvSIA == S_OFF )
        NeedGEM = true; 
     else {
    	 C1[ii]->IterDone = 0;
    	 NeedGEM = false; 
     }
    	 // Here we compare this node for current time and for previous time - works for AIA and PIA
	 for( ic=0; ic < CH->nICb; ic++)    // do we check charge here?
	 {     // It has to be checked on minimal allowed c0 value
		 if( C1[ii]->bIC[ic] < mtp->cez )
	     { // to prevent loss of Independent Component
			 C1[ii]->bIC[ic] = mtp->cez;
		 }
		 dc = C0[ii]->bIC[ic] - C1[ii]->bIC[ic];
		 if( fabs( dc ) > min( mtp->cdv, (C1[ii]->bIC[ic] * 1e-3 )))
			 NeedGEM = true;  // we still need to recalculate equilibrium
                          // in this node because its vector b has changed
	 }
     C1[ii]->bIC[CH->nICb-1] = 0.;   // zeroing charge off in bulk composition
//     NeedGEM = true;
     if( mode == NEED_GEM_PIA )
     {   // smart algorithm
    	 if( iaN[ii] == true )
    	 {
    		 Mode = NEED_GEM_AIA;
       	 }
    	 else {
    		 Mode = NEED_GEM_PIA;
    		 if( mtp->PvSIA == S_ON )   // force loading of primal solution into GEMIPM
    			 Mode *= -1;            // othervise use internal (old) primal solution
    	 }
     }
     else Mode = NEED_GEM_AIA;
// NeedGEM = true;  // Mode = NEED_GEM_AIA;   // debugging - calculating all nodes with AIA!
     if( NeedGEM )
     {
        RetCode = na->RunGEM( ii, Mode );
        // Returns GEMIPM2 calculation time in sec after the last call to GEM_run()
        mtp->TimeGEM +=	na->GEM_CalcTime();
        // checking RetCode from GEM IPM calculation
        if( !(RetCode==OK_GEM_AIA || RetCode == OK_GEM_PIA ))
        {
//          cout << "CalcIPM4 " << RetCode << endl;
          char buf[200];
          gstring err_msg;
          iRet = false;

          sprintf( buf, " Node= %-8d  Step= %-8d\n", ii, mtp->ct );
          err_msg = buf;
          switch( RetCode )
          {
                case BAD_GEM_AIA:
                      err_msg += "Bad GEM result using simplex IA";
                      break;
                case  ERR_GEM_AIA:
                      err_msg += "GEM calculation error using simplex IA";
                      break;
                case  BAD_GEM_PIA:
                      err_msg += "Bad GEM result using previous solution IA";
                      break;
                case  ERR_GEM_PIA:
                      err_msg += "GEM calculation error using previous solution IA";
                      break;
               case  TERROR_GEM:  err_msg +=  "Terminal error in GEMIPM2 module";
          }
          if( mtp->PvMO != S_OFF && diffile )
          {  fprintf( diffile, "\nError reported from GEMIPM2 module\n%s\n",
                    err_msg.c_str() );
          }
#ifndef IPMGEMPLUGIN
          else
          {  err_msg += "\n Continue?";
             if( !vfQuestion( window(),
                 "Error reported from GEMIPM2 module",err_msg.c_str() ))
                     Error("Error reported from GEMIPM2 module",
                     "Process stopped by the user");
          }
#endif
        }
     }
     else { // GEM calculation for this node not needed
//       C0[ii]->IterDone = 0;
       C1[ii]->IterDone = 0; // number of GEMIPM iterations is set to 0
     }
   }  // ii   end of node iteration loop
   return iRet;
}

// The mass transport iteration time step
void TGEM2MT::MassTransAdvecStart()
{
    mtp->dx = mtp->cLen/mtp->nC;
    mtp->dTau = 0.5*(mtp->dx/mtp->fVel)*1/mtp->tf;
    mtp->oTau = 0;
    mtp->cTau = mtp->Tau[START_];
    // mtp->cTau = 0;
    mtp->ct = 0;
}

// The mass transport iteration time step
void TGEM2MT::MassTransParticleStart()
{
    double dt_adv, dt_dif;
    mtp->dx = mtp->cLen/mtp->nC;
    dt_adv = (mtp->dx / mtp->fVel);    // Courant - Neumann criteria
    dt_dif = (mtp->dx*mtp->dx)/2./(mtp->al_in*mtp->fVel + mtp->Dif_in);
    mtp->dTau = min(dt_adv, dt_dif) / mtp->tf;
//    mtp->dTau = 0.5*(mtp->dx/mtp->fVel)*1./mtp->tf;  // Courant criterion
    mtp->oTau = 0;
    mtp->cTau = mtp->Tau[START_];
    // mtp->cTau = 0;
    mtp->ct = 0;
    pa->setUpCounters();
}


// The mass transport iteration time step
void TGEM2MT::MassTransParticleStep( bool CompMode )
{
   mtp->ct += 1;
   mtp->oTau = mtp->cTau;
   mtp->cTau += mtp->dTau;
   pa->GEMPARTRACK( mtp->PsMode, CompMode, mtp->oTau, mtp->cTau );
}

// The mass transport iteration time step
// CompMode: true: Transport via dependent components in Aq phase (except H2O)
//           false: Transport via dissolved indepedent components (with H2O)
//
void TGEM2MT::MassTransAdvecStep( bool CompMode )
{
 double c0, c1, cm1,  cm2, cmax, c12, cm12, 
        charge, c0new, dc, cr, aji, fmolal;      // some help variables
 int ii, ic, jc;

 //  Getting direct access to TNodeArray class data
 DATACH* CH = na->pCSD();       // DataCH structure
 // DATABRPTR* C0 = na->pNodT0();  // nodes at current time point
 // DATABRPTR* C1 = na->pNodT1();  // nodes at previous time point

   mtp->ct += 1;
   mtp->cTau += mtp->dTau;
   cr = mtp->fVel * mtp->dTau/mtp->dx;

   for( ii = 2; ii< mtp->nC-1; ii++) // node iteration, -1 the right boundary is open ....
   {
     mtp->qc = (short)ii;
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
            		 aji = nodeCH_A( jc, ic );
            		 if( aji )
            		     node1_bIC(ii, ic) -= aji * dc;
            	 }
    	     }
    	     else {  // setting the DC amount to minimum
    	    	 node1_xDC( ii, jc ) = mtp->cez;
    	    	 dc = c0/fmolal - mtp->cez;
    	    	 for( ic=0; ic<CH->nICb; ic++)  // incrementing independent components
            	 {
            		 aji = nodeCH_A( jc, ic );
            		 if( aji )
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
    		 niw = node1_xDC( ii, CH->nDCinPH[0]-1 )* nodeCH_A( CH->nDCinPH[0]-1, ic );
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
    		 if( (c0/fmolal - dc) > mtp->cez ) // min( mtp->cez, node1_bIC(ii, ic) * 1e-4 ))
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
	 charge = node1_bIC(ii, CH->nICb-1 );
node1_bIC(ii, CH->nICb-1 ) = 0.0;		// debugging
   } // end of loop over nodes
}

// Main call for the mass transport iterations in 1D case
//
// mode is NEED_GEM_AIA or NEED_GEM_PIA (see DATABR.H) mode of GEM initial approximation
// If NEED_GEM_PIA then the program will try smart initial approximation for the nodes
// (PIA when possible, AIA in the vicinity of fronts, pH and pe barriers).
// returns true, if cancelled/interrupted by the user; false if finished Ok
bool TGEM2MT::Trans1D( char mode )
{
  int evrt =10;
  bool iRet = false;
  bool CompMode;   // Component transport mode: true: DC; false: IC
  int nStart = 0, nEnd = mtp->nC;
  int NodesSetToAIA;

FILE* logfile;
FILE* ph_file;
FILE* diffile = NULL;

if( mtp->PsDDc == S_ON && mtp->PsDIc == S_OFF )  // Set of DC transport using record switches
	CompMode = true; 
if( mtp->PsDDc == S_OFF && mtp->PsDIc == S_ON )  // Set of IC transport using record switches
	CompMode = false; 

if( mtp->PvMO != S_OFF )
{
// Preparations: opening output files for monitoring 1D profiles
logfile = fopen( "ICaq-log.dat", "w+" );    // Total dissolved element molarities
if( !logfile)
  return iRet;
ph_file = fopen( "Ph-log.dat", "w+" );   // Mole amounts of phases
if( !ph_file)
  return iRet;
diffile = fopen( "ICdif-log.dat", "w+" );   //  Element amount diffs for t and t-1
if( !diffile)
  return iRet;
}

// time scales testing
clock_t t_start, t_end, t_out, t_out2;
clock_t outp_time = (clock_t)0;
t_start = clock();
mtp->TimeGEM = 0.0; 

   if( mtp->iStat != AS_RUN  )
   {  switch( mtp->PsMode )
     {
       case GMT_MODE_A:   // A: 1D advection (numerical) coupled FMT scoping model
    	                MassTransAdvecStart();
                        nStart = 1; nEnd = mtp->nC-1;
                 break;
       case GMT_MODE_V:   // V: const volume coupled FMT scoping model
       case GMT_MODE_W:   // W: random-walk advection-diffusion coupled FMT scoping model
    	                MassTransParticleStart();
                 break;
      default: // more mass transport models here
               break;
     }
//
   // Calculation of chemical equilibria in all nodes at the beginning
   // with the Simplex initial approximation
     CalcIPM( NEED_GEM_AIA, nStart, nEnd, diffile );
   }
   mtp->iStat = AS_READY;


if( mtp->PvMO != S_OFF )
{
// Data collection for monitoring: Initial state (at t=0)
t_out = clock();
na->logDiffsIC( diffile, mtp->ct, mtp->cTau/(365*86400), mtp->nC, 1 );
na->logProfileAqIC( logfile, mtp->ct, mtp->cTau/(365*86400), mtp->nC, 1 );
na->logProfilePhMol( ph_file, mtp->ct, mtp->cTau/(365*86400), mtp->nC, 1 );
t_out2 = clock();
outp_time += ( t_out2 - t_out);
}

#ifndef IPMGEMPLUGIN
      if(  mtp->PvMSg != S_OFF && vfQuestion(window(),
             GetName(), "Use graphic monitoring?") )
        {
            RecordPlot( 0 );
        }
#endif

//  This loop contains the mass transport iteration time step
     do {   // time iteration step

#ifndef IPMGEMPLUGIN
       iRet = pVisor->Message( window(), GetName(),
           "Calculating Reactive Mass Transport (RMT)\n"
           "Please, wait (may take long)...", mtp->ct, mtp->ntM );
#endif
       if( iRet )
         break;

        //  the mass transport iteration time step
         switch( mtp->PsMode )
         {
          case GMT_MODE_A: MassTransAdvecStep( CompMode );
                  break;
          case GMT_MODE_V:
          case GMT_MODE_W: MassTransParticleStep( CompMode );
                                      break;
          default: // more mass transport models here
                  break;
        }

        // The analysis of GEM IA modes in nodes - optional
//        NodesSetToAIA = CheckPIAinNodes1D( mode, nStart, nEnd );
         
        //   Here we call a loop on GEM calculations over nodes
        //   parallelization should affect this loop only
        CalcIPM( mode, nStart, nEnd, diffile );

if( mtp->PvMO != S_OFF )
{
t_out = clock();
na->logDiffsIC( diffile, mtp->ct, mtp->cTau/(365*86400), mtp->nC, evrt );
    // logging differences after the MT iteration loop
t_out2 = clock();
outp_time += ( t_out2 -  t_out);
}
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
          if( mtp->PsMode == GMT_MODE_V ||
                mtp->PsMode == GMT_MODE_W )
            pa->CopyfromT1toT0();

if( mtp->PvMO != S_OFF )
{
// Data collection for monitoring: Current state
t_out = clock();
na->logProfileAqIC( logfile, mtp->ct, mtp->cTau/(365*86400), mtp->nC, evrt );
na->logProfilePhMol( ph_file, mtp->ct, mtp->cTau/(365*86400), mtp->nC, evrt );
t_out2 = clock();
outp_time += ( t_out2 - t_out);
}

     } while ( mtp->cTau < mtp->Tau[STOP_] && mtp->ct < mtp->ntM );


t_end = clock();
double dtime = ( t_end- t_start );
double clc_sec = CLOCKS_PER_SEC;

if( mtp->PvMO != S_OFF )
{
fprintf( diffile,
  "\nTotal time of calculation %lg s;  Time of output %lg s;  Whole run time %lg s;  Pure GEM run time %lg s\n",
    (dtime-outp_time)/clc_sec,  outp_time/clc_sec, dtime/clc_sec, mtp->TimeGEM );
fclose( logfile );
fclose( ph_file );
fclose( diffile );
}

#ifndef IPMGEMPLUGIN
pVisor->CloseMessage();
#endif

  return iRet;
}



// --------------------- end of m_gem2mtt.cpp ---------------------------

