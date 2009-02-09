//--------------------------------------------------------------------
// $Id: main.cpp 717 2006-07-06 08:06:21Z gems $
//
// Demo test of usage of the TNode class for implementing a simple
// direct coupling scheme between FMT and GEM in a single-GEM-call
// fashion, assuming that the chemical speciation and all dynamic
// parameter data are kept in the FMT part, which calls GEMIPM
// calculation once per node.

// TNode class implements a  simple C/C++ interface between GEM IPM
// and FMT codes. Works with DATACH and work DATABR structures
// without using the TNodearray class
//
// Copyright (C) 2006 S.Dmytriyeva, D.Kulik
//
// This file is part of GEMIPM2K code for thermodynamic modelling
// by Gibbs energy minimization

// This file may be distributed under the licence terms defined
// in GEMIPM2K.QAL
//
// See also http://les.web.psi.ch/Software/GEMS-PSI
// mailto://gems2.support@psi.ch
//-------------------------------------------------------------------

#include <time.h>
#include <math.h>
#include <string.h>

#include "node.h"

int outTest = 0; 

#define nNodes  2 // set here how many nodes you need

int main( int argc, char* argv[] )
 {
   long int nTimes = 1e3;   // Maximum number of time iteration steps

	// Analyzing command line arguments
     // Default arguments
     char ipm_input_file_list_name[256] = "chemsys.lst";
     char dbr_input_file_name[256] = "chemsys-dbr.dat";
     char fmt_input_file_name[256] = "fmtparam.dat";

     if (argc >= 2 )
       strncpy( ipm_input_file_list_name, argv[1], 256);
         // list of files needed as input for initializing GEMIPM2K
     if (argc >= 3 )
           strncpy( dbr_input_file_name, argv[2], 256);
             // input file for boundary conditions
     if (argc >= 4 )
       strncpy( fmt_input_file_name, argv[3], 256);
         // your optional file with FMT input parameters

   // Creating TNode structure accessible trough node pointer
   TNode* node  = new TNode();

   // Here we read the files needed as input for initializing GEMIPM2K
   // The easiest way to prepare them is to use GEMS-PSI code (GEM2MT module)
   if( node->GEM_init( ipm_input_file_list_name ) )
   {
	   cout << "error occured during reading the files" ;
	   return 1;  // error occured during reading the files
   }

   // allocations and defaults for other FMT parameters can be added here

   // Here you can read your file with some FMT parameters and initial data
   // if( my_fmt_input(fmt_input_file_name) )
   //   return 2;

   // Number of ICs, DCs, Phases and Phases-solutions kept in the node
   // DATABR structure for exchange with GEMIPM - for your convenience
   long int nIC, nDC, nPH, nPS;
   //long int i,   j,   k,   ks;    // indices for direct access to components
                             // and phases data in the DataCH framework

   // Getting direct access to DataCH structure in GEMIPM2K memory
   DATACH* dCH = node->pCSD();
   if( !dCH  )
       return 3;

   // Getting direct access to work node DATABR structure which
   // exchanges data between GEMIPM and FMT parts
   DATABR* dBR = node->pCNode();
   if( !dBR  )
       return 4;

   // Extracting data bridge array sizes
   nIC = dCH->nICb;
   nDC = dCH->nDCb;
   nPH = dCH->nPHb;
   nPS = dCH->nPSb;

   // Allocating work memory for FMT part (here only chemical variables)
   // for one node only (real FMT problems consider many nodes)
   // Names are consistent with the DataBridge structure (see "databr.h")

   long int m_NodeHandle[nNodes], m_NodeStatusCH[nNodes], m_IterDone[nNodes];

   double m_T[nNodes], m_P[nNodes], m_Vs[nNodes], m_Ms[nNodes],
          m_Gs[nNodes], m_Hs[nNodes], m_IC[nNodes], m_pH[nNodes], m_pe[nNodes],
          m_Eh[nNodes];

   double *m_xDC, *m_gam, *m_xPH, *m_aPH, *m_vPS, *m_mPS,*m_bPS,
         *m_xPA, *m_dul, *m_dll, *m_bIC, *m_rMB, *m_uIC;

   m_bIC = (double*)malloc( nNodes*nIC*sizeof(double) );
   m_rMB = (double*)malloc( nNodes*nIC*sizeof(double) );
   m_uIC = (double*)malloc( nNodes*nIC*sizeof(double) );
   m_xDC = (double*)malloc( nNodes*nDC*sizeof(double) );
   m_gam = (double*)malloc( nNodes*nDC*sizeof(double) );
   m_dul = (double*)malloc( nNodes*nDC*sizeof(double) );
   m_dll = (double*)malloc( nNodes*nDC*sizeof(double) );
   m_aPH = (double*)malloc( nNodes*nPH*sizeof(double) );
   m_xPH = (double*)malloc( nNodes*nPH*sizeof(double) );
   m_vPS = (double*)malloc( nNodes*nPS*sizeof(double) );
   m_mPS = (double*)malloc( nNodes*nPS*sizeof(double) );
   m_bPS = (double*)malloc( nNodes*nIC*nPS*sizeof(double) );
   m_xPA = (double*)malloc( nNodes*nPS*sizeof(double) );

   // (1) ---------------------------------------------
   // Initialization of GEMIPM and chemical data kept in the FMT part
   // Can be done in a loop over nodes if there are many nodes

   cout << "Begin Initialiation part" << endl;

   long int in;
   for(  in=1; in<nNodes; in++ )
   {
     dBR->NodeStatusCH = NEED_GEM_AIA; // direct access to node DATABR structure

     // re-calculating equilibrium by calling GEMIPM
     m_NodeStatusCH[in] = node->GEM_run(1., false);


     if( !( m_NodeStatusCH[in] == OK_GEM_AIA || m_NodeStatusCH[in] == OK_GEM_SIA ) )
     {
    	 cout << "Bad task in Initialiation part";
    	 return 5;
     }
     // Extracting chemical data into FMT part
     node->GEM_restore_MT( m_NodeHandle[in], m_NodeStatusCH[in], m_T[in],
       m_P[in], m_Vs[in], m_Ms[in],
       m_bIC+in*nIC, m_dul+in*nDC, m_dll+in*nDC, m_aPH+in*nPH );
        // Extracting GEMIPM output data to FMT part
     node->GEM_to_MT( m_NodeHandle[in], m_NodeStatusCH[in], m_IterDone[in],
       m_Vs[in], m_Ms[in], m_Gs[in], m_Hs[in], m_IC[in], m_pH[in], m_pe[in],
       m_Eh[in], m_rMB+in*nIC, m_uIC+in*nIC, m_xDC+in*nDC, m_gam+in*nDC,
       m_xPH+in*nPH, m_vPS+in*nPS, m_mPS+in*nPS,
       m_bPS+in*nIC*nPS, m_xPA+in*nPS );

     // Here the file output for the initial conditions can be implemented
   }

  // Initialization of GEMIPM and chemical data kept in the FMT part
  // Can be done in a loop over boundary nodes

  // Read DATABR structure from text file (read boundary condition)
      TNode::na->GEM_read_dbr( dbr_input_file_name );

  for(  in=0; in<1; in++ )
  {
   dBR->NodeStatusCH = NEED_GEM_AIA; // direct access to node DATABR structure

  // re-calculating equilibrium by calling GEMIPM
   m_NodeStatusCH[in] = node->GEM_run( 1., false );

   if( !( m_NodeStatusCH[in] == OK_GEM_AIA || m_NodeStatusCH[in] == OK_GEM_SIA ) )
  {
 	 cout << "Bad task in Initialiation part";
     return 5;
  }
  // Extracting chemical data into FMT part
   node->GEM_restore_MT( m_NodeHandle[in], m_NodeStatusCH[in], m_T[in],
    m_P[in], m_Vs[in], m_Ms[in],
    m_bIC+in*nIC, m_dul+in*nDC, m_dll+in*nDC, m_aPH+in*nPH );
     // Extracting GEMIPM output data to FMT part
   node->GEM_to_MT( m_NodeHandle[in], m_NodeStatusCH[in], m_IterDone[in],
    m_Vs[in], m_Ms[in], m_Gs[in], m_Hs[in], m_IC[in], m_pH[in], m_pe[in],
    m_Eh[in], m_rMB+in*nIC, m_uIC+in*nIC, m_xDC+in*nDC, m_gam+in*nDC,
    m_xPH+in*nPH, m_vPS+in*nPS, m_mPS+in*nPS,
    m_bPS+in*nIC*nPS, m_xPA+in*nPS );

  // Here the file output for the initial conditions can be implemented
 }

   cout << "End Initialiation part" << endl;

   clock_t t_start11, t_end11;
   t_start11 = clock();

   // (2) ----------------------------------------------
   // Work loop for the coupled FMT-GEM modelling

   cout << "Begin Coupled Modelling part" << endl;
   long int xCa = node->IC_name_to_xDB("Ca");
   long int xSi = node->IC_name_to_xDB("Si");
   long int xO = node->IC_name_to_xDB("O");
   long int xH = node->IC_name_to_xDB("H");
   long int xCalcite = node->Ph_name_to_xDB("Calcite");
   long int xSiO2 = node->Ph_name_to_xDB("Silica-amorph");

   double dC, dS, iC[nNodes], iS[nNodes], mass[nNodes];
   double nC = 100, nS = 200; // Number of different point by Calcite and SiO2
   dC = (1. - 0./*1.e-7*/)/nC; 
   dS = (1. - 0./*1.e-7*/)/nS; 
   for( in=0; in<nNodes; in++ )
   {
	   iC[in] =1.;  
	   iS[in] =-1.;
	   mass[in] = 0.;
   }
   
   cout << "Start Tnode test: " << ipm_input_file_list_name << " "
          << dbr_input_file_name << endl;
    cout << " nNodes = " << nNodes << "  nTimes = " << nTimes
          << "  deltaC = " << dC << " deltaS = " << dS << endl;

   // Checking indexes
   cout << "xCa= " << xCa << " xSi=" << xSi << " xO=" << xO << " xH=" << xH
        << " xCalcite=" << xCalcite << " xSiO2=" << xSiO2 << endl;


   for( long int it=0; it<nTimes; it++ )  // iterations over time
   {

     // Loop over nodes for calculating the mass transport step
     for(  in=0; in<nNodes; in++ )
     {
       if( it > 0 )
       {
           if( m_bIC[in*nIC+xCa] <= (1e-7+dC) )
        	   iC[in] = 1.;
           else
               if( m_bIC[in*nIC+xCa] > 1 )
            	   iC[in] = -1.;
        	   
    	   m_bIC[in*nIC+xCa]+= dC*iC[in];
           m_bIC[in*nIC+xH] += 2.*dC*iC[in];
           m_bIC[in*nIC+xO] += 2.*dC*iC[in];
           
           mass[in] += dC*iC[in] * dCH->ICmm[/*IC_xDB_to_xCH*/(xCa)];
           mass[in] += 2. * dC*iC[in] * dCH->ICmm[/*IC_xDB_to_xCH*/(xH)];
           mass[in] += 2. * dC*iC[in] * dCH->ICmm[/*IC_xDB_to_xCH*/(xO)];
           
           if( m_bIC[in*nIC+xSi] <= (1e-7+dS) )
        	   iS[in] = 1.;
           else
               if( m_bIC[in*nIC+xSi] > 1 )
            	   iS[in] = -1.;

           m_bIC[in*nIC+xSi] += dS*iS[in];
           m_bIC[in*nIC+xO]  += 2.*dS*iS[in];
           mass[in] += dS*iS[in] * dCH->ICmm[/*IC_xDB_to_xCH*/(xSi)];
           mass[in] += 2. * dS*iS[in] * dCH->ICmm[/*IC_xDB_to_xCH*/(xO)];
       }
       else
       {
   		for (long int i=0;i<dCH->nICb;i++){
    			mass[in] += m_bIC[in*nIC+i]*dCH->ICmm[i];
       }
      }
     } 

     //cout << " it = " << it << "  dt = " << dt << "  tc = " << tc << endl;
     
     //     cout << " Chemical loop begins: " << endl;
     // Loop over nodes for calculating the chemical equilibration step
     for( in=0; in<nNodes; in++ )
     {
        m_NodeHandle[in] = in;
        m_NodeStatusCH[in] = NEED_GEM_SIA; // or NEED_GEM_SIA

        // Setting input data for GEMIPM
        node->GEM_from_MT( m_NodeHandle[in], m_NodeStatusCH[in],
             m_T[in], m_P[in], m_Vs[in], m_Ms[in],
             m_bIC+in*nIC, m_dul+in*nDC, m_dll+in*nDC, m_aPH+in*nPH );
        node->GEM_set_MT( (double)it, 1. );

if( in == 1 && it == 363 )
{	   
    //     sprintf(NextRecipeOutFileName , "%s.out", NextRecipeFileName );
    //     cout << "  See output in the " << NextRecipeOutFileName << " file" << endl;
        node->GEM_write_dbr( "dbr_363_before.out", false, false );
    //      sprintf(NextRecipeOutFileName , "%s.Dump.out", NextRecipeFileName );
    //      cout << "See dump output in the " << NextRecipeOutFileName << " file" << endl;
       node->GEM_print_ipm( "ipm_363_before.out"  );
       outTest = 1;
     
}     
        // Calling GEMIPM calculation
        m_NodeStatusCH[in] = node->GEM_run( mass[in]*1e-3, true );

if( in == 1 && it == 363 )
{	   
            //     sprintf(NextRecipeOutFileName , "%s.out", NextRecipeFileName );
            //     cout << "  See output in the " << NextRecipeOutFileName << " file" << endl;
                node->GEM_write_dbr( "dbr_363_after.out", false, false );
            //      sprintf(NextRecipeOutFileName , "%s.Dump.out", NextRecipeFileName );
            //      cout << "See dump output in the " << NextRecipeOutFileName << " file" << endl;
               node->GEM_print_ipm( "ipm_363_after.out"  );
               outTest = 0;
}     
        if( ( m_NodeStatusCH[in] == ERR_GEM_AIA ||
               m_NodeStatusCH[in] == ERR_GEM_SIA ) )
        {
        	cout << "Error Calling GEMIPM calculation it = " << it << " Node " << in 
            << " Ca= " << m_bIC[in*nIC+xCa] << " Si= " << m_bIC[in*nIC+xSi] << endl;
        }
        else
        {	
            if( ( m_NodeStatusCH[in] == BAD_GEM_AIA ||
                   m_NodeStatusCH[in] == BAD_GEM_SIA ) )
            {
            	cout << "Warning Calling GEMIPM calculation it = " << it << " Node " << in 
                << " Ca= " << m_bIC[in*nIC+xCa] << " Si= " << m_bIC[in*nIC+xSi] << endl;
            }	
        // Extracting GEMIPM output data to FMT part
          node->GEM_to_MT( m_NodeHandle[in], m_NodeStatusCH[in], m_IterDone[in],
          m_Vs[in], m_Ms[in], m_Gs[in], m_Hs[in], m_IC[in], m_pH[in], m_pe[in],
          m_Eh[in],m_rMB+in*nIC, m_uIC+in*nIC, m_xDC+in*nDC, m_gam+in*nDC,
          m_xPH+in*nPH, m_vPS+in*nPS, m_mPS+in*nPS,
          m_bPS+in*nIC*nPS, m_xPA+in*nPS  );

        // Here the debug print for each node in can be implemented
        if( !(it % (nTimes/100)))
          cout << "  in = " << in << "  it = " << it 
               << " Cal= " << m_xPH[in*nPH+xCalcite] 
               << " Silica= " << m_xPH[in*nPH+xSiO2] << endl;
        }  
     }
    
     // Here the output for the current state at tc can be implemented
  }

   t_end11 = clock();
  double dtime = ( t_end11- t_start11 );
  double clc_sec = CLOCKS_PER_SEC;
  cout <<  "Total time of calculation  s; " <<  (dtime)/clc_sec << endl;
  // cout << " End Coupled Modelling part" << endl;

  // (3) ----------------------------------------------
  // Calculations finished - t_end reached

  // freeing dynamic arrays
  free( m_xDC );
  free( m_gam );
  free( m_xPH );
  free( m_vPS );
  free( m_mPS );
  free( m_bPS );
  free( m_xPA );
  free( m_dul );
  free( m_dll );
  free( m_bIC );
  free( m_rMB );
  free( m_uIC );

  // deleting GEMIPM and data exchange memory structures
  delete node;

  cout << endl << "Finished Ok" << endl;

  return 0;
}

//---------------------------------------------------------------------------
// end of main.cpp for node class usage example
