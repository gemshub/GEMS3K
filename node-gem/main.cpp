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
// See also http://gems.web.psi.ch/
// mailto://gems2.support@psi.ch
//-------------------------------------------------------------------

#include <time.h>
#include <math.h>
#include <string.h>

#include "node.h"

//The case of data exchange in computer memory
int main( int argc, char* argv[] )
 {
   // Analyzing command line arguments ( Default arguments)
     char ipm_input_file_list_name[256] = "system-dat.lst";
//     char dbr_input_file_name[256] = "system-dbr.dat";

     if (argc >= 2 )  // list of files needed as input for initializing GEMIPM2K
       strncpy( ipm_input_file_list_name, argv[1], 256);
//     if (argc >= 3 ) // input file for boundary conditions
//       strncpy( dbr_input_file_name, argv[2], 256);

     // Creates TNode structure instance accessible trough the "node" pointer
      TNode* node  = new TNode();

      // (1) Initialization of GEMIPM2K internal data by reading  files
      //     whose names are given in the input_system_file_list_name
      if( node->GEM_init( ipm_input_file_list_name ) )
      {
          cout << "Error occured during reading the files" ;
          return 1;
      }

      // Getting direct access to work node DATABR structure which exchanges the
      // data with GEM IPM2 (already filled out by reading the DBR input file)
      DATABR* dBR = node->pCNode();

      // Getting direct access to DataCH structure in GEMIPM2K instance memory
      DATACH* dCH = node->pCSD();
      // Extracting data bridge array sizes
      long int nIC = dCH->nICb;
      long int nDC = dCH->nDCb;
      long int nPH = dCH->nPHb;
      long int nPS = dCH->nPSb;

      long int nNodes = 11; // Number of nodes in the transport part
                              // (can be much greater in a real transport code)    
     
      // Allocating memory for keeping node data (real transport codes may do it differently)   
      long int *m_NodeHandle, *m_NodeStatusCH, *m_IterDone;
      m_NodeHandle = new long int [nNodes];
      m_NodeStatusCH = new long int [nNodes];
      m_IterDone = new long int [nNodes];

      double *m_T, *m_P, *m_Vs, *m_Ms, *m_Gs, *m_Hs, *m_IC, *m_pH, *m_pe, *m_Eh;
      m_T = new double [nNodes];
      m_P = new double [nNodes];
      m_Vs = new double [nNodes];
      m_Ms = new double [nNodes];
      m_Gs = new double [nNodes];
      m_Hs = new double [nNodes];
      m_IC = new double [nNodes];
      m_pH = new double [nNodes];
      m_pe = new double [nNodes];
      m_Eh = new double [nNodes];

      double **m_xDC, **m_gam, **m_xPH, **m_aPH, **m_vPS, **m_mPS,**m_bPS,
             **m_xPA, **m_dul, **m_dll, **m_bIC, **m_rMB, **m_uIC;
      m_xDC = new double *[nNodes];
      m_gam = new double *[nNodes];
      m_xPH = new double *[nNodes];
      m_aPH = new double *[nNodes];
      m_vPS = new double *[nNodes];
      m_mPS = new double *[nNodes];
      m_bPS = new double *[nNodes];
      m_xPA = new double *[nNodes];
      m_dul = new double *[nNodes];
      m_dll = new double *[nNodes];
      m_bIC = new double *[nNodes];
      m_rMB = new double *[nNodes];
      m_uIC = new double *[nNodes];

      for (long int in=0; in<nNodes; in++)
      {
         m_bIC[in] = new double [nIC];
         m_rMB[in] = new double [nIC];
         m_uIC[in] = new double [nIC];
         m_xDC[in] = new double [nDC];
         m_gam[in] = new double [nDC];
         m_dul[in] = new double [nDC];
         m_dll[in] = new double [nDC];
         m_aPH[in] = new double [nPH];
         m_xPH[in] = new double [nPH];
         m_vPS[in] = new double [nPS];
         m_mPS[in] = new double [nPS];
         m_xPA[in] = new double [nPS];
         m_bPS[in] = new double [nIC*nPS];
      }

cout << "Begin Initialiation part" << endl;
      // Initialization of GEMIPM2K and chemical information for nodes kept in the FMT part
      long int in;
      for(  in=0; in<nNodes; in++ )
      {
    	  // Asking GEM to run with automatic initial approximation
    	  dBR->NodeStatusCH = NEED_GEM_AIA;
    	  // (2) re-calculating equilibrium by calling GEMIPM2K, getting the status back
          m_NodeStatusCH[in] = node->GEM_run(1., false);
          if( !( m_NodeStatusCH[in] == OK_GEM_AIA || m_NodeStatusCH[in] == OK_GEM_SIA ) )
          {
              cout << "Error occured during re-calculating equilibrium" ;
              return 5;
          }
          // (6) Extracting GEMIPM input data to mass-transport program arrays
          node->GEM_restore_MT( m_NodeHandle[in], m_NodeStatusCH[in], m_T[in], m_P[in],
            m_Vs[in], m_Ms[in], m_bIC[in], m_dul[in], m_dll[in], m_aPH[in] );
          
          // (7) Extracting GEMIPM output data to mass-transport program arrays
          node->GEM_to_MT( m_NodeHandle[in], m_NodeStatusCH[in], m_IterDone[in],
            m_Vs[in], m_Ms[in], m_Gs[in], m_Hs[in], m_IC[in], m_pH[in], m_pe[in],
            m_Eh[in], m_rMB[in], m_uIC[in], m_xDC[in], m_gam[in], m_xPH[in],
            m_vPS[in], m_mPS[in], m_bPS[in], m_xPA[in] );

          // Here the setup of initial differences between node compositions,
          //    temperatures, etc. can be implemented
          // Here the file output for the initial conditions can be implemented
        }

cout << "End Initialiation part" << endl;
        int xCalcite = node->Ph_name_to_xDB("Calcite");
        int xDolomite = node->Ph_name_to_xDB("Dolomite-dis");
      
      // Main loop - iterations over nTimes time steps
      long int it, nTimes = 99;
      for( it=0; it<nTimes; it++ ) 
      {
         // Mass transport loop over nodes (here not a real transport model)
         double parcel[3];
         for (int i=0; i<3; i++)
         {
            parcel[i] = 0.01*m_bIC[0][i];
         }
         for(  in=1; in<nNodes; in++ )
         {
              // some operators that change in some nodes some amounts of some migrating chemical species (m_xDC array)
              // (or  some amounts of migrating chemical elements in m_bIC array), possibly using data from other arrays
              // in such a way that the mass conservation within the whole array of nodes  is retained
           for (int i=0; i<3; i++)
           {
              m_bIC[in][i] += parcel[i];
              m_bIC[in-1][i] -= parcel[i];
           }
           // The above example loop implements a zero-order flux of the first three
           //   chemical elements in one direction
           // real advective/diffusive transport models are much more complex, but essentially
           //   do similar things
        }

        // Chemical equilibration loop over nodes
        for( in=0; in<nNodes; in++ )
        {
          m_NodeHandle[in] = in;
          m_NodeStatusCH[in] = NEED_GEM_SIA;
          // (8) Setting input data for GEMIPM to use available node speciation as
          // initial approximation
          node->GEM_from_MT( m_NodeHandle[in], m_NodeStatusCH[in],
                  m_T[in], m_P[in], m_Vs[in], m_Ms[in],
                  m_bIC[in], m_dul[in], m_dll[in], m_aPH[in], m_xDC[in], m_gam[in] );
          // (9)   Passing current FMT iteration information into the work DATABR structure
          node->GEM_set_MT( (double)it, 1. );
 
          // Calling GEMIPM calculation
          m_NodeStatusCH[in] = node->GEM_run( 1., true );

       	  if( ( m_NodeStatusCH[in] == ERR_GEM_AIA || m_NodeStatusCH[in] == ERR_GEM_SIA ) )
          {
              cout << "Error: GEM calculation results are not retrieved. Time step"
                 << it << " node " << in << endl;
          }
           else
           {   
            if( ( m_NodeStatusCH[in] == BAD_GEM_AIA || m_NodeStatusCH[in] == BAD_GEM_SIA ) )
            {
               cout << "Warning about insufficient quality of GEM solution, but GEM results are retrieved"
               << it << " node " << in << endl;
            }              
            else // (7) Extracting GEMIPM output data to FMT part
              node->GEM_to_MT( m_NodeHandle[in], m_NodeStatusCH[in], m_IterDone[in],
                m_Vs[in], m_Ms[in], m_Gs[in], m_Hs[in], m_IC[in], m_pH[in], m_pe[in],
                m_Eh[in],m_rMB[in], m_uIC[in], m_xDC[in], m_gam[in], m_xPH[in],
                m_vPS[in], m_mPS[in], m_bPS[in], m_xPA[in] );
           }
     // Here, the output upon completion of the time step is usually implemented
     //  to monitor the coupled simulation or collect results
        cout << "Time step" << it << " node " << in ;
        cout << " Cal= " << m_xPH[in][xCalcite] <<
                " Dol= " << m_xPH[in][xDolomite] <<
                " pH= " << m_pH[in] << endl;
        }
   }

   // Calculations finished - end time reached

   // Final output e.g. of total simulation time or of the final distribution of
   //  components and phases in all nodes can be implemented here
cout << " End Coupled Modelling part" << endl;

   // Deleting chemical data arrays for nodes
   for (long int in=0; in<nNodes; in++)  
   {
      delete[]m_bIC[in];
      delete[]m_rMB[in];
      delete[]m_uIC[in];
      delete[]m_xDC[in];
      delete[]m_gam[in];
      delete[]m_dul[in];
      delete[]m_dll[in];
      delete[]m_aPH[in];
      delete[]m_xPH[in];
      delete[]m_vPS[in];
      delete[]m_mPS[in];
      delete[]m_xPA[in];
      delete[]m_bPS[in];
   }
   delete[]m_xDC;
   delete[]m_gam;
   delete[]m_xPH;
   delete[]m_aPH;
   delete[]m_vPS;
   delete[]m_mPS;
   delete[]m_bPS;
   delete[]m_xPA;
   delete[]m_dul;
   delete[]m_dll;
   delete[]m_bIC;
   delete[]m_rMB;
   delete[]m_uIC;

   delete[]m_NodeHandle;
   delete[]m_NodeStatusCH;
   delete[]m_IterDone;
   delete[]m_T;
   delete[]m_P;
   delete[]m_Vs;
   delete[]m_Ms;
   delete[]m_Gs;
   delete[] m_Hs;
   delete[]m_IC;
   delete[]m_pH;
   delete[]m_pe;
   delete[]m_Eh;

   // deleting GEMIPM and data exchange memory structures
   delete node;

  // end of example
  return 0;
}

//---------------------------------------------------------------------------
// end of main.cpp for node class usage example
