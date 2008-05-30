//--------------------------------------------------------------------
// $Id: main.cpp 182 2008-05-27 08:06:21Z gems $
//
// gemnode
// Demo test of usage of the TNode class for implementing a simple
// direct coupling scheme between FMT and GEM in a single-GEM-call
// fashion, assuming that the chemical speciation and all dynamic
// parameter data are kept in the FMT part, which calls GEMIPM
// calculation once per node.

// TNode class implements a  simple C/C++ interface between GEM IPM
// and FMT codes. Works with DATACH and work DATABR structures
// without using the TNodearray class
//
// Copyright (C) 2006,2008 S.Dmytrieva, D.Kulik
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
const int  nNodes =  5;   // set here how many nodes you need

int main( int argc, char* argv[] )
 {
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

// Creating TNode structure accessible trough the pointer node
   TNode* node  = new TNode();

// Here we read the files needed as input for initializing GEMIPM2K
  // The easiest way to prepare them is to use GEMS-PSI code
   if( node->GEM_init( ipm_input_file_list_name ) )
       return 1;  // error occured when reading the files (unix2dos?)

// int nNodes = 5;     // number of local equilibrium nodes, 1 or more
   int nTimes = 100;      // Maximum number of time iteration steps
   unsigned int nTotIt;   // Number of GEM iterations per time step
   double t_start = 0., t_end = 10000., dt = 100., tc = 1.;
   double CalcTime = 0.0;
   double MeanIt, MeanFIA, MeanIPM, MeanPRL;  // mean numbers of GEM iterations
   int nPrecLoops =0;     // Number of IPM-2 precision enhancement loops
   int nIterFIA =0;       // Number of EnterFeasibleDomain() iterations
   int nIterIPM =0;       // Number of main IPM iterations
   int nIterTotal=0;      // Total FIA+IPM iterations done

   cout << "Start of gemnode PIA test: " << ipm_input_file_list_name << " "
         << dbr_input_file_name << endl;
   cout << " nNodes = " << nNodes << "  nTimes = " << nTimes
         << "  t_start = " << t_start << " t_end = " << t_end
         << "  dt = " << dt << "  tc = " << tc << endl;

// allocations and defaults for other FMT parameters can be added here
// . . . . . . . . . . . . . . . .
// Here you can read your file with some FMT parameters and initial data
   // if( my_fmt_input(fmt_input_file_name) )
   //   return 2;

// Number of ICs, DCs, Phases and Phases-solutions kept in the node DATABR
// structure for exchange with GEMIPM - for your convenience
   int nIC, nDC, nPH, nPS;
//   int i,   j,   k,   ks;    // indices for direct access to components
                               // and phases data in the DataCH framework

// Getting direct access to DATACH structure in GEMIPM2K memory
   DATACH* dCH = node->pCSD();
   if( !dCH  )
       return 3;

// Getting direct access to the work node DATABR structure which exchanges
// the data between GEM and FMT parts
   DATABR* dBR = node->pCNode();
   if( !dBR  )
       return 4;

// Extracting data bridge array sizes
   nIC = dCH->nICb;
   nDC = dCH->nDCb;
   nPH = dCH->nPHb;
   nPS = dCH->nPSb;

// Allocating work memory for FMT part (here only chemical variables)
   // for nNodes (a few) nodes (real FMT problems consider many more nodes).
   // Names of arrays are taken consistent with the DATABR structure
   //  (see "databr.h") for better readability

   short m_NodeHandle[nNodes], m_NodeStatusCH[nNodes], m_IterDone[nNodes];
   int nPrecL[nNodes], nFIA[nNodes], nIPM[nNodes];

   double m_T[nNodes], m_P[nNodes], m_Vs[nNodes], m_Ms[nNodes],
          m_Gs[nNodes], m_Hs[nNodes], m_IC[nNodes], m_pH[nNodes],
	  m_pe[nNodes], m_Eh[nNodes];

   double *m_xDC, *m_gam, *m_xPH, *m_aPH, *m_vPS, *m_mPS,*m_bPS,
         *m_xPA, *m_dul, *m_dll, *m_bIC, *m_rMB, *m_uIC;

   m_bIC = new double[ nNodes*nIC ];
   m_rMB = new double[ nNodes*nIC ];
   m_uIC = new double[ nNodes*nIC ];
   m_xDC = new double[ nNodes*nDC ];
   m_gam = new double[ nNodes*nDC ];
   m_dul = new double[ nNodes*nDC ];
   m_dll = new double[ nNodes*nDC ];
   m_aPH = new double[ nNodes*nPH ];
   m_xPH = new double[ nNodes*nPH ];
   m_vPS = new double[ nNodes*nPS ];
   m_mPS = new double[ nNodes*nPS ];
   m_xPA = new double[ nNodes*nPS ];
   m_bPS = new double[ nNodes*nIC*nPS ];

// (1) ---------------------------------------------
// Initialization of chemical data kept in the FMT part.
   // Can be done in a loop over nodes, if there are many nodes
   //   cout << "Begin Initialisation Part" << endl;
   int in;
   for(  in=1; in<nNodes; in++ )  // Node with in=0 will be initialized later
   {
     dBR->NodeStatusCH = NEED_GEM_AIA; // direct access to node DATABR structure

// re-calculating equilibrium by calling GEMIPM2K, getting the status
     m_NodeStatusCH[in] = node->GEM_run( false );

     if( !( m_NodeStatusCH[in] == OK_GEM_AIA || m_NodeStatusCH[in] == OK_GEM_SIA ) )
        return 5; // GEM IPM did not converge properly

// Extracting GEM IPM input chemical data into FMT part
     node->GEM_restore_MT( m_NodeHandle[in], m_NodeStatusCH[in],
       m_T[in], m_P[in], m_Vs[in], m_Ms[in],
       m_bIC+in*nIC, m_dul+in*nDC, m_dll+in*nDC, m_aPH+in*nPH );

// Extracting GEM IPM output chemical data into FMT part
     node->GEM_to_MT( m_NodeHandle[in], m_NodeStatusCH[in], m_IterDone[in],
       m_Vs[in], m_Ms[in], m_Gs[in], m_Hs[in], m_IC[in], m_pH[in], m_pe[in],
       m_Eh[in], m_rMB+in*nIC, m_uIC+in*nIC, m_xDC+in*nDC, m_gam+in*nDC,
       m_xPH+in*nPH, m_vPS+in*nPS, m_mPS+in*nPS,
       m_bPS+in*nIC*nPS, m_xPA+in*nPS );

//  Uncomment this to test variable pressures and temperatures
//         m_T[in] += in*5;
//         m_P[in] += (in-1)*20;
//         m_T[in] += in*7;
//         m_P[in] += (in-1)*20;

// Here the file output for the initial conditions can be implemented
// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
   }

// Read DATABR structure from text file separately (boundary condition)
      TNode::na->GEM_read_dbr( dbr_input_file_name );

  for(  in=0; in<1; in++ ) // Initialising boundary condition node(s)
  {
   dBR->NodeStatusCH = NEED_GEM_AIA; // activating GEM IPM for automatic initial
                                     // approximation
// re-calculating equilibrium by calling GEMIPM using previous primal solution in this node
//   m_NodeStatusCH[in] = node->GEM_run( true );
   // re-calculating equilibrium by calling GEMIPM using previous content of GEMIPM structure
      m_NodeStatusCH[in] = node->GEM_run( false );
  if( !( m_NodeStatusCH[in] == OK_GEM_AIA || m_NodeStatusCH[in] == OK_GEM_SIA ) )
     return 5;
// Extracting GEM IPM input chemical data into FMT part
   node->GEM_restore_MT( m_NodeHandle[in], m_NodeStatusCH[in], m_T[in],
    m_P[in], m_Vs[in], m_Ms[in],
    m_bIC+in*nIC, m_dul+in*nDC, m_dll+in*nDC, m_aPH+in*nPH );
// Extracting GEM IPM output chemical data to FMT part
   node->GEM_to_MT( m_NodeHandle[in], m_NodeStatusCH[in], m_IterDone[in],
    m_Vs[in], m_Ms[in], m_Gs[in], m_Hs[in], m_IC[in], m_pH[in], m_pe[in],
    m_Eh[in], m_rMB+in*nIC, m_uIC+in*nIC, m_xDC+in*nDC, m_gam+in*nDC,
    m_xPH+in*nPH, m_vPS+in*nPS, m_mPS+in*nPS,
    m_bPS+in*nIC*nPS, m_xPA+in*nPS );

// Here the file output for the initial boundary conditions can be implemented
// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
   }

   cout << "End of gem_node test initialisation" << endl;
   clock_t t_start11, t_end11;
   t_start11 = clock();  // Setting time when the coupled modelling starts

// (2) ----------------------------------------------
// Work loop for the coupled FMT-GEM modelling

   cout << "Begin Coupled Modelling Part" << endl;
      // Getting DATABR indexes for chemical species to be monitored

   int xiC = node->IC_name_to_xDB("C");
   int xiCa = node->IC_name_to_xDB("Ca");
   int xiH = node->IC_name_to_xDB("H");
   int xiO = node->IC_name_to_xDB("O");
   int xiSr = node->IC_name_to_xDB("Sr");
   int xiZz = node->IC_name_to_xDB("Zz");
   int xiCl = node->IC_name_to_xDB("Cl");
   int xCal = node->Ph_name_to_xDB("(Sr,Ca)CO3(reg)");
   int xStr = node->Ph_name_to_xDB("Strontianite");
   // Checking indexes
   cout << " xiC= " << xiC << "  xiCa= " << xiCa << "  xiH= " << xiH
        << "  xiO= " << xiO << "  xiSr=" << xiSr << "  xiZz=" << xiZz
        << "  xCalcite=" << xCal << "  xStrontianite=" << xStr << endl;
   nTotIt = 0;

   for( int it=0; it<nTimes; it++ )  // iterations over time
   {
     int in;
//   cout << " FMT loop begins: " << endl;
// Loop over nodes for calculating the mass transport step
   // ( this whole loop may be replaced by call(s) to FMT subroutines)
     for(  in=1; in<nNodes; in++ )
     {
       ; // Add here some operators as function of tc and dt
       // that transfer mass between nodes (i.e. change some elements of
       // m_bIC array). Take care about mass conservation and consistency
       // of chemical compositions of the nodes.
       // in this example, simply adding SrCl2 to m_bIC vector
       // in order to cause the conversion of calcite to strontianite
       if( it > 0 )
       {
         if( xiSr >= 0 )
           m_bIC[in*nIC+xiSr] += dt*in*1e-6;
         if( xiCl >= 0 )
           m_bIC[in*nIC+xiCl] += dt*in*2e-6;
       }
       // Alternatively, the transport may be simulated using m_xDC arrays.
       // In this case, an alternative (overloaded) call to GEM_from_MT()
       // should be used (see below)
     }
//     cout << " FMT loop ends: ";
     cout << " it = " << it << "  tc = " << tc << endl;

//     cout << " Chemical loop begins: " << endl;

// Loop over nodes for calculating the chemical equilibration step
     for( in=0; in<nNodes; in++ )
     {
//if( !in)
        cout << "  in = " << in << "  T = " << m_T[in];
        m_NodeHandle[in] = in;

// Below you can switch between AIA and PIA initial approximation modes
//         m_NodeStatusCH[in] = NEED_GEM_AIA;    // tests are marked *.out2A 
        m_NodeStatusCH[in] = NEED_GEM_SIA;      // tests are marked *.out2P

// Setting input data for GEM IPM
        node->GEM_from_MT( m_NodeHandle[in], m_NodeStatusCH[in],
            m_T[in], m_P[in], m_Vs[in], m_Ms[in],
            m_bIC+in*nIC, m_dul+in*nDC, m_dll+in*nDC, m_aPH+in*nPH // );
	    ,m_xDC+in*nDC, m_gam+in*nDC );   // this overload call works also
                               // in NEED_GEM_AIA mode, but only in v.2.2!
//if(!in)
//    for( int j=0; j<nDC; j++ )
//      cout << "  " << m_xDC[in*nDC+j];
// if(!in)
//    for( int j=0; j<nDC; j++ )
//      cout << "  " << m_gam[in*nDC+j];
//  Alternative call to correct bulk chemical composition using changed
//     m_xDC (chemical species amounts) data. Take care that this function
//     actually compresses the xDC values into mole amounts of elements
//     and adds them to the elements of bIC vector.
//   node->GEM_from_MT( m_NodeHandle[in], m_NodeStatusCH[in],
//     m_T[in], m_P[in], m_Vs[in], m_Ms[in],
//     m_bIC+in*nIC, m_dul+in*nDC, m_dll+in*nDC, m_aPH+in*nPH, m_xDC+in*nDC );

// Calling GEM IPM2 calculation
// re-calculating equilibrium by calling GEMIPM using previous primal solution in this node
    m_NodeStatusCH[in] = node->GEM_run( true );
// re-calculating equilibrium by calling GEMIPM using previous content of GEMIPM structure
//       m_NodeStatusCH[in] = node->GEM_run( false );
       if( !( m_NodeStatusCH[in] == OK_GEM_AIA || m_NodeStatusCH[in] == OK_GEM_SIA ) )
            return 5;
        CalcTime += node->GEM_CalcTime();  // Incrementing calculation time - only v.2.2.0
        nIterTotal += node->GEM_Iterations( nPrecL[in], nFIA[in], nIPM[in] );
        nPrecLoops += nPrecL[in];
        nIterFIA += nFIA[in];
        nIterIPM += nIPM[in];

// Extracting GEM IPM output data to FMT part
        node->GEM_to_MT( m_NodeHandle[in], m_NodeStatusCH[in], m_IterDone[in],
          m_Vs[in], m_Ms[in], m_Gs[in], m_Hs[in], m_IC[in], m_pH[in], m_pe[in],
          m_Eh[in],m_rMB+in*nIC, m_uIC+in*nIC, m_xDC+in*nDC, m_gam+in*nDC,
          m_xPH+in*nPH, m_vPS+in*nPS, m_mPS+in*nPS,
          m_bPS+in*nIC*nPS, m_xPA+in*nPS  );

        nTotIt += m_IterDone[in];

// Here the debug print for each node can be implemented
// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

// if( !in )
//{ //        cout << " Gem run ends: ";
/*
        cout << "  bC= " << m_bIC[in*nIC+xiC] << "  bCa= " << m_bIC[in*nIC+xiCa]
             << "  bH= " << m_bIC[in*nIC+xiH] << "  bO= " << m_bIC[in*nIC+xiO]
             << "  bSr= " << m_bIC[in*nIC+xiSr] << "  bZz= " << m_bIC[in*nIC+xiZz];
//      cout << endl;  */
        cout << "     Cal= " << m_xPH[in*nPH+xCal] <<
                "  Str= " << m_xPH[in*nPH+xStr];
        cout << "  [Ca]= " << m_bPS[in*nIC*nPS+xiCa] <<
                "  [Sr]= " << m_bPS[in*nIC*nPS+xiSr];
        cout << "  pH= " << m_pH[in] << "  It= " <<  m_IterDone[in]; // << endl;
        cout << "  nPRL= " << nPrecL[in] << "  nFIA= "  << nFIA[in] << "  nIPM= " << nIPM[in] << endl;
//}
   }
//    cout << " Chemical loop ends: " << endl;

// Here the output for the current transport state at tc can be implemented
// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    tc += dt;  // time increment
  }

  t_end11 = clock(); // getting end time of coupled calculations
  double dtime = ( t_end11- t_start11 );
  double clc_sec = CLOCKS_PER_SEC;
  MeanIt = double(nIterTotal)/500.; // per 1 GEMIPM2 call ( nTotIt double(nNodes*nTimes);
  MeanFIA = double(nIterFIA)/500.;
  MeanIPM = double(nIterIPM)/500.;
  MeanPRL = double(nPrecLoops)/500.;

  cout <<  "Pure GEM IPM2 runtime , s: " <<  CalcTime << endl; // Only v. 2.2.0
  cout <<  "Total time of calculation, s: " <<  (dtime)/clc_sec << endl;
  cout << "    Mean GEMIPM2 iterations per node: " << MeanPRL << " " <<
               MeanFIA << " " << MeanIPM << " " << MeanIt << endl;

  cout << " This gem_node test ";

// (3) ----------------------------------------------
// Calculations finished - t_end reached

// freeing dynamic arrays
  delete[] m_xDC;
  delete[] m_gam;
  delete[] m_xPH;
  delete[] m_vPS;
  delete[] m_mPS;
  delete[] m_bPS;
  delete[] m_xPA;
  delete[] m_dul;
  delete[] m_dll;
  delete[] m_bIC;
  delete[] m_rMB;
  delete[] m_uIC;
  delete[] m_aPH;

// deleting GEMIPM and data exchange memory structures
  delete node;

  cout << "has finished Ok" << endl;

  return 0;
}

//---------------------------------------------------------------------------
// end of main.cpp for gemnode - the TNode class usage example
