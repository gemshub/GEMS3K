//--------------------------------------------------------------------
// $Id$
//
// Demo test of usage of the TNode class for implementing a simple
// direct coupling scheme between FMT and GEM in a single-GEM-call
// fashion, assuming that the chemical speciation and all dynamic
// parameter data are kept in the FMT part, which calls GEM IPM
// calculation once per node.

// TNode class implements a  simple C/C++ interface between GEMS3K
// and FMT codes. Works with DATACH and work DATABR structures
//
// Copyright (C) 2006,2012 S.Dmytriyeva, D.Kulik, G.Kosakowski
//
// This file is part of the GEMS3K code for thermodynamic modelling
// by Gibbs energy minimization
//
// See also http://gems.web.psi.ch/GEMS3K/
// mailto://gems2.support@psi.ch
//-------------------------------------------------------------------

#include "main.h"

//The case of data exchange in computer memory
int main( int argc, char* argv[] )
{
   // Analyzing command line arguments ( Default arguments)
   char ipm_input_file_list_name[256] = "system-dat.lst";
   char dbr_input_file_name[256] = "system-dbr.dat";

   if (argc >= 2 )  // list of files needed as input for initializing GEMIPM2K
       strncpy( ipm_input_file_list_name, argv[1], 256);
   if (argc >= 3 ) // input file for boundary conditions
       strncpy( dbr_input_file_name, argv[2], 256);

    // Creates TNode structure instance accessible trough the "node" pointer
    TNode* node  = new TNode();

    // (1) Initialization of GEMS3K internal data by reading  files
    //     whose names are given in the ipm_input_system_file_list_name
    if( node->GEM_init( ipm_input_file_list_name ) )
    {
          cout << "Error occured during reading the files" ;
          return 1;
    }

    // Getting direct access to work node DATABR structure which exchanges the
    // data with GEMS3K (already filled out by reading the DBR input file)
    DATABR* dBR = node->pCNode();

    // Getting direct access to DataCH structure in GEMS3K instance memory
    DATACH* dCH = node->pCSD();

    // Creating memory for mass transport nodes
    // 11 nodes, 99 time steps
    TMyTransport mt( 11, 100, dCH->nICb, dCH->nDCb, dCH->nPHb, dCH->nPSb );

    // Initialization of GEMS3K and chemical information for nodes kept in the MT part
    long int in;
    for(  in=1; in< mt.nNodes; in++ )
    {
        // Asking GEM IPM to run with automatic initial approximation
        dBR->NodeStatusCH = NEED_GEM_AIA;
        // (2) re-calculating equilibrium by calling GEMIPM2K, getting the status back
        mt.aNodeStatusCH[in] = node->GEM_run( false);
        if( !( mt.aNodeStatusCH[in] == OK_GEM_AIA || mt.aNodeStatusCH[in] == OK_GEM_SIA ) )
        {
              cout << "Error occured during re-calculating equilibrium" ;
              return 5;
        }

        // Extracting GEM IPM input data to mass-transport program arrays
        node->GEM_restore_MT( mt.aNodeHandle[in], mt.aNodeStatusCH[in], mt.aT[in], mt.aP[in],
            mt.aVs[in], mt.aMs[in], mt.abIC[in], mt.adul[in], mt.adll[in], mt.aaPH[in] );
          
        // Extracting GEM IPM output data to mass-transport program arrays
        node->GEM_to_MT( mt.aNodeHandle[in], mt.aNodeStatusCH[in], mt.aIterDone[in],
            mt.aVs[in], mt.aMs[in], mt.aGs[in], mt.aHs[in], mt.aIC[in], mt.apH[in], mt.ape[in],
            mt.aEh[in], mt.arMB[in], mt.auIC[in], mt.axDC[in], mt.agam[in], mt.axPH[in],
            mt.avPS[in], mt.amPS[in], mt.abPS[in], mt.axPA[in], mt.aaPh[in], mt.abSP[in] );

        // Here the setup of initial differences between node compositions,
        //    temperatures, etc. can be implemented
        //
        // Here the file output for the initial conditions can be implemented
    }

    // Read DATABR structure from text file (read boundary condition on the left)
    node->GEM_read_dbr( dbr_input_file_name );

    for(  in=0; in<1; in++ )
    {
        // Asking GEM IPM to run with automatic initial approximation
        dBR->NodeStatusCH = NEED_GEM_AIA;
        // (2) Re-calculating chemical equilibrium by calling GEM
        mt.aNodeStatusCH[in] = node->GEM_run( false );
        if( !( mt.aNodeStatusCH[in] == OK_GEM_AIA || mt.aNodeStatusCH[in] == OK_GEM_SIA ) )
        {
              cout << "Error occured during re-calculating chemical equilibrium" ;
              return 5;
        }

        // (6) Extracting GEMIPM input data to mass-transport program arrays
        node->GEM_restore_MT( mt.aNodeHandle[in], mt.aNodeStatusCH[in], mt.aT[in], mt.aP[in],
            mt.aVs[in], mt.aMs[in], mt.abIC[in], mt.adul[in], mt.adll[in], mt.aaPH[in] );
          
        // (7) Extracting GEMIPM output data to mass-transport program arrays
        node->GEM_to_MT( mt.aNodeHandle[in], mt.aNodeStatusCH[in], mt.aIterDone[in],
            mt.aVs[in], mt.aMs[in], mt.aGs[in], mt.aHs[in], mt.aIC[in], mt.apH[in], mt.ape[in],
            mt.aEh[in], mt.arMB[in], mt.auIC[in], mt.axDC[in], mt.agam[in], mt.axPH[in],
            mt.avPS[in], mt.amPS[in], mt.abPS[in], mt.axPA[in], mt.aaPh[in], mt.abSP[in] );

        // Here the setup of initial differences between node compositions,
        //    temperatures, etc. can be implemented
        //
        // Here the file output for the initial conditions can be implemented
    }
      
    // Main loop - iterations over nTimes time steps
    int xCalcite = node->Ph_name_to_xDB("Calcite");
    int xDolomite = node->Ph_name_to_xDB("Dolomite-dis");
    int xAq_gen = node->Ph_name_to_xDB("aq_gen");
    long int ICndx[5];
    ICndx[0] = node->IC_name_to_xDB("Ca");
    ICndx[1] = node->IC_name_to_xDB("C");
    ICndx[2] = node->IC_name_to_xDB("O");
    ICndx[3] = node->IC_name_to_xDB("Mg");
    ICndx[4] = node->IC_name_to_xDB("Cl");
    // Checking indexes
    cout << "xCa= " << ICndx[0] << " xC=" << ICndx[1] << " xO=" << ICndx[2] << " xMg="
         << ICndx[3] << " xCl=" << ICndx[4] << endl << " xCalcite=" << xCalcite
         << " xDolomite=" << xDolomite << " xAq_gen=" << xAq_gen << endl << endl;
    double stoich[5] = { 0., 0., 0., 1., 2. }; // defines what is 'transported'

    long int it;
    for( it=0; it< mt.nTimes; it++ )
    {
       cout << "Time step  " << it << endl;
       // Mass transport loop over nodes (not a real transport model)
       mt.OneTimeStepRun( stoich, ICndx, 5 );

       // Chemical equilibration loop over nodes
       for( in=0; in< mt.nNodes; in++ )
       {
          mt.aNodeHandle[in] = in;
          mt.aNodeStatusCH[in] = NEED_GEM_SIA;
          // (8) Setting input data for GEM IPM to use available node speciation as
          // initial approximation
          node->GEM_from_MT( mt.aNodeHandle[in], mt.aNodeStatusCH[in],
                  mt.aT[in], mt.aP[in], mt.aVs[in], mt.aMs[in],
                  mt.abIC[in], mt.adul[in], mt.adll[in], mt.aaPH[in], mt.axDC[in], mt.agam[in] );
          // (9)   Passing current FMT iteration information into the work DATABR structure
          node->GEM_set_MT( (double)it, 1. );
 
          // Calling GEMIPM calculation
          mt.aNodeStatusCH[in] = node->GEM_run( true );

          if( ( mt.aNodeStatusCH[in] == ERR_GEM_AIA || mt.aNodeStatusCH[in] == ERR_GEM_SIA ||
                        mt.aNodeStatusCH[in] ==  T_ERROR_GEM ) )
          {
              cout << "Error: GEM calculation results are not retrieved. Time step"
                 << it << " node " << in << endl;
          }
          else
          {
            if( ( mt.aNodeStatusCH[in] == BAD_GEM_AIA || mt.aNodeStatusCH[in] == BAD_GEM_SIA  ) )
            {
               cout << "Insufficient quality of GEM solution, but GEM results are retrieved"
               << it << " node " << in << endl;
            }              
            else // (7) Extracting GEMIPM output data to FMT part
              node->GEM_to_MT( mt.aNodeHandle[in], mt.aNodeStatusCH[in], mt.aIterDone[in],
                mt.aVs[in], mt.aMs[in], mt.aGs[in], mt.aHs[in], mt.aIC[in], mt.apH[in], mt.ape[in],
                mt.aEh[in],mt.arMB[in], mt.auIC[in], mt.axDC[in], mt.agam[in], mt.axPH[in],
                mt.avPS[in], mt.amPS[in], mt.abPS[in], mt.axPA[in], mt.aaPh[in], mt.abSP[in] );
          }
          // Here, the output upon completion of the time step is usually implemented
          //  to monitor the coupled simulation or collect results
          cout << "  Node " << in ;
          cout << ": Aq= " << mt.axPH[in][xAq_gen] << " pH= " << mt.apH[in] <<
                  "  Calcite= " << mt.axPH[in][xCalcite] << endl;
      }
  }
  // Calculations finished - end time reached

  // Final output e.g. of total simulation time or of the final distribution of
  //  components and phases in all nodes can be implemented here

  // deleting GEM IPM and data exchange memory structures
  delete node;
  // end of example
  return 0;
}


TMyTransport::TMyTransport( long int p_nNod, long int p_nTim, long int p_nIC, long int p_nDC,
              long int p_nPH, long int p_nPS )
{

    nNodes = p_nNod;
    nTimes = p_nTim;
    nIC = p_nIC;
    nDC = p_nDC;
    nPH = p_nPH;
    nPS = p_nPS;

    aNodeHandle = new long int [nNodes];
    aNodeStatusCH = new long int [nNodes];
    aIterDone = new long int [nNodes];

    aT = new double [nNodes];
    aP = new double [nNodes];
    aVs = new double [nNodes];
    aMs = new double [nNodes];
    aGs = new double [nNodes];
    aHs = new double [nNodes];
    aIC = new double [nNodes];
    apH = new double [nNodes];
    ape = new double [nNodes];
    aEh = new double [nNodes];

    axDC = new double *[nNodes];
    agam = new double *[nNodes];
    axPH = new double *[nNodes];
    aaPH = new double *[nNodes];
    avPS = new double *[nNodes];
    amPS = new double *[nNodes];
    abPS = new double *[nNodes];
    axPA = new double *[nNodes];
    aaPh = new double *[nNodes];
    adul = new double *[nNodes];
    adll = new double *[nNodes];
    abIC = new double *[nNodes];
    arMB = new double *[nNodes];
    auIC = new double *[nNodes];
    abSP = new double *[nNodes];

    for (long int in=0; in<nNodes; in++)
    {
         abIC[in] = new double [nIC];
         arMB[in] = new double [nIC];
         auIC[in] = new double [nIC];
         axDC[in] = new double [nDC];
         agam[in] = new double [nDC];
         adul[in] = new double [nDC];
         adll[in] = new double [nDC];
         aaPH[in] = new double [nPH];
         axPH[in] = new double [nPH];
         avPS[in] = new double [nPS];
         amPS[in] = new double [nPS];
         axPA[in] = new double [nPS];
         aaPh[in] = new double [nPH];
         abPS[in] = new double [nIC*nPS];
         abSP[in] = new double [nIC];
    }
}

TMyTransport::~TMyTransport()
{
    // Deleting chemical data arrays for nodes
    for (long int in=0; in<nNodes; in++)
    {
        delete[]abIC[in];
        delete[]arMB[in];
        delete[]auIC[in];
        delete[]axDC[in];
        delete[]agam[in];

        delete[]adul[in];
        delete[]adll[in];

        delete[]aaPH[in];
        delete[]axPH[in];
        delete[]avPS[in];
        delete[]amPS[in];
        delete[]axPA[in];
        delete[]aaPh[in];
        delete[]abPS[in];
        delete[]abSP[in];
    }
// return;
    delete[]axDC;
    delete[]agam;
    delete[]axPH;
    delete[]aaPH;
    delete[]avPS;
    delete[]amPS;
    delete[]abPS;
    delete[]axPA;
    delete[]aaPh;
    delete[]adul;
    delete[]adll;
    delete[]abIC;
    delete[]arMB;
    delete[]auIC;
    delete[]abSP;

    delete[]aNodeHandle;
    delete[]aNodeStatusCH;
    delete[]aIterDone;
    delete[]aT;
    delete[]aP;
    delete[]aVs;
    delete[]aMs;
    delete[]aGs;
    delete[]aHs;
    delete[]aIC;
    delete[]apH;
    delete[]ape;
    delete[]aEh;
}

// A very simple example of transport algorithm
void TMyTransport::OneTimeStepRun( double *stoicf, long int *ICndx, long int nICndx )
{
    double parcel[nICndx];
    long int in;
    for(  in=1; in< nNodes; in++ )
    {
    // some operators that change in some nodes some amounts of some migrating
    // chemical species (axDC array or  some amounts of migrating chemical
    // elements in abIC array), possibly using data from other arrays in such
    // a way that the mass conservation within the whole array of nodes is retained
        for (int i=0; i<nICndx; i++)
        {
            parcel[i] = 0.01* stoicf[i]*abIC[in-1][ICndx[i]];
            abIC[in][ICndx[i]] += parcel[i];
//            if( in > 1 )
//                abIC[in-1][ICndx[i]] -= parcel[i];
        }
    // The above example loop implements a zero-order flux of MgCl2 in one direction.
    // Real advective/diffusive transport models are much more complex, but essentially
    //   do similar things
    }
}

//---------------------------------------------------------------------------
// end of main.cpp for node class usage example
