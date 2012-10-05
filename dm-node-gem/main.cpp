//--------------------------------------------------------------------
// Dan Miron
// Test program using the TNode class for calculating solubilities for
// different input recepies (dbr) files
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
//   if (argc >= 3 ) // input file for boundary conditions
//       strncpy( dbr_input_file_name, argv[2], 256);

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
    DATABR* dBR = node->pCNode(); // DM pointer to work with dBR data

    // Getting direct access to DataCH structure in GEMS3K instance memory
    DATACH* dCH = node->pCSD(); // DM pointer to work with dCH data

    // Creating memory for mass transport nodes
    // 5 experiments
    TMyExperiments mt( 5, dCH->nICb, dCH->nDCb, dCH->nPHb, dCH->nPSb );
    long int in=0;
    if (argc >= 3 )
     {
            char NextRecipeFileName[256];
            char NextRecipeOutFileName[300];
            char input_recipes_file_list_path[256-fileNameLength] = "/tp_test";
            char input_recipes_file_list_name[256];
            char (*recipes)[64];


            strncpy( input_recipes_file_list_name, argv[2], 256);
        // list of additional recipes (dbr.dat files) e.g. for simulation
        // of a titration or another irreversible process
            // Reading list of recipes names from file
            cout<< input_recipes_file_list_name << "aalaa" << endl;
            recipes = f_getfiles( input_recipes_file_list_name, input_recipes_file_list_path, in, ',' );
            cout << recipes[in]<<endl;

           for(int in=0; in < mt.nNodes; in++ )
        {
           // Trying to read the next file name
              sprintf(NextRecipeFileName , "%s%s", input_recipes_file_list_path, recipes[in] );

              cout<<recipes[in]<< " " <<endl;

          // (5) Reading the next DBR file with different input composition or temperature
          node->GEM_read_dbr( NextRecipeFileName );

          // Asking GEM IPM2 to run (faster) with smart initial approximation
          dBR->NodeStatusCH = NEED_GEM_SIA;
          mt.aNodeStatusCH[in] = node->GEM_run( false );
          if( mt.aNodeStatusCH[in] == OK_GEM_AIA || mt.aNodeStatusCH[in] == OK_GEM_SIA  )
          {
              sprintf(NextRecipeOutFileName , "%s.out", NextRecipeFileName );
              cout<<NextRecipeOutFileName<< " " <<endl;
              node->GEM_write_dbr( NextRecipeOutFileName, false, true );
              sprintf(NextRecipeOutFileName , "%s.Dump.out", NextRecipeFileName );
              node->GEM_print_ipm( NextRecipeOutFileName );
          }
          else {
                 // error message, debugging printout
                sprintf(NextRecipeOutFileName , "%s.Dump.out", NextRecipeFileName );
                node->GEM_print_ipm( NextRecipeOutFileName );
  //??              return 5; // GEM IPM did not converge properly
                }


          // Extracting GEM IPM input data to mass-transport program arrays
          node->GEM_restore_MT( mt.aNodeHandle[in], mt.aNodeStatusCH[in], mt.aT[in], mt.aP[in], ///// DM (6) Passes (copies) the GEMS3K input data from the work instance of DATABR structure.
              mt.aVs[in], mt.aMs[in], mt.abIC[in], mt.adul[in], mt.adll[in], mt.aaPH[in] );

          // Diferences between node temperatures and pressures
          switch (in)
          {
          case 0:
              mt.aT[in]=873.15;
              mt.aP[in]=1e+08;
              break;
          case 1:
              mt.aT[in]=773.15;
              mt.aP[in]=1e+08;
              break;
          case 2:
              mt.aT[in]=873.15;
              mt.aP[in]=1e+08;
              break;
          case 3:
              mt.aT[in]=973.15;
              mt.aP[in]=2e+08;
              break;
          case 4:
              mt.aT[in]=773.15;
              mt.aP[in]=1e+08;
              break;
           }

        }
     }
           // end of possible loop on input recipes


  // *********************** DM part of old node-gem modified ********************* //


//node->GEM_read_dbr( "test-dbr2.dat" );

//    // Initialization of GEMS3K and chemical information for nodes kept in the MT part
//    long int in=0;
//    for(  in=0; in<mt.nNodes; in++ ) // DM mt.nodes number of experiments
//    {

////        if (in==2) {
//////            node->GEM_init( ipm_input_file_list_name );
////        node->GEM_read_dbr( "test-dbr2.dat" );
////        }

//        // Asking GEM IPM to run with automatic initial approximation
//        dBR->NodeStatusCH = NEED_GEM_AIA;
//        // (2) re-calculating equilibrium by calling GEMIPM2K, getting the status back
//        mt.aNodeStatusCH[in] = node->GEM_run( false );
//        if( !( mt.aNodeStatusCH[in] == OK_GEM_AIA || mt.aNodeStatusCH[in] == OK_GEM_SIA ) )
//        {
//              cout << "Error occured during re-calculating equilibrium" ;
//              return 5;
//        }
//        // Extracting GEM IPM input data to mass-transport program arrays
//        node->GEM_restore_MT( mt.aNodeHandle[in], mt.aNodeStatusCH[in], mt.aT[in], mt.aP[in], ///// DM (6) Passes (copies) the GEMS3K input data from the work instance of DATABR structure.
//            mt.aVs[in], mt.aMs[in], mt.abIC[in], mt.adul[in], mt.adll[in], mt.aaPH[in] );

//        // Diferences between node temperatures and pressures
//        switch (in)
//        {
//        case 1:
//            mt.aT[in]=873.15;
//            mt.aP[in]=1e+08;
//            break;
//        case 2:
//            mt.aT[in]=773.15;
//            mt.aP[in]=1e+08;
//            break;
//        case 3:
//            mt.aT[in]=973.15;
//            mt.aP[in]=2e+08;
//            break;
//        case 4:
//            mt.aT[in]=773.15;
//            mt.aP[in]=75e+06;
//            break;
//         }

//                  std::strstream filename;
//                  filename << "test-dbr" << (int)in << ".dat";
//                  cout << filename.str();
////                  node->GEM_write_dbr(filename.str(), false, false, false);
////        node->GEM_write_dbr("test-dbr", false, false, false);

//        // Extracting GEM IPM output data to mass-transport program arrays // DM here we could extract the solubilities
//        node->GEM_to_MT( mt.aNodeHandle[in], mt.aNodeStatusCH[in], mt.aIterDone[in], //  DM Retrieves the GEMIPM2 chemical speciation calculation results from the work DATABR structure instance
//            mt.aVs[in], mt.aMs[in], mt.aGs[in], mt.aHs[in], mt.aIC[in], mt.apH[in], mt.ape[in],
//            mt.aEh[in], mt.arMB[in], mt.auIC[in], mt.axDC[in], mt.agam[in], mt.axPH[in],
//            mt.avPS[in], mt.amPS[in], mt.abPS[in], mt.axPA[in], mt.aaPh[in], mt.abSP[in] );


        // Here the setup of initial differences between node compositions,
        //    temperatures, etc. can be implemented
        //
        // Here the file output for the initial conditions can be implemented
//    }

//    Read DATABR structure from text file (read boundary condition on the left)
//    node->GEM_read_dbr( dbr_input_file_name );

//    for(  in=0; in<mt.nNodes; in++ )
//    {
//        // Asking GEM IPM to run with automatic initial approximation
//        dBR->NodeStatusCH = NEED_GEM_AIA;
//        // (2) Re-calculating chemical equilibrium by calling GEM
//        mt.aNodeStatusCH[in] = node->GEM_run( false );
//        if( !( mt.aNodeStatusCH[in] == OK_GEM_AIA || mt.aNodeStatusCH[in] == OK_GEM_SIA ) )
//        {
//              cout << "Error occured during re-calculating chemical equilibrium" ;
//              return 5;
//        }


////        // (6) Extracting GEMIPM input data to mass-transport program arrays
////        node->GEM_restore_MT( mt.aNodeHandle[in], mt.aNodeStatusCH[in], mt.aT[in], mt.aP[in],
////            mt.aVs[in], mt.aMs[in], mt.abIC[in], mt.adul[in], mt.adll[in], mt.aaPH[in] );

          
////        // (7) Extracting GEMIPM output data to mass-transport program arrays
////        node->GEM_to_MT( mt.aNodeHandle[in], mt.aNodeStatusCH[in], mt.aIterDone[in],
////            mt.aVs[in], mt.aMs[in], mt.aGs[in], mt.aHs[in], mt.aIC[in], mt.apH[in], mt.ape[in],
////            mt.aEh[in], mt.arMB[in], mt.auIC[in], mt.axDC[in], mt.agam[in], mt.axPH[in],
////            mt.avPS[in], mt.amPS[in], mt.abPS[in], mt.axPA[in], mt.aaPh[in], mt.abSP[in] );

//        // Here the setup of initial differences between node compositions,
//        //    temperatures, etc. can be implemented
//        //
//        // Here the file output for the initial conditions can be implemented
//    }


      // *********************** DM part of old node-gem modified ********************* //



      
    // Main loop - iterations over nTimes time steps
    int xCorundum = node->Ph_name_to_xDB("Corundum");
    int xAq_gen = node->Ph_name_to_xDB("aq_gen");
    long int ICndx[3]; // DM dCH->nICb - number of independent components kept in DBR file and DATABR memory structure * exchange with 7
    ICndx[0] = node->IC_name_to_xDB("Si");
    ICndx[1] = node->IC_name_to_xDB("Al");
    ICndx[2] = node->IC_name_to_xDB("Na");
    ICndx[3] = node->IC_name_to_xDB("Cl");
    // Checking indexes
    cout << "xSi= " << ICndx[0] << " xAl=" << ICndx[1] << " xNa=" << ICndx[2] << " xCl="
         << ICndx[3] << endl << " xCorundum=" << xCorundum
         << " xAq_gen=" << xAq_gen << endl << endl;

    // double stoich[5] = { 0., 0., 0., 1., 2. }; // defines what is 'transported'

//    long int it;
//    for( it=0; it< mt.nTimes; it++ )
//       {
//       cout << "Time step  " << it << endl;
//       // Mass transport loop over nodes (not a real transport model)
//       mt.OneTimeStepRun( mt );  // DM creates all nodes - somehow I have to create all experiments - read the class experiment,

       // Chemical equilibration loop over nodes
       for( in=0; in< mt.nNodes; in++ )
       {

//           if (in==2) {
////               node->GEM_init( ipm_input_file_list_name );
//           node->GEM_read_dbr( "test-dbr2.dat" );
//           }

//           std::strstream filename;
//           filename << "test-dbr" << (int)in << ".dat";
//          node->GEM_read_dbr( filename.str() );
          mt.aNodeHandle[in] = in;
          mt.aNodeStatusCH[in] = NEED_GEM_SIA;
          // (8) Setting input data for GEM IPM to use available node speciation as
          // initial approximation
          node->GEM_from_MT( mt.aNodeHandle[in], mt.aNodeStatusCH[in],
                  mt.aT[in], mt.aP[in], mt.aVs[in], mt.aMs[in],
                  mt.abIC[in], mt.adul[in], mt.adll[in], mt.aaPH[in], mt.axDC[in], mt.agam[in] );

          cout << " P " << mt.aP[in] << " T " << mt.aT[in] << endl;
         // cout << "T= " << node->CNode->TK << " P= " << node->CNode->P << endl;
          // (9)   Passing current FMT iteration information into the work DATABR structure
//          node->GEM_set_MT( (double)in, 1. );

//          std::strstream filename;
//          filename << "test-dbr" << (int)in << ".dat";
//          cout << filename.str();
//          node->GEM_write_dbr(filename.str(), false, false, false);
 
          // Calling GEMIPM calculation
          mt.aNodeStatusCH[in] = node->GEM_run( true );

          if( ( mt.aNodeStatusCH[in] == ERR_GEM_AIA || mt.aNodeStatusCH[in] == ERR_GEM_SIA ||
                        mt.aNodeStatusCH[in] ==  T_ERROR_GEM ) )
          {
              cout << "Error: GEM calculation results are not retrieved. Time step";
//                 << it << " node " << in << endl;
          }
          else
          {
            if( ( mt.aNodeStatusCH[in] == BAD_GEM_AIA || mt.aNodeStatusCH[in] == BAD_GEM_SIA  ) )
            {
               cout << "Insufficient quality of GEM solution, but GEM results are retrieved";
//               << it << " node " << in << endl;
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
                  "  Corundum= " << mt.axPH[in][xCorundum] << endl << "moles of Al= " << node->Get_mIC(ICndx[1]) << " P " << mt.aP[in] << " T " << mt.aT[in] << endl;
      }
//  }
  // Calculations finished - end time reached
  // Final output e.g. of total simulation time or of the final distribution of
  //  components and phases in all nodes can be implemented here


  // *********************** DM test of different functions ********************* //

    int alabala = node->DC_name_to_xDB( "AlHSiO3+2"); // gets the index of calcite
    cout << " AlHSiO3+2 G0= " << node->DC_G0(alabala, 150000000, 673.15, false); // prints the G0 of specie in the output file



  // *********************** DM test of different functions ********************* //


  // deleting GEM IPM and data exchange memory structures
  delete node;
  // end of example
  return 0;
}


TMyExperiments::TMyExperiments( long int p_nNod, long int p_nIC, long int p_nDC,
              long int p_nPH, long int p_nPS )
{

    nNodes = p_nNod;
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

TMyExperiments::~TMyExperiments()
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


//void TMyExperiments::OneTimeStepRun(TMyExperiments &mt )
//{

//        long int in;
//        for(  in=1; in< nNodes; in++ )
//        {
//        // some operators that change in some nodes some amounts of some migrating
//        // chemical species (axDC array or  some amounts of migrating chemical
//        // elements in abIC array), possibly using data from other arrays in such
//        // a way that the mass conservation within the whole array of nodes is retained




////            for (int i=0; i<nICndx; i++)
////            {
////                parcel[i] = 0.01* stoicf[i]*abIC[in-1][ICndx[i]];
////                abIC[in][ICndx[i]] += parcel[i];
////    //            if( in > 1 )
////    //                abIC[in-1][ICndx[i]] -= parcel[i];
////            }
//        // The above example loop implements a zero-order flux of MgCl2 in one direction.
//        // Real advective/diffusive transport models are much more complex, but essentially
//        //   do similar things
//        }
//}



//---------------------------------------------------------------------------
// end of main.cpp for node class usage example
