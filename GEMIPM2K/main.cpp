//--------------------------------------------------------------------
// $Id: main.cpp 1004 2008-01-23 14:06:44Z gems $
//
// Demo test of usage of the TNode class for implementing a simple
// batch-like calculation of equilibria using text file input and
// GEMIPM2 numerical kernel

// TNode class implements a  simple C/C++ interface of GEMIPM2K.
// It works with DATACH and work DATABR structures and respective
// .dch.dat (chemical system definition) and .dbr.dat (recipe or data bridge)
// files. In addition, the program reads a .ipm.dat file which can be
// used for tuning up numerical controls of GEM IPM2 algorithm and for
// setting up the parameters of non-ideal mixing models.
//
// Copyright (C) 2007 D.Kulik, S.Dmitrieva
//
// This file is part of GEMIPM2K code for thermodynamic modelling
// by Gibbs energy minimization

// This file may be distributed under the licence terms defined
// in GEMIPM2K.QAL
//
// See also http://gems.web.psi.ch
// mailto://gems2.support@psi.ch
//-------------------------------------------------------------------

#include <time.h>
#include <math.h>
#include <string.h>

#include "node.h"

int main( int argc, char* argv[] )
 {
   double PureTime = 0.0;
   int nNodes = 1, nRecipes = 0;
   char  (*recipes)[fileNameLength] = 0;

   //  char *recipes = NULL;
 //  char *ThisRecipe = NULL;

   // Analyzing command line arguments
     // Default arguments
     char input_system_file_list_name[256] = "system.lst";
     char input_recipes_file_list_name[256] = "more_recipes.lst";
   
     if (argc >= 2 )
       strncpy( input_system_file_list_name, argv[1], 256);
         // list of CSD files needed as input for initializing GEMIPM2K

   // Creating TNode structure accessible trough the "node" pointer
   TNode* node  = new TNode();

   cout << "Welcome to GEMIPM2K v. 2.2.0 solver of (geo)chemical equilibria! "
        << endl;

   // Here we read the files needed as input for initializing GEMIPM2K
   // The easiest way to prepare them is to use GEMS-PSI code (GEM2MT module)
   if( node->GEM_init( input_system_file_list_name ) )
   {
     cout << "Error: Chemical system definition files listed in " << endl
          << input_system_file_list_name << endl
          << " cannot be found or are corrupt" << endl;
      return 1;  // error occured during reading the files
   }

  // Number of ICs, DCs, Phases and Phases-solutions kept in the node
  // DATABR structure for exchange with GEMIPM - for your convenience
  int nIC, nDC, nPH, nPS;
  int i,   j,   k,   ks;    // indices for direct access to components
                            // and phases data in the DataCH framework
  short nodeHandle, NodeStatusCH, IterDone;

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
  
  //
  // Calculate start DATABR structure reading from files needed as input for initializing GEMIPM2K
  //
  dBR->NodeStatusFMT = No_transport; 
  dBR->NodeStatusCH = NEED_GEM_AIA; // direct access to node DATABR structure
  dBR->NodeHandle = -1;
  // re-calculating equilibrium by calling GEMIPM2K, getting the status
  NodeStatusCH = node->GEM_run( false );
  PureTime += node->GEM_CalcTime();

  if( !( NodeStatusCH == OK_GEM_AIA || NodeStatusCH == OK_GEM_PIA ) )
      {
          cout << "Error: GEMIPM2 did not converge well " << NodeStatusCH << endl;
          cout << "See *.Dump.out file" << endl;
          node->GEM_print_ipm( NULL );
          cout << "Bye! " << endl;
          return 5; // GEM IPM did not converge properly
      }

  IterDone = dBR->IterDone;
  cout << "GEMIPM2 did " << IterDone << " iterations in " << node->GEM_CalcTime() 
        << " sec in mode " << NodeStatusCH << "." << endl;
  cout << "  See output in the *.out file" << endl;
  node->GEM_write_dbr( NULL, false, true );
//  cout << "See details in the *.Dump.out file" << endl;
//  node->GEM_print_ipm( NULL );


  // Working with list of additional recipes
  if (argc >= 3 )
  {  
	  char NextRecipeFileName[256];
	  char NextRecipeOutFileName[300];
	  char input_recipes_file_list_path[256-fileNameLength] = "";
    
	  strncpy( input_recipes_file_list_name, argv[2], 256);
      // list of additional recipes (dbr.dat files) e.g. for simulation
      // of a titration or another irreversible process
 	 // Reading list of recipes names from file 
 	  recipes = f_getfiles(  input_recipes_file_list_name, 
 	       		 input_recipes_file_list_path, nRecipes, ',');			 
     if( nRecipes < 0 || nRecipes > 1000 )
    	 goto FINISH;
     
 	  
 	 for(int cRecipe=0; cRecipe < nRecipes; cRecipe++ )
     { 
        // Trying to read the next file name 
 	    sprintf(NextRecipeFileName , "%s%s", input_recipes_file_list_path, recipes[cRecipe] );
 	    TNode::na->GEM_read_dbr( NextRecipeFileName );
 	    
 	    dBR->NodeStatusFMT = No_transport; 
 		dBR->NodeStatusCH = NEED_GEM_PIA; // direct access to node DATABR structure
 		dBR->NodeHandle = cRecipe;
        // re-calculating equilibrium by calling GEMIPM2K, getting the status
        NodeStatusCH = node->GEM_run( false );
        PureTime += node->GEM_CalcTime();

        if( !( NodeStatusCH == OK_GEM_AIA || NodeStatusCH == OK_GEM_PIA ) )
        {
          cout << "Error: GEMIPM2 did not converge well " << NodeStatusCH << endl;
   	      sprintf(NextRecipeOutFileName , "%s.Dump.out", NextRecipeFileName );
          cout << "See " << NextRecipeOutFileName << " file" << endl;
          node->GEM_print_ipm( NextRecipeOutFileName );
          cout << "Bye! " << endl;
          return 5; // GEM IPM did not converge properly
        }

      IterDone = dBR->IterDone;
      cout << "GEMIPM2 did " << IterDone << " iterations in " << node->GEM_CalcTime() 
            << " sec in mode " << NodeStatusCH << "." << endl;
      sprintf(NextRecipeOutFileName , "%s.out", NextRecipeFileName );
      cout << "  See output in the " << NextRecipeOutFileName << " file" << endl;
      node->GEM_write_dbr( NextRecipeOutFileName, false, false );
//      sprintf(NextRecipeOutFileName , "%s.Dump.out", NextRecipeFileName );
//      cout << "See dump output in the " << NextRecipeOutFileName << " file" << endl;
//      node->GEM_print_ipm( NextRecipeOutFileName  );
    
   }  // end for

  } // end if
  
FINISH:
  // Finished! 
  // deleting GEMIPM and data exchange memory structures
   delete node;
   if( recipes ) delete recipes;
   cout <<  "Pure time of calculation, s: " <<  PureTime << endl;
   cout << "Bye!" << endl;

   return 0;

 }

//---------------------------------------------------------------------------
// end of main.cpp for node class usage - GEM single calculation example
