//--------------------------------------------------------------------
// $Id: main.cpp 1386 2009-08-06 12:43:51Z gems $
//
// Demo test of usage of the TNode class for implementing a simple
// batch-like calculation of equilibria using text file input and
// GEMIPM2 numerical kernel

// TNode class implements a  simple C/C++ interface of GEMIPM2K.
// It works with DATACH and work DATABR structures and respective
// DCH (chemical system definition) and DBR (recipe or data bridge)
// data files. In addition, the program reads an IPM inlut file which 
// can be used for tuning up numerical controls of GEM IPM2 algorithm 
// and for setting up the parameters of non-ideal mixing models.
//
// Copyright (C) 2007,2008 D.Kulik, S.Dmitrieva
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
#include <iomanip>


//The simplest case: data exchange using disk files only
int main( int argc, char* argv[] )
 {
   long nRecipes = 0;
   char (*recipes)[fileNameLength] = 0;
   
   // Analyzing command line arguments
   // Default arguments
   char input_system_file_list_name[256] = "system-dat.lst";
   char input_recipes_file_list_name[256] = "more_recipes.lst";
   
   if (argc >= 2 )
       strncpy( input_system_file_list_name, argv[1], 256);
   // list of DCH, IPM and DBR input files for initializing GEMIPM2K
   
   // Creates TNode structure instance accessible trough the "node" pointer
   TNode* node  = new TNode();
  
   // (1) Initialization of GEMIPM2K internal data by reading  files
   //     whose names are given in the input_system_file_list_name
  if( node->GEM_init( input_system_file_list_name ) )
  {
      // error occured during reading the files
      return 1;
  }

   // Getting direct access to work node DATABR structure which exchanges the
   // data with GEM IPM2 (already filled out by reading the DBR input file)
   DATABR* dBR = node->pCNode(); 

   
   // Asking GEM to run with automatic initial approximation 
   dBR->NodeStatusCH = NEED_GEM_AIA;

   // (2) re-calculating equilibrium by calling GEMIPM2K, getting the status back
   int NodeStatusCH = node->GEM_run( false );

   if( NodeStatusCH == OK_GEM_AIA || NodeStatusCH == OK_GEM_SIA  )
   {    // (3) Writing results in default DBR file
       node->GEM_write_dbr( NULL, false, true );
//       node->GEM_print_ipm( NULL );   // possible debugging printout
   }
   else {
      // (4) possible return status analysis, error message
       node->GEM_print_ipm( NULL );   // possible debugging printout
       return 5; // GEM IPM did not converge properly      
        }

   // Here a possible loop on input recipes begins
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
  	  
  	 for(int cRecipe=0; cRecipe < nRecipes; cRecipe++ )
      { 
         // Trying to read the next file name 
  	    sprintf(NextRecipeFileName , "%s%s", input_recipes_file_list_path, recipes[cRecipe] );
  
        // (5) Reading the next DBR file with different input composition or temperature
        node->GEM_read_dbr( NextRecipeFileName );

        // Asking GEM IPM2 to run (faster) with smart initial approximation 
        dBR->NodeStatusCH = NEED_GEM_SIA;       

        NodeStatusCH = node->GEM_run( false );

        if( NodeStatusCH == OK_GEM_AIA || NodeStatusCH == OK_GEM_SIA  )
        {    
            sprintf(NextRecipeOutFileName , "%s.out", NextRecipeFileName );
            node->GEM_write_dbr( NextRecipeOutFileName, false, true );

        }
        else {
               // error message, debugging printout
     	      sprintf(NextRecipeOutFileName , "%s.Dump.out", NextRecipeFileName );
              node->GEM_print_ipm( NextRecipeOutFileName );
               // GEM IPM did not converge properly
             }

      }
   }	 
  	 // end of possible loop on input recipes
   delete node;
   if( recipes ) delete recipes;

 // End of example  
   return 0; 
 }  
   

//---------------------------------------------------------------------------
// end of main.cpp for TNode class usage - GEM single calculation example
