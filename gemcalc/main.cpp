//--------------------------------------------------------------------
// $Id: main.cpp 686 2012-06-01 14:10:22Z kulik $
//
// Demo test of usage of the TNode class for implementing a simple
// batch-like calculation of equilibria using text file input and
// GEM-IPM-3 numerical kernel.

// TNode class implements a simple C/C++ interface of GEMS3K code.
// It works with DATACH and work DATABR structures and respective
// DCH (chemical system definition) and DBR (recipe or data bridge)
// data files. In addition, the program reads an IPM input file which
// can be used for tuning up numerical controls of GEM IPM-3 algorithm
// and for setting up the parameters of non-ideal mixing models.
//
// Copyright (C) 2007-2012 D.Kulik, S.Dmytriyeva
// <GEMS Development Team, mailto:gems2.support@psi.ch>
//
// This file is part of the GEMS3K code for thermodynamic modelling
// by Gibbs energy minimization <http://gems.web.psi.ch/GEMS3K/>
//
// GEMS3K is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation, either version 3 of
// the License, or (at your option) any later version.

// GEMS3K is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with GEMS3K code. If not, see <http://www.gnu.org/licenses/>
//-------------------------------------------------------------------

#include <time.h>
#include <math.h>
#include <string>
#include <iomanip>
#include "GEMS3K/node.h"
#include "GEMS3K/v_detail.h"
using namespace std;

#define fileNameLength 64
/// Get Path of file and Reading list of file names from it, return number of files
char  (* f_getfiles(const char *f_name, char *Path,
        long int& nElem, char delim ))[fileNameLength];

//The simplest case: data exchange using disk files only
int main( int argc, char* argv[] )
 {
   long nRecipes = 0;
   char (*recipes)[fileNameLength] = nullptr;
   
   // Analyzing command line arguments
   // Default arguments
   char input_system_file_list_name[256] = "system-dat.lst";
   char input_recipes_file_list_name[256] = "more_recipes.lst";
   
   if (argc >= 2 )
       strncpy( input_system_file_list_name, argv[1], 256);
   // list of DCH, IPM and DBR input files for initializing GEMS3K
   
   // Creates TNode structure instance accessible through the "node" pointer
   TNode* node  = new TNode();
  
   // (1) Initialization of GEMS3K internal data by reading  files
   //     whose names are given in the input_system_file_list_name
  if( node->GEM_init( input_system_file_list_name ) )
  {
      // error occured during reading the files
      cout << "error occured during reading the files" << endl;
      return 1;
  }

   // Getting direct access to work node DATABR structure which exchanges the
   // data with GEM IPM3 (already filled out by reading the DBR input file)
   DATABR* dBR = node->pCNode(); 

  // test internal functions
   long int xCa = node->IC_name_to_xCH("Ca");
   long int xCa_ion = node->DC_name_to_xCH("Ca+2");
   long int xCal = node->DC_name_to_xCH("Cal");
   long int xaq = node->Ph_name_to_xCH("aq_gen");
   long int xCalcite = node->Ph_name_to_xCH("Calcite");
   long int xbCa = node->IC_xCH_to_xDB(xCa);
   long int xbCa_ion = node->DC_xCH_to_xDB(xCa_ion);
   long int xbCal = node->DC_xCH_to_xDB(xCal);
   long int xbaq = node->Ph_xCH_to_xDB(xaq);
   long int xbCalcite = node->Ph_xCH_to_xDB(xCalcite);
   
   cout << "          CH  BR" << endl;
   cout << " Ca       " << xCa << "   " << xbCa << endl;
   cout << " Ca+2     " << xCa_ion << "   " << xbCa_ion << endl;
   cout << " Cal     " << xCal << "  " << xbCal << endl;
   cout << " aq_gen   " << xaq << "   " << xbaq << endl;
   cout << " Calcite  " << xCalcite << "   " << xbCalcite << endl;
   cout << setprecision(7) << setw(10) << endl;
   
   // Asking GEM to run with automatic initial approximation 
   dBR->NodeStatusCH = NEED_GEM_AIA;

//node->Ph_Enthalpy(0);
//node->Ph_Enthalpy(1);
//node->Ph_Enthalpy(2);
//node->Ph_Enthalpy(3);

   node->GEM_print_ipm( "BeforeCalcPhase.txt" );   // possible debugging printout

// (2) re-calculating equilibrium by calling GEMS3K, getting the status back
   long NodeStatusCH = node->GEM_run( false );

   if( NodeStatusCH == OK_GEM_AIA || NodeStatusCH == OK_GEM_SIA  )
   {    // (3) Writing results in default DBR file
       node->GEM_write_dbr( nullptr, false, true, false );
node->GEM_print_ipm( "AfterCalcPhase.txt" );   // possible debugging printout
   }
   else {
      // (4) possible return status analysis, error message
       node->GEM_print_ipm( nullptr );   // possible debugging printout
       return 5; // GEM IPM did not converge properly - error message needed
   }
   cout << "SatIndx aq" << node->Ph_SatInd(0) << " gas " << node->Ph_SatInd(1) << " s1 " << node->Ph_SatInd(2)
        << " s2 " << node->Ph_SatInd(3) << " s3 " << node->Ph_SatInd(4) << " s4 " << node->Ph_SatInd(5) << endl;
return 0;
   // test internal functions
  cout << "Ph_Volume   Aq: " << node->Ph_Volume(xbaq) <<  " Calcite: " << node->Ph_Volume(xbCalcite) << endl;   
  cout << "Ph_Mass     Aq: " << node->Ph_Mass(xbaq) <<  " Calcite: " << node->Ph_Mass(xbCalcite) << endl;   
  cout << "Ph_SatInd   Aq: " << node->Ph_SatInd(xbaq) <<  " Calcite: " << node->Ph_SatInd(xbCalcite) << endl;   
   
  cout << endl;
  cout << "Ca+2    Get_nDC  " << node->Get_nDC(xbCa_ion) <<  " DC_n  " << node->DC_n(xCa_ion) << endl;   
  cout << "Cal     Get_nDC  " << node->Get_nDC(xbCal) <<  " DC_n  " << node->DC_n(xCal) << endl;   
  cout << "Ca+2    Get_muDC " << node->Get_muDC(xbCa_ion) <<  " DC_mu " << node->DC_mu(xCa_ion) << endl;   
  cout << "Cal     Get_muDC " << node->Get_muDC(xbCal) <<  " DC_mu " << node->DC_mu(xCal) << endl;   
  cout << "Ca+2    Get_aDC  " << node->Get_aDC(xbCa_ion) <<  " DC_a  " << node->DC_a(xCa_ion) << endl;   
  cout << "Cal     Get_aDC  " << node->Get_aDC(xbCal) <<  " DC_a  " << node->DC_a(xCal) << endl;   
  cout << "Ca+2    Get_cDC  " << node->Get_cDC(xbCa_ion) <<  " DC_c  " << node->DC_c(xCa_ion) << endl;   
  cout << "Cal     Get_cDC  " << node->Get_cDC(xbCal) <<  " DC_c  " << node->DC_c(xCal) << endl;   
  cout << "Ca+2    Get_gDC  " << node->Get_gDC(xbCa_ion) <<  " DC_g  " << node->DC_g(xCa_ion) << endl;   
  cout << "Cal     Get_gDC  " << node->Get_gDC(xbCal) <<  " DC_g  " << node->DC_g(xCal) << endl;   
  cout << endl;
  cout << "G0   Ca+2: " << node->DC_G0( xCa_ion, node->cP(), node->cTK(), false ) <<  " Cal: " << node->DC_G0( xCal, node->cP(), node->cTK(), false ) << endl;   
  cout << "V0   Ca+2: " << node->DC_V0( xCa_ion, node->cP(), node->cTK() ) <<  " Cal: " << node->DC_V0( xCal, node->cP(), node->cTK() ) << endl;   
  cout << "H0   Ca+2: " << node->DC_H0( xCa_ion, node->cP(), node->cTK() ) <<  " Cal: " << node->DC_H0( xCal, node->cP(), node->cTK() ) << endl;   
  cout << "S0   Ca+2: " << node->DC_S0( xCa_ion, node->cP(), node->cTK() ) <<  " Cal: " << node->DC_S0( xCal, node->cP(), node->cTK() ) << endl;   
  cout << "Cp0  Ca+2: " << node->DC_Cp0( xCa_ion, node->cP(), node->cTK() ) <<  " Cal: " << node->DC_Cp0( xCal, node->cP(), node->cTK() ) << endl;   
  cout << "A0   Ca+2: " << node->DC_A0( xCa_ion, node->cP(), node->cTK() ) <<  " Cal: " << node->DC_A0( xCal, node->cP(), node->cTK() ) << endl;   
  cout << "U0   Ca+2: " << node->DC_U0( xCa_ion, node->cP(), node->cTK() ) <<  " Cal: " << node->DC_U0( xCal, node->cP(), node->cTK() ) << endl;   
   
   // Here a possible loop on more input recipes begins
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
        sprintf(NextRecipeFileName , "%s\%s", input_recipes_file_list_path, recipes[cRecipe] );

        // (5) Reading the next DBR file with different input composition or temperature
        node->GEM_read_dbr( NextRecipeFileName );

        // Asking GEM IPM2 to run (faster) with smart initial approximation 
        dBR->NodeStatusCH = NEED_GEM_SIA;       

        NodeStatusCH = node->GEM_run( false );

        if( NodeStatusCH == OK_GEM_AIA || NodeStatusCH == OK_GEM_SIA  )
        {    sprintf(NextRecipeOutFileName , "%s.nc.out", NextRecipeFileName );
             node->GEM_write_dbr( NextRecipeOutFileName, false, false, false );
             sprintf(NextRecipeOutFileName , "%s.nc.Dump.out", NextRecipeFileName );
             node->GEM_print_ipm( NextRecipeOutFileName );
        }
        else {
               // error message, debugging printout
     	      sprintf(NextRecipeOutFileName , "%s.Dump.out", NextRecipeFileName );
              node->GEM_print_ipm( NextRecipeOutFileName );
//              return 5; // GEM IPM did not converge properly - error message needed
              }
      }
   }	 
   // end of possible loop on input recipes
   delete node;
   if( recipes ) delete[] recipes;

 // End of example  
   return 0; 
}
   
const int bGRAN = 20;

// Get Path of file and Reading list of file names from it, return number of files
char  (* f_getfiles(const char *f_name, char *Path,
        long int& nElem, char delim ))[fileNameLength]
{
  int ii, bSize = bGRAN;
  char  (*filesList)[fileNameLength];
  char  (*filesListNew)[fileNameLength];
  filesList = new char[bSize][fileNameLength];
  string name;

// Get path
   string path_;
   string flst_name = f_name;
   size_t pos = flst_name.rfind("/");
   path_ = "";
   if( pos < string::npos )
      path_ = flst_name.substr(0, pos+1);
   strncpy( Path, path_.c_str(), 256-fileNameLength);
   Path[255] = '\0';

//  open file stream for the file names list file
   fstream f_lst( f_name/*flst_name.c_str()*/, ios::in );
   ErrorIf( !f_lst.good(), f_name, "Fileopen error");

// Reading list of names from file
  nElem = 0;
  while( !f_lst.eof() )
  {
    f_getline( f_lst, name, delim);
    if( nElem >= bSize )
    {    bSize = bSize+bGRAN;
         filesListNew = new char[bSize][fileNameLength];
         for( ii=0; ii<nElem-1; ii++ )
           strncpy( filesListNew[ii], filesList[ii], fileNameLength);
         delete[] filesList;
         filesList =  filesListNew;
    }
    strncpy( filesList[nElem], name.c_str(), fileNameLength);
    filesList[nElem][fileNameLength-1] = '\0';
    nElem++;
  }

  // Realloc memory for reading size
  if( nElem != bSize )
  {
    filesListNew = new char[nElem][fileNameLength];
    for(  ii=0; ii<nElem; ii++ )
      strncpy( filesListNew[ii], filesList[ii], fileNameLength);
    delete[] filesList;
    filesList =  filesListNew;
  }

  return filesList;
}


//---------------------------------------------------------------------------
// end of main.cpp for TNode class usage - GEMS3K single calculation example
