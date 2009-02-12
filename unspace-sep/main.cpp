//-------------------------------------------------------------------
// $Id: main.cpp 792 2006-09-19 08:10:41Z gems $
//
// Debugging version of a standalone variant of UnSpace module 
//   extracted from GEMS-PSI v.2.2.2 code
//
//  This standalone program is intended for high-performance
//    computing of large UnSpace tasks (in perspective)
//
// Reads GEMIPM2K input files with setup of the parent chemical 
//    thermodynamic system using TNode class functions from 
//    GEMIPM2K program
//
// Also reads in the file with the UnSpace task setup  
//
// Copyright (C) 2007,2008 S.Dmytrieva, D.Kulik
//
//-------------------------------------------------------------------

#include <time.h>
#include <math.h>
#include "io_arrays.h"
#include "ms_unspace.h"
#include "gstring.h"

istream&
f_getline(istream& is, gstring& str, char delim);

//---------------------------------------------------------------------------

int main( int argc, char* argv[] )
 {
//     int       RetC = 0;
     gstring unspace_in1 = "unspace_init.txt";
     gstring chemsys_in1 = "ipmfiles-dat.lst";
     gstring unspace_out = "UnSpaceTest.out";

// from argv
      if (argc >= 2 )
       chemsys_in1 = argv[1];
      if (argc >= 3 )
    	  unspace_in1 = argv[2];

   try{
cout << "Welcome to Standalone UnSpace program!" << endl; 
    // Allocate the main module here
	 TUnSpace::pu = new TUnSpace();

    // Here we read the MULTI structure, DATACH and DATABR files prepared from GEMS
cout << " * Reading chemical system definition from file name list " << chemsys_in1.c_str() << endl;
	 if( TUnSpace::pu->TaskSystemInit( chemsys_in1.c_str() ) )
	       return 1;  // error reading files

    // Here we read the UnSpace task definition, prepared from GEMS or by hand
cout << " * Reading UnSpace task definition from file " << unspace_in1.c_str() << endl;
	 if( TUnSpace::pu->ReadTask( unspace_in1.c_str() ))
        return 1;  // error reading files

   // here we call the UnSpace sensitivity analysis calculation
cout << " * Starting the UnSpace sensitivity analysis calculations (may take time...)" << endl;
	 TUnSpace::pu->CalcTask( unspace_in1.c_str() );
    
    // Here we write the UnSPace output into text file 
cout << " * Writing UnSpace task calculation results to *.res file "  << unspace_out.c_str() << endl;
	 if( TUnSpace::pu->WriteTask( unspace_out.c_str() ))
        return 1;  // error writing files
cout << "UnSpace success! Bye! " << endl << endl;  
   }
   catch(TError& err)
       {
        fstream f_log("unsplog.txt", ios::out|ios::app );
        f_log << err.title.c_str() << ": " << err.mess.c_str() << endl;
        return 1;
       }

   return 0; // RetC;
}

//---------------------------------------------------------------------------

