//--------------------------------------------------------------------
// $Id$
//
// description
//
// Copyright (c) 2016 A.Yapparova
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
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with GEMS3K code. If not, see <http://www.gnu.org/licenses/>.
//-------------------------------------------------------------------

#include "main.h"

//The case of data exchange in computer memory
int main( int argc, char* argv[] )
{
     // Analyzing command line arguments ( Default arguments)
    char input_system_file_list_name[256] = "system-dat.lst";

    if (argc >= 2 )
        strncpy( input_system_file_list_name, argv[1], 256);

     // Creates TNode structure instance accessible trough the "node" pointer
     TNode* node  = new TNode();

     // (1) Initialization of GEMS3K internal data by reading  files
     //     whose names are given in the ipm_input_system_file_list_name
     if( node->GEM_init( input_system_file_list_name ) )
     {
           cout << "Error occured during reading the files" ;
           return 1;
     }

     double p, T, G0, V0, h, S, cp;
     long int xCH = node->DC_name_to_xCH("H2O@");

     double pmin, pmax, pstep, Tmin, Tmax, Tstep;

     cout<<"\n Enter pmin ";
     cin>>pmin;
     cout<<"\n Enter pmax ";
     cin>>pmax;
     cout<<"\n Enter pstep ";
     cin>>pstep;
     cout<<"\n Enter Tmin ";
     cin>>Tmin;
     cout<<"\n Enter Tmax ";
     cin>>Tmax;
     cout<<"\n Enter Tstep ";
     cin>>Tstep;


     if (pstep <= 0)
         pstep = 1.;
     if (Tstep <= 0)
         Tstep = 1.;

     if (pmin>pmax)
     {
         cout<<"\n pmin>pmax"<<endl;
         return -1;
     }

     if (Tmin>Tmax)
     {
         cout<<"\n Tmin>Tmax"<<endl;
         return -1;
     }

     //pmin = 1.5e5, pmax = 9.5e5, pstep = 1.e5;
     //T = 473.15;

     cout << setprecision(10);
     cout<<"\nT\tp\tG0\tV0\th\tS\tcp";

     for (p=pmin;p<=pmax;p+=pstep)
         for (T=Tmin;T<=Tmax;T+=Tstep)
         {
            G0 = node->DC_G0(xCH, p, T, false);
            V0 = node->DC_V0(xCH, p, T);
            h  = node->DC_H0(xCH, p, T);
            S  = node->DC_S0(xCH, p, T);
            cp = node->DC_Cp0(xCH, p, T);
            cout<<"\n"<<T<<"\t"<<p<<"\t"<<G0<<"\t"<<V0
               <<"\t"<<h<<"\t"<<S<<"\t"<<cp;
         }

     cout<<endl;

  // deleting GEM IPM and data exchange memory structures
  delete node;

  return 0;
}
