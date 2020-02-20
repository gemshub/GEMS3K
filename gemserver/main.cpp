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

//#include <time.h>
//#include <math.h>
//#include <string.h>

#include "node.h"
#include <iomanip>

#include <zmq.hpp>
#include <string>
#include <iostream>
#ifndef _WIN32
#include <sys/stat.h>
#include <unistd.h>
#else
#include <windows.h>

#define sleep(n)    Sleep(n)
#endif

gstring replace( const gstring str, const char* old_part, const char* new_part);
double  CalculateEquilibriumServer( const gstring& lst_f_name );


int main () {

    //  Prepare our context and socket
    zmq::context_t context (1);
    zmq::socket_t socket (context, ZMQ_REP);
    socket.bind ("tcp://*:5555");

    std::cout << "ZeroMQ server start..." << "\nSocket :  " << socket.handle()  << std::endl;

    // create server_data directory
    mkdir("server_data", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    while (true) {
        zmq::recv_flags rsv_flag = zmq::recv_flags::none;
        zmq::message_t request;

        //  Wait for next request from client
        socket.recv (request, rsv_flag);
        string path =  request.to_string();//string( static_cast<const char*>(request.data()));
        std::cout << "Received: " << path << std::endl;

        //  Do some 'work'
        auto time = CalculateEquilibriumServer( path );
        auto stime = std::to_string(time);

        //  Send reply back to client
        zmq::send_flags snd_flags=zmq::send_flags::none;
        zmq::message_t reply(stime.begin(), stime.end());
        //memcpy( reply.data(), "World", 5 );
        socket.send( reply, snd_flags );
    }
    return 0;
}
   

// Run process of calculate equilibria into the GEMS3K side
double  CalculateEquilibriumServer( const gstring& lst_f_name )
{
    double ret=0.;
    try
    {
        auto dbr_lst_f_name  = replace( lst_f_name,"-dat.lst","-dbr-out.dat");

        // Creates TNode structure instance accessible through the "node" pointer
        std::shared_ptr<TNode> node(new TNode());

        // (1) Initialization of GEMS3K internal data by reading  files
        //     whose names are given in the lst_f_name
        if( node->GEM_init( lst_f_name.c_str() ) )
        {
            // error occured during reading the files
            cout << "error occured during reading the files" << endl;
            return ret;
        }

        // (2) re-calculating equilibrium by calling GEMS3K, getting the status back
        long NodeStatusCH = node->GEM_run( false );
        ret  = node->GEM_CalcTime();

        if( NodeStatusCH == OK_GEM_AIA || NodeStatusCH == OK_GEM_SIA  ){
            // (3) Writing results in default DBR file
            node->GEM_write_dbr( nullptr/*dbr_lst_f_name.c_str()*/, false, false, false );
            node->GEM_print_ipm( "GEMipmOK.txt" );   // possible debugging printout
        }
        else {
            // (4) possible return status analysis, error message
            node->GEM_print_ipm( "GEMipmError.txt" );   // possible debugging printout
            return ret; // GEM IPM did not converge properly - error message needed
        }

    }catch(TError& err)
    {


    }
    catch(...)
    {

    }
    return ret;
}

// The method replace() returns a copy of the string
// in which the first occurrence of old have been replaced with new
gstring replace( const gstring str, const char* old_part, const char* new_part)
{
    size_t pos = str.find( old_part );
    if( pos == gstring::npos )
      return str;

    gstring res = str.substr(0, pos);
            res += new_part;
            res += str.substr( pos+strlen(old_part));
    return res;
}


//---------------------------------------------------------------------------
// end of main.cpp for TNode class usage - GEMS3K single calculation example
