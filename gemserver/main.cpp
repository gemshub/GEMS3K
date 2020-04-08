
#define OLD


#include <iomanip>

#include <zmq.hpp>
#include <zmq_addon.hpp>
#include <string>
#include <iostream>
#ifndef _WIN32
#include <sys/stat.h>
#include <unistd.h>
#else
#include <windows.h>

#define sleep(n)    Sleep(n)
#endif

double  CalculateEquilibriumServer( const std::string& lst_f_name );
/// Run process of calculate equilibria into the GEMS3K side (read from strings)
double  CalculateEquilibriumServer( const std::vector<std::string>& dch_ipm_dbr, std::string& dbr_result );



int main () {

    //  Prepare our context and socket
    zmq::context_t context (1);
    zmq::socket_t socket (context, ZMQ_REP);
    socket.bind ("tcp://*:5555");

    std::cout << "ZeroMQ server start... "  << socket.handle()  << std::endl;

#ifdef OLD
    // create server_data directory
    mkdir("server_data", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
#endif

    while (true) {

#ifdef OLD
        zmq::recv_flags rsv_flag = zmq::recv_flags::none;
        zmq::message_t request;

        //  Wait for next request from client
        socket.recv (request, rsv_flag);
        std::string path =  request.to_string();//string( static_cast<const char*>(request.data()));
        std::cout << "Received: " << path << std::endl;

        //  Do some 'work'
        auto time = CalculateEquilibriumServer( path );
        auto stime = std::to_string(time);

        //  Send reply back to client
        zmq::send_flags snd_flags=zmq::send_flags::none;
        zmq::message_t reply(stime.begin(), stime.end());
        //memcpy( reply.data(), "World", 5 );
        socket.send( reply, snd_flags );

#else
       std::vector<zmq::message_t> omsgs;
       auto oret = zmq::recv_multipart(socket, std::back_inserter(omsgs));

       auto stime = std::string("10");
       std::string dbr_result;
       std::vector<std::string> msgs_data;
       for( const auto& msg: omsgs )
           msgs_data.push_back(  msg.to_string() );

       if(  msgs_data.size() >= 4 and msgs_data[0] == "system")
       {
         // one system calculation
           auto time = CalculateEquilibriumServer( msgs_data, dbr_result );
           stime = std::to_string(time);
       }

       std::vector<zmq::message_t> msgs_vec;
       msgs_vec.push_back( zmq::message_t(stime.begin(), stime.end()));
       msgs_vec.push_back( zmq::message_t(dbr_result.begin(), dbr_result.end()));
       auto iret = zmq::send_multipart(socket, msgs_vec);

#endif

    }
    return 0;
}
   


//---------------------------------------------------------------------------
