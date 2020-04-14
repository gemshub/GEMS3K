//#define OLD

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
std::vector<std::string>  CalculateEquilibriumServer( const std::vector<std::string>& dch_ipm_dbr );




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
       std::vector<std::string> ret_msgs;

       // get input data
       std::vector<zmq::message_t> omsgs;
       /*auto oret =*/ zmq::recv_multipart(socket, std::back_inserter(omsgs));
       std::vector<std::string> msgs_data;
       for( const auto& msg: omsgs )
           msgs_data.push_back(  msg.to_string() );

       // execute command
       if(  msgs_data.size() >= 4 and msgs_data[0] == "system")
       {
         // one system calculation
           ret_msgs = CalculateEquilibriumServer( msgs_data );
       }

       // send result
       std::vector<zmq::message_t> msgs_vec;
       for( const auto& msg: ret_msgs )
           msgs_vec.push_back(  zmq::message_t(msg.begin(), msg.end()) );
       /*auto iret =*/ zmq::send_multipart(socket, msgs_vec);

#endif

    }
    return 0;
}
   


//---------------------------------------------------------------------------
