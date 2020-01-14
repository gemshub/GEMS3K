#include "node.h"
#include <iomanip>



//
//  Hello World client in C++
//  Connects REQ socket to tcp://localhost:5555
//  Sends "Hello" to server, expects "World" back
//
#include <zmq.hpp>
#include <string>
#include <iostream>

int main ()
{
    //  Prepare our context and socket
    zmq::context_t context (1);
    zmq::socket_t socket (context, ZMQ_REQ);

    std::cout << "Connecting to hello world server…" << std::endl;
    socket.connect ("tcp://localhost:5555");
    gstring path = "server_data/toServer-dat.lst";

    //  Do 5 requests, waiting each time for a response
    for (int request_nbr = 0; request_nbr < 1; request_nbr++)
    {
        zmq::send_flags snd_flags=zmq::send_flags::none;
        zmq::message_t request (path.begin(), path.end());
        std::cout << request_nbr << "Sending:" << request.str() << "…" << std::endl;
        socket.send (request, snd_flags);

        //  Get the reply.
        zmq::recv_flags rsv_flag = zmq::recv_flags::none;
        zmq::message_t reply;
        socket.recv (reply, rsv_flag);
        std::cout << "Received:" << reply.str() << std::endl;
    }
    return 0;
}


//---------------------------------------------------------------------------
// end of main.cpp for TNode class usage - GEMS3K single calculation example
