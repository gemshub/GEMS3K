#include <zmq.hpp>
#include <zmq_addon.hpp>
#include <string>
#include <iostream>
#include "tnodetask.h"

const char* one_system_task = "system";
const char* only_dbr_task = "dbr";
const char* nodearray_task = "nodearray";

int main ()
{
    NodeGEMSTask task_data;

    //  Prepare our context and socket
    zmq::context_t context (1);
    zmq::socket_t socket (context, ZMQ_REP);
    socket.bind ("tcp://*:5555");

    std::cout << "ZeroMQ server start... "  << socket.handle()  << std::endl;

    while (true) {

        std::vector<std::string> ret_msgs;

        // get input data
        std::vector<zmq::message_t> omsgs;
        /*auto oret =*/ zmq::recv_multipart(socket, std::back_inserter(omsgs));
        std::vector<std::string> msgs_data;
        for( const auto& msg: omsgs )
            msgs_data.push_back(  msg.to_string() );

        // execute command
        if(  msgs_data.size() >= 4 && msgs_data[0] == one_system_task)
        {
            //std::cout << "Init... "  << msgs_data.size() << std::endl;
            ret_msgs =  task_data.initData( msgs_data[1], msgs_data[2], msgs_data[3] );
            if( ret_msgs.empty() )
                ret_msgs = task_data.calculateEquilibrium("");
        }
        else if(  msgs_data.size() >= 4 && msgs_data[0] == nodearray_task )
        {
            //std::cout << "gem2mt... "  << msgs_data.size() << std::endl;
            ret_msgs =  task_data.initData( msgs_data[1], msgs_data[2], msgs_data[3] );
            if( ret_msgs.empty() )
            {
                //task_data.calculateEquilibrium("");
                ret_msgs.push_back("ok");
            }
        }
        else if(  msgs_data.size() >= 2 && msgs_data[0] == only_dbr_task )
        {
            //std::cout << "dbr... "  << msgs_data.size() << std::endl;
            ret_msgs = task_data.calculateEquilibrium(msgs_data[1]);
        }

        // send result
        std::vector<zmq::message_t> msgs_vec;
        for( const auto& msg: ret_msgs )
            msgs_vec.push_back(  zmq::message_t(msg.begin(), msg.end()) );
        /*auto iret =*/ zmq::send_multipart(socket, msgs_vec);

    }
    return 0;
}



//---------------------------------------------------------------------------
