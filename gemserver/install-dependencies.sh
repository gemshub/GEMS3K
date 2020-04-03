#!/bin/bash
# Installing dependencies needed to build GEMS-Reactoro for linux 16.04 or 18.04

#sudo rm -f /usr/local/lib/libzmq.a
#sudo rm -f /usr/local/include/zmq.hpp

threads=3
workfolder=${PWD}

# ZeroMQ core engine in C++, implements ZMTP/3.1 (https://github.com/zeromq/libzmq)
# if no zeromq installed in /usr/local/lib/zeromq.a (/usr/local/include/zeromq)
test -f /usr/local/lib/libzmq.a || {

        # Building zeromq library
        mkdir -p ~/code && \
                cd ~/code && \
                git clone https://github.com/zeromq/libzmq.git && \
                cd libzmq && \
                mkdir -p build && \
                cd build && \
                cmake .. -DCMAKE_CXX_FLAGS=-fPIC -DCMAKE_BUILD_TYPE=Release && \
                make -j $threads && \
                sudo make install

        # Removing generated build files
        cd ~ && \
                 rm -rf ~/code
}

# Header-only C++ binding for libzmq (https://github.com/zeromq/cppzmq)
# zeromq/cppzmq
test -f /usr/local/include/zmq.hpp || {

	# Building zmq library
	mkdir -p ~/code && \
		cd ~/code && \
                git clone https://github.com/zeromq/cppzmq.git && \
                cd cppzmq && \
		mkdir -p build && \
		cd build && \
                cmake ..  && \
		make && \
		sudo make install

	# Removing generated build files
	cd ~ && \
		 rm -rf ~/code
}

sudo ldconfig
