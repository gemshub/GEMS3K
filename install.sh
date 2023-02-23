#!/bin/bash

#./install-dependencies.sh

threads=1
BRANCH_GEMS3K=master
BuildType=Release
#BuildType=Debug
InstallPrefix=/usr/local
workfolder=${PWD}

mkdir -p build
cd build
cmake .. -DCMAKE_CXX_FLAGS=-fPIC -DCMAKE_BUILD_TYPE=$BuildType -DCMAKE_INSTALL_PREFIX=$InstallPrefix
make -j $threads 
sudo make install

if [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
   sudo ldconfig
fi
