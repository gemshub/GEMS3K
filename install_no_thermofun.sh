#!/bin/bash

#./install-dependencies.sh NO_THERMOFUN

threads=1
BRANCH_GEMS3K=master
BuildType=Release
#BuildType=Debug
InstallPrefix=/usr/local
#InstallPrefix=/home/sveta/devGEMS/gemshub/local
workfolder=${PWD}


mkdir -p build
cd build
cmake .. -DCMAKE_CXX_FLAGS=-fPIC -DCMAKE_BUILD_TYPE=$BuildType -DCMAKE_INSTALL_PREFIX=$InstallPrefix -DUSE_THERMOFUN=OFF
make -j $threads 
sudo make install

if [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
   sudo ldconfig
fi

# cmake ..  -DXGEMS_USE_THERMOFUN=OFF
