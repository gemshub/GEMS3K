#!/bin/bash
# Installing dependencies needed to build thermofun on (k)ubuntu linux 20.04

if [ "$(uname)" == "Darwin" ]; then

    # Do under Mac OS X platform
    #Needs Xcode and ArangoDB server locally installed
    brew upgrade
    brew install cmake
    brew install spdlog

    EXTN=dylib

elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then

    #Needs gcc v.5 or higher and ArangoDB server locally installed
    sudo apt-get update
    sudo apt-get install libspdlog-dev
    #sudo apt-get install libeigen3-dev
    #sudo ln -s /usr/include/eigen3/Eigen /usr/include/Eigen
    #sudo apt-get install pybind11-dev
    #sudo apt-get install libspdlog-dev

    EXTN=so
fi

# Uncomment what is necessary to reinstall by force 
#sudo rm -f /usr/local/include/nlohmann/json.hpp
#sudo rm -rf /usr/local/include/eigen3/Eigen
#sudo rm -rf /usr/local/include/pybind11
#sudo rm -f /usr/local/lib/libChemicalFun.$EXTN
#sudo rm -f /usr/local/lib/libThermoFun.$EXTN

threads=3
BRANCH_TFUN=master
git status

# nlohmann/json
test -f /usr/local/include/nlohmann/json.hpp || {

	# Building yaml-cpp library
	mkdir -p ~/code && \
        cd ~/code && \
        git clone https://github.com/nlohmann/json.git && \
        cd json && \
        mkdir -p build && \
        cd build && \
        cmake .. -DJSON_BuildTests=OFF && \
        make && \
        sudo make install

	# Removing generated build files
	cd ~ && \
        rm -rf ~/code
}


# Eigen3 math library (added for building and installing xGEMS)
# if not installed in /usr/local/include/eigen3)
test -d /usr/local/include/eigen3/Eigen || {

       # Downloading and unpacking eigen3 source code into ~/code/eigen
       mkdir -p ~/code && \
       cd ~/code && \
       git clone https://gitlab.com/libeigen/eigen.git -b before-git-migration && \
       cd eigen && \
       mkdir -p build && \
       cd build && \
       cmake .. && \
       sudo make install

       # Removing generated build files
       cd ~ && \
       rm -rf ~/code
}

#Pybind11
test -d /usr/local/include/pybind11 || {

        # Building yaml-cpp library
        mkdir -p ~/code && \
        cd ~/code && \
        git clone https://github.com/pybind/pybind11.git && \
        cd pybind11 && \
        mkdir -p build && \
        cd build && \
        cmake .. -DPYBIND11_TEST=OFF && \
        make && \
        sudo make install

        # Removing generated build files
        cd ~ && \
        rm -rf ~/code
}

# ChemicalFun library
# if no ChemicalFun installed in /usr/local/lib/ (/usr/local/include/ChemicalFun)
test -f /usr/local/lib/libChemicalFun.$EXTN || {

        # Building thermofun library
        mkdir -p ~/code && \
        cd ~/code && \
        git clone https://bitbucket.org/gems4/chemicalfun.git -b $BRANCH_TFUN  && \
        cd chemicalfun && \
        mkdir -p build && \
        cd build && \
        cmake .. -DCMAKE_CXX_FLAGS=-fPIC -DCMAKE_BUILD_TYPE=$BUILD_TYPE -DTFUN_BUILD_PYTHON=OFF && \
        make -j $threads && \
        sudo make install

        # Removing generated build files
        cd ~ && \
        rm -rf ~/code
}


# ThermoFun library 
# if no ThermoFun installed in /usr/local/lib/libThermoFun.a (/usr/local/include/ThermoFun)
test -f /usr/local/lib/libThermoFun.$EXTN || {

	# Building thermofun library
	mkdir -p ~/code && \
        cd ~/code && \
        git clone https://bitbucket.org/gems4/thermofun.git -b $BRANCH_TFUN && \
        cd thermofun && \
        mkdir -p build && \
        cd build && \
        cmake .. -DCMAKE_CXX_FLAGS=-fPIC -DCMAKE_BUILD_TYPE=Release && \
        make -j $threads && \
        sudo make install

	# Removing generated build files
	cd ~ && \
        rm -rf ~/code
}

if [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
   sudo ldconfig
fi

