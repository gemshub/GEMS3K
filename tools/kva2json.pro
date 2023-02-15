#  qmake project file for the kva2json example (part of GEMS3K standalone code)
# (c) 2020 GEMS Developer Team
 
TEMPLATE = app
LANGUAGE = C++
TARGET = kva2json
VERSION = 3.4.6

CONFIG -= qt
CONFIG += warn_on
CONFIG += thread console
CONFIG += c++17
CONFIG += sanitaze sanitaze_thread

DEFINES += NODEARRAYLEVEL
#DEFINES += USE_NLOHMANNJSON
DEFINES += USE_THERMOFUN
DEFINES += USE_THERMO_LOG
DEFINES += OVERFLOW_EXCEPT  #compile with nan inf exceptions

!win32 {

DEFINES += __unix
QMAKE_CFLAGS += -pedantic -Wall -Wextra -Wwrite-strings -Werror

QMAKE_CXXFLAGS += -fPIC -Wall -Wextra -Wformat-nonliteral -Wcast-align -Wpointer-arith \
 -Wmissing-declarations -Winline \ # -Wundef \ #-Weffc++ \
 -Wcast-qual -Wshadow -Wwrite-strings -Wno-unused-parameter \
 -Wfloat-equal -pedantic -ansi

}

GEMS3K_CPP = ../GEMS3K
GEMS3K_H   = $$GEMS3K_CPP

DEPENDPATH += .
DEPENDPATH += $$GEMS3K_H

INCLUDEPATH += .
INCLUDEPATH += $$GEMS3K_H

contains(DEFINES, USE_THERMOFUN) {

#ThermoFun_CPP   =  ../ThermoFun
#ThermoFun_H     =   $$ThermoFun_CPP
#DEPENDPATH += $$ThermoFun_H
#INCLUDEPATH += $$ThermoFun_H
#include($$ThermoFun_CPP/ThermoFun.pri)
LIBS += -lThermoFun -lChemicalFun

} ## end USE_THERMOFUN


QMAKE_LFLAGS +=
#QMAKE_CXXFLAGS += -Wall -Wno-unused
OBJECTS_DIR = obj

include($$GEMS3K_CPP/gems3k.pri) 

#HEADERS   +=   args_tool.h
#SOURCES   +=   kva2json.cpp

SOURCES   +=   thread_test.cpp
