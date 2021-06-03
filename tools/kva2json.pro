#  qmake project file for the kva2json example (part of GEMS3K standalone code)
# (c) 2020 GEMS Developer Team
 
TEMPLATE = app
LANGUAGE = C++
TARGET = kva2json
VERSION = 3.4.6

CONFIG -= qt
CONFIG -= warn_on
CONFIG += debug
CONFIG += console
CONFIG += c++17


#DEFINES += IPMGEMPLUGIN
DEFINES += NODEARRAYLEVEL
DEFINES += NOPARTICLEARRAY
#DEFINES += USE_NLOHMANNJSON
DEFINES += OVERFLOW_EXCEPT  #compile with nan inf exceptions

!win32 {

DEFINES += __unix
QMAKE_CFLAGS += -pedantic -Wall -Wextra -Wwrite-strings -Werror

QMAKE_CXXFLAGS += -fPIC -Wall -Wextra -Wformat-nonliteral -Wcast-align -Wpointer-arith \
 -Wmissing-declarations -Winline \ # -Wundef \ #-Weffc++ \
 -Wcast-qual -Wshadow -Wwrite-strings -Wno-unused-parameter \
 -Wfloat-equal -pedantic -ansi


#QMAKE_CXXFLAGS += -fvisibility-inlines-hidden -std=c++17 -fmessage-length=0 -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong \
#-fno-plt -O2 -ffunction-sections -pipe -isystem -O3 -DNDEBUG -fPIC -Wall -Wno-misleading-indentation -Wno-ignored-attributes -Wno-pedantic \
#-Wno-variadic-macros -Wno-deprecated -std=gnu++1z -MD -MT

}

GEMS3K_CPP = ../GEMS3K
GEMS3K_H   = $$GEMS3K_CPP

DEPENDPATH += .
DEPENDPATH += $$GEMS3K_H

INCLUDEPATH += .
INCLUDEPATH += $$GEMS3K_H


QMAKE_LFLAGS +=
#QMAKE_CXXFLAGS += -Wall -Wno-unused
OBJECTS_DIR = obj

include($$GEMS3K_CPP/gems3k.pri) 

HEADERS   +=   args_tool.h
SOURCES   +=   kva2json.cpp

