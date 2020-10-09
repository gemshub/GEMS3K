#  qmake project file for the gemcalc example (part of GEMS3K standalone code)
# (c) 2012-2019 GEMS Developer Team
 
TEMPLATE = app
LANGUAGE = C++
TARGET = tools
VERSION = 3.4.6

CONFIG -= qt
CONFIG -= warn_on
CONFIG += debug
#CONFIG += windows
CONFIG += console
CONFIG += c++17


#DEFINES += IPMGEMPLUGIN
DEFINES += NODEARRAYLEVEL
DEFINES += NOPARTICLEARRAY
#DEFINES += USE_OLD_KV_IO_FILES

!win32 {
  DEFINES += __unix
QMAKE_CFLAGS += pedantic -Wall -Wextra -Wwrite-strings -Werror

QMAKE_CXXFLAGS += -Wall -Wextra -Wformat-nonliteral -Wcast-align -Wpointer-arith \
 -Wmissing-declarations -Winline -Wundef \ #-Weffc++ \
 -Wcast-qual -Wshadow -Wwrite-strings -Wno-unused-parameter \
 -Wfloat-equal -pedantic -ansi

}

GEMS3K_CPP = ../GEMS3K
GEMS3K_H   = $$GEMS3K_CPP

DEPENDPATH +=
DEPENDPATH += .
DEPENDPATH += $$GEMS3K_H

INCLUDEPATH += 
INCLUDEPATH += .
INCLUDEPATH += $$GEMS3K_H

QMAKE_LFLAGS +=
QMAKE_CXXFLAGS += -Wall -Wno-unused
OBJECTS_DIR = obj

SOURCES    +=   calc2convert.cpp

include($$GEMS3K_CPP/gems3k.pri) 

HEADERS += \
    args_impex.h


