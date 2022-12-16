#  qmake project file for the nodearray-gem example (part of GEMS3K standalone code)
# (c) 2012-2019 GEMS Developer Team

TEMPLATE = app
LANGUAGE = C++
TARGET = nodearrs
VERSION = 3.4.6

CONFIG -= qt
CONFIG += warn_on
CONFIG += debug
CONFIG += thread
CONFIG += console
CONFIG += c++17
CONFIG += sanitaze sanitaze_thread

DEFINES += NODEARRAYLEVEL
#DEFINES += useOMP
#DEFINES += USE_NLOHMANNJSON

!win32 {
  DEFINES += __unix
}

GEMS3K_CPP = ../GEMS3K
GEMS3K_H   = $$GEMS3K_CPP

DEPENDPATH +=
DEPENDPATH += .
DEPENDPATH += $$GEMS3K_H

INCLUDEPATH += /usr/lib/gcc/x86_64-linux-gnu/9/include
INCLUDEPATH += .
INCLUDEPATH += $$GEMS3K_H

#QMAKE_CXXFLAGS += -fopenmp -pg
QMAKE_CXXFLAGS += -fopenmp
QMAKE_LFLAGS +=  -fopenmp
LIBS += -lgomp -lpthread

#QMAKE_LFLAGS += -pg
OBJECTS_DIR = obj

include($$GEMS3K_CPP/gems3k.pri)

HEADERS	 +=  m_gem2mt.h \
             particlearray.h

SOURCES  +=   main.cpp \
              particlearray.cpp \
              m_gem2mtsep.cpp \
              m_gem2mtt.cpp \
              m_gem2mtbox.cpp \
              m_gem2mtfor.cpp
