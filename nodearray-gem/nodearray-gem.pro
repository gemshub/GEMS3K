#  qmake project file for the nodearray-gem example (part of GEMS3K standalone code)
# (c) 2012-2019 GEMS Developer Team

TEMPLATE = app
LANGUAGE = C++
TARGET = nodearrs
VERSION = 3.4.6

CONFIG -= qt
CONFIG += warn_on
CONFIG += thread console
CONFIG += c++17
CONFIG += sanitaze sanitaze_thread

#DEFINES += NO_NODEARRAYLEVEL
#DEFINES += USE_NLOHMANNJSON
#DEFINES += useOMP
#DEFINES += NO_USE_THERMOFUN
#DEFINES += USE_THERMO_LOG
DEFINES += OVERFLOW_EXCEPT  #compile with nan inf exceptions

!win32 {
  DEFINES += __unix
}

GEMS3K_CPP = ../GEMS3K
GEMS3K_H   = $$GEMS3K_CPP

DEPENDPATH +=
DEPENDPATH += .
DEPENDPATH += $$GEMS3K_H

INCLUDEPATH += .
INCLUDEPATH += $$GEMS3K_H

#QMAKE_CXXFLAGS += -fopenmp -pg
#QMAKE_CXXFLAGS += -fopenmp
#QMAKE_LFLAGS +=  -fopenmp
#LIBS += -lgomp -lpthread

#QMAKE_LFLAGS += -pg
OBJECTS_DIR = obj

!contains(DEFINES, NO_USE_THERMOFUN) {
LIBS += -lThermoFun -lChemicalFun
} ## end NO_USE_THERMOFUN

include($$GEMS3K_CPP/gems3k.pri)

HEADERS	 +=  m_gem2mt.h \
             particlearray.h

SOURCES  +=   main.cpp \
              particlearray.cpp \
              m_gem2mtsep.cpp \
              m_gem2mtt.cpp \
              m_gem2mtbox.cpp \
              m_gem2mtfor.cpp
