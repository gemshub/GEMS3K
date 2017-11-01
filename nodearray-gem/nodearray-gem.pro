#  qmake project file for the nodearray-gem example (part of GEMS3K standalone code)
# (c) 2012 GEMS Developer Team

TEMPLATE = app
LANGUAGE = C++
TARGET = nodearrs
VERSION = 3.0.0

CONFIG -= qt
CONFIG += warn_on
CONFIG += debug
#CONFIG += windows
CONFIG += console

DEFINES += IPMGEMPLUGIN
DEFINES += NODEARRAYLEVEL
#DEFINES += NOPARTICLEARRAY
#DEFINES += SEPGEM2MTMODE

!win32 {
  DEFINES += __unix
}

GEMS3K_CPP = ../GEMS3K
GEMS3K_H   = $$GEMS3K_CPP

DEPENDPATH +=
DEPENDPATH += .
DEPENDPATH += $$GEMS3K_H

INCLUDEPATH +=
INCLUDEPATH += .
INCLUDEPATH += $$GEMS3K_H

QMAKE_CXXFLAGS += -fopenmp
QMAKE_LFLAGS +=  -fopenmp
LIBS += -lgomp -lpthread

QMAKE_LFLAGS +=
OBJECTS_DIR = obj


HEADERS	 +=  m_gem2mt.h \
             particlearray.h

SOURCES  +=   main.cpp \
              particlearray.cpp \
              m_gem2mtt.cpp \
              m_gem2mtbox.cpp \
              m_gem2mtfor.cpp \
              m_gem2mtsep.cpp

include($$GEMS3K_CPP/gems3k.pri)
