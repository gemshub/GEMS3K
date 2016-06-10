#  qmake project file for the test-interp example (part of GEMS3K standalone code)
# © 2016 GEMS Developer Team
 
TEMPLATE = app
LANGUAGE = C++
TARGET = test-interp
VERSION = 3.1.0

CONFIG -= qt
CONFIG += warn_on
CONFIG += debug
#CONFIG += windows
CONFIG += console

DEFINES += IPMGEMPLUGIN
#DEFINES += NODEARRAYLEVEL
DEFINES += NOPARTICLEARRAY

!win32:DEFINES += __unix

GEMS3K_CPP = ../GEMS3K
GEMS3K_H   = $$GEMS3K_CPP

TEST_INTERP_CPP = .
TEST_INTERP_H   = $$TEST_INTERP_CPP

DEPENDPATH +=
DEPENDPATH += $$TEST_INTERP_H
DEPENDPATH += $$GEMS3K_H

INCLUDEPATH += 
INCLUDEPATH += $$TEST_INTERP_H
INCLUDEPATH += $$GEMS3K_H

OBJECTS_DIR = obj

include($$TEST_INTERP_CPP/test-interp.pri)
include($$GEMS3K_CPP/gems3k.pri) 
