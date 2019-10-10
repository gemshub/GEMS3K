#  qmake project file for the node-gem example (part of GEMS3K standalone code)
#  2012-2019 GEMS Developer Team
 
TEMPLATE = app
LANGUAGE = C++
TARGET = node-gem
VERSION = 3.4.6

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

NODE_GEM_CPP = .
NODE_GEM_H   = $$NODE_GEM_CPP

DEPENDPATH +=
DEPENDPATH += $$NODE_GEM_H
DEPENDPATH += $$GEMS3K_H

INCLUDEPATH += 
INCLUDEPATH += $$NODE_GEM_H
INCLUDEPATH += $$GEMS3K_H

OBJECTS_DIR = obj

include($$NODE_GEM_CPP/node-gem.pri)
include($$GEMS3K_CPP/gems3k.pri) 
