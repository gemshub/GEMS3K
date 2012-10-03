#  qmake project file for the node-gem example (part of GEMS3K standalone code)
# © 2012 GEMS Developer Team
 
TEMPLATE = app
LANGUAGE = C++
TARGET = dm-node-gem
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

DM_NODE_GEM_CPP = .
DM_NODE_GEM_H   = $$DM_NODE_GEM_CPP

DEPENDPATH +=
DEPENDPATH += $$DM_NODE_GEM_H
DEPENDPATH += $$GEMS3K_H

INCLUDEPATH += 
INCLUDEPATH += $$DM_NODE_GEM_H
INCLUDEPATH += $$GEMS3K_H

OBJECTS_DIR = obj

include($$DM_NODE_GEM_CPP/dm-node-gem.pri)
include($$GEMS3K_CPP/gems3k.pri) 
