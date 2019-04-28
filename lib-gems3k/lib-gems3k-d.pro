##########################################################################################
#      qmake project file for generation of a dynamically linked GEMS3K library          #
#  configure release version:	            qmake "CONFIG += release" LIB_GEMS3K.pro     #
#  configure debug version:	            qmake "CONFIG += debug" LIB_GEMS3K.pro       #
#   for optional debug output of the ELVIS model add DEFINES += ELVIS_DEBUG              #
##########################################################################################

TEMPLATE	= lib
LANGUAGE	= C++
TARGET		= gems3k
VERSION		= 3.4.6

CONFIG		-= qt
CONFIG		+= warn_on
CONFIG		+= console

QMAKE_CC	= gcc
QMAKE_CXX	= g++

#CONFIG( release,  debug|release ) {
#	message( "Configuring for release build ..." )
#	QMAKE_CFLAGS_RELEASE = -O2
#	QMAKE_CXXFLAGS_RELEASE = -O2
#}

#CONFIG( debug,  debug|release ) {
#	message( "Configuring for debug build ..." )
#	QMAKE_CFLAGS_DEBUG   = -g -Wall -pedantic -fexceptions
#	QMAKE_CXXFLAGS_DEBUG = -g -Wall -pedantic -fexceptions
#}

DEFINES += IPMGEMPLUGIN
#DEFINES += NODEARRAYLEVEL
DEFINES += NOPARTICLEARRAY

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

QMAKE_LFLAGS += 
OBJECTS_DIR = obj

include($$GEMS3K_CPP/gems3k.pri) 
