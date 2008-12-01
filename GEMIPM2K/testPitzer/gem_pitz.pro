#qmake -spec  win32-borland gemsipm2k.pro > a
#qmake -spec win32-msvc.net gemsipm2k.pro > a

TEMPLATE	= app
LANGUAGE        = C++
TARGET		= gempitz
VERSION         = 2.2.0

CONFIG		-= qt
CONFIG		+=  warn_on debug windows
CONFIG		+= console

DEFINES         += IPMGEMPLUGIN


!win32 {
  DEFINES += __unix
}

win32-borland {
       	DEFINES += __win32_borland
        #  Debug, RTTI, exceptions, Visual C - compatible
        QMAKE_CFLAGS += -x -xd -xp -VM -RT
        QMAKE_CXXFLAGS += -x -xd -xp -VM -RT
}

DEPENDPATH += ;.;../ 
INCLUDEPATH += ;.;../ 
   
   HEADERS	 += s_pitzer.h \
                    ../verror.h   \
                    ../array.h \
                    ../v_user.h  \
                    ../gstring.h 
   
   SOURCES	 +=  s_pitzer.cpp  \
                     ../gstring.cpp \ 
                 main.cpp 
