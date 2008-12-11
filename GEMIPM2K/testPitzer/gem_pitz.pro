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
   
   HEADERS	 += ../s_fgl.h \
                    ../verror.h   \
                    ../array.h \
                    ../v_user.h  \
                    ../gstring.h 
   
   SOURCES	 +=   ../s_fgl.cpp \
                      ../s_fgl1.cpp \
                      ../s_fgl2.cpp \
                      ../s_fgl3.cpp \
                      ../s_fgl4.cpp  \
                      ../gstring.cpp \ 
                      main.cpp 
