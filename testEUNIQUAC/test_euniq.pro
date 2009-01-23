#qmake -spec  win32-borland gemsipm2k.pro > a
#qmake -spec win32-msvc.net gemsipm2k.pro > a

TEMPLATE	= app
LANGUAGE        = C++
TARGET		= euniquac
VERSION         = 2.3.0

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

DEPENDPATH += ;.;../GEMIPM2K 
INCLUDEPATH += ;.;../GEMIPM2K 
   
   HEADERS	 += 	../GEMIPM2K/s_fgl.h \
                    ../GEMIPM2K/verror.h   \
                    ../GEMIPM2K/array.h \
                    ../GEMIPM2K/v_user.h  \
                    ../GEMIPM2K/gstring.h 
   
   SOURCES	 +=     ../GEMIPM2K/s_fgl.cpp \
          			../GEMIPM2K/s_fgl1.cpp  \
          			../GEMIPM2K/s_fgl2.cpp \
          			../GEMIPM2K/s_fgl3.cpp  \
          			../GEMIPM2K/s_fgl4.cpp  \
                    ../GEMIPM2K/gstring.cpp \ 
                    main.cpp 
