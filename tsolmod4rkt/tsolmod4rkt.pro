#  qmake project file for the kva2json example (part of GEMS3K standalone code)
# (c) 2020 GEMS Developer Team
 
TEMPLATE = app
LANGUAGE = C++
TARGET = kva2json
VERSION = 3.4.6

CONFIG -= qt
CONFIG += warn_on
CONFIG += thread console
CONFIG += c++17
CONFIG += sanitaze sanitaze_thread

DEFINES += NODEARRAYLEVEL
#DEFINES += USE_NLOHMANNJSON
DEFINES += USE_THERMOFUN
DEFINES += USE_THERMO_LOG
DEFINES += OVERFLOW_EXCEPT  #compile with nan inf exceptions
DEFINES += SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_OFF

!win32 {

DEFINES += __unix
QMAKE_CFLAGS += -pedantic -Wall -Wextra -Wwrite-strings -Werror

QMAKE_CXXFLAGS += -fPIC -Wall -Wextra -Wformat-nonliteral -Wcast-align -Wpointer-arith \
 -Wmissing-declarations -Winline \ # -Wundef \ #-Weffc++ \
 -Wcast-qual -Wshadow -Wwrite-strings -Wno-unused-parameter \
 -Wfloat-equal -pedantic -ansi

}

GEMS3K_CPP = ../GEMS3K
GEMS3K_H   = $$GEMS3K_CPP

DEPENDPATH += .
DEPENDPATH += $$GEMS3K_H

INCLUDEPATH += .
INCLUDEPATH += $$GEMS3K_H

contains(DEFINES, USE_THERMOFUN) {
#ThermoFun_CPP   =  ../ThermoFun
#ThermoFun_H     =   $$ThermoFun_CPP
#DEPENDPATH += $$ThermoFun_H
#INCLUDEPATH += $$ThermoFun_H
#include($$ThermoFun_CPP/ThermoFun.pri)
LIBS += -lThermoFun -lChemicalFun

} ## end USE_THERMOFUN


OBJECTS_DIR = obj

!contains(DEFINES, USE_NLOHMANNJSON) {
  SJSON_CPP = $$GEMS3K_H/..
  SJSON_H   = $$SJSON_CPP
  DEPENDPATH += $$SJSON_H
  INCLUDEPATH += $$SJSON_H
  HEADERS	 += $$SJSON_H/simdjson/simdjson.h
}

HEADERS	 +=         $$GEMS3K_H/verror.h  \
                    $$GEMS3K_H/gdatastream.h  \
                    $$GEMS3K_H/num_methods.h \
                    $$GEMS3K_H/s_solmod.h \
                    $$GEMS3K_H/m_const_base.h  \
                    $$GEMS3K_H/databr.h \
                    $$GEMS3K_H/datach.h \
                    $$GEMS3K_H/io_template.h \
                    $$GEMS3K_H/io_nlohmann.h \
                    $$GEMS3K_H/io_keyvalue.h \
                    $$GEMS3K_H/io_simdjson.h \
                    $$GEMS3K_H/gems3k_impex.h \
                    $$GEMS3K_H/v_detail.h \
                    $$GEMS3K_H/v_service.h \
                    $$GEMS3K_H/jsonconfig.h \
                    $$GEMS3K_H/datach_api.h \
                    solmodcalc.h \
                    tsolmod_multi.h

SOURCES	  +=          $$GEMS3K_CPP/gdatastream.cpp  \
                      $$GEMS3K_CPP/num_methods.cpp \
                      $$GEMS3K_CPP/s_solmod.cpp \
                      $$GEMS3K_CPP/s_solmod2.cpp \
                      $$GEMS3K_CPP/s_solmod3.cpp \
                      $$GEMS3K_CPP/s_solmod4.cpp \
                      $$GEMS3K_CPP/s_solmod5.cpp \
                      $$GEMS3K_CPP/s_solmod6.cpp \
                      $$GEMS3K_CPP/s_sorpmod.cpp \
                      $$GEMS3K_CPP/io_template.cpp \
                      $$GEMS3K_CPP/io_nlohmann.cpp \
                      $$GEMS3K_CPP/io_keyvalue.cpp \
                      $$GEMS3K_CPP/io_simdjson.cpp \
                      $$GEMS3K_CPP/gems3k_impex.cpp \
                      $$GEMS3K_CPP/v_detail.cpp \
                      $$GEMS3K_CPP/v_service.cpp \
                      $$GEMS3K_CPP/jsonconfig.cpp \
                      $$GEMS3K_CPP/datach_api.cpp \
                      $$GEMS3K_CPP/datach_formats.cpp \
                      solmodcalc.cpp \
                      tsolmod4rkt.cpp \
                      tsolmod_multi_add.cpp \
                      tsolmod_multi_alloc.cpp \
                      tsolmod_multi_file.cpp \
                      tsolmod_multi_format.cpp

