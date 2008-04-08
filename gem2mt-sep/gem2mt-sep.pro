#qmake -spec win32-borland gemsipm2k.pro > a
#qmake -spec win32-msvc.net gemsipm2k.pro > a

TEMPLATE	= app
LANGUAGE        = C++
TARGET		= gem2mts
VERSION         = 2.2.3

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

LIBS_CPP       =  ../GEMIPM2K
KERNEL_CPP     =  ../GEMIPM2K
MODS_CPP       =  ../GEMIPM2K
SUBMOD_CPP     =  ../GEMIPM2K
NUMERICS_CPP     =  ../GEMIPM2K
GEM2MTS_CPP	= .

#LIBS_CPP       =  .
#KERNEL_CPP     =  .
#MODS_CPP       =  .
#SUBMOD_CPP     =  .
#NUMERICS_CPP     =  .

KERNEL_H     =  $$KERNEL_CPP
LIBS_H       =  $$LIBS_CPP
MODS_H       =  $$MODS_CPP
SUBMOD_H     =  $$SUBMOD_CPP
NUMERICS_H     =  $$NUMERICS_CPP
GEM2MTS_H    = $$GEM2MTS_CPP

DEPENDPATH   += ;.;$$KERNEL_H;$$LIBS_H;$$MODS_H;$$SUBMOD_H;$$NUMERICS_H;GEM2MTS_H
INCLUDEPATH   +=;.;$$KERNEL_H;$$LIBS_H;$$MODS_H;$$SUBMOD_H;$$NUMERICS_H;GEM2MTS_H

OBJECTS_DIR       = obj

       HEADERS	 += $$LIBS_H/array.h \
                    $$LIBS_H/verror.h  \
                    $$LIBS_H/gstring.h  \
                    $$LIBS_H/gdatastream.h  \
                    $$KERNEL_H/v_user.h \
                    $$MODS_H/s_fgl.h \
                    $$SUBMOD_H/databr.h \
                    $$SUBMOD_H/datach.h \
                    $$SUBMOD_H/nodearray.h \
                    $$SUBMOD_H/node.h \
                    $$GEM2MTS_H/particlearray.h \
                    $$SUBMOD_H/ms_multi.h \
	            $$SUBMOD_H/io_arrays.h \
		    $$GEM2MTS_H/ms_gem2mt.h \
		    $$MODS_H/m_const.h \
		    $$MODS_H/m_param.h  \
                     $$NUMERICS_H/tnt_i_refvec.h \
                     $$NUMERICS_H/tnt_array1d.h \
                     $$NUMERICS_H/tnt_array2d.h \
                     $$NUMERICS_H/tnt.h \
                     $$NUMERICS_H/jama_cholesky.h \
                     $$NUMERICS_H/jama_lu.h \
                     $$NUMERICS_H/num_methods.h

        SOURCES	  +=  $$LIBS_CPP/gstring.cpp  \
                      $$LIBS_CPP/gdatastream.cpp  \
                    $$NUMERICS_CPP/num_methods.cpp \
                    $$MODS_CPP/s_fgl.cpp \
                    $$MODS_CPP/s_fgl2.cpp \
                    $$GEM2MTS_CPP/m_gem2mtt.cpp \
		    $$GEM2MTS_CPP/m_gem2mtbox.cpp \
		    $$GEM2MTS_CPP/ms_gem2mt.cpp \
                    $$SUBMOD_CPP/nodearray.cpp \
                    $$SUBMOD_CPP/node.cpp \
                    $$SUBMOD_CPP/node_format.cpp \
                    $$GEM2MTS_CPP/particlearray.cpp \
                    $$SUBMOD_CPP/ms_multi_file.cpp \
                    $$SUBMOD_CPP/ms_multi_format.cpp \
                    $$SUBMOD_CPP/ipm_chemical.cpp \
                    $$SUBMOD_CPP/ipm_chemical2.cpp \
                    $$SUBMOD_CPP/ipm_chemical3.cpp \
                    $$SUBMOD_CPP/ipm_main.cpp \
                    $$SUBMOD_CPP/ipm_simplex.cpp \
		    $$SUBMOD_CPP/io_arrays.cpp \
                    $$GEM2MTS_CPP/main.cpp \
		    $$MODS_CPP/ms_param.cpp
