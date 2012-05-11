#qmake -spec win32-borland nodearray-gem.pro > a
#qmake -spec win32-msvc.net nodearray-gem.pro > a

TEMPLATE	= app
LANGUAGE        = C++
TARGET		= nodearrs
VERSION         = 2.2.3

CONFIG		-= qt
CONFIG		+=  warn_on debug windows
CONFIG		+= console

DEFINES         += IPMGEMPLUGIN
DEFINES         += NOPARTICLEARRAY

!win32 {
  DEFINES += __unix
}


win32-borland {
    	DEFINES += __win32_borland
        #  Debug, RTTI, exceptions, Visual C - compatible
        QMAKE_CFLAGS += -x -xd -xp -VM -RT
        QMAKE_CXXFLAGS += -x -xd -xp -VM -RT
}

LIBS_CPP       =  ../GEMS3K
KERNEL_CPP     =  ../GEMS3K
MODS_CPP       =  ../GEMS3K
SUBMOD_CPP     =  ../GEMS3K
NUMERICS_CPP   =  ../GEMS3K
GEMIPM_CPP     =  ../GEMS3K

#LIBS_CPP       =  .
#KERNEL_CPP     =  .
#MODS_CPP       =  .
#SUBMOD_CPP     =  .
#NUMERICS_CPP   =  .
#GEMIPM_CPP     =  .

KERNEL_H     =  $$KERNEL_CPP
LIBS_H       =  $$LIBS_CPP
MODS_H       =  $$MODS_CPP
SUBMOD_H     =  $$SUBMOD_CPP
NUMERICS_H   =  $$NUMERICS_CPP
GEMIPM_H     =  $$GEMIPM_CPP


DEPENDPATH   += ;.;$$GEMIPM_H;$$KERNEL_H;$$LIBS_H;$$MODS_H;$$SUBMOD_H;$$NUMERICS_H
INCLUDEPATH   +=;.;$$GEMIPM_H;$$KERNEL_H;$$LIBS_H;$$MODS_H;$$SUBMOD_H;$$NUMERICS_H

OBJECTS_DIR       = obj

       HEADERS	 += $$LIBS_H/verror.h  \
                    $$LIBS_H/gdatastream.h  \
                    $$KERNEL_H/v_user.h \
                    $$MODS_H/s_fgl.h \
                    $$SUBMOD_H/databr.h \
                    $$SUBMOD_H/datach.h \
                    $$SUBMOD_H/nodearray.h \
                    $$SUBMOD_H/node.h \
                    $$SUBMOD_H/ms_multi.h \
	                $$SUBMOD_H/io_arrays.h \
		            $$GEMIPM_H/m_const.h \
		            $$GEMIPM_H/m_param.h  \
                    $$NUMERICS_H/num_methods.h \
                    $$NUMERICS_H/tnt_i_refvec.h \
                    $$NUMERICS_H/tnt_array1d.h \
                    $$NUMERICS_H/tnt_array2d.h \
                    $$NUMERICS_H/tnt.h \
                    $$NUMERICS_H/jama_cholesky.h \
                    $$NUMERICS_H/jama_lu.h

        SOURCES	  +=  $$LIBS_CPP/gdatastream.cpp  \
                     $$MODS_CPP/s_fgl.cpp \
                     $$MODS_CPP/s_fgl1.cpp \
                     $$MODS_CPP/s_fgl2.cpp \
                     $$MODS_CPP/s_fgl3.cpp \
                     $$MODS_CPP/s_fgl4.cpp \
                    $$SUBMOD_CPP/node.cpp \
                    $$SUBMOD_CPP/nodearray.cpp \
                    $$SUBMOD_CPP/node_format.cpp \
                    $$SUBMOD_CPP/ms_multi_file.cpp \
                    $$SUBMOD_CPP/ms_multi_format.cpp \
                    $$SUBMOD_CPP/ipm_chemical.cpp \
                    $$SUBMOD_CPP/ipm_chemical2.cpp \
                    $$SUBMOD_CPP/ipm_chemical3.cpp \
                    $$SUBMOD_CPP/ipm_main.cpp \
                    $$SUBMOD_CPP/ipm_simplex.cpp \
		            $$SUBMOD_CPP/io_arrays.cpp \
                    main.cpp \
		            $$GEMIPM_CPP/ms_param.cpp \
                    $$NUMERICS_CPP/num_methods.cpp
