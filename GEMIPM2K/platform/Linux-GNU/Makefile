#############################################################################
# GEMIPM2K - standalone GEMS version 
# Makefile for building: gemcalc
#############################################################################

####### Compiler, tools and options

CC       = g++
CXX      = g++
LEX      = flex
YACC     = yacc
CFLAGS   = -pipe -Wall -W -g  -DIPMGEMPLUGIN -D__unix
CXXFLAGS = -pipe -Wall -W -g  -DIPMGEMPLUGIN -D__unix
LEXFLAGS = 
YACCFLAGS= -d
INCPATH  = -I.
LINK     = gcc
LFLAGS   = -lstdc++
LIBS     = $(SUBLIBS)  
AR       = ar cqs
RANLIB   = 
TAR      = tar -cf
GZIP     = gzip -9f
COPY     = cp -f
COPY_FILE= $(COPY)
COPY_DIR = $(COPY) -r
INSTALL_FILE= $(COPY_FILE)
INSTALL_DIR = $(COPY_DIR)
DEL_FILE = rm -f
SYMLINK  = ln -sf
DEL_DIR  = rmdir
MOVE     = mv -f
CHK_DIR_EXISTS= test -d
MKDIR    = mkdir -p

####### Output directory

OBJECTS_DIR = obj/

####### Files

HEADERS = array.h \
		verror.h \
		gstring.h \
		gdatastream.h \
		v_user.h \
		s_fgl.h \
		databr.h \
		datach.h \
		node.h \
		ms_multi.h \
		io_arrays.h \
		m_const.h \
		m_param.h \
		num_methods.h \
		tnt_i_refvec.h \
		tnt_array1d.h \
		tnt_array2d.h \
		tnt.h \
		jama_cholesky.h \
		jama_lu.h
SOURCES = gstring.cpp \
		gdatastream.cpp \
		s_fgl.cpp \
		s_fgl1.cpp \
		s_fgl2.cpp \
		s_fgl3.cpp \
		s_fgl4.cpp \
		node.cpp \
		node_format.cpp \
		ms_multi_file.cpp \
		ms_multi_format.cpp \
		ipm_chemical.cpp \
		ipm_chemical2.cpp \
		ipm_chemical3.cpp \
		ipm_main.cpp \
		ipm_simplex.cpp \
		io_arrays.cpp \
		main.cpp \
		ms_param.cpp \
		num_methods.cpp
OBJECTS = obj/gstring.o \
		obj/gdatastream.o \
		obj/s_fgl.o \
		obj/s_fgl1.o \
		obj/s_fgl2.o \
		obj/s_fgl3.o \
		obj/s_fgl4.o \
		obj/node.o \
		obj/node_format.o \
		obj/ms_multi_file.o \
		obj/ms_multi_format.o \
		obj/ipm_chemical.o \
		obj/ipm_chemical2.o \
		obj/ipm_chemical3.o \
		obj/ipm_main.o \
		obj/ipm_simplex.o \
		obj/io_arrays.o \
		obj/main.o \
		obj/ms_param.o \
		obj/num_methods.o


DESTDIR  = 
TARGET   = gemcalc

first: all
####### Implicit rules

.SUFFIXES: .c .o .cpp .cc .cxx .C

.cpp.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o $@ $<

.cc.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o $@ $<

.cxx.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o $@ $<

.C.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o $@ $<

.c.o:
	$(CC) -c $(CFLAGS) $(INCPATH) -o $@ $<

####### Build rules

all: Makefile $(TARGET)

$(TARGET):  $(OBJECTS)  
	$(LINK) $(LFLAGS) -o $(TARGET) $(OBJECTS)  $(OBJCOMP) $(LIBS)


clean:
	-$(DEL_FILE) $(OBJECTS)
	-$(DEL_FILE) *~ core *.core


####### Sub-libraries

distclean: clean
	-$(DEL_FILE) $(TARGET)


FORCE:

####### Compile

obj/gstring.o: gstring.cpp gstring.h \
		v_user.h \
		array.h \
		verror.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/gstring.o gstring.cpp

obj/gdatastream.o: gdatastream.cpp gdatastream.h \
		verror.h \
		gstring.h \
		array.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/gdatastream.o gdatastream.cpp

obj/s_fgl.o: s_fgl.cpp s_fgl.h \
		verror.h \
		m_const.h \
		gstring.h \
		array.h \
		v_user.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/s_fgl.o s_fgl.cpp

obj/s_fgl1.o: s_fgl1.cpp verror.h \
		s_fgl.h \
		gstring.h \
		array.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/s_fgl1.o s_fgl1.cpp

obj/s_fgl2.o: s_fgl2.cpp s_fgl.h \
		m_const.h \
		verror.h \
		v_user.h \
		gstring.h \
		array.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/s_fgl2.o s_fgl2.cpp

obj/s_fgl3.o: s_fgl3.cpp verror.h \
		s_fgl.h \
		gstring.h \
		array.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/s_fgl3.o s_fgl3.cpp

obj/s_fgl4.o: s_fgl4.cpp verror.h \
		s_fgl.h \
		gstring.h \
		array.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/s_fgl4.o s_fgl4.cpp

obj/node.o: node.cpp node.h \
		gdatastream.h \
		num_methods.h \
		m_param.h \
		datach.h \
		databr.h \
		ms_multi.h \
		verror.h \
		gstring.h \
		array.h \
		m_const.h \
		s_fgl.h \
		v_user.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/node.o node.cpp

obj/node_format.o: node_format.cpp io_arrays.h \
		node.h \
		gdatastream.h \
		gstring.h \
		array.h \
		m_param.h \
		datach.h \
		databr.h \
		ms_multi.h \
		verror.h \
		m_const.h \
		s_fgl.h \
		v_user.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/node_format.o node_format.cpp

obj/ms_multi_file.o: ms_multi_file.cpp io_arrays.h \
		m_param.h \
		gdatastream.h \
		gstring.h \
		array.h \
		ms_multi.h \
		verror.h \
		m_const.h \
		s_fgl.h \
		v_user.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/ms_multi_file.o ms_multi_file.cpp

obj/ms_multi_format.o: ms_multi_format.cpp io_arrays.h \
		m_param.h \
		node.h \
		gstring.h \
		array.h \
		gdatastream.h \
		ms_multi.h \
		verror.h \
		m_const.h \
		s_fgl.h \
		v_user.h \
		datach.h \
		databr.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/ms_multi_format.o ms_multi_format.cpp

obj/ipm_chemical.o: ipm_chemical.cpp m_param.h \
		gdatastream.h \
		ms_multi.h \
		verror.h \
		gstring.h \
		array.h \
		m_const.h \
		s_fgl.h \
		v_user.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/ipm_chemical.o ipm_chemical.cpp

obj/ipm_chemical2.o: ipm_chemical2.cpp m_param.h \
		gdatastream.h \
		ms_multi.h \
		verror.h \
		gstring.h \
		array.h \
		m_const.h \
		s_fgl.h \
		v_user.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/ipm_chemical2.o ipm_chemical2.cpp

obj/ipm_chemical3.o: ipm_chemical3.cpp m_param.h \
		s_fgl.h \
		gdatastream.h \
		ms_multi.h \
		verror.h \
		gstring.h \
		array.h \
		m_const.h \
		v_user.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/ipm_chemical3.o ipm_chemical3.cpp

obj/ipm_main.o: ipm_main.cpp m_param.h \
		jama_lu.h \
		jama_cholesky.h \
		node.h \
		gdatastream.h \
		ms_multi.h \
		verror.h \
		gstring.h \
		array.h \
		m_const.h \
		s_fgl.h \
		v_user.h \
		tnt.h \
		tnt_array1d.h \
		tnt_array2d.h \
		tnt_i_refvec.h \
		datach.h \
		databr.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/ipm_main.o ipm_main.cpp

obj/ipm_simplex.o: ipm_simplex.cpp m_param.h \
		gdatastream.h \
		ms_multi.h \
		verror.h \
		gstring.h \
		array.h \
		m_const.h \
		s_fgl.h \
		v_user.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/ipm_simplex.o ipm_simplex.cpp

obj/io_arrays.o: io_arrays.cpp io_arrays.h \
		verror.h \
		v_user.h \
		gstring.h \
		array.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/io_arrays.o io_arrays.cpp

obj/main.o: main.cpp node.h \
		m_param.h \
		datach.h \
		databr.h \
		gdatastream.h \
		ms_multi.h \
		verror.h \
		gstring.h \
		array.h \
		m_const.h \
		s_fgl.h \
		v_user.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/main.o main.cpp

obj/ms_param.o: ms_param.cpp m_param.h \
		num_methods.h \
		gdatastream.h \
		node.h \
		ms_multi.h \
		verror.h \
		gstring.h \
		array.h \
		m_const.h \
		s_fgl.h \
		v_user.h \
		datach.h \
		databr.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/ms_param.o ms_param.cpp

obj/num_methods.o: num_methods.cpp v_user.h \
		num_methods.h \
		verror.h \
		gstring.h \
		array.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/num_methods.o num_methods.cpp

####### Install

install:  

uninstall:  

