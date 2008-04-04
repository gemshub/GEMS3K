CC = gcc
CXX = g++
FORT = gfortran-4.2
CFLAGS = -c  -g -O2 -pedantic -fno-nonansi-builtins 
# -lstd -lC -lpgc -lrt -Msecond_underscore
CCFLAGS =   -c -g -O2 -DIPMGEMPLUGIN -D__unix
CCFLAGS1 =   -c -g -O2 -DIPMGEMPLUGIN -D__unix
FFLAGS = -c -g -O2 -fbounds-check -fno-underscoring
LINK = gfortran-4.2
LFLAGS = -lstdc++
EXEC = mcotacgem

${EXEC}: mcotac1d.o datio.o erech.o inistat2.o inparf.o setpar_.o walk2_.o wegdat1d_.o hydro1d_.o \
       gstring.o gdatastream.o s_fgl.o s_fgl2.o node.o node_format.o ms_multi_file.o \
       ms_multi_format.o ipm_chemical.o ipm_chemical2.o ipm_chemical3.o \
       ipm_main.o ipm_simplex.o io_arrays.o f_gem_node.o ms_param.o num_methods.o
       
	${LINK} -o ${EXEC}  ${LFLAGS} mcotac1d.o datio.o \
                 erech.o inistat2.o inparf.o setpar_.o walk2_.o wegdat1d_.o hydro1d_.o \
       gstring.o gdatastream.o s_fgl.o s_fgl2.o node.o node_format.o ms_multi_file.o \
       ms_multi_format.o ipm_chemical.o ipm_chemical2.o ipm_chemical3.o \
       ipm_main.o ipm_simplex.o io_arrays.o f_gem_node.o ms_param.o num_methods.o

mcotac1d.o: mcotac1d.f gwheader.inc kinetics.inc
	cp mcotac1d.f mcotac1d.fpp
	$(FORT)  $(FFLAGS) mcotac1d.fpp
datio.o: datio.f gwheader.inc kinetics.inc
	cp datio.f datio.fpp
	$(FORT)  $(FFLAGS) datio.fpp
erech.o: erech.f gwheader.inc kinetics.inc
	cp erech.f erech.fpp
	$(FORT)  $(FFLAGS) erech.fpp
inistat2.o: inistat2.f gwheader.inc kinetics.inc
	cp inistat2.f inistat2.fpp
	$(FORT)  $(FFLAGS)  inistat2.fpp
inparf.o: inparf.f gwheader.inc
	cp inparf.f inparf.fpp
	$(FORT)  $(FFLAGS)  inparf.fpp
setpar_.o: setpar_.c gwheader.h
	$(CC)  $(CFLAGS) setpar_.c
walk2_.o: walk2_.c gwheader.h
	$(CC)  $(CFLAGS) walk2_.c
wegdat1d_.o: wegdat1d_.c gwheader.h
	$(CC)  $(CFLAGS) wegdat1d_.c
hydro1d_.o: hydro1d_.c  gwheader.h
	$(CC)  $(CFLAGS)  hydro1d_.c

gstring.o: gstring.cpp gstring.h \
		array.h
	$(CXX)  $(CCFLAGS) gstring.cpp
gdatastream.o: gdatastream.cpp gdatastream.h \
		verror.h \
		gstring.h \
		array.h
	$(CXX)  $(CCFLAGS) gdatastream.cpp
s_fgl.o: s_fgl.cpp s_fgl.h \
		m_const.h \
		verror.h \
		v_user.h \
		gstring.h \
		array.h
	$(CXX)  $(CCFLAGS) s_fgl.cpp
s_fgl2.o: s_fgl2.cpp s_fgl.h \
		m_const.h \
		verror.h \
		v_user.h \
		gstring.h \
		array.h
	$(CXX)  $(CCFLAGS) s_fgl2.cpp
node.o: node.cpp node.h \
		gdatastream.h \
		m_param.h \
		datach.h \
		databr.h \
		ms_multi.h \
		verror.h \
		gstring.h \
		array.h \
		m_const.h \
		v_user.h
	$(CXX)  $(CCFLAGS) node.cpp
node_format.o: node_format.cpp io_arrays.h \
                node.h \
		gdatastream.h \
		m_param.h \
		datach.h \
		databr.h \
		ms_multi.h \
		verror.h \
		gstring.h \
		array.h \
		m_const.h \
		v_user.h
	$(CXX)  $(CCFLAGS) node_format.cpp
ms_multi_file.o: ms_multi_file.cpp  io_arrays.h \
                m_param.h \
		gdatastream.h \
		ms_multi.h \
		verror.h \
		gstring.h \
		array.h \
		m_const.h \
		v_user.h
	$(CXX)  $(CCFLAGS) ms_multi_file.cpp
ms_multi_format.o: ms_multi_format.cpp  io_arrays.h \
                m_param.h \
		node.h \
		gdatastream.h \
		ms_multi.h \
		verror.h \
		gstring.h \
		array.h \
		m_const.h \
		v_user.h \
		datach.h \
		databr.h
	$(CXX)  $(CCFLAGS) ms_multi_format.cpp
ipm_chemical.o: ipm_chemical.cpp m_param.h \
		gdatastream.h \
		ms_multi.h \
		verror.h \
		gstring.h \
		array.h \
		m_const.h \
		v_user.h
	$(CXX)  $(CCFLAGS) ipm_chemical.cpp
ipm_chemical2.o: ipm_chemical2.cpp m_param.h \
		gdatastream.h \
		ms_multi.h \
		verror.h \
		gstring.h \
		array.h \
		m_const.h \
		v_user.h
	$(CXX)  $(CCFLAGS) ipm_chemical2.cpp
ipm_chemical3.o: ipm_chemical3.cpp m_param.h \
		s_fgl.h \
		gdatastream.h \
		ms_multi.h \
		verror.h \
		gstring.h \
		array.h \
		m_const.h \
		v_user.h
	$(CXX)  $(CCFLAGS) ipm_chemical3.cpp
ipm_main.o: ipm_main.cpp m_param.h \
		jama_lu.h \
		jama_cholesky.h \
		node.h \
		gdatastream.h \
		ms_multi.h \
		verror.h \
		gstring.h \
		array.h \
		m_const.h \
		v_user.h \
		tnt.h \
		tnt_array1d.h \
		tnt_array2d.h \
		tnt_i_refvec.h \
		datach.h \
		databr.h
	$(CXX)  $(CCFLAGS) -DPGf90 ipm_main.cpp
ipm_simplex.o: ipm_simplex.cpp m_param.h \
		gdatastream.h \
		ms_multi.h \
		verror.h \
		gstring.h \
		array.h \
		m_const.h \
		v_user.h
	$(CXX)  $(CCFLAGS) ipm_simplex.cpp
io_arrays.o: io_arrays.cpp io_arrays.h \
		verror.h \
		gstring.h \
		array.h
	$(CXX)  $(CCFLAGS1) io_arrays.cpp
f_gem_node.o: f_gem_node.cpp node.h \
		m_param.h \
		datach.h \
		databr.h \
		gdatastream.h \
		ms_multi.h \
		verror.h \
		gstring.h \
		array.h \
		m_const.h \
		v_user.h
	$(CXX)  $(CCFLAGS1) f_gem_node.cpp
ms_param.o: ms_param.cpp m_param.h \
		num_methods.h \
		gdatastream.h \
		node.h \
		ms_multi.h \
		verror.h \
		gstring.h \
		array.h \
		m_const.h \
		v_user.h \
		datach.h \
		databr.h
	$(CXX)  $(CCFLAGS) ms_param.cpp
num_methods.o: num_methods.cpp v_user.h \
		num_methods.h \
		verror.h \
		gstring.h \
		array.h
	$(CXX)  $(CCFLAGS) num_methods.cpp

#gwheader.h: gwheader.h
#f_gem_node.h: f_gem_node.h
#gwheader.inc: gwheader.inc
#kinetics.inc: kinetics.inc

clean: 
	rm -f *.o *.ti