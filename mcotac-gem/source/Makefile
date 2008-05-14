CC = mpicc
CXX = mpiCC
FORT = mpif90
CFLAGS = -c  -g -O -DIPMGEMPLUGIN -pedantic -D__unix -D__GNU  -D__MPI -fbounds-check 
# -lstd -lC -lpgc -lrt -Msecond_underscore
CCFLAGS =   -c -g -O -DIPMGEMPLUGIN -D__unix -D__MPI -D__GNU -fbounds-check 
CCFLAGS1 =   -c -g -O -DIPMGEMPLUGIN -D__unix -D__MPI -D__GNU -fbounds-check 
FFLAGS = -c -g -O -fno-underscoring -fbounds-check -DIPMGEMPLUGIN -D__unix -D__GNU -D__MPI 
LINK = mpif90
LFLAGS = -fbounds-check -lstdc++ -lgfortran
GEMINCLUDE = ../../GEMIPM2K
GEMLIB     = ../../GEMIPM2K
EXEC = mcotacgem-mpi

${EXEC}: mcotac1d.o datio.o erech.o inistat2.o inparf.o setpar_.o walk2_.o wegdat1d_.o hydro1d_.o \
       f_gem_node.o
       
	${LINK} -o ${EXEC}  ${LFLAGS} mcotac1d.o datio.o \
                 erech.o inistat2.o inparf.o setpar_.o walk2_.o wegdat1d_.o hydro1d_.o f_gem_node.o -L${GEMLIB} -lgemipm2k

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

f_gem_node.o: f_gem_node.cpp ${GEMINCLUDE}/node.h \
		${GEMINCLUDE}/m_param.h \
		${GEMINCLUDE}/datach.h \
		${GEMINCLUDE}/databr.h \
		${GEMINCLUDE}/gdatastream.h \
		${GEMINCLUDE}/ms_multi.h \
		${GEMINCLUDE}/verror.h \
		${GEMINCLUDE}/gstring.h \
		${GEMINCLUDE}/array.h \
		${GEMINCLUDE}/m_const.h \
		${GEMINCLUDE}/v_user.h
	$(CXX)  $(CCFLAGS1) -I${GEMINCLUDE}  f_gem_node.cpp

#gwheader.h: gwheader.h
#f_gem_node.h: f_gem_node.h
#gwheader.inc: gwheader.inc
#kinetics.inc: kinetics.inc

clean: 
	rm -f *.o *.ti *.fpp