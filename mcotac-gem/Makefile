CC = pgcc
CXX = pgCC
FORT = pgf90
CFLAGS = -c  -g -O2 -DIPMGEMPLUGIN -D__unix -D__PGI -Mpreprocess
# -lstd -lC -lpgc -lrt -Msecond_underscore
CCFLAGS =   -c -g -O2 -DIPMGEMPLUGIN -D__unix -D__PGI -Mpreprocess
CCFLAGS1 =   -c -g -O2 -DIPMGEMPLUGIN -D__unix -D__PGI -Mpreprocess
FFLAGS = -c -g -O2 -DIPMGEMPLUGIN -D__unix -D__GNU -D__PGI -Mnosecond_underscore -Mpreprocess
LINK = pgf90
LFLAGS = 
LIBS = -lrt -lstd -lC -lpgc 
GEMINCLUDE = ../GEMIPM2K
GEMLIB     = -L../GEMIPM2K -lgemipm2k-pgi
EXEC = mcotacgem-pgi

${EXEC}: mcotac1d.o datio.o erech.o inistat2.o inparf.o setpar_.o walk2_.o wegdat1d_.o hydro1d_.o \
       f_gem_node.o
       
	${LINK} -o ${EXEC}  ${LFLAGS} mcotac1d.o datio.o \
          erech.o inistat2.o inparf.o setpar_.o walk2_.o wegdat1d_.o hydro1d_.o f_gem_node.o  ${GEMLIB}  ${LIBS}

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


f_gem_node.o: f_gem_node.cpp ../GEMIPM2K/node.h \
		../GEMIPM2K/m_param.h \
		../GEMIPM2K/datach.h \
		../GEMIPM2K/databr.h \
		../GEMIPM2K/gdatastream.h \
		../GEMIPM2K/ms_multi.h \
		../GEMIPM2K/verror.h \
		../GEMIPM2K/gstring.h \
		../GEMIPM2K/array.h \
		../GEMIPM2K/m_const.h \
		../GEMIPM2K/v_user.h
	$(CXX)  $(CCFLAGS1) -I${GEMINCLUDE}  f_gem_node.cpp

#gwheader.h: gwheader.h
#f_gem_node.h: f_gem_node.h
#gwheader.inc: gwheader.inc
#kinetics.inc: kinetics.inc

clean: 
	rm -f *.o *.ti
