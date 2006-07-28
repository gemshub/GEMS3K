# Microsoft Developer Studio Generated NMAKE File, Based on test.dsp
!IF "$(CFG)" == ""
CFG=test - Win32 Debug
!MESSAGE No configuration specified. Defaulting to test - Win32 Debug.
!ENDIF 

!IF "$(CFG)" != "test - Win32 Release" && "$(CFG)" != "test - Win32 Debug"
!MESSAGE Invalid configuration "$(CFG)" specified.
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "test.mak" CFG="test - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "test - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "test - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 
!ERROR An invalid configuration is specified.
!ENDIF 

!IF "$(OS)" == "Windows_NT"
NULL=
!ELSE 
NULL=nul
!ENDIF 

!IF  "$(CFG)" == "test - Win32 Release"

OUTDIR=.\Release
INTDIR=.\Release
# Begin Custom Macros
OutDir=.\Release
# End Custom Macros

ALL : "$(OUTDIR)\test.exe"


CLEAN :
	-@erase "$(INTDIR)\datio.obj"
	-@erase "$(INTDIR)\erech.obj"
	-@erase "$(INTDIR)\gdatastream.obj"
	-@erase "$(INTDIR)\gstring.obj"
	-@erase "$(INTDIR)\hydro1d_.obj"
	-@erase "$(INTDIR)\inistat2.obj"
	-@erase "$(INTDIR)\inparf.obj"
	-@erase "$(INTDIR)\ipm_chemical.obj"
	-@erase "$(INTDIR)\ipm_chemical2.obj"
	-@erase "$(INTDIR)\ipm_chemical3.obj"
	-@erase "$(INTDIR)\ipm_main.obj"
	-@erase "$(INTDIR)\ipm_simplex.obj"
	-@erase "$(INTDIR)\mainfromf.obj"
	-@erase "$(INTDIR)\mcotac1d.obj"
	-@erase "$(INTDIR)\ms_multi_file.obj"
	-@erase "$(INTDIR)\ms_multi_format.obj"
	-@erase "$(INTDIR)\ms_param.obj"
	-@erase "$(INTDIR)\node.obj"
	-@erase "$(INTDIR)\node_format.obj"
	-@erase "$(INTDIR)\s_fgl.obj"
	-@erase "$(INTDIR)\setpar.obj"
	-@erase "$(INTDIR)\setpar_.obj"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(INTDIR)\walk2_.obj"
	-@erase "$(INTDIR)\wegdat1d_.obj"
	-@erase "$(OUTDIR)\test.exe"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

F90=pgf90
F90_PROJ=/compile_only /include:"$(INTDIR)\\" /nologo /warn:nofileopt /module:"Release/" /object:"Release/" 
F90_OBJS=.\Release/

.SUFFIXES: .fpp

.for{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.f{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.f90{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.fpp{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

CPP=pgcc.exe
CPP_PROJ=/nologo /ML /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /Fp"$(INTDIR)\test.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /c 

.c{$(INTDIR)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cpp{$(INTDIR)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cxx{$(INTDIR)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.c{$(INTDIR)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cpp{$(INTDIR)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cxx{$(INTDIR)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

RSC=rc.exe
BSC32=bscmake.exe
BSC32_FLAGS=/nologo /o"$(OUTDIR)\test.bsc" 
BSC32_SBRS= \
	
LINK32=pgf90
LINK32_FLAGS=kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /incremental:no /pdb:"$(OUTDIR)\test.pdb" /machine:I386 /out:"$(OUTDIR)\test.exe" 
LINK32_OBJS= \
	"$(INTDIR)\datio.obj" \
	"$(INTDIR)\erech.obj" \
	"$(INTDIR)\gdatastream.obj" \
	"$(INTDIR)\gstring.obj" \
	"$(INTDIR)\hydro1d_.obj" \
	"$(INTDIR)\inistat2.obj" \
	"$(INTDIR)\inparf.obj" \
	"$(INTDIR)\ipm_chemical.obj" \
	"$(INTDIR)\ipm_chemical2.obj" \
	"$(INTDIR)\ipm_chemical3.obj" \
	"$(INTDIR)\ipm_main.obj" \
	"$(INTDIR)\ipm_simplex.obj" \
	"$(INTDIR)\mainfromf.obj" \
	"$(INTDIR)\mcotac1d.obj" \
	"$(INTDIR)\ms_multi_file.obj" \
	"$(INTDIR)\ms_multi_format.obj" \
	"$(INTDIR)\ms_param.obj" \
	"$(INTDIR)\node.obj" \
	"$(INTDIR)\node_format.obj" \
	"$(INTDIR)\s_fgl.obj" \
	"$(INTDIR)\setpar.obj" \
	"$(INTDIR)\setpar_.obj" \
	"$(INTDIR)\walk2_.obj" \
	"$(INTDIR)\wegdat1d_.obj"

"$(OUTDIR)\test.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ELSEIF  "$(CFG)" == "test - Win32 Debug"

OUTDIR=.\Debug
INTDIR=.\Debug
# Begin Custom Macros
OutDir=.\Debug
# End Custom Macros

ALL : "$(OUTDIR)\test.exe"


CLEAN :
	-@erase "$(INTDIR)\datio.obj"
	-@erase "$(INTDIR)\erech.obj"
	-@erase "$(INTDIR)\gdatastream.obj"
	-@erase "$(INTDIR)\gstring.obj"
	-@erase "$(INTDIR)\hydro1d_.obj"
	-@erase "$(INTDIR)\inistat2.obj"
	-@erase "$(INTDIR)\inparf.obj"
	-@erase "$(INTDIR)\ipm_chemical.obj"
	-@erase "$(INTDIR)\ipm_chemical2.obj"
	-@erase "$(INTDIR)\ipm_chemical3.obj"
	-@erase "$(INTDIR)\ipm_main.obj"
	-@erase "$(INTDIR)\ipm_simplex.obj"
	-@erase "$(INTDIR)\mainfromf.obj"
	-@erase "$(INTDIR)\mcotac1d.obj"
	-@erase "$(INTDIR)\ms_multi_file.obj"
	-@erase "$(INTDIR)\ms_multi_format.obj"
	-@erase "$(INTDIR)\ms_param.obj"
	-@erase "$(INTDIR)\node.obj"
	-@erase "$(INTDIR)\node_format.obj"
	-@erase "$(INTDIR)\s_fgl.obj"
	-@erase "$(INTDIR)\setpar.obj"
	-@erase "$(INTDIR)\setpar_.obj"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(INTDIR)\vc60.pdb"
	-@erase "$(INTDIR)\walk2_.obj"
	-@erase "$(INTDIR)\wegdat1d_.obj"
	-@erase "$(OUTDIR)\test.exe"
	-@erase "$(OUTDIR)\test.ilk"
	-@erase "$(OUTDIR)\test.pdb"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

F90=pgf90.exe
F90_PROJ=/check:bounds /compile_only /debug:full /include:"$(INTDIR)\\" /nologo /traceback /warn:argument_checking /warn:nofileopt /module:"Debug/" /object:"Debug/" /pdbfile:"Debug/DF60.PDB" 
F90_OBJS=.\Debug/

.SUFFIXES: .fpp

.for{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.f{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.f90{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.fpp{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

CPP=pgcc
CPP_PROJ=/nologo /MLd /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /D "IPMGEMPLUGIN" /Fp"$(INTDIR)\test.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /GZ /c 

.c{$(INTDIR)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cpp{$(INTDIR)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cxx{$(INTDIR)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.c{$(INTDIR)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cpp{$(INTDIR)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cxx{$(INTDIR)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

RSC=rc.exe
BSC32=bscmake.exe
BSC32_FLAGS=/nologo /o"$(OUTDIR)\test.bsc" 
BSC32_SBRS= \
	
LINK32=link.exe
LINK32_FLAGS=kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /incremental:yes /pdb:"$(OUTDIR)\test.pdb" /debug /machine:I386 /nodefaultlib:"libcd libcp" /force /out:"$(OUTDIR)\test.exe" /pdbtype:sept 
LINK32_OBJS= \
	"$(INTDIR)\datio.obj" \
	"$(INTDIR)\erech.obj" \
	"$(INTDIR)\gdatastream.obj" \
	"$(INTDIR)\gstring.obj" \
	"$(INTDIR)\hydro1d_.obj" \
	"$(INTDIR)\inistat2.obj" \
	"$(INTDIR)\inparf.obj" \
	"$(INTDIR)\ipm_chemical.obj" \
	"$(INTDIR)\ipm_chemical2.obj" \
	"$(INTDIR)\ipm_chemical3.obj" \
	"$(INTDIR)\ipm_main.obj" \
	"$(INTDIR)\ipm_simplex.obj" \
	"$(INTDIR)\mainfromf.obj" \
	"$(INTDIR)\mcotac1d.obj" \
	"$(INTDIR)\ms_multi_file.obj" \
	"$(INTDIR)\ms_multi_format.obj" \
	"$(INTDIR)\ms_param.obj" \
	"$(INTDIR)\node.obj" \
	"$(INTDIR)\node_format.obj" \
	"$(INTDIR)\s_fgl.obj" \
	"$(INTDIR)\setpar.obj" \
	"$(INTDIR)\setpar_.obj" \
	"$(INTDIR)\walk2_.obj" \
	"$(INTDIR)\wegdat1d_.obj"

"$(OUTDIR)\test.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ENDIF 


!IF "$(NO_EXTERNAL_DEPS)" != "1"
!IF EXISTS("test.dep")
!INCLUDE "test.dep"
!ELSE 
!MESSAGE Warning: cannot find "test.dep"
!ENDIF 
!ENDIF 


!IF "$(CFG)" == "test - Win32 Release" || "$(CFG)" == "test - Win32 Debug"
SOURCE=..\datio.f

"$(INTDIR)\datio.obj" : $(SOURCE) "$(INTDIR)"
	$(F90) $(F90_PROJ) $(SOURCE)


SOURCE=..\erech.f

"$(INTDIR)\erech.obj" : $(SOURCE) "$(INTDIR)"
	$(F90) $(F90_PROJ) $(SOURCE)


SOURCE=..\gdatastream.cpp

"$(INTDIR)\gdatastream.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\gstring.cpp

"$(INTDIR)\gstring.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\hydro1d_.c

"$(INTDIR)\hydro1d_.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\inistat2.f

"$(INTDIR)\inistat2.obj" : $(SOURCE) "$(INTDIR)"
	$(F90) $(F90_PROJ) $(SOURCE)


SOURCE=..\inparf.f

"$(INTDIR)\inparf.obj" : $(SOURCE) "$(INTDIR)"
	$(F90) $(F90_PROJ) $(SOURCE)


SOURCE=..\ipm_chemical.cpp

"$(INTDIR)\ipm_chemical.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\ipm_chemical2.cpp

"$(INTDIR)\ipm_chemical2.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\ipm_chemical3.cpp

"$(INTDIR)\ipm_chemical3.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\ipm_main.cpp

"$(INTDIR)\ipm_main.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\ipm_simplex.cpp

"$(INTDIR)\ipm_simplex.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\mainfromf.cpp

"$(INTDIR)\mainfromf.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\mcotac1d.f

"$(INTDIR)\mcotac1d.obj" : $(SOURCE) "$(INTDIR)"
	$(F90) $(F90_PROJ) $(SOURCE)


SOURCE=..\ms_multi_file.cpp

"$(INTDIR)\ms_multi_file.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\ms_multi_format.cpp

"$(INTDIR)\ms_multi_format.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\ms_param.cpp

"$(INTDIR)\ms_param.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\node.cpp

"$(INTDIR)\node.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\node_format.cpp

"$(INTDIR)\node_format.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\s_fgl.cpp

"$(INTDIR)\s_fgl.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\setpar.cpp

"$(INTDIR)\setpar.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\setpar_.c

"$(INTDIR)\setpar_.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\walk2_.c

"$(INTDIR)\walk2_.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\wegdat1d_.c

"$(INTDIR)\wegdat1d_.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)



!ENDIF 

