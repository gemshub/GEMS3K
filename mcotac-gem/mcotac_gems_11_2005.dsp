# Microsoft Developer Studio Project File - Name="mcotac_gems_11_2005" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=mcotac_gems_11_2005 - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "mcotac_gems_11_2005.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "mcotac_gems_11_2005.mak" CFG="mcotac_gems_11_2005 - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "mcotac_gems_11_2005 - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "mcotac_gems_11_2005 - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "mcotac_gems_11_2005 - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
# ADD BASE F90 /compile_only /include:"Release/" /nologo /warn:nofileopt
# ADD F90 /compile_only /include:"Release/" /nologo /warn:nofileopt
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD BASE RSC /l 0x807 /d "NDEBUG"
# ADD RSC /l 0x807 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386

!ELSEIF  "$(CFG)" == "mcotac_gems_11_2005 - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE F90 /check:bounds /compile_only /debug:full /include:"Debug/" /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD F90 /check:bounds /compile_only /debug:full /include:"Debug/" /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /D "IPMGEMPLUGIN" /YX /FD /GZ /c
# ADD BASE RSC /l 0x807 /d "_DEBUG"
# ADD RSC /l 0x807 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /incremental:no /debug /machine:I386 /nodefaultlib:"libcd libcpd libc" /force /pdbtype:sept

!ENDIF 

# Begin Target

# Name "mcotac_gems_11_2005 - Win32 Release"
# Name "mcotac_gems_11_2005 - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat;f90;for;f;fpp"
# Begin Source File

SOURCE=.\datio.f
DEP_F90_DATIO=\
	".\gwheader.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\erech.f
DEP_F90_ERECH=\
	".\gwheader.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\gdatastream.cpp
# End Source File
# Begin Source File

SOURCE=.\gstring.cpp
# End Source File
# Begin Source File

SOURCE=.\hydro1d_.c
# End Source File
# Begin Source File

SOURCE=.\inistat2.f
DEP_F90_INIST=\
	".\gwheader.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\inparf.f
DEP_F90_INPAR=\
	".\gwheader.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\ipm_chemical.cpp
# End Source File
# Begin Source File

SOURCE=.\ipm_chemical2.cpp
# End Source File
# Begin Source File

SOURCE=.\ipm_chemical3.cpp
# End Source File
# Begin Source File

SOURCE=.\ipm_main.cpp
# End Source File
# Begin Source File

SOURCE=.\ipm_simplex.cpp
# End Source File
# Begin Source File

SOURCE=.\mainfromf.cpp
# End Source File
# Begin Source File

SOURCE=.\mcotac1d.f
DEP_F90_MCOTA=\
	".\gwheader.inc"\
	
NODEP_F90_MCOTA=\
	".\Debug\datach.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\ms_multi_file.cpp
# End Source File
# Begin Source File

SOURCE=.\ms_multi_format.cpp
# End Source File
# Begin Source File

SOURCE=.\ms_param.cpp
# End Source File
# Begin Source File

SOURCE=.\node.cpp
# End Source File
# Begin Source File

SOURCE=.\node_format.cpp
# End Source File
# Begin Source File

SOURCE=.\s_fgl.cpp
# End Source File
# Begin Source File

SOURCE=.\setpar.cpp
# End Source File
# Begin Source File

SOURCE=.\setpar_.c
# End Source File
# Begin Source File

SOURCE=.\walk2_.c
# End Source File
# Begin Source File

SOURCE=.\wegdat1d_.c
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl;fi;fd"
# Begin Source File

SOURCE=.\array.h
# End Source File
# Begin Source File

SOURCE=.\databr.h
# End Source File
# Begin Source File

SOURCE=.\datach.h
# End Source File
# Begin Source File

SOURCE=.\datach.inc
# End Source File
# Begin Source File

SOURCE=.\gdatastream.h
# End Source File
# Begin Source File

SOURCE=.\gstring.h
# End Source File
# Begin Source File

SOURCE=.\gwheader.h
# End Source File
# Begin Source File

SOURCE=.\jama_cholesky.h
# End Source File
# Begin Source File

SOURCE=.\jama_lu.h
# End Source File
# Begin Source File

SOURCE=.\m_const.h
# End Source File
# Begin Source File

SOURCE=.\m_param.h
# End Source File
# Begin Source File

SOURCE=.\ms_multi.h
# End Source File
# Begin Source File

SOURCE=.\node.h
# End Source File
# Begin Source File

SOURCE=.\s_fgl.h
# End Source File
# Begin Source File

SOURCE=.\tnt.h
# End Source File
# Begin Source File

SOURCE=.\tnt_array1d.h
# End Source File
# Begin Source File

SOURCE=.\tnt_array2d.h
# End Source File
# Begin Source File

SOURCE=.\tnt_i_refvec.h
# End Source File
# Begin Source File

SOURCE=.\v_user.h
# End Source File
# Begin Source File

SOURCE=.\verror.h
# End Source File
# End Group
# Begin Group "Resource Files"

# PROP Default_Filter "ico;cur;bmp;dlg;rc2;rct;bin;rgs;gif;jpg;jpeg;jpe"
# End Group
# Begin Source File

SOURCE=.\gwheader.inc
# End Source File
# Begin Source File

SOURCE=.\kinetics.inc
# End Source File
# End Target
# End Project
