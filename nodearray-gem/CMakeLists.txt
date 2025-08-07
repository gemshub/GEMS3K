# Recursively collect all header files from the current directory
#file(GLOB HEADER_FILES *.h)
set(HEADER_FILES
     "m_gem2mt.h"
     "particlearray.h")

# Recursively collect all source files from the current directory
#file(GLOB SOURCE_FILES *.cpp)
set(SOURCE_FILES
     "particlearray.cpp"
     "m_gem2mtsep.cpp"
     "m_gem2mtt.cpp"
     "m_gem2mtbox.cpp"
     "m_gem2mtfor.cpp")

# Compile the C++ source code into object files
add_library(GEM2MT_OBJECT OBJECT ${HEADER_FILES} ${SOURCE_FILES})

# Check if a shared library is to be built
if(BUILD_SHARED_LIBS)
    # Enable automatic creation of a module definition (.def) file for a SHARED library on Windows.
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
    add_library(GEM2MT SHARED $<TARGET_OBJECTS:GEM2MT_OBJECT>)
    add_library(GEMS3K::GEM2MT ALIAS GEM2MT)
    target_link_libraries(GEM2MT PUBLIC GEMS3K::GEMS3K)
    install( TARGETS GEM2MT
       EXPORT GEMS3KTargets
       ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT libraries
       LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT libraries
       RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT libraries)
    # Install debug symbols
    if(MSVC)
       install(
          FILES $<TARGET_PDB_FILE:GEM2MT>
          DESTINATION ${CMAKE_INSTALL_BINDIR}
          COMPONENT libraries
          OPTIONAL )
    endif()
endif()


# Check if a static library is to be built
if(BUILD_STATIC_LIBS)
   add_library(GEM2MT-static STATIC $<TARGET_OBJECTS:GEM2MT_OBJECT>)
   add_library(GEMS3K::GEM2MT-static ALIAS GEM2MT-static)
   target_link_libraries(GEM2MT-static PUBLIC GEMS3K::GEMS3K-static)
   install( TARGETS GEM2MT-static
      EXPORT GEMS3KTargets
      ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT libraries
      LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT libraries
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT libraries)
endif()

# Create an install target for the header files
install(DIRECTORY ${GEM2MT_SOURCE_DIR}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} COMPONENT headers
    FILES_MATCHING PATTERN "*.h")

