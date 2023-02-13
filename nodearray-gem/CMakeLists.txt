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
    add_library(GEM2MT_SHARED SHARED $<TARGET_OBJECTS:GEM2MT_OBJECT>)
    set_target_properties(GEM2MT_SHARED PROPERTIES OUTPUT_NAME GEM2MT)
    target_link_libraries(GEM2MT_SHARED PUBLIC GEMS3K)
    #if(USE_THERMOFUN MATCHES ON)
    #  target_link_libraries(GEM2MT_SHARED PUBLIC ThermoFun::ThermoFun)
    #endif()
    install( TARGETS GEM2MT_SHARED
             EXPORT GEM2MTTargets DESTINATION "lib" COMPONENT libraries)
endif()


# Check if a static library is to be built
if(BUILD_STATIC_LIBS)
   add_library(GEM2MT_STATIC STATIC $<TARGET_OBJECTS:GEM2MT_OBJECT>)
   set_target_properties(GEM2MT_STATIC PROPERTIES OUTPUT_NAME GEM2MT-static)
   target_link_libraries(GEM2MT_STATIC PUBLIC GEMS3K-static)
   #if(USE_THERMOFUN MATCHES ON)
   #  target_link_libraries(GEM2MT_STATIC PUBLIC ThermoFun::ThermoFun)
   #endif()
   install( TARGETS GEM2MT_STATIC
           EXPORT GEM2MTTargets DESTINATION "lib" COMPONENT libraries)
endif()


# Install debug symbols
if(MSVC)
    install(
        FILES $<TARGET_PDB_FILE:GEM2MT_SHARED>
        DESTINATION ${CMAKE_INSTALL_BINDIR}
        COMPONENT libraries
        OPTIONAL
    )
endif()

# Create an install target for the header files
install(DIRECTORY ${GEM2MT_SOURCE_DIR}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} COMPONENT headers
    FILES_MATCHING PATTERN "*.h")

