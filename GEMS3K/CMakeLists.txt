# Recursively collect all header files from the current directory
file(GLOB HEADER_FILES *.h)

# Recursively collect all source files from the current directory
file(GLOB SOURCE_FILES *.cpp)

# Recursively collect all source files from the current directory
#file(GLOB SIMDJSON_FILES ${CMAKE_SOURCE_DIR}/simdjson/simdjson.cpp)

# Compile the C++ source code into object files
add_library(GEMS3K_OBJECT OBJECT ${HEADER_FILES} ${SOURCE_FILES})

if(USE_SPDLOG_PRECOMPILED)
   if(NOT TARGET spdlog::spdlog)
       find_package(spdlog CONFIG REQUIRED)
   endif()
endif()

# Check if a shared library is to be built
if(BUILD_SHARED_LIBS)
    # Enable automatic creation of a module definition (.def) file for a SHARED library on Windows.
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
    add_library(GEMS3K SHARED $<TARGET_OBJECTS:GEMS3K_OBJECT>)
    add_library(GEMS3K::GEMS3K ALIAS GEMS3K)
    if(USE_THERMOFUN MATCHES ON)
      target_link_libraries(GEMS3K ThermoFun::ThermoFun)
    else()
      if(USE_SPDLOG_PRECOMPILED)
         target_compile_definitions(GEMS3K PUBLIC USE_SPDLOG_PRECOMPILED)
         target_link_libraries(GEMS3K PUBLIC spdlog::spdlog)
      endif()
    endif()
    install( TARGETS GEMS3K
             EXPORT GEMS3KTargets DESTINATION "lib" COMPONENT libraries)

    # Install debug symbols
    if(MSVC)
    install(
        FILES $<TARGET_PDB_FILE:GEMS3K>
        DESTINATION ${CMAKE_INSTALL_BINDIR}
        COMPONENT libraries
        OPTIONAL)
    endif()
endif()

# Check if a static library is to be built
if(BUILD_STATIC_LIBS)
   add_library(GEMS3K-static STATIC $<TARGET_OBJECTS:GEMS3K_OBJECT>)
    add_library(GEMS3K::GEMS3K-static ALIAS GEMS3K-static)
   if(USE_THERMOFUN MATCHES ON)
     target_link_libraries(GEMS3K-static PUBLIC ThermoFun::ThermoFun)
   else()
       if(USE_SPDLOG_PRECOMPILED)
          target_compile_definitions(GEMS3K-static PUBLIC USE_SPDLOG_PRECOMPILED)
          target_link_libraries(GEMS3K-static PUBLIC spdlog::spdlog)
       endif()
   endif()
   install( TARGETS GEMS3K-static
           EXPORT GEMS3KTargets DESTINATION "lib" COMPONENT libraries)
endif()



# Install GEMS3K header files
install(DIRECTORY ${PROJECT_SOURCE_DIR}/GEMS3K
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} COMPONENT headers
    FILES_MATCHING PATTERN "*.hpp" PATTERN "*.h")



