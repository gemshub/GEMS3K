# Recursively collect all header files from the current directory
file(GLOB HEADER_FILES *.h)

# Recursively collect all source files from the current directory
file(GLOB SOURCE_FILES *.cpp)

# Recursively collect all source files from the current directory
#file(GLOB SIMDJSON_FILES ${CMAKE_SOURCE_DIR}/simdjson/simdjson.cpp)

# Compile the C++ source code into object files
#add_library(GEMS3K_OBJECT OBJECT ${HEADER_FILES} ${SOURCE_FILES})

if(USE_SPDLOG_PRECOMPILED)
   if(NOT TARGET spdlog::spdlog)
       find_package(spdlog CONFIG REQUIRED)
   endif()
endif()

# General gate for benchmark/diagnostics-only instrumentation (condition-number
# estimation, per-phase timing, and future additions in the same spirit —
# structured failure counters, stall/convergence-trajectory export, etc., per
# GEMS3K/CLAUDE.md's sequencing plan). This code adds real overhead per call
# (e.g. condition-number estimation costs ~44% of the linear-solve step itself
# — see ipm_main.cpp's MakeAndSolveSystemOfLinearEquations()), so it's OFF by
# default: normal/production use of GEMS3K pays nothing extra. gems-benchmark's
# build turns this ON to get per-system diagnostics. New instrumentation of
# this kind should reuse GEMS3K_BENCHMARK_DIAGNOSTICS rather than adding a new
# macro per feature.
option(ENABLE_BENCHMARK_DIAGNOSTICS "Enable benchmark/diagnostics-only instrumentation (condition-number estimation, per-phase timing, etc.) in the IPM solver (adds overhead; for benchmarking/diagnostics only)" OFF)

# Check if a shared library is to be built
if(BUILD_SHARED_LIBS)
   # Enable automatic creation of a module definition (.def) file for a SHARED library on Windows.
   set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS TRUE)

   add_library(GEMS3K ${HEADER_FILES} ${SOURCE_FILES})
   add_library(GEMS3K::GEMS3K ALIAS GEMS3K)
   if(ENABLE_BENCHMARK_DIAGNOSTICS)
       target_compile_definitions(GEMS3K PRIVATE GEMS3K_BENCHMARK_DIAGNOSTICS)
   endif()
   if(USE_THERMOFUN MATCHES ON)
       target_compile_definitions(GEMS3K PUBLIC USE_SPDLOG_PRECOMPILED)
       target_link_libraries(GEMS3K PUBLIC ThermoFun::ThermoFun)
   else()
       if(USE_SPDLOG_PRECOMPILED)
           target_compile_definitions(GEMS3K PUBLIC USE_SPDLOG_PRECOMPILED)
           target_link_libraries(GEMS3K PUBLIC spdlog::spdlog)
       endif()
   endif()

   install( TARGETS GEMS3K
             EXPORT GEMS3KTargets
             ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT libraries
             LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT libraries
             RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT libraries)

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
   add_library(GEMS3K-static STATIC ${HEADER_FILES} ${SOURCE_FILES})
    add_library(GEMS3K::GEMS3K-static ALIAS GEMS3K-static)
   if(ENABLE_BENCHMARK_DIAGNOSTICS)
       target_compile_definitions(GEMS3K-static PRIVATE GEMS3K_BENCHMARK_DIAGNOSTICS)
   endif()
   if(USE_THERMOFUN MATCHES ON)
     target_link_libraries(GEMS3K-static PUBLIC ThermoFun::ThermoFun)
   endif()
   #if(USE_SPDLOG_PRECOMPILED)
   #   target_compile_definitions(GEMS3K-static PUBLIC USE_SPDLOG_PRECOMPILED)
   #   target_link_libraries(GEMS3K-static PUBLIC spdlog::spdlog)
   #endif()
   install( TARGETS GEMS3K-static
           EXPORT GEMS3KTargets DESTINATION "lib" COMPONENT libraries)
endif()

# Install GEMS3K header files
install(DIRECTORY ${PROJECT_SOURCE_DIR}/GEMS3K
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} COMPONENT headers
    FILES_MATCHING PATTERN "*.hpp" PATTERN "*.h")



