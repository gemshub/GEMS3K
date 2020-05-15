# Require a certain version of cmake
cmake_minimum_required(VERSION 3.7)

# Set the name of the project
project(GEMS3K VERSION 4.0.0 LANGUAGES CXX C)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_MACOSX_RPATH ON)

# Include the cmake variables with values for installation directories
include(GNUInstallDirs)

# Set the cmake module path of the project
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules")

# Use ccache to speed up repeated compilations
include(CCache)

# Ensure proper configuration if in a conda environment
include(CondaAware)

# Set the output directories of the built libraries and binaries
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Set libraries to be compiled with position independent code mode (i.e., fPIC option in GNU compilers)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Set the default build type to Release
if(NOT CMAKE_BUILD_TYPE)
    message(STATUS "Setting build type to 'Release' as none was specified.")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# Define which types of libraries to build
option(BUILD_SHARED_LIBS "Build shared libraries." ON)
option(BUILD_STATIC_LIBS "Build static libraries." ON)

# Set the list of compiler flags for GNU compiler
if(${CMAKE_CXX_COMPILER_ID} MATCHES "GNU")
    add_compile_options(-std=c++11 -Wall -Wno-pedantic -Wno-variadic-macros -Wno-deprecated)
endif()

# Set the list of compiler flags for Clang compiler
if(${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
    add_compile_options(-std=c++11 -Wall -Wno-ignored-attributes -Wno-pedantic -Wno-variadic-macros -Wno-deprecated)
endif()

# Set the list of compiler flags for Intel compiler
if(${CMAKE_CXX_COMPILER_ID} MATCHES "Intel")
    add_compile_options(-std=c++11 -Wall -Wno-variadic-macros -Wno-deprecated)
endif()

# Set the list of compiler flags for MSVC compiler
if(${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
    add_compile_options(
        /D_SCL_SECURE_NO_WARNINGS
        /D_CRT_SECURE_NO_WARNINGS=1
        /MP4
        /EHsc
        /D_SILENCE_TR1_NAMESPACE_DEPRECATION_WARNING
        /DNOMINMAX
    )
endif()


# Set the jsonio14 source directory path
set(GEMS3K_SOURCE_DIR ${CMAKE_SOURCE_DIR}/GEMS3K)
set(GEMS3K_HEADER_DIR ${CMAKE_SOURCE_DIR}/GEMS3K)

# Set the include directories
include_directories(${GEMS3K_HEADER_DIR})
if(${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
    include_directories("/usr/local/include")
    link_directories("/usr/local/lib")
endif()

# Find all GEMS3K dependencies
include(GEMS3KFindDeps)

# Compile GEMS3K
add_subdirectory(GEMS3K)
link_directories(${CMAKE_BINARY_DIR}/GEMS3K)
link_directories(${CMAKE_BINARY_DIR}/lib)

# Compile GEMS3K
add_subdirectory(gemserver)

# Install the cmake config files that permit users to use find_package(gems3k)
#include(GEMS3KInstallCMakeConfigFiles)
