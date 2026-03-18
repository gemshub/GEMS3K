# Require a certain version of cmake
cmake_minimum_required(VERSION 3.16..4.2)

# Set the name of the project
project(GEMS3K VERSION 4.5.2 LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
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

# Define which types of libraries to build
option(BUILD_SHARED_LIBS "Build shared libraries." ON)
option(BUILD_STATIC_LIBS "Build static libraries." ON)
option(USE_THERMOFUN "Use ThermoFun" ON)
option(BUILD_GEM2MT "Install GEM2MT lib." OFF)
option(BUILD_NLOHMANNJSON "Use nlohmann json." OFF)
option(BUILD_TOOLS "Build tools" ON)
option(BUILD_SOLMOD "Build isolate TSolMod " ON)
option(BUILD_SOLMOD_PYTHON "Build isolate TSolMod " ON)

# Currently is not setup to produce a dynamic library using MSVC, only static
#if(${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
#    set(BUILD_SHARED_LIBS OFF)
#endif()

# Used into conda only
if(DEFINED ENV{CONDA_PREFIX})
   option(USE_SPDLOG_PRECOMPILED "Use spdlog in compiled version" ON)
else()
   option(USE_SPDLOG_PRECOMPILED "Use spdlog in compiled version" OFF)
endif()

# Set the default build type to Release
if(NOT CMAKE_BUILD_TYPE)
    message(STATUS "Setting build type to 'Release' as none was specified.")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# Set the list of compiler flags for GNU compiler
if(${CMAKE_CXX_COMPILER_ID} MATCHES "GNU")
    add_compile_options( -pthread -Wall -Wno-misleading-indentation -Wno-ignored-attributes -Wno-variadic-macros -Wno-deprecated)
endif()

# Set the list of compiler flags for Clang compiler
if(${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
    add_compile_options(-Wall -Wno-ignored-attributes -Wno-pedantic -Wno-variadic-macros -Wno-deprecated)
endif()

# Set the list of compiler flags for Intel compiler
if(${CMAKE_CXX_COMPILER_ID} MATCHES "Intel")
    add_compile_options(-Wall -Wno-variadic-macros -Wno-deprecated)
endif()

# Set the list of compiler flags for MSVC compiler
if(${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
    add_compile_options(
        /D_SCL_SECURE_NO_WARNINGS
        /D_CRT_SECURE_NO_WARNINGS=1
        /MP4
        /utf-8
        /EHsc
        /D_SILENCE_TR1_NAMESPACE_DEPRECATION_WARNING
        /DNOMINMAX
    )
endif()

# test thread safe
#set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=thread")
#set (CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG} -fsanitize=thread")

# Set the source directory path
set(GEMS3K_SOURCE_DIR ${CMAKE_SOURCE_DIR}/GEMS3K)
set(GEMS3K_HEADER_DIR ${CMAKE_SOURCE_DIR}/GEMS3K)

if(${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
    include_directories("/usr/local/include")
    link_directories("/usr/local/lib")
endif()

# Set the include directories
include_directories(${GEMS3K_HEADER_DIR})
include_directories(${CMAKE_SOURCE_DIR})

# Find all GEMS3K dependencies
include(GEMS3KFindDeps)

# Versioning
set(GEMS3K_VERSION "${PROJECT_VERSION}")
find_package(Git)
if(GIT_FOUND)
  execute_process(
    COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    OUTPUT_VARIABLE GEMS3K_VERSION_HASH
    ERROR_QUIET
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )
  execute_process(
    COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    OUTPUT_VARIABLE GEMS3K_GIT_BRANCH
    ERROR_QUIET
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if(NOT GEMS3K_VERSION_HASH)
        set(GEMS3K_VERSION_HASH "---")
    endif()
    if(NOT GEMS3K_GIT_BRANCH)
        set(GEMS3K_GIT_BRANCH "---")
    endif()
else()
  set(GEMS3K_VERSION_HASH "<unknown>")
  set(GEMS3K_GIT_BRANCH "<unknown>")
endif()

if(USE_THERMOFUN MATCHES ON)
   set(GEMS3K_ChemicalFun "${ChemicalFun_VERSION}")
   set(GEMS3K_ThermoFun "${ThermoFun_VERSION}")
else()
   set(GEMS3K_ChemicalFun "-")
   set(GEMS3K_ThermoFun "-")
endif()

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/gems3k_version.h.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/GEMS3K/gems3k_version.h"
  @ONLY
  )

# Set some necessary definitions
if(USE_THERMOFUN MATCHES OFF)
  add_definitions(-DNO_THERMOFUN)
else()
  #add_definitions(-DUSE_THERMO_LOG)
endif()

if(BUILD_NLOHMANNJSON MATCHES ON)
  add_definitions(-DUSE_NLOHMANNJSON)
endif()

if(BUILD_GEM2MT MATCHES ON)
  set(GEM2MT_SOURCE_DIR ${CMAKE_SOURCE_DIR}/nodearray-gem)
  set(GEM2MT_HEADER_DIR ${CMAKE_SOURCE_DIR}/nodearray-gem)
  include_directories(${GEM2MT_HEADER_DIR})
endif()


# Compile GEMS3K
add_subdirectory(GEMS3K)
link_directories(${CMAKE_BINARY_DIR}/GEMS3K)
link_directories(${CMAKE_BINARY_DIR}/lib)

if(BUILD_GEM2MT MATCHES ON)
    add_subdirectory(nodearray-gem)
    link_directories(${CMAKE_BINARY_DIR}/nodearray-gem)
endif()

if(BUILD_TOOLS)
   add_subdirectory(tools)
endif()

if(BUILD_SOLMOD)
   add_subdirectory(tsolmod4rkt)
endif()

# Install the cmake config files that permit users to use find_package(gems3k)
include(GEMS3KInstallCMakeConfigFiles)
