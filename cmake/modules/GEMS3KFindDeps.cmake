#find_package(ZeroMQ REQUIRED)
#find_library(ZEROMQ_LIB zeromq)
#if(NOT ZEROMQ_LIB)
#  message(FATAL_ERROR "zeromq library not found")
#endif()

# Find pybind11 library (if needed)
if(BUILD_SOLMOD_PYTHON)
    find_package (Python3 COMPONENTS Interpreter Development)
    find_package(pybind11 REQUIRED)
    if(NOT pybind11_FOUND)
        message(WARNING "Could not find pybind11 - the Python module `thermofun` will not be built.")
        set(BUILD_SOLMOD_PYTHON OFF)
    else()
        message(STATUS "Found pybind11 v${pybind11_VERSION}: ${pybind11_INCLUDE_DIRS}")
    endif()
endif()

# Eigen (header-only) is needed by the native solver's linear algebra.
# With Optima, Optima is found first and brings the Eigen it was built with.
if(USE_OPTIMA_SOLVER)
  find_package(Optima REQUIRED)
  # GEMS3K needs the modified Optima: its Trace.hpp header and the jumpahead option.
  include(CheckCXXSourceCompiles)
  set(CMAKE_REQUIRED_LIBRARIES Optima::Optima)
  unset(GEMS3K_OPTIMA_IS_MODIFIED CACHE)   # checked again on every configure: Optima_DIR can change
  check_cxx_source_compiles("
    #include <Optima/Optima.hpp>
    #include <Optima/Trace.hpp>
    int main() { Optima::Options o; return o.jumpahead ? 1 : 0; }" GEMS3K_OPTIMA_IS_MODIFIED)
  unset(CMAKE_REQUIRED_LIBRARIES)
  if(NOT GEMS3K_OPTIMA_IS_MODIFIED)
    message(FATAL_ERROR "The Optima found in ${Optima_DIR} does not have Optima/Trace.hpp or Options::jumpahead. "
                        "Use the modified Optima of gemshub: https://github.com/gemshub/optima, branch optima_develop.")
  endif()
endif()
find_package(Eigen3 REQUIRED NO_MODULE)
if(Eigen3_VERSION VERSION_LESS 3.4)
  message(FATAL_ERROR "Eigen 3.4 or newer is required (found ${Eigen3_VERSION})")
endif()
message(STATUS "Found Eigen3 v${Eigen3_VERSION}")

if(USE_THERMOFUN MATCHES ON)
  find_package(ThermoFun REQUIRED)
  if(NOT ThermoFun_FOUND)
    message(FATAL_ERROR "ThermoFun library not found")
  else()
    message(STATUS "Found ThermoFun v${ThermoFun_VERSION}")
    # ThermoFun before 0.7.0 keeps a reference to the pressure in the key of its property cache, so a later
    # pressure changes the old keys. A ThermoFun before 0.7.0 cannot switch that cache off after it is made.
    if(ThermoFun_VERSION VERSION_LESS 0.7.0)
      message(WARNING "ThermoFun ${ThermoFun_VERSION} has a property cache that can return values of an earlier pressure. Use ThermoFun 0.7.0 or newer (see environment.yml).")
    endif()
    # ThermoFun 1.0 changes the default gas properties and the volume unit. This GEMS3K adds the pressure and
    # fugacity terms of gases itself and reads volumes in J/bar, so it needs a ThermoFun 0.x.
    if(NOT ThermoFun_VERSION VERSION_LESS 1.0)
      message(FATAL_ERROR "ThermoFun ${ThermoFun_VERSION} is not supported: this GEMS3K needs ThermoFun 0.7 or a later 0.x version.")
    endif()
  endif()

  find_package(ChemicalFun REQUIRED)
  if(NOT ChemicalFun_FOUND)
    message(FATAL_ERROR "ChemicalFun library not found")
  else()
    message(STATUS "Found ChemicalFun v${ChemicalFun_VERSION}")
  endif()

# Find the autodiff library (header-only); fetch it if it is not installed (not needed with -DTFUN_USE_AUTODIFF=OFF)
#if(NOT TFUN_USE_AUTODIFF)
#    message(STATUS "Building without autodiff: the models use plain numbers, all derivatives are 0")
#else()
find_package(autodiff 1.1.1 CONFIG QUIET)
if(autodiff_FOUND)
    message(STATUS "Found autodiff: ${autodiff_DIR} (found version \"${autodiff_VERSION}\")")
else()
    message(STATUS "autodiff not found - fetching it from GitHub")
    include(FetchContent)
    set(AUTODIFF_BUILD_TESTS    OFF CACHE BOOL "" FORCE)
    set(AUTODIFF_BUILD_PYTHON   OFF CACHE BOOL "" FORCE)
    set(AUTODIFF_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
    set(AUTODIFF_BUILD_DOCS     OFF CACHE BOOL "" FORCE)
    FetchContent_Declare(autodiff
        GIT_REPOSITORY https://github.com/autodiff/autodiff.git
        # v1.1.2 plus the commit that makes Eigen optional (ThermoFun does not need it)
        GIT_TAG        b0a4feff5b2a61262e94305452ac53369fe35e75)
    FetchContent_MakeAvailable(autodiff)
endif()
#endif() # TFUN_USE_AUTODIFF

  if(USE_SPDLOG_PRECOMPILED)
    if(NOT TARGET spdlog::spdlog)
       find_package(spdlog CONFIG REQUIRED)
       if(NOT spdlog_FOUND)
           message(FATAL_ERROR "spdlog library not found")
       else()
           message(STATUS "Found spdlog v${spdlog_VERSION}")
       endif()
   endif()
 endif()
endif()

