# options
message(STATUS "Building isolate tsolmod4rkt library ")
file(COPY "${CMAKE_SOURCE_DIR}/tsolmod4rkt-build/Thermo-time-in" DESTINATION "${CMAKE_BINARY_DIR}/bin")
file(COPY "${CMAKE_SOURCE_DIR}/tsolmod4rkt-build/Thermo-time-all" DESTINATION "${CMAKE_BINARY_DIR}/bin")
file(COPY "${CMAKE_SOURCE_DIR}/tsolmod4rkt/solmod4rkt_demo.py" DESTINATION "${CMAKE_BINARY_DIR}/bin")
file(COPY "${CMAKE_SOURCE_DIR}/tsolmod4rkt/solmod4rkt_test01.py" DESTINATION "${CMAKE_BINARY_DIR}/bin")
file(COPY "${CMAKE_SOURCE_DIR}/tsolmod4rkt/solmod4rkt_test01.ipynb" DESTINATION "${CMAKE_BINARY_DIR}/bin")


include_directories(${CMAKE_SOURCE_DIR})

set(SOLMOD_HEADER_FILES
                    ${CMAKE_SOURCE_DIR}/GEMS3K/verror.h
                    ${CMAKE_SOURCE_DIR}/GEMS3K/gdatastream.h
                    ${CMAKE_SOURCE_DIR}/GEMS3K/num_methods.h
                    ${CMAKE_SOURCE_DIR}/GEMS3K/s_solmod.h
                    ${CMAKE_SOURCE_DIR}/GEMS3K/m_const_base.h
                    ${CMAKE_SOURCE_DIR}/GEMS3K/databr.h
                    ${CMAKE_SOURCE_DIR}/GEMS3K/datach.h
                    ${CMAKE_SOURCE_DIR}/GEMS3K/io_template.h
                    ${CMAKE_SOURCE_DIR}/GEMS3K/io_nlohmann.h
                    ${CMAKE_SOURCE_DIR}/GEMS3K/io_keyvalue.h
                    ${CMAKE_SOURCE_DIR}/GEMS3K/io_simdjson.h
                    ${CMAKE_SOURCE_DIR}/GEMS3K/gems3k_impex.h
                    ${CMAKE_SOURCE_DIR}/GEMS3K/v_detail.h
                    ${CMAKE_SOURCE_DIR}/GEMS3K/v_service.h
                    ${CMAKE_SOURCE_DIR}/GEMS3K/jsonconfig.h
                    ${CMAKE_SOURCE_DIR}/GEMS3K/datach_api.h
                    solmodengine.h
                    solmodfactory.h )


set(SOLMOD_SOURCE_FILES
                      ${CMAKE_SOURCE_DIR}/GEMS3K/gdatastream.cpp
                      ${CMAKE_SOURCE_DIR}/GEMS3K/num_methods.cpp
                      ${CMAKE_SOURCE_DIR}/GEMS3K/s_solmod.cpp
                      ${CMAKE_SOURCE_DIR}/GEMS3K/s_solmod2.cpp
                      ${CMAKE_SOURCE_DIR}/GEMS3K/s_solmod3.cpp
                      ${CMAKE_SOURCE_DIR}/GEMS3K/s_solmod4.cpp
                      ${CMAKE_SOURCE_DIR}/GEMS3K/s_solmod5.cpp
                      ${CMAKE_SOURCE_DIR}/GEMS3K/s_solmod6.cpp
                      ${CMAKE_SOURCE_DIR}/GEMS3K/s_sorpmod.cpp
                      ${CMAKE_SOURCE_DIR}/GEMS3K/io_template.cpp
                      ${CMAKE_SOURCE_DIR}/GEMS3K/io_nlohmann.cpp
                      ${CMAKE_SOURCE_DIR}/GEMS3K/io_keyvalue.cpp
                      ${CMAKE_SOURCE_DIR}/GEMS3K/io_simdjson.cpp
                      ${CMAKE_SOURCE_DIR}/GEMS3K/gems3k_impex.cpp
                      ${CMAKE_SOURCE_DIR}/GEMS3K/v_detail.cpp
                      ${CMAKE_SOURCE_DIR}/GEMS3K/v_service.cpp
                      ${CMAKE_SOURCE_DIR}/GEMS3K/jsonconfig.cpp
                      ${CMAKE_SOURCE_DIR}/GEMS3K/datach_api.cpp
                      ${CMAKE_SOURCE_DIR}/GEMS3K/datach_formats.cpp
                      solmodengine.cpp
                      solmodfactory.cpp
                      solmodfactory_alloc.cpp
                      solmodfactory_format.cpp )



# Compile the C++ source code into object files
add_library(SOLMOD_OBJECT OBJECT ${SOLMOD_HEADER_FILES} ${SOLMOD_SOURCE_FILES})

#if(USE_SPDLOG_PRECOMPILED)
#   if(NOT TARGET spdlog::spdlog)
#       find_package(spdlog CONFIG REQUIRED)
#   endif()
#endif()


# Check if a shared library is to be built
if(BUILD_SHARED_LIBS)
    # Enable automatic creation of a module definition (.def) file for a SHARED library on Windows.
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
    add_library(tsolmod4rkt SHARED $<TARGET_OBJECTS:SOLMOD_OBJECT>)
    add_library(GEMS3K::tsolmod4rkt ALIAS tsolmod4rkt)
    if(USE_THERMOFUN MATCHES ON)
      target_compile_definitions(tsolmod4rkt PUBLIC USE_SPDLOG_PRECOMPILED)
      target_link_libraries(tsolmod4rkt PUBLIC ThermoFun::ThermoFun)
    else()
      if(USE_SPDLOG_PRECOMPILED)
         target_compile_definitions(tsolmod4rkt PUBLIC USE_SPDLOG_PRECOMPILED)
         target_link_libraries(tsolmod4rkt PUBLIC spdlog::spdlog)
      endif()
    endif()
    install( TARGETS tsolmod4rkt
             EXPORT GEMS3KTargets DESTINATION "lib" COMPONENT libraries)

    # Install debug symbols
    if(MSVC)
    install(
        FILES $<TARGET_PDB_FILE:solmod>
        DESTINATION ${CMAKE_INSTALL_BINDIR}
        COMPONENT libraries
        OPTIONAL)
    endif()
endif()

# Check if a static library is to be built
if(BUILD_STATIC_LIBS)
   add_library(tsolmod4rkt-static STATIC $<TARGET_OBJECTS:SOLMOD_OBJECT>)
    add_library(GEMS3K::tsolmod4rkt-static ALIAS tsolmod4rkt-static)
   if(USE_THERMOFUN MATCHES ON)
     target_link_libraries(tsolmod4rkt-static PUBLIC ThermoFun::ThermoFun)
   else()
       if(USE_SPDLOG_PRECOMPILED)
          target_compile_definitions(tsolmod4rkt-static PUBLIC USE_SPDLOG_PRECOMPILED)
          target_link_libraries(tsolmod4rkt-static PUBLIC spdlog::spdlog)
       endif()
   endif()
   install( TARGETS tsolmod4rkt-static
           EXPORT GEMS3KTargets DESTINATION "lib" COMPONENT libraries)
endif()

# Install solmod header files
install(FILES ${SOLMOD_HEADER_FILES}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/tsolmod4rkt )

# Test application
find_package(Threads REQUIRED)
set(Test_apps
       solmod4rkt_test01
       solmodfactory_demo
       solmodengine_demo )

foreach (appName ${Test_apps})
    add_executable(${appName} ${appName}.cpp)
    if(${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
        target_link_libraries(${appName} PUBLIC GEMS3K::tsolmod4rkt-static Threads::Threads)
    else()
        target_link_libraries(${appName} PUBLIC GEMS3K::tsolmod4rkt Threads::Threads)
    endif()
    #install(TARGETS ${appName}} DESTINATION bin)
endforeach()

# Find the Python interpreter
#if(THERMOMATCH_BUILD_PYTHON)
    find_package(pybind11 REQUIRED)
    message(STATUS "Found pybind11 v${pybind11_VERSION}: ${pybind11_INCLUDE_DIRS}")
    if(NOT pybind11_FOUND)
        set(CHEMICALFUN_BUILD_PYTHON OFF)
        message(WARNING "pybind11 has not been found. "
            "Cannot build python bindings for thermomatch")
    endif()

    add_subdirectory(python)
#endif()



