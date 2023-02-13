# options
message(STATUS "Building json <-> key-value converter ")

include_directories(${CMAKE_SOURCE_DIR})

add_executable(kva2json kva2json.cpp)

if(BUILD_GEM2MT MATCHES ON)
  if(${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
    target_link_libraries(kva2json PUBLIC GEMS3K-static GEM2MT-static)
  else()
    target_link_libraries(kva2json PUBLIC GEMS3K GEM2MT pthread)
  endif()
else()
  if(${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
    target_link_libraries( kva2json PUBLIC GEMS3K-static)
  else()
    target_link_libraries(kva2json PUBLIC GEMS3K pthread)
  endif()
endif()

#if(USE_THERMOFUN MATCHES ON)
#   target_link_libraries(kva2json PUBLIC ThermoFun::ThermoFun)
#endif()

install(TARGETS kva2json DESTINATION bin)

