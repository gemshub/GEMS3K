# options
message(STATUS "Building json <-> key-value converter ")

include_directories(${CMAKE_SOURCE_DIR})

set(Tools
       kva2json
       thread_test )

foreach (toolName ${Tools})

    add_executable(${toolName} ${toolName}.cpp)

    if(BUILD_GEM2MT MATCHES ON)
      if(${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
          target_link_libraries(${toolName} PUBLIC GEMS3K-static GEM2MT-static pthread)
      else()
          target_link_libraries(${toolName} PUBLIC GEMS3K GEM2MT pthread)
      endif()
    else()
      if(${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
          target_link_libraries(${toolName} PUBLIC GEMS3K-static pthread)
      else()
          target_link_libraries(${toolName} PUBLIC GEMS3K pthread)
      endif()
    endif()

    install(TARGETS ${toolName} DESTINATION bin)

endforeach()

