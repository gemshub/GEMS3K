# options
message(STATUS "Building json <-> key-value converter ")

include_directories(${CMAKE_SOURCE_DIR})

find_package(Threads REQUIRED)

set(Tools
       kva2json
       thread_test
       thermofun_test)

foreach (toolName ${Tools})

    add_executable(${toolName} ${toolName}.cpp)

    if(BUILD_GEM2MT MATCHES ON)
      if(${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
          target_link_libraries(${toolName} PUBLIC GEMS3K::GEMS3K-static GEMS3K::GEM2MT-static Threads::Threads)
      else()
          target_link_libraries(${toolName} PUBLIC GEMS3K::GEMS3K GEMS3K::GEM2MT Threads::Threads)
      endif()
    else()
      if(${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
          target_link_libraries(${toolName} PUBLIC GEMS3K::GEMS3K-static Threads::Threads)
      else()
          target_link_libraries(${toolName} PUBLIC GEMS3K::GEMS3K Threads::Threads)
      endif()
    endif()

    install(TARGETS ${toolName} DESTINATION bin)

endforeach()

