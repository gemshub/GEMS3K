# options
message(STATUS "Building GEMS3K-server")

include_directories(${CMAKE_SOURCE_DIR})

add_executable(gems3k-server tnodetask.cpp main.cpp)

if(${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
    target_link_libraries(gems3k-server 
        PUBLIC GEMS3K-static
        PUBLIC libzmq )
else()
    target_link_libraries(gems3k-server 
        PUBLIC GEMS3K
        PUBLIC zmq )
endif()

install(TARGETS gems3k-server DESTINATION bin)

