# options
message(STATUS "Building GEMS3K-server")

include_directories(${CMAKE_SOURCE_DIR})

add_executable(gems3k-server tnodetask.cpp main.cpp)
target_link_libraries(gems3k-server 
PUBLIC GEMS3K
PUBLIC zmq )

install(TARGETS gems3k-server DESTINATION bin)

