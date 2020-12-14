# options
message(STATUS "Building json <-> key-value converter ")

include_directories(${CMAKE_SOURCE_DIR})

add_executable(kva2json kva2json.cpp)

if(${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
    target_link_libraries( kva2json PUBLIC GEMS3K-static )
else()
    target_link_libraries(kva2json PUBLIC GEMS3K )
endif()

install(TARGETS kva2json DESTINATION bin)

