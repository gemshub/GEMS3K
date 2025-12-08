message(STATUS "Building python library ")
file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/__init__.py" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")

# Collect the source files of PyThermoFun
file(GLOB_RECURSE HEADER_FILES *.hpp)
file(GLOB_RECURSE SOURCE_FILES *.cpp)

# Include the python directory
include_directories(${PROJECT_SOURCE_DIR}/tsolmod4rkt)

# Create the PySolmod4rkt library
pybind11_add_module(PySolmod4rkt ${SOURCE_FILES} NO_EXTRAS) # NO_EXTRAS used to decrease linking time, but causes larger module file size

if(${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
    target_link_libraries(PySolmod4rkt LINK_PRIVATE GEMS3K::tsolmod4rkt-static)
else()
    target_link_libraries(PySolmod4rkt LINK_PRIVATE GEMS3K::tsolmod4rkt)
endif()
