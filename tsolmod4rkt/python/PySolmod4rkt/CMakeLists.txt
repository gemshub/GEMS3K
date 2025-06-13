# Collect the source files of PySolmod4rkt
file(GLOB_RECURSE HEADER_FILES *.hpp)
file(GLOB_RECURSE SOURCE_FILES *.cpp)

# Include the python directory
include_directories(${PROJECT_SOURCE_DIR}/tsolmod4rkt)
include_directories(${PROJECT_SOURCE_DIR}/tsolmod4rkt/python)

# Create the PySolmod4rkt library
pybind11_add_module(PySolmod4rkt ${SOURCE_FILES} NO_EXTRAS) # NO_EXTRAS used to decrease linking time, but causes larger module file size

target_link_libraries(PySolmod4rkt LINK_PRIVATE GEMS3K::tsolmod4rkt)
