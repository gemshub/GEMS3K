# Find pybind11
#find_package(pybind11 REQUIRED)

# Build the Python wrappers of tsolmod4rkt library
add_subdirectory(PySolmod4rkt)

# Build the Python package for built-in mixing models
add_subdirectory(solmod4rkt)
