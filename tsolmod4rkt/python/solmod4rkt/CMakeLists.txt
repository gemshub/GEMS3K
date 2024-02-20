# Copy the thermo impex directory to the binary directory
file(COPY . DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Define a variable with the name of the compiled python module file as defined by pybind11 (which is used in the file setup.py.in)
set(SOLMOD4RKT_PYTHON_MODULE_FILENAME PySolmod4rkt${PYTHON_MODULE_EXTENSION})

# Configure the setup.py file
set(SETUP_PY_IN ${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in)
set(SETUP_PY ${CMAKE_CURRENT_BINARY_DIR}/setup.py)

# Configure the setup.py script with the provided CMake variables
configure_file(${SETUP_PY_IN} ${SETUP_PY})

# https://www.scivision.dev/cmake-install-python-package/
# https://pip.pypa.io/en/stable/topics/local-project-installs/#editable-installs
# detect virtualenv and set Pip args accordingly
if(DEFINED ENV{VIRTUAL_ENV} OR DEFINED ENV{CONDA_PREFIX})
  set(_PIP_ARG)
else()
  set(_PIP_ARG "--user")
endif()
#message(STATUS "_PIP_ARG: ${_PIP_ARG}")

# Create a custom target to build solmod4rkt python package during build stage
add_custom_target(solmod4rkt ALL
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:PySolmod4rkt>
        ${CMAKE_CURRENT_BINARY_DIR}/$<TARGET_FILE_NAME:PySolmod4rkt>
        COMMAND ${PYTHON_EXECUTABLE}  -m pip install --prefix ${CMAKE_BINARY_DIR} --ignore-installed .
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# Set dependencies of PySolmod4rkt target
add_dependencies(solmod4rkt PySolmod4rkt)

# Set the path where the python package is installed to CMAKE_INSTALL_PREFIX if not given
if(NOT DEFINED SOLMOD4RKT_PYTHON_INSTALL_PREFIX)

# Install the SOLMOD4RKT python package
install(CODE
"
    if(EXISTS ${CMAKE_CURRENT_BINARY_DIR}/../../../lib/PySolmod4rkt.pdb)
        string(REPLACE .pyd .pdb SOLMOD4RKT_PDB_FILENAME \"${SOLMOD4RKT_PYTHON_MODULE_FILENAME}\")

        execute_process(
            COMMAND \${CMAKE_COMMAND} -E copy ../../../lib/PySolmod4rkt.pdb \${SOLMOD4RKT_PDB_FILENAME}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    endif()

    execute_process(
        COMMAND ${PYTHON_EXECUTABLE} -m pip install ${CMAKE_CURRENT_BINARY_DIR}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
")
else()

# If the path is already in Windows format (with backslashes), it can't be added directly
# to the string below, otherwise CMake will later complain about "Invalid escape sequence".
file(TO_CMAKE_PATH "${SOLMOD4RKT_PYTHON_INSTALL_PREFIX}" SOLMOD4RKT_PYTHON_INSTALL_PREFIX)

# Install the SOLMOD4RKT python package using setuptools
install(CODE
"
    file(TO_NATIVE_PATH \"${SOLMOD4RKT_PYTHON_INSTALL_PREFIX}\" SOLMOD4RKT_PYTHON_INSTALL_PREFIX_NATIVE)

    if(EXISTS ${CMAKE_CURRENT_BINARY_DIR}/../../../lib/PySolmod4rkt.pdb)
        string(REPLACE .pyd .pdb SOLMOD4RKT_PDB_FILENAME \"${SOLMOD4RKT_PYTHON_MODULE_FILENAME}\")

        execute_process(
            COMMAND \${CMAKE_COMMAND} -E copy ../../../lib/PySolmod4rkt.pdb \${SOLMOD4RKT_PDB_FILENAME}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    endif()

    execute_process(
    COMMAND ${PYTHON_EXECUTABLE} -m pip install --prefix \${SOLMOD4RKT_PYTHON_INSTALL_PREFIX_NATIVE} ${CMAKE_CURRENT_BINARY_DIR}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
")
endif()
