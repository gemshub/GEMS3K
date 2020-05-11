# Recursively collect all header files from the current directory
file(GLOB HEADER_FILES *.h)

# Recursively collect all source files from the current directory
file(GLOB SOURCE_FILES *.cpp)

# Set some necessary definitions
#add_definitions(-DIPMGEMPLUGIN)
add_definitions(-DNOPARTICLEARRAY)
add_definitions(-DNODEARRAYLEVEL)

# Compile the C++ source code into object files
add_library(GEMS3K_OBJECT OBJECT ${HEADER_FILES} ${SOURCE_FILES})

# Check if a shared library is to be built
if(BUILD_SHARED_LIBS)
    add_library(GEMS3K_SHARED SHARED $<TARGET_OBJECTS:GEMS3K_OBJECT>)
#if( ${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
    target_link_libraries(GEMS3K_SHARED zmq )
#endif()
    set_target_properties(GEMS3K_SHARED PROPERTIES OUTPUT_NAME GEMS3K)
    install(TARGETS GEMS3K_SHARED DESTINATION "lib" COMPONENT libraries)
endif()

# Check if a static library is to be built
if(BUILD_STATIC_LIBS)
   add_library(GEMS3K_STATIC STATIC $<TARGET_OBJECTS:GEMS3K_OBJECT>)
   set_target_properties(GEMS3K_STATIC PROPERTIES OUTPUT_NAME GEMS3K-static)
   install(TARGETS GEMS3K_STATIC DESTINATION "lib" COMPONENT libraries)
endif()

# Create an install target for the header files
# Install the header files preserving the directory hierarchy
install(FILES ${HEADER_FILES}
        DESTINATION include/GEMS3K )
#install(DIRECTORY ${GEMS3K_SOURCE_DIR}/GEMS3K
#    DESTINATION "include" COMPONENT headers
#    FILES_MATCHING PATTERN "*.h")
