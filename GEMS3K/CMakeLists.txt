# Recursively collect all header files from the current directory
file(GLOB HEADER_FILES *.h)

# Recursively collect all source files from the current directory
file(GLOB SOURCE_FILES *.cpp)

# Recursively collect all source files from the current directory
#file(GLOB SIMDJSON_FILES ${CMAKE_SOURCE_DIR}/simdjson/simdjson.cpp)

# Compile the C++ source code into object files
add_library(GEMS3K_OBJECT OBJECT ${HEADER_FILES} ${SOURCE_FILES}  )

# Check if a shared library is to be built
if(BUILD_SHARED_LIBS)
    # Enable automatic creation of a module definition (.def) file for a SHARED library on Windows.
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
    add_library(GEMS3K_SHARED SHARED $<TARGET_OBJECTS:GEMS3K_OBJECT>)
    set_target_properties(GEMS3K_SHARED PROPERTIES OUTPUT_NAME GEMS3K)
    if(USE_THERMOFUN MATCHES ON)
      target_link_libraries(GEMS3K_SHARED PUBLIC ThermoFun::ThermoFun)
    endif()
    install( TARGETS GEMS3K_SHARED
             EXPORT GEMS3KTargets DESTINATION "lib" COMPONENT libraries)
endif()


# Check if a static library is to be built
if(BUILD_STATIC_LIBS)
   add_library(GEMS3K_STATIC STATIC $<TARGET_OBJECTS:GEMS3K_OBJECT>)
   set_target_properties(GEMS3K_STATIC PROPERTIES OUTPUT_NAME GEMS3K-static)
   if(USE_THERMOFUN MATCHES ON)
     target_link_libraries(GEMS3K_STATIC PUBLIC ThermoFun::ThermoFun)
   endif()
   install( TARGETS GEMS3K_STATIC
           EXPORT GEMS3KTargets DESTINATION "lib" COMPONENT libraries)
endif()


# Install debug symbols
if(MSVC)
    install(
        FILES $<TARGET_PDB_FILE:GEMS3K_SHARED>
        DESTINATION ${CMAKE_INSTALL_BINDIR}
        COMPONENT libraries
        OPTIONAL
    )
endif()

# Install GEMS3K header files
install(DIRECTORY ${PROJECT_SOURCE_DIR}/GEMS3K
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} COMPONENT headers
    FILES_MATCHING PATTERN "*.hpp" PATTERN "*.h")



