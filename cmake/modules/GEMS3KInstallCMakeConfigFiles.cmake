# The path where cmake config files are installed
set(GEMS3K_INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/GEMS3K)

install(EXPORT GEMS3KTargets
    FILE GEMS3KTargets.cmake
    NAMESPACE GEMS3K::
    DESTINATION ${GEMS3K_INSTALL_CONFIGDIR}
    COMPONENT cmake)

include(CMakePackageConfigHelpers)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/GEMS3KConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion)

configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules/GEMS3KConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/GEMS3KConfig.cmake
    INSTALL_DESTINATION ${GEMS3K_INSTALL_CONFIGDIR}
    PATH_VARS GEMS3K_INSTALL_CONFIGDIR)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/GEMS3KConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/GEMS3KConfigVersion.cmake
    DESTINATION ${GEMS3K_INSTALL_CONFIGDIR})
