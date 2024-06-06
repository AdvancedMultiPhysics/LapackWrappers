# Set some CMake properties
CMAKE_MINIMUM_REQUIRED(VERSION 3.21)
CMAKE_POLICY( SET CMP0069 NEW )


# Print a banner
MESSAGE("================")
MESSAGE("Lapack Wrappers ")
MESSAGE("================")


# Set the project name
SET( PROJ LapackWrappers )                  # Set the project name for CMake
SET( LapackWrappers_LIB LapackWrappers )    # Set the final library name


# Initialize the project (should disable languages that the TPL builder will enable)
PROJECT( ${PROJ} LANGUAGES )


# Set the install path
IF( PREFIX )
    SET( CMAKE_INSTALL_PREFIX ${PREFIX} )
ELSEIF( INSTALL_DIR )
    SET( CMAKE_INSTALL_PREFIX ${INSTALL_DIR} )
ELSEIF( TPL_INSTALL_DIR )
    SET( CMAKE_INSTALL_PREFIX ${TPL_INSTALL_DIR} )
ELSE()
    SET( CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_BINARY_DIR} )
ENDIF()
SET( ${PROJ}_INSTALL_DIR ${CMAKE_INSTALL_PREFIX} )
SET( ${PROJ}_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR} )
MESSAGE("Installing LapackWrappers in: ${CMAKE_INSTALL_PREFIX}")


# Set sared/static library
IF ( ENABLE_SHARED )
    SET( LIB_TYPE SHARED )
ELSE()
    SET( LIB_TYPE STATIC )
ENDIF()


# Find blas/lapack
IF ( DISABLE_BLAS OR DISABLE_LAPACK )
    # Configuring without lapack/blas
    ENABLE_LANGUAGE( CXX )
    INCLUDE( cmake/macros.cmake )
    INCLUDE( cmake/WriteRepoVersion.cmake )
ELSEIF ( TPL_DIRECTORY )
    # Use the TPL package to specify blas/lapack
    FIND_PACKAGE( TPLs REQUIRED COMPONENTS LAPACK PATHS "${TPL_DIRECTORY}" )
    SET( CMAKE_INSTALL_RPATH "${BLAS_DIR}" "${LAPACK_DIR}" ${CMAKE_INSTALL_RPATH} )
    INCLUDE( cmake/WriteLapackHeader.cmake )
ELSE()
    # Find Blas/Lapack
    ENABLE_LANGUAGE( CXX )
    SET( CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH} )
    FIND_PACKAGE( BlasLapack )
    INCLUDE( cmake/macros.cmake )
    INCLUDE( cmake/WriteRepoVersion.cmake )
    INCLUDE( cmake/WriteLapackHeader.cmake )
ENDIF()
INCLUDE_DIRECTORIES( "${CMAKE_INSTALL_PREFIX}/include" )


# Update rpath
SET( CMAKE_INSTALL_RPATH "${${PROJ}_INSTALL_DIR}/lib" ${CMAKE_INSTALL_RPATH} )


# Enable tests
ENABLE_TESTING()
INCLUDE( CTest )


# Create custom targets for distclean
ADD_DISTCLEAN( libLapackWrappers.a libLapackWrappers.so blas_lapack.h LapackWrappers-*.*.*.tar.gz )


# Create a release target
WRITE_REPO_VERSION()
CREATE_RELEASE()


# Install the library
FILE( COPY "${${PROJ}_SOURCE_DIR}/LapackWrappers.h" DESTINATION "${${PROJ}_INSTALL_DIR}/include" )
CHECK_ENABLE_FLAG( DISABLE_BLAS 0 )
CHECK_ENABLE_FLAG( DISABLE_LAPACK 0 )
IF ( DISABLE_BLAS OR DISABLE_LAPACK )
    SET( LapackWrappersSource LapackWrappers.cpp LapackWrappersInfo.cpp LapackWrappers.default.cpp )
    ADD_DEFINITIONS( -DDISABLE_LAPACK )
ELSE()
    SET( LapackWrappersSource LapackWrappers.cpp LapackWrappers.double.cpp LapackWrappers.float.cpp LapackWrappers.complex.cpp LapackWrappersInfo.cpp LapackWrappersMachineParams.cpp )
ENDIF()
ADD_LIBRARY( LapackWrappers ${LIB_TYPE} ${LapackWrappersSource} )
TARGET_LINK_LIBRARIES( LapackWrappers ${BLAS_LAPACK_LINK} )
INSTALL( TARGETS LapackWrappers EXPORT LapackWrappersTargets DESTINATION "${${PROJ}_INSTALL_DIR}/lib" )
INSTALL( EXPORT LapackWrappersTargets
         FILE LapackWrappersTargets.cmake
         NAMESPACE LapackWrappers::
         DESTINATION ${${PROJ}_INSTALL_DIR}/lib/cmake
)


# Generating a Package Configuration File
INCLUDE( CMakePackageConfigHelpers )
SET( INCLUDE_INSTALL_DIR  "${CMAKE_INSTALL_PREFIX}/include" CACHE PATH "Location of header files" )
SET( SYSCONFIG_INSTALL_DIR ${${PROJ}_INSTALL_DIR}/lib/cmake CACHE PATH "Location of configuration files" )
CONFIGURE_PACKAGE_CONFIG_FILE( cmake/LapackWrappers.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/cmake/LapackWrappersConfig.cmake
    INSTALL_DESTINATION "${${PROJ}_INSTALL_DIR}/lib/cmake/LapackWrappers"
    PATH_VARS INCLUDE_INSTALL_DIR SYSCONFIG_INSTALL_DIR )
WRITE_BASIC_PACKAGE_VERSION_FILE(
    ${CMAKE_CURRENT_BINARY_DIR}/cmake/LapackWrappersVersion.cmake
    VERSION ${${PROJ}_VERSION}
    COMPATIBILITY SameMajorVersion )
INSTALL( FILES ${CMAKE_CURRENT_BINARY_DIR}/cmake/LapackWrappersConfig.cmake
               ${CMAKE_CURRENT_BINARY_DIR}/cmake/LapackWrappersVersion.cmake
         DESTINATION ${${PROJ}_INSTALL_DIR}/lib/cmake )

# Add the tests
ADD_SUBDIRECTORY( test )

