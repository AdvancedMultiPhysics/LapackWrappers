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
SET( CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH} )
INCLUDE_DIRECTORIES( "${CMAKE_INSTALL_PREFIX}/include" )


# Find blas/lapack
IF ( TPL_DIRECTORY )
    # Use the TPL package to specify blas/lapack
    FIND_PACKAGE( TPLs REQUIRED OPTIONAL_COMPONENTS LAPACK PATHS "${TPL_DIRECTORY}" )
    IF ( NOT DEFINED LAPACK_FOUND )
        SET( DISABLE_BLAS TRUE )
        SET( DISABLE_LAPACK TRUE )
    ELSEIF ( LAPACK_FOUND AND NOT ( DISABLE_BLAS OR DISABLE_LAPACK ) )
        SET( CMAKE_INSTALL_RPATH "${BLAS_DIR}" "${LAPACK_DIR}" ${CMAKE_INSTALL_RPATH} )
        INCLUDE( cmake/WriteLapackHeader.cmake )
    ELSE()
        SET( DISABLE_BLAS TRUE )
        SET( DISABLE_LAPACK TRUE )
    ENDIF()
ELSE()
    # Standalone configure
    ENABLE_LANGUAGE( CXX )
    INCLUDE( cmake/macros.cmake )
    INCLUDE( cmake/WriteRepoVersion.cmake )
    SET_COMPILER_FLAGS()
    IF ( NOT ( DISABLE_BLAS OR DISABLE_LAPACK ) )
        FIND_PACKAGE( BlasLapack )
    ENDIF()
ENDIF()


# Set shared/static library
IF ( NOT DEFINED LIB_TYPE )
    IF ( ENABLE_SHARED )
        SET( LIB_TYPE SHARED )
    ELSE()
        SET( LIB_TYPE STATIC )
    ENDIF()
ENDIF()


# Set system specific libraries
# This should move inside the finds for specific lapack libraries
IF ( ${CMAKE_SYSTEM_NAME} STREQUAL "Windows" )
    # Windows specific system libraries
ELSEIF( ${CMAKE_SYSTEM_NAME} STREQUAL "Linux" )
    # Linux specific system libraries
    SET( SYSTEM_LIBS -lpthread -ldl )
ELSEIF( ${CMAKE_SYSTEM_NAME} STREQUAL "Darwin" )
    # Max specific system libraries
    SET( SYSTEM_LIBS -lpthread -ldl )
ELSEIF( ${CMAKE_SYSTEM_NAME} STREQUAL "Generic" )
    # Generic system libraries
ELSE()
    MESSAGE( FATAL_ERROR "OS not detected" )
ENDIF()


# Update rpath
SET( CMAKE_INSTALL_RPATH "${${PROJ}_INSTALL_DIR}/lib" ${CMAKE_INSTALL_RPATH} )


# Enable tests
INCLUDE( CTest )


# Add distclean target
ADD_DISTCLEAN( libLapackWrappers.a libLapackWrappers.so blas_lapack.h LapackWrappers-*.*.*.tar.gz )


# Check if we want to enable coverage
IF ( ENABLE_GCOV OR ENABLE_COVERAGE )
    FIND_PACKAGE( Coverage )
ENDIF()
IF ( COVERAGE_FOUND )
    SET( ENABLE_COVERAGE true )
    SET( COVERAGE_FLAGS ${COVERAGE_FLAGS} )
    SET( COVERAGE_LIBS  ${COVERAGE_LIBS} )
ENDIF()


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


# Add library
ADD_LIBRARY( LapackWrappers ${LIB_TYPE} ${LapackWrappersSource} )
TARGET_LINK_LIBRARIES( LapackWrappers ${BLAS_LAPACK_LINK} ${COVERAGE_LIBS} ${SYSTEM_LIBS} )
INSTALL( TARGETS LapackWrappers EXPORT LapackWrappersTargets DESTINATION "${${PROJ}_INSTALL_DIR}/lib" )
INSTALL( EXPORT LapackWrappersTargets
         FILE LapackWrappersTargets.cmake
         NAMESPACE LapackWrappers::
         DESTINATION ${${PROJ}_INSTALL_DIR}/lib/cmake/LapackWrappers
)
TARGET_COMPILE_DEFINITIONS( LapackWrappers PUBLIC ${COVERAGE_FLAGS} )


# Generating a Package Configuration File
INCLUDE( CMakePackageConfigHelpers )
SET( INCLUDE_INSTALL_DIR  "${CMAKE_INSTALL_PREFIX}/include" CACHE PATH "Location of header files" )
SET( SYSCONFIG_INSTALL_DIR ${${PROJ}_INSTALL_DIR}/lib/cmake/LapackWrappers CACHE PATH "Location of configuration files" )
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
         DESTINATION ${${PROJ}_INSTALL_DIR}/lib/cmake/LapackWrappers )


# Add the tests
IF (BUILD_TESTING)
    ADD_SUBDIRECTORY( test )
ENDIF()
