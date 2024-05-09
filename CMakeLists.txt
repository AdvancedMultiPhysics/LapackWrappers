# Set some CMake properties
CMAKE_MINIMUM_REQUIRED(VERSION 3.1)


# Print a banner
MESSAGE("================")
MESSAGE("Lapack Wrappers ")
MESSAGE("================")


# Set the project name
SET( PROJ LapackWrappers )


# Initialize the project (should disable langagues that the TPL builder will enable)
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


# Load the TPL package
IF ( NOT TPL_DIRECTORY )
    MESSAGE( FATAL_ERROR "TPL_DIRECTORY must be specified")
ENDIF()
SET( CMAKE_MODULE_PATH "${TPL_DIRECTORY}" ${CMAKE_MODULE_PATH} )
FIND_PACKAGE( TPLs REQUIRED COMPONENTS LAPACK )
SET( CMAKE_INSTALL_RPATH "${${PROJ}_INSTALL_DIR}/lib" "${BLAS_DIR}" "${LAPACK_DIR}" ${CMAKE_INSTALL_RPATH} )


# Include the TPL inclde paths and libraries
INCLUDE_DIRECTORIES( ${TPL_INCLUDE_DIRS} )
SET( TPL_LIBS ${TPL_LIBRARIES} )


# Enable tests
ENABLE_TESTING()
INCLUDE( CTest )


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
TARGET_LINK_LIBRARIES( LapackWrappers ${TPL_LIBRARIES} )
INSTALL( TARGETS LapackWrappers DESTINATION "${${PROJ}_INSTALL_DIR}/lib" )



