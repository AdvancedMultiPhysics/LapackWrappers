# ctest script for building, running, and submitting the test results 
# Usage:  ctest -S script,build
#   build = debug / optimized / weekly / valgrind
# Note: this test will use use the number of processors defined in the variable N_PROCS,
#   the environmental variable N_PROCS, or the number of processors available (if not specified)


# Set the Project variables
SET( PROJ LapackWrappers )


# Set platform specific variables
SITE_NAME( HOSTNAME )
STRING( REGEX REPLACE  "-(ext|login)(..|.)"  ""  HOSTNAME  "${HOSTNAME}" )
SET( USE_MPI             $ENV{USE_MPI}             )
SET( CC                  $ENV{CC}                  )
SET( CXX                 $ENV{CXX}                 )
SET( CFLAGS              $ENV{CFLAGS}              ) 
SET( CXXFLAGS            $ENV{CXXFLAGS}            )
SET( CXX_STD             $ENV{CXX_STD}             )
SET( MPIEXEC             $ENV{MPIEXEC}             )
SET( TIMER_DIRECTORY    "$ENV{TIMER_DIRECTORY}"    )
SET( COVERAGE_COMMAND    $ENV{COVERAGE_COMMAND}    )
SET( VALGRIND_COMMAND    $ENV{VALGRIND_COMMAND}    )
SET( CMAKE_MAKE_PROGRAM "$ENV{CMAKE_MAKE_PROGRAM}" )
SET( CTEST_CMAKE_GENERATOR $ENV{CTEST_CMAKE_GENERATOR} )
SET( LDLIBS              $ENV{LDLIBS}              )
SET( LDFLAGS             $ENV{LDFLAGS}             )
SET( MPI_DIRECTORY       $ENV{MPI_DIRECTORY}       )
SET( MPI_INCLUDE         $ENV{MPI_INCLUDE}         )
SET( MPI_LINK_FLAGS     "$ENV{MPI_LINK_FLAGS}"     )
SET( MPI_LIBRARIES       $ENV{MPI_LIBRARIES}       )
SET( MPIEXEC             $ENV{MPIEXEC}             )
SET( BUILD_SERIAL        $ENV{BUILD_SERIAL}        )
SET( SKIP_TESTS          $ENV{SKIP_TESTS}          )
SET( BUILDNAME_POSTFIX  "$ENV{BUILDNAME_POSTFIX}"  )
SET( CTEST_URL          "$ENV{CTEST_URL}"          )
SET( TPL_DIRECTORY      "$ENV{TPL_DIRECTORY}"      )
SET( MKL_INCLUDE_DIR    "$ENV{MKL_INCLUDE_DIR}"    )
SET( MKL_LIB_DIR        "$ENV{MKL_LIB_DIR}"        )
SET( LAPACK_INSTALL_DIR "$ENV{LAPACK_INSTALL_DIR}" )

# Get the source directory based on the current directory
IF ( NOT ${PROJ}_SOURCE_DIR )
    SET( ${PROJ}_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/.." )
ENDIF()
IF ( NOT CMAKE_MAKE_PROGRAM )
    SET( CMAKE_MAKE_PROGRAM make )
ENDIF()


# Check that we specified the build type to run
SET( RUN_WEEKLY FALSE )
IF( NOT CTEST_SCRIPT_ARG )
    MESSAGE(FATAL_ERROR "No build specified: ctest -S /path/to/script,build (debug/optimized/valgrind")
ELSEIF( ${CTEST_SCRIPT_ARG} STREQUAL "debug" )
    SET( CTEST_BUILD_NAME "${PROJ}-debug" )
    SET( CMAKE_BUILD_TYPE "Debug" )
    SET( CTEST_COVERAGE_COMMAND ${COVERAGE_COMMAND} )
    SET( ENABLE_GCOV "true" )
    SET( USE_VALGRIND FALSE )
ELSEIF( (${CTEST_SCRIPT_ARG} STREQUAL "optimized") OR (${CTEST_SCRIPT_ARG} STREQUAL "opt") )
    SET( CTEST_BUILD_NAME "${PROJ}-opt" )
    SET( CMAKE_BUILD_TYPE "Release" )
    SET( CTEST_COVERAGE_COMMAND )
    SET( ENABLE_GCOV "false" )
    SET( USE_VALGRIND FALSE )
ELSEIF( (${CTEST_SCRIPT_ARG} STREQUAL "weekly") )
    SET( CTEST_BUILD_NAME "${PROJ}-Weekly" )
    SET( CMAKE_BUILD_TYPE "Release" )
    SET( CTEST_COVERAGE_COMMAND )
    SET( ENABLE_GCOV "false" )
    SET( USE_VALGRIND FALSE )
    SET( RUN_WEEKLY TRUE )
ELSEIF( ${CTEST_SCRIPT_ARG} STREQUAL "valgrind" )
    SET( CTEST_BUILD_NAME "${PROJ}-valgrind" )
    SET( CMAKE_BUILD_TYPE "Debug" )
    SET( CTEST_COVERAGE_COMMAND )
    SET( ENABLE_GCOV "false" )
    SET( USE_VALGRIND TRUE )
ELSEIF( ${CTEST_SCRIPT_ARG} STREQUAL "doc" )
    SET( CTEST_BUILD_NAME "${PROJ}-doc" )
    SET( CMAKE_BUILD_TYPE "Release" )
    SET( BUILD_ONLY_DOCS "true" )
ELSE()
    MESSAGE(FATAL_ERROR "Invalid build (${CTEST_SCRIPT_ARG}): ctest -S /path/to/script,build (debug/opt/valgrind")
ENDIF()
IF ( BUILDNAME_POSTFIX )
    SET( CTEST_BUILD_NAME "${CTEST_BUILD_NAME}-${BUILDNAME_POSTFIX}" )
ENDIF()
IF ( NOT COVERAGE_COMMAND )
    SET( ENABLE_GCOV "false" )
ENDIF()


# Set the number of processors
IF( NOT DEFINED N_PROCS )
    SET( N_PROCS $ENV{N_PROCS} )
ENDIF()
IF ( NOT DEFINED N_PROCS )
    SET(N_PROCS 1)
    # Linux:
    SET(cpuinfo_file "/proc/cpuinfo")
    IF(EXISTS "${cpuinfo_file}")
        FILE(STRINGS "${cpuinfo_file}" procs REGEX "^processor.: [0-9]+$")
        list(LENGTH procs N_PROCS)
    ENDIF()
    # Mac:
    IF(APPLE)
        find_program(cmd_sys_pro "sysctl")
        if(cmd_sys_pro)
            execute_process(COMMAND ${cmd_sys_pro} hw.physicalcpu OUTPUT_VARIABLE info)
            STRING(REGEX REPLACE "^.*hw.physicalcpu: ([0-9]+).*$" "\\1" N_PROCS "${info}")
        ENDIF()
    ENDIF()
    # Windows:
    IF(WIN32)
        SET(N_PROCS "$ENV{NUMBER_OF_PROCESSORS}")
    ENDIF()
ENDIF()


# Set the nightly start time
# This controls the version of a checkout from cvs/svn (ignored for mecurial/git)
# This does not control the start of the day displayed on CDash, that is controled by the CDash project settings
SET( NIGHTLY_START_TIME "$ENV{NIGHTLY_START_TIME}" )
IF ( NOT NIGHTLY_START_TIME )
    SET( NIGHTLY_START_TIME "18:00:00 EST" )
ENDIF()
SET( CTEST_NIGHTLY_START_TIME ${NIGHTLY_START_TIME} )


# Set basic variables
SET( CTEST_PROJECT_NAME "${PROJ}" )
SET( CTEST_SOURCE_DIRECTORY "${${PROJ}_SOURCE_DIR}" )
SET( CTEST_BINARY_DIRECTORY "." )
SET( CTEST_DASHBOARD "Nightly" )
SET( CTEST_CUSTOM_MAXIMUM_NUMBER_OF_ERRORS 500 )
SET( CTEST_CUSTOM_MAXIMUM_NUMBER_OF_WARNINGS 500 )
SET( CTEST_CUSTOM_MAXIMUM_PASSED_TEST_OUTPUT_SIZE 10000 )
SET( CTEST_CUSTOM_MAXIMUM_FAILED_TEST_OUTPUT_SIZE 10000 )
SET( CTEST_COMMAND "\"${CTEST_EXECUTABLE_NAME}\" -D ${CTEST_DASHBOARD}" )
IF ( BUILD_SERIAL )
    SET( CTEST_BUILD_COMMAND "${CMAKE_MAKE_PROGRAM} -i install" )
ELSE()
    SET( CTEST_BUILD_COMMAND "${CMAKE_MAKE_PROGRAM} -i -j ${N_PROCS} install" )
ENDIF()
SET( CTEST_CUSTOM_WARNING_EXCEPTION )
SET( CTEST_CUSTOM_ERROR_EXCEPTION )


# Set timeouts: 5 minutes for debug, 2 for opt, and 30 minutes for valgrind/weekly
IF ( USE_VALGRIND )
    SET( CTEST_TEST_TIMEOUT 1800 )
ELSEIF ( RUN_WEEKLY )
    SET( CTEST_TEST_TIMEOUT 1800 )
ELSEIF( ${CMAKE_BUILD_TYPE} STREQUAL "Debug" )
    SET( CTEST_TEST_TIMEOUT 300 )
ELSE()
    SET( CTEST_TEST_TIMEOUT 120 )
ENDIF()


# Set valgrind options
#SET (VALGRIND_COMMAND_OPTIONS "--tool=memcheck --leak-check=yes --track-fds=yes --num-callers=50 --show-reachable=yes --track-origins=yes --malloc-fill=0xff --free-fill=0xfe --suppressions=${${PROJ}_SOURCE_DIR}/ValgrindSuppresionFile" )
SET( VALGRIND_COMMAND_OPTIONS  "--tool=memcheck --leak-check=yes --track-fds=yes --num-callers=50 --show-reachable=yes --suppressions=${${PROJ}_SOURCE_DIR}/ValgrindSuppresionFile" )
IF ( USE_VALGRIND )
    SET( MEMORYCHECK_COMMAND ${VALGRIND_COMMAND} )
    SET( MEMORYCHECKCOMMAND ${VALGRIND_COMMAND} )
    SET( CTEST_MEMORYCHECK_COMMAND ${VALGRIND_COMMAND} )
    SET( CTEST_MEMORYCHECKCOMMAND ${VALGRIND_COMMAND} )
    SET( CTEST_MEMORYCHECK_COMMAND_OPTIONS ${VALGRIND_COMMAND_OPTIONS} )
    SET( CTEST_MEMORYCHECKCOMMAND_OPTIONS  ${VALGRIND_COMMAND_OPTIONS} )
ENDIF()


# Clear the binary directory and create an initial cache
EXECUTE_PROCESS( COMMAND ${CMAKE_COMMAND} -E remove -f CMakeCache.txt )
EXECUTE_PROCESS( COMMAND ${CMAKE_COMMAND} -E remove_directory CMakeFiles )
FILE(WRITE "${CTEST_BINARY_DIRECTORY}/CMakeCache.txt" "CTEST_TEST_CTEST:BOOL=1")


# Set the configure options
SET( CTEST_OPTIONS "-DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}" )
IF ( TPL_DIRECTORY )
    SET( CTEST_OPTIONS "${CTEST_OPTIONS};-DTPL_DIRECTORY=${TPL_DIRECTORY}")
ELSE()
    SET( CTEST_OPTIONS "${CTEST_OPTIONS};-DENABLE_SHARED:BOOL=TRUE" )
    SET( CTEST_OPTIONS "${CTEST_OPTIONS};-DENABLE_GCOV:BOOL=${ENABLE_GCOV}" )
    IF ( MKL_INCLUDE_DIR )
        SET( CTEST_OPTIONS "${CTEST_OPTIONS};-DLAPACK_INSTALL_DIR=${MKL_INCLUDE_DIR}" )
        SET( CTEST_OPTIONS "${CTEST_OPTIONS};-DMKL_INCLUDE_DIR=${MKL_INCLUDE_DIR}" )
        SET( CTEST_OPTIONS "${CTEST_OPTIONS};-DMKL_LIB_DIR=${MKL_LIB_DIR}" )
    ELSEIF ( LAPACK_INSTALL_DIR )
        SET( CTEST_OPTIONS "${CTEST_OPTIONS};-DLAPACK_INSTALL_DIR=${LAPACK_INSTALL_DIR}" )
    ELSE()
        SET( CTEST_OPTIONS "${CTEST_OPTIONS};-DDISABLE_LAPACK:BOOL=TRUE" )
    ENDIF()
ENDIF()
MESSAGE("Configure options:")
MESSAGE("   ${CTEST_OPTIONS}")


# Configure the drop site
IF ( NOT CTEST_SITE )
    SET( CTEST_SITE ${HOSTNAME} )
ENDIF()
IF ( NOT CTEST_URL )
    SET( CTEST_DROP_METHOD "http" )
    SET( CTEST_DROP_LOCATION "/CDash/submit.php?project=LapackWrappers" )
    SET( CTEST_DROP_SITE_CDASH TRUE )
    SET( DROP_SITE_CDASH TRUE )
    SET( CTEST_DROP_SITE ${CTEST_SITE} )
ELSE()
    STRING( REPLACE "PROJECT" "LapackWrappers" CTEST_URL "${CTEST_URL}" )
    SET( CTEST_SUBMIT_URL "${CTEST_URL}" )
ENDIF()


# Configure and run the tests
CTEST_START( "${CTEST_DASHBOARD}" )
CTEST_UPDATE()
CTEST_SUBMIT( PARTS Update )
CTEST_CONFIGURE(
    BUILD   ${CTEST_BINARY_DIRECTORY}
    SOURCE  ${CTEST_SOURCE_DIRECTORY}
    OPTIONS "${CTEST_OPTIONS}"
)
CTEST_SUBMIT( PARTS Configure )


# Run the configure/build/test
CTEST_BUILD()
CTEST_SUBMIT( PARTS Build )
EXECUTE_PROCESS( COMMAND ${CMAKE_MAKE_PROGRAM} install )
IF ( SKIP_TESTS )
    # Do not run tests
    SET( CTEST_COVERAGE_COMMAND )
ELSEIF ( USE_VALGRIND )
    CTEST_MEMCHECK( EXCLUDE "(procs|WEEKLY|cppcheck|cppclean|test_crash)"  PARALLEL_LEVEL ${N_PROCS} )
ELSEIF ( EXCLUDE_WEEKLY )
    CTEST_TEST( EXCLUDE WEEKLY  PARALLEL_LEVEL ${N_PROCS} )
ELSE()
    CTEST_TEST( PARALLEL_LEVEL ${N_PROCS} )
ENDIF()
IF( CTEST_COVERAGE_COMMAND )
    CTEST_COVERAGE()
ENDIF()
CTEST_SUBMIT( PARTS Test )
CTEST_SUBMIT( PARTS Coverage )
CTEST_SUBMIT( PARTS MemCheck )
CTEST_SUBMIT( PARTS Done )


# Write a message to test for success in the ctest-builder
MESSAGE( "ctest_script ran to completion" )


