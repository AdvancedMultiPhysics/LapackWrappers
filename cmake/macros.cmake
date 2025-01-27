# Append a list to a file
FUNCTION( APPEND_LIST FILENAME VARS PREFIX POSTFIX )
    FOREACH ( tmp ${VARS} )
        FILE( APPEND "${FILENAME}" "${PREFIX}" )
        FILE( APPEND "${FILENAME}" "${tmp}" )
        FILE( APPEND "${FILENAME}" "${POSTFIX}" )
    ENDFOREACH ()
ENDFUNCTION()


# Check if a flag is enabled
MACRO( CHECK_ENABLE_FLAG FLAG DEFAULT )
    IF ( NOT DEFINED ${FLAG} )
        SET( ${FLAG} ${DEFAULT} )
    ELSEIF ( ${FLAG} STREQUAL "" )
        SET( ${FLAG} ${DEFAULT} )
    ELSEIF ( ( ${${FLAG}} STREQUAL "FALSE" ) OR ( ${${FLAG}} STREQUAL "false" ) OR ( ${${FLAG}} STREQUAL "0" ) OR ( ${${FLAG}} STREQUAL "OFF" ) )
        SET( ${FLAG} 0 )
    ELSEIF ( ( ${${FLAG}} STREQUAL "TRUE" ) OR ( ${${FLAG}} STREQUAL "true" ) OR ( ${${FLAG}} STREQUAL "1" ) OR ( ${${FLAG}} STREQUAL "ON" ) )
        SET( ${FLAG} 1 )
    ELSE()
        MESSAGE( "Bad value for ${FLAG} (${${FLAG}}); use true or false" )
    ENDIF()
ENDMACRO()


# Set compiler flags
MACRO( SET_COMPILER_FLAGS )
    IF ( CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX OR (${CMAKE_CXX_COMPILER_ID} MATCHES "GNU") )
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -Woverloaded-virtual -Wsign-compare -Wformat-security -Wformat-overflow=2 -Wno-aggressive-loop-optimizations")
    ELSEIF ( MSVC OR MSVC_IDE OR MSVC60 OR MSVC70 OR MSVC71 OR MSVC80 OR CMAKE_COMPILER_2005 OR MSVC90 OR MSVC10 )
        # Add Microsoft specifc compiler options
        SET( CMAKE_C_FLAGS     "${CMAKE_C_FLAGS} /D _SCL_SECURE_NO_WARNINGS /D _CRT_SECURE_NO_WARNINGS /D _ITERATOR_DEBUG_LEVEL=0 /wd4267 /Zc:preprocessor" )
        SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /D _SCL_SECURE_NO_WARNINGS /D _CRT_SECURE_NO_WARNINGS /D _ITERATOR_DEBUG_LEVEL=0 /wd4267 /Zc:preprocessor" )
    ELSEIF ( ${CMAKE_CXX_COMPILER_ID} MATCHES "Intel" ) 
        SET(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} -Wall" )
    ELSEIF ( ${CMAKE_CXX_COMPILER_ID} MATCHES "PGI" )
        SET(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} -Minform=inform -Mlist --display_error_number --diag_suppress 111,128,185")
    ELSEIF ( (${CMAKE_CXX_COMPILER_ID} MATCHES "CRAY") OR (${CMAKE_CXX_COMPILER_ID} MATCHES "Cray") )
        SET(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS}")
    ELSEIF ( (${CMAKE_CXX_COMPILER_ID} MATCHES "CLANG") OR (${CMAKE_CXX_COMPILER_ID} MATCHES "Clang") )
        SET(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic -Wno-missing-braces -Wmissing-field-initializers -ftemplate-depth=1024")
    ELSEIF ( ${CMAKE_CXX_COMPILER_ID} MATCHES "XL" )
        SET(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} -Wall -ftemplate-depth=512")
    ELSEIF ( (${CMAKE_C_COMPILER_ID} MATCHES "NVHPC") OR (${CMAKE_CXX_COMPILER_ID} MATCHES "NVHPC") )
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -Woverloaded-virtual -Wsign-compare -Wformat-security")
    ELSE()
        MESSAGE( "CMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}")
        MESSAGE( "CMAKE_CXX_COMPILER_ID=${CMAKE_CXX_COMPILER_ID}")
        MESSAGE(FATAL_ERROR "Unknown C/C++ compiler")
    ENDIF()
    IF ( "${CMAKE_BUILD_TYPE}" STREQUAL "Debug" AND NOT ("${CMAKE_CXX_FLAGS_DEBUG}" MATCHES "-D_DEBUG") )
        SET( CMAKE_C_FLAGS_DEBUG   " ${CMAKE_C_FLAGS_DEBUG}   -DDEBUG -D_DEBUG" )
        SET( CMAKE_CXX_FLAGS_DEBUG " ${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG -D_DEBUG" )        
    ENDIF()
    # Enable GLIBCXX_DEBUG flags
    CHECK_ENABLE_FLAG( ENABLE_GXX_DEBUG 0 )
    IF ( ENABLE_GXX_DEBUG AND NOT ("${CMAKE_CXX_FLAGS_DEBUG}" MATCHES "-D_GLIBCXX_DEBUG") ) 
        SET( CMAKE_CXX_FLAGS_DEBUG " ${CMAKE_CXX_FLAGS_DEBUG} -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC" )
    ENDIF()
    MESSAGE( "CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")
ENDMACRO()


# add custom target distclean
# cleans and removes cmake generated files etc.
FUNCTION( ADD_DISTCLEAN ${ARGN} )
    SET(DISTCLEANED
        assembly
        cmake.depends
        cmake.check_depends
        CMakeCache.txt
        CMakeFiles
        CMakeTmp
        CMakeDoxy*
        cmake.check_cache
        *.cmake
        compile.log
        cppcheck-build
        Doxyfile
        Makefile
        core core.*
        DartConfiguration.tcl
        install_manifest.txt
        Testing
        include
        doc
        docs
        examples
        latex_docs
        lib
        Makefile.config
        install_manifest.txt
        test
        matlab
        Matlab
        mex
        tmp
        #tmp#
        bin
        cmake
        cppclean
        compile_commands.json
        TPLs.h
        ${ARGN}
    )
    ADD_CUSTOM_TARGET(distclean @echo cleaning for source distribution)
    IF (UNIX)
        ADD_CUSTOM_COMMAND(
            DEPENDS clean
            COMMENT "distribution clean"
            COMMAND rm
            ARGS    -Rf ${DISTCLEANED}
            TARGET  distclean
        )
    ELSE()
        SET( DISTCLEANED
            ${DISTCLEANED}
            *.vcxproj*
            ipch
            x64
            Debug
        )
        SET( DISTCLEAN_FILE "${CMAKE_CURRENT_BINARY_DIR}/distclean.bat" )
        FILE( WRITE  "${DISTCLEAN_FILE}" "del /s /q /f " )
        APPEND_LIST( "${DISTCLEAN_FILE}" "${DISTCLEANED}" " " " " )
        FILE( APPEND "${DISTCLEAN_FILE}" "\n" )
        APPEND_LIST( "${DISTCLEAN_FILE}" "${DISTCLEANED}" "for /d %%x in ("   ") do rd /s /q \"%%x\"\n" )
        ADD_CUSTOM_COMMAND(
            DEPENDS clean
            COMMENT "distribution clean"
            COMMAND distclean.bat & del /s/q/f distclean.bat
            TARGET  distclean
        )
    ENDIF()
ENDFUNCTION()

