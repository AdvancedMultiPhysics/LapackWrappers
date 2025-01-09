#include "LapackWrappers.h"

#ifndef DISABLE_LAPACK
#include "blas_lapack.h"
#endif

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdio>
#include <functional>
#include <iostream>
#include <limits>
#include <mutex>
#include <random>
#include <string.h>
#include <string>
#include <thread>


// Choose the OS
#if defined( WIN32 ) || defined( _WIN32 ) || defined( WIN64 ) || defined( _WIN64 )
#define WINDOWS
#include <process.h>
#include <stdlib.h>
#include <windows.h>
#else
#define LINUX
#include <pthread.h>
#include <unistd.h>
#endif


// Function to replace all instances of a string with another
static inline std::string strrep(
    const std::string &in, const std::string &s, const std::string &r )
{
    std::string str( in );
    size_t i = 0;
    while ( i < str.length() ) {
        i = str.find( s, i );
        if ( i == std::string::npos ) {
            break;
        }
        str.replace( i, s.length(), r );
        i += r.length();
    }
    return str;
}


// Function to run a command and capture stdout
[[maybe_unused]] static inline std::string runCommand(
    std::function<void( void )> fun, const std::string &prefix )
{
    fflush( stdout ); // clean everything first
    char buffer[2048];
    memset( buffer, 0, sizeof( buffer ) );
    auto out                  = dup( STDOUT_FILENO );
    [[maybe_unused]] auto tmp = freopen( "NUL", "a", stdout );
    setvbuf( stdout, buffer, _IOFBF, 2048 );
    fun();
    tmp = freopen( "NUL", "a", stdout );
    dup2( out, STDOUT_FILENO );
    setvbuf( stdout, NULL, _IONBF, 2048 );
    std::string str = prefix + strrep( buffer, "\n", "\n  " );
    while ( !str.empty() ) {
        char tmp = str.back();
        if ( tmp > 32 && tmp != ' ' )
            break;
        str.pop_back();
    }
    str += "\n";
    return str;
}


/****************************************************************************
 *  Function to set an environemental variable                               *
 ****************************************************************************/
#if defined( USE_MKL ) || defined( USE_MATLAB_LAPACK )
static void setenv( const char *name, const char *value )
{
    static std::mutex lock;
    lock.lock();
#if defined( WIN32 ) || defined( _WIN32 ) || defined( WIN64 ) || defined( _WIN64 ) || \
    defined( _MSC_VER )
    SetEnvironmentVariable( name, value ) != 0;
#else
    if ( value != nullptr )
        ::setenv( name, value, 1 );
    else
        ::unsetenv( name );
#endif
    lock.unlock();
}
#endif


/******************************************************************
 * Set the number of threads to use                                *
 ******************************************************************/
static int setThreads( [[maybe_unused]] int N )
{
    int N2 = 0;
#if defined( USE_MKL )
    setenv( "MKL_NUM_THREADS", std::to_string( N ).c_str() );
    N2 = N;
#elif defined( USE_OPENBLAS )
    openblas_set_num_threads( N );
    N2 = openblas_get_num_threads();
#elif defined( USE_MATLAB_LAPACK )
    setenv( "OMP_NUM_THREADS", std::to_string( N ).c_str() );
    setenv( "MKL_NUM_THREADS", std::to_string( N ).c_str() );
    N2 = N;
#endif
    return N2;
}
template<>
int Lapack<float>::set_num_threads( int N )
{
    return setThreads( N );
}
template<>
int Lapack<double>::set_num_threads( int N )
{
    return setThreads( N );
}
static bool disable_threads()
{
    setThreads( 1 );
    return true;
}
bool global_lapack_threads_disabled = disable_threads();


/******************************************************************
 * Set the vendor string                                           *
 ******************************************************************/
// clang-format off
#ifdef DISABLE_LAPACK
    static constexpr char LapackVendor[] = "Built-in lapack routines";
#elif defined( USE_ATLAS )
    static constexpr char LapackVendor[] = "ATLAS";
#elif defined( USE_ACML )
    static constexpr char LapackVendor[] = "ACML";
#elif defined( USE_MKL )
    static constexpr char LapackVendor[] = "MKL";
#elif defined( USE_MATLAB_LAPACK )
    static constexpr char LapackVendor[] = "MATLAB LAPACK";
#elif defined( USE_VECLIB )
    static constexpr char LapackVendor[] = "VECLIB";
#elif defined( USE_OPENBLAS )
    static constexpr char LapackVendor[] = "OpenBLAS";
#else
    static constexpr char LapackVendor[] = "Unknown";
#endif
// clang-format on


/******************************************************************
 * Print the lapack information                                    *
 ******************************************************************/
template<typename TYPE>
std::string Lapack<TYPE>::info()
{
    // Print the vendor info
    std::string msg( LapackVendor );
    msg += "\n";
    // Use vendor-specific utility functions
#ifdef USE_OPENBLAS
    msg += "  " + std::string( openblas_get_config() );
#endif
    // Get vendor specific output (capture stdout)
#ifdef USE_ACML
    msg += runCommand( acmlinfo() );
#endif
    return msg;
}
template std::string Lapack<double>::info();
template std::string Lapack<float>::info();
