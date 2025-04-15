#include "LapackWrappers.h"

#include <chrono>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <thread>


using Complex = std::complex<double>;


// Get the time difference in us
static inline int64_t diff( std::chrono::time_point<std::chrono::system_clock> t1,
    std::chrono::time_point<std::chrono::system_clock> t2 )
{
    return static_cast<int64_t>( 1e6 * std::chrono::duration<double>( t2 - t1 ).count() );
}


// Call Lapack::run_test
template<typename TYPE>
void run_test( const std::string &routine, int N, int &us, double &error, int &err )
{
    auto t1   = std::chrono::system_clock::now();
    bool pass = Lapack<TYPE>::run_test( routine, N, error );
    auto t2   = std::chrono::system_clock::now();
    us        = static_cast<int>( diff( t1, t2 ) / N );
    err       = pass ? 0 : 1;
}


// Print test result
void printTest(
    bool pass, bool print, const std::string &name, int us, double error, int &N_errors )
{
    if ( pass && print ) {
        printf( "%7s:  passed:  %5i us  (%e)\n", name.c_str(), us, error );
    } else if ( !pass ) {
        printf( "%7s:  failed:  %5i us  (%e)\n", name.c_str(), us, error );
        N_errors++;
    }
}


// Run basic tests
template<class TYPE>
char getPrefix();
template<class TYPE>
const char *getTypeString();
template<>
char getPrefix<float>()
{
    return 's';
}
template<>
char getPrefix<double>()
{
    return 'd';
}
template<>
char getPrefix<Complex>()
{
    return 'z';
}
template<>
const char *getTypeString<float>()
{
    return "single";
}
template<>
const char *getTypeString<double>()
{
    return "double";
}
template<>
const char *getTypeString<Complex>()
{
    return "complex";
}
template<class TYPE>
int runBasic( bool print_all )
{
    if ( print_all ) {
        printf( "\nMachine parameters (%s)\n", getTypeString<TYPE>() );
        printf( "%s\n", Lapack<TYPE>::machineParams().print().data() );
        printf( "Running %s precision basic tests\n", getTypeString<TYPE>() );
        auto tests = Lapack<TYPE>::list_all_tests();
        if ( tests.empty() ) {
            printf( "   No tests detected\n" );
            return 0;
        }
        printf( "   %s", tests[0].data() );
        for ( size_t i = 1; i < tests.size(); i++ )
            printf( ", %s", tests[i].data() );
        printf( "\n" );
    }
    int N_errors = Lapack<TYPE>::run_all_test();
    if ( N_errors == 0 && print_all ) {
        printf( "  passed %cp tests\n", getPrefix<TYPE>() );
    } else if ( N_errors != 0 ) {
        printf( "  failed %i %cp tests\n", N_errors, getPrefix<TYPE>() );
    }
    printf( "\n" );
    return N_errors;
}


// Run performance tests
template<class TYPE>
int runPerformance( bool print_all, int N_test )
{
    if ( print_all )
        printf( "\nGetting %s precision test times\n", getTypeString<TYPE>() );
    auto tests = Lapack<TYPE>::list_all_tests();
    std::vector<double> error( tests.size() );
    std::vector<int> time( tests.size() );
    std::vector<int> err( tests.size() );
    int N_errors = 0;
    for ( size_t i = 0; i < tests.size(); i++ ) {
        run_test<TYPE>( tests[i], N_test, time[i], error[i], err[i] );
        printTest( err[i] == 0, print_all, tests[i], time[i], error[i], N_errors );
    }
    if ( N_errors == 0 && print_all ) {
        printf( "  passed %cp timing tests\n", getPrefix<TYPE>() );
    } else if ( N_errors != 0 ) {
        printf( "  failed %d %cp timing tests\n", N_errors, getPrefix<TYPE>() );
    }
    return N_errors;
}


// The main function to run all tests
int runAll( bool print_all )
{
    // Number of times to run the tests for timing results
    // Note: the tests are designed to take ~ the same time/test
    const int N_test = 20;

    // Store the number of errors
    int N_errors = 0;

    // Run the basic tests
    N_errors += runBasic<double>( print_all );
    N_errors += runBasic<float>( print_all );
    N_errors += runBasic<Complex>( print_all );

    // Get the times for the tests
    N_errors += runPerformance<double>( print_all, N_test );
    N_errors += runPerformance<float>( print_all, N_test );

    // Run the tests in parallel to check for parallel bugs
    {
        if ( print_all )
            printf( "\nRunning parallel tests\n" );
        int N_threads = 8;
        auto tests    = Lapack<double>::list_all_tests();
        for ( auto &test : tests ) {
            if ( test == "drand" )
                continue; // random number generation is blocking (not parallel)
            auto t1 = std::chrono::system_clock::now();
            std::thread threads[128];
            int time2[128];
            int err2[128];
            double error2[128];
            for ( int j = 0; j < N_threads; j++ )
                threads[j] = std::thread( run_test<double>, test, N_test, std::ref( time2[j] ),
                    std::ref( error2[j] ), std::ref( err2[j] ) );
            for ( int j = 0; j < N_threads; j++ )
                threads[j].join();
            auto t2      = std::chrono::system_clock::now();
            bool pass    = true;
            double error = 0;
            for ( int j = 0; j < N_threads; j++ ) {
                pass  = pass && err2[j] == 0;
                error = std::max( error, error2[j] );
            }
            auto us = static_cast<int>( diff( t1, t2 ) / ( N_test * N_threads ) );
            printTest( pass, print_all, test, us, error, N_errors );
        }
    }

    return N_errors;
}


// main
int main( int argc, char *argv[] )
{
    // Read inputs
    int N      = 1;
    bool print = true;
    if ( argc == 2 ) {
        print = false;
        N     = std::stoi( argv[1] );
    } else if ( argc > 2 ) {
        std::cerr << "Invald call";
        return 1;
    }

    // Set the number of threads
    Lapack<double>::set_num_threads( 1 );

    // Print the lapack version
    if ( print )
        std::cout << Lapack<double>::info();

    // Run the tests
    int N_errors = 0;
    for ( int i = 0; i < N; i++ )
        N_errors += runAll( print );

    // Finished
    if ( N_errors == 0 )
        std::cout << "\nAll tests passed\n";
    else
        std::cout << "\nSome tests failed\n";
    return N_errors == 0 ? 0 : 1;
}
