#define NOMINMAX
#include "LapackWrappers.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdio>
#include <iostream>
#include <limits>
#include <mutex>
#include <random>
#include <string.h>
#include <string>
#include <thread>
#include <type_traits>


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


// Define some sizes of the problems (try to normalize the time/test)
#define TEST_SIZE_VEC 100000   // Vector tests O(N)
#define TEST_SIZE_MATVEC 1000  // Matrix-vector tests O(N^2)
#define TEST_SIZE_MAT 150      // Matrix-matrix / Dense solves tests O(N^3)
#define TEST_SIZE_TRI 5000     // Tridiagonal/banded tests
#define TEST_SIZE_TRI_MAT 2000 // Tridiagonal/banded solve tests


// Choose precision to perfom calculations
template<class TYPE>
class extended
{
};
template<>
class extended<float>
{
public:
    typedef double TYPE2;
};
template<>
class extended<double>
{
public:
    typedef long double TYPE2;
};
template<>
class extended<std::complex<float>>
{
public:
    typedef std::complex<double> TYPE2;
};
template<>
class extended<std::complex<double>>
{
public:
    typedef std::complex<long double> TYPE2;
};


// Get epsilon
template<class TYPE>
constexpr double epsilon()
{
    if constexpr ( std::is_same_v<TYPE, float> ) {
        return std::numeric_limits<float>::epsilon();
    } else if constexpr ( std::is_same_v<TYPE, double> ) {
        return std::numeric_limits<double>::epsilon();
    } else if constexpr ( std::is_same_v<TYPE, std::complex<double>> ) {
        return 2 * std::numeric_limits<double>::epsilon();
    } else {
        throw std::logic_error( "Not finished" );
    }
}


// Lists the tests availible to run
#ifndef DISABLE_LAPACK
template<>
std::vector<std::string> Lapack<double>::list_all_tests()
{
    return { "dcopy", "dscal", "dnrm2", "daxpy", "dgemv", "dgemm", "dasum", "ddot", "dgesv",
        "dgtsv", "dgbsv", "dgetrf", "dgttrf", "dgbtrf", "dgetrs", "dgttrs", "dgetri", "drand" };
}
template<>
std::vector<std::string> Lapack<float>::list_all_tests()
{
    return { "scopy", "sscal", "snrm2", "saxpy", "sgemv", "sgemm", "sasum", "sdot", "sgesv",
        "sgtsv", "sgbsv", "sgetrf", "sgttrf", "sgbtrf", "sgetrs", "sgttrs", "sgetri", "srand" };
}
template<>
std::vector<std::string> Lapack<std::complex<double>>::list_all_tests()
{
    return {};
}
#else
template<>
std::vector<std::string> Lapack<double>::list_all_tests()
{
    return { "dcopy", "dscal", "dnrm2", "daxpy", "dasum", "ddot", "drand", "dgemv", "dgemm" };
}
template<>
std::vector<std::string> Lapack<float>::list_all_tests()
{
    return { "scopy", "sscal", "snrm2", "saxpy", "sasum", "sdot", "srand", "sgemv", "sgemm" };
}
template<>
std::vector<std::string> Lapack<std::complex<double>>::list_all_tests()
{
    return { "zcopy", "zscal", "znrm2", "zaxpy", "zasum", "zdot", "zrand", "zgemv", "zgemm" };
}
#endif


// Run all the tests
template<typename TYPE>
int run_basic_tests();
template<typename TYPE>
int Lapack<TYPE>::run_all_test()
{
    int N_errors = 0;
    int N        = 2; // We want two iterations to ensure the test works for N>1
    auto tests   = Lapack<TYPE>::list_all_tests();
    for ( auto test : tests ) {
        double error = 0;
        bool pass    = Lapack<TYPE>::run_test( test, N, error );
        if ( !pass ) {
            printf( "test_%s failed (%e)\n", test.c_str(), error );
            N_errors++;
        }
    }
    if ( run_basic_tests<TYPE>() != 0 ) {
        printf( "basic tests failed\n" );
        N_errors++;
    }
    return N_errors;
}
template int Lapack<double>::run_all_test();
template int Lapack<float>::run_all_test();
template int Lapack<std::complex<double>>::run_all_test();


// Fill a vector with random TYPE precision data
template<>
void Lapack<float>::random( int N, float *data )
{
    static std::mutex lock;
    static std::random_device rd;
    static std::mt19937 gen( rd() );
    static std::uniform_real_distribution<float> dis( 0, 1 );
    lock.lock();
    for ( int i = 0; i < N; i++ )
        data[i] = dis( gen );
    lock.unlock();
}
template<>
void Lapack<double>::random( int N, double *data )
{
    static std::mutex lock;
    static std::random_device rd;
    static std::mt19937_64 gen( rd() );
    static std::uniform_real_distribution<double> dis( 0, 1 );
    lock.lock();
    for ( int i = 0; i < N; i++ )
        data[i] = dis( gen );
    lock.unlock();
}
template<>
void Lapack<std::complex<double>>::random( int N, std::complex<double> *data )
{
    static std::mutex lock;
    static std::random_device rd;
    static std::mt19937_64 gen( rd() );
    static std::uniform_real_distribution<double> dis( 0, 1 );
    lock.lock();
    for ( int i = 0; i < N; i++ )
        data[i] = std::complex<double>( dis( gen ), dis( gen ) );
    lock.unlock();
}


// Check if two vectors are approximately equal
template<typename TYPE>
static inline bool approx_equal( int N, const TYPE *x1, const TYPE *x2, const double tol = 1e-12 )
{
    bool pass = true;
    for ( int i = 0; i < N; i++ )
        pass = pass && std::abs( x1[i] - x2[i] ) <= tol * 0.5 * std::abs( x1[i] + x2[i] );
    return pass;
}


// Return the L2 error
template<typename TYPE>
static inline double L2Error( int N, const TYPE *x1, const TYPE *x2 )
{
    double norm = 0.0;
    for ( int i = 0; i < N; i++ )
        norm += std::norm( x1[i] - x2[i] );
    return sqrt( norm );
}


// Generate a matrix/rhs
template<typename TYPE>
void generateRhs( int N, TYPE *x )
{
    Lapack<TYPE>::random( N, x );
}
template<typename TYPE>
void generateMatrix( int N, TYPE *A )
{
    Lapack<TYPE>::random( N * N, A );
    for ( int i = 0; i < N; i++ ) {
        TYPE d       = A[i + i * N];
        A[i + i * N] = 0;
        A[i + i * N] = d - Lapack<TYPE>::asum( N, &A[i * N], 1 );
    }
}
template<typename TYPE>
void generateMatrixBanded( int N, int KL, int KU, TYPE *A )
{
    for ( int i = 0; i < N * N; i++ )
        A[i] = 0;
    for ( int i = 0; i < N; i++ ) {
        int K1 = std::max( i - KL, 0 );
        int K2 = std::min( i + KU, N - 1 );
        Lapack<TYPE>::random( K2 - K1 + 1, &A[K1 + i * N] );
        TYPE d       = A[i + i * N];
        A[i + i * N] = 0;
        A[i + i * N] = d - Lapack<TYPE>::asum( K2 - K1 + 1, &A[K1 + i * N], 1 );
    }
}
template<typename TYPE>
void generateTriDiag( int N, TYPE *DL, TYPE *D, TYPE *DU )
{
    Lapack<TYPE>::random( N - 1, DL );
    Lapack<TYPE>::random( N - 1, DU );
    Lapack<TYPE>::random( N, D );
    D[0] -= DL[0];
    for ( int i = 1; i < N - 1; i++ )
        D[i] -= DL[i] + DU[i - 1];
    D[N - 1] -= DU[N - 2];
}
template<typename TYPE>
void extractBanded( int N, int KL, int KU, const TYPE *A, TYPE *AB )
{
    const int K2 = 2 * KL + KU + 1;
    for ( int k = 0; k < N; k++ ) {
        for ( int k2 = -KL; k2 <= KU; k2++ ) {
            if ( k + k2 < 0 || k + k2 >= N )
                continue;
            AB[k2 + 2 * KL + k * K2] = A[k + k2 + k * N];
        }
    }
}
template<typename TYPE>
void extractTriDiag( int N, const TYPE *A, TYPE *DL, TYPE *D, TYPE *DU )
{
    D[0] = A[0];
    for ( int k = 1; k < N; k++ ) {
        D[k]      = A[k + k * N];
        DU[k - 1] = A[( k - 1 ) + k * N];
        DL[k - 1] = A[k + ( k - 1 ) * N];
    }
}


// Some basic tests (mostly to improve coverage)
template<typename TYPE>
int run_basic_tests()
{
    int N    = 50;
    int KL   = 2;
    int KU   = 3;
    int K2   = 2 * KL + KU + 1;
    auto x   = new TYPE[N];
    auto A   = new TYPE[N * N];
    auto D   = new TYPE[N];
    auto D2  = new TYPE[N];
    auto DL  = new TYPE[N - 1];
    auto DL2 = new TYPE[N - 1];
    auto DU  = new TYPE[N - 1];
    auto DU2 = new TYPE[N - 1];
    auto AB  = new TYPE[N * K2];
    generateRhs( N, x );
    generateMatrix( N, A );
    generateMatrixBanded( N, KL, KU, A );
    generateTriDiag( N, DL, D, DU );
    extractBanded( N, KL, KU, A, AB );
    extractTriDiag( N, A, DL2, D2, DU2 );
    delete[] x;
    delete[] A;
    delete[] D;
    delete[] DL;
    delete[] DU;
    delete[] D2;
    delete[] DL2;
    delete[] DU2;
    delete[] AB;
    return 0;
}


// Test random
template<typename TYPE>
static bool test_random( int N, double &error )
{
    if constexpr ( std::is_floating_point<TYPE>::value ) {
        constexpr int Nb = 25; // Number of bins
        int K            = TEST_SIZE_VEC;
        TYPE *x          = new TYPE[K];
        int count[Nb]    = { 0 };
        error            = 0;
        for ( int it = 0; it < N; it++ ) {
            Lapack<TYPE>::random( K, x );
            TYPE sum = 0;
            for ( int i = 0; i < K; i++ ) {
                error = ( error == 0 && x[i] < 0.0 ) ? -1 : error;
                error = ( error == 0 && x[i] > 1.0 ) ? -2 : error;
                sum += x[i];
                int j = std::min<int>( floor( x[i] * Nb ), Nb - 1 );
                count[j]++;
            }
            TYPE avg = sum / K;
            if ( error == 0 && ( avg < 0.35 || avg > 0.65 ) )
                error = -3;
        }
        delete[] x;
        if ( error != 0 )
            return false;
        // If we did not encounter an error, check the Chi-Square Distribution
        double E  = ( K * N ) / static_cast<double>( Nb );
        double X2 = 0.0;
        for ( int i = 0; i < Nb; i++ )
            X2 += ( count[i] - E ) * ( count[i] - E ) / E;
        double X2c = 52.6; // Critical value for 0.999 (we will fail 0.1% of the time)
        error      = X2 / X2c;
        return error <= 1.0;
    } else if ( std::is_same_v<TYPE, std::complex<double>> ) {
        // Need to improve complex tests
        return true;
    } else {
        throw std::logic_error( "Not finihsed" );
    }
}


// Test copy
template<typename TYPE>
static bool test_copy( int N, double &error )
{
    const int K = TEST_SIZE_VEC;
    TYPE *x1    = new TYPE[K];
    TYPE *x2    = new TYPE[K];
    Lapack<TYPE>::random( K, x1 );
    int N_errors = 0;
    error        = 0;
    for ( int i = 0; i < N; i++ ) {
        memset( (void *) x2, 0xB6, K * sizeof( TYPE ) );
        Lapack<TYPE>::copy( K, x1, 1, x2, 1 );
        if ( !approx_equal( K, x1, x2 ) )
            N_errors++;
    }
    delete[] x1;
    delete[] x2;
    return N_errors == 0;
}


// Test scal
template<typename TYPE>
static bool test_scal( int N, double &error )
{
    const int K = TEST_SIZE_VEC;
    TYPE *x0    = new TYPE[K];
    TYPE *x1    = new TYPE[K];
    TYPE *x2    = new TYPE[K];
    Lapack<TYPE>::random( K, x0 );
    const TYPE pi = static_cast<TYPE>( 3.141592653589793 );
    for ( int j = 0; j < K; j++ )
        x1[j] = pi * x0[j];
    int N_errors = 0;
    error        = 0;
    double tol   = 10 * epsilon<TYPE>();
    for ( int i = 0; i < N; i++ ) {
        memcpy( x2, x1, K * sizeof( TYPE ) );
        Lapack<TYPE>::scal( K, pi, x0, 1 );
        if ( !approx_equal( K, x1, x2, tol ) )
            N_errors++;
    }
    delete[] x0;
    delete[] x1;
    delete[] x2;
    return N_errors == 0;
}

// Test nrm2
template<typename TYPE>
static bool test_nrm2( int N, double &error )
{
    const int K = TEST_SIZE_VEC;
    TYPE *x     = new TYPE[K];
    Lapack<TYPE>::random( K, x );
    double ans1 = 0.0;
    for ( int j = 0; j < K; j++ )
        ans1 += std::norm( x[j] );
    ans1         = sqrt( ans1 );
    int N_errors = 0;
    error        = 0;
    double tol   = sqrt( K ) * epsilon<TYPE>();
    for ( int i = 0; i < N; i++ ) {
        TYPE ans2  = Lapack<TYPE>::nrm2( K, x, 1 );
        double err = std::abs( ans1 - ans2 ) / sqrt( K );
        error      = std::max( error, err );
        if ( err > tol )
            N_errors++;
    }
    delete[] x;
    return N_errors == 0;
}

// Test asum
template<typename TYPE>
static bool test_asum( int N, double &error )
{
    const int K = 2 * TEST_SIZE_VEC;
    // Create a random set of numbers and the sum (simple method)
    TYPE *x = new TYPE[K];
    Lapack<TYPE>::random( K, x );
    long double sum = 0.0;
    for ( int j = 0; j < K; j++ )
        sum += std::abs( x[j] );
    double ans1 = sum;
    // Check asum
    int N_errors = 0;
    error        = 0;
    double tol   = 100 * epsilon<TYPE>();
    for ( int i = 0; i < N; i++ ) {
        double ans2 = Lapack<TYPE>::asum( K, x, 1 );
        double err  = std::abs( ans1 - ans2 ) / K;
        error       = std::max( error, err );
        if ( err > tol )
            N_errors++;
    }
    delete[] x;
    return N_errors == 0;
}

// Test dot
template<typename TYPE>
static bool test_dot( int N, double &error )
{
    typedef typename extended<TYPE>::TYPE2 TYPE2;
    const int K = TEST_SIZE_VEC;
    TYPE *x1    = new TYPE[K];
    TYPE *x2    = new TYPE[K];
    Lapack<TYPE>::random( K, x1 );
    Lapack<TYPE>::random( K, x2 );
    TYPE2 ans1 = 0.0;
    for ( int j = 0; j < K; j++ )
        ans1 += static_cast<TYPE2>( x1[j] ) * static_cast<TYPE2>( x2[j] );
    int N_errors = 0;
    error        = 0;
    double tol   = 50 * epsilon<TYPE>();
    for ( int i = 0; i < N; i++ ) {
        TYPE2 ans2 = Lapack<TYPE>::dot( K, x1, 1, x2, 1 );
        double err = std::abs( ans1 - ans2 ) / K;
        error      = std::max( error, err );
        if ( err > tol )
            N_errors++;
    }
    delete[] x1;
    delete[] x2;
    return N_errors == 0;
}

// Test axpy
template<typename TYPE>
static bool test_axpy( int N, double &error )
{
    const int K = TEST_SIZE_VEC;
    TYPE *x     = new TYPE[K];
    TYPE *y0    = new TYPE[K];
    TYPE *y1    = new TYPE[K];
    TYPE *y2    = new TYPE[K];
    Lapack<TYPE>::random( K, x );
    Lapack<TYPE>::random( K, y0 );
    const TYPE pi = static_cast<TYPE>( 3.141592653589793 );
    for ( int j = 0; j < K; j++ )
        y1[j] = y0[j] + pi * x[j];
    error = 0;
    for ( int i = 0; i < N; i++ ) {
        memcpy( y2, y0, K * sizeof( TYPE ) );
        Lapack<TYPE>::axpy( K, pi, x, 1, y2, 1 );
        double err = L2Error( K, y1, y2 );
        error      = std::max( error, err );
    }
    double tol = sqrt( K ) * epsilon<TYPE>();
    bool pass  = error <= tol;
    delete[] x;
    delete[] y0;
    delete[] y1;
    delete[] y2;
    return pass;
}

// Test gemv
template<typename TYPE>
static bool test_gemv( int N, double &error )
{
    int K1   = TEST_SIZE_MATVEC + 1;
    int K2   = TEST_SIZE_MATVEC - 1;
    TYPE *A  = new TYPE[K1 * K2];
    TYPE *x  = new TYPE[K2];
    TYPE *y  = new TYPE[K1];
    TYPE *y1 = new TYPE[K1];
    TYPE *y2 = new TYPE[K1];
    Lapack<TYPE>::random( K1 * K2, A );
    Lapack<TYPE>::random( K2, x );
    Lapack<TYPE>::random( K1, y );
    const TYPE alpha = static_cast<TYPE>( 3.141592653589793 );
    const TYPE beta  = static_cast<TYPE>( 1.414213562373095 );
    for ( int j = 0; j < K1; j++ ) {
        y1[j] = beta * y[j];
        for ( int k = 0; k < K2; k++ )
            y1[j] += alpha * A[j + k * K1] * x[k];
    }
    int N_errors = 0;
    error        = 0;
    double norm  = Lapack<TYPE>::nrm2( K1, y1, 1 );
    double tol   = K1 * epsilon<TYPE>();
    for ( int i = 0; i < N; i++ ) {
        memcpy( y2, y, K1 * sizeof( TYPE ) );
        Lapack<TYPE>::gemv( 'N', K1, K2, alpha, A, K1, x, 1, beta, y2, 1 );
        error = std::max( error, L2Error( K1, y1, y2 ) / norm );
        if ( !approx_equal( K1, y1, y2, tol ) )
            N_errors++;
    }
    delete[] A;
    delete[] x;
    delete[] y;
    delete[] y1;
    delete[] y2;
    return N_errors == 0;
}

// Test gemm
template<typename TYPE>
static bool test_gemm( int N, double &error )
{
    const int K = TEST_SIZE_MAT;
    TYPE *A     = new TYPE[K * K];
    TYPE *B     = new TYPE[K * K];
    TYPE *C     = new TYPE[K * K];
    TYPE *C1    = new TYPE[K * K];
    TYPE *C2    = new TYPE[K * K];
    Lapack<TYPE>::random( K * K, A );
    Lapack<TYPE>::random( K * K, B );
    Lapack<TYPE>::random( K * K, C );
    const TYPE alpha = static_cast<TYPE>( 3.141592653589793 );
    const TYPE beta  = static_cast<TYPE>( 1.414213562373095 );
    for ( int i = 0; i < K; i++ ) {
        for ( int j = 0; j < K; j++ ) {
            C1[i + j * K] = beta * C[i + j * K];
            for ( int k = 0; k < K; k++ )
                C1[i + j * K] += alpha * A[i + k * K] * B[k + j * K];
        }
    }
    int N_errors = 0;
    error        = 0;
    double norm  = Lapack<TYPE>::nrm2( K * K, C1, 1 );
    double tol   = K * 10 * epsilon<TYPE>();
    for ( int i = 0; i < N; i++ ) {
        memcpy( C2, C, K * K * sizeof( TYPE ) );
        Lapack<TYPE>::gemm( 'N', 'N', K, K, K, alpha, A, K, B, K, beta, C2, K );
        error = std::max( error, L2Error( K * K, C1, C2 ) / norm );
        if ( !approx_equal( K * K, C1, C2, tol ) )
            N_errors++;
    }
    delete[] A;
    delete[] B;
    delete[] C;
    delete[] C1;
    delete[] C2;
    return N_errors == 0;
}


#ifndef DISABLE_LAPACK


// Test gesv
template<typename TYPE>
static bool test_gesv( int N, double &error )
{
    // Test solving a diagonal matrix
    const int K = TEST_SIZE_MAT;
    TYPE *A     = new TYPE[K * K];
    TYPE *x1    = new TYPE[K];
    TYPE *x2    = new TYPE[K];
    TYPE *b     = new TYPE[K];
    int *IPIV   = new int[K];
    memset( (void *) A, 0, K * K * sizeof( TYPE ) );
    Lapack<TYPE>::random( K, x2 );
    Lapack<TYPE>::random( K, b );
    for ( int k = 0; k < K; k++ ) {
        A[k + k * K] = x2[k] + epsilon<TYPE>();
        x1[k]        = b[k] / ( x2[k] + epsilon<TYPE>() );
    }
    int N_errors = 0;
    error        = 0;
    double norm  = Lapack<TYPE>::nrm2( K, x1, 1 );
    double tol   = K * 10 * epsilon<TYPE>();
    for ( int i = 0; i < N; i++ ) {
        memcpy( x2, b, K * sizeof( TYPE ) );
        int err = 0;
        Lapack<TYPE>::gesv( K, 1, A, K, IPIV, x2, K, err );
        N_errors += err == 0 ? 0 : 1;
        error = std::max( error, L2Error( K, x1, x2 ) / norm );
        if ( !approx_equal( K, x1, x2, tol ) )
            N_errors++;
    }
    delete[] A;
    delete[] x1;
    delete[] x2;
    delete[] b;
    delete[] IPIV;
    return N_errors == 0;
}

// Test gtsv
template<typename TYPE>
static bool test_gtsv( int N, double &error )
{
    // Test solving a tri-diagonal matrix by comparing to dgtsv
    const int K = TEST_SIZE_TRI_MAT / 2;
    TYPE *A     = new TYPE[K * K];
    TYPE *D     = new TYPE[K];
    TYPE *D2    = new TYPE[K];
    TYPE *DL    = new TYPE[K - 1];
    TYPE *DL2   = new TYPE[K - 1];
    TYPE *DU    = new TYPE[K - 1];
    TYPE *DU2   = new TYPE[K - 1];
    TYPE *x1    = new TYPE[K];
    TYPE *x2    = new TYPE[K];
    TYPE *b     = new TYPE[K];
    int *IPIV   = new int[K];
    generateRhs( K, b );
    generateMatrixBanded( K, 1, 1, A );
    extractTriDiag( K, A, DL, D, DU );
    memcpy( x1, b, K * sizeof( TYPE ) );
    int err = 0;
    Lapack<TYPE>::gesv( K, 1, A, K, IPIV, x1, K, err );
    int N_errors = 0;
    error        = 0;
    for ( int i = 0; i < N; i++ ) {
        memcpy( x2, b, K * sizeof( TYPE ) );
        memcpy( D2, D, K * sizeof( TYPE ) );
        memcpy( DL2, DL, ( K - 1 ) * sizeof( TYPE ) );
        memcpy( DU2, DU, ( K - 1 ) * sizeof( TYPE ) );
        Lapack<TYPE>::gtsv( K, 1, DL2, D2, DU2, x2, K, err );
        N_errors += err == 0 ? 0 : 1;
        if ( err != 0 )
            printf( "Error calling gtsv (%i)\n", err );
        double err2 = L2Error( N, x1, x2 );
        double norm = Lapack<TYPE>::nrm2( N, x1, 1 );
        error       = std::max( error, err2 / norm );
    }
    double tol = 2e4 * epsilon<TYPE>();
    if ( error > tol ) {
        printf( "test_gtsv error (%e) exceeded tolerance (%e)\n", error, tol );
        N_errors++;
    }
    delete[] A;
    delete[] D;
    delete[] D2;
    delete[] DL;
    delete[] DL2;
    delete[] DU;
    delete[] DU2;
    delete[] x1;
    delete[] x2;
    delete[] b;
    delete[] IPIV;
    return N_errors == 0;
}
// Test gbsv
template<typename TYPE>
static bool test_gbsv( int N, double &error )
{
    // Test solving a banded-diagonal matrix by comparing to dgtsv
    //    N = 6, KL = 2, KU = 1:
    //        *    *    *    +    +    +
    //        *    *    +    +    +    +
    //        *   a12  a23  a34  a45  a56
    //       a11  a22  a33  a44  a55  a66
    //       a21  a32  a43  a54  a65   *
    //       a31  a42  a53  a64   *    *
    const int K  = TEST_SIZE_TRI_MAT / 2;
    const int KL = 2;
    const int KU = 2;
    const int K2 = 2 * KL + KU + 1;
    TYPE *A      = new TYPE[K * K];
    TYPE *AB     = new TYPE[K * K2];
    TYPE *AB2    = new TYPE[K * K2];
    TYPE *x1     = new TYPE[K];
    TYPE *x2     = new TYPE[K];
    TYPE *b      = new TYPE[K];
    int *IPIV    = new int[K];
    generateRhs( K, b );
    generateMatrixBanded( K, KL, KU, A );
    extractBanded( K, KL, KU, A, AB );
    memcpy( x1, b, K * sizeof( TYPE ) );
    int err = 0;
    error   = 0;
    Lapack<TYPE>::gesv( K, 1, A, K, IPIV, x1, K, err );
    int N_errors = 0;
    for ( int i = 0; i < N; i++ ) {
        memcpy( x2, b, K * sizeof( TYPE ) );
        memcpy( AB2, AB, K * K2 * sizeof( TYPE ) );
        Lapack<TYPE>::gbsv( K, KL, KU, 1, AB2, K2, IPIV, x2, K, err );
        N_errors += err == 0 ? 0 : 1;
        double norm = Lapack<TYPE>::nrm2( K, x1, 1 );
        double err2 = L2Error( K, x1, x2 );
        error       = std::max( error, err2 / norm );
    }
    const double tol = 4 * K * sqrt( K ) * epsilon<TYPE>();
    if ( error > tol ) {
        printf( "   test_gbsv error (%e) exceeded tolerance (%e)\n", error, tol );
        N_errors++;
    }
    delete[] A;
    delete[] AB;
    delete[] AB2;
    delete[] x1;
    delete[] x2;
    delete[] b;
    delete[] IPIV;
    return N_errors == 0;
}

// Test getrf
template<typename TYPE>
static bool test_getrf( int N, double &error )
{
    // Check dgetrf by performing a factorization and solve and comparing to dgesv
    const int K = TEST_SIZE_MAT;
    TYPE *A     = new TYPE[K * K];
    TYPE *A2    = new TYPE[K * K];
    TYPE *x1    = new TYPE[K];
    TYPE *x2    = new TYPE[K];
    TYPE *b     = new TYPE[K];
    int *IPIV   = new int[K];
    generateRhs( K, b );
    generateMatrix( K, A );
    memcpy( A2, A, K * K * sizeof( TYPE ) );
    memcpy( x1, b, K * sizeof( TYPE ) );
    int err = 0;
    Lapack<TYPE>::gesv( K, 1, A2, K, IPIV, x1, K, err );
    int N_errors = 0;
    for ( int i = 0; i < N; i++ ) {
        memcpy( A2, A, K * K * sizeof( TYPE ) );
        Lapack<TYPE>::getrf( K, K, A2, K, IPIV, err );
        N_errors += err == 0 ? 0 : 1;
    }
    memcpy( x2, b, K * sizeof( TYPE ) );
    Lapack<TYPE>::getrs( 'N', K, 1, A2, K, IPIV, x2, K, err );
    double norm = Lapack<TYPE>::nrm2( K, x1, 1 );
    double err2 = L2Error( K, x1, x2 );
    double tol  = 10.0 * norm * epsilon<TYPE>();
    if ( err2 > tol )
        N_errors++;
    error = err2 / norm;
    delete[] A;
    delete[] A2;
    delete[] x1;
    delete[] x2;
    delete[] b;
    delete[] IPIV;
    return N_errors == 0;
}

// Test gttrf
template<typename TYPE>
static bool test_gttrf( int N, double &error )
{
    // Check dgttrf by performing a factorization and solve and comparing to dgtsv
    const int K = 5 * TEST_SIZE_TRI;
    TYPE *D     = new TYPE[K];
    TYPE *D2    = new TYPE[K];
    TYPE *DL    = new TYPE[K - 1];
    TYPE *DL2   = new TYPE[K - 1];
    TYPE *DU    = new TYPE[K - 1];
    TYPE *DU2   = new TYPE[K - 1];
    TYPE *DU3   = new TYPE[K - 2];
    TYPE *x1    = new TYPE[K];
    TYPE *x2    = new TYPE[K];
    TYPE *b     = new TYPE[K];
    int *IPIV   = new int[K];
    generateRhs( K, b );
    generateTriDiag( K, DL, D, DU );
    memcpy( x1, b, K * sizeof( TYPE ) );
    memcpy( D2, D, K * sizeof( TYPE ) );
    memcpy( DL2, DL, ( K - 1 ) * sizeof( TYPE ) );
    memcpy( DU2, DU, ( K - 1 ) * sizeof( TYPE ) );
    int err = 0;
    Lapack<TYPE>::gtsv( K, 1, DL2, D2, DU2, x1, K, err );
    int N_errors = 0;
    for ( int i = 0; i < N; i++ ) {
        memcpy( D2, D, K * sizeof( TYPE ) );
        memcpy( DL2, DL, ( K - 1 ) * sizeof( TYPE ) );
        memcpy( DU2, DU, ( K - 1 ) * sizeof( TYPE ) );
        Lapack<TYPE>::gttrf( K, DL2, D2, DU2, DU3, IPIV, err );
        N_errors += err == 0 ? 0 : 1;
    }
    memcpy( x2, b, K * sizeof( TYPE ) );
    Lapack<TYPE>::gttrs( 'N', K, 1, DL2, D2, DU2, DU3, IPIV, x2, K, err );
    double norm = Lapack<TYPE>::nrm2( K, x1, 1 );
    double err2 = L2Error( K, x1, x2 ) / norm;
    double tol  = 100 * K * epsilon<TYPE>();
    if ( err2 > tol ) {
        printf( "   gttrf exceeded tolerance: error = %e, tol = %e\n", err2, tol );
        N_errors++;
    }
    error = err2 / norm;
    delete[] D;
    delete[] D2;
    delete[] DL;
    delete[] DL2;
    delete[] DU;
    delete[] DU2;
    delete[] DU3;
    delete[] x1;
    delete[] x2;
    delete[] b;
    delete[] IPIV;
    return N_errors == 0;
}

// Test gbtrf
template<typename TYPE>
static bool test_gbtrf( int N, double &error )
{
    // Check dgbtrf by performing a factorization and solve and comparing to dgbsv
    const int K  = TEST_SIZE_TRI;
    const int KL = 2;
    const int KU = 2;
    const int K2 = 2 * KL + KU + 1;
    TYPE *AB     = new TYPE[K * K2];
    TYPE *AB2    = new TYPE[K * K2];
    TYPE *x1     = new TYPE[K];
    TYPE *x2     = new TYPE[K];
    TYPE *b      = new TYPE[K];
    int *IPIV    = new int[K];
    Lapack<TYPE>::random( K * K2, AB );
    Lapack<TYPE>::random( K, b );
    int err = 0;
    memcpy( x1, b, K * sizeof( TYPE ) );
    memcpy( AB2, AB, K * K2 * sizeof( TYPE ) );
    Lapack<TYPE>::gbsv( K, KL, KU, 1, AB2, K2, IPIV, x1, K, err );
    int N_errors = 0;
    for ( int i = 0; i < N; i++ ) {
        memcpy( AB2, AB, K * K2 * sizeof( TYPE ) );
        Lapack<TYPE>::gbtrf( K, K, KL, KU, AB2, K2, IPIV, err );
        N_errors += err == 0 ? 0 : 1;
    }
    memcpy( x2, b, K * sizeof( TYPE ) );
    Lapack<TYPE>::gbtrs( 'N', K, KL, KU, 1, AB2, K2, IPIV, x2, K, err );
    double norm = Lapack<TYPE>::nrm2( K, x1, 1 );
    double err2 = L2Error( K, x1, x2 );
    double tol  = 10 * norm * epsilon<TYPE>();
    if ( err2 > tol )
        N_errors++;
    error = err2 / norm;
    delete[] AB;
    delete[] AB2;
    delete[] x1;
    delete[] x2;
    delete[] b;
    delete[] IPIV;
    return N_errors == 0;
}

// Test getrs
template<typename TYPE>
static bool test_getrs( int N, double &error )
{
    // Check dgetrs by performing a factorization and solve and comparing to dgesv
    const int K = 2 * TEST_SIZE_MAT;
    TYPE *A     = new TYPE[K * K];
    TYPE *A2    = new TYPE[K * K];
    TYPE *x1    = new TYPE[K];
    TYPE *x2    = new TYPE[K];
    TYPE *b     = new TYPE[K];
    int *IPIV   = new int[K];
    Lapack<TYPE>::random( K * K, A );
    Lapack<TYPE>::random( K, b );
    int err = 0;
    error   = 0;
    memcpy( A2, A, K * K * sizeof( TYPE ) );
    memcpy( x1, b, K * sizeof( TYPE ) );
    Lapack<TYPE>::gesv( K, 1, A2, K, IPIV, x1, K, err );
    int N_errors = 0;
    Lapack<TYPE>::getrf( K, K, A, K, IPIV, err );
    for ( int i = 0; i < N; i++ ) {
        memcpy( A2, A, K * K * sizeof( TYPE ) );
        memcpy( x2, b, K * sizeof( TYPE ) );
        Lapack<TYPE>::getrs( 'N', K, 1, A2, K, IPIV, x2, K, err );
        N_errors += err == 0 ? 0 : 1;
        double norm = Lapack<TYPE>::nrm2( K, x1, 1 );
        double err2 = L2Error( K, x1, x2 );
        double tol  = 10.0 * norm * epsilon<TYPE>();
        if ( err > tol )
            N_errors++;
        error = std::max( error, err2 / norm );
    }
    delete[] A;
    delete[] A2;
    delete[] x1;
    delete[] x2;
    delete[] b;
    delete[] IPIV;
    return N_errors == 0;
}

// Test gttrs
template<typename TYPE>
static bool test_gttrs( int N, double &error )
{
    // Check dgttrs by performing a factorization and solve and comparing to dgtsv
    const int K = 5 * TEST_SIZE_TRI;
    TYPE *D     = new TYPE[K];
    TYPE *D2    = new TYPE[K];
    TYPE *DL    = new TYPE[K - 1];
    TYPE *DL2   = new TYPE[K - 1];
    TYPE *DU    = new TYPE[K - 1];
    TYPE *DU2   = new TYPE[K - 1];
    TYPE *DU3   = new TYPE[K - 2];
    TYPE *DU4   = new TYPE[K - 2];
    TYPE *x1    = new TYPE[K];
    TYPE *x2    = new TYPE[K];
    TYPE *b     = new TYPE[K];
    int *IPIV   = new int[K];
    generateRhs( K, b );
    generateTriDiag( K, DL, D, DU );
    int err = 0;
    error   = 0;
    memcpy( x1, b, K * sizeof( TYPE ) );
    memcpy( D2, D, K * sizeof( TYPE ) );
    memcpy( DL2, DL, ( K - 1 ) * sizeof( TYPE ) );
    memcpy( DU2, DU, ( K - 1 ) * sizeof( TYPE ) );
    Lapack<TYPE>::gtsv( K, 1, DL2, D2, DU2, x1, K, err );
    Lapack<TYPE>::gttrf( K, DL, D, DU, DU3, IPIV, err );
    int N_errors = 0;
    for ( int i = 0; i < N; i++ ) {
        memcpy( D2, D, K * sizeof( TYPE ) );
        memcpy( DL2, DL, ( K - 1 ) * sizeof( TYPE ) );
        memcpy( DU2, DU, ( K - 1 ) * sizeof( TYPE ) );
        memcpy( DU4, DU3, ( K - 2 ) * sizeof( TYPE ) );
        memcpy( x2, b, K * sizeof( TYPE ) );
        Lapack<TYPE>::gttrs( 'N', K, 1, DL2, D2, DU2, DU4, IPIV, x2, K, err );
        N_errors += err == 0 ? 0 : 1;
        double norm = Lapack<TYPE>::nrm2( K, x1, 1 );
        double err2 = L2Error( K, x1, x2 );
        double tol  = 300 * K * norm * epsilon<TYPE>();
        if ( err2 > tol )
            N_errors++;
        error = std::max( error, err2 / norm );
    }
    delete[] D;
    delete[] D2;
    delete[] DL;
    delete[] DL2;
    delete[] DU;
    delete[] DU2;
    delete[] DU3;
    delete[] DU4;
    delete[] x1;
    delete[] x2;
    delete[] b;
    delete[] IPIV;
    return N_errors == 0;
}

// Test gbtrs
template<typename TYPE>
static bool test_gbtrs( int N, double &error )
{
    // Check dgbtrs by performing a factorization and solve and comparing to dgbsv
    const int K  = TEST_SIZE_TRI;
    const int KL = 2;
    const int KU = 2;
    const int K2 = 2 * KL + KU + 1;
    TYPE *AB     = new TYPE[K * K2];
    TYPE *AB2    = new TYPE[K * K2];
    TYPE *x1     = new TYPE[K];
    TYPE *x2     = new TYPE[K];
    TYPE *b      = new TYPE[K];
    int *IPIV    = new int[K];
    Lapack<TYPE>::random( K * K2, AB );
    Lapack<TYPE>::random( K, b );
    int err = 0;
    error   = 0;
    memcpy( x1, b, K * sizeof( TYPE ) );
    memcpy( AB2, AB, K * K2 * sizeof( TYPE ) );
    Lapack<TYPE>::gbsv( K, KL, KU, 1, AB2, K2, IPIV, x1, K, err );
    Lapack<TYPE>::gbtrf( K, K, KL, KU, AB, K2, IPIV, err );
    int N_errors = 0;
    for ( int i = 0; i < N; i++ ) {
        memcpy( AB2, AB, K * K2 * sizeof( TYPE ) );
        memcpy( x2, b, K * sizeof( TYPE ) );
        Lapack<TYPE>::gbtrs( 'N', K, KL, KU, 1, AB2, K2, IPIV, x2, K, err );
        N_errors += err == 0 ? 0 : 1;
        double norm = Lapack<TYPE>::nrm2( K, x1, 1 );
        double err2 = L2Error( K, x1, x2 );
        double tol  = sqrt( K ) * norm * epsilon<TYPE>();
        if ( err2 > tol )
            N_errors++;
        error = std::max( error, err2 / norm );
    }
    delete[] AB;
    delete[] AB2;
    delete[] x1;
    delete[] x2;
    delete[] b;
    delete[] IPIV;
    return N_errors == 0;
}

// Test getri
template<typename TYPE>
static bool test_getri( int N, double &error )
{
    // Check getri by performing a factorization, calculating the inverse,
    //   multiplying the rhs, and comparing to gesv
    const int K     = TEST_SIZE_MAT;
    const int LWORK = 8 * K;
    TYPE *A         = new TYPE[K * K];
    TYPE *A2        = new TYPE[K * K];
    TYPE *x1        = new TYPE[K];
    TYPE *x2        = new TYPE[K];
    TYPE *b         = new TYPE[K];
    int *IPIV       = new int[K];
    TYPE *WORK      = new TYPE[LWORK];
    Lapack<TYPE>::random( K * K, A );
    Lapack<TYPE>::random( K, b );
    int err      = 0;
    error        = 0;
    int N_errors = 0;
    memcpy( A2, A, K * K * sizeof( TYPE ) );
    memcpy( x1, b, K * sizeof( TYPE ) );
    Lapack<TYPE>::gesv( K, 1, A2, K, IPIV, x1, K, err );
    if ( err != 0 )
        printf( "Error in gesv within test_getri\n" );
    double norm = Lapack<TYPE>::nrm2( K, x1, 1 );
    // Compute LU factorization
    Lapack<TYPE>::getrf( K, K, A, K, IPIV, err );
    if ( err != 0 )
        printf( "Error in getrf within test_getri\n" );
    for ( int i = 0; i < N; i++ ) {
        // Compute the inverse
        memcpy( A2, A, K * K * sizeof( TYPE ) );
        Lapack<TYPE>::getri( K, A2, K, IPIV, WORK, LWORK, err );
        if ( err != 0 ) {
            printf( "Error in getri within test_getri\n" );
            N_errors++;
        }
        // Perform the mat-vec
        memset( (void *) x2, 0xB6, K * sizeof( TYPE ) );
        Lapack<TYPE>::gemv( 'N', K, K, 1, A2, K, b, 1, 0, x2, 1 );
        double err2 = L2Error( K, x1, x2 ) / norm;
        error       = std::max( error, err2 );
    }
    // Check the result
    double tol = 10 * K * epsilon<TYPE>();
    if ( error > tol ) {
        printf( "   getri exceeded tolerance: error = %e, tol = %e\n", error, tol );
        N_errors++;
    }
    delete[] A;
    delete[] A2;
    delete[] x1;
    delete[] x2;
    delete[] b;
    delete[] IPIV;
    delete[] WORK;
    return N_errors == 0;
}


#else


// Tests that are not valid without lapack
template<typename TYPE>
static bool test_gesv( int, double & )
{
    printf( "test is disabled without external Lapack" );
    return false;
}
template<typename TYPE>
static bool test_gtsv( int, double & )
{
    printf( "test is disabled without external Lapack" );
    return false;
}
template<typename TYPE>
static bool test_gbsv( int, double & )
{
    printf( "test is disabled without external Lapack" );
    return false;
}
template<typename TYPE>
static bool test_getrf( int, double & )
{
    printf( "test is disabled without external Lapack" );
    return false;
}
template<typename TYPE>
static bool test_gttrf( int, double & )
{
    printf( "test is disabled without external Lapack" );
    return false;
}
template<typename TYPE>
static bool test_gbtrf( int, double & )
{
    printf( "test is disabled without external Lapack" );
    return false;
}
template<typename TYPE>
static bool test_getrs( int, double & )
{
    printf( "test is disabled without external Lapack" );
    return false;
}
template<typename TYPE>
static bool test_gttrs( int, double & )
{
    printf( "test is disabled without external Lapack" );
    return false;
}
template<typename TYPE>
static bool test_gbtrs( int, double & )
{
    printf( "test is disabled without external Lapack" );
    return false;
}
template<typename TYPE>
static bool test_getri( int, double & )
{
    printf( "test is disabled without external Lapack" );
    return false;
}


#endif


/******************************************************************
 * Run a given test                                                *
 ******************************************************************/
template<class TYPE>
bool Lapack<TYPE>::run_test( const std::string &routine, int N, double &error )
{
    auto name = routine.substr( 1 );
    std::transform( name.begin(), name.end(), name.begin(), ::tolower );
    bool pass = true;
    if ( name == "copy" ) {
        pass = test_copy<TYPE>( N, error );
    } else if ( name == "scal" ) {
        pass = test_scal<TYPE>( N, error );
    } else if ( name == "nrm2" ) {
        pass = test_nrm2<TYPE>( N, error );
    } else if ( name == "axpy" ) {
        pass = test_axpy<TYPE>( N, error );
    } else if ( name == "gemv" ) {
        pass = test_gemv<TYPE>( N, error );
    } else if ( name == "gemm" ) {
        pass = test_gemm<TYPE>( N, error );
    } else if ( name == "asum" ) {
        pass = test_asum<TYPE>( N, error );
    } else if ( name == "dot" ) {
        pass = test_dot<TYPE>( N, error );
    } else if ( name == "gesv" ) {
        pass = test_gesv<TYPE>( N, error );
    } else if ( name == "gtsv" ) {
        pass = test_gtsv<TYPE>( N, error );
    } else if ( name == "gbsv" ) {
        pass = test_gbsv<TYPE>( N, error );
    } else if ( name == "getrf" ) {
        pass = test_getrf<TYPE>( N, error );
    } else if ( name == "gttrf" ) {
        pass = test_gttrf<TYPE>( N, error );
    } else if ( name == "gbtrf" ) {
        pass = test_gbtrf<TYPE>( N, error );
    } else if ( name == "getrs" ) {
        pass = test_getrs<TYPE>( N, error );
    } else if ( name == "gttrs" ) {
        pass = test_gttrs<TYPE>( N, error );
    } else if ( name == "gbtrs" ) {
        pass = test_gbtrs<TYPE>( N, error );
    } else if ( name == "getri" ) {
        pass = test_getri<TYPE>( N, error );
    } else if ( name == "rand" ) {
        pass = test_random<TYPE>( N, error );
    } else {
        std::cerr << "Unknown test: " << name << std::endl;
        return false;
    }
    return pass;
}
template bool Lapack<float>::run_test( const std::string &, int, double & );
template bool Lapack<double>::run_test( const std::string &, int, double & );
template bool Lapack<std::complex<double>>::run_test( const std::string &, int, double & );
