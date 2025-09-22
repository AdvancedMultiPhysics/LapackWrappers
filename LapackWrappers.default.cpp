#include "LapackWrappers.h"

#include <complex>
#include <math.h>
#include <stdexcept>


#define ASSERT( EXP )                                            \
    do {                                                         \
        if ( !( EXP ) ) {                                        \
            throw std::logic_error( "Failed assertion: " #EXP ); \
        }                                                        \
    } while ( 0 )


// Helper functions
static inline bool getTrans( char TRANS )
{
    if ( TRANS == 'N' || TRANS == 'n' ) {
        return false;
    } else if ( TRANS == 'T' || TRANS == 'T' || TRANS == 'C' || TRANS == 'C' ) {
        return true;
    } else {
        throw std::logic_error( "invalid value for TRANS" );
    }
}


// lamch
template<typename TYPE>
TYPE Lapack<TYPE>::lamch( char cmach )
{
    if ( cmach == 'E' ) {
        static_assert( std::numeric_limits<TYPE>::radix == 2 );
        return 0.5 * std::numeric_limits<TYPE>::epsilon();
    } else if ( cmach == 'S' ) {
        return std::numeric_limits<TYPE>::min();
    } else if ( cmach == 'B' ) {
        return std::numeric_limits<TYPE>::radix;
    } else if ( cmach == 'P' ) {
        return std::numeric_limits<TYPE>::epsilon();
    } else if ( cmach == 'N' ) {
        return std::numeric_limits<TYPE>::digits;
    } else if ( cmach == 'R' ) {
        return 1;
    } else if ( cmach == 'M' ) {
        return std::numeric_limits<TYPE>::min_exponent;
    } else if ( cmach == 'U' ) {
        return std::numeric_limits<TYPE>::min();
    } else if ( cmach == 'L' ) {
        return std::numeric_limits<TYPE>::max_exponent;
    } else if ( cmach == 'O' ) {
        return std::numeric_limits<TYPE>::max();
    } else {
        return 0;
    }
}
template float Lapack<float>::lamch( char );
template double Lapack<double>::lamch( char );


// Define the member functions
template<class TYPE>
void Lapack<TYPE>::copy( int N, const TYPE *DX, int INCX, TYPE *DY, int INCY )
{
    auto X = DX;
    auto Y = DY;
    for ( int i = 0; i < N; i++, X += INCX, Y += INCY )
        *Y = *X;
}
// Define the member functions
template<class TYPE>
void Lapack<TYPE>::swap( int N, TYPE *DX, int INCX, TYPE *DY, int INCY )
{
    auto X = DX;
    auto Y = DY;
    for ( int i = 0; i < N; i++, X += INCX, Y += INCY ) {
        TYPE Z = *Y;
        *Y     = *X;
        *X     = Z;
    }
}
template<class TYPE>
void Lapack<TYPE>::scal( int N, TYPE DA, TYPE *DX, int INCX )
{
    auto X = DX;
    for ( int i = 0; i < N; i++, X += INCX )
        *X *= DA;
}
template<class TYPE>
double Lapack<TYPE>::nrm2( int N, const TYPE *DX, int INCX )
{
    double s = 0;
    auto X   = DX;
    for ( int i = 0; i < N; i++, X += INCX )
        s += std::norm( *X );
    return sqrt( s );
}
template<class TYPE>
int Lapack<TYPE>::iamax( int N, const TYPE *DX, int INCX )
{
    if ( N <= 1 )
        return N;
    int k  = 0;
    auto X = DX + INCX;
    auto y = std::abs( *DX );
    for ( int i = 1; i < N; i++, X += INCX ) {
        auto x = std::abs( *X );
        if ( x > y ) {
            k = i;
            y = x;
        }
    }
    return k;
}
template<class TYPE>
void Lapack<TYPE>::axpy( int N, TYPE DA, const TYPE *DX, int INCX, TYPE *DY, int INCY )
{
    auto X = DX;
    auto Y = DY;
    for ( int i = 0; i < N; i++, X += INCX, Y += INCY )
        *Y += DA * ( *X );
}
template<class TYPE>
double Lapack<TYPE>::asum( int N, const TYPE *DX, int INCX )
{
    if ( N <= 1 )
        return N;
    auto X   = DX;
    double s = 0;
    for ( int i = 0; i < N; i++, X += INCX )
        s += std::abs( *X );
    return s;
}
template<class TYPE>
TYPE Lapack<TYPE>::dot( int N, const TYPE *DX, int INCX, const TYPE *DY, int INCY )
{
    TYPE s = 0;
    auto X = DX;
    auto Y = DY;
    for ( int i = 0; i < N; i++, X += INCX, Y += INCY )
        s += ( *X ) * ( *Y );
    return s;
}
template<class TYPE>
void Lapack<TYPE>::gemv( char TRANS, int M, int N, TYPE alpha, const TYPE *A, int LDA,
    const TYPE *DX, int INCX, TYPE beta, TYPE *DY, int INCY )
{
    ASSERT( M >= 0 );
    ASSERT( N >= 0 );
    ASSERT( LDA >= 0 );
    ASSERT( INCX >= 1 );
    ASSERT( INCY >= 1 );
    bool trans = getTrans( TRANS );
    int Nx     = trans ? M : N;
    int Ny     = trans ? N : M;
    for ( int i = 0, iy = 0; i < Ny; i++, iy += INCY )
        DY[iy] = beta * DY[iy];
    constexpr TYPE zero( 0 );
    if ( alpha == zero )
        return;
    if ( TRANS == 'N' || TRANS == 'n' ) {
        for ( int i = 0, iy = 0; i < Ny; i++, iy += INCY ) {
            TYPE Ax = zero;
            for ( int j = 0, jx = 0; j < Nx; j++, jx += INCX )
                Ax += A[i + j * LDA] * DX[jx];
            DY[iy] += alpha * Ax;
        }
    } else {
        for ( int i = 0, iy = 0; i < Ny; i++, iy += INCY ) {
            TYPE Ax = zero;
            for ( int j = 0, jx = 0; j < Nx; j++, jx += INCX )
                Ax += A[j + i * LDA] * DX[jx];
            DY[iy] += alpha * Ax;
        }
    }
}
template<class TYPE>
void Lapack<TYPE>::gemm( char TRANSA, char TRANSB, int M, int N, int K, TYPE alpha, const TYPE *A,
    int LDA, const TYPE *B, int LDB, TYPE beta, TYPE *C, int LDC )
{
    bool transa = getTrans( TRANSA );
    bool transb = getTrans( TRANSB );
    int nrowa   = transa ? K : M;
    int nrowb   = transb ? N : K;
    ASSERT( M >= 0 );
    ASSERT( N >= 0 );
    ASSERT( K >= 0 );
    ASSERT( LDA >= std::max( 1, nrowa ) );
    ASSERT( LDB >= std::max( 1, nrowb ) );
    ASSERT( LDC >= std::max( 1, M ) );
    for ( int i = 0; i < N * M; i++ )
        C[i] = beta * C[i];
    constexpr TYPE zero( 0 );
    if ( alpha == zero )
        return;
    if ( !transb ) {
        if ( !transa ) {
            // C := alpha*A*B + beta*C.
            for ( int j = 0; j < N; j++ ) {
                for ( int l = 0; l < K; l++ ) {
                    TYPE temp = alpha * B[l + j * LDB];
                    for ( int i = 0; i < M; i++ ) {
                        C[i + j * LDC] += temp * A[i + l * LDA];
                    }
                }
            }
        } else {
            // C := alpha*A**T*B + beta*C
            for ( int j = 0; j < N; j++ ) {
                for ( int i = 0; i < M; i++ ) {
                    TYPE temp = zero;
                    for ( int l = 0; l < K; l++ ) {
                        temp = temp + A[l + i * LDA] * B[l + j * LDB];
                    }
                    C[i + j * LDC] += alpha * temp;
                }
            }
        }
    } else {
        if ( !transa ) {
            // C := alpha*A*B**T + beta*C
            for ( int j = 0; j < N; j++ ) {
                for ( int l = 0; l < K; l++ ) {
                    TYPE temp = alpha * B[j + l * LDB];
                    for ( int i = 0; i < M; i++ ) {
                        C[i + j * LDC] += temp * A[i + l * LDA];
                    }
                }
            }
        } else {
            // C := alpha*A**T*B**T + beta*C
            for ( int j = 0; j < N; j++ ) {
                for ( int i = 0; i < M; i++ ) {
                    TYPE temp = zero;
                    for ( int l = 0; l < K; l++ ) {
                        temp = temp + A[l + i * LDA] * B[j + l * LDB];
                    }
                    C[i + j * LDC] += alpha * temp;
                }
            }
        }
    }
}
template<class TYPE>
void Lapack<TYPE>::ger( int, int, TYPE, const TYPE *, int, const TYPE *, int, TYPE *, int )
{
    throw std::logic_error( "ger is not currently supported without blas/lapack" );
}
template<class TYPE>
void Lapack<TYPE>::gesv( int, int, TYPE *, int, int *, TYPE *, int, int & )
{
    throw std::logic_error( "gesv is not currently supported without blas/lapack" );
}
template<class TYPE>
void Lapack<TYPE>::gtsv( int, int, TYPE *, TYPE *, TYPE *, TYPE *, int, int & )
{
    throw std::logic_error( "gtsv is not currently supported without blas/lapack" );
}
template<class TYPE>
void Lapack<TYPE>::gbsv( int, int, int, int, TYPE *, int, int *, TYPE *, int, int & )
{
    throw std::logic_error( "gbsv is not currently supported without blas/lapack" );
}
template<class TYPE>
void Lapack<TYPE>::getrf( int, int, TYPE *, int, int *, int & )
{
    throw std::logic_error( "getrf is not currently supported without blas/lapack" );
}
template<class TYPE>
void Lapack<TYPE>::gttrf( int, TYPE *, TYPE *, TYPE *, TYPE *, int *, int & )
{
    throw std::logic_error( "gttrf is not currently supported without blas/lapack" );
}
template<class TYPE>
void Lapack<TYPE>::gbtrf( int, int, int, int, TYPE *, int, int *, int & )
{
    throw std::logic_error( "gbtrf is not currently supported without blas/lapack" );
}
template<class TYPE>
void Lapack<TYPE>::getrs( char, int, int, const TYPE *, int, const int *, TYPE *, int, int & )
{
    throw std::logic_error( "getrs is not currently supported without blas/lapack" );
}
template<class TYPE>
void Lapack<TYPE>::gttrs( char, int, int, const TYPE *, const TYPE *, const TYPE *, const TYPE *,
    const int *, TYPE *, int, int & )
{
    throw std::logic_error( "gttrs is not currently supported without blas/lapack" );
}
template<class TYPE>
void Lapack<TYPE>::gbtrs(
    char, int, int, int, int, const TYPE *, int, const int *, TYPE *, int, int & )
{
    throw std::logic_error( "gbtrs is not currently supported without blas/lapack" );
}
template<class TYPE>
void Lapack<TYPE>::getri( int, TYPE *, int, const int *, TYPE *, int, int & )
{
    throw std::logic_error( "getri is not currently supported without blas/lapack" );
}
template<class TYPE>
void Lapack<TYPE>::trsm( char, char, char, char, int, int, TYPE, const TYPE *, int, TYPE *, int )
{
    throw std::logic_error( "trsm is not currently supported without blas/lapack" );
}


// Explicit instatiations
#define INSTANTIATE( TYPE )                                                                       \
    template void Lapack<TYPE>::copy( int, const TYPE *, int, TYPE *, int );                      \
    template void Lapack<TYPE>::swap( int, TYPE *, int, TYPE *, int );                            \
    template void Lapack<TYPE>::scal( int, TYPE, TYPE *, int );                                   \
    template double Lapack<TYPE>::nrm2( int, const TYPE *, int );                                 \
    template int Lapack<TYPE>::iamax( int, const TYPE *, int );                                   \
    template void Lapack<TYPE>::axpy( int, TYPE, const TYPE *, int, TYPE *, int );                \
    template void Lapack<TYPE>::gemv(                                                             \
        char, int, int, TYPE, const TYPE *, int, const TYPE *, int, TYPE, TYPE *, int );          \
    template void Lapack<TYPE>::gemm( char, char, int, int, int K, TYPE ALPHA, const TYPE *, int, \
        const TYPE *, int, TYPE, TYPE *, int );                                                   \
    template double Lapack<TYPE>::asum( int, const TYPE *, int );                                 \
    template TYPE Lapack<TYPE>::dot( int, const TYPE *, int, const TYPE *, int );                 \
    template void Lapack<TYPE>::ger(                                                              \
        int, int, TYPE, const TYPE *, int, const TYPE *, int, TYPE *, int );                      \
    template void Lapack<TYPE>::gesv( int, int, TYPE *, int, int *, TYPE *, int, int & );         \
    template void Lapack<TYPE>::gtsv( int, int, TYPE *, TYPE *, TYPE *, TYPE *, int, int & );     \
    template void Lapack<TYPE>::gbsv(                                                             \
        int, int, int, int, TYPE *, int, int *, TYPE *, int, int & );                             \
    template void Lapack<TYPE>::getrf( int, int, TYPE *, int, int *, int & );                     \
    template void Lapack<TYPE>::gttrf( int, TYPE *, TYPE *, TYPE *, TYPE *, int *, int & );       \
    template void Lapack<TYPE>::gbtrf( int, int, int, int, TYPE *, int, int *, int & );           \
    template void Lapack<TYPE>::getrs(                                                            \
        char, int, int, const TYPE *, int, const int *, TYPE *, int, int & );                     \
    template void Lapack<TYPE>::gttrs( char, int, int, const TYPE *, const TYPE *, const TYPE *,  \
        const TYPE *, const int *, TYPE *, int, int & );                                          \
    template void Lapack<TYPE>::gbtrs(                                                            \
        char, int, int, int, int, const TYPE *, int, const int *, TYPE *, int, int & );           \
    template void Lapack<TYPE>::getri( int, TYPE *, int, const int *, TYPE *, int, int & );       \
    template void Lapack<TYPE>::trsm(                                                             \
        char, char, char, char, int, int, TYPE, const TYPE *, int, TYPE *, int )

INSTANTIATE( float );
INSTANTIATE( double );
INSTANTIATE( std::complex<double> );
