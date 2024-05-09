#include "LapackWrappers.h"

#include <complex>
#include <math.h>
#include <stdexcept>


#define NULL_USE( variable )                 \
    do {                                     \
        if ( 0 ) {                           \
            char *temp = (char *) &variable; \
            temp++;                          \
        }                                    \
    } while ( 0 )
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
    for ( int i = 0; i < N; i++, X += INCX, Y += INCY )
        *Y = *X;
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
    auto X = DX;
    auto y = std::abs( *X );
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
void Lapack<TYPE>::ger(
    int N, int M, TYPE alpha, const TYPE *x, int INCX, const TYPE *y, int INCY, TYPE *A, int LDA )
{
    NULL_USE( N );
    NULL_USE( M );
    NULL_USE( alpha );
    NULL_USE( x );
    NULL_USE( INCX );
    NULL_USE( y );
    NULL_USE( INCY );
    NULL_USE( A );
    NULL_USE( LDA );
    throw std::logic_error( "ger is not currently supported without blas/lapack" );
}
template<class TYPE>
void Lapack<TYPE>::gesv( int N, int NRHS, TYPE *A, int LDA, int *IPIV, TYPE *B, int LDB, int &INFO )
{
    NULL_USE( N );
    NULL_USE( NRHS );
    NULL_USE( A );
    NULL_USE( LDA );
    NULL_USE( IPIV );
    NULL_USE( B );
    NULL_USE( LDB );
    NULL_USE( INFO );
    throw std::logic_error( "gesv is not currently supported without blas/lapack" );
}
template<class TYPE>
void Lapack<TYPE>::gtsv( int N, int NRHS, TYPE *DL, TYPE *D, TYPE *DU, TYPE *B, int LDB, int &INFO )
{
    NULL_USE( N );
    NULL_USE( NRHS );
    NULL_USE( DL );
    NULL_USE( D );
    NULL_USE( DU );
    NULL_USE( B );
    NULL_USE( LDB );
    NULL_USE( INFO );
    throw std::logic_error( "gtsv is not currently supported without blas/lapack" );
}
template<class TYPE>
void Lapack<TYPE>::gbsv(
    int N, int KL, int KU, int NRHS, TYPE *AB, int LDAB, int *IPIV, TYPE *B, int LDB, int &INFO )
{
    NULL_USE( N );
    NULL_USE( KL );
    NULL_USE( KU );
    NULL_USE( NRHS );
    NULL_USE( AB );
    NULL_USE( LDAB );
    NULL_USE( IPIV );
    NULL_USE( B );
    NULL_USE( LDB );
    NULL_USE( INFO );
    throw std::logic_error( "gbsv is not currently supported without blas/lapack" );
}
template<class TYPE>
void Lapack<TYPE>::getrf( int M, int N, TYPE *A, int LDA, int *IPIV, int &INFO )
{
    NULL_USE( M );
    NULL_USE( N );
    NULL_USE( A );
    NULL_USE( LDA );
    NULL_USE( IPIV );
    NULL_USE( INFO );
    throw std::logic_error( "getrf is not currently supported without blas/lapack" );
}
template<class TYPE>
void Lapack<TYPE>::gttrf( int N, TYPE *DL, TYPE *D, TYPE *DU, TYPE *DU2, int *IPIV, int &INFO )
{
    NULL_USE( N );
    NULL_USE( DL );
    NULL_USE( D );
    NULL_USE( DU );
    NULL_USE( DU2 );
    NULL_USE( IPIV );
    NULL_USE( INFO );
    throw std::logic_error( "gttrf is not currently supported without blas/lapack" );
}
template<class TYPE>
void Lapack<TYPE>::gbtrf( int M, int N, int KL, int KU, TYPE *AB, int LDAB, int *IPIV, int &INFO )
{
    NULL_USE( M );
    NULL_USE( N );
    NULL_USE( KL );
    NULL_USE( KU );
    NULL_USE( AB );
    NULL_USE( LDAB );
    NULL_USE( IPIV );
    NULL_USE( INFO );
    throw std::logic_error( "gbtrf is not currently supported without blas/lapack" );
}
template<class TYPE>
void Lapack<TYPE>::getrs( char TRANS, int N, int NRHS, const TYPE *A, int LDA, const int *IPIV,
    TYPE *B, int LDB, int &INFO )
{
    NULL_USE( TRANS );
    NULL_USE( N );
    NULL_USE( NRHS );
    NULL_USE( A );
    NULL_USE( LDA );
    NULL_USE( IPIV );
    NULL_USE( B );
    NULL_USE( LDB );
    NULL_USE( INFO );
    throw std::logic_error( "getrs is not currently supported without blas/lapack" );
}
template<class TYPE>
void Lapack<TYPE>::gttrs( char TRANS, int N, int NRHS, const TYPE *DL, const TYPE *D,
    const TYPE *DU, const TYPE *DU2, const int *IPIV, TYPE *B, int LDB, int &INFO )
{
    NULL_USE( TRANS );
    NULL_USE( N );
    NULL_USE( NRHS );
    NULL_USE( DL );
    NULL_USE( D );
    NULL_USE( DU );
    NULL_USE( DU2 );
    NULL_USE( IPIV );
    NULL_USE( B );
    NULL_USE( LDB );
    NULL_USE( INFO );
    throw std::logic_error( "gttrs is not currently supported without blas/lapack" );
}
template<class TYPE>
void Lapack<TYPE>::gbtrs( char TRANS, int N, int KL, int KU, int NRHS, const TYPE *AB, int LDAB,
    const int *IPIV, TYPE *B, int LDB, int &INFO )
{
    NULL_USE( TRANS );
    NULL_USE( N );
    NULL_USE( KL );
    NULL_USE( KU );
    NULL_USE( NRHS );
    NULL_USE( AB );
    NULL_USE( LDAB );
    NULL_USE( IPIV );
    NULL_USE( B );
    NULL_USE( LDB );
    NULL_USE( INFO );
    throw std::logic_error( "gbtrs is not currently supported without blas/lapack" );
}
template<class TYPE>
void Lapack<TYPE>::getri(
    int N, TYPE *A, int LDA, const int *IPIV, TYPE *WORK, int LWORK, int &INFO )
{
    NULL_USE( N );
    NULL_USE( A );
    NULL_USE( LDA );
    NULL_USE( IPIV );
    NULL_USE( WORK );
    NULL_USE( LWORK );
    NULL_USE( INFO );
    throw std::logic_error( "getri is not currently supported without blas/lapack" );
}
template<class TYPE>
void Lapack<TYPE>::trsm( char SIDE, char UPLO, char TRANS, char DIAG, int M, int N, TYPE ALPHA,
    const TYPE *A, int LDA, TYPE *B, int LDB )
{
    NULL_USE( SIDE );
    NULL_USE( UPLO );
    NULL_USE( TRANS );
    NULL_USE( DIAG );
    NULL_USE( M );
    NULL_USE( N );
    NULL_USE( ALPHA );
    NULL_USE( A );
    NULL_USE( LDA );
    NULL_USE( B );
    NULL_USE( LDB );
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
