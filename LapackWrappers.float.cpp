#include "LapackWrappers.h"
#include "blas_lapack.h"

#include <mutex>
#include <stdexcept>


// Define macro to handle name mangling
#ifndef FORTRAN_WRAPPER
#if defined( USE_ACML )
#define FORTRAN_WRAPPER( x ) x##_
#elif defined( _WIN32 ) || defined( __hpux ) || defined( USE_MKL )
#define FORTRAN_WRAPPER( x ) x
#elif defined( USE_VECLIB ) || defined( USE_OPENBLAS ) || defined( USE_CBLAS )
#define FORTRAN_WRAPPER( x ) x##_
inline CBLAS_SIDE SIDE2( char SIDE )
{
    return ( SIDE = 'L' || SIDE == 'l' ) ? CblasLeft : CblasRight;
}
inline CBLAS_UPLO UPLO2( char UPLO )
{
    return ( UPLO = 'U' || UPLO == 'u' ) ? CblasUpper : CblasLower;
}
inline CBLAS_DIAG DIAG2( char DIAG )
{
    return ( DIAG = 'U' || DIAG == 'u' ) ? CblasUnit : CblasNonUnit;
}
inline CBLAS_TRANSPOSE TRANS2( char TRANS )
{
    CBLAS_TRANSPOSE ans = CblasNoTrans;
    if ( TRANS == 'N' || TRANS == 'n' ) {
        ans = CblasNoTrans;
    } else if ( TRANS == 'T' || TRANS == 't' ) {
        ans = CblasTrans;
    } else if ( TRANS == 'C' || TRANS == 'c' ) {
        ans = CblasConjTrans;
    }
    return ans;
}
#else
#define FORTRAN_WRAPPER( x ) x##_
#endif
#endif


// Lock for thread safety
[[maybe_unused]] static std::mutex d_mutex;


// lamch
#undef slamch
template<>
float Lapack<float>::lamch( char cmach )
{
#ifdef USE_ATLAS
    return clapack_slamch( cmach );
#elif defined( USE_ACML )
    return ::slamch( cmach );
#elif defined( USE_VECLIB )
    return FORTRAN_WRAPPER( ::slamch )( &cmach );
#elif defined( USE_OPENBLAS )
    return FORTRAN_WRAPPER( ::slamch )( &cmach, 1 );
#elif defined( USE_LAPACKE )
    return ::LAPACKE_slamch( cmach );
#else
    return FORTRAN_WRAPPER( ::slamch )( &cmach );
#endif
}


// Define the member functions
#undef scopy
template<>
void Lapack<float>::copy( int N, const float *DX, int INCX, float *DY, int INCY )
{
#ifdef USE_ATLAS
    cblas_scopy( N, DX, INCX, DY, INCY );
#elif defined( USE_VECLIB )
    cblas_scopy( N, DX, INCX, DY, INCY );
#elif defined( USE_OPENBLAS )
    cblas_scopy( N, DX, INCX, DY, INCY );
#elif defined( USE_CBLAS )
    cblas_scopy( N, DX, INCX, DY, INCY );
#elif defined( USE_MATLAB_LAPACK )
    ptrdiff_t Nl = N, INCXl = INCX, INCYl = INCY;
    FORTRAN_WRAPPER( ::scopy )( &Nl, (float *) DX, &INCXl, DY, &INCYl );
#else
    FORTRAN_WRAPPER( ::scopy )( &N, (float *) DX, &INCX, DY, &INCY );
#endif
}
// Define the member functions
#undef sswap
template<>
void Lapack<float>::swap( int N, float *DX, int INCX, float *DY, int INCY )
{
#ifdef USE_ATLAS
    cblas_sswap( N, DX, INCX, DY, INCY );
#elif defined( USE_VECLIB )
    cblas_sswap( N, DX, INCX, DY, INCY );
#elif defined( USE_OPENBLAS )
    cblas_sswap( N, DX, INCX, DY, INCY );
#elif defined( USE_CBLAS )
    cblas_sswap( N, DX, INCX, DY, INCY );
#elif defined( USE_MATLAB_LAPACK )
    ptrdiff_t Nl = N, INCXl = INCX, INCYl = INCY;
    FORTRAN_WRAPPER( ::sswap )( &Nl, DX, &INCXl, DY, &INCYl );
#else
    FORTRAN_WRAPPER( ::sswap )( &N, DX, &INCX, DY, &INCY );
#endif
}
#undef sscal
template<>
void Lapack<float>::scal( int N, float DA, float *DX, int INCX )
{
#ifdef USE_ATLAS
    cblas_sscal( N, DA, DX, INCX );
#elif defined( USE_VECLIB )
    cblas_sscal( N, DA, DX, INCX );
#elif defined( USE_OPENBLAS )
    cblas_sscal( N, DA, DX, INCX );
#elif defined( USE_CBLAS )
    cblas_sscal( N, DA, DX, INCX );
#elif defined( USE_MATLAB_LAPACK )
    ptrdiff_t Np = N, INCXp = INCX;
    FORTRAN_WRAPPER( ::sscal )( &Np, &DA, DX, &INCXp );
#else
    FORTRAN_WRAPPER( ::sscal )( &N, &DA, DX, &INCX );
#endif
}
#undef snrm2
template<>
double Lapack<float>::nrm2( int N, const float *DX, int INCX )
{
#ifdef USE_ATLAS
    return cblas_snrm2( N, DX, INCX );
#elif defined( USE_VECLIB )
    return cblas_snrm2( N, DX, INCX );
#elif defined( USE_OPENBLAS )
    return cblas_snrm2( N, DX, INCX );
#elif defined( USE_CBLAS )
    return cblas_snrm2( N, DX, INCX );
#elif defined( USE_MATLAB_LAPACK )
    ptrdiff_t Np = N, INCXp = INCX;
    return FORTRAN_WRAPPER( ::snrm2 )( &Np, (float *) DX, &INCXp );
#else
    return FORTRAN_WRAPPER( ::snrm2 )( &N, (float *) DX, &INCX );
#endif
}
#undef isamax
template<>
int Lapack<float>::iamax( int N, const float *DX, int INCX )
{
#ifdef USE_ATLAS
    return cblas_isamax( N, DX, INCX ) - 1;
#elif defined( USE_VECLIB )
    return cblas_isamax( N, DX, INCX ) - 1;
#elif defined( USE_OPENBLAS )
    return cblas_isamax( N, DX, INCX ) - 1;
#elif defined( USE_CBLAS )
    return cblas_isamax( N, DX, INCX ) - 1;
#elif defined( USE_MATLAB_LAPACK )
    ptrdiff_t Np = N, INCXp = INCX;
    return FORTRAN_WRAPPER( ::isamax )( &Np, (float *) DX, &INCXp ) - 1;
#else
    return FORTRAN_WRAPPER( ::isamax )( &N, (float *) DX, &INCX ) - 1;
#endif
}
#undef saxpy
template<>
void Lapack<float>::axpy( int N, float DA, const float *DX, int INCX, float *DY, int INCY )
{
#ifdef USE_ATLAS
    cblas_saxpy( N, DA, DX, INCX, DY, INCY );
#elif defined( USE_VECLIB )
    cblas_saxpy( N, DA, DX, INCX, DY, INCY );
#elif defined( USE_OPENBLAS )
    cblas_saxpy( N, DA, DX, INCX, DY, INCY );
#elif defined( USE_CBLAS )
    cblas_saxpy( N, DA, DX, INCX, DY, INCY );
#elif defined( USE_MATLAB_LAPACK )
    ptrdiff_t Np = N, INCXp = INCX, INCYp = INCY;
    FORTRAN_WRAPPER( ::saxpy )( &Np, &DA, (float *) DX, &INCXp, DY, &INCYp );
#else
    FORTRAN_WRAPPER( ::saxpy )( &N, &DA, (float *) DX, &INCX, DY, &INCY );
#endif
}
#undef sgemv
template<>
void Lapack<float>::gemv( char TRANS, int M, int N, float ALPHA, const float *A, int LDA,
    const float *DX, int INCX, float BETA, float *DY, int INCY )
{
#ifdef USE_ATLAS
    cblas_sgemv(
        CblasColMajor, (CBLAS_TRANSPOSE) TRANS, M, N, ALPHA, A, LDA, DX, INCX, BETA, DY, INCY );
#elif defined( USE_ACML )
    // FORTRAN_WRAPPER(::sgemv)(&TRANS,&M,&N,&ALPHA,(float*)A,&LDA,(float*)DX,&INCX,&BETA,DY,&INCY,1);
    ::sgemv( TRANS, M, N, ALPHA, (float *) A, LDA, (float *) DX, INCX, BETA, DY, INCY );
#elif defined( USE_VECLIB )
    cblas_sgemv( CblasColMajor, TRANS2( TRANS ), M, N, ALPHA, A, LDA, DX, INCX, BETA, DY, INCY );
#elif defined( USE_OPENBLAS )
    cblas_sgemv( CblasColMajor, TRANS2( TRANS ), M, N, ALPHA, A, LDA, DX, INCX, BETA, DY, INCY );
#elif defined( USE_CBLAS )
    cblas_sgemv( CblasColMajor, TRANS2( TRANS ), M, N, ALPHA, A, LDA, DX, INCX, BETA, DY, INCY );
#elif defined( USE_MATLAB_LAPACK )
    ptrdiff_t Mp = M, Np = N, LDAp = LDA, INCXp = INCX, INCYp = INCY;
    FORTRAN_WRAPPER( ::sgemv )
    ( &TRANS, &Mp, &Np, &ALPHA, (float *) A, &LDAp, (float *) DX, &INCXp, &BETA, DY, &INCYp );
#else
    FORTRAN_WRAPPER( ::sgemv )
    ( &TRANS, &M, &N, &ALPHA, (float *) A, &LDA, (float *) DX, &INCX, &BETA, DY, &INCY );
#endif
}
#undef sgemm
template<>
void Lapack<float>::gemm( char TRANSA, char TRANSB, int M, int N, int K, float ALPHA,
    const float *A, int LDA, const float *B, int LDB, float BETA, float *C, int LDC )
{
#ifdef USE_ATLAS
    cblas_sgemm( CblasColMajor, (CBLAS_TRANSPOSE) TRANSA, (CBLAS_TRANSPOSE) TRANSB, M, N, K, ALPHA,
        A, LDA, B, LDB, BETA, C, LDC );
#elif defined( USE_ACML )
    FORTRAN_WRAPPER( ::sgemm )
    ( &TRANSA, &TRANSB, &M, &N, &K, &ALPHA, (float *) A, &LDA, (float *) B, &LDB, &BETA, C, &LDC, 1,
        1 );
//::sgemm(TRANSA,TRANSA,M,N,K,ALPHA,(float*)A,LDA,(float*)B,LDB,BETA,C,LDC);
#elif defined( USE_VECLIB )
    cblas_sgemm( CblasColMajor, TRANS2( TRANSA ), TRANS2( TRANSB ), M, N, K, ALPHA, A, LDA, B, LDB,
        BETA, C, LDC );
#elif defined( USE_OPENBLAS )
    cblas_sgemm( CblasColMajor, TRANS2( TRANSA ), TRANS2( TRANSB ), M, N, K, ALPHA, A, LDA, B, LDB,
        BETA, C, LDC );
#elif defined( USE_CBLAS )
    cblas_sgemm( CblasColMajor, TRANS2( TRANSA ), TRANS2( TRANSB ), M, N, K, ALPHA, A, LDA, B, LDB,
        BETA, C, LDC );
#elif defined( USE_MATLAB_LAPACK )
    ptrdiff_t Mp = M, Np = N, Kp = K, LDAp = LDA, LDBp = LDB, LDCp = LDC;
    FORTRAN_WRAPPER( ::sgemm )
    ( &TRANSA, &TRANSB, &Mp, &Np, &Kp, &ALPHA, (float *) A, &LDAp, (float *) B, &LDBp, &BETA, C,
        &LDCp );
#else
    FORTRAN_WRAPPER( ::sgemm )
    ( &TRANSA, &TRANSB, &M, &N, &K, &ALPHA, (float *) A, &LDA, (float *) B, &LDB, &BETA, C, &LDC );
#endif
}
#undef sasum
template<>
double Lapack<float>::asum( int N, const float *DX, int INCX )
{
#ifdef USE_ATLAS
    return cblas_sasum( N, DX, INCX );
#elif defined( USE_VECLIB )
    return cblas_sasum( N, DX, INCX );
#elif defined( USE_OPENBLAS )
    return cblas_sasum( N, DX, INCX );
#elif defined( USE_CBLAS )
    return cblas_sasum( N, DX, INCX );
#elif defined( USE_MATLAB_LAPACK )
    ptrdiff_t Np = N, INCXp = INCX;
    return FORTRAN_WRAPPER( ::sasum )( &Np, (float *) DX, &INCXp );
#else
    return FORTRAN_WRAPPER( ::sasum )( &N, (float *) DX, &INCX );
#endif
}
#undef sdot
template<>
float Lapack<float>::dot( int N, const float *DX, int INCX, const float *DY, int INCY )
{
#ifdef USE_ATLAS
    return cblas_sdot( N, DX, INCX, DY, INCY );
#elif defined( USE_VECLIB )
    return cblas_sdot( N, DX, INCX, DY, INCY );
#elif defined( USE_OPENBLAS )
    return cblas_sdot( N, DX, INCX, DY, INCY );
#elif defined( USE_CBLAS )
    return cblas_sdot( N, DX, INCX, DY, INCY );
#elif defined( USE_MATLAB_LAPACK )
    ptrdiff_t Np = N, INCXp = INCX, INCYp = INCY;
    return FORTRAN_WRAPPER( ::sdot )( &Np, (float *) DX, &INCXp, (float *) DY, &INCYp );
#else
    return FORTRAN_WRAPPER( ::sdot )( &N, (float *) DX, &INCX, (float *) DY, &INCY );
#endif
}
#undef sger
template<>
void Lapack<float>::ger( int N, int M, float alpha, const float *x, int INCX, const float *y,
    int INCY, float *A, int LDA )
{
#ifdef USE_ATLAS
    cblas_sger( N, M, alpha, x, INCX, y, INCY, A, LDA );
#elif defined( USE_VECLIB )
    cblas_sger( CblasColMajor, N, M, alpha, x, INCX, y, INCY, A, LDA );
#elif defined( USE_OPENBLAS )
    cblas_sger( CblasColMajor, N, M, alpha, x, INCX, y, INCY, A, LDA );
#elif defined( USE_CBLAS )
    cblas_sger( CblasColMajor, N, M, alpha, x, INCX, y, INCY, A, LDA );
#elif defined( USE_MATLAB_LAPACK )
    ptrdiff_t Np = N, Mp = M, INCXp = INCX, INCYp = INCY, LDAp = LDA;
    FORTRAN_WRAPPER( ::sger )
    ( &Np, &Mp, &alpha, (float *) x, &INCXp, (float *) y, &INCYp, A, &LDAp );
#else
    FORTRAN_WRAPPER( ::sger )( &N, &M, &alpha, (float *) x, &INCX, (float *) y, &INCY, A, &LDA );
#endif
}
#undef sgesv
template<>
void Lapack<float>::gesv(
    int N, int NRHS, float *A, int LDA, int *IPIV, float *B, int LDB, int &INFO )
{
#ifdef USE_ATLAS
    INFO = clapack_sgesv( CblasColMajor, N, NRHS, A, LDA, IPIV, B, LDB );
#elif defined( USE_VECLIB )
    sgesv_( &N, &NRHS, A, &LDA, IPIV, B, &LDB, &INFO );
#elif defined( USE_OPENBLAS )
    sgesv_( &N, &NRHS, A, &LDA, IPIV, B, &LDB, &INFO );
#elif defined( USE_MATLAB_LAPACK )
    ptrdiff_t Np = N, NRHSp = NRHS, LDAp = LDA, LDBp = LDB, INFOp;
    ptrdiff_t *IPIVp = new ptrdiff_t[N];
    FORTRAN_WRAPPER( ::sgesv )( &Np, &NRHSp, A, &LDAp, IPIVp, B, &LDBp, &INFOp );
    for ( int i = 0; i < N; i++ ) {
        IPIV[i] = static_cast<int>( IPIVp[i] );
    }
    delete[] IPIVp;
    INFO = static_cast<int>( INFOp );
#else
    FORTRAN_WRAPPER( ::sgesv )( &N, &NRHS, (float *) A, &LDA, IPIV, B, &LDB, &INFO );
#endif
}
#undef sgtsv
template<>
void Lapack<float>::gtsv(
    int N, int NRHS, float *DL, float *D, float *DU, float *B, int LDB, int &INFO )
{
#ifdef USE_ATLAS
    throw std::logic_error( "ATLAS does not support sgtsv" );
#elif defined( USE_VECLIB )
    FORTRAN_WRAPPER( ::sgtsv )( &N, &NRHS, DL, D, DU, B, &LDB, &INFO );
#elif defined( USE_OPENBLAS )
    FORTRAN_WRAPPER( ::sgtsv )( &N, &NRHS, DL, D, DU, B, &LDB, &INFO );
#elif defined( USE_MATLAB_LAPACK )
    ptrdiff_t N2 = N, NRHS2 = NRHS, LDB2 = LDB, INFOp;
    FORTRAN_WRAPPER( ::sgtsv )( &N2, &NRHS2, DL, D, DU, B, &LDB2, &INFOp );
    INFO = static_cast<int>( INFOp );
#else
    FORTRAN_WRAPPER( ::sgtsv )( &N, &NRHS, DL, D, DU, B, &LDB, &INFO );
#endif
}
#undef sgbsv
template<>
void Lapack<float>::gbsv(
    int N, int KL, int KU, int NRHS, float *AB, int LDAB, int *IPIV, float *B, int LDB, int &INFO )
{
#ifdef USE_ATLAS
    throw std::logic_error( "ATLAS does not support sgbsv" );
#elif defined( USE_VECLIB )
    FORTRAN_WRAPPER( ::sgbsv )( &N, &KL, &KU, &NRHS, AB, &LDAB, IPIV, B, &LDB, &INFO );
#elif defined( USE_OPENBLAS )
    FORTRAN_WRAPPER( ::sgbsv )( &N, &KL, &KU, &NRHS, AB, &LDAB, IPIV, B, &LDB, &INFO );
#elif defined( USE_MATLAB_LAPACK )
    ptrdiff_t Np = N, KLp = KL, KUp = KU, NRHSp = NRHS, LDABp = LDAB, LDBp = LDB, INFOp;
    ptrdiff_t *IPIVp = new ptrdiff_t[N];
    FORTRAN_WRAPPER( ::sgbsv )( &Np, &KLp, &KUp, &NRHSp, AB, &LDABp, IPIVp, B, &LDBp, &INFOp );
    for ( int i = 0; i < N; i++ ) {
        IPIV[i] = static_cast<int>( IPIVp[i] );
    }
    delete[] IPIVp;
    INFO = static_cast<int>( INFOp );
#else
    FORTRAN_WRAPPER( ::sgbsv )( &N, &KL, &KU, &NRHS, AB, &LDAB, IPIV, B, &LDB, &INFO );
#endif
}
#undef sgetrf
template<>
void Lapack<float>::getrf( int M, int N, float *A, int LDA, int *IPIV, int &INFO )
{
#ifdef USE_ATLAS
    INFO = clapack_sgetrf( CblasColMajor, M, N, A, LDA, IPIV );
#elif defined( USE_VECLIB )
    FORTRAN_WRAPPER( ::sgetrf )( &M, &N, A, &LDA, IPIV, &INFO );
#elif defined( USE_OPENBLAS )
    FORTRAN_WRAPPER( ::sgetrf )( &M, &N, A, &LDA, IPIV, &INFO );
#elif defined( USE_MATLAB_LAPACK )
    ptrdiff_t Np = N, Mp = M, LDAp = LDA, INFOp;
    ptrdiff_t *IPIVp = new ptrdiff_t[N];
    FORTRAN_WRAPPER( ::sgetrf )( &Mp, &Np, A, &LDAp, IPIVp, &INFOp );
    for ( int i = 0; i < N; i++ ) {
        IPIV[i] = static_cast<int>( IPIVp[i] );
    }
    delete[] IPIVp;
    INFO = static_cast<int>( INFOp );
#else
    FORTRAN_WRAPPER( ::sgetrf )( &M, &N, A, &LDA, IPIV, &INFO );
#endif
}
#undef sgttrf
template<>
void Lapack<float>::gttrf( int N, float *DL, float *D, float *DU, float *DU2, int *IPIV, int &INFO )
{
#ifdef USE_ATLAS
    throw std::logic_error( "ATLAS does not support sgttrf" );
#elif defined( USE_VECLIB )
    FORTRAN_WRAPPER( ::sgttrf )( &N, DL, D, DU, DU2, IPIV, &INFO );
#elif defined( USE_OPENBLAS )
    FORTRAN_WRAPPER( ::sgttrf )( &N, DL, D, DU, DU2, IPIV, &INFO );
#elif defined( USE_MATLAB_LAPACK )
    ptrdiff_t Np     = N, INFOp;
    ptrdiff_t *IPIVp = new ptrdiff_t[N];
    FORTRAN_WRAPPER( ::sgttrf )( &Np, DL, D, DU, DU2, IPIVp, &INFOp );
    for ( int i = 0; i < N; i++ ) {
        IPIV[i] = static_cast<int>( IPIVp[i] );
    }
    delete[] IPIVp;
    INFO = static_cast<int>( INFOp );
#else
    FORTRAN_WRAPPER( ::sgttrf )( &N, DL, D, DU, DU2, IPIV, &INFO );
#endif
}
#undef sgbtrf
template<>
void Lapack<float>::gbtrf( int M, int N, int KL, int KU, float *AB, int LDAB, int *IPIV, int &INFO )
{
#ifdef USE_ATLAS
    throw std::logic_error( "ATLAS does not support sgbtrf" );
#elif defined( USE_VECLIB )
    FORTRAN_WRAPPER( ::sgbtrf )( &M, &N, &KL, &KU, AB, &LDAB, IPIV, &INFO );
#elif defined( USE_OPENBLAS )
    FORTRAN_WRAPPER( ::sgbtrf )( &M, &N, &KL, &KU, AB, &LDAB, IPIV, &INFO );
#elif defined( USE_MATLAB_LAPACK )
    ptrdiff_t Mp = M, Np = N, KLp = KL, KUp = KU, LDABp = LDAB, INFOp;
    ptrdiff_t *IPIVp = new ptrdiff_t[N];
    FORTRAN_WRAPPER( ::sgbtrf )( &Mp, &Np, &KLp, &KUp, AB, &LDABp, IPIVp, &INFOp );
    for ( int i = 0; i < N; i++ ) {
        IPIV[i] = static_cast<int>( IPIVp[i] );
    }
    delete[] IPIVp;
    INFO = static_cast<int>( INFOp );
#elif defined( USE_ACML )
    d_mutex.lock();
    FORTRAN_WRAPPER( ::sgbtrf )( &M, &N, &KL, &KU, AB, &LDAB, IPIV, &INFO );
    d_mutex.unlock();
#else
    FORTRAN_WRAPPER( ::sgbtrf )( &M, &N, &KL, &KU, AB, &LDAB, IPIV, &INFO );
#endif
}
#undef sgetrs
template<>
void Lapack<float>::getrs( char TRANS, int N, int NRHS, const float *A, int LDA, const int *IPIV,
    float *B, int LDB, int &INFO )
{
#ifdef USE_ATLAS
    INFO = clapack_sgetrs( CblasColMajor, (CBLAS_TRANSPOSE) TRANS, N, NRHS, A, LDA, IPIV, B, LDB );
#elif defined( USE_ACML )
    ::sgetrs( TRANS, N, NRHS, (float *) A, LDA, (int *) IPIV, B, LDB, &INFO );
#elif defined( USE_VECLIB )
    FORTRAN_WRAPPER( ::sgetrs )
    ( &TRANS, &N, &NRHS, (float *) A, &LDA, (int *) IPIV, B, &LDB, &INFO );
#elif defined( USE_OPENBLAS )
    FORTRAN_WRAPPER( ::sgetrs )
    ( &TRANS, &N, &NRHS, (float *) A, &LDA, (int *) IPIV, B, &LDB, &INFO, 1 );
#elif defined( USE_LAPACKE )
    FORTRAN_WRAPPER( ::sgetrs )
    ( &TRANS, &N, &NRHS, (float *) A, &LDA, (int *) IPIV, B, &LDB, &INFO );
#elif defined( USE_MATLAB_LAPACK )
    ptrdiff_t Np = N, NRHSp = NRHS, LDAp = LDA, LDBp = LDB, INFOp;
    ptrdiff_t *IPIVp = new ptrdiff_t[N];
    for ( int i = 0; i < N; i++ ) {
        IPIVp[i] = IPIV[i];
    }
    FORTRAN_WRAPPER( ::sgetrs )( &TRANS, &Np, &NRHSp, (float *) A, &LDAp, IPIVp, B, &LDBp, &INFOp );
    delete[] IPIVp;
    INFO = static_cast<int>( INFOp );
#else
    FORTRAN_WRAPPER( ::sgetrs )
    ( &TRANS, &N, &NRHS, (float *) A, &LDA, (int *) IPIV, B, &LDB, &INFO );
#endif
}
#undef sgttrs
template<>
void Lapack<float>::gttrs( char TRANS, int N, int NRHS, const float *DL, const float *D,
    const float *DU, const float *DU2, const int *IPIV, float *B, int LDB, int &INFO )
{
#ifdef USE_ATLAS
    throw std::logic_error( "ATLAS does not support sgttrs" );
#elif defined( USE_ACML )
    ::sgttrs( TRANS, N, NRHS, (float *) DL, (float *) D, (float *) DU, (float *) DU2, (int *) IPIV,
        B, LDB, &INFO );
#elif defined( USE_VECLIB )
    FORTRAN_WRAPPER( ::sgttrs )
    ( &TRANS, &N, &NRHS, (float *) DL, (float *) D, (float *) DU, (float *) DU2, (int *) IPIV, B,
        &LDB, &INFO );
#elif defined( USE_OPENBLAS )
    FORTRAN_WRAPPER( ::sgttrs )
    ( &TRANS, &N, &NRHS, (float *) DL, (float *) D, (float *) DU, (float *) DU2, (int *) IPIV, B,
        &LDB, &INFO, 1 );
#elif defined( USE_MATLAB_LAPACK )
    ptrdiff_t Np = N, NRHSp = NRHS, LDBp = LDB, INFOp;
    ptrdiff_t *IPIVp = new ptrdiff_t[N];
    for ( int i = 0; i < N; i++ ) {
        IPIVp[i] = IPIV[i];
    }
    FORTRAN_WRAPPER( ::sgttrs )
    ( &TRANS, &Np, &NRHSp, (float *) DL, (float *) D, (float *) DU, (float *) DU2, IPIVp, B, &LDBp,
        &INFOp );
    delete[] IPIVp;
    INFO = static_cast<int>( INFOp );
#else
    FORTRAN_WRAPPER( ::sgttrs )
    ( &TRANS, &N, &NRHS, (float *) DL, (float *) D, (float *) DU, (float *) DU2, (int *) IPIV, B,
        &LDB, &INFO );
#endif
}
#undef sgbtrs
template<>
void Lapack<float>::gbtrs( char TRANS, int N, int KL, int KU, int NRHS, const float *AB, int LDAB,
    const int *IPIV, float *B, int LDB, int &INFO )
{
#ifdef USE_ATLAS
    throw std::logic_error( "ATLAS does not support sgbtrs" );
#elif defined( USE_ACML )
    ::sgbtrs( TRANS, N, KL, KU, NRHS, (float *) AB, LDAB, (int *) IPIV, B, LDB, &INFO );
#elif defined( USE_VECLIB )
    FORTRAN_WRAPPER( ::sgbtrs )
    ( &TRANS, &N, &KL, &KU, &NRHS, (float *) AB, &LDAB, (int *) IPIV, B, &LDB, &INFO );
#elif defined( USE_OPENBLAS )
    FORTRAN_WRAPPER( ::sgbtrs )
    ( &TRANS, &N, &KL, &KU, &NRHS, (float *) AB, &LDAB, (int *) IPIV, B, &LDB, &INFO, 1 );
#elif defined( USE_MATLAB_LAPACK )
    ptrdiff_t Np = N, KLp = KL, KUp = KU, NRHSp = NRHS, LDABp = LDAB, LDBp = LDB, INFOp;
    ptrdiff_t *IPIVp = new ptrdiff_t[N];
    for ( int i = 0; i < N; i++ ) {
        IPIVp[i] = IPIV[i];
    }
    FORTRAN_WRAPPER( ::sgbtrs )
    ( &TRANS, &Np, &KLp, &KUp, &NRHSp, (float *) AB, &LDABp, IPIVp, B, &LDBp, &INFOp );
    INFO = static_cast<int>( INFOp );
    delete[] IPIVp;
#else
    FORTRAN_WRAPPER( ::sgbtrs )
    ( &TRANS, &N, &KL, &KU, &NRHS, (float *) AB, &LDAB, (int *) IPIV, B, &LDB, &INFO );
#endif
}
#undef sgetri
template<>
void Lapack<float>::getri(
    int N, float *A, int LDA, const int *IPIV, float *WORK, int LWORK, int &INFO )
{
#ifdef USE_ATLAS
    INFO = clapack_sgetri( CblasColMajor, N, A, LDA, IPIV );
#elif defined( USE_ACML )
    ::sgetri_( &N, A, &LDA, (int *) IPIV, WORK, &LWORK, &INFO );
#elif defined( USE_VECLIB )
    FORTRAN_WRAPPER( ::sgetri )( &N, A, &LDA, (int *) IPIV, WORK, &LWORK, &INFO );
#elif defined( USE_OPENBLAS )
    FORTRAN_WRAPPER( ::sgetri )( &N, A, &LDA, (int *) IPIV, WORK, &LWORK, &INFO );
#elif defined( USE_MATLAB_LAPACK )
    ptrdiff_t Np = N, LDAp = LDA, LWORKp = LWORK, INFOp;
    ptrdiff_t *IPIVp = new ptrdiff_t[N];
    for ( int i = 0; i < N; i++ ) {
        IPIVp[i] = IPIV[i];
    }
    FORTRAN_WRAPPER( ::sgetri )( &Np, A, &LDAp, IPIVp, WORK, &LWORKp, &INFOp );
    INFO = static_cast<int>( INFOp );
    delete[] IPIVp;
#else
    FORTRAN_WRAPPER( ::sgetri )( &N, A, &LDA, (int *) IPIV, WORK, &LWORK, &INFO );
#endif
}
#undef strsm
template<>
void Lapack<float>::trsm( char SIDE, char UPLO, char TRANS, char DIAG, int M, int N, float ALPHA,
    const float *A, int LDA, float *B, int LDB )
{
#ifdef USE_ATLAS
    throw std::logic_error( "strsm not implemented for ATLAS" );
#elif defined( USE_ACML )
    char SIDE2[2] = { SIDE, 0 }, UPLO2[2] = { UPLO, 0 }, TRANS2[2] = { TRANS, 0 },
         DIAG2[2] = { DIAG, 0 };
    ::strsm_( SIDE2, UPLO2, TRANS2, DIAG2, &M, &N, &ALPHA, (float *) A, &LDA, B, &LDB, 1, 1, 1, 1 );
#elif defined( USE_VECLIB ) || defined( USE_OPENBLAS ) || defined( USE_CBLAS )
    cblas_strsm( CblasColMajor, SIDE2( SIDE ), UPLO2( UPLO ), TRANS2( TRANS ), DIAG2( DIAG ), M, N,
        ALPHA, (float *) A, LDA, B, LDB );
#elif defined( USE_MATLAB_LAPACK )
    ptrdiff_t Mp = M, Np = N, LDAp = LDA, LDBp = LDB;
    FORTRAN_WRAPPER( ::strsm )
    ( &SIDE, &UPLO, &TRANS, &DIAG, &Mp, &Np, &ALPHA, (float *) A, &LDAp, B, &LDBp );
#else
    FORTRAN_WRAPPER( ::strsm )
    ( &SIDE, &UPLO, &TRANS, &DIAG, &M, &N, &ALPHA, (float *) A, &LDA, B, &LDB );
#endif
}
