#include "LapackWrappers.h"
#include "blas_lapack.h"

#include <complex>
#include <cstring>
#include <mutex>
#include <stdexcept>


// Define complex types for lapack interface
#if defined( USE_MKL )
typedef MKL_Complex16 Complex;
static_assert( sizeof( MKL_Complex16 ) == sizeof( std::complex<double> ) );
#else
typedef double _Complex Complex;
#endif
static inline Complex convert( const std::complex<double> &x )
{
    static_assert( sizeof( std::complex<double> ) == sizeof( Complex ) );
    Complex y;
    memcpy( &y, &x, sizeof( Complex ) );
    return y;
}


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


// Helper functions
#if defined( USE_MATLAB_LAPACK )
template<class T1, class T2>
void convert( int N, const T1 *x, T2 *y )
{
    for ( int i = 0; i < N; i++ )
        y[i] = x[i];
}
#endif


// Define the member functions
#undef scopy
template<>
void Lapack<std::complex<double>>::copy(
    int N, const std::complex<double> *DX_, int INCX, std::complex<double> *DY_, int INCY )
{
    auto DX = reinterpret_cast<const Complex *>( DX_ );
    auto DY = reinterpret_cast<Complex *>( DY_ );
#ifdef USE_ATLAS
    cblas_zcopy( N, DX, INCX, DY, INCY );
#elif defined( USE_VECLIB )
    cblas_zcopy( N, DX, INCX, DY, INCY );
#elif defined( USE_OPENBLAS )
    cblas_zcopy( N, DX, INCX, DY, INCY );
#elif defined( USE_CBLAS )
    cblas_zcopy( N, DX, INCX, DY, INCY );
#elif defined( USE_MATLAB_LAPACK )
    ptrdiff_t Nl = N, INCXl = INCX, INCYl = INCY;
    FORTRAN_WRAPPER( ::zcopy )( &Nl, (double *) DX, &INCXl, (double *) DY, &INCYl );
#else
    FORTRAN_WRAPPER( ::zcopy )( &N, const_cast<Complex *>( DX ), &INCX, DY, &INCY );
#endif
}
// Define the member functions
#undef sswap
template<>
void Lapack<std::complex<double>>::swap(
    int N, std::complex<double> *DX_, int INCX, std::complex<double> *DY_, int INCY )
{
    auto DX = reinterpret_cast<Complex *>( DX_ );
    auto DY = reinterpret_cast<Complex *>( DY_ );
#ifdef USE_ATLAS
    cblas_zswap( N, DX, INCX, DY, INCY );
#elif defined( USE_VECLIB )
    cblas_zswap( N, DX, INCX, DY, INCY );
#elif defined( USE_OPENBLAS )
    cblas_zswap( N, DX, INCX, DY, INCY );
#elif defined( USE_CBLAS )
    cblas_zswap( N, DX, INCX, DY, INCY );
#elif defined( USE_MATLAB_LAPACK )
    ptrdiff_t Nl = N, INCXl = INCX, INCYl = INCY;
    FORTRAN_WRAPPER( ::zswap )( &Nl, (double *) DX, &INCXl, (double *) DY, &INCYl );
#else
    FORTRAN_WRAPPER( ::zswap )( &N, const_cast<Complex *>( DX ), &INCX, DY, &INCY );
#endif
}
#undef sscal
template<>
void Lapack<std::complex<double>>::scal(
    int N, std::complex<double> DA_, std::complex<double> *DX_, int INCX )
{
    auto DA = reinterpret_cast<Complex *>( &DA_ );
    auto DX = reinterpret_cast<Complex *>( DX_ );
#ifdef USE_ATLAS
    cblas_zscal( N, DA, DX, INCX );
#elif defined( USE_VECLIB )
    cblas_zscal( N, DA, DX, INCX );
#elif defined( USE_OPENBLAS )
    cblas_zscal( N, DA, DX, INCX );
#elif defined( USE_CBLAS )
    cblas_zscal( N, DA, DX, INCX );
#elif defined( USE_MATLAB_LAPACK )
    ptrdiff_t Np = N, INCXp = INCX;
    FORTRAN_WRAPPER( ::zscal )( &Np, (double *) DA, (double *) DX, &INCXp );
#else
    FORTRAN_WRAPPER( ::zscal )( &N, DA, DX, &INCX );
#endif
}
#undef snrm2
template<>
double Lapack<std::complex<double>>::nrm2( int N, const std::complex<double> *DX_, int INCX )
{
    auto DX = reinterpret_cast<const Complex *>( DX_ );
#ifdef USE_ATLAS
    return cblas_dznrm2( N, DX, INCX );
#elif defined( USE_VECLIB )
    return cblas_dznrm2( N, DX, INCX );
#elif defined( USE_OPENBLAS )
    return cblas_dznrm2( N, DX, INCX );
#elif defined( USE_CBLAS )
    return cblas_dznrm2( N, DX, INCX );
#elif defined( USE_MATLAB_LAPACK )
    ptrdiff_t Np = N, INCXp = INCX;
    return FORTRAN_WRAPPER( ::dznrm2 )( &Np, (double *) DX, &INCXp );
#else
    return FORTRAN_WRAPPER( ::dznrm2 )( &N, const_cast<Complex *>( DX ), &INCX );
#endif
}
#undef isamax
template<>
int Lapack<std::complex<double>>::iamax( int N, const std::complex<double> *DX_, int INCX )
{
    auto DX = reinterpret_cast<const Complex *>( DX_ );
#ifdef USE_ATLAS
    return cblas_izamax( N, DX, INCX ) - 1;
#elif defined( USE_VECLIB )
    return cblas_izamax( N, DX, INCX ) - 1;
#elif defined( USE_OPENBLAS )
    return cblas_izamax( N, DX, INCX ) - 1;
#elif defined( USE_CBLAS )
    return cblas_izamax( N, DX, INCX ) - 1;
#elif defined( USE_MATLAB_LAPACK )
    ptrdiff_t Np = N, INCXp = INCX;
    return FORTRAN_WRAPPER( ::izamax )( &Np, (double *) DX, &INCXp ) - 1;
#else
    return FORTRAN_WRAPPER( ::izamax )( &N, const_cast<Complex *>( DX ), &INCX ) - 1;
#endif
}
#undef saxpy
template<>
void Lapack<std::complex<double>>::axpy( int N, std::complex<double> DA_,
    const std::complex<double> *DX_, int INCX, std::complex<double> *DY_, int INCY )
{
    auto DA = convert( DA_ );
    auto DX = reinterpret_cast<const Complex *>( DX_ );
    auto DY = reinterpret_cast<Complex *>( DY_ );
#ifdef USE_ATLAS
    cblas_zaxpy( N, &DA, DX, INCX, DY, INCY );
#elif defined( USE_VECLIB )
    cblas_zaxpy( N, &DA, DX, INCX, DY, INCY );
#elif defined( USE_OPENBLAS )
    cblas_zaxpy( N, &DA, DX, INCX, DY, INCY );
#elif defined( USE_CBLAS )
    cblas_zaxpy( N, &DA, DX, INCX, DY, INCY );
#elif defined( USE_MATLAB_LAPACK )
    ptrdiff_t Np = N, INCXp = INCX, INCYp = INCY;
    FORTRAN_WRAPPER( ::zaxpy )( &Np, (double *) &DA, (double *) DX, &INCXp, (double *) DY, &INCYp );
#else
    FORTRAN_WRAPPER( ::zaxpy )( &N, &DA, const_cast<Complex *>( DX ), &INCX, DY, &INCY );
#endif
}
#undef sgemv
template<>
void Lapack<std::complex<double>>::gemv( char TRANS, int M, int N, std::complex<double> ALPHA_,
    const std::complex<double> *A_, int LDA, const std::complex<double> *DX_, int INCX,
    std::complex<double> BETA_, std::complex<double> *DY_, int INCY )
{
    auto ALPHA = convert( ALPHA_ );
    auto BETA  = convert( BETA_ );
    auto A     = reinterpret_cast<const Complex *>( A_ );
    auto DX    = reinterpret_cast<const Complex *>( DX_ );
    auto DY    = reinterpret_cast<Complex *>( DY_ );
#ifdef USE_ATLAS
    cblas_zgemv(
        CblasColMajor, (CBLAS_TRANSPOSE) TRANS, M, N, &ALPHA, A, LDA, DX, INCX, &BETA, DY, INCY );
#elif defined( USE_ACML )
    ::zgemv( TRANS, M, N, &ALPHA, const_cast<Complex *>( A ), LDA, const_cast<Complex *>( DX ),
        INCX, &BETA, DY, INCY );
#elif defined( USE_VECLIB )
    cblas_zgemv( CblasColMajor, TRANS2( TRANS ), M, N, &ALPHA, A, LDA, DX, INCX, &BETA, DY, INCY );
#elif defined( USE_OPENBLAS )
    cblas_zgemv( CblasColMajor, TRANS2( TRANS ), M, N, &ALPHA, A, LDA, DX, INCX, &BETA, DY, INCY );
#elif defined( USE_CBLAS )
    cblas_zgemv( CblasColMajor, TRANS2( TRANS ), M, N, &ALPHA, A, LDA, DX, INCX, &BETA, DY, INCY );
#elif defined( USE_MATLAB_LAPACK )
    ptrdiff_t Mp = M, Np = N, LDAp = LDA, INCXp = INCX, INCYp = INCY;
    FORTRAN_WRAPPER( ::zgemv )
    ( &TRANS, &Mp, &Np, (double *) &ALPHA, (double *) A, &LDAp, (double *) DX, &INCXp,
        (double *) &BETA, (double *) DY, &INCYp );
#else
    FORTRAN_WRAPPER( ::zgemv )
    ( &TRANS, &M, &N, &ALPHA, const_cast<Complex *>( A ), &LDA, const_cast<Complex *>( DX ), &INCX,
        &BETA, DY, &INCY );
#endif
}
#undef sgemm
template<>
void Lapack<std::complex<double>>::gemm( char TRANSA, char TRANSB, int M, int N, int K,
    std::complex<double> ALPHA_, const std::complex<double> *A_, int LDA,
    const std::complex<double> *B_, int LDB, std::complex<double> BETA_, std::complex<double> *C_,
    int LDC )
{
    auto ALPHA = convert( ALPHA_ );
    auto BETA  = convert( BETA_ );
    auto A     = reinterpret_cast<const Complex *>( A_ );
    auto B     = reinterpret_cast<const Complex *>( B_ );
    auto C     = reinterpret_cast<Complex *>( C_ );
#ifdef USE_ATLAS
    cblas_zgemm( CblasColMajor, (CBLAS_TRANSPOSE) TRANSA, (CBLAS_TRANSPOSE) TRANSB, M, N, K, &ALPHA,
        A, LDA, B, LDB, &BETA, C, LDC );
#elif defined( USE_ACML )
    FORTRAN_WRAPPER( ::zgemm )
    ( &TRANSA, &TRANSB, &M, &N, &K, &ALPHA, const_cast<Complex *>( A ), &LDA,
        const_cast<Complex *>( B ), &LDB, &BETA, C, &LDC, 1, 1 );
//::zgemm(TRANSA,TRANSA,M,N,K,ALPHA,(std::complex<double>*)A,LDA,(std::complex<double>*)B,LDB,BETA,C,LDC);
#elif defined( USE_VECLIB )
    cblas_zgemm( CblasColMajor, TRANS2( TRANSA ), TRANS2( TRANSB ), M, N, K, &ALPHA, A, LDA, B, LDB,
        &BETA, C, LDC );
#elif defined( USE_OPENBLAS )
    cblas_zgemm( CblasColMajor, TRANS2( TRANSA ), TRANS2( TRANSB ), M, N, K, &ALPHA, A, LDA, B, LDB,
        &BETA, C, LDC );
#elif defined( USE_CBLAS )
    cblas_zgemm( CblasColMajor, TRANS2( TRANSA ), TRANS2( TRANSB ), M, N, K, &ALPHA, A, LDA, B, LDB,
        &BETA, C, LDC );
#elif defined( USE_MATLAB_LAPACK )
    ptrdiff_t Mp = M, Np = N, Kp = K, LDAp = LDA, LDBp = LDB, LDCp = LDC;
    FORTRAN_WRAPPER( ::zgemm )
    ( &TRANSA, &TRANSB, &Mp, &Np, &Kp, (double *) &ALPHA, (double *) A, &LDAp, (double *) B, &LDBp,
        (double *) &BETA, (double *) C, &LDCp );
#else
    FORTRAN_WRAPPER( ::zgemm )
    ( &TRANSA, &TRANSB, &M, &N, &K, &ALPHA, const_cast<Complex *>( A ), &LDA,
        const_cast<Complex *>( B ), &LDB, &BETA, const_cast<Complex *>( C ), &LDC );
#endif
}
#undef sasum
template<>
double Lapack<std::complex<double>>::asum( int N, const std::complex<double> *DX_, int INCX )
{
    if ( N == 0 )
        return 0;
    auto DX = reinterpret_cast<const Complex *>( DX_ );
#ifdef USE_ATLAS
    return cblas_dzasum( N, DX, INCX );
#elif defined( USE_VECLIB )
    return cblas_dzasum( N, DX, INCX );
#elif defined( USE_OPENBLAS )
    return cblas_dzasum( N, DX, INCX );
#elif defined( USE_CBLAS )
    return cblas_dzasum( N, DX, INCX );
#elif defined( USE_MATLAB_LAPACK )
    ptrdiff_t Np = N, INCXp = INCX;
    return FORTRAN_WRAPPER( ::dzasum )( &Np, (double *) DX, &INCXp );
#else
    return FORTRAN_WRAPPER( ::dzasum )( &N, const_cast<Complex *>( DX ), &INCX );
#endif
}
#undef sdot
template<>
std::complex<double> Lapack<std::complex<double>>::dot(
    int N, const std::complex<double> *DX_, int INCX, const std::complex<double> *DY_, int INCY )
{
    auto DX = reinterpret_cast<const Complex *>( DX_ );
    auto DY = reinterpret_cast<const Complex *>( DY_ );
#ifdef USE_ATLAS
    return cblas_zdotc( N, DX, INCX, DY, INCY );
#elif defined( USE_VECLIB )
    return cblas_zdotc( N, DX, INCX, DY, INCY );
#elif defined( USE_OPENBLAS )
    return cblas_zdotc( N, DX, INCX, DY, INCY );
#elif defined( USE_CBLAS )
    return cblas_zdotc( N, DX, INCX, DY, INCY );
#elif defined( USE_MATLAB_LAPACK )
    ptrdiff_t Np = N, INCXp = INCX, INCYp = INCY;
    auto rtn = FORTRAN_WRAPPER( ::zdotc )( &Np, (double *) DX, &INCXp, (double *) DY, &INCYp );
    return std::complex<double>( rtn.r, rtn.i );
#elif defined( USE_MKL )
    MKL_Complex16 c;
    zdotc( &c, &N, const_cast<Complex *>( DX ), &INCX, const_cast<Complex *>( DY ), &INCY );
    return { c.real, c.imag };
#else
    return FORTRAN_WRAPPER( ::zdotc )(
        &N, const_cast<Complex *>( DX ), &INCX, const_cast<Complex *>( DY ), &INCY );
#endif
}
#undef sger
template<>
void Lapack<std::complex<double>>::ger( int N, int M, std::complex<double> ALPHA_,
    const std::complex<double> *x_, int INCX, const std::complex<double> *y_, int INCY,
    std::complex<double> *A_, int LDA )
{
    auto ALPHA = convert( ALPHA_ );
    auto x     = reinterpret_cast<const Complex *>( x_ );
    auto y     = reinterpret_cast<const Complex *>( y_ );
    auto A     = reinterpret_cast<Complex *>( A_ );
#ifdef USE_ATLAS
    cblas_zgerc( N, M, ALPHA, x, INCX, y, INCY, A, LDA );
#elif defined( USE_VECLIB )
    cblas_zgerc( CblasColMajor, N, M, &ALPHA, x, INCX, y, INCY, A, LDA );
#elif defined( USE_OPENBLAS )
    cblas_zgerc( CblasColMajor, N, M, &ALPHA, x, INCX, y, INCY, A, LDA );
#elif defined( USE_CBLAS )
    cblas_zgerc( CblasColMajor, N, M, &ALPHA, x, INCX, y, INCY, A, LDA );
#elif defined( USE_MATLAB_LAPACK )
    ptrdiff_t Np = N, Mp = M, INCXp = INCX, INCYp = INCY, LDAp = LDA;
    FORTRAN_WRAPPER( ::zgerc )
    ( &Np, &Mp, (double *) &ALPHA, (double *) x, &INCXp, (double *) y, &INCYp, (double *) A,
        &LDAp );
#else
    FORTRAN_WRAPPER( ::zgerc )
    ( &N, &M, &ALPHA, const_cast<Complex *>( x ), &INCX, const_cast<Complex *>( y ), &INCY,
        const_cast<Complex *>( A ), &LDA );
#endif
}
#undef sgesv
template<>
void Lapack<std::complex<double>>::gesv( int N, int NRHS, std::complex<double> *A_, int LDA,
    int *IPIV, std::complex<double> *B_, int LDB, int &INFO )
{
    auto A = reinterpret_cast<Complex *>( A_ );
    auto B = reinterpret_cast<Complex *>( B_ );
#ifdef USE_ATLAS
    INFO = clapack_zgesv( CblasColMajor, N, NRHS, A, LDA, IPIV, B, LDB );
#elif defined( USE_VECLIB )
    zgesv_( &N, &NRHS, A, &LDA, IPIV, B, &LDB, &INFO );
#elif defined( USE_OPENBLAS )
    zgesv_( &N, &NRHS, A, &LDA, IPIV, B, &LDB, &INFO );
#elif defined( USE_MATLAB_LAPACK )
    ptrdiff_t Np = N, NRHSp = NRHS, LDAp = LDA, LDBp = LDB, INFOp;
    ptrdiff_t *IPIVp = new ptrdiff_t[N];
    FORTRAN_WRAPPER( ::zgesv )
    ( &Np, &NRHSp, (double *) A, &LDAp, IPIVp, (double *) B, &LDBp, &INFOp );
    convert( N, IPIVp, IPIV );
    delete[] IPIVp;
    INFO = static_cast<int>( INFOp );
#else
    FORTRAN_WRAPPER( ::zgesv )
    ( &N, &NRHS, const_cast<Complex *>( A ), &LDA, IPIV, const_cast<Complex *>( B ), &LDB, &INFO );
#endif
}
#undef sgtsv
template<>
void Lapack<std::complex<double>>::gtsv( int N, int NRHS, std::complex<double> *DL_,
    std::complex<double> *D_, std::complex<double> *DU_, std::complex<double> *B_, int LDB,
    int &INFO )
{
    auto DL = reinterpret_cast<Complex *>( DL_ );
    auto D  = reinterpret_cast<Complex *>( D_ );
    auto DU = reinterpret_cast<Complex *>( DU_ );
    auto B  = reinterpret_cast<Complex *>( B_ );
#ifdef USE_ATLAS
    throw std::logic_error( "ATLAS does not support sgtsv" );
#elif defined( USE_VECLIB )
    FORTRAN_WRAPPER( ::zgtsv )( &N, &NRHS, DL, D, DU, B, &LDB, &INFO );
#elif defined( USE_OPENBLAS )
    FORTRAN_WRAPPER( ::zgtsv )( &N, &NRHS, DL, D, DU, B, &LDB, &INFO );
#elif defined( USE_MATLAB_LAPACK )
    ptrdiff_t N2 = N, NRHS2 = NRHS, LDB2 = LDB, INFOp;
    FORTRAN_WRAPPER( ::zgtsv )
    ( &N2, &NRHS2, (double *) DL, (double *) D, (double *) DU, (double *) B, &LDB2, &INFOp );
    INFO = static_cast<int>( INFOp );
#else
    FORTRAN_WRAPPER( ::zgtsv )
    ( &N, &NRHS, DL, D, DU, const_cast<Complex *>( B ), &LDB, &INFO );
#endif
}
#undef sgbsv
template<>
void Lapack<std::complex<double>>::gbsv( int N, int KL, int KU, int NRHS, std::complex<double> *AB_,
    int LDAB, int *IPIV, std::complex<double> *B_, int LDB, int &INFO )
{
    auto AB = reinterpret_cast<Complex *>( AB_ );
    auto B  = reinterpret_cast<Complex *>( B_ );
#ifdef USE_ATLAS
    throw std::logic_error( "ATLAS does not support sgbsv" );
#elif defined( USE_VECLIB )
    FORTRAN_WRAPPER( ::zgbsv )( &N, &KL, &KU, &NRHS, AB, &LDAB, IPIV, B, &LDB, &INFO );
#elif defined( USE_OPENBLAS )
    FORTRAN_WRAPPER( ::zgbsv )( &N, &KL, &KU, &NRHS, AB, &LDAB, IPIV, B, &LDB, &INFO );
#elif defined( USE_MATLAB_LAPACK )
    ptrdiff_t Np = N, KLp = KL, KUp = KU, NRHSp = NRHS, LDABp = LDAB, LDBp = LDB, INFOp;
    ptrdiff_t *IPIVp = new ptrdiff_t[N];
    FORTRAN_WRAPPER( ::zgbsv )
    ( &Np, &KLp, &KUp, &NRHSp, (double *) AB, &LDABp, IPIVp, (double *) B, &LDBp, &INFOp );
    convert( N, IPIVp, IPIV );
    delete[] IPIVp;
    INFO = static_cast<int>( INFOp );
#else
    FORTRAN_WRAPPER( ::zgbsv )
    ( &N, &KL, &KU, &NRHS, const_cast<Complex *>( AB ), &LDAB, IPIV, const_cast<Complex *>( B ),
        &LDB, &INFO );
#endif
}
#undef sgetrf
template<>
void Lapack<std::complex<double>>::getrf(
    int M, int N, std::complex<double> *A_, int LDA, int *IPIV, int &INFO )
{
    auto A = reinterpret_cast<Complex *>( A_ );
#ifdef USE_ATLAS
    INFO = clapack_sgetrf( CblasColMajor, M, N, A, LDA, IPIV );
#elif defined( USE_VECLIB )
    FORTRAN_WRAPPER( ::zgetrf )( &M, &N, A, &LDA, IPIV, &INFO );
#elif defined( USE_OPENBLAS )
    FORTRAN_WRAPPER( ::zgetrf )( &M, &N, A, &LDA, IPIV, &INFO );
#elif defined( USE_MATLAB_LAPACK )
    ptrdiff_t Np = N, Mp = M, LDAp = LDA, INFOp;
    ptrdiff_t *IPIVp = new ptrdiff_t[N];
    FORTRAN_WRAPPER( ::zgetrf )( &Mp, &Np, (double *) A, &LDAp, IPIVp, &INFOp );
    convert( N, IPIVp, IPIV );
    delete[] IPIVp;
    INFO = static_cast<int>( INFOp );
#else
    FORTRAN_WRAPPER( ::zgetrf )( &M, &N, A, &LDA, IPIV, &INFO );
#endif
}
#undef sgttrf
template<>
void Lapack<std::complex<double>>::gttrf( int N, std::complex<double> *DL_,
    std::complex<double> *D_, std::complex<double> *DU_, std::complex<double> *DU2_, int *IPIV,
    int &INFO )
{
    auto DL  = reinterpret_cast<Complex *>( DL_ );
    auto D   = reinterpret_cast<Complex *>( D_ );
    auto DU  = reinterpret_cast<Complex *>( DU_ );
    auto DU2 = reinterpret_cast<Complex *>( DU2_ );
#ifdef USE_ATLAS
    throw std::logic_error( "ATLAS does not support sgttrf" );
#elif defined( USE_VECLIB )
    FORTRAN_WRAPPER( ::zgttrf )( &N, DL, D, DU, DU2, IPIV, &INFO );
#elif defined( USE_OPENBLAS )
    FORTRAN_WRAPPER( ::zgttrf )( &N, DL, D, DU, DU2, IPIV, &INFO );
#elif defined( USE_MATLAB_LAPACK )
    ptrdiff_t Np     = N, INFOp;
    ptrdiff_t *IPIVp = new ptrdiff_t[N];
    FORTRAN_WRAPPER( ::zgttrf )
    ( &Np, (double *) DL, (double *) D, (double *) DU, (double *) DU2, IPIVp, &INFOp );
    convert( N, IPIVp, IPIV );
    delete[] IPIVp;
    INFO = static_cast<int>( INFOp );
#else
    FORTRAN_WRAPPER( ::zgttrf )( &N, DL, D, DU, DU2, IPIV, &INFO );
#endif
}
#undef sgbtrf
template<>
void Lapack<std::complex<double>>::gbtrf(
    int M, int N, int KL, int KU, std::complex<double> *AB_, int LDAB, int *IPIV, int &INFO )
{
    auto AB = reinterpret_cast<Complex *>( AB_ );
#ifdef USE_ATLAS
    throw std::logic_error( "ATLAS does not support sgbtrf" );
#elif defined( USE_VECLIB )
    FORTRAN_WRAPPER( ::zgbtrf )( &M, &N, &KL, &KU, AB, &LDAB, IPIV, &INFO );
#elif defined( USE_OPENBLAS )
    FORTRAN_WRAPPER( ::zgbtrf )( &M, &N, &KL, &KU, AB, &LDAB, IPIV, &INFO );
#elif defined( USE_MATLAB_LAPACK )
    ptrdiff_t Mp = M, Np = N, KLp = KL, KUp = KU, LDABp = LDAB, INFOp;
    ptrdiff_t *IPIVp = new ptrdiff_t[N];
    FORTRAN_WRAPPER( ::zgbtrf )( &Mp, &Np, &KLp, &KUp, (double *) AB, &LDABp, IPIVp, &INFOp );
    convert( N, IPIVp, IPIV );
    delete[] IPIVp;
    INFO = static_cast<int>( INFOp );
#elif defined( USE_ACML )
    d_mutex.lock();
    FORTRAN_WRAPPER( ::zgbtrf )( &M, &N, &KL, &KU, AB, &LDAB, IPIV, &INFO );
    d_mutex.unlock();
#else
    FORTRAN_WRAPPER( ::zgbtrf )( &M, &N, &KL, &KU, AB, &LDAB, IPIV, &INFO );
#endif
}
#undef sgetrs
template<>
void Lapack<std::complex<double>>::getrs( char TRANS, int N, int NRHS,
    const std::complex<double> *A_, int LDA, const int *IPIV, std::complex<double> *B_, int LDB,
    int &INFO )
{
    auto A = reinterpret_cast<const Complex *>( A_ );
    auto B = reinterpret_cast<Complex *>( B_ );
#ifdef USE_ATLAS
    INFO = clapack_sgetrs( CblasColMajor, (CBLAS_TRANSPOSE) TRANS, N, NRHS, A, LDA, IPIV, B, LDB );
#elif defined( USE_ACML )
    ::zgetrs( TRANS, N, NRHS, const_cast<Complex *>( A ), LDA, (int *) IPIV, B, LDB, &INFO );
#elif defined( USE_VECLIB )
    FORTRAN_WRAPPER( ::zgetrs )
    ( &TRANS, &N, &NRHS, const_cast<Complex *>( A ), &LDA, (int *) IPIV, B, &LDB, &INFO );
#elif defined( USE_OPENBLAS )
    FORTRAN_WRAPPER( ::zgetrs )
    ( &TRANS, &N, &NRHS, const_cast<Complex *>( A ), &LDA, (int *) IPIV, B, &LDB, &INFO, 1 );
#elif defined( USE_MATLAB_LAPACK )
    ptrdiff_t Np = N, NRHSp = NRHS, LDAp = LDA, LDBp = LDB, INFOp;
    ptrdiff_t *IPIVp = new ptrdiff_t[N];
    convert( N, IPIV, IPIVp );
    FORTRAN_WRAPPER( ::zgetrs )
    ( &TRANS, &Np, &NRHSp, (double *) A, &LDAp, IPIVp, (double *) B, &LDBp, &INFOp );
    delete[] IPIVp;
    INFO = static_cast<int>( INFOp );
#else
    FORTRAN_WRAPPER( ::zgetrs )
    ( &TRANS, &N, &NRHS, const_cast<Complex *>( A ), &LDA, (int *) IPIV, B, &LDB, &INFO );
#endif
}
#undef sgttrs
template<>
void Lapack<std::complex<double>>::gttrs( char TRANS, int N, int NRHS,
    const std::complex<double> *DL_, const std::complex<double> *D_,
    const std::complex<double> *DU_, const std::complex<double> *DU2_, const int *IPIV,
    std::complex<double> *B_, int LDB, int &INFO )
{
    auto DL  = reinterpret_cast<const Complex *>( DL_ );
    auto D   = reinterpret_cast<const Complex *>( D_ );
    auto DU  = reinterpret_cast<const Complex *>( DU_ );
    auto DU2 = reinterpret_cast<const Complex *>( DU2_ );
    auto B   = reinterpret_cast<Complex *>( B_ );
#ifdef USE_ATLAS
    throw std::logic_error( "ATLAS does not support sgttrs" );
#elif defined( USE_ACML )
    ::zgttrs( TRANS, N, NRHS, DL, D, DU, DU2, (int *) IPIV, B, LDB, &INFO );
#elif defined( USE_VECLIB )
    FORTRAN_WRAPPER( ::zgttrs )
    ( &TRANS, &N, &NRHS, DL, D, DU, DU2, (int *) IPIV, B, &LDB, &INFO );
#elif defined( USE_OPENBLAS )
    FORTRAN_WRAPPER( ::zgttrs )
    ( &TRANS, &N, &NRHS, DL, D, DU, DU2, (int *) IPIV, B, &LDB, &INFO, 1 );
#elif defined( USE_MATLAB_LAPACK )
    ptrdiff_t Np = N, NRHSp = NRHS, LDBp = LDB, INFOp;
    ptrdiff_t *IPIVp = new ptrdiff_t[N];
    convert( N, IPIV, IPIVp );
    FORTRAN_WRAPPER( ::zgttrs )
    ( &TRANS, &Np, &NRHSp, (double *) DL, (double *) D, (double *) DU, (double *) DU2, IPIVp,
        (double *) B, &LDBp, &INFOp );
    delete[] IPIVp;
    INFO = static_cast<int>( INFOp );
#else
    FORTRAN_WRAPPER( ::zgttrs )
    ( &TRANS, &N, &NRHS, const_cast<Complex *>( DL ), const_cast<Complex *>( D ),
        const_cast<Complex *>( DU ), const_cast<Complex *>( DU2 ), (int *) IPIV, B, &LDB, &INFO );
#endif
}
#undef sgbtrs
template<>
void Lapack<std::complex<double>>::gbtrs( char TRANS, int N, int KL, int KU, int NRHS,
    const std::complex<double> *AB_, int LDAB, const int *IPIV, std::complex<double> *B_, int LDB,
    int &INFO )
{
    auto AB = reinterpret_cast<const Complex *>( AB_ );
    auto B  = reinterpret_cast<Complex *>( B_ );
#ifdef USE_ATLAS
    throw std::logic_error( "ATLAS does not support sgbtrs" );
#elif defined( USE_ACML )
    ::zgbtrs(
        TRANS, N, KL, KU, NRHS, const_cast<Complex *>( AB ), LDAB, (int *) IPIV, B, LDB, &INFO );
#elif defined( USE_VECLIB )
    FORTRAN_WRAPPER( ::zgbtrs )
    ( &TRANS, &N, &KL, &KU, &NRHS, const_cast<Complex *>( AB ), &LDAB, (int *) IPIV, B, &LDB,
        &INFO );
#elif defined( USE_OPENBLAS )
    FORTRAN_WRAPPER( ::zgbtrs )
    ( &TRANS, &N, &KL, &KU, &NRHS, const_cast<Complex *>( AB ), &LDAB, (int *) IPIV, B, &LDB, &INFO,
        1 );
#elif defined( USE_MATLAB_LAPACK )
    ptrdiff_t Np = N, KLp = KL, KUp = KU, NRHSp = NRHS, LDABp = LDAB, LDBp = LDB, INFOp;
    ptrdiff_t *IPIVp = new ptrdiff_t[N];
    convert( N, IPIV, IPIVp );
    FORTRAN_WRAPPER( ::zgbtrs )
    ( &TRANS, &Np, &KLp, &KUp, &NRHSp, (double *) AB, &LDABp, IPIVp, (double *) B, &LDBp, &INFOp );
    INFO = static_cast<int>( INFOp );
    delete[] IPIVp;
#else
    FORTRAN_WRAPPER( ::zgbtrs )
    ( &TRANS, &N, &KL, &KU, &NRHS, const_cast<Complex *>( AB ), &LDAB, (int *) IPIV, B, &LDB,
        &INFO );
#endif
}
#undef sgetri
template<>
void Lapack<std::complex<double>>::getri( int N, std::complex<double> *A_, int LDA, const int *IPIV,
    std::complex<double> *WORK_, int LWORK, int &INFO )
{
    auto A    = reinterpret_cast<Complex *>( A_ );
    auto WORK = reinterpret_cast<Complex *>( WORK_ );
#ifdef USE_ATLAS
    INFO = clapack_sgetri( CblasColMajor, N, A, LDA, IPIV );
#elif defined( USE_ACML )
    ::zgetri_( &N, A, &LDA, (int *) IPIV, WORK, &LWORK, &INFO );
#elif defined( USE_VECLIB )
    FORTRAN_WRAPPER( ::zgetri )( &N, A, &LDA, (int *) IPIV, WORK, &LWORK, &INFO );
#elif defined( USE_OPENBLAS )
    FORTRAN_WRAPPER( ::zgetri )( &N, A, &LDA, (int *) IPIV, WORK, &LWORK, &INFO );
#elif defined( USE_MATLAB_LAPACK )
    ptrdiff_t Np = N, LDAp = LDA, LWORKp = LWORK, INFOp;
    ptrdiff_t *IPIVp = new ptrdiff_t[N];
    convert( N, IPIV, IPIVp );
    FORTRAN_WRAPPER( ::zgetri )
    ( &Np, (double *) A, &LDAp, IPIVp, (double *) WORK, &LWORKp, &INFOp );
    INFO = static_cast<int>( INFOp );
    delete[] IPIVp;
#else
    FORTRAN_WRAPPER( ::zgetri )( &N, A, &LDA, (int *) IPIV, WORK, &LWORK, &INFO );
#endif
}
#undef strsm
template<>
void Lapack<std::complex<double>>::trsm( char SIDE, char UPLO, char TRANS, char DIAG, int M, int N,
    std::complex<double> ALPHA_, const std::complex<double> *A_, int LDA, std::complex<double> *B_,
    int LDB )
{
    auto ALPHA = convert( ALPHA_ );
    auto A     = reinterpret_cast<const Complex *>( A_ );
    auto B     = reinterpret_cast<Complex *>( B_ );
#ifdef USE_ATLAS
    throw std::logic_error( "ztrsm not implemented for ATLAS" );
#elif defined( USE_ACML )
    char SIDE2[2] = { SIDE, 0 }, UPLO2[2] = { UPLO, 0 }, TRANS2[2] = { TRANS, 0 },
         DIAG2[2] = { DIAG, 0 };
    ::ztrsm_( SIDE2, UPLO2, TRANS2, DIAG2, &M, &N, &ALPHA, const_cast<Complex *>( A ), &LDA, B,
        &LDB, 1, 1, 1, 1 );
#elif defined( USE_VECLIB )
    cblas_ztrsm( CblasColMajor, SIDE2( SIDE ), UPLO2( UPLO ), TRANS2( TRANS ), DIAG2( DIAG ), M, N,
        &ALPHA, const_cast<Complex *>( A ), LDA, B, LDB );
#elif defined( USE_OPENBLAS )
    cblas_ztrsm( CblasColMajor, SIDE2( SIDE ), UPLO2( UPLO ), TRANS2( TRANS ), DIAG2( DIAG ), M, N,
        &ALPHA, const_cast<Complex *>( A ), LDA, B, LDB );
#elif defined( USE_MATLAB_LAPACK )
    ptrdiff_t Mp = M, Np = N, LDAp = LDA, LDBp = LDB;
    FORTRAN_WRAPPER( ::ztrsm )
    ( &SIDE, &UPLO, &TRANS, &DIAG, &Mp, &Np, (double *) &ALPHA, (double *) A, &LDAp, (double *) B,
        &LDBp );
#else
    FORTRAN_WRAPPER( ::ztrsm )
    ( &SIDE, &UPLO, &TRANS, &DIAG, &M, &N, &ALPHA, const_cast<Complex *>( A ), &LDA, B, &LDB );
#endif
}
