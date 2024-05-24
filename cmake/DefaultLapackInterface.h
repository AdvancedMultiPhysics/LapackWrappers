#ifndef FORTRAN_LAPACK_CALLS
#define FORTRAN_LAPACK_CALLS

#include <complex>

extern "C" {

#if defined( _WIN32 ) || defined( __hpux )
#define FORTRAN_WRAPPER( x ) x
#else
#define FORTRAN_WRAPPER( x ) x##_
#endif

// double precision functions
// Misc
extern void FORTRAN_WRAPPER( dfill )( int*, double*, double* );
extern int FORTRAN_WRAPPER( idamax )( int*, double*, int* );
extern double FORTRAN_WRAPPER( dlamch )( char* );
// Level 1 BLAS Routines
extern void FORTRAN_WRAPPER( dswap )( int*, double*, int*, double*, int* );
extern void FORTRAN_WRAPPER( dscal )( int*, double*, double*, int* );
extern void FORTRAN_WRAPPER( dcopy )( int*, double*, int*, double*, int* );
extern void FORTRAN_WRAPPER( daxpy )( int*, double*, double*, int*, double*, int* );
extern double FORTRAN_WRAPPER( ddot )( int*, double*, int*, double*, int* );
extern double FORTRAN_WRAPPER( dasum )( int*, double*, int* );
extern double FORTRAN_WRAPPER( damax )( int*, double*, int* );
extern double FORTRAN_WRAPPER( dnrm2 )( int*, double*, int* );
// Level 2 BLAS Routines
extern void FORTRAN_WRAPPER( dgemv )(
    char*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int* );
extern void FORTRAN_WRAPPER( dger )(
    int*, int*, double*, double*, int*, double*, int*, double*, int* );
// Level 3 BLAS Routines
extern void FORTRAN_WRAPPER( dgemm )(
    char*, char*, int*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int* );
// LAPACK Routines
extern void FORTRAN_WRAPPER( dgesv )( int*, int*, double*, int*, int*, double*, int*, int* );
extern void FORTRAN_WRAPPER( dgtsv )( int*, int*, double*, double*, double*, double*, int*, int* );
extern void FORTRAN_WRAPPER( dgbsv )(
    int*, int*, int*, int*, double*, int*, int*, double*, int*, int* );
extern void FORTRAN_WRAPPER( dgbmv )(
    char*, int*, int*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int* );
extern void FORTRAN_WRAPPER( dgetrf )( int*, int*, double*, int*, int*, int* );
extern void FORTRAN_WRAPPER( dgttrf )( int*, double*, double*, double*, double*, int*, int* );
extern void FORTRAN_WRAPPER( dgbtrf )( int*, int*, int*, int*, double*, int*, int*, int* );
extern void FORTRAN_WRAPPER( dgetrs )(
    char*, int*, int*, double*, int*, int*, double*, int*, int* );
extern void FORTRAN_WRAPPER( dgttrs )(
    char*, int*, int*, double*, double*, double*, double*, int*, double*, int*, int* );
extern void FORTRAN_WRAPPER( dgbtrs )(
    char*, int*, int*, int*, int*, double*, int*, int*, double*, int*, int* );
extern void FORTRAN_WRAPPER( dgetri )( int*, double*, int*, int*, double*, int*, int* );
extern void FORTRAN_WRAPPER( dtrsm )(
    char*, char*, char*, char*, int*, int*, double*, double*, int*, double*, int* );

// single precision functions
// Misc
extern void FORTRAN_WRAPPER( sfill )( int*, float*, float* );
extern int FORTRAN_WRAPPER( isamax )( int*, float*, int* );
extern float FORTRAN_WRAPPER( slamch )( char* );
// Level 1 BLAS Routines
extern void FORTRAN_WRAPPER( sswap )( int*, float*, int*, float*, int* );
extern void FORTRAN_WRAPPER( sscal )( int*, float*, float*, int* );
extern void FORTRAN_WRAPPER( scopy )( int*, float*, int*, float*, int* );
extern void FORTRAN_WRAPPER( saxpy )( int*, float*, float*, int*, float*, int* );
extern float FORTRAN_WRAPPER( sdot )( int*, float*, int*, float*, int* );
extern float FORTRAN_WRAPPER( sasum )( int*, float*, int* );
extern float FORTRAN_WRAPPER( samax )( int*, float*, int* );
extern float FORTRAN_WRAPPER( snrm2 )( int*, float*, int* );
// Level 2 BLAS Routines
extern void FORTRAN_WRAPPER( sgemv )(
    char*, int*, int*, float*, float*, int*, float*, int*, float*, float*, int* );
extern void FORTRAN_WRAPPER( sger )( int*, int*, float*, float*, int*, float*, int*, float*, int* );
// Level 3 BLAS Routines
extern void FORTRAN_WRAPPER( sgemm )(
    char*, char*, int*, int*, int*, float*, float*, int*, float*, int*, float*, float*, int* );

// LAPACK Routines
extern void FORTRAN_WRAPPER( sgesv )( int*, int*, float*, int*, int*, float*, int*, int* );
extern void FORTRAN_WRAPPER( sgtsv )( int*, int*, float*, float*, float*, float*, int*, int* );
extern void FORTRAN_WRAPPER( sgbsv )(
    int*, int*, int*, int*, float*, int*, int*, float*, int*, int* );
extern void FORTRAN_WRAPPER( sgbmv )(
    char*, int*, int*, int*, int*, float*, float*, int*, float*, int*, float*, float*, int* );
extern void FORTRAN_WRAPPER( sgetrf )( int*, int*, float*, int*, int*, int* );
extern void FORTRAN_WRAPPER( sgttrf )( int*, float*, float*, float*, float*, int*, int* );
extern void FORTRAN_WRAPPER( sgbtrf )( int*, int*, int*, int*, float*, int*, int*, int* );
extern void FORTRAN_WRAPPER( sgetrs )( char*, int*, int*, float*, int*, int*, float*, int*, int* );
extern void FORTRAN_WRAPPER( sgttrs )(
    char*, int*, int*, float*, float*, float*, float*, int*, float*, int*, int* );
extern void FORTRAN_WRAPPER( sgbtrs )(
    char*, int*, int*, int*, int*, float*, int*, int*, float*, int*, int* );
extern void FORTRAN_WRAPPER( sgetri )( int*, float*, int*, int*, float*, int*, int* );
extern void FORTRAN_WRAPPER( strsm )(
    char*, char*, char*, char*, int*, int*, float*, float*, int*, float*, int* );

// complex double precision functions
typedef double _Complex Complex;
// Misc
extern void FORTRAN_WRAPPER( zfill )( int*, Complex*, Complex* );
extern int FORTRAN_WRAPPER( izamax )( int*, Complex*, int* );
extern Complex FORTRAN_WRAPPER( zlamch )( char* );
// Level 1 BLAS Routines
extern void FORTRAN_WRAPPER( zswap )( int*, Complex*, int*, Complex*, int* );
extern void FORTRAN_WRAPPER( zscal )( int*, Complex*, Complex*, int* );
extern void FORTRAN_WRAPPER( zcopy )( int*, Complex*, int*, Complex*, int* );
extern void FORTRAN_WRAPPER( zaxpy )( int*, Complex*, Complex*, int*, Complex*, int* );
extern double FORTRAN_WRAPPER( zdotc )( int*, Complex*, int*, Complex*, int* );
extern double FORTRAN_WRAPPER( dzasum )( int*, Complex*, int* );
extern double FORTRAN_WRAPPER( zamax )( int*, Complex*, int* );
extern double FORTRAN_WRAPPER( dznrm2 )( int*, Complex*, int* );
// Level 2 BLAS Routines
extern void FORTRAN_WRAPPER( zgemv )(
    char*, int*, int*, Complex*, Complex*, int*, Complex*, int*, Complex* eta, Complex*, int* );
extern void FORTRAN_WRAPPER( zgerc )(
    int*, int*, Complex*, Complex*, int*, Complex*, int*, Complex*, int* );
// Level 3 BLAS Routines
extern void FORTRAN_WRAPPER( zgemm )( char*, char*, int*, int*, int*, Complex*, Complex*, int*,
    Complex*, int*, Complex*, Complex*, int* );
// LAPACK Routines
extern void FORTRAN_WRAPPER( zgesv )( int*, int*, Complex*, int*, int*, Complex*, int*, int* );
extern void FORTRAN_WRAPPER( zgtsv )(
    int*, int*, Complex*, Complex*, Complex*, Complex*, int*, int* );
extern void FORTRAN_WRAPPER( zgbsv )(
    int*, int*, int*, int*, Complex*, int*, int*, Complex*, int*, int* );
extern void FORTRAN_WRAPPER( zgbmv )( char*, int*, int*, int*, int*, Complex*, Complex*, int*,
    Complex*, int*, Complex*, Complex*, int* );
extern void FORTRAN_WRAPPER( zgetrf )( int*, int*, Complex*, int*, int*, int* );
extern void FORTRAN_WRAPPER( zgttrf )( int*, Complex*, Complex*, Complex*, Complex*, int*, int* );
extern void FORTRAN_WRAPPER( zgbtrf )( int*, int*, int*, int*, Complex*, int*, int*, int* );
extern void FORTRAN_WRAPPER( zgetrs )(
    char*, int*, int*, Complex*, int*, int*, Complex*, int*, int* );
extern void FORTRAN_WRAPPER( zgttrs )(
    char*, int*, int*, Complex*, Complex*, Complex*, Complex*, int*, Complex*, int*, int* );
extern void FORTRAN_WRAPPER( zgbtrs )(
    char*, int*, int*, int*, int*, Complex*, int*, int*, Complex*, int*, int* );
extern void FORTRAN_WRAPPER( zgetri )( int*, Complex*, int*, int*, Complex*, int*, int* );
extern void FORTRAN_WRAPPER( ztrsm )(
    char*, char*, char*, char*, int*, int*, Complex*, Complex*, int*, Complex*, int* );
}

#endif
