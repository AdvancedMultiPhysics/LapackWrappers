#include "LapackWrappers.h"
#include "blas_lapack.h"


// Define macro to handle name mangling
#ifndef FORTRAN_WRAPPER
#if defined( USE_ACML )
#define FORTRAN_WRAPPER( x ) x##_
#elif defined( _WIN32 ) || defined( __hpux ) || defined( USE_MKL )
#define FORTRAN_WRAPPER( x ) x
#elif defined( USE_VECLIB )
#define FORTRAN_WRAPPER( x ) x##_
#elif defined( USE_OPENBLAS )
#define FORTRAN_WRAPPER( x ) x##_
#else
#define FORTRAN_WRAPPER( x ) x##_
#endif
#endif


/******************************************************************
 * Get the machine parameters by lamch                             *
 ******************************************************************/
#undef dlamch
template<>
float Lapack<float>::lamch( char cmach )
{
#ifdef USE_ATLAS
    return clapack_dlamch( cmach );
#elif defined( USE_ACML )
    return ::dlamch( cmach );
#elif defined( USE_VECLIB )
    return FORTRAN_WRAPPER( ::dlamch )( &cmach );
#elif defined( USE_OPENBLAS )
    return FORTRAN_WRAPPER( ::dlamch )( &cmach, 1 );
#elif defined( USE_LAPACKE )
    return ::LAPACKE_dlamch( cmach );
#else
    return FORTRAN_WRAPPER( ::dlamch )( &cmach );
#endif
}
template<>
double Lapack<double>::lamch( char cmach )
{
#ifdef USE_ATLAS
    return clapack_dlamch( cmach );
#elif defined( USE_ACML )
    return ::dlamch( cmach );
#elif defined( USE_VECLIB )
    return FORTRAN_WRAPPER( ::dlamch )( &cmach );
#elif defined( USE_OPENBLAS )
    return FORTRAN_WRAPPER( ::dlamch )( &cmach, 1 );
#elif defined( USE_LAPACKE )
    return ::LAPACKE_dlamch( cmach );
#else
    return FORTRAN_WRAPPER( ::dlamch )( &cmach );
#endif
}


/******************************************************************
 * Get the machine parameters by lamch                             *
 ******************************************************************/
template<typename TYPE>
LapackMachineParams Lapack<TYPE>::machineParams()
{
    LapackMachineParams data;
    data.eps   = Lapack<TYPE>::lamch( 'E' );
    data.sfmin = Lapack<TYPE>::lamch( 'S' );
    data.base  = Lapack<TYPE>::lamch( 'B' );
    data.prec  = Lapack<TYPE>::lamch( 'P' );
    data.t     = Lapack<TYPE>::lamch( 'N' );
    data.rnd   = Lapack<TYPE>::lamch( 'R' );
    data.emin  = Lapack<TYPE>::lamch( 'M' );
    data.rmin  = Lapack<TYPE>::lamch( 'U' );
    data.emax  = Lapack<TYPE>::lamch( 'L' );
    data.rmax  = Lapack<TYPE>::lamch( 'O' );
    return data;
}
template LapackMachineParams Lapack<double>::machineParams();
template LapackMachineParams Lapack<float>::machineParams();
std::string LapackMachineParams::print() const
{
    char msg[1000];
    char *tmp = msg;
    tmp += sprintf( tmp, "  eps   = %-13.6e    relative machine precision\n", eps );
    tmp += sprintf( tmp, "  sfmin = %-13.6e    safe minimum\n", sfmin );
    tmp += sprintf( tmp, "  base  = %-11i      base of the machine\n", base );
    tmp += sprintf( tmp, "  prec  = %-13.6e    eps*base\n", prec );
    tmp += sprintf( tmp, "  t     = %-11i      number of digits in the mantissa\n", t );
    tmp += sprintf( tmp, "  rnd   = %-11i      1 when rounding occurs in addition, 0 otherwise\n",
        rnd ? 1 : 0 );
    tmp += sprintf( tmp, "  emin  = %-11i      minimum exponent before underflow\n", emin );
    tmp += sprintf( tmp, "  rmin  = %-13.6e    underflow threshold - base**(emin-1)\n", rmin );
    tmp += sprintf( tmp, "  emax  = %-11i      largest exponent before overflow\n", emax );
    tmp += sprintf( tmp, "  rmax  = %-13.6e    overflow threshold - (base**emax)*(1-eps)\n", rmax );
    return std::string( msg );
}
