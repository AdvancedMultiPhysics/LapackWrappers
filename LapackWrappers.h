#ifndef USE_LAPACK_WRAPPER
#define USE_LAPACK_WRAPPER

#include <cstdint>
#include <string>
#include <vector>


struct LapackMachineParams {
    bool rnd;                  // true if rounding occurs in addition, false otherwise
    uint8_t base;              // base of the machine
    uint8_t t;                 // number of digits in the mantissa
    short emin;                // minimum exponent before underflow
    short emax;                // largest exponent before overflow
    double eps;                // relative machine precision
    double sfmin;              // safe minimum
    double prec;               // eps*base
    double rmin;               // underflow threshold - base**(emin-1)
    double rmax;               // overflow threshold - (base**emax)*(1-eps)
    std::string print() const; // Print the machine parameters
};


/*! \class Lapack
    \brief A wrapper class for BLAS/LAPACK routines
    \details  This class is a wrapper class to provide a common interface to
      BLAS and LAPACK routines.
*/
template<typename TYPE>
class Lapack
{
public:
    /*!
     * \brief   copy
     * \details  copy copies a vector x, to a vector y
     * @param[in]  N        The number of values to copy
     * @param[in]  x        The source vector
     * @param[in]  INCX     The spacing between points in x
     * @param[in]  y        The destination vector
     * @param[out] INCY     The spacing between points in y
     */
    static void copy( int N, const TYPE *x, int INCX, TYPE *y, int INCY );

    /*!
     * \brief   swap
     * \details  swap swaps two vectors
     * @param[in]     N     The number of values to copy
     * @param[in,out] x     The first vector
     * @param[in]  INCX     The spacing between points in x
     * @param[in,out] y     The second vector
     * @param[in]  INCY     The spacing between points in x
     */
    static void swap( int N, TYPE *x, int INCX, TYPE *y, int INCY );

    /*!
     * \brief   scal
     * \details  scal scales a vector by a constant.  x = a*x
     * @param[in]  N        The number of values to copy
     * @param[in]  a        The scale factor
     * @param[in,out] x     The vector
     * @param[in]  INCX     The spacing between points in x
     */
    static void scal( int N, TYPE a, TYPE *x, int INCX );

    /*!
     * \brief   nrm2
     * \details  nrm2 returns the euclidean norm of a vector via the function
     * @param[in]  N        The number of values to copy
     * @param[in]  x        The input vector
     * @param[in]  INCX     The spacing between points in x
     */
    static double nrm2( int N, const TYPE *x, int INCX );

    /*!
     * \brief   iamax
     * \details  iamax finds the index of element having maximum absolute value.
     *    Note: the returned index is 0 (C++) based.
     * @param[in]  N        The number of values to copy
     * @param[in]  x        The input vector
     * @param[in]  INCX     The spacing between points in x
     */
    static int iamax( int N, const TYPE *x, int INCX );

    /*!
     * \brief   axpy
     * \details  axpy scales a vector by a constant plus a vector.  y = a*x + y
     * @param[in] N         The number of values to copy
     * @param[in] a         The scale factor
     * @param[in] x         The source vector
     * @param[in] INCX      The spacing between points in x
     * @param[in,out] y     The destination vector
     * @param[in] INCY      The spacing between points in x
     */
    static void axpy( int N, TYPE a, const TYPE *x, int INCX, TYPE *y, int INCY );

    /*!
     * \brief   gemv
     * \details  gemv performs one of the matrix-vector operations
     *      y := alpha*A*x + beta*y, or y := alpha*A'*x + beta*y,
     *   where alpha and beta are scalars, x and y are vectors and A
     *   is an m by n matrix.
     * @param[in] TRANS     On entry, TRANS specifies the operation to be performed as follows:
     *                      TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
     *                      TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
     *                      TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
     * @param[in] M         The number of rows of the matrix A.  M >= 0
     * @param[in] N         The number of columns of the matrix A.  N >= 0
     * @param[in] alpha     The scale factor alpha
     * @param[in] A         The leading m by n part of the array A must contain the matrix of
     * coefficients.
     * @param[in] LDA       LDA specifies the first dimension of A as declared in the
     *                      calling (sub) program. LDA must be at least max( 1, m ).
     * @param[in] x         The source vector
     * @param[in] INCX      The spacing between points in x
     * @param[in] beta      The scale factor beta
     * @param[in,out] y     The destination vector
     * @param[in] INCY      The spacing between points in y
     */
    static void gemv( char TRANS, int M, int N, TYPE alpha, const TYPE *A, int LDA, const TYPE *x,
        int INCX, TYPE beta, TYPE *y, int INCY );

    /*!
     * \brief   gemm
     * \details  gemm performs one of the matrix-matrix operations
     *       C := alpha*op( A )*op( B ) + beta*C,
     *    where  op( X ) is one of
     *       op( X ) = X   or   op( X ) = X',
     *    alpha and beta are scalars, and A, B and C are matrices,
     *    with op( A ) an m by k matrix,  op( B )  a  k by n matrix
     *    and  C an m by n matrix.
     * @param[in] TRANSA    TRANSA specifies the form of op( A ) to be
     *                      used in the matrix multiplication as follows:
     *                      TRANSA = 'N' or 'n',  op( A ) = A.
     *                      TRANSA = 'T' or 't',  op( A ) = A'.
     *                      TRANSA = 'C' or 'c',  op( A ) = A'.
     * @param[in] TRANSB    TRANSA specifies the form of op( B ) to be
     *                      used in the matrix multiplication as follows:
     *                      TRANSA = 'N' or 'n',  op( B ) = B.
     *                      TRANSA = 'T' or 't',  op( B ) = B'.
     *                      TRANSA = 'C' or 'c',  op( B ) = B'.
     * @param[in] M         M specifies  the number of rows of the matrix op( A )
     *                      and of the matrix C.  M must be at least zero.
     * @param[in] N         N specifies the number of columns of the matrix op( B )
     *                      and of the matrix C.  N must be at least zero.
     * @param[in] K         K specifies  the number of columns of the matrix op( A )
     *                      and the number of rows of the matrix op( B ).
     *                      K must be at least  zero.
     * @param[in] alpha     The scalar alpha.
     * @param[in] A         array of DIMENSION ( LDA, ka ), where
     *                      k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
     *                      Before entry with  TRANSA = 'N' or 'n',  the leading
     *                      m by k part of the array  A  must contain the matrix
     *                      A,  otherwise the leading  k by m  part of the array
     *                      A  must contain  the matrix A.
     * @param[in] LDA       LDA specifies the first dimension of A as
     *                      declared in the calling (sub) program. When  TRANSA =
     *                      'N' or 'n' then LDA must be at least  max( 1, m ),
     *                      otherwise  LDA must be at least  max( 1, k ).
     * @param[in] B         array of DIMENSION ( LDB, kb ), where
     *                      n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
     *                      Before entry with  TRANSB = 'N' or 'n',  the leading
     *                      k by n part of the array  B  must contain the matrix
     *                      B,  otherwise the leading  n by k  part of the array
     *                      B  must contain  the matrix B.  Unchanged on exit.
     * @param[in] LDB       LDB specifies the first dimension of B as
     *                      declared in the calling (sub) program. When  TRANSB =
     *                      'N' or 'n' then LDB must be at least  max( 1, k ),
     *                      otherwise  LDB must be at least  max( 1, n ).
     * @param[in] beta      The scalar  beta.  When beta is supplied as zero then
     *                      C need not be set on input.  Unchanged on exit.
     * @param[in,out] C     Array of DIMENSION ( LDC, n ).
     *                      Before entry, the leading  m by n  part of the array
     *                      C must contain the matrix  C,  except when  beta  is
     *                      zero, in which case C need not be set on entry.  On
     *                      exit, the array  C  is overwritten by the  m by n
     *                      matrix ( alpha*op( A )*op( B ) + beta*C ).
     * @param[in] LDC    - INTEGER.
     *                      On entry, LDC specifies the first dimension of C as
     *                      declared in  the  calling  (sub)  program.   LDC
     *                      must  be  at  least max( 1, m ).  Unchanged on exit.
     */
    static void gemm( char TRANSA, char TRANSB, int M, int N, int K, TYPE alpha, const TYPE *A,
        int LDA, const TYPE *B, int LDB, TYPE beta, TYPE *C, int LDC );

    /*!
     * \brief   asum
     * \details  asum sums a vector
     * @param[in]  N        The number of values to copy
     * @param[in]  x        The source vector
     * @param[in]  INCX     The spacing between points in x
     */
    static double asum( int N, const TYPE *x, int INCX );

    /*!
     * \brief   dot
     * \details  dot computes the dot product between two vectors
     * @param[in]  N        The number of values to copy
     * @param[in]  x        The source vector x
     * @param[in]  INCX     The spacing between points in x
     * @param[in]  y        The source vector y
     * @param[in]  INCY     The spacing between points in y
     */
    static TYPE dot( int N, const TYPE *x, int INCX, const TYPE *y, int INCY );

    /*!
     * \brief   ger
     * \details  ger performs the rank 1 operation
     *     A := alpha*x*y' + A,
     *  where alpha is a scalar, x is an m element vector, y is an n element
     *  vector and A is an m by n matrix.
     * @param[in]  M        The number of rows of the matrix A
     * @param[in]  N        The number of columns of the matrix A
     * @param[in]  alpha    The scalar alpha
     * @param[in]  x        The source vector x
     * @param[in]  INCX     The spacing between points in x
     * @param[in]  y        The source vector y
     * @param[in]  INCY     The spacing between points in y
     * @param[in,out]  A    On entry, the N-by-N coefficient matrix A.
     *                      On exit, the updated matrix.
     * @param[in]  LDA      The leading dimension of the array A.  LDA >= max(1,M).
     */
    static void ger( int N, int M, TYPE alpha, const TYPE *x, int INCX, const TYPE *y, int INCY,
        TYPE *A, int LDA );

    /*!
     * \brief   gesv
     * \details  gesv computes the solution to a real system of linear equations
     *       A * X = B,
     *    where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
     *    The LU decomposition with partial pivoting and row interchanges is
     *    used to factor A as
     *       A = P * L * U,
     *    where P is a permutation matrix, L is unit lower triangular, and U is
     *    upper triangular.  The factored form of A is then used to solve the
     *    system of equations A * X = B.
     * @param[in]  N        The number of linear equations, i.e., the order of the
     *                      matrix A.  N >= 0.
     * @param[in]  NRHS     The number of right hand sides, i.e., the number of columns
     *                      of the matrix B.  NRHS >= 0.
     * @param[in]  A        On entry, the N-by-N coefficient matrix A.
     *                      On exit, the factors L and U from the factorization
     *                      A = P*L*U; the unit diagonal elements of L are not stored.
     * @param[in]  LDA      The leading dimension of the array A.  LDA >= max(1,N).
     * @param[out] IPIV     The pivot indices that define the permutation matrix P;
     *                      row i of the matrix was interchanged with row IPIV(i).
     * @param[in,out] B     On entry, the N-by-NRHS matrix of right hand side matrix B.
     *                      On exit, if INFO = 0, the N-by-NRHS solution matrix X.
     * @param[in]  LDB      The leading dimension of the array B.  LDB >= max(1,N).
     * @param[out] INFO     Exit code
     *                      = 0:  successful exit
     *                      < 0:  if INFO = -i, the i-th argument had an illegal value
     *                      > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
     *                            has been completed, but the factor U is exactly
     *                            singular, so the solution could not be computed.
     */
    static void gesv( int N, int NRHS, TYPE *A, int LDA, int *IPIV, TYPE *B, int LDB, int &INFO );

    /*!
     * \brief   gtsv
     * \details  gtsv computes the solution to a real system of linear equations
     *       A * X = B,
     *    where A is an n by n tridiagonal matrix, by Gaussian elimination with
     *    partial pivoting.  Note that the equation  A**T*X = B  may be solved by interchanging the
     *    order of the arguments DU and DL.
     * @param[in]  N        The number of linear equations, i.e., the order of the
     *                      matrix A.  N >= 0.
     * @param[in]  NRHS     The number of right hand sides, i.e., the number of columns
     *                      of the matrix B.  NRHS >= 0.
     * @param[in,out] DL    On entry, DL must contain the (n-1) subdiagonal elements of A.
     *                      On exit, DL is overwritten by the (n-2) elements of the second
     *                      superdiagonal of the upper triangular matrix U from the LU
     *                      factorization of A, in DL(1), ..., DL(n-2).
     * @param[in,out] D     On entry, D must contain the diagonal elements of A.
     *                      On exit, D is overwritten by the n diagonal elements of U.
     * @param[in,out] DU    On entry, DU must contain the (n-1) superdiagonal
     *                      elements of A.  On exit, DU is overwritten by the
     *                      (n-1) elements of the first superdiagonal of U.
     * @param[in,out] B     On entry, the N-by-NRHS matrix of right hand side matrix B.
     *                      On exit, if INFO = 0, the N-by-NRHS solution matrix X.
     * @param[in]  LDB      The leading dimension of the array B.  LDB >= max(1,N).
     * @param[out] INFO     Exit code
     *                      = 0:  successful exit
     *                      < 0:  if INFO = -i, the i-th argument had an illegal value
     *                      > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
     *                            has been completed, but the factor U is exactly
     *                            singular, so the solution could not be computed.
     */
    static void gtsv( int N, int NRHS, TYPE *DL, TYPE *D, TYPE *DU, TYPE *B, int LDB, int &INFO );

    /*!
     * \brief   gbsv
     * \details  gbsv computes the solution to a real system of linear equations
     *       A * X = B,
     *    where A is a band matrix of order N with KL subdiagonals
     *    and KU superdiagonals, and X and B are N-by-NRHS matrices.
     *    The LU decomposition with partial pivoting and row interchanges is
     *    used to factor A as A = L * U, where L is a product of permutation
     *    and unit lower triangular matrices with KL subdiagonals, and U is
     *    upper triangular with KL+KU superdiagonals.  The factored form of A
     *    is then used to solve the system of equations A * X = B.
     *
     *    Further Details
     *    ===============
     *
     *    The band storage scheme is illustrated by the following example, when
     *    M = N = 6, KL = 2, KU = 1:
     *
     *    On entry:                       On exit:
     *
     *        *    *    *    +    +    +       *    *    *   u14  u25  u36
     *        *    *    +    +    +    +       *    *   u13  u24  u35  u46
     *        *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
     *       a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
     *       a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   *
     *       a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    *
     *
     *    Array elements marked * are not used by the routine; elements marked
     *    + need not be set on entry, but are required by the routine to store
     *    elements of U because of fill-in resulting from the row interchanges.
     *
     * @param[in]  N        The number of linear equations, i.e., the order of the
     *                      matrix A.  N >= 0.
     * @param[in]  KL       The number of subdiagonals within the band of A.  KL >= 0.
     * @param[in]  KU       The number of superdiagonals within the band of A.  KU >= 0.
     * @param[in]  NRHS     The number of right hand sides, i.e., the number of columns
     *                      of the matrix B.  NRHS >= 0.
     * @param[in,out] AB    On entry, the matrix A in band storage, in rows KL+1 to
     *                      2*KL+KU+1; rows 1 to KL of the array need not be set.
     *                      The j-th column of A is stored in the j-th column of the
     *                      array AB as follows:
     *                      AB(KL+KU+1+i-j,j) = A(i,j) for max(1,j-KU)<=i<=min(N,j+KL)
     *                      On exit, details of the factorization: U is stored as an
     *                      upper triangular band matrix with KL+KU superdiagonals in
     *                      rows 1 to KL+KU+1, and the multipliers used during the
     *                      factorization are stored in rows KL+KU+2 to 2*KL+KU+1.
     * @param[out] LDAB     The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
     * @param[out] IPIV     The pivot indices that define the permutation matrix P;
     *                      row i of the matrix was interchanged with row IPIV(i).
     * @param[in,out] B     On entry, the N-by-NRHS matrix of right hand side matrix B.
     *                      On exit, if INFO = 0, the N-by-NRHS solution matrix X.
     * @param[in]  LDB      The leading dimension of the array B.  LDB >= max(1,N).
     * @param[out] INFO     Exit code
     *                      = 0:  successful exit
     *                      < 0:  if INFO = -i, the i-th argument had an illegal value
     *                      > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
     *                            has been completed, but the factor U is exactly
     *                            singular, so the solution could not be computed.
     */
    static void gbsv( int N, int KL, int KU, int NRHS, TYPE *AB, int LDAB, int *IPIV, TYPE *B,
        int LDB, int &INFO );


    /*!
     * \brief   getrf
     * \details  getrf computes an LU factorization of a real matrix A
     *    using elimination with partial pivoting and row interchanges.
     *    The factorization has the form:
     *       A = L * U
     *    where L is a product of permutation and unit lower bi-diagonal
     *    matrices and U is upper triangular with nonzeros in only the
     *    main diagonal and first two superdiagonals.
     * @param[in]  M        The number of rows of the matrix A.  M >= 0.
     * @param[in]  N        The number of columns of the matrix A.  N >= 0.
     * @param[in,out] A     On entry, the M-by-N matrix to be factored.
     *                      On exit, the factors L and U from the factorization A =
     *                      P*L*U; the unit diagonal elements of L are not stored.
     * @param[in]  LDA      The leading dimension of the array A.  LDA >= max(1,M).
     * @param[out] IPIV     The pivot indices that define the permutation matrix P;
     *                      row i of the matrix was interchanged with row IPIV(i).
     * @param[out] INFO     Exit code
     *                      = 0:  successful exit
     *                      < 0:  if INFO = -i, the i-th argument had an illegal value
     *                      > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
     *                            has been completed, but the factor U is exactly
     *                            singular, so the solution could not be computed.
     */
    static void getrf( int M, int N, TYPE *A, int LDA, int *IPIV, int &INFO );

    /*!
     * \brief   gttrf
     * \details  gttrf computes an LU factorization of a real tridiagonal
     *    matrix A using elimination with partial pivoting and row interchanges.
     *    The factorization has the form
     *       A = L * U
     *    where L is a product of permutation and unit lower bidiago-
     *    nal matrices and U is upper triangular with nonzeros in only
     *    the main diagonal and first two superdiagonals.
     * @param[in]  N        The number of linear equations, i.e., the order of the
     *                      matrix A.  N >= 0.
     * @param[in,out] DL    On entry, DL must contain the (n-1) subdiagonal elements of A.
     *                      On exit, DL is overwritten by the (n-2) elements of the second
     *                      superdiagonal of the upper triangular matrix U from the LU
     *                      factorization of A, in DL(1), ..., DL(n-2).
     *
     * @param[in,out] D     On entry, D must contain the diagonal elements of A.
     *                      On exit, D is overwritten by the n diagonal elements of U.
     * @param[in,out] DU    On entry, DU must contain the (n-1) superdiagonal
     *                      elements of A.  On exit, DU is overwritten by the
     *                      (n-1) elements of the first superdiagonal of U.
     * @param[out] DU2      On exit, DU2 is overwritten by the (n-2) elements of
     *                      the second superdiagonal of U.
     * @param[out] IPIV     The pivot indices that define the permutation matrix P;
     *                      row i of the matrix was interchanged with row IPIV(i).
     * @param[out] INFO     Exit code
     *                      = 0:  successful exit
     *                      < 0:  if INFO = -i, the i-th argument had an illegal value
     *                      > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
     *                            has been completed, but the factor U is exactly
     *                            singular, so the solution could not be computed.
     */
    static void gttrf( int N, TYPE *DL, TYPE *D, TYPE *DU, TYPE *DU2, int *IPIV, int &INFO );

    /*!
     * \brief   gbtrf
     * \details  gbtrf computes an LU factorization of a real m-by-n band
     *    matrix A using partial pivoting with row interchanges.
     * @param[in]  M        The number of rows of the matrix A.  M >= 0.
     * @param[in]  N        The number of columns of the matrix A.  N >= 0.
     * @param[in] KL        The number of subdiagonals within the band of A.  KL >= 0.
     * @param[in] KU        The number of superdiagonals within the band of A.  KU >= 0.
     * @param[in,out] AB    On entry, the matrix A in band storage, in rows KL+1 to
     *                         2*KL+KU+1; rows 1 to KL of the array need not be set.
     *                         The j-th column of A is stored in the j-th column of the array AB as
     *                         follows: AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
     *                      On exit, details of the factorization: U is stored as an upper
     *                          triangular band matrix with KL+KU super-diagonals in rows
     *                          1 to KL+KU+1, and the multipliers used during the factorization
     *                          are stored in row KL+KU+2 to 2*KL+KU+1.
     * @param[in] LDAB      The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
     * @param[out] IPIV     The pivot indices; for 1 <= i <= min(M,N), row i of
     *                      the matrix was interchanged with row IPIV(i).
     * @param[out] INFO     Exit code
     *                      = 0:  successful exit
     *                      < 0:  if INFO = -i, the i-th argument had an illegal value
     *                      > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
     *                            has been completed, but the factor U is exactly
     *                            singular, and division by zero will occur if
     *                            it is used to solve a system of equations.
     */
    static void gbtrf( int M, int N, int KL, int KU, TYPE *AB, int LDAB, int *IPIV, int &INFO );

    /*!
     * \brief   getrs
     * \details  getrs solves a system of linear equations
     *       A*X = B  or  A'*X = B
     *    with a general N-by-N matrix A using the LU factorization computed by DGETRF.
     * @param[in] TRANS     Specifies the form of the system of equations:
     *                      = 'N':  A * X = B  (No transpose)
     *                      = 'T':  A'* X = B  (Transpose)
     *                      = 'C':  A'* X = B  (Conjugate transpose = Transpose)
     * @param[in] N         The order of the matrix A.  N >= 0.
     * @param[in] NRHS      The number of right hand sides, i.e., the number of
     *                      columns of the matrix B.  NRHS >= 0.
     * @param[in] A         The factors L and U from the factorization A = P*L*U as computed by
     * DGETRF.
     * @param[in] LDA       The leading dimension of the array A.  LDA >= max(1,N).
     * @param[in] IPIV      The pivot indices from DGETRF; for 1<=i<=N, row i of
     *                      the matrix was interchanged with row IPIV(i).  (LDB,NRHS)
     * @param[in,out] B     On entry, the right hand side matrix B.  On exit, the solution matrix X.
     * @param[in] LDB       The leading dimension of the array B.  LDB >= max(1,N).
     * @param[out] INFO     Exit code
     *                      = 0:  successful exit
     *                      < 0:  if INFO = -i, the i-th argument had an illegal value
     */
    static void getrs( char TRANS, int N, int NRHS, const TYPE *A, int LDA, const int *IPIV,
        TYPE *B, int LDB, int &INFO );

    /*!
     * \brief   gttrs
     * \details  gttrs solves one of the systems of equations
     *       A*X = B  or  A'*X = B, with a tridiagonal matrix A using
     *    the LU factorization computed by DGTTRF.
     * @param[in] TRANS     Specifies the form of the system of equations:
     *                      = 'N':  A * X = B  (No transpose)
     *                      = 'T':  A'* X = B  (Transpose)
     *                      = 'C':  A'* X = B  (Conjugate transpose = Transpose)
     * @param[in] N         The order of the matrix A.  N >= 0.
     * @param[in] NRHS      The number of right hand sides, i.e., the number of
     *                      columns of the matrix B.  NRHS >= 0.
     * @param[in] DL        The (n-1) multipliers that define the matrix L from
     *                      the LU factorization of A.
     * @param[in] D         The n diagonal elements of the upper triangular
     *                      matrix U from the LU factorization of A.
     * @param[in] DU        The (n-1) elements of the first superdiagonal of U.
     * @param[in] DU2       The (n-2) elements of the second superdiagonal of U.
     * @param[in] IPIV      The pivot indices; for 1 <= i <= n, row i of the
     *                      matrix was interchanged with row IPIV(i).  IPIV(i)
     *                      will always be either i or i+1; IPIV(i) = i indi-
     *                      cates a row interchange was not required.
     * @param[in,out] B     On entry, the right hand side matrix B.  On exit, the solution matrix X.
     * @param[in] LDB       The leading dimension of the array B.  LDB >= max(1,N).
     * @param[out] INFO     Exit code
     *                      = 0:  successful exit
     *                      < 0:  if INFO = -i, the i-th argument had an illegal value
     */
    static void gttrs( char TRANS, int N, int NRHS, const TYPE *DL, const TYPE *D, const TYPE *DU,
        const TYPE *DU2, const int *IPIV, TYPE *B, int LDB, int &INFO );

    /*!
     * \brief   gbtrs
     * \details  gbtrs solves a system of linear equations
     *       A * X = B  or  A' * X = B with a general band matrix A
     *    using the LU factorization computed by DGBTRF.
     * @param[in] TRANS     Specifies the form of the system of equations:
     *                      = 'N':  A * X = B  (No transpose)
     *                      = 'T':  A'* X = B  (Transpose)
     *                      = 'C':  A'* X = B  (Conjugate transpose = Transpose)
     * @param[in] N         The order of the matrix A.  N >= 0.
     * @param[in] KL        The number of subdiagonals within the band of A.  KL >= 0.
     * @param[in] KU        The number of superdiagonals within the band of A.  KU >= 0.
     * @param[in] NRHS      The number of right hand sides, i.e., the number of
     *                      columns of the matrix B.  NRHS >= 0.
     * @param[in] AB        Details of the LU factorization of the band matrix
     *                      A, as computed by DGBTRF.  U is stored as an upper
     *                      triangular band matrix with KL+KU superdiagonals in
     *                      rows 1 to KL+KU+1, and the multipliers used during
     *                      the factorization are stored in rows KL+KU+2 to
     *                      2*KL+KU+1.
     * @param[in] LDAB       The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
     * @param[in] IPIV      The pivot indices; for 1 <= i <= n, row i of the
     *                      matrix was interchanged with row IPIV(i).  IPIV(i)
     *                      will always be either i or i+1; IPIV(i) = i indi-
     *                      cates a row interchange was not required.
     * @param[in,out] B     On entry, the right hand side matrix B.  On exit, the solution matrix X.
     * @param[in] LDB       The leading dimension of the array B.  LDB >= max(1,N).
     * @param[out] INFO     Exit code
     *                      = 0:  successful exit
     *                      < 0:  if INFO = -i, the i-th argument had an illegal value
     */
    static void gbtrs( char TRANS, int N, int KL, int KU, int NRHS, const TYPE *AB, int LDAB,
        const int *IPIV, TYPE *B, int LDB, int &INFO );

    /*!
     * \brief   getri
     * \details  getri computes the inverse of a matrix using the LU factorization
     *    computed by DGETRF.
     * @param[in]  N        The order of the matrix A.  N >= 0.
     * @param[in,out] A     On entry, the factors L and U from the factorization
     *                      A = P*L*U as computed by DGETRF.
     *                      On exit, if INFO = 0, the inverse of the original matrix A.
     * @param[in]  LDA      The leading dimension of the array A.  LDA >= max(1,N).
     * @param[in]  IPIV     The pivot indices from DGETRF; for 1<=i<=N, row i of the
     *                      matrix was interchanged with row IPIV(i).
     * @param[in,out] WORK  Array, dimension (MAX(1,LWORK))
     *                      On exit, if INFO=0, then WORK(1) returns the optimal LWORK.
     * @param[in] LWORK     The dimension of the array WORK.  LWORK >= max(1,N).
     *                      For optimal performance LWORK >= N*NB, where NB is
     *                      the optimal blocksize returned by ILAENV.
     *                      If LWORK = -1, then a workspace query is assumed; the routine
     *                      only calculates the optimal size of the WORK array, returns
     *                      this value as the first entry of the WORK array, and no error
     *                      message related to LWORK is issued by XERBLA.
     * @param[out] INFO     Exit code
     *                      = 0:  successful exit
     *                      < 0:  if INFO = -i, the i-th argument had an illegal value
     *                      > 0:  if INFO = i, U(i,i) is exactly zero; the matrix is
     */
    static void getri( int N, TYPE *A, int LDA, const int *IPIV, TYPE *WORK, int LWORK, int &INFO );

    /*!
     * \brief   trsm
     * \details  trsm solves one of the matrix equations
     *     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
     *  where alpha is a scalar, X and B are m by n matrices, A is a unit, or
     *  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
     *     op( A ) = A   or   op( A ) = A'.
     *  The matrix X is overwritten on B.
     * @param[in] SIDE      Specifies whether op( A ) appears on the left or right of X as follows:
     *                      'L' or 'l'   op( A )*X = alpha*B
     *                      'R' or 'r'   X*op( A ) = alpha*B
     * @param[in] UPLO      Specifies whether the matrix A is an upper or lower triangular matrix as
     * follows:
     *                      'U' or 'u'   A is an upper triangular matrix
     *                      'L' or 'l'   A is a lower triangular matrix
     * @param[in] TRANS     The operation to be performed as follows:
     *                      'N' or 'n'   y := alpha*A*x + beta*y
     *                      'T' or 't'   y := alpha*A'*x + beta*y
     *                      'C' or 'c'   y := alpha*A'*x + beta*y
     * @param[in] DIAG      Is A unit triangular as follows
     *                      'U' or 'u'   A is assumed to be unit triangular
     *                      'N' or 'n'   A is not assumed to be unit triangular
     * @param[in] M         Specifies the number of rows of B. M must be at least zero.
     * @param[in] N         Specifies the number of columns of B.  N must be at least zero.
     * @param[in] ALPHA     Specifies the scalar alpha.  When alpha is zero then A is not
     *                      referenced and B need not be set before entry.
     * @param[in] A         Input array of DIMENSION ( LDA, k ), where k is m when SIDE = 'L' or 'l'
     *                      and is n when SIDE = 'R' or 'r'.
     *                      When UPLO = 'U' or 'u', the leading k by k upper triangular part
     *                          of the array A must contain the upper triangular matrix and the
     *                          strictly lower triangular part of A is not referenced.
     *                      When UPLO = 'L' or 'l', the leading k by k lower triangular part
     *                          of the array A must contain the lower triangular matrix and the
     *                          strictly upper triangular part of A is not referenced.
     *                      Note that when DIAG = 'U' or 'u', the diagonal elements of
     *                          A are not referenced either, but are assumed to be unity.
     * @param[in] LDA       The first dimension of A as declared in the calling (sub) program.
     *                      When SIDE = 'L' or 'l' then LDA must be at least max(1,m),
     *                      when SIDE = 'R' or 'r' then LDA must be at least max(1,n).
     * @param[in,out] B     On entry, the leading m by n part of the array B must contain the
     * right-hand
     *                      side matrix B, and on exit is overwritten by the solution matrix X.
     * @param[in] LDB       The first dimension of B as declared in the calling (sub) program.
     *                      LDB must be at least max(1,m).
     */
    static void trsm( char SIDE, char UPLO, char TRANS, char DIAG, int M, int N, TYPE ALPHA,
        const TYPE *A, int LDA, TYPE *B, int LDB );

    /*!
     * \brief   lamch
     * \details  lamch determines TYPE precision machine parameters.
     *
     * @param[in] cmach     Specifies the value to be returned by DLAMCH:
     *                      'E' or 'e':   eps   = relative machine precision
     *                      'S' or 's :   sfmin = safe minimum, such that 1/sfmin does not overflow
     *                      'B' or 'b':   base  = base of the machine
     *                      'P' or 'p':   prec  = eps*base
     *                      'N' or 'n':   TYPE     = number of (base) digits in the mantissa
     *                      'R' or 'r':   rnd   = 1.0 when rounding occurs in addition, 0.0
     * otherwise
     *                      'M' or 'm':   emin  = minimum exponent before (gradual) underflow
     *                      'U' or 'u':   rmin  = underflow threshold - base**(emin-1)
     *                      'L' or 'l':   emax  = largest exponent before overflow
     *                      'O' or 'o':   rmax  = overflow threshold  - (base**emax)*(1-eps)
     * @return              Return the requested value
     */
    static TYPE lamch( char cmach );

    /*!
     * \brief  Fill a vector with random numbers
     * \details  This will fill a vector with random numbers in the range (0,1)
     * @param[in] N         Specifies the number of entries to fill
     * @param[out] X        Output vector
     */
    static void random( int N, TYPE *X );


    /*!
     * \brief   Run a test for a given LAPACK/BLAS routine
     * \details  This will run a simple test for a given LAPACK/BLAS routine.
     *      While this only runs some simple tests, it should detect basic errors.
     * @param[in] routine   The routine to test
     * @param[in] N         The number of times to repeat the test (useful for thread-safety
     * testing)
     * @param[out] error    The largest error detected
     * @return              Return true if test passed
     */
    static bool run_test( const std::string &routine, int N, double &error );

    //! List all of the availible tests to run
    static std::vector<std::string> list_all_tests();

    /*!
     * \brief   Run the basic test suite
     * \details  This will run all the simple tests
     * @return              The number of failures detected
     */
    static int run_all_test();

    //! Return machine parameter info
    static LapackMachineParams machineParams();

    //! Print lapack info
    static std::string info();

    //! Set the number of threads to use for BLAS/LAPACK
    static int set_num_threads( int N );
};


#endif
