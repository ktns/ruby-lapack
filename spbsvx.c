#include "rb_lapack.h"

static VALUE
rb_spbsvx(int argc, VALUE *argv, VALUE self){
  VALUE rb_fact;
  char fact; 
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_kd;
  integer kd; 
  VALUE rb_ab;
  real *ab; 
  VALUE rb_afb;
  real *afb; 
  VALUE rb_equed;
  char equed; 
  VALUE rb_s;
  real *s; 
  VALUE rb_b;
  real *b; 
  VALUE rb_x;
  real *x; 
  VALUE rb_rcond;
  real rcond; 
  VALUE rb_ferr;
  real *ferr; 
  VALUE rb_berr;
  real *berr; 
  VALUE rb_info;
  integer info; 
  VALUE rb_ab_out__;
  real *ab_out__;
  VALUE rb_afb_out__;
  real *afb_out__;
  VALUE rb_s_out__;
  real *s_out__;
  VALUE rb_b_out__;
  real *b_out__;
  real *work;
  integer *iwork;

  integer ldab;
  integer n;
  integer ldafb;
  integer ldb;
  integer nrhs;
  integer ldx;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  x, rcond, ferr, berr, info, ab, afb, equed, s, b = NumRu::Lapack.spbsvx( fact, uplo, kd, ab, afb, equed, s, b)\n    or\n  NumRu::Lapack.spbsvx  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SPBSVX( FACT, UPLO, N, KD, NRHS, AB, LDAB, AFB, LDAFB, EQUED, S, B, LDB, X, LDX, RCOND, FERR, BERR, WORK, IWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  SPBSVX uses the Cholesky factorization A = U**T*U or A = L*L**T to\n*  compute the solution to a real system of linear equations\n*     A * X = B,\n*  where A is an N-by-N symmetric positive definite band matrix and X\n*  and B are N-by-NRHS matrices.\n*\n*  Error bounds on the solution and a condition estimate are also\n*  provided.\n*\n*  Description\n*  ===========\n*\n*  The following steps are performed:\n*\n*  1. If FACT = 'E', real scaling factors are computed to equilibrate\n*     the system:\n*        diag(S) * A * diag(S) * inv(diag(S)) * X = diag(S) * B\n*     Whether or not the system will be equilibrated depends on the\n*     scaling of the matrix A, but if equilibration is used, A is\n*     overwritten by diag(S)*A*diag(S) and B by diag(S)*B.\n*\n*  2. If FACT = 'N' or 'E', the Cholesky decomposition is used to\n*     factor the matrix A (after equilibration if FACT = 'E') as\n*        A = U**T * U,  if UPLO = 'U', or\n*        A = L * L**T,  if UPLO = 'L',\n*     where U is an upper triangular band matrix, and L is a lower\n*     triangular band matrix.\n*\n*  3. If the leading i-by-i principal minor is not positive definite,\n*     then the routine returns with INFO = i. Otherwise, the factored\n*     form of A is used to estimate the condition number of the matrix\n*     A.  If the reciprocal of the condition number is less than machine\n*     precision, INFO = N+1 is returned as a warning, but the routine\n*     still goes on to solve for X and compute error bounds as\n*     described below.\n*\n*  4. The system of equations is solved for X using the factored form\n*     of A.\n*\n*  5. Iterative refinement is applied to improve the computed solution\n*     matrix and calculate error bounds and backward error estimates\n*     for it.\n*\n*  6. If equilibration was used, the matrix X is premultiplied by\n*     diag(S) so that it solves the original system before\n*     equilibration.\n*\n\n*  Arguments\n*  =========\n*\n*  FACT    (input) CHARACTER*1\n*          Specifies whether or not the factored form of the matrix A is\n*          supplied on entry, and if not, whether the matrix A should be\n*          equilibrated before it is factored.\n*          = 'F':  On entry, AFB contains the factored form of A.\n*                  If EQUED = 'Y', the matrix A has been equilibrated\n*                  with scaling factors given by S.  AB and AFB will not\n*                  be modified.\n*          = 'N':  The matrix A will be copied to AFB and factored.\n*          = 'E':  The matrix A will be equilibrated if necessary, then\n*                  copied to AFB and factored.\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U':  Upper triangle of A is stored;\n*          = 'L':  Lower triangle of A is stored.\n*\n*  N       (input) INTEGER\n*          The number of linear equations, i.e., the order of the\n*          matrix A.  N >= 0.\n*\n*  KD      (input) INTEGER\n*          The number of superdiagonals of the matrix A if UPLO = 'U',\n*          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.\n*\n*  NRHS    (input) INTEGER\n*          The number of right-hand sides, i.e., the number of columns\n*          of the matrices B and X.  NRHS >= 0.\n*\n*  AB      (input/output) REAL array, dimension (LDAB,N)\n*          On entry, the upper or lower triangle of the symmetric band\n*          matrix A, stored in the first KD+1 rows of the array, except\n*          if FACT = 'F' and EQUED = 'Y', then A must contain the\n*          equilibrated matrix diag(S)*A*diag(S).  The j-th column of A\n*          is stored in the j-th column of the array AB as follows:\n*          if UPLO = 'U', AB(KD+1+i-j,j) = A(i,j) for max(1,j-KD)<=i<=j;\n*          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(N,j+KD).\n*          See below for further details.\n*\n*          On exit, if FACT = 'E' and EQUED = 'Y', A is overwritten by\n*          diag(S)*A*diag(S).\n*\n*  LDAB    (input) INTEGER\n*          The leading dimension of the array A.  LDAB >= KD+1.\n*\n*  AFB     (input or output) REAL array, dimension (LDAFB,N)\n*          If FACT = 'F', then AFB is an input argument and on entry\n*          contains the triangular factor U or L from the Cholesky\n*          factorization A = U**T*U or A = L*L**T of the band matrix\n*          A, in the same storage format as A (see AB).  If EQUED = 'Y',\n*          then AFB is the factored form of the equilibrated matrix A.\n*\n*          If FACT = 'N', then AFB is an output argument and on exit\n*          returns the triangular factor U or L from the Cholesky\n*          factorization A = U**T*U or A = L*L**T.\n*\n*          If FACT = 'E', then AFB is an output argument and on exit\n*          returns the triangular factor U or L from the Cholesky\n*          factorization A = U**T*U or A = L*L**T of the equilibrated\n*          matrix A (see the description of A for the form of the\n*          equilibrated matrix).\n*\n*  LDAFB   (input) INTEGER\n*          The leading dimension of the array AFB.  LDAFB >= KD+1.\n*\n*  EQUED   (input or output) CHARACTER*1\n*          Specifies the form of equilibration that was done.\n*          = 'N':  No equilibration (always true if FACT = 'N').\n*          = 'Y':  Equilibration was done, i.e., A has been replaced by\n*                  diag(S) * A * diag(S).\n*          EQUED is an input argument if FACT = 'F'; otherwise, it is an\n*          output argument.\n*\n*  S       (input or output) REAL array, dimension (N)\n*          The scale factors for A; not accessed if EQUED = 'N'.  S is\n*          an input argument if FACT = 'F'; otherwise, S is an output\n*          argument.  If FACT = 'F' and EQUED = 'Y', each element of S\n*          must be positive.\n*\n*  B       (input/output) REAL array, dimension (LDB,NRHS)\n*          On entry, the N-by-NRHS right hand side matrix B.\n*          On exit, if EQUED = 'N', B is not modified; if EQUED = 'Y',\n*          B is overwritten by diag(S) * B.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B.  LDB >= max(1,N).\n*\n*  X       (output) REAL array, dimension (LDX,NRHS)\n*          If INFO = 0 or INFO = N+1, the N-by-NRHS solution matrix X to\n*          the original system of equations.  Note that if EQUED = 'Y',\n*          A and B are modified on exit, and the solution to the\n*          equilibrated system is inv(diag(S))*X.\n*\n*  LDX     (input) INTEGER\n*          The leading dimension of the array X.  LDX >= max(1,N).\n*\n*  RCOND   (output) REAL\n*          The estimate of the reciprocal condition number of the matrix\n*          A after equilibration (if done).  If RCOND is less than the\n*          machine precision (in particular, if RCOND = 0), the matrix\n*          is singular to working precision.  This condition is\n*          indicated by a return code of INFO > 0.\n*\n*  FERR    (output) REAL array, dimension (NRHS)\n*          The estimated forward error bound for each solution vector\n*          X(j) (the j-th column of the solution matrix X).\n*          If XTRUE is the true solution corresponding to X(j), FERR(j)\n*          is an estimated upper bound for the magnitude of the largest\n*          element in (X(j) - XTRUE) divided by the magnitude of the\n*          largest element in X(j).  The estimate is as reliable as\n*          the estimate for RCOND, and is almost always a slight\n*          overestimate of the true error.\n*\n*  BERR    (output) REAL array, dimension (NRHS)\n*          The componentwise relative backward error of each solution\n*          vector X(j) (i.e., the smallest relative change in\n*          any element of A or B that makes X(j) an exact solution).\n*\n*  WORK    (workspace) REAL array, dimension (3*N)\n*\n*  IWORK   (workspace) INTEGER array, dimension (N)\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*          > 0:  if INFO = i, and i is\n*                <= N:  the leading minor of order i of A is\n*                       not positive definite, so the factorization\n*                       could not be completed, and the solution has not\n*                       been computed. RCOND = 0 is returned.\n*                = N+1: U is nonsingular, but RCOND is less than machine\n*                       precision, meaning that the matrix is singular\n*                       to working precision.  Nevertheless, the\n*                       solution and error bounds are computed because\n*                       there are a number of situations where the\n*                       computed solution can be more accurate than the\n*                       value of RCOND would suggest.\n*\n\n*  Further Details\n*  ===============\n*\n*  The band storage scheme is illustrated by the following example, when\n*  N = 6, KD = 2, and UPLO = 'U':\n*\n*  Two-dimensional storage of the symmetric matrix A:\n*\n*     a11  a12  a13\n*          a22  a23  a24\n*               a33  a34  a35\n*                    a44  a45  a46\n*                         a55  a56\n*     (aij=conjg(aji))         a66\n*\n*  Band storage of the upper triangle of A:\n*\n*      *    *   a13  a24  a35  a46\n*      *   a12  a23  a34  a45  a56\n*     a11  a22  a33  a44  a55  a66\n*\n*  Similarly, if UPLO = 'L' the format of A is as follows:\n*\n*     a11  a22  a33  a44  a55  a66\n*     a21  a32  a43  a54  a65   *\n*     a31  a42  a53  a64   *    *\n*\n*  Array elements marked * are not used by the routine.\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 8)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 8)", argc);
  rb_fact = argv[0];
  rb_uplo = argv[1];
  rb_kd = argv[2];
  rb_ab = argv[3];
  rb_afb = argv[4];
  rb_equed = argv[5];
  rb_s = argv[6];
  rb_b = argv[7];

  fact = StringValueCStr(rb_fact)[0];
  uplo = StringValueCStr(rb_uplo)[0];
  kd = NUM2INT(rb_kd);
  equed = StringValueCStr(rb_equed)[0];
  if (!NA_IsNArray(rb_ab))
    rb_raise(rb_eArgError, "ab (4th argument) must be NArray");
  if (NA_RANK(rb_ab) != 2)
    rb_raise(rb_eArgError, "rank of ab (4th argument) must be %d", 2);
  ldab = NA_SHAPE0(rb_ab);
  n = NA_SHAPE1(rb_ab);
  if (NA_TYPE(rb_ab) != NA_SFLOAT)
    rb_ab = na_change_type(rb_ab, NA_SFLOAT);
  ab = NA_PTR_TYPE(rb_ab, real*);
  if (!NA_IsNArray(rb_afb))
    rb_raise(rb_eArgError, "afb (5th argument) must be NArray");
  if (NA_RANK(rb_afb) != 2)
    rb_raise(rb_eArgError, "rank of afb (5th argument) must be %d", 2);
  ldafb = NA_SHAPE0(rb_afb);
  if (NA_SHAPE1(rb_afb) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of afb must be the same as shape 1 of ab");
  if (NA_TYPE(rb_afb) != NA_SFLOAT)
    rb_afb = na_change_type(rb_afb, NA_SFLOAT);
  afb = NA_PTR_TYPE(rb_afb, real*);
  if (!NA_IsNArray(rb_s))
    rb_raise(rb_eArgError, "s (7th argument) must be NArray");
  if (NA_RANK(rb_s) != 1)
    rb_raise(rb_eArgError, "rank of s (7th argument) must be %d", 1);
  if (NA_SHAPE0(rb_s) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of s must be the same as shape 1 of ab");
  if (NA_TYPE(rb_s) != NA_SFLOAT)
    rb_s = na_change_type(rb_s, NA_SFLOAT);
  s = NA_PTR_TYPE(rb_s, real*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (8th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (8th argument) must be %d", 2);
  ldb = NA_SHAPE0(rb_b);
  nrhs = NA_SHAPE1(rb_b);
  if (NA_TYPE(rb_b) != NA_SFLOAT)
    rb_b = na_change_type(rb_b, NA_SFLOAT);
  b = NA_PTR_TYPE(rb_b, real*);
  ldx = MAX(1,n);
  {
    int shape[2];
    shape[0] = ldx;
    shape[1] = nrhs;
    rb_x = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  x = NA_PTR_TYPE(rb_x, real*);
  {
    int shape[1];
    shape[0] = nrhs;
    rb_ferr = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  ferr = NA_PTR_TYPE(rb_ferr, real*);
  {
    int shape[1];
    shape[0] = nrhs;
    rb_berr = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  berr = NA_PTR_TYPE(rb_berr, real*);
  {
    int shape[2];
    shape[0] = ldab;
    shape[1] = n;
    rb_ab_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  ab_out__ = NA_PTR_TYPE(rb_ab_out__, real*);
  MEMCPY(ab_out__, ab, real, NA_TOTAL(rb_ab));
  rb_ab = rb_ab_out__;
  ab = ab_out__;
  {
    int shape[2];
    shape[0] = ldafb;
    shape[1] = n;
    rb_afb_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  afb_out__ = NA_PTR_TYPE(rb_afb_out__, real*);
  MEMCPY(afb_out__, afb, real, NA_TOTAL(rb_afb));
  rb_afb = rb_afb_out__;
  afb = afb_out__;
  {
    int shape[1];
    shape[0] = n;
    rb_s_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  s_out__ = NA_PTR_TYPE(rb_s_out__, real*);
  MEMCPY(s_out__, s, real, NA_TOTAL(rb_s));
  rb_s = rb_s_out__;
  s = s_out__;
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = nrhs;
    rb_b_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rb_b_out__, real*);
  MEMCPY(b_out__, b, real, NA_TOTAL(rb_b));
  rb_b = rb_b_out__;
  b = b_out__;
  work = ALLOC_N(real, (3*n));
  iwork = ALLOC_N(integer, (n));

  spbsvx_(&fact, &uplo, &n, &kd, &nrhs, ab, &ldab, afb, &ldafb, &equed, s, b, &ldb, x, &ldx, &rcond, ferr, berr, work, iwork, &info);

  free(work);
  free(iwork);
  rb_rcond = rb_float_new((double)rcond);
  rb_info = INT2NUM(info);
  rb_equed = rb_str_new(&equed,1);
  return rb_ary_new3(10, rb_x, rb_rcond, rb_ferr, rb_berr, rb_info, rb_ab, rb_afb, rb_equed, rb_s, rb_b);
}

void
init_lapack_spbsvx(VALUE mLapack){
  rb_define_module_function(mLapack, "spbsvx", rb_spbsvx, -1);
}
