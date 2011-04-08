#include "rb_lapack.h"

extern VOID cppsvx_(char *fact, char *uplo, integer *n, integer *nrhs, complex *ap, complex *afp, char *equed, real *s, complex *b, integer *ldb, complex *x, integer *ldx, real *rcond, real *ferr, real *berr, complex *work, real *rwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_cppsvx(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_fact;
  char fact; 
  VALUE rblapack_uplo;
  char uplo; 
  VALUE rblapack_ap;
  complex *ap; 
  VALUE rblapack_afp;
  complex *afp; 
  VALUE rblapack_equed;
  char equed; 
  VALUE rblapack_s;
  real *s; 
  VALUE rblapack_b;
  complex *b; 
  VALUE rblapack_x;
  complex *x; 
  VALUE rblapack_rcond;
  real rcond; 
  VALUE rblapack_ferr;
  real *ferr; 
  VALUE rblapack_berr;
  real *berr; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_ap_out__;
  complex *ap_out__;
  VALUE rblapack_afp_out__;
  complex *afp_out__;
  VALUE rblapack_s_out__;
  real *s_out__;
  VALUE rblapack_b_out__;
  complex *b_out__;
  complex *work;
  real *rwork;

  integer n;
  integer ldb;
  integer nrhs;
  integer ldx;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  x, rcond, ferr, berr, info, ap, afp, equed, s, b = NumRu::Lapack.cppsvx( fact, uplo, ap, afp, equed, s, b, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE CPPSVX( FACT, UPLO, N, NRHS, AP, AFP, EQUED, S, B, LDB, X, LDX, RCOND, FERR, BERR, WORK, RWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  CPPSVX uses the Cholesky factorization A = U**H*U or A = L*L**H to\n*  compute the solution to a complex system of linear equations\n*     A * X = B,\n*  where A is an N-by-N Hermitian positive definite matrix stored in\n*  packed format and X and B are N-by-NRHS matrices.\n*\n*  Error bounds on the solution and a condition estimate are also\n*  provided.\n*\n*  Description\n*  ===========\n*\n*  The following steps are performed:\n*\n*  1. If FACT = 'E', real scaling factors are computed to equilibrate\n*     the system:\n*        diag(S) * A * diag(S) * inv(diag(S)) * X = diag(S) * B\n*     Whether or not the system will be equilibrated depends on the\n*     scaling of the matrix A, but if equilibration is used, A is\n*     overwritten by diag(S)*A*diag(S) and B by diag(S)*B.\n*\n*  2. If FACT = 'N' or 'E', the Cholesky decomposition is used to\n*     factor the matrix A (after equilibration if FACT = 'E') as\n*        A = U'* U ,  if UPLO = 'U', or\n*        A = L * L',  if UPLO = 'L',\n*     where U is an upper triangular matrix, L is a lower triangular\n*     matrix, and ' indicates conjugate transpose.\n*\n*  3. If the leading i-by-i principal minor is not positive definite,\n*     then the routine returns with INFO = i. Otherwise, the factored\n*     form of A is used to estimate the condition number of the matrix\n*     A.  If the reciprocal of the condition number is less than machine\n*     precision, INFO = N+1 is returned as a warning, but the routine\n*     still goes on to solve for X and compute error bounds as\n*     described below.\n*\n*  4. The system of equations is solved for X using the factored form\n*     of A.\n*\n*  5. Iterative refinement is applied to improve the computed solution\n*     matrix and calculate error bounds and backward error estimates\n*     for it.\n*\n*  6. If equilibration was used, the matrix X is premultiplied by\n*     diag(S) so that it solves the original system before\n*     equilibration.\n*\n\n*  Arguments\n*  =========\n*\n*  FACT    (input) CHARACTER*1\n*          Specifies whether or not the factored form of the matrix A is\n*          supplied on entry, and if not, whether the matrix A should be\n*          equilibrated before it is factored.\n*          = 'F':  On entry, AFP contains the factored form of A.\n*                  If EQUED = 'Y', the matrix A has been equilibrated\n*                  with scaling factors given by S.  AP and AFP will not\n*                  be modified.\n*          = 'N':  The matrix A will be copied to AFP and factored.\n*          = 'E':  The matrix A will be equilibrated if necessary, then\n*                  copied to AFP and factored.\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U':  Upper triangle of A is stored;\n*          = 'L':  Lower triangle of A is stored.\n*\n*  N       (input) INTEGER\n*          The number of linear equations, i.e., the order of the\n*          matrix A.  N >= 0.\n*\n*  NRHS    (input) INTEGER\n*          The number of right hand sides, i.e., the number of columns\n*          of the matrices B and X.  NRHS >= 0.\n*\n*  AP      (input/output) COMPLEX array, dimension (N*(N+1)/2)\n*          On entry, the upper or lower triangle of the Hermitian matrix\n*          A, packed columnwise in a linear array, except if FACT = 'F'\n*          and EQUED = 'Y', then A must contain the equilibrated matrix\n*          diag(S)*A*diag(S).  The j-th column of A is stored in the\n*          array AP as follows:\n*          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;\n*          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.\n*          See below for further details.  A is not modified if\n*          FACT = 'F' or 'N', or if FACT = 'E' and EQUED = 'N' on exit.\n*\n*          On exit, if FACT = 'E' and EQUED = 'Y', A is overwritten by\n*          diag(S)*A*diag(S).\n*\n*  AFP     (input or output) COMPLEX array, dimension (N*(N+1)/2)\n*          If FACT = 'F', then AFP is an input argument and on entry\n*          contains the triangular factor U or L from the Cholesky\n*          factorization A = U**H*U or A = L*L**H, in the same storage\n*          format as A.  If EQUED .ne. 'N', then AFP is the factored\n*          form of the equilibrated matrix A.\n*\n*          If FACT = 'N', then AFP is an output argument and on exit\n*          returns the triangular factor U or L from the Cholesky\n*          factorization A = U**H*U or A = L*L**H of the original\n*          matrix A.\n*\n*          If FACT = 'E', then AFP is an output argument and on exit\n*          returns the triangular factor U or L from the Cholesky\n*          factorization A = U**H*U or A = L*L**H of the equilibrated\n*          matrix A (see the description of AP for the form of the\n*          equilibrated matrix).\n*\n*  EQUED   (input or output) CHARACTER*1\n*          Specifies the form of equilibration that was done.\n*          = 'N':  No equilibration (always true if FACT = 'N').\n*          = 'Y':  Equilibration was done, i.e., A has been replaced by\n*                  diag(S) * A * diag(S).\n*          EQUED is an input argument if FACT = 'F'; otherwise, it is an\n*          output argument.\n*\n*  S       (input or output) REAL array, dimension (N)\n*          The scale factors for A; not accessed if EQUED = 'N'.  S is\n*          an input argument if FACT = 'F'; otherwise, S is an output\n*          argument.  If FACT = 'F' and EQUED = 'Y', each element of S\n*          must be positive.\n*\n*  B       (input/output) COMPLEX array, dimension (LDB,NRHS)\n*          On entry, the N-by-NRHS right hand side matrix B.\n*          On exit, if EQUED = 'N', B is not modified; if EQUED = 'Y',\n*          B is overwritten by diag(S) * B.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B.  LDB >= max(1,N).\n*\n*  X       (output) COMPLEX array, dimension (LDX,NRHS)\n*          If INFO = 0 or INFO = N+1, the N-by-NRHS solution matrix X to\n*          the original system of equations.  Note that if EQUED = 'Y',\n*          A and B are modified on exit, and the solution to the\n*          equilibrated system is inv(diag(S))*X.\n*\n*  LDX     (input) INTEGER\n*          The leading dimension of the array X.  LDX >= max(1,N).\n*\n*  RCOND   (output) REAL\n*          The estimate of the reciprocal condition number of the matrix\n*          A after equilibration (if done).  If RCOND is less than the\n*          machine precision (in particular, if RCOND = 0), the matrix\n*          is singular to working precision.  This condition is\n*          indicated by a return code of INFO > 0.\n*\n*  FERR    (output) REAL array, dimension (NRHS)\n*          The estimated forward error bound for each solution vector\n*          X(j) (the j-th column of the solution matrix X).\n*          If XTRUE is the true solution corresponding to X(j), FERR(j)\n*          is an estimated upper bound for the magnitude of the largest\n*          element in (X(j) - XTRUE) divided by the magnitude of the\n*          largest element in X(j).  The estimate is as reliable as\n*          the estimate for RCOND, and is almost always a slight\n*          overestimate of the true error.\n*\n*  BERR    (output) REAL array, dimension (NRHS)\n*          The componentwise relative backward error of each solution\n*          vector X(j) (i.e., the smallest relative change in\n*          any element of A or B that makes X(j) an exact solution).\n*\n*  WORK    (workspace) COMPLEX array, dimension (2*N)\n*\n*  RWORK   (workspace) REAL array, dimension (N)\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*          > 0:  if INFO = i, and i is\n*                <= N:  the leading minor of order i of A is\n*                       not positive definite, so the factorization\n*                       could not be completed, and the solution has not\n*                       been computed. RCOND = 0 is returned.\n*                = N+1: U is nonsingular, but RCOND is less than machine\n*                       precision, meaning that the matrix is singular\n*                       to working precision.  Nevertheless, the\n*                       solution and error bounds are computed because\n*                       there are a number of situations where the\n*                       computed solution can be more accurate than the\n*                       value of RCOND would suggest.\n*\n\n*  Further Details\n*  ===============\n*\n*  The packed storage scheme is illustrated by the following example\n*  when N = 4, UPLO = 'U':\n*\n*  Two-dimensional storage of the Hermitian matrix A:\n*\n*     a11 a12 a13 a14\n*         a22 a23 a24\n*             a33 a34     (aij = conjg(aji))\n*                 a44\n*\n*  Packed storage of the upper triangle of A:\n*\n*  AP = [ a11, a12, a22, a13, a23, a33, a14, a24, a34, a44 ]\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  x, rcond, ferr, berr, info, ap, afp, equed, s, b = NumRu::Lapack.cppsvx( fact, uplo, ap, afp, equed, s, b, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 7)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 7)", argc);
  rblapack_fact = argv[0];
  rblapack_uplo = argv[1];
  rblapack_ap = argv[2];
  rblapack_afp = argv[3];
  rblapack_equed = argv[4];
  rblapack_s = argv[5];
  rblapack_b = argv[6];
  if (rb_options != Qnil) {
  }

  uplo = StringValueCStr(rblapack_uplo)[0];
  equed = StringValueCStr(rblapack_equed)[0];
  if (!NA_IsNArray(rblapack_b))
    rb_raise(rb_eArgError, "b (7th argument) must be NArray");
  if (NA_RANK(rblapack_b) != 2)
    rb_raise(rb_eArgError, "rank of b (7th argument) must be %d", 2);
  nrhs = NA_SHAPE1(rblapack_b);
  ldb = NA_SHAPE0(rblapack_b);
  if (NA_TYPE(rblapack_b) != NA_SCOMPLEX)
    rblapack_b = na_change_type(rblapack_b, NA_SCOMPLEX);
  b = NA_PTR_TYPE(rblapack_b, complex*);
  fact = StringValueCStr(rblapack_fact)[0];
  if (!NA_IsNArray(rblapack_s))
    rb_raise(rb_eArgError, "s (6th argument) must be NArray");
  if (NA_RANK(rblapack_s) != 1)
    rb_raise(rb_eArgError, "rank of s (6th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_s);
  if (NA_TYPE(rblapack_s) != NA_SFLOAT)
    rblapack_s = na_change_type(rblapack_s, NA_SFLOAT);
  s = NA_PTR_TYPE(rblapack_s, real*);
  if (!NA_IsNArray(rblapack_ap))
    rb_raise(rb_eArgError, "ap (3th argument) must be NArray");
  if (NA_RANK(rblapack_ap) != 1)
    rb_raise(rb_eArgError, "rank of ap (3th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_ap) != (n*(n+1)/2))
    rb_raise(rb_eRuntimeError, "shape 0 of ap must be %d", n*(n+1)/2);
  if (NA_TYPE(rblapack_ap) != NA_SCOMPLEX)
    rblapack_ap = na_change_type(rblapack_ap, NA_SCOMPLEX);
  ap = NA_PTR_TYPE(rblapack_ap, complex*);
  if (!NA_IsNArray(rblapack_afp))
    rb_raise(rb_eArgError, "afp (4th argument) must be NArray");
  if (NA_RANK(rblapack_afp) != 1)
    rb_raise(rb_eArgError, "rank of afp (4th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_afp) != (n*(n+1)/2))
    rb_raise(rb_eRuntimeError, "shape 0 of afp must be %d", n*(n+1)/2);
  if (NA_TYPE(rblapack_afp) != NA_SCOMPLEX)
    rblapack_afp = na_change_type(rblapack_afp, NA_SCOMPLEX);
  afp = NA_PTR_TYPE(rblapack_afp, complex*);
  ldx = MAX(1,n);
  {
    int shape[2];
    shape[0] = ldx;
    shape[1] = nrhs;
    rblapack_x = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  x = NA_PTR_TYPE(rblapack_x, complex*);
  {
    int shape[1];
    shape[0] = nrhs;
    rblapack_ferr = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  ferr = NA_PTR_TYPE(rblapack_ferr, real*);
  {
    int shape[1];
    shape[0] = nrhs;
    rblapack_berr = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  berr = NA_PTR_TYPE(rblapack_berr, real*);
  {
    int shape[1];
    shape[0] = n*(n+1)/2;
    rblapack_ap_out__ = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  ap_out__ = NA_PTR_TYPE(rblapack_ap_out__, complex*);
  MEMCPY(ap_out__, ap, complex, NA_TOTAL(rblapack_ap));
  rblapack_ap = rblapack_ap_out__;
  ap = ap_out__;
  {
    int shape[1];
    shape[0] = n*(n+1)/2;
    rblapack_afp_out__ = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  afp_out__ = NA_PTR_TYPE(rblapack_afp_out__, complex*);
  MEMCPY(afp_out__, afp, complex, NA_TOTAL(rblapack_afp));
  rblapack_afp = rblapack_afp_out__;
  afp = afp_out__;
  {
    int shape[1];
    shape[0] = n;
    rblapack_s_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  s_out__ = NA_PTR_TYPE(rblapack_s_out__, real*);
  MEMCPY(s_out__, s, real, NA_TOTAL(rblapack_s));
  rblapack_s = rblapack_s_out__;
  s = s_out__;
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = nrhs;
    rblapack_b_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rblapack_b_out__, complex*);
  MEMCPY(b_out__, b, complex, NA_TOTAL(rblapack_b));
  rblapack_b = rblapack_b_out__;
  b = b_out__;
  work = ALLOC_N(complex, (2*n));
  rwork = ALLOC_N(real, (n));

  cppsvx_(&fact, &uplo, &n, &nrhs, ap, afp, &equed, s, b, &ldb, x, &ldx, &rcond, ferr, berr, work, rwork, &info);

  free(work);
  free(rwork);
  rblapack_rcond = rb_float_new((double)rcond);
  rblapack_info = INT2NUM(info);
  rblapack_equed = rb_str_new(&equed,1);
  return rb_ary_new3(10, rblapack_x, rblapack_rcond, rblapack_ferr, rblapack_berr, rblapack_info, rblapack_ap, rblapack_afp, rblapack_equed, rblapack_s, rblapack_b);
}

void
init_lapack_cppsvx(VALUE mLapack){
  rb_define_module_function(mLapack, "cppsvx", rblapack_cppsvx, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
