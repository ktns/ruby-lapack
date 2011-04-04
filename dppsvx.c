#include "rb_lapack.h"

extern VOID dppsvx_(char *fact, char *uplo, integer *n, integer *nrhs, doublereal *ap, doublereal *afp, char *equed, doublereal *s, doublereal *b, integer *ldb, doublereal *x, integer *ldx, doublereal *rcond, doublereal *ferr, doublereal *berr, doublereal *work, integer *iwork, integer *info);

static VALUE
rb_dppsvx(int argc, VALUE *argv, VALUE self){
  VALUE rb_fact;
  char fact; 
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_ap;
  doublereal *ap; 
  VALUE rb_afp;
  doublereal *afp; 
  VALUE rb_equed;
  char equed; 
  VALUE rb_s;
  doublereal *s; 
  VALUE rb_b;
  doublereal *b; 
  VALUE rb_x;
  doublereal *x; 
  VALUE rb_rcond;
  doublereal rcond; 
  VALUE rb_ferr;
  doublereal *ferr; 
  VALUE rb_berr;
  doublereal *berr; 
  VALUE rb_info;
  integer info; 
  VALUE rb_ap_out__;
  doublereal *ap_out__;
  VALUE rb_afp_out__;
  doublereal *afp_out__;
  VALUE rb_s_out__;
  doublereal *s_out__;
  VALUE rb_b_out__;
  doublereal *b_out__;
  doublereal *work;
  integer *iwork;

  integer n;
  integer ldb;
  integer nrhs;
  integer ldx;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  x, rcond, ferr, berr, info, ap, afp, equed, s, b = NumRu::Lapack.dppsvx( fact, uplo, ap, afp, equed, s, b)\n    or\n  NumRu::Lapack.dppsvx  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DPPSVX( FACT, UPLO, N, NRHS, AP, AFP, EQUED, S, B, LDB, X, LDX, RCOND, FERR, BERR, WORK, IWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  DPPSVX uses the Cholesky factorization A = U**T*U or A = L*L**T to\n*  compute the solution to a real system of linear equations\n*     A * X = B,\n*  where A is an N-by-N symmetric positive definite matrix stored in\n*  packed format and X and B are N-by-NRHS matrices.\n*\n*  Error bounds on the solution and a condition estimate are also\n*  provided.\n*\n*  Description\n*  ===========\n*\n*  The following steps are performed:\n*\n*  1. If FACT = 'E', real scaling factors are computed to equilibrate\n*     the system:\n*        diag(S) * A * diag(S) * inv(diag(S)) * X = diag(S) * B\n*     Whether or not the system will be equilibrated depends on the\n*     scaling of the matrix A, but if equilibration is used, A is\n*     overwritten by diag(S)*A*diag(S) and B by diag(S)*B.\n*\n*  2. If FACT = 'N' or 'E', the Cholesky decomposition is used to\n*     factor the matrix A (after equilibration if FACT = 'E') as\n*        A = U**T* U,  if UPLO = 'U', or\n*        A = L * L**T,  if UPLO = 'L',\n*     where U is an upper triangular matrix and L is a lower triangular\n*     matrix.\n*\n*  3. If the leading i-by-i principal minor is not positive definite,\n*     then the routine returns with INFO = i. Otherwise, the factored\n*     form of A is used to estimate the condition number of the matrix\n*     A.  If the reciprocal of the condition number is less than machine\n*     precision, INFO = N+1 is returned as a warning, but the routine\n*     still goes on to solve for X and compute error bounds as\n*     described below.\n*\n*  4. The system of equations is solved for X using the factored form\n*     of A.\n*\n*  5. Iterative refinement is applied to improve the computed solution\n*     matrix and calculate error bounds and backward error estimates\n*     for it.\n*\n*  6. If equilibration was used, the matrix X is premultiplied by\n*     diag(S) so that it solves the original system before\n*     equilibration.\n*\n\n*  Arguments\n*  =========\n*\n*  FACT    (input) CHARACTER*1\n*          Specifies whether or not the factored form of the matrix A is\n*          supplied on entry, and if not, whether the matrix A should be\n*          equilibrated before it is factored.\n*          = 'F':  On entry, AFP contains the factored form of A.\n*                  If EQUED = 'Y', the matrix A has been equilibrated\n*                  with scaling factors given by S.  AP and AFP will not\n*                  be modified.\n*          = 'N':  The matrix A will be copied to AFP and factored.\n*          = 'E':  The matrix A will be equilibrated if necessary, then\n*                  copied to AFP and factored.\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U':  Upper triangle of A is stored;\n*          = 'L':  Lower triangle of A is stored.\n*\n*  N       (input) INTEGER\n*          The number of linear equations, i.e., the order of the\n*          matrix A.  N >= 0.\n*\n*  NRHS    (input) INTEGER\n*          The number of right hand sides, i.e., the number of columns\n*          of the matrices B and X.  NRHS >= 0.\n*\n*  AP      (input/output) DOUBLE PRECISION array, dimension (N*(N+1)/2)\n*          On entry, the upper or lower triangle of the symmetric matrix\n*          A, packed columnwise in a linear array, except if FACT = 'F'\n*          and EQUED = 'Y', then A must contain the equilibrated matrix\n*          diag(S)*A*diag(S).  The j-th column of A is stored in the\n*          array AP as follows:\n*          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;\n*          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.\n*          See below for further details.  A is not modified if\n*          FACT = 'F' or 'N', or if FACT = 'E' and EQUED = 'N' on exit.\n*\n*          On exit, if FACT = 'E' and EQUED = 'Y', A is overwritten by\n*          diag(S)*A*diag(S).\n*\n*  AFP     (input or output) DOUBLE PRECISION array, dimension\n*                            (N*(N+1)/2)\n*          If FACT = 'F', then AFP is an input argument and on entry\n*          contains the triangular factor U or L from the Cholesky\n*          factorization A = U'*U or A = L*L', in the same storage\n*          format as A.  If EQUED .ne. 'N', then AFP is the factored\n*          form of the equilibrated matrix A.\n*\n*          If FACT = 'N', then AFP is an output argument and on exit\n*          returns the triangular factor U or L from the Cholesky\n*          factorization A = U'*U or A = L*L' of the original matrix A.\n*\n*          If FACT = 'E', then AFP is an output argument and on exit\n*          returns the triangular factor U or L from the Cholesky\n*          factorization A = U'*U or A = L*L' of the equilibrated\n*          matrix A (see the description of AP for the form of the\n*          equilibrated matrix).\n*\n*  EQUED   (input or output) CHARACTER*1\n*          Specifies the form of equilibration that was done.\n*          = 'N':  No equilibration (always true if FACT = 'N').\n*          = 'Y':  Equilibration was done, i.e., A has been replaced by\n*                  diag(S) * A * diag(S).\n*          EQUED is an input argument if FACT = 'F'; otherwise, it is an\n*          output argument.\n*\n*  S       (input or output) DOUBLE PRECISION array, dimension (N)\n*          The scale factors for A; not accessed if EQUED = 'N'.  S is\n*          an input argument if FACT = 'F'; otherwise, S is an output\n*          argument.  If FACT = 'F' and EQUED = 'Y', each element of S\n*          must be positive.\n*\n*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)\n*          On entry, the N-by-NRHS right hand side matrix B.\n*          On exit, if EQUED = 'N', B is not modified; if EQUED = 'Y',\n*          B is overwritten by diag(S) * B.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B.  LDB >= max(1,N).\n*\n*  X       (output) DOUBLE PRECISION array, dimension (LDX,NRHS)\n*          If INFO = 0 or INFO = N+1, the N-by-NRHS solution matrix X to\n*          the original system of equations.  Note that if EQUED = 'Y',\n*          A and B are modified on exit, and the solution to the\n*          equilibrated system is inv(diag(S))*X.\n*\n*  LDX     (input) INTEGER\n*          The leading dimension of the array X.  LDX >= max(1,N).\n*\n*  RCOND   (output) DOUBLE PRECISION\n*          The estimate of the reciprocal condition number of the matrix\n*          A after equilibration (if done).  If RCOND is less than the\n*          machine precision (in particular, if RCOND = 0), the matrix\n*          is singular to working precision.  This condition is\n*          indicated by a return code of INFO > 0.\n*\n*  FERR    (output) DOUBLE PRECISION array, dimension (NRHS)\n*          The estimated forward error bound for each solution vector\n*          X(j) (the j-th column of the solution matrix X).\n*          If XTRUE is the true solution corresponding to X(j), FERR(j)\n*          is an estimated upper bound for the magnitude of the largest\n*          element in (X(j) - XTRUE) divided by the magnitude of the\n*          largest element in X(j).  The estimate is as reliable as\n*          the estimate for RCOND, and is almost always a slight\n*          overestimate of the true error.\n*\n*  BERR    (output) DOUBLE PRECISION array, dimension (NRHS)\n*          The componentwise relative backward error of each solution\n*          vector X(j) (i.e., the smallest relative change in\n*          any element of A or B that makes X(j) an exact solution).\n*\n*  WORK    (workspace) DOUBLE PRECISION array, dimension (3*N)\n*\n*  IWORK   (workspace) INTEGER array, dimension (N)\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*          > 0:  if INFO = i, and i is\n*                <= N:  the leading minor of order i of A is\n*                       not positive definite, so the factorization\n*                       could not be completed, and the solution has not\n*                       been computed. RCOND = 0 is returned.\n*                = N+1: U is nonsingular, but RCOND is less than machine\n*                       precision, meaning that the matrix is singular\n*                       to working precision.  Nevertheless, the\n*                       solution and error bounds are computed because\n*                       there are a number of situations where the\n*                       computed solution can be more accurate than the\n*                       value of RCOND would suggest.\n*\n\n*  Further Details\n*  ===============\n*\n*  The packed storage scheme is illustrated by the following example\n*  when N = 4, UPLO = 'U':\n*\n*  Two-dimensional storage of the symmetric matrix A:\n*\n*     a11 a12 a13 a14\n*         a22 a23 a24\n*             a33 a34     (aij = conjg(aji))\n*                 a44\n*\n*  Packed storage of the upper triangle of A:\n*\n*  AP = [ a11, a12, a22, a13, a23, a33, a14, a24, a34, a44 ]\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 7)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 7)", argc);
  rb_fact = argv[0];
  rb_uplo = argv[1];
  rb_ap = argv[2];
  rb_afp = argv[3];
  rb_equed = argv[4];
  rb_s = argv[5];
  rb_b = argv[6];

  uplo = StringValueCStr(rb_uplo)[0];
  equed = StringValueCStr(rb_equed)[0];
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (7th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (7th argument) must be %d", 2);
  nrhs = NA_SHAPE1(rb_b);
  ldb = NA_SHAPE0(rb_b);
  if (NA_TYPE(rb_b) != NA_DFLOAT)
    rb_b = na_change_type(rb_b, NA_DFLOAT);
  b = NA_PTR_TYPE(rb_b, doublereal*);
  fact = StringValueCStr(rb_fact)[0];
  if (!NA_IsNArray(rb_s))
    rb_raise(rb_eArgError, "s (6th argument) must be NArray");
  if (NA_RANK(rb_s) != 1)
    rb_raise(rb_eArgError, "rank of s (6th argument) must be %d", 1);
  n = NA_SHAPE0(rb_s);
  if (NA_TYPE(rb_s) != NA_DFLOAT)
    rb_s = na_change_type(rb_s, NA_DFLOAT);
  s = NA_PTR_TYPE(rb_s, doublereal*);
  if (!NA_IsNArray(rb_ap))
    rb_raise(rb_eArgError, "ap (3th argument) must be NArray");
  if (NA_RANK(rb_ap) != 1)
    rb_raise(rb_eArgError, "rank of ap (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_ap) != (n*(n+1)/2))
    rb_raise(rb_eRuntimeError, "shape 0 of ap must be %d", n*(n+1)/2);
  if (NA_TYPE(rb_ap) != NA_DFLOAT)
    rb_ap = na_change_type(rb_ap, NA_DFLOAT);
  ap = NA_PTR_TYPE(rb_ap, doublereal*);
  if (!NA_IsNArray(rb_afp))
    rb_raise(rb_eArgError, "afp (4th argument) must be NArray");
  if (NA_RANK(rb_afp) != 1)
    rb_raise(rb_eArgError, "rank of afp (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_afp) != (n*(n+1)/2))
    rb_raise(rb_eRuntimeError, "shape 0 of afp must be %d", n*(n+1)/2);
  if (NA_TYPE(rb_afp) != NA_DFLOAT)
    rb_afp = na_change_type(rb_afp, NA_DFLOAT);
  afp = NA_PTR_TYPE(rb_afp, doublereal*);
  ldx = MAX(1,n);
  {
    int shape[2];
    shape[0] = ldx;
    shape[1] = nrhs;
    rb_x = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  x = NA_PTR_TYPE(rb_x, doublereal*);
  {
    int shape[1];
    shape[0] = nrhs;
    rb_ferr = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  ferr = NA_PTR_TYPE(rb_ferr, doublereal*);
  {
    int shape[1];
    shape[0] = nrhs;
    rb_berr = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  berr = NA_PTR_TYPE(rb_berr, doublereal*);
  {
    int shape[1];
    shape[0] = n*(n+1)/2;
    rb_ap_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  ap_out__ = NA_PTR_TYPE(rb_ap_out__, doublereal*);
  MEMCPY(ap_out__, ap, doublereal, NA_TOTAL(rb_ap));
  rb_ap = rb_ap_out__;
  ap = ap_out__;
  {
    int shape[1];
    shape[0] = n*(n+1)/2;
    rb_afp_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  afp_out__ = NA_PTR_TYPE(rb_afp_out__, doublereal*);
  MEMCPY(afp_out__, afp, doublereal, NA_TOTAL(rb_afp));
  rb_afp = rb_afp_out__;
  afp = afp_out__;
  {
    int shape[1];
    shape[0] = n;
    rb_s_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  s_out__ = NA_PTR_TYPE(rb_s_out__, doublereal*);
  MEMCPY(s_out__, s, doublereal, NA_TOTAL(rb_s));
  rb_s = rb_s_out__;
  s = s_out__;
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = nrhs;
    rb_b_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rb_b_out__, doublereal*);
  MEMCPY(b_out__, b, doublereal, NA_TOTAL(rb_b));
  rb_b = rb_b_out__;
  b = b_out__;
  work = ALLOC_N(doublereal, (3*n));
  iwork = ALLOC_N(integer, (n));

  dppsvx_(&fact, &uplo, &n, &nrhs, ap, afp, &equed, s, b, &ldb, x, &ldx, &rcond, ferr, berr, work, iwork, &info);

  free(work);
  free(iwork);
  rb_rcond = rb_float_new((double)rcond);
  rb_info = INT2NUM(info);
  rb_equed = rb_str_new(&equed,1);
  return rb_ary_new3(10, rb_x, rb_rcond, rb_ferr, rb_berr, rb_info, rb_ap, rb_afp, rb_equed, rb_s, rb_b);
}

void
init_lapack_dppsvx(VALUE mLapack){
  rb_define_module_function(mLapack, "dppsvx", rb_dppsvx, -1);
}
