#include "rb_lapack.h"

static VALUE
rb_dsgesv(int argc, VALUE *argv, VALUE self){
  VALUE rb_a;
  doublereal *a; 
  VALUE rb_b;
  doublereal *b; 
  VALUE rb_ipiv;
  integer *ipiv; 
  VALUE rb_x;
  doublereal *x; 
  VALUE rb_iter;
  integer iter; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  doublereal *a_out__;
  doublereal *work;
  real *swork;

  integer lda;
  integer n;
  integer ldb;
  integer nrhs;
  integer ldx;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  ipiv, x, iter, info, a = NumRu::Lapack.dsgesv( a, b)\n    or\n  NumRu::Lapack.dsgesv  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DSGESV( N, NRHS, A, LDA, IPIV, B, LDB, X, LDX, WORK, SWORK, ITER, INFO )\n\n*  Purpose\n*  =======\n*\n*  DSGESV computes the solution to a real system of linear equations\n*     A * X = B,\n*  where A is an N-by-N matrix and X and B are N-by-NRHS matrices.\n*\n*  DSGESV first attempts to factorize the matrix in SINGLE PRECISION\n*  and use this factorization within an iterative refinement procedure\n*  to produce a solution with DOUBLE PRECISION normwise backward error\n*  quality (see below). If the approach fails the method switches to a\n*  DOUBLE PRECISION factorization and solve.\n*\n*  The iterative refinement is not going to be a winning strategy if\n*  the ratio SINGLE PRECISION performance over DOUBLE PRECISION\n*  performance is too small. A reasonable strategy should take the\n*  number of right-hand sides and the size of the matrix into account.\n*  This might be done with a call to ILAENV in the future. Up to now, we\n*  always try iterative refinement.\n*\n*  The iterative refinement process is stopped if\n*      ITER > ITERMAX\n*  or for all the RHS we have:\n*      RNRM < SQRT(N)*XNRM*ANRM*EPS*BWDMAX\n*  where\n*      o ITER is the number of the current iteration in the iterative\n*        refinement process\n*      o RNRM is the infinity-norm of the residual\n*      o XNRM is the infinity-norm of the solution\n*      o ANRM is the infinity-operator-norm of the matrix A\n*      o EPS is the machine epsilon returned by DLAMCH('Epsilon')\n*  The value ITERMAX and BWDMAX are fixed to 30 and 1.0D+00\n*  respectively.\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The number of linear equations, i.e., the order of the\n*          matrix A.  N >= 0.\n*\n*  NRHS    (input) INTEGER\n*          The number of right hand sides, i.e., the number of columns\n*          of the matrix B.  NRHS >= 0.\n*\n*  A       (input or input/ouptut) DOUBLE PRECISION array,\n*          dimension (LDA,N)\n*          On entry, the N-by-N coefficient matrix A.\n*          On exit, if iterative refinement has been successfully used\n*          (INFO.EQ.0 and ITER.GE.0, see description below), then A is\n*          unchanged, if double precision factorization has been used\n*          (INFO.EQ.0 and ITER.LT.0, see description below), then the\n*          array A contains the factors L and U from the factorization\n*          A = P*L*U; the unit diagonal elements of L are not stored.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,N).\n*\n*  IPIV    (output) INTEGER array, dimension (N)\n*          The pivot indices that define the permutation matrix P;\n*          row i of the matrix was interchanged with row IPIV(i).\n*          Corresponds either to the single precision factorization\n*          (if INFO.EQ.0 and ITER.GE.0) or the double precision\n*          factorization (if INFO.EQ.0 and ITER.LT.0).\n*\n*  B       (input) DOUBLE PRECISION array, dimension (LDB,NRHS)\n*          The N-by-NRHS right hand side matrix B.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B.  LDB >= max(1,N).\n*\n*  X       (output) DOUBLE PRECISION array, dimension (LDX,NRHS)\n*          If INFO = 0, the N-by-NRHS solution matrix X.\n*\n*  LDX     (input) INTEGER\n*          The leading dimension of the array X.  LDX >= max(1,N).\n*\n*  WORK    (workspace) DOUBLE PRECISION array, dimension (N*NRHS)\n*          This array is used to hold the residual vectors.\n*\n*  SWORK   (workspace) REAL array, dimension (N*(N+NRHS))\n*          This array is used to use the single precision matrix and the\n*          right-hand sides or solutions in single precision.\n*\n*  ITER    (output) INTEGER\n*          < 0: iterative refinement has failed, double precision\n*               factorization has been performed\n*               -1 : the routine fell back to full precision for\n*                    implementation- or machine-specific reasons\n*               -2 : narrowing the precision induced an overflow,\n*                    the routine fell back to full precision\n*               -3 : failure of SGETRF\n*               -31: stop the iterative refinement after the 30th\n*                    iterations\n*          > 0: iterative refinement has been sucessfully used.\n*               Returns the number of iterations\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*          > 0:  if INFO = i, U(i,i) computed in DOUBLE PRECISION is\n*                exactly zero.  The factorization has been completed,\n*                but the factor U is exactly singular, so the solution\n*                could not be computed.\n*\n*  =========\n*\n\n");
    return Qnil;
  }
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rb_a = argv[0];
  rb_b = argv[1];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (1th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (1th argument) must be %d", 2);
  lda = NA_SHAPE0(rb_a);
  n = NA_SHAPE1(rb_a);
  if (NA_TYPE(rb_a) != NA_DFLOAT)
    rb_a = na_change_type(rb_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rb_a, doublereal*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (2th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (2th argument) must be %d", 2);
  ldb = NA_SHAPE0(rb_b);
  nrhs = NA_SHAPE1(rb_b);
  if (NA_TYPE(rb_b) != NA_DFLOAT)
    rb_b = na_change_type(rb_b, NA_DFLOAT);
  b = NA_PTR_TYPE(rb_b, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_ipiv = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  ipiv = NA_PTR_TYPE(rb_ipiv, integer*);
  ldx = MAX(1,n);
  {
    int shape[2];
    shape[0] = ldx;
    shape[1] = nrhs;
    rb_x = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  x = NA_PTR_TYPE(rb_x, doublereal*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rb_a_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rb_a_out__, doublereal*);
  MEMCPY(a_out__, a, doublereal, NA_TOTAL(rb_a));
  rb_a = rb_a_out__;
  a = a_out__;
  work = ALLOC_N(doublereal, (n*nrhs));
  swork = ALLOC_N(real, (n*(n+nrhs)));

  dsgesv_(&n, &nrhs, a, &lda, ipiv, b, &ldb, x, &ldx, work, swork, &iter, &info);

  free(work);
  free(swork);
  rb_iter = INT2NUM(iter);
  rb_info = INT2NUM(info);
  return rb_ary_new3(5, rb_ipiv, rb_x, rb_iter, rb_info, rb_a);
}

void
init_lapack_dsgesv(VALUE mLapack){
  rb_define_module_function(mLapack, "dsgesv", rb_dsgesv, -1);
}
