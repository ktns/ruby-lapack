#include "rb_lapack.h"

static VALUE
rb_sggglm(int argc, VALUE *argv, VALUE self){
  VALUE rb_a;
  real *a; 
  VALUE rb_b;
  real *b; 
  VALUE rb_d;
  real *d; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_x;
  real *x; 
  VALUE rb_y;
  real *y; 
  VALUE rb_work;
  real *work; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  real *a_out__;
  VALUE rb_b_out__;
  real *b_out__;
  VALUE rb_d_out__;
  real *d_out__;

  integer lda;
  integer m;
  integer ldb;
  integer p;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  x, y, work, info, a, b, d = NumRu::Lapack.sggglm( a, b, d, lwork)\n    or\n  NumRu::Lapack.sggglm  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SGGGLM( N, M, P, A, LDA, B, LDB, D, X, Y, WORK, LWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  SGGGLM solves a general Gauss-Markov linear model (GLM) problem:\n*\n*          minimize || y ||_2   subject to   d = A*x + B*y\n*              x\n*\n*  where A is an N-by-M matrix, B is an N-by-P matrix, and d is a\n*  given N-vector. It is assumed that M <= N <= M+P, and\n*\n*             rank(A) = M    and    rank( A B ) = N.\n*\n*  Under these assumptions, the constrained equation is always\n*  consistent, and there is a unique solution x and a minimal 2-norm\n*  solution y, which is obtained using a generalized QR factorization\n*  of the matrices (A, B) given by\n*\n*     A = Q*(R),   B = Q*T*Z.\n*           (0)\n*\n*  In particular, if matrix B is square nonsingular, then the problem\n*  GLM is equivalent to the following weighted linear least squares\n*  problem\n*\n*               minimize || inv(B)*(d-A*x) ||_2\n*                   x\n*\n*  where inv(B) denotes the inverse of B.\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The number of rows of the matrices A and B.  N >= 0.\n*\n*  M       (input) INTEGER\n*          The number of columns of the matrix A.  0 <= M <= N.\n*\n*  P       (input) INTEGER\n*          The number of columns of the matrix B.  P >= N-M.\n*\n*  A       (input/output) REAL array, dimension (LDA,M)\n*          On entry, the N-by-M matrix A.\n*          On exit, the upper triangular part of the array A contains\n*          the M-by-M upper triangular matrix R.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A. LDA >= max(1,N).\n*\n*  B       (input/output) REAL array, dimension (LDB,P)\n*          On entry, the N-by-P matrix B.\n*          On exit, if N <= P, the upper triangle of the subarray\n*          B(1:N,P-N+1:P) contains the N-by-N upper triangular matrix T;\n*          if N > P, the elements on and above the (N-P)th subdiagonal\n*          contain the N-by-P upper trapezoidal matrix T.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B. LDB >= max(1,N).\n*\n*  D       (input/output) REAL array, dimension (N)\n*          On entry, D is the left hand side of the GLM equation.\n*          On exit, D is destroyed.\n*\n*  X       (output) REAL array, dimension (M)\n*  Y       (output) REAL array, dimension (P)\n*          On exit, X and Y are the solutions of the GLM problem.\n*\n*  WORK    (workspace/output) REAL array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK. LWORK >= max(1,N+M+P).\n*          For optimum performance, LWORK >= M+min(N,P)+max(N,P)*NB,\n*          where NB is an upper bound for the optimal blocksizes for\n*          SGEQRF, SGERQF, SORMQR and SORMRQ.\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*          = 1:  the upper triangular factor R associated with A in the\n*                generalized QR factorization of the pair (A, B) is\n*                singular, so that rank(A) < M; the least squares\n*                solution could not be computed.\n*          = 2:  the bottom (N-M) by (N-M) part of the upper trapezoidal\n*                factor T associated with B in the generalized QR\n*                factorization of the pair (A, B) is singular, so that\n*                rank( A B ) < N; the least squares solution could not\n*                be computed.\n*\n\n*  ===================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_a = argv[0];
  rb_b = argv[1];
  rb_d = argv[2];
  rb_lwork = argv[3];

  lwork = NUM2INT(rb_lwork);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (1th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (1th argument) must be %d", 2);
  lda = NA_SHAPE0(rb_a);
  m = NA_SHAPE1(rb_a);
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (2th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (2th argument) must be %d", 2);
  ldb = NA_SHAPE0(rb_b);
  p = NA_SHAPE1(rb_b);
  if (NA_TYPE(rb_b) != NA_SFLOAT)
    rb_b = na_change_type(rb_b, NA_SFLOAT);
  b = NA_PTR_TYPE(rb_b, real*);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (3th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (3th argument) must be %d", 1);
  n = NA_SHAPE0(rb_d);
  if (NA_TYPE(rb_d) != NA_SFLOAT)
    rb_d = na_change_type(rb_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rb_d, real*);
  {
    int shape[1];
    shape[0] = m;
    rb_x = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  x = NA_PTR_TYPE(rb_x, real*);
  {
    int shape[1];
    shape[0] = p;
    rb_y = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  y = NA_PTR_TYPE(rb_y, real*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rb_work = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rb_work, real*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = m;
    rb_a_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rb_a_out__, real*);
  MEMCPY(a_out__, a, real, NA_TOTAL(rb_a));
  rb_a = rb_a_out__;
  a = a_out__;
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = p;
    rb_b_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rb_b_out__, real*);
  MEMCPY(b_out__, b, real, NA_TOTAL(rb_b));
  rb_b = rb_b_out__;
  b = b_out__;
  {
    int shape[1];
    shape[0] = n;
    rb_d_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  d_out__ = NA_PTR_TYPE(rb_d_out__, real*);
  MEMCPY(d_out__, d, real, NA_TOTAL(rb_d));
  rb_d = rb_d_out__;
  d = d_out__;

  sggglm_(&n, &m, &p, a, &lda, b, &ldb, d, x, y, work, &lwork, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(7, rb_x, rb_y, rb_work, rb_info, rb_a, rb_b, rb_d);
}

void
init_lapack_sggglm(VALUE mLapack){
  rb_define_module_function(mLapack, "sggglm", rb_sggglm, -1);
}
