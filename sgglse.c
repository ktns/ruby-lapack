#include "rb_lapack.h"

static VALUE
rb_sgglse(int argc, VALUE *argv, VALUE self){
  VALUE rb_a;
  real *a; 
  VALUE rb_b;
  real *b; 
  VALUE rb_c;
  real *c; 
  VALUE rb_d;
  real *d; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_x;
  real *x; 
  VALUE rb_work;
  real *work; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  real *a_out__;
  VALUE rb_b_out__;
  real *b_out__;
  VALUE rb_c_out__;
  real *c_out__;
  VALUE rb_d_out__;
  real *d_out__;

  integer lda;
  integer n;
  integer ldb;
  integer m;
  integer p;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  x, work, info, a, b, c, d = NumRu::Lapack.sgglse( a, b, c, d, lwork)\n    or\n  NumRu::Lapack.sgglse  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SGGLSE( M, N, P, A, LDA, B, LDB, C, D, X, WORK, LWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  SGGLSE solves the linear equality-constrained least squares (LSE)\n*  problem:\n*\n*          minimize || c - A*x ||_2   subject to   B*x = d\n*\n*  where A is an M-by-N matrix, B is a P-by-N matrix, c is a given\n*  M-vector, and d is a given P-vector. It is assumed that\n*  P <= N <= M+P, and\n*\n*           rank(B) = P and  rank( (A) ) = N.\n*                                ( (B) )\n*\n*  These conditions ensure that the LSE problem has a unique solution,\n*  which is obtained using a generalized RQ factorization of the\n*  matrices (B, A) given by\n*\n*     B = (0 R)*Q,   A = Z*T*Q.\n*\n\n*  Arguments\n*  =========\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix A.  M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrices A and B. N >= 0.\n*\n*  P       (input) INTEGER\n*          The number of rows of the matrix B. 0 <= P <= N <= M+P.\n*\n*  A       (input/output) REAL array, dimension (LDA,N)\n*          On entry, the M-by-N matrix A.\n*          On exit, the elements on and above the diagonal of the array\n*          contain the min(M,N)-by-N upper trapezoidal matrix T.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A. LDA >= max(1,M).\n*\n*  B       (input/output) REAL array, dimension (LDB,N)\n*          On entry, the P-by-N matrix B.\n*          On exit, the upper triangle of the subarray B(1:P,N-P+1:N)\n*          contains the P-by-P upper triangular matrix R.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B. LDB >= max(1,P).\n*\n*  C       (input/output) REAL array, dimension (M)\n*          On entry, C contains the right hand side vector for the\n*          least squares part of the LSE problem.\n*          On exit, the residual sum of squares for the solution\n*          is given by the sum of squares of elements N-P+1 to M of\n*          vector C.\n*\n*  D       (input/output) REAL array, dimension (P)\n*          On entry, D contains the right hand side vector for the\n*          constrained equation.\n*          On exit, D is destroyed.\n*\n*  X       (output) REAL array, dimension (N)\n*          On exit, X is the solution of the LSE problem.\n*\n*  WORK    (workspace/output) REAL array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK. LWORK >= max(1,M+N+P).\n*          For optimum performance LWORK >= P+min(M,N)+max(M,N)*NB,\n*          where NB is an upper bound for the optimal blocksizes for\n*          SGEQRF, SGERQF, SORMQR and SORMRQ.\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*          = 1:  the upper triangular factor R associated with B in the\n*                generalized RQ factorization of the pair (B, A) is\n*                singular, so that rank(B) < P; the least squares\n*                solution could not be computed.\n*          = 2:  the (N-P) by (N-P) part of the upper trapezoidal factor\n*                T associated with A in the generalized RQ factorization\n*                of the pair (B, A) is singular, so that\n*                rank( (A) ) < N; the least squares solution could not\n*                    ( (B) )\n*                be computed.\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_a = argv[0];
  rb_b = argv[1];
  rb_c = argv[2];
  rb_d = argv[3];
  rb_lwork = argv[4];

  lwork = NUM2INT(rb_lwork);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (1th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (1th argument) must be %d", 2);
  lda = NA_SHAPE0(rb_a);
  n = NA_SHAPE1(rb_a);
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (2th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (2th argument) must be %d", 2);
  ldb = NA_SHAPE0(rb_b);
  if (NA_SHAPE1(rb_b) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of b must be the same as shape 1 of a");
  if (NA_TYPE(rb_b) != NA_SFLOAT)
    rb_b = na_change_type(rb_b, NA_SFLOAT);
  b = NA_PTR_TYPE(rb_b, real*);
  if (!NA_IsNArray(rb_c))
    rb_raise(rb_eArgError, "c (3th argument) must be NArray");
  if (NA_RANK(rb_c) != 1)
    rb_raise(rb_eArgError, "rank of c (3th argument) must be %d", 1);
  m = NA_SHAPE0(rb_c);
  if (NA_TYPE(rb_c) != NA_SFLOAT)
    rb_c = na_change_type(rb_c, NA_SFLOAT);
  c = NA_PTR_TYPE(rb_c, real*);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (4th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (4th argument) must be %d", 1);
  p = NA_SHAPE0(rb_d);
  if (NA_TYPE(rb_d) != NA_SFLOAT)
    rb_d = na_change_type(rb_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rb_d, real*);
  {
    int shape[1];
    shape[0] = n;
    rb_x = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  x = NA_PTR_TYPE(rb_x, real*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rb_work = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rb_work, real*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rb_a_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rb_a_out__, real*);
  MEMCPY(a_out__, a, real, NA_TOTAL(rb_a));
  rb_a = rb_a_out__;
  a = a_out__;
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = n;
    rb_b_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rb_b_out__, real*);
  MEMCPY(b_out__, b, real, NA_TOTAL(rb_b));
  rb_b = rb_b_out__;
  b = b_out__;
  {
    int shape[1];
    shape[0] = m;
    rb_c_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  c_out__ = NA_PTR_TYPE(rb_c_out__, real*);
  MEMCPY(c_out__, c, real, NA_TOTAL(rb_c));
  rb_c = rb_c_out__;
  c = c_out__;
  {
    int shape[1];
    shape[0] = p;
    rb_d_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  d_out__ = NA_PTR_TYPE(rb_d_out__, real*);
  MEMCPY(d_out__, d, real, NA_TOTAL(rb_d));
  rb_d = rb_d_out__;
  d = d_out__;

  sgglse_(&m, &n, &p, a, &lda, b, &ldb, c, d, x, work, &lwork, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(7, rb_x, rb_work, rb_info, rb_a, rb_b, rb_c, rb_d);
}

void
init_lapack_sgglse(VALUE mLapack){
  rb_define_module_function(mLapack, "sgglse", rb_sgglse, -1);
}
