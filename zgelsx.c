#include "rb_lapack.h"

extern VOID zgelsx_(integer *m, integer *n, integer *nrhs, doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, integer *jpvt, doublereal *rcond, integer *rank, doublecomplex *work, doublereal *rwork, integer *info);

static VALUE
rb_zgelsx(int argc, VALUE *argv, VALUE self){
  VALUE rb_m;
  integer m; 
  VALUE rb_a;
  doublecomplex *a; 
  VALUE rb_b;
  doublecomplex *b; 
  VALUE rb_jpvt;
  integer *jpvt; 
  VALUE rb_rcond;
  doublereal rcond; 
  VALUE rb_rank;
  integer rank; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  doublecomplex *a_out__;
  VALUE rb_b_out__;
  doublecomplex *b_out__;
  VALUE rb_jpvt_out__;
  integer *jpvt_out__;
  doublecomplex *work;
  doublereal *rwork;

  integer lda;
  integer n;
  integer ldb;
  integer nrhs;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  rank, info, a, b, jpvt = NumRu::Lapack.zgelsx( m, a, b, jpvt, rcond)\n    or\n  NumRu::Lapack.zgelsx  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZGELSX( M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND, RANK, WORK, RWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  This routine is deprecated and has been replaced by routine ZGELSY.\n*\n*  ZGELSX computes the minimum-norm solution to a complex linear least\n*  squares problem:\n*      minimize || A * X - B ||\n*  using a complete orthogonal factorization of A.  A is an M-by-N\n*  matrix which may be rank-deficient.\n*\n*  Several right hand side vectors b and solution vectors x can be\n*  handled in a single call; they are stored as the columns of the\n*  M-by-NRHS right hand side matrix B and the N-by-NRHS solution\n*  matrix X.\n*\n*  The routine first computes a QR factorization with column pivoting:\n*      A * P = Q * [ R11 R12 ]\n*                  [  0  R22 ]\n*  with R11 defined as the largest leading submatrix whose estimated\n*  condition number is less than 1/RCOND.  The order of R11, RANK,\n*  is the effective rank of A.\n*\n*  Then, R22 is considered to be negligible, and R12 is annihilated\n*  by unitary transformations from the right, arriving at the\n*  complete orthogonal factorization:\n*     A * P = Q * [ T11 0 ] * Z\n*                 [  0  0 ]\n*  The minimum-norm solution is then\n*     X = P * Z' [ inv(T11)*Q1'*B ]\n*                [        0       ]\n*  where Q1 consists of the first RANK columns of Q.\n*\n\n*  Arguments\n*  =========\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix A.  M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix A.  N >= 0.\n*\n*  NRHS    (input) INTEGER\n*          The number of right hand sides, i.e., the number of\n*          columns of matrices B and X. NRHS >= 0.\n*\n*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)\n*          On entry, the M-by-N matrix A.\n*          On exit, A has been overwritten by details of its\n*          complete orthogonal factorization.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,M).\n*\n*  B       (input/output) COMPLEX*16 array, dimension (LDB,NRHS)\n*          On entry, the M-by-NRHS right hand side matrix B.\n*          On exit, the N-by-NRHS solution matrix X.\n*          If m >= n and RANK = n, the residual sum-of-squares for\n*          the solution in the i-th column is given by the sum of\n*          squares of elements N+1:M in that column.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B. LDB >= max(1,M,N).\n*\n*  JPVT    (input/output) INTEGER array, dimension (N)\n*          On entry, if JPVT(i) .ne. 0, the i-th column of A is an\n*          initial column, otherwise it is a free column.  Before\n*          the QR factorization of A, all initial columns are\n*          permuted to the leading positions; only the remaining\n*          free columns are moved as a result of column pivoting\n*          during the factorization.\n*          On exit, if JPVT(i) = k, then the i-th column of A*P\n*          was the k-th column of A.\n*\n*  RCOND   (input) DOUBLE PRECISION\n*          RCOND is used to determine the effective rank of A, which\n*          is defined as the order of the largest leading triangular\n*          submatrix R11 in the QR factorization with pivoting of A,\n*          whose estimated condition number < 1/RCOND.\n*\n*  RANK    (output) INTEGER\n*          The effective rank of A, i.e., the order of the submatrix\n*          R11.  This is the same as the order of the submatrix T11\n*          in the complete orthogonal factorization of A.\n*\n*  WORK    (workspace) COMPLEX*16 array, dimension\n*                      (min(M,N) + max( N, 2*min(M,N)+NRHS )),\n*\n*  RWORK   (workspace) DOUBLE PRECISION array, dimension (2*N)\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_m = argv[0];
  rb_a = argv[1];
  rb_b = argv[2];
  rb_jpvt = argv[3];
  rb_rcond = argv[4];

  rcond = NUM2DBL(rb_rcond);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_DCOMPLEX)
    rb_a = na_change_type(rb_a, NA_DCOMPLEX);
  a = NA_PTR_TYPE(rb_a, doublecomplex*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (3th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (3th argument) must be %d", 2);
  nrhs = NA_SHAPE1(rb_b);
  ldb = NA_SHAPE0(rb_b);
  if (NA_TYPE(rb_b) != NA_DCOMPLEX)
    rb_b = na_change_type(rb_b, NA_DCOMPLEX);
  b = NA_PTR_TYPE(rb_b, doublecomplex*);
  m = NUM2INT(rb_m);
  if (!NA_IsNArray(rb_jpvt))
    rb_raise(rb_eArgError, "jpvt (4th argument) must be NArray");
  if (NA_RANK(rb_jpvt) != 1)
    rb_raise(rb_eArgError, "rank of jpvt (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_jpvt) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of jpvt must be the same as shape 1 of a");
  if (NA_TYPE(rb_jpvt) != NA_LINT)
    rb_jpvt = na_change_type(rb_jpvt, NA_LINT);
  jpvt = NA_PTR_TYPE(rb_jpvt, integer*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rb_a_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rb_a_out__, doublecomplex*);
  MEMCPY(a_out__, a, doublecomplex, NA_TOTAL(rb_a));
  rb_a = rb_a_out__;
  a = a_out__;
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = nrhs;
    rb_b_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rb_b_out__, doublecomplex*);
  MEMCPY(b_out__, b, doublecomplex, NA_TOTAL(rb_b));
  rb_b = rb_b_out__;
  b = b_out__;
  {
    int shape[1];
    shape[0] = n;
    rb_jpvt_out__ = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  jpvt_out__ = NA_PTR_TYPE(rb_jpvt_out__, integer*);
  MEMCPY(jpvt_out__, jpvt, integer, NA_TOTAL(rb_jpvt));
  rb_jpvt = rb_jpvt_out__;
  jpvt = jpvt_out__;
  work = ALLOC_N(doublecomplex, (MIN(m,n) + MAX(n,2*(MIN(m,n))+nrhs)));
  rwork = ALLOC_N(doublereal, (2*n));

  zgelsx_(&m, &n, &nrhs, a, &lda, b, &ldb, jpvt, &rcond, &rank, work, rwork, &info);

  free(work);
  free(rwork);
  rb_rank = INT2NUM(rank);
  rb_info = INT2NUM(info);
  return rb_ary_new3(5, rb_rank, rb_info, rb_a, rb_b, rb_jpvt);
}

void
init_lapack_zgelsx(VALUE mLapack){
  rb_define_module_function(mLapack, "zgelsx", rb_zgelsx, -1);
}
