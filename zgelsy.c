#include "rb_lapack.h"

static VALUE
rb_zgelsy(int argc, VALUE *argv, VALUE self){
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
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_rank;
  integer rank; 
  VALUE rb_work;
  doublecomplex *work; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  doublecomplex *a_out__;
  VALUE rb_b_out__;
  doublecomplex *b_out__;
  VALUE rb_jpvt_out__;
  integer *jpvt_out__;
  doublereal *rwork;

  integer lda;
  integer n;
  integer ldb;
  integer nrhs;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  rank, work, info, a, b, jpvt = NumRu::Lapack.zgelsy( m, a, b, jpvt, rcond, lwork)\n    or\n  NumRu::Lapack.zgelsy  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZGELSY( M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND, RANK, WORK, LWORK, RWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  ZGELSY computes the minimum-norm solution to a complex linear least\n*  squares problem:\n*      minimize || A * X - B ||\n*  using a complete orthogonal factorization of A.  A is an M-by-N\n*  matrix which may be rank-deficient.\n*\n*  Several right hand side vectors b and solution vectors x can be\n*  handled in a single call; they are stored as the columns of the\n*  M-by-NRHS right hand side matrix B and the N-by-NRHS solution\n*  matrix X.\n*\n*  The routine first computes a QR factorization with column pivoting:\n*      A * P = Q * [ R11 R12 ]\n*                  [  0  R22 ]\n*  with R11 defined as the largest leading submatrix whose estimated\n*  condition number is less than 1/RCOND.  The order of R11, RANK,\n*  is the effective rank of A.\n*\n*  Then, R22 is considered to be negligible, and R12 is annihilated\n*  by unitary transformations from the right, arriving at the\n*  complete orthogonal factorization:\n*     A * P = Q * [ T11 0 ] * Z\n*                 [  0  0 ]\n*  The minimum-norm solution is then\n*     X = P * Z' [ inv(T11)*Q1'*B ]\n*                [        0       ]\n*  where Q1 consists of the first RANK columns of Q.\n*\n*  This routine is basically identical to the original xGELSX except\n*  three differences:\n*    o The permutation of matrix B (the right hand side) is faster and\n*      more simple.\n*    o The call to the subroutine xGEQPF has been substituted by the\n*      the call to the subroutine xGEQP3. This subroutine is a Blas-3\n*      version of the QR factorization with column pivoting.\n*    o Matrix B (the right hand side) is updated with Blas-3.\n*\n\n*  Arguments\n*  =========\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix A.  M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix A.  N >= 0.\n*\n*  NRHS    (input) INTEGER\n*          The number of right hand sides, i.e., the number of\n*          columns of matrices B and X. NRHS >= 0.\n*\n*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)\n*          On entry, the M-by-N matrix A.\n*          On exit, A has been overwritten by details of its\n*          complete orthogonal factorization.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,M).\n*\n*  B       (input/output) COMPLEX*16 array, dimension (LDB,NRHS)\n*          On entry, the M-by-NRHS right hand side matrix B.\n*          On exit, the N-by-NRHS solution matrix X.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B. LDB >= max(1,M,N).\n*\n*  JPVT    (input/output) INTEGER array, dimension (N)\n*          On entry, if JPVT(i) .ne. 0, the i-th column of A is permuted\n*          to the front of AP, otherwise column i is a free column.\n*          On exit, if JPVT(i) = k, then the i-th column of A*P\n*          was the k-th column of A.\n*\n*  RCOND   (input) DOUBLE PRECISION\n*          RCOND is used to determine the effective rank of A, which\n*          is defined as the order of the largest leading triangular\n*          submatrix R11 in the QR factorization with pivoting of A,\n*          whose estimated condition number < 1/RCOND.\n*\n*  RANK    (output) INTEGER\n*          The effective rank of A, i.e., the order of the submatrix\n*          R11.  This is the same as the order of the submatrix T11\n*          in the complete orthogonal factorization of A.\n*\n*  WORK    (workspace/output) COMPLEX*16 array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK.\n*          The unblocked strategy requires that:\n*            LWORK >= MN + MAX( 2*MN, N+1, MN+NRHS )\n*          where MN = min(M,N).\n*          The block algorithm requires that:\n*            LWORK >= MN + MAX( 2*MN, NB*(N+1), MN+MN*NB, MN+NB*NRHS )\n*          where NB is an upper bound on the blocksize returned\n*          by ILAENV for the routines ZGEQP3, ZTZRZF, CTZRQF, ZUNMQR,\n*          and ZUNMRZ.\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  RWORK   (workspace) DOUBLE PRECISION array, dimension (2*N)\n*\n*  INFO    (output) INTEGER\n*          = 0: successful exit\n*          < 0: if INFO = -i, the i-th argument had an illegal value\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*    A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville, USA\n*    E. Quintana-Orti, Depto. de Informatica, Universidad Jaime I, Spain\n*    G. Quintana-Orti, Depto. de Informatica, Universidad Jaime I, Spain\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rb_m = argv[0];
  rb_a = argv[1];
  rb_b = argv[2];
  rb_jpvt = argv[3];
  rb_rcond = argv[4];
  rb_lwork = argv[5];

  m = NUM2INT(rb_m);
  rcond = NUM2DBL(rb_rcond);
  lwork = NUM2INT(rb_lwork);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  lda = NA_SHAPE0(rb_a);
  n = NA_SHAPE1(rb_a);
  if (NA_TYPE(rb_a) != NA_DCOMPLEX)
    rb_a = na_change_type(rb_a, NA_DCOMPLEX);
  a = NA_PTR_TYPE(rb_a, doublecomplex*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (3th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (3th argument) must be %d", 2);
  ldb = NA_SHAPE0(rb_b);
  nrhs = NA_SHAPE1(rb_b);
  if (NA_TYPE(rb_b) != NA_DCOMPLEX)
    rb_b = na_change_type(rb_b, NA_DCOMPLEX);
  b = NA_PTR_TYPE(rb_b, doublecomplex*);
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
    int shape[1];
    shape[0] = MAX(1,lwork);
    rb_work = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rb_work, doublecomplex*);
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
    shape[0] = DIM_LEN(n);
    rb_jpvt_out__ = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  jpvt_out__ = NA_PTR_TYPE(rb_jpvt_out__, integer*);
  MEMCPY(jpvt_out__, jpvt, integer, NA_TOTAL(rb_jpvt));
  rb_jpvt = rb_jpvt_out__;
  jpvt = jpvt_out__;
  rwork = ALLOC_N(doublereal, (2*n));

  zgelsy_(&m, &n, &nrhs, a, &lda, b, &ldb, jpvt, &rcond, &rank, work, &lwork, rwork, &info);

  free(rwork);
  rb_rank = INT2NUM(rank);
  rb_info = INT2NUM(info);
  return rb_ary_new3(6, rb_rank, rb_work, rb_info, rb_a, rb_b, rb_jpvt);
}

void
init_lapack_zgelsy(VALUE mLapack){
  rb_define_module_function(mLapack, "zgelsy", rb_zgelsy, -1);
}
