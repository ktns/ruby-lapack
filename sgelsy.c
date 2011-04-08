#include "rb_lapack.h"

extern VOID sgelsy_(integer *m, integer *n, integer *nrhs, real *a, integer *lda, real *b, integer *ldb, integer *jpvt, real *rcond, integer *rank, real *work, integer *lwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_sgelsy(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_m;
  integer m; 
  VALUE rblapack_a;
  real *a; 
  VALUE rblapack_b;
  real *b; 
  VALUE rblapack_jpvt;
  integer *jpvt; 
  VALUE rblapack_rcond;
  real rcond; 
  VALUE rblapack_lwork;
  integer lwork; 
  VALUE rblapack_rank;
  integer rank; 
  VALUE rblapack_work;
  real *work; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_a_out__;
  real *a_out__;
  VALUE rblapack_b_out__;
  real *b_out__;
  VALUE rblapack_jpvt_out__;
  integer *jpvt_out__;

  integer lda;
  integer n;
  integer ldb;
  integer nrhs;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  rank, work, info, a, b, jpvt = NumRu::Lapack.sgelsy( m, a, b, jpvt, rcond, lwork, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE SGELSY( M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND, RANK, WORK, LWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  SGELSY computes the minimum-norm solution to a real linear least\n*  squares problem:\n*      minimize || A * X - B ||\n*  using a complete orthogonal factorization of A.  A is an M-by-N\n*  matrix which may be rank-deficient.\n*\n*  Several right hand side vectors b and solution vectors x can be\n*  handled in a single call; they are stored as the columns of the\n*  M-by-NRHS right hand side matrix B and the N-by-NRHS solution\n*  matrix X.\n*\n*  The routine first computes a QR factorization with column pivoting:\n*      A * P = Q * [ R11 R12 ]\n*                  [  0  R22 ]\n*  with R11 defined as the largest leading submatrix whose estimated\n*  condition number is less than 1/RCOND.  The order of R11, RANK,\n*  is the effective rank of A.\n*\n*  Then, R22 is considered to be negligible, and R12 is annihilated\n*  by orthogonal transformations from the right, arriving at the\n*  complete orthogonal factorization:\n*     A * P = Q * [ T11 0 ] * Z\n*                 [  0  0 ]\n*  The minimum-norm solution is then\n*     X = P * Z' [ inv(T11)*Q1'*B ]\n*                [        0       ]\n*  where Q1 consists of the first RANK columns of Q.\n*\n*  This routine is basically identical to the original xGELSX except\n*  three differences:\n*    o The call to the subroutine xGEQPF has been substituted by the\n*      the call to the subroutine xGEQP3. This subroutine is a Blas-3\n*      version of the QR factorization with column pivoting.\n*    o Matrix B (the right hand side) is updated with Blas-3.\n*    o The permutation of matrix B (the right hand side) is faster and\n*      more simple.\n*\n\n*  Arguments\n*  =========\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix A.  M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix A.  N >= 0.\n*\n*  NRHS    (input) INTEGER\n*          The number of right hand sides, i.e., the number of\n*          columns of matrices B and X. NRHS >= 0.\n*\n*  A       (input/output) REAL array, dimension (LDA,N)\n*          On entry, the M-by-N matrix A.\n*          On exit, A has been overwritten by details of its\n*          complete orthogonal factorization.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,M).\n*\n*  B       (input/output) REAL array, dimension (LDB,NRHS)\n*          On entry, the M-by-NRHS right hand side matrix B.\n*          On exit, the N-by-NRHS solution matrix X.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B. LDB >= max(1,M,N).\n*\n*  JPVT    (input/output) INTEGER array, dimension (N)\n*          On entry, if JPVT(i) .ne. 0, the i-th column of A is permuted\n*          to the front of AP, otherwise column i is a free column.\n*          On exit, if JPVT(i) = k, then the i-th column of AP\n*          was the k-th column of A.\n*\n*  RCOND   (input) REAL\n*          RCOND is used to determine the effective rank of A, which\n*          is defined as the order of the largest leading triangular\n*          submatrix R11 in the QR factorization with pivoting of A,\n*          whose estimated condition number < 1/RCOND.\n*\n*  RANK    (output) INTEGER\n*          The effective rank of A, i.e., the order of the submatrix\n*          R11.  This is the same as the order of the submatrix T11\n*          in the complete orthogonal factorization of A.\n*\n*  WORK    (workspace/output) REAL array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK.\n*          The unblocked strategy requires that:\n*             LWORK >= MAX( MN+3*N+1, 2*MN+NRHS ),\n*          where MN = min( M, N ).\n*          The block algorithm requires that:\n*             LWORK >= MAX( MN+2*N+NB*(N+1), 2*MN+NB*NRHS ),\n*          where NB is an upper bound on the blocksize returned\n*          by ILAENV for the routines SGEQP3, STZRZF, STZRQF, SORMQR,\n*          and SORMRZ.\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  INFO    (output) INTEGER\n*          = 0: successful exit\n*          < 0: If INFO = -i, the i-th argument had an illegal value.\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*    A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville, USA\n*    E. Quintana-Orti, Depto. de Informatica, Universidad Jaime I, Spain\n*    G. Quintana-Orti, Depto. de Informatica, Universidad Jaime I, Spain\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  rank, work, info, a, b, jpvt = NumRu::Lapack.sgelsy( m, a, b, jpvt, rcond, lwork, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rblapack_m = argv[0];
  rblapack_a = argv[1];
  rblapack_b = argv[2];
  rblapack_jpvt = argv[3];
  rblapack_rcond = argv[4];
  rblapack_lwork = argv[5];
  if (rb_options != Qnil) {
  }

  rcond = (real)NUM2DBL(rblapack_rcond);
  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_a);
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_SFLOAT)
    rblapack_a = na_change_type(rblapack_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rblapack_a, real*);
  if (!NA_IsNArray(rblapack_b))
    rb_raise(rb_eArgError, "b (3th argument) must be NArray");
  if (NA_RANK(rblapack_b) != 2)
    rb_raise(rb_eArgError, "rank of b (3th argument) must be %d", 2);
  nrhs = NA_SHAPE1(rblapack_b);
  ldb = NA_SHAPE0(rblapack_b);
  if (NA_TYPE(rblapack_b) != NA_SFLOAT)
    rblapack_b = na_change_type(rblapack_b, NA_SFLOAT);
  b = NA_PTR_TYPE(rblapack_b, real*);
  m = NUM2INT(rblapack_m);
  lwork = NUM2INT(rblapack_lwork);
  if (!NA_IsNArray(rblapack_jpvt))
    rb_raise(rb_eArgError, "jpvt (4th argument) must be NArray");
  if (NA_RANK(rblapack_jpvt) != 1)
    rb_raise(rb_eArgError, "rank of jpvt (4th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_jpvt) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of jpvt must be the same as shape 1 of a");
  if (NA_TYPE(rblapack_jpvt) != NA_LINT)
    rblapack_jpvt = na_change_type(rblapack_jpvt, NA_LINT);
  jpvt = NA_PTR_TYPE(rblapack_jpvt, integer*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rblapack_work = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rblapack_work, real*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rblapack_a_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rblapack_a_out__, real*);
  MEMCPY(a_out__, a, real, NA_TOTAL(rblapack_a));
  rblapack_a = rblapack_a_out__;
  a = a_out__;
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = nrhs;
    rblapack_b_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rblapack_b_out__, real*);
  MEMCPY(b_out__, b, real, NA_TOTAL(rblapack_b));
  rblapack_b = rblapack_b_out__;
  b = b_out__;
  {
    int shape[1];
    shape[0] = n;
    rblapack_jpvt_out__ = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  jpvt_out__ = NA_PTR_TYPE(rblapack_jpvt_out__, integer*);
  MEMCPY(jpvt_out__, jpvt, integer, NA_TOTAL(rblapack_jpvt));
  rblapack_jpvt = rblapack_jpvt_out__;
  jpvt = jpvt_out__;

  sgelsy_(&m, &n, &nrhs, a, &lda, b, &ldb, jpvt, &rcond, &rank, work, &lwork, &info);

  rblapack_rank = INT2NUM(rank);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(6, rblapack_rank, rblapack_work, rblapack_info, rblapack_a, rblapack_b, rblapack_jpvt);
}

void
init_lapack_sgelsy(VALUE mLapack){
  rb_define_module_function(mLapack, "sgelsy", rblapack_sgelsy, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
