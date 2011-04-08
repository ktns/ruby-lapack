#include "rb_lapack.h"

extern VOID dgelsx_(integer *m, integer *n, integer *nrhs, doublereal *a, integer *lda, doublereal *b, integer *ldb, integer *jpvt, doublereal *rcond, integer *rank, doublereal *work, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dgelsx(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_m;
  integer m; 
  VALUE rblapack_a;
  doublereal *a; 
  VALUE rblapack_b;
  doublereal *b; 
  VALUE rblapack_jpvt;
  integer *jpvt; 
  VALUE rblapack_rcond;
  doublereal rcond; 
  VALUE rblapack_rank;
  integer rank; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_a_out__;
  doublereal *a_out__;
  VALUE rblapack_b_out__;
  doublereal *b_out__;
  VALUE rblapack_jpvt_out__;
  integer *jpvt_out__;
  doublereal *work;

  integer lda;
  integer n;
  integer ldb;
  integer nrhs;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  rank, info, a, b, jpvt = NumRu::Lapack.dgelsx( m, a, b, jpvt, rcond, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DGELSX( M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND, RANK, WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  This routine is deprecated and has been replaced by routine DGELSY.\n*\n*  DGELSX computes the minimum-norm solution to a real linear least\n*  squares problem:\n*      minimize || A * X - B ||\n*  using a complete orthogonal factorization of A.  A is an M-by-N\n*  matrix which may be rank-deficient.\n*\n*  Several right hand side vectors b and solution vectors x can be\n*  handled in a single call; they are stored as the columns of the\n*  M-by-NRHS right hand side matrix B and the N-by-NRHS solution\n*  matrix X.\n*\n*  The routine first computes a QR factorization with column pivoting:\n*      A * P = Q * [ R11 R12 ]\n*                  [  0  R22 ]\n*  with R11 defined as the largest leading submatrix whose estimated\n*  condition number is less than 1/RCOND.  The order of R11, RANK,\n*  is the effective rank of A.\n*\n*  Then, R22 is considered to be negligible, and R12 is annihilated\n*  by orthogonal transformations from the right, arriving at the\n*  complete orthogonal factorization:\n*     A * P = Q * [ T11 0 ] * Z\n*                 [  0  0 ]\n*  The minimum-norm solution is then\n*     X = P * Z' [ inv(T11)*Q1'*B ]\n*                [        0       ]\n*  where Q1 consists of the first RANK columns of Q.\n*\n\n*  Arguments\n*  =========\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix A.  M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix A.  N >= 0.\n*\n*  NRHS    (input) INTEGER\n*          The number of right hand sides, i.e., the number of\n*          columns of matrices B and X. NRHS >= 0.\n*\n*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)\n*          On entry, the M-by-N matrix A.\n*          On exit, A has been overwritten by details of its\n*          complete orthogonal factorization.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,M).\n*\n*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)\n*          On entry, the M-by-NRHS right hand side matrix B.\n*          On exit, the N-by-NRHS solution matrix X.\n*          If m >= n and RANK = n, the residual sum-of-squares for\n*          the solution in the i-th column is given by the sum of\n*          squares of elements N+1:M in that column.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B. LDB >= max(1,M,N).\n*\n*  JPVT    (input/output) INTEGER array, dimension (N)\n*          On entry, if JPVT(i) .ne. 0, the i-th column of A is an\n*          initial column, otherwise it is a free column.  Before\n*          the QR factorization of A, all initial columns are\n*          permuted to the leading positions; only the remaining\n*          free columns are moved as a result of column pivoting\n*          during the factorization.\n*          On exit, if JPVT(i) = k, then the i-th column of A*P\n*          was the k-th column of A.\n*\n*  RCOND   (input) DOUBLE PRECISION\n*          RCOND is used to determine the effective rank of A, which\n*          is defined as the order of the largest leading triangular\n*          submatrix R11 in the QR factorization with pivoting of A,\n*          whose estimated condition number < 1/RCOND.\n*\n*  RANK    (output) INTEGER\n*          The effective rank of A, i.e., the order of the submatrix\n*          R11.  This is the same as the order of the submatrix T11\n*          in the complete orthogonal factorization of A.\n*\n*  WORK    (workspace) DOUBLE PRECISION array, dimension\n*                      (max( min(M,N)+3*N, 2*min(M,N)+NRHS )),\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  rank, info, a, b, jpvt = NumRu::Lapack.dgelsx( m, a, b, jpvt, rcond, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rblapack_m = argv[0];
  rblapack_a = argv[1];
  rblapack_b = argv[2];
  rblapack_jpvt = argv[3];
  rblapack_rcond = argv[4];
  if (rb_options != Qnil) {
  }

  rcond = NUM2DBL(rblapack_rcond);
  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_a);
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_DFLOAT)
    rblapack_a = na_change_type(rblapack_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rblapack_a, doublereal*);
  if (!NA_IsNArray(rblapack_b))
    rb_raise(rb_eArgError, "b (3th argument) must be NArray");
  if (NA_RANK(rblapack_b) != 2)
    rb_raise(rb_eArgError, "rank of b (3th argument) must be %d", 2);
  nrhs = NA_SHAPE1(rblapack_b);
  ldb = NA_SHAPE0(rblapack_b);
  if (NA_TYPE(rblapack_b) != NA_DFLOAT)
    rblapack_b = na_change_type(rblapack_b, NA_DFLOAT);
  b = NA_PTR_TYPE(rblapack_b, doublereal*);
  m = NUM2INT(rblapack_m);
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
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rblapack_a_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rblapack_a_out__, doublereal*);
  MEMCPY(a_out__, a, doublereal, NA_TOTAL(rblapack_a));
  rblapack_a = rblapack_a_out__;
  a = a_out__;
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = nrhs;
    rblapack_b_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rblapack_b_out__, doublereal*);
  MEMCPY(b_out__, b, doublereal, NA_TOTAL(rblapack_b));
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
  work = ALLOC_N(doublereal, (MAX((MIN(m,n))+3*n,2*(MIN(m,n))*nrhs)));

  dgelsx_(&m, &n, &nrhs, a, &lda, b, &ldb, jpvt, &rcond, &rank, work, &info);

  free(work);
  rblapack_rank = INT2NUM(rank);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(5, rblapack_rank, rblapack_info, rblapack_a, rblapack_b, rblapack_jpvt);
}

void
init_lapack_dgelsx(VALUE mLapack){
  rb_define_module_function(mLapack, "dgelsx", rblapack_dgelsx, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
