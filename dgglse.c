#include "rb_lapack.h"

extern VOID dgglse_(integer *m, integer *n, integer *p, doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *c, doublereal *d, doublereal *x, doublereal *work, integer *lwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dgglse(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_a;
  doublereal *a; 
  VALUE rblapack_b;
  doublereal *b; 
  VALUE rblapack_c;
  doublereal *c; 
  VALUE rblapack_d;
  doublereal *d; 
  VALUE rblapack_lwork;
  integer lwork; 
  VALUE rblapack_x;
  doublereal *x; 
  VALUE rblapack_work;
  doublereal *work; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_a_out__;
  doublereal *a_out__;
  VALUE rblapack_b_out__;
  doublereal *b_out__;
  VALUE rblapack_c_out__;
  doublereal *c_out__;
  VALUE rblapack_d_out__;
  doublereal *d_out__;

  integer lda;
  integer n;
  integer ldb;
  integer m;
  integer p;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  x, work, info, a, b, c, d = NumRu::Lapack.dgglse( a, b, c, d, lwork, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DGGLSE( M, N, P, A, LDA, B, LDB, C, D, X, WORK, LWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  DGGLSE solves the linear equality-constrained least squares (LSE)\n*  problem:\n*\n*          minimize || c - A*x ||_2   subject to   B*x = d\n*\n*  where A is an M-by-N matrix, B is a P-by-N matrix, c is a given\n*  M-vector, and d is a given P-vector. It is assumed that\n*  P <= N <= M+P, and\n*\n*           rank(B) = P and  rank( (A) ) = N.\n*                                ( (B) )\n*\n*  These conditions ensure that the LSE problem has a unique solution,\n*  which is obtained using a generalized RQ factorization of the\n*  matrices (B, A) given by\n*\n*     B = (0 R)*Q,   A = Z*T*Q.\n*\n\n*  Arguments\n*  =========\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix A.  M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrices A and B. N >= 0.\n*\n*  P       (input) INTEGER\n*          The number of rows of the matrix B. 0 <= P <= N <= M+P.\n*\n*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)\n*          On entry, the M-by-N matrix A.\n*          On exit, the elements on and above the diagonal of the array\n*          contain the min(M,N)-by-N upper trapezoidal matrix T.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A. LDA >= max(1,M).\n*\n*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,N)\n*          On entry, the P-by-N matrix B.\n*          On exit, the upper triangle of the subarray B(1:P,N-P+1:N)\n*          contains the P-by-P upper triangular matrix R.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B. LDB >= max(1,P).\n*\n*  C       (input/output) DOUBLE PRECISION array, dimension (M)\n*          On entry, C contains the right hand side vector for the\n*          least squares part of the LSE problem.\n*          On exit, the residual sum of squares for the solution\n*          is given by the sum of squares of elements N-P+1 to M of\n*          vector C.\n*\n*  D       (input/output) DOUBLE PRECISION array, dimension (P)\n*          On entry, D contains the right hand side vector for the\n*          constrained equation.\n*          On exit, D is destroyed.\n*\n*  X       (output) DOUBLE PRECISION array, dimension (N)\n*          On exit, X is the solution of the LSE problem.\n*\n*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK. LWORK >= max(1,M+N+P).\n*          For optimum performance LWORK >= P+min(M,N)+max(M,N)*NB,\n*          where NB is an upper bound for the optimal blocksizes for\n*          DGEQRF, SGERQF, DORMQR and SORMRQ.\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*          = 1:  the upper triangular factor R associated with B in the\n*                generalized RQ factorization of the pair (B, A) is\n*                singular, so that rank(B) < P; the least squares\n*                solution could not be computed.\n*          = 2:  the (N-P) by (N-P) part of the upper trapezoidal factor\n*                T associated with A in the generalized RQ factorization\n*                of the pair (B, A) is singular, so that\n*                rank( (A) ) < N; the least squares solution could not\n*                    ( (B) )\n*                be computed.\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  x, work, info, a, b, c, d = NumRu::Lapack.dgglse( a, b, c, d, lwork, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rblapack_a = argv[0];
  rblapack_b = argv[1];
  rblapack_c = argv[2];
  rblapack_d = argv[3];
  rblapack_lwork = argv[4];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (1th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (1th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_a);
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_DFLOAT)
    rblapack_a = na_change_type(rblapack_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rblapack_a, doublereal*);
  if (!NA_IsNArray(rblapack_b))
    rb_raise(rb_eArgError, "b (2th argument) must be NArray");
  if (NA_RANK(rblapack_b) != 2)
    rb_raise(rb_eArgError, "rank of b (2th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_b) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of b must be the same as shape 1 of a");
  ldb = NA_SHAPE0(rblapack_b);
  if (NA_TYPE(rblapack_b) != NA_DFLOAT)
    rblapack_b = na_change_type(rblapack_b, NA_DFLOAT);
  b = NA_PTR_TYPE(rblapack_b, doublereal*);
  if (!NA_IsNArray(rblapack_c))
    rb_raise(rb_eArgError, "c (3th argument) must be NArray");
  if (NA_RANK(rblapack_c) != 1)
    rb_raise(rb_eArgError, "rank of c (3th argument) must be %d", 1);
  m = NA_SHAPE0(rblapack_c);
  if (NA_TYPE(rblapack_c) != NA_DFLOAT)
    rblapack_c = na_change_type(rblapack_c, NA_DFLOAT);
  c = NA_PTR_TYPE(rblapack_c, doublereal*);
  if (!NA_IsNArray(rblapack_d))
    rb_raise(rb_eArgError, "d (4th argument) must be NArray");
  if (NA_RANK(rblapack_d) != 1)
    rb_raise(rb_eArgError, "rank of d (4th argument) must be %d", 1);
  p = NA_SHAPE0(rblapack_d);
  if (NA_TYPE(rblapack_d) != NA_DFLOAT)
    rblapack_d = na_change_type(rblapack_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rblapack_d, doublereal*);
  lwork = NUM2INT(rblapack_lwork);
  {
    int shape[1];
    shape[0] = n;
    rblapack_x = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  x = NA_PTR_TYPE(rblapack_x, doublereal*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rblapack_work = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rblapack_work, doublereal*);
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
    shape[1] = n;
    rblapack_b_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rblapack_b_out__, doublereal*);
  MEMCPY(b_out__, b, doublereal, NA_TOTAL(rblapack_b));
  rblapack_b = rblapack_b_out__;
  b = b_out__;
  {
    int shape[1];
    shape[0] = m;
    rblapack_c_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  c_out__ = NA_PTR_TYPE(rblapack_c_out__, doublereal*);
  MEMCPY(c_out__, c, doublereal, NA_TOTAL(rblapack_c));
  rblapack_c = rblapack_c_out__;
  c = c_out__;
  {
    int shape[1];
    shape[0] = p;
    rblapack_d_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  d_out__ = NA_PTR_TYPE(rblapack_d_out__, doublereal*);
  MEMCPY(d_out__, d, doublereal, NA_TOTAL(rblapack_d));
  rblapack_d = rblapack_d_out__;
  d = d_out__;

  dgglse_(&m, &n, &p, a, &lda, b, &ldb, c, d, x, work, &lwork, &info);

  rblapack_info = INT2NUM(info);
  return rb_ary_new3(7, rblapack_x, rblapack_work, rblapack_info, rblapack_a, rblapack_b, rblapack_c, rblapack_d);
}

void
init_lapack_dgglse(VALUE mLapack){
  rb_define_module_function(mLapack, "dgglse", rblapack_dgglse, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
