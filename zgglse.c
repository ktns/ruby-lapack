#include "rb_lapack.h"

extern VOID zgglse_(integer *m, integer *n, integer *p, doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, doublecomplex *c, doublecomplex *d, doublecomplex *x, doublecomplex *work, integer *lwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_zgglse(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_a;
  doublecomplex *a; 
  VALUE rblapack_b;
  doublecomplex *b; 
  VALUE rblapack_c;
  doublecomplex *c; 
  VALUE rblapack_d;
  doublecomplex *d; 
  VALUE rblapack_lwork;
  integer lwork; 
  VALUE rblapack_x;
  doublecomplex *x; 
  VALUE rblapack_work;
  doublecomplex *work; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_a_out__;
  doublecomplex *a_out__;
  VALUE rblapack_b_out__;
  doublecomplex *b_out__;
  VALUE rblapack_c_out__;
  doublecomplex *c_out__;
  VALUE rblapack_d_out__;
  doublecomplex *d_out__;

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
      printf("%s\n", "USAGE:\n  x, work, info, a, b, c, d = NumRu::Lapack.zgglse( a, b, c, d, lwork, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZGGLSE( M, N, P, A, LDA, B, LDB, C, D, X, WORK, LWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  ZGGLSE solves the linear equality-constrained least squares (LSE)\n*  problem:\n*\n*          minimize || c - A*x ||_2   subject to   B*x = d\n*\n*  where A is an M-by-N matrix, B is a P-by-N matrix, c is a given\n*  M-vector, and d is a given P-vector. It is assumed that\n*  P <= N <= M+P, and\n*\n*           rank(B) = P and  rank( ( A ) ) = N.\n*                                ( ( B ) )\n*\n*  These conditions ensure that the LSE problem has a unique solution,\n*  which is obtained using a generalized RQ factorization of the\n*  matrices (B, A) given by\n*\n*     B = (0 R)*Q,   A = Z*T*Q.\n*\n\n*  Arguments\n*  =========\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix A.  M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrices A and B. N >= 0.\n*\n*  P       (input) INTEGER\n*          The number of rows of the matrix B. 0 <= P <= N <= M+P.\n*\n*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)\n*          On entry, the M-by-N matrix A.\n*          On exit, the elements on and above the diagonal of the array\n*          contain the min(M,N)-by-N upper trapezoidal matrix T.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A. LDA >= max(1,M).\n*\n*  B       (input/output) COMPLEX*16 array, dimension (LDB,N)\n*          On entry, the P-by-N matrix B.\n*          On exit, the upper triangle of the subarray B(1:P,N-P+1:N)\n*          contains the P-by-P upper triangular matrix R.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B. LDB >= max(1,P).\n*\n*  C       (input/output) COMPLEX*16 array, dimension (M)\n*          On entry, C contains the right hand side vector for the\n*          least squares part of the LSE problem.\n*          On exit, the residual sum of squares for the solution\n*          is given by the sum of squares of elements N-P+1 to M of\n*          vector C.\n*\n*  D       (input/output) COMPLEX*16 array, dimension (P)\n*          On entry, D contains the right hand side vector for the\n*          constrained equation.\n*          On exit, D is destroyed.\n*\n*  X       (output) COMPLEX*16 array, dimension (N)\n*          On exit, X is the solution of the LSE problem.\n*\n*  WORK    (workspace/output) COMPLEX*16 array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK. LWORK >= max(1,M+N+P).\n*          For optimum performance LWORK >= P+min(M,N)+max(M,N)*NB,\n*          where NB is an upper bound for the optimal blocksizes for\n*          ZGEQRF, CGERQF, ZUNMQR and CUNMRQ.\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*          = 1:  the upper triangular factor R associated with B in the\n*                generalized RQ factorization of the pair (B, A) is\n*                singular, so that rank(B) < P; the least squares\n*                solution could not be computed.\n*          = 2:  the (N-P) by (N-P) part of the upper trapezoidal factor\n*                T associated with A in the generalized RQ factorization\n*                of the pair (B, A) is singular, so that\n*                rank( (A) ) < N; the least squares solution could not\n*                    ( (B) )\n*                be computed.\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  x, work, info, a, b, c, d = NumRu::Lapack.zgglse( a, b, c, d, lwork, [:usage => usage, :help => help])\n");
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
  if (NA_TYPE(rblapack_a) != NA_DCOMPLEX)
    rblapack_a = na_change_type(rblapack_a, NA_DCOMPLEX);
  a = NA_PTR_TYPE(rblapack_a, doublecomplex*);
  if (!NA_IsNArray(rblapack_b))
    rb_raise(rb_eArgError, "b (2th argument) must be NArray");
  if (NA_RANK(rblapack_b) != 2)
    rb_raise(rb_eArgError, "rank of b (2th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_b) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of b must be the same as shape 1 of a");
  ldb = NA_SHAPE0(rblapack_b);
  if (NA_TYPE(rblapack_b) != NA_DCOMPLEX)
    rblapack_b = na_change_type(rblapack_b, NA_DCOMPLEX);
  b = NA_PTR_TYPE(rblapack_b, doublecomplex*);
  if (!NA_IsNArray(rblapack_c))
    rb_raise(rb_eArgError, "c (3th argument) must be NArray");
  if (NA_RANK(rblapack_c) != 1)
    rb_raise(rb_eArgError, "rank of c (3th argument) must be %d", 1);
  m = NA_SHAPE0(rblapack_c);
  if (NA_TYPE(rblapack_c) != NA_DCOMPLEX)
    rblapack_c = na_change_type(rblapack_c, NA_DCOMPLEX);
  c = NA_PTR_TYPE(rblapack_c, doublecomplex*);
  if (!NA_IsNArray(rblapack_d))
    rb_raise(rb_eArgError, "d (4th argument) must be NArray");
  if (NA_RANK(rblapack_d) != 1)
    rb_raise(rb_eArgError, "rank of d (4th argument) must be %d", 1);
  p = NA_SHAPE0(rblapack_d);
  if (NA_TYPE(rblapack_d) != NA_DCOMPLEX)
    rblapack_d = na_change_type(rblapack_d, NA_DCOMPLEX);
  d = NA_PTR_TYPE(rblapack_d, doublecomplex*);
  lwork = NUM2INT(rblapack_lwork);
  {
    int shape[1];
    shape[0] = n;
    rblapack_x = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  x = NA_PTR_TYPE(rblapack_x, doublecomplex*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rblapack_work = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rblapack_work, doublecomplex*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rblapack_a_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rblapack_a_out__, doublecomplex*);
  MEMCPY(a_out__, a, doublecomplex, NA_TOTAL(rblapack_a));
  rblapack_a = rblapack_a_out__;
  a = a_out__;
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = n;
    rblapack_b_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rblapack_b_out__, doublecomplex*);
  MEMCPY(b_out__, b, doublecomplex, NA_TOTAL(rblapack_b));
  rblapack_b = rblapack_b_out__;
  b = b_out__;
  {
    int shape[1];
    shape[0] = m;
    rblapack_c_out__ = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  c_out__ = NA_PTR_TYPE(rblapack_c_out__, doublecomplex*);
  MEMCPY(c_out__, c, doublecomplex, NA_TOTAL(rblapack_c));
  rblapack_c = rblapack_c_out__;
  c = c_out__;
  {
    int shape[1];
    shape[0] = p;
    rblapack_d_out__ = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  d_out__ = NA_PTR_TYPE(rblapack_d_out__, doublecomplex*);
  MEMCPY(d_out__, d, doublecomplex, NA_TOTAL(rblapack_d));
  rblapack_d = rblapack_d_out__;
  d = d_out__;

  zgglse_(&m, &n, &p, a, &lda, b, &ldb, c, d, x, work, &lwork, &info);

  rblapack_info = INT2NUM(info);
  return rb_ary_new3(7, rblapack_x, rblapack_work, rblapack_info, rblapack_a, rblapack_b, rblapack_c, rblapack_d);
}

void
init_lapack_zgglse(VALUE mLapack){
  rb_define_module_function(mLapack, "zgglse", rblapack_zgglse, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
