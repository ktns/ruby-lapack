#include "rb_lapack.h"

extern VOID zgels_(char *trans, integer *m, integer *n, integer *nrhs, doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, doublecomplex *work, integer *lwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_zgels(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_trans;
  char trans; 
  VALUE rblapack_m;
  integer m; 
  VALUE rblapack_a;
  doublecomplex *a; 
  VALUE rblapack_b;
  doublecomplex *b; 
  VALUE rblapack_lwork;
  integer lwork; 
  VALUE rblapack_work;
  doublecomplex *work; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_a_out__;
  doublecomplex *a_out__;
  VALUE rblapack_b_out__;
  doublecomplex *b_out__;

  integer lda;
  integer n;
  integer ldb;
  integer nrhs;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  work, info, a, b = NumRu::Lapack.zgels( trans, m, a, b, lwork, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZGELS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  ZGELS solves overdetermined or underdetermined complex linear systems\n*  involving an M-by-N matrix A, or its conjugate-transpose, using a QR\n*  or LQ factorization of A.  It is assumed that A has full rank.\n*\n*  The following options are provided:\n*\n*  1. If TRANS = 'N' and m >= n:  find the least squares solution of\n*     an overdetermined system, i.e., solve the least squares problem\n*                  minimize || B - A*X ||.\n*\n*  2. If TRANS = 'N' and m < n:  find the minimum norm solution of\n*     an underdetermined system A * X = B.\n*\n*  3. If TRANS = 'C' and m >= n:  find the minimum norm solution of\n*     an undetermined system A**H * X = B.\n*\n*  4. If TRANS = 'C' and m < n:  find the least squares solution of\n*     an overdetermined system, i.e., solve the least squares problem\n*                  minimize || B - A**H * X ||.\n*\n*  Several right hand side vectors b and solution vectors x can be\n*  handled in a single call; they are stored as the columns of the\n*  M-by-NRHS right hand side matrix B and the N-by-NRHS solution\n*  matrix X.\n*\n\n*  Arguments\n*  =========\n*\n*  TRANS   (input) CHARACTER*1\n*          = 'N': the linear system involves A;\n*          = 'C': the linear system involves A**H.\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix A.  M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix A.  N >= 0.\n*\n*  NRHS    (input) INTEGER\n*          The number of right hand sides, i.e., the number of\n*          columns of the matrices B and X. NRHS >= 0.\n*\n*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)\n*          On entry, the M-by-N matrix A.\n*            if M >= N, A is overwritten by details of its QR\n*                       factorization as returned by ZGEQRF;\n*            if M <  N, A is overwritten by details of its LQ\n*                       factorization as returned by ZGELQF.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,M).\n*\n*  B       (input/output) COMPLEX*16 array, dimension (LDB,NRHS)\n*          On entry, the matrix B of right hand side vectors, stored\n*          columnwise; B is M-by-NRHS if TRANS = 'N', or N-by-NRHS\n*          if TRANS = 'C'.\n*          On exit, if INFO = 0, B is overwritten by the solution\n*          vectors, stored columnwise:\n*          if TRANS = 'N' and m >= n, rows 1 to n of B contain the least\n*          squares solution vectors; the residual sum of squares for the\n*          solution in each column is given by the sum of squares of the\n*          modulus of elements N+1 to M in that column;\n*          if TRANS = 'N' and m < n, rows 1 to N of B contain the\n*          minimum norm solution vectors;\n*          if TRANS = 'C' and m >= n, rows 1 to M of B contain the\n*          minimum norm solution vectors;\n*          if TRANS = 'C' and m < n, rows 1 to M of B contain the\n*          least squares solution vectors; the residual sum of squares\n*          for the solution in each column is given by the sum of\n*          squares of the modulus of elements M+1 to N in that column.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B. LDB >= MAX(1,M,N).\n*\n*  WORK    (workspace/output) COMPLEX*16 array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK.\n*          LWORK >= max( 1, MN + max( MN, NRHS ) ).\n*          For optimal performance,\n*          LWORK >= max( 1, MN + max( MN, NRHS )*NB ).\n*          where MN = min(M,N) and NB is the optimum block size.\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*          > 0:  if INFO =  i, the i-th diagonal element of the\n*                triangular factor of A is zero, so that A does not have\n*                full rank; the least squares solution could not be\n*                computed.\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  work, info, a, b = NumRu::Lapack.zgels( trans, m, a, b, lwork, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rblapack_trans = argv[0];
  rblapack_m = argv[1];
  rblapack_a = argv[2];
  rblapack_b = argv[3];
  rblapack_lwork = argv[4];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_a);
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_DCOMPLEX)
    rblapack_a = na_change_type(rblapack_a, NA_DCOMPLEX);
  a = NA_PTR_TYPE(rblapack_a, doublecomplex*);
  if (!NA_IsNArray(rblapack_b))
    rb_raise(rb_eArgError, "b (4th argument) must be NArray");
  if (NA_RANK(rblapack_b) != 2)
    rb_raise(rb_eArgError, "rank of b (4th argument) must be %d", 2);
  nrhs = NA_SHAPE1(rblapack_b);
  ldb = NA_SHAPE0(rblapack_b);
  if (NA_TYPE(rblapack_b) != NA_DCOMPLEX)
    rblapack_b = na_change_type(rblapack_b, NA_DCOMPLEX);
  b = NA_PTR_TYPE(rblapack_b, doublecomplex*);
  m = NUM2INT(rblapack_m);
  lwork = NUM2INT(rblapack_lwork);
  trans = StringValueCStr(rblapack_trans)[0];
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
    shape[1] = nrhs;
    rblapack_b_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rblapack_b_out__, doublecomplex*);
  MEMCPY(b_out__, b, doublecomplex, NA_TOTAL(rblapack_b));
  rblapack_b = rblapack_b_out__;
  b = b_out__;

  zgels_(&trans, &m, &n, &nrhs, a, &lda, b, &ldb, work, &lwork, &info);

  rblapack_info = INT2NUM(info);
  return rb_ary_new3(4, rblapack_work, rblapack_info, rblapack_a, rblapack_b);
}

void
init_lapack_zgels(VALUE mLapack){
  rb_define_module_function(mLapack, "zgels", rblapack_zgels, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
