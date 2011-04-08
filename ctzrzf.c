#include "rb_lapack.h"

extern VOID ctzrzf_(integer *m, integer *n, complex *a, integer *lda, complex *tau, complex *work, integer *lwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_ctzrzf(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_a;
  complex *a; 
  VALUE rblapack_lwork;
  integer lwork; 
  VALUE rblapack_tau;
  complex *tau; 
  VALUE rblapack_work;
  complex *work; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_a_out__;
  complex *a_out__;

  integer lda;
  integer n;
  integer m;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  tau, work, info, a = NumRu::Lapack.ctzrzf( a, lwork, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE CTZRZF( M, N, A, LDA, TAU, WORK, LWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  CTZRZF reduces the M-by-N ( M<=N ) complex upper trapezoidal matrix A\n*  to upper triangular form by means of unitary transformations.\n*\n*  The upper trapezoidal matrix A is factored as\n*\n*     A = ( R  0 ) * Z,\n*\n*  where Z is an N-by-N unitary matrix and R is an M-by-M upper\n*  triangular matrix.\n*\n\n*  Arguments\n*  =========\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix A.  M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix A.  N >= M.\n*\n*  A       (input/output) COMPLEX array, dimension (LDA,N)\n*          On entry, the leading M-by-N upper trapezoidal part of the\n*          array A must contain the matrix to be factorized.\n*          On exit, the leading M-by-M upper triangular part of A\n*          contains the upper triangular matrix R, and elements M+1 to\n*          N of the first M rows of A, with the array TAU, represent the\n*          unitary matrix Z as a product of M elementary reflectors.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,M).\n*\n*  TAU     (output) COMPLEX array, dimension (M)\n*          The scalar factors of the elementary reflectors.\n*\n*  WORK    (workspace/output) COMPLEX array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK.  LWORK >= max(1,M).\n*          For optimum performance LWORK >= M*NB, where NB is\n*          the optimal blocksize.\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*    A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville, USA\n*\n*  The factorization is obtained by Householder's method.  The kth\n*  transformation matrix, Z( k ), which is used to introduce zeros into\n*  the ( m - k + 1 )th row of A, is given in the form\n*\n*     Z( k ) = ( I     0   ),\n*              ( 0  T( k ) )\n*\n*  where\n*\n*     T( k ) = I - tau*u( k )*u( k )',   u( k ) = (   1    ),\n*                                                 (   0    )\n*                                                 ( z( k ) )\n*\n*  tau is a scalar and z( k ) is an ( n - m ) element vector.\n*  tau and z( k ) are chosen to annihilate the elements of the kth row\n*  of X.\n*\n*  The scalar tau is returned in the kth element of TAU and the vector\n*  u( k ) in the kth row of A, such that the elements of z( k ) are\n*  in  a( k, m + 1 ), ..., a( k, n ). The elements of R are returned in\n*  the upper triangular part of A.\n*\n*  Z is given by\n*\n*     Z =  Z( 1 ) * Z( 2 ) * ... * Z( m ).\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  tau, work, info, a = NumRu::Lapack.ctzrzf( a, lwork, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rblapack_a = argv[0];
  rblapack_lwork = argv[1];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (1th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (1th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_a);
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_SCOMPLEX)
    rblapack_a = na_change_type(rblapack_a, NA_SCOMPLEX);
  a = NA_PTR_TYPE(rblapack_a, complex*);
  lwork = NUM2INT(rblapack_lwork);
  m = lda;
  {
    int shape[1];
    shape[0] = m;
    rblapack_tau = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  tau = NA_PTR_TYPE(rblapack_tau, complex*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rblapack_work = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rblapack_work, complex*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rblapack_a_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rblapack_a_out__, complex*);
  MEMCPY(a_out__, a, complex, NA_TOTAL(rblapack_a));
  rblapack_a = rblapack_a_out__;
  a = a_out__;

  ctzrzf_(&m, &n, a, &lda, tau, work, &lwork, &info);

  rblapack_info = INT2NUM(info);
  return rb_ary_new3(4, rblapack_tau, rblapack_work, rblapack_info, rblapack_a);
}

void
init_lapack_ctzrzf(VALUE mLapack){
  rb_define_module_function(mLapack, "ctzrzf", rblapack_ctzrzf, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
