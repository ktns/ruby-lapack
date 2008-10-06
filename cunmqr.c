#include "rb_lapack.h"

static VALUE
rb_cunmqr(int argc, VALUE *argv, VALUE self){
  VALUE rb_side;
  char side; 
  VALUE rb_trans;
  char trans; 
  VALUE rb_m;
  integer m; 
  VALUE rb_a;
  complex *a; 
  VALUE rb_tau;
  complex *tau; 
  VALUE rb_c;
  complex *c; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_work;
  complex *work; 
  VALUE rb_info;
  integer info; 
  VALUE rb_c_out__;
  complex *c_out__;

  integer lda;
  integer k;
  integer ldc;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  work, info, c = NumRu::Lapack.cunmqr( side, trans, m, a, tau, c, lwork)\n    or\n  NumRu::Lapack.cunmqr  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE CUNMQR( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, LWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  CUNMQR overwrites the general complex M-by-N matrix C with\n*\n*                  SIDE = 'L'     SIDE = 'R'\n*  TRANS = 'N':      Q * C          C * Q\n*  TRANS = 'C':      Q**H * C       C * Q**H\n*\n*  where Q is a complex unitary matrix defined as the product of k\n*  elementary reflectors\n*\n*        Q = H(1) H(2) . . . H(k)\n*\n*  as returned by CGEQRF. Q is of order M if SIDE = 'L' and of order N\n*  if SIDE = 'R'.\n*\n\n*  Arguments\n*  =========\n*\n*  SIDE    (input) CHARACTER*1\n*          = 'L': apply Q or Q**H from the Left;\n*          = 'R': apply Q or Q**H from the Right.\n*\n*  TRANS   (input) CHARACTER*1\n*          = 'N':  No transpose, apply Q;\n*          = 'C':  Conjugate transpose, apply Q**H.\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix C. M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix C. N >= 0.\n*\n*  K       (input) INTEGER\n*          The number of elementary reflectors whose product defines\n*          the matrix Q.\n*          If SIDE = 'L', M >= K >= 0;\n*          if SIDE = 'R', N >= K >= 0.\n*\n*  A       (input) COMPLEX array, dimension (LDA,K)\n*          The i-th column must contain the vector which defines the\n*          elementary reflector H(i), for i = 1,2,...,k, as returned by\n*          CGEQRF in the first k columns of its array argument A.\n*          A is modified by the routine but restored on exit.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.\n*          If SIDE = 'L', LDA >= max(1,M);\n*          if SIDE = 'R', LDA >= max(1,N).\n*\n*  TAU     (input) COMPLEX array, dimension (K)\n*          TAU(i) must contain the scalar factor of the elementary\n*          reflector H(i), as returned by CGEQRF.\n*\n*  C       (input/output) COMPLEX array, dimension (LDC,N)\n*          On entry, the M-by-N matrix C.\n*          On exit, C is overwritten by Q*C or Q**H*C or C*Q**H or C*Q.\n*\n*  LDC     (input) INTEGER\n*          The leading dimension of the array C. LDC >= max(1,M).\n*\n*  WORK    (workspace/output) COMPLEX array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK.\n*          If SIDE = 'L', LWORK >= max(1,N);\n*          if SIDE = 'R', LWORK >= max(1,M).\n*          For optimum performance LWORK >= N*NB if SIDE = 'L', and\n*          LWORK >= M*NB if SIDE = 'R', where NB is the optimal\n*          blocksize.\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 7)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 7)", argc);
  rb_side = argv[0];
  rb_trans = argv[1];
  rb_m = argv[2];
  rb_a = argv[3];
  rb_tau = argv[4];
  rb_c = argv[5];
  rb_lwork = argv[6];

  side = StringValueCStr(rb_side)[0];
  trans = StringValueCStr(rb_trans)[0];
  m = NUM2INT(rb_m);
  lwork = NUM2INT(rb_lwork);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (4th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (4th argument) must be %d", 2);
  lda = NA_SHAPE0(rb_a);
  k = NA_SHAPE1(rb_a);
  if (NA_TYPE(rb_a) != NA_SCOMPLEX)
    rb_a = na_change_type(rb_a, NA_SCOMPLEX);
  a = NA_PTR_TYPE(rb_a, complex*);
  if (!NA_IsNArray(rb_tau))
    rb_raise(rb_eArgError, "tau (5th argument) must be NArray");
  if (NA_RANK(rb_tau) != 1)
    rb_raise(rb_eArgError, "rank of tau (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_tau) != k)
    rb_raise(rb_eRuntimeError, "shape 0 of tau must be the same as shape 1 of a");
  if (NA_TYPE(rb_tau) != NA_SCOMPLEX)
    rb_tau = na_change_type(rb_tau, NA_SCOMPLEX);
  tau = NA_PTR_TYPE(rb_tau, complex*);
  if (!NA_IsNArray(rb_c))
    rb_raise(rb_eArgError, "c (6th argument) must be NArray");
  if (NA_RANK(rb_c) != 2)
    rb_raise(rb_eArgError, "rank of c (6th argument) must be %d", 2);
  ldc = NA_SHAPE0(rb_c);
  n = NA_SHAPE1(rb_c);
  if (NA_TYPE(rb_c) != NA_SCOMPLEX)
    rb_c = na_change_type(rb_c, NA_SCOMPLEX);
  c = NA_PTR_TYPE(rb_c, complex*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rb_work = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rb_work, complex*);
  {
    int shape[2];
    shape[0] = ldc;
    shape[1] = n;
    rb_c_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  c_out__ = NA_PTR_TYPE(rb_c_out__, complex*);
  MEMCPY(c_out__, c, complex, NA_TOTAL(rb_c));
  rb_c = rb_c_out__;
  c = c_out__;

  cunmqr_(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(3, rb_work, rb_info, rb_c);
}

void
init_lapack_cunmqr(VALUE mLapack){
  rb_define_module_function(mLapack, "cunmqr", rb_cunmqr, -1);
}
