#include "rb_lapack.h"

extern VOID zunml2_(char *side, char *trans, integer *m, integer *n, integer *k, doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *c, integer *ldc, doublecomplex *work, integer *info);

static VALUE
rb_zunml2(int argc, VALUE *argv, VALUE self){
  VALUE rb_side;
  char side; 
  VALUE rb_trans;
  char trans; 
  VALUE rb_a;
  doublecomplex *a; 
  VALUE rb_tau;
  doublecomplex *tau; 
  VALUE rb_c;
  doublecomplex *c; 
  VALUE rb_info;
  integer info; 
  VALUE rb_c_out__;
  doublecomplex *c_out__;
  doublecomplex *work;

  integer lda;
  integer m;
  integer k;
  integer ldc;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, c = NumRu::Lapack.zunml2( side, trans, a, tau, c)\n    or\n  NumRu::Lapack.zunml2  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZUNML2( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  ZUNML2 overwrites the general complex m-by-n matrix C with\n*\n*        Q * C  if SIDE = 'L' and TRANS = 'N', or\n*\n*        Q'* C  if SIDE = 'L' and TRANS = 'C', or\n*\n*        C * Q  if SIDE = 'R' and TRANS = 'N', or\n*\n*        C * Q' if SIDE = 'R' and TRANS = 'C',\n*\n*  where Q is a complex unitary matrix defined as the product of k\n*  elementary reflectors\n*\n*        Q = H(k)' . . . H(2)' H(1)'\n*\n*  as returned by ZGELQF. Q is of order m if SIDE = 'L' and of order n\n*  if SIDE = 'R'.\n*\n\n*  Arguments\n*  =========\n*\n*  SIDE    (input) CHARACTER*1\n*          = 'L': apply Q or Q' from the Left\n*          = 'R': apply Q or Q' from the Right\n*\n*  TRANS   (input) CHARACTER*1\n*          = 'N': apply Q  (No transpose)\n*          = 'C': apply Q' (Conjugate transpose)\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix C. M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix C. N >= 0.\n*\n*  K       (input) INTEGER\n*          The number of elementary reflectors whose product defines\n*          the matrix Q.\n*          If SIDE = 'L', M >= K >= 0;\n*          if SIDE = 'R', N >= K >= 0.\n*\n*  A       (input) COMPLEX*16 array, dimension\n*                               (LDA,M) if SIDE = 'L',\n*                               (LDA,N) if SIDE = 'R'\n*          The i-th row must contain the vector which defines the\n*          elementary reflector H(i), for i = 1,2,...,k, as returned by\n*          ZGELQF in the first k rows of its array argument A.\n*          A is modified by the routine but restored on exit.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A. LDA >= max(1,K).\n*\n*  TAU     (input) COMPLEX*16 array, dimension (K)\n*          TAU(i) must contain the scalar factor of the elementary\n*          reflector H(i), as returned by ZGELQF.\n*\n*  C       (input/output) COMPLEX*16 array, dimension (LDC,N)\n*          On entry, the m-by-n matrix C.\n*          On exit, C is overwritten by Q*C or Q'*C or C*Q' or C*Q.\n*\n*  LDC     (input) INTEGER\n*          The leading dimension of the array C. LDC >= max(1,M).\n*\n*  WORK    (workspace) COMPLEX*16 array, dimension\n*                                   (N) if SIDE = 'L',\n*                                   (M) if SIDE = 'R'\n*\n*  INFO    (output) INTEGER\n*          = 0: successful exit\n*          < 0: if INFO = -i, the i-th argument had an illegal value\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_side = argv[0];
  rb_trans = argv[1];
  rb_a = argv[2];
  rb_tau = argv[3];
  rb_c = argv[4];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  m = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_DCOMPLEX)
    rb_a = na_change_type(rb_a, NA_DCOMPLEX);
  a = NA_PTR_TYPE(rb_a, doublecomplex*);
  side = StringValueCStr(rb_side)[0];
  if (!NA_IsNArray(rb_c))
    rb_raise(rb_eArgError, "c (5th argument) must be NArray");
  if (NA_RANK(rb_c) != 2)
    rb_raise(rb_eArgError, "rank of c (5th argument) must be %d", 2);
  n = NA_SHAPE1(rb_c);
  ldc = NA_SHAPE0(rb_c);
  if (NA_TYPE(rb_c) != NA_DCOMPLEX)
    rb_c = na_change_type(rb_c, NA_DCOMPLEX);
  c = NA_PTR_TYPE(rb_c, doublecomplex*);
  if (!NA_IsNArray(rb_tau))
    rb_raise(rb_eArgError, "tau (4th argument) must be NArray");
  if (NA_RANK(rb_tau) != 1)
    rb_raise(rb_eArgError, "rank of tau (4th argument) must be %d", 1);
  k = NA_SHAPE0(rb_tau);
  if (NA_TYPE(rb_tau) != NA_DCOMPLEX)
    rb_tau = na_change_type(rb_tau, NA_DCOMPLEX);
  tau = NA_PTR_TYPE(rb_tau, doublecomplex*);
  trans = StringValueCStr(rb_trans)[0];
  {
    int shape[2];
    shape[0] = ldc;
    shape[1] = n;
    rb_c_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  c_out__ = NA_PTR_TYPE(rb_c_out__, doublecomplex*);
  MEMCPY(c_out__, c, doublecomplex, NA_TOTAL(rb_c));
  rb_c = rb_c_out__;
  c = c_out__;
  work = ALLOC_N(doublecomplex, (lsame_(&side,"L") ? n : lsame_(&side,"R") ? m : 0));

  zunml2_(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &info);

  free(work);
  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_info, rb_c);
}

void
init_lapack_zunml2(VALUE mLapack){
  rb_define_module_function(mLapack, "zunml2", rb_zunml2, -1);
}
