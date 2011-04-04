#include "rb_lapack.h"

extern VOID zupmtr_(char *side, char *uplo, char *trans, integer *m, integer *n, doublecomplex *ap, doublecomplex *tau, doublecomplex *c, integer *ldc, doublecomplex *work, integer *info);

static VALUE
rb_zupmtr(int argc, VALUE *argv, VALUE self){
  VALUE rb_side;
  char side; 
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_trans;
  char trans; 
  VALUE rb_m;
  integer m; 
  VALUE rb_ap;
  doublecomplex *ap; 
  VALUE rb_tau;
  doublecomplex *tau; 
  VALUE rb_c;
  doublecomplex *c; 
  VALUE rb_info;
  integer info; 
  VALUE rb_c_out__;
  doublecomplex *c_out__;
  doublecomplex *work;

  integer ldc;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, c = NumRu::Lapack.zupmtr( side, uplo, trans, m, ap, tau, c)\n    or\n  NumRu::Lapack.zupmtr  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZUPMTR( SIDE, UPLO, TRANS, M, N, AP, TAU, C, LDC, WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  ZUPMTR overwrites the general complex M-by-N matrix C with\n*\n*                  SIDE = 'L'     SIDE = 'R'\n*  TRANS = 'N':      Q * C          C * Q\n*  TRANS = 'C':      Q**H * C       C * Q**H\n*\n*  where Q is a complex unitary matrix of order nq, with nq = m if\n*  SIDE = 'L' and nq = n if SIDE = 'R'. Q is defined as the product of\n*  nq-1 elementary reflectors, as returned by ZHPTRD using packed\n*  storage:\n*\n*  if UPLO = 'U', Q = H(nq-1) . . . H(2) H(1);\n*\n*  if UPLO = 'L', Q = H(1) H(2) . . . H(nq-1).\n*\n\n*  Arguments\n*  =========\n*\n*  SIDE    (input) CHARACTER*1\n*          = 'L': apply Q or Q**H from the Left;\n*          = 'R': apply Q or Q**H from the Right.\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U': Upper triangular packed storage used in previous\n*                 call to ZHPTRD;\n*          = 'L': Lower triangular packed storage used in previous\n*                 call to ZHPTRD.\n*\n*  TRANS   (input) CHARACTER*1\n*          = 'N':  No transpose, apply Q;\n*          = 'C':  Conjugate transpose, apply Q**H.\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix C. M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix C. N >= 0.\n*\n*  AP      (input) COMPLEX*16 array, dimension\n*                               (M*(M+1)/2) if SIDE = 'L'\n*                               (N*(N+1)/2) if SIDE = 'R'\n*          The vectors which define the elementary reflectors, as\n*          returned by ZHPTRD.  AP is modified by the routine but\n*          restored on exit.\n*\n*  TAU     (input) COMPLEX*16 array, dimension (M-1) if SIDE = 'L'\n*                                     or (N-1) if SIDE = 'R'\n*          TAU(i) must contain the scalar factor of the elementary\n*          reflector H(i), as returned by ZHPTRD.\n*\n*  C       (input/output) COMPLEX*16 array, dimension (LDC,N)\n*          On entry, the M-by-N matrix C.\n*          On exit, C is overwritten by Q*C or Q**H*C or C*Q**H or C*Q.\n*\n*  LDC     (input) INTEGER\n*          The leading dimension of the array C. LDC >= max(1,M).\n*\n*  WORK    (workspace) COMPLEX*16 array, dimension\n*                                   (N) if SIDE = 'L'\n*                                   (M) if SIDE = 'R'\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 7)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 7)", argc);
  rb_side = argv[0];
  rb_uplo = argv[1];
  rb_trans = argv[2];
  rb_m = argv[3];
  rb_ap = argv[4];
  rb_tau = argv[5];
  rb_c = argv[6];

  side = StringValueCStr(rb_side)[0];
  m = NUM2INT(rb_m);
  if (!NA_IsNArray(rb_c))
    rb_raise(rb_eArgError, "c (7th argument) must be NArray");
  if (NA_RANK(rb_c) != 2)
    rb_raise(rb_eArgError, "rank of c (7th argument) must be %d", 2);
  n = NA_SHAPE1(rb_c);
  ldc = NA_SHAPE0(rb_c);
  if (NA_TYPE(rb_c) != NA_DCOMPLEX)
    rb_c = na_change_type(rb_c, NA_DCOMPLEX);
  c = NA_PTR_TYPE(rb_c, doublecomplex*);
  trans = StringValueCStr(rb_trans)[0];
  uplo = StringValueCStr(rb_uplo)[0];
  if (!NA_IsNArray(rb_tau))
    rb_raise(rb_eArgError, "tau (6th argument) must be NArray");
  if (NA_RANK(rb_tau) != 1)
    rb_raise(rb_eArgError, "rank of tau (6th argument) must be %d", 1);
  if (NA_SHAPE0(rb_tau) != (m-1))
    rb_raise(rb_eRuntimeError, "shape 0 of tau must be %d", m-1);
  if (NA_TYPE(rb_tau) != NA_DCOMPLEX)
    rb_tau = na_change_type(rb_tau, NA_DCOMPLEX);
  tau = NA_PTR_TYPE(rb_tau, doublecomplex*);
  if (!NA_IsNArray(rb_ap))
    rb_raise(rb_eArgError, "ap (5th argument) must be NArray");
  if (NA_RANK(rb_ap) != 1)
    rb_raise(rb_eArgError, "rank of ap (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_ap) != (m*(m+1)/2))
    rb_raise(rb_eRuntimeError, "shape 0 of ap must be %d", m*(m+1)/2);
  if (NA_TYPE(rb_ap) != NA_DCOMPLEX)
    rb_ap = na_change_type(rb_ap, NA_DCOMPLEX);
  ap = NA_PTR_TYPE(rb_ap, doublecomplex*);
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

  zupmtr_(&side, &uplo, &trans, &m, &n, ap, tau, c, &ldc, work, &info);

  free(work);
  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_info, rb_c);
}

void
init_lapack_zupmtr(VALUE mLapack){
  rb_define_module_function(mLapack, "zupmtr", rb_zupmtr, -1);
}
