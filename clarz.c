#include "rb_lapack.h"

static VALUE
rb_clarz(int argc, VALUE *argv, VALUE self){
  VALUE rb_side;
  char side; 
  VALUE rb_m;
  integer m; 
  VALUE rb_l;
  integer l; 
  VALUE rb_v;
  complex *v; 
  VALUE rb_incv;
  integer incv; 
  VALUE rb_tau;
  complex tau; 
  VALUE rb_c;
  complex *c; 
  VALUE rb_c_out__;
  complex *c_out__;
  complex *work;

  integer ldc;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  c = NumRu::Lapack.clarz( side, m, l, v, incv, tau, c)\n    or\n  NumRu::Lapack.clarz  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE CLARZ( SIDE, M, N, L, V, INCV, TAU, C, LDC, WORK )\n\n*  Purpose\n*  =======\n*\n*  CLARZ applies a complex elementary reflector H to a complex\n*  M-by-N matrix C, from either the left or the right. H is represented\n*  in the form\n*\n*        H = I - tau * v * v'\n*\n*  where tau is a complex scalar and v is a complex vector.\n*\n*  If tau = 0, then H is taken to be the unit matrix.\n*\n*  To apply H' (the conjugate transpose of H), supply conjg(tau) instead\n*  tau.\n*\n*  H is a product of k elementary reflectors as returned by CTZRZF.\n*\n\n*  Arguments\n*  =========\n*\n*  SIDE    (input) CHARACTER*1\n*          = 'L': form  H * C\n*          = 'R': form  C * H\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix C.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix C.\n*\n*  L       (input) INTEGER\n*          The number of entries of the vector V containing\n*          the meaningful part of the Householder vectors.\n*          If SIDE = 'L', M >= L >= 0, if SIDE = 'R', N >= L >= 0.\n*\n*  V       (input) COMPLEX array, dimension (1+(L-1)*abs(INCV))\n*          The vector v in the representation of H as returned by\n*          CTZRZF. V is not used if TAU = 0.\n*\n*  INCV    (input) INTEGER\n*          The increment between elements of v. INCV <> 0.\n*\n*  TAU     (input) COMPLEX\n*          The value tau in the representation of H.\n*\n*  C       (input/output) COMPLEX array, dimension (LDC,N)\n*          On entry, the M-by-N matrix C.\n*          On exit, C is overwritten by the matrix H * C if SIDE = 'L',\n*          or C * H if SIDE = 'R'.\n*\n*  LDC     (input) INTEGER\n*          The leading dimension of the array C. LDC >= max(1,M).\n*\n*  WORK    (workspace) COMPLEX array, dimension\n*                         (N) if SIDE = 'L'\n*                      or (M) if SIDE = 'R'\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*    A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville, USA\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 7)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 7)", argc);
  rb_side = argv[0];
  rb_m = argv[1];
  rb_l = argv[2];
  rb_v = argv[3];
  rb_incv = argv[4];
  rb_tau = argv[5];
  rb_c = argv[6];

  side = StringValueCStr(rb_side)[0];
  m = NUM2INT(rb_m);
  l = NUM2INT(rb_l);
  incv = NUM2INT(rb_incv);
  tau.r = (real)NUM2DBL(rb_funcall(rb_tau, rb_intern("real"), 0));
  tau.i = (real)NUM2DBL(rb_funcall(rb_tau, rb_intern("imag"), 0));
  if (!NA_IsNArray(rb_v))
    rb_raise(rb_eArgError, "v (4th argument) must be NArray");
  if (NA_RANK(rb_v) != 1)
    rb_raise(rb_eArgError, "rank of v (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_v) != (1+(l-1)*abs(incv)))
    rb_raise(rb_eRuntimeError, "shape 0 of v must be %d", 1+(l-1)*abs(incv));
  if (NA_TYPE(rb_v) != NA_SCOMPLEX)
    rb_v = na_change_type(rb_v, NA_SCOMPLEX);
  v = NA_PTR_TYPE(rb_v, complex*);
  if (!NA_IsNArray(rb_c))
    rb_raise(rb_eArgError, "c (7th argument) must be NArray");
  if (NA_RANK(rb_c) != 2)
    rb_raise(rb_eArgError, "rank of c (7th argument) must be %d", 2);
  ldc = NA_SHAPE0(rb_c);
  n = NA_SHAPE1(rb_c);
  if (NA_TYPE(rb_c) != NA_SCOMPLEX)
    rb_c = na_change_type(rb_c, NA_SCOMPLEX);
  c = NA_PTR_TYPE(rb_c, complex*);
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
  work = ALLOC_N(complex, (lsame_(&side,"L") ? n : lsame_(&side,"R") ? m : 0));

  clarz_(&side, &m, &n, &l, v, &incv, &tau, c, &ldc, work);

  free(work);
  return rb_c;
}

void
init_lapack_clarz(VALUE mLapack){
  rb_define_module_function(mLapack, "clarz", rb_clarz, -1);
}
