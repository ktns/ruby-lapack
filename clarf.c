#include "rb_lapack.h"

extern VOID clarf_(char *side, integer *m, integer *n, complex *v, integer *incv, complex *tau, complex *c, integer *ldc, complex *work);

static VALUE sHelp, sUsage;

static VALUE
rblapack_clarf(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_side;
  char side; 
  VALUE rblapack_m;
  integer m; 
  VALUE rblapack_v;
  complex *v; 
  VALUE rblapack_incv;
  integer incv; 
  VALUE rblapack_tau;
  complex tau; 
  VALUE rblapack_c;
  complex *c; 
  VALUE rblapack_c_out__;
  complex *c_out__;
  complex *work;

  integer ldc;
  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  c = NumRu::Lapack.clarf( side, m, v, incv, tau, c, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE CLARF( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )\n\n*  Purpose\n*  =======\n*\n*  CLARF applies a complex elementary reflector H to a complex M-by-N\n*  matrix C, from either the left or the right. H is represented in the\n*  form\n*\n*        H = I - tau * v * v'\n*\n*  where tau is a complex scalar and v is a complex vector.\n*\n*  If tau = 0, then H is taken to be the unit matrix.\n*\n*  To apply H' (the conjugate transpose of H), supply conjg(tau) instead\n*  tau.\n*\n\n*  Arguments\n*  =========\n*\n*  SIDE    (input) CHARACTER*1\n*          = 'L': form  H * C\n*          = 'R': form  C * H\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix C.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix C.\n*\n*  V       (input) COMPLEX array, dimension\n*                     (1 + (M-1)*abs(INCV)) if SIDE = 'L'\n*                  or (1 + (N-1)*abs(INCV)) if SIDE = 'R'\n*          The vector v in the representation of H. V is not used if\n*          TAU = 0.\n*\n*  INCV    (input) INTEGER\n*          The increment between elements of v. INCV <> 0.\n*\n*  TAU     (input) COMPLEX\n*          The value tau in the representation of H.\n*\n*  C       (input/output) COMPLEX array, dimension (LDC,N)\n*          On entry, the M-by-N matrix C.\n*          On exit, C is overwritten by the matrix H * C if SIDE = 'L',\n*          or C * H if SIDE = 'R'.\n*\n*  LDC     (input) INTEGER\n*          The leading dimension of the array C. LDC >= max(1,M).\n*\n*  WORK    (workspace) COMPLEX array, dimension\n*                         (N) if SIDE = 'L'\n*                      or (M) if SIDE = 'R'\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  c = NumRu::Lapack.clarf( side, m, v, incv, tau, c, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rblapack_side = argv[0];
  rblapack_m = argv[1];
  rblapack_v = argv[2];
  rblapack_incv = argv[3];
  rblapack_tau = argv[4];
  rblapack_c = argv[5];
  if (rb_options != Qnil) {
  }

  tau.r = (real)NUM2DBL(rb_funcall(rblapack_tau, rb_intern("real"), 0));
  tau.i = (real)NUM2DBL(rb_funcall(rblapack_tau, rb_intern("imag"), 0));
  side = StringValueCStr(rblapack_side)[0];
  m = NUM2INT(rblapack_m);
  incv = NUM2INT(rblapack_incv);
  if (!NA_IsNArray(rblapack_c))
    rb_raise(rb_eArgError, "c (6th argument) must be NArray");
  if (NA_RANK(rblapack_c) != 2)
    rb_raise(rb_eArgError, "rank of c (6th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_c);
  ldc = NA_SHAPE0(rblapack_c);
  if (NA_TYPE(rblapack_c) != NA_SCOMPLEX)
    rblapack_c = na_change_type(rblapack_c, NA_SCOMPLEX);
  c = NA_PTR_TYPE(rblapack_c, complex*);
  if (!NA_IsNArray(rblapack_v))
    rb_raise(rb_eArgError, "v (3th argument) must be NArray");
  if (NA_RANK(rblapack_v) != 1)
    rb_raise(rb_eArgError, "rank of v (3th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_v) != (1 + (m-1)*abs(incv)))
    rb_raise(rb_eRuntimeError, "shape 0 of v must be %d", 1 + (m-1)*abs(incv));
  if (NA_TYPE(rblapack_v) != NA_SCOMPLEX)
    rblapack_v = na_change_type(rblapack_v, NA_SCOMPLEX);
  v = NA_PTR_TYPE(rblapack_v, complex*);
  {
    int shape[2];
    shape[0] = ldc;
    shape[1] = n;
    rblapack_c_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  c_out__ = NA_PTR_TYPE(rblapack_c_out__, complex*);
  MEMCPY(c_out__, c, complex, NA_TOTAL(rblapack_c));
  rblapack_c = rblapack_c_out__;
  c = c_out__;
  work = ALLOC_N(complex, (lsame_(&side,"L") ? n : lsame_(&side,"R") ? m : 0));

  clarf_(&side, &m, &n, v, &incv, &tau, c, &ldc, work);

  free(work);
  return rblapack_c;
}

void
init_lapack_clarf(VALUE mLapack){
  rb_define_module_function(mLapack, "clarf", rblapack_clarf, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
