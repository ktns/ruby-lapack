#include "rb_lapack.h"

extern VOID zlaic1_(integer *job, integer *j, doublecomplex *x, doublereal *sest, doublecomplex *w, doublecomplex *gamma, doublereal *sestpr, doublecomplex *s, doublecomplex *c);

static VALUE sHelp, sUsage;

static VALUE
rblapack_zlaic1(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_job;
  integer job; 
  VALUE rblapack_x;
  doublecomplex *x; 
  VALUE rblapack_sest;
  doublereal sest; 
  VALUE rblapack_w;
  doublecomplex *w; 
  VALUE rblapack_gamma;
  doublecomplex gamma; 
  VALUE rblapack_sestpr;
  doublereal sestpr; 
  VALUE rblapack_s;
  doublecomplex s; 
  VALUE rblapack_c;
  doublecomplex c; 

  integer j;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  sestpr, s, c = NumRu::Lapack.zlaic1( job, x, sest, w, gamma, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZLAIC1( JOB, J, X, SEST, W, GAMMA, SESTPR, S, C )\n\n*  Purpose\n*  =======\n*\n*  ZLAIC1 applies one step of incremental condition estimation in\n*  its simplest version:\n*\n*  Let x, twonorm(x) = 1, be an approximate singular vector of an j-by-j\n*  lower triangular matrix L, such that\n*           twonorm(L*x) = sest\n*  Then ZLAIC1 computes sestpr, s, c such that\n*  the vector\n*                  [ s*x ]\n*           xhat = [  c  ]\n*  is an approximate singular vector of\n*                  [ L     0  ]\n*           Lhat = [ w' gamma ]\n*  in the sense that\n*           twonorm(Lhat*xhat) = sestpr.\n*\n*  Depending on JOB, an estimate for the largest or smallest singular\n*  value is computed.\n*\n*  Note that [s c]' and sestpr**2 is an eigenpair of the system\n*\n*      diag(sest*sest, 0) + [alpha  gamma] * [ conjg(alpha) ]\n*                                            [ conjg(gamma) ]\n*\n*  where  alpha =  conjg(x)'*w.\n*\n\n*  Arguments\n*  =========\n*\n*  JOB     (input) INTEGER\n*          = 1: an estimate for the largest singular value is computed.\n*          = 2: an estimate for the smallest singular value is computed.\n*\n*  J       (input) INTEGER\n*          Length of X and W\n*\n*  X       (input) COMPLEX*16 array, dimension (J)\n*          The j-vector x.\n*\n*  SEST    (input) DOUBLE PRECISION\n*          Estimated singular value of j by j matrix L\n*\n*  W       (input) COMPLEX*16 array, dimension (J)\n*          The j-vector w.\n*\n*  GAMMA   (input) COMPLEX*16\n*          The diagonal element gamma.\n*\n*  SESTPR  (output) DOUBLE PRECISION\n*          Estimated singular value of (j+1) by (j+1) matrix Lhat.\n*\n*  S       (output) COMPLEX*16\n*          Sine needed in forming xhat.\n*\n*  C       (output) COMPLEX*16\n*          Cosine needed in forming xhat.\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  sestpr, s, c = NumRu::Lapack.zlaic1( job, x, sest, w, gamma, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rblapack_job = argv[0];
  rblapack_x = argv[1];
  rblapack_sest = argv[2];
  rblapack_w = argv[3];
  rblapack_gamma = argv[4];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_w))
    rb_raise(rb_eArgError, "w (4th argument) must be NArray");
  if (NA_RANK(rblapack_w) != 1)
    rb_raise(rb_eArgError, "rank of w (4th argument) must be %d", 1);
  j = NA_SHAPE0(rblapack_w);
  if (NA_TYPE(rblapack_w) != NA_DCOMPLEX)
    rblapack_w = na_change_type(rblapack_w, NA_DCOMPLEX);
  w = NA_PTR_TYPE(rblapack_w, doublecomplex*);
  if (!NA_IsNArray(rblapack_x))
    rb_raise(rb_eArgError, "x (2th argument) must be NArray");
  if (NA_RANK(rblapack_x) != 1)
    rb_raise(rb_eArgError, "rank of x (2th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_x) != j)
    rb_raise(rb_eRuntimeError, "shape 0 of x must be the same as shape 0 of w");
  if (NA_TYPE(rblapack_x) != NA_DCOMPLEX)
    rblapack_x = na_change_type(rblapack_x, NA_DCOMPLEX);
  x = NA_PTR_TYPE(rblapack_x, doublecomplex*);
  gamma.r = NUM2DBL(rb_funcall(rblapack_gamma, rb_intern("real"), 0));
  gamma.i = NUM2DBL(rb_funcall(rblapack_gamma, rb_intern("imag"), 0));
  job = NUM2INT(rblapack_job);
  sest = NUM2DBL(rblapack_sest);

  zlaic1_(&job, &j, x, &sest, w, &gamma, &sestpr, &s, &c);

  rblapack_sestpr = rb_float_new((double)sestpr);
  rblapack_s = rb_funcall(rb_gv_get("Complex"), rb_intern("new"), 2, rb_float_new((double)(s.r)), rb_float_new((double)(s.i)));
  rblapack_c = rb_funcall(rb_gv_get("Complex"), rb_intern("new"), 2, rb_float_new((double)(c.r)), rb_float_new((double)(c.i)));
  return rb_ary_new3(3, rblapack_sestpr, rblapack_s, rblapack_c);
}

void
init_lapack_zlaic1(VALUE mLapack){
  rb_define_module_function(mLapack, "zlaic1", rblapack_zlaic1, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
