#include "rb_lapack.h"

extern VOID slaic1_(integer *job, integer *j, real *x, real *sest, real *w, real *gamma, real *sestpr, real *s, real *c);

static VALUE
rb_slaic1(int argc, VALUE *argv, VALUE self){
  VALUE rb_job;
  integer job; 
  VALUE rb_x;
  real *x; 
  VALUE rb_sest;
  real sest; 
  VALUE rb_w;
  real *w; 
  VALUE rb_gamma;
  real gamma; 
  VALUE rb_sestpr;
  real sestpr; 
  VALUE rb_s;
  real s; 
  VALUE rb_c;
  real c; 

  integer j;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  sestpr, s, c = NumRu::Lapack.slaic1( job, x, sest, w, gamma)\n    or\n  NumRu::Lapack.slaic1  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLAIC1( JOB, J, X, SEST, W, GAMMA, SESTPR, S, C )\n\n*  Purpose\n*  =======\n*\n*  SLAIC1 applies one step of incremental condition estimation in\n*  its simplest version:\n*\n*  Let x, twonorm(x) = 1, be an approximate singular vector of an j-by-j\n*  lower triangular matrix L, such that\n*           twonorm(L*x) = sest\n*  Then SLAIC1 computes sestpr, s, c such that\n*  the vector\n*                  [ s*x ]\n*           xhat = [  c  ]\n*  is an approximate singular vector of\n*                  [ L     0  ]\n*           Lhat = [ w' gamma ]\n*  in the sense that\n*           twonorm(Lhat*xhat) = sestpr.\n*\n*  Depending on JOB, an estimate for the largest or smallest singular\n*  value is computed.\n*\n*  Note that [s c]' and sestpr**2 is an eigenpair of the system\n*\n*      diag(sest*sest, 0) + [alpha  gamma] * [ alpha ]\n*                                            [ gamma ]\n*\n*  where  alpha =  x'*w.\n*\n\n*  Arguments\n*  =========\n*\n*  JOB     (input) INTEGER\n*          = 1: an estimate for the largest singular value is computed.\n*          = 2: an estimate for the smallest singular value is computed.\n*\n*  J       (input) INTEGER\n*          Length of X and W\n*\n*  X       (input) REAL array, dimension (J)\n*          The j-vector x.\n*\n*  SEST    (input) REAL\n*          Estimated singular value of j by j matrix L\n*\n*  W       (input) REAL array, dimension (J)\n*          The j-vector w.\n*\n*  GAMMA   (input) REAL\n*          The diagonal element gamma.\n*\n*  SESTPR  (output) REAL\n*          Estimated singular value of (j+1) by (j+1) matrix Lhat.\n*\n*  S       (output) REAL\n*          Sine needed in forming xhat.\n*\n*  C       (output) REAL\n*          Cosine needed in forming xhat.\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_job = argv[0];
  rb_x = argv[1];
  rb_sest = argv[2];
  rb_w = argv[3];
  rb_gamma = argv[4];

  if (!NA_IsNArray(rb_w))
    rb_raise(rb_eArgError, "w (4th argument) must be NArray");
  if (NA_RANK(rb_w) != 1)
    rb_raise(rb_eArgError, "rank of w (4th argument) must be %d", 1);
  j = NA_SHAPE0(rb_w);
  if (NA_TYPE(rb_w) != NA_SFLOAT)
    rb_w = na_change_type(rb_w, NA_SFLOAT);
  w = NA_PTR_TYPE(rb_w, real*);
  if (!NA_IsNArray(rb_x))
    rb_raise(rb_eArgError, "x (2th argument) must be NArray");
  if (NA_RANK(rb_x) != 1)
    rb_raise(rb_eArgError, "rank of x (2th argument) must be %d", 1);
  if (NA_SHAPE0(rb_x) != j)
    rb_raise(rb_eRuntimeError, "shape 0 of x must be the same as shape 0 of w");
  if (NA_TYPE(rb_x) != NA_SFLOAT)
    rb_x = na_change_type(rb_x, NA_SFLOAT);
  x = NA_PTR_TYPE(rb_x, real*);
  gamma = (real)NUM2DBL(rb_gamma);
  job = NUM2INT(rb_job);
  sest = (real)NUM2DBL(rb_sest);

  slaic1_(&job, &j, x, &sest, w, &gamma, &sestpr, &s, &c);

  rb_sestpr = rb_float_new((double)sestpr);
  rb_s = rb_float_new((double)s);
  rb_c = rb_float_new((double)c);
  return rb_ary_new3(3, rb_sestpr, rb_s, rb_c);
}

void
init_lapack_slaic1(VALUE mLapack){
  rb_define_module_function(mLapack, "slaic1", rb_slaic1, -1);
}
