#include "rb_lapack.h"

extern VOID slasd5_(integer *i, real *d, real *z, real *delta, real *rho, real *dsigma, real *work);

static VALUE
rb_slasd5(int argc, VALUE *argv, VALUE self){
  VALUE rb_i;
  integer i; 
  VALUE rb_d;
  real *d; 
  VALUE rb_z;
  real *z; 
  VALUE rb_rho;
  real rho; 
  VALUE rb_delta;
  real *delta; 
  VALUE rb_dsigma;
  real dsigma; 
  real *work;


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  delta, dsigma = NumRu::Lapack.slasd5( i, d, z, rho)\n    or\n  NumRu::Lapack.slasd5  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_i = argv[0];
  rb_d = argv[1];
  rb_z = argv[2];
  rb_rho = argv[3];

  rho = (real)NUM2DBL(rb_rho);
  if (!NA_IsNArray(rb_z))
    rb_raise(rb_eArgError, "z (3th argument) must be NArray");
  if (NA_RANK(rb_z) != 1)
    rb_raise(rb_eArgError, "rank of z (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_z) != (2))
    rb_raise(rb_eRuntimeError, "shape 0 of z must be %d", 2);
  if (NA_TYPE(rb_z) != NA_SFLOAT)
    rb_z = na_change_type(rb_z, NA_SFLOAT);
  z = NA_PTR_TYPE(rb_z, real*);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (2th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (2th argument) must be %d", 1);
  if (NA_SHAPE0(rb_d) != (2))
    rb_raise(rb_eRuntimeError, "shape 0 of d must be %d", 2);
  if (NA_TYPE(rb_d) != NA_SFLOAT)
    rb_d = na_change_type(rb_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rb_d, real*);
  i = NUM2INT(rb_i);
  {
    int shape[1];
    shape[0] = 2;
    rb_delta = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  delta = NA_PTR_TYPE(rb_delta, real*);
  work = ALLOC_N(real, (2));

  slasd5_(&i, d, z, delta, &rho, &dsigma, work);

  free(work);
  rb_dsigma = rb_float_new((double)dsigma);
  return rb_ary_new3(2, rb_delta, rb_dsigma);
}

void
init_lapack_slasd5(VALUE mLapack){
  rb_define_module_function(mLapack, "slasd5", rb_slasd5, -1);
}
