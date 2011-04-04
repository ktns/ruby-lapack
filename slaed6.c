#include "rb_lapack.h"

extern VOID slaed6_(integer *kniter, logical *orgati, real *rho, real *d, real *z, real *finit, real *tau, integer *info);

static VALUE
rb_slaed6(int argc, VALUE *argv, VALUE self){
  VALUE rb_kniter;
  integer kniter; 
  VALUE rb_orgati;
  logical orgati; 
  VALUE rb_rho;
  real rho; 
  VALUE rb_d;
  real *d; 
  VALUE rb_z;
  real *z; 
  VALUE rb_finit;
  real finit; 
  VALUE rb_tau;
  real tau; 
  VALUE rb_info;
  integer info; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  tau, info = NumRu::Lapack.slaed6( kniter, orgati, rho, d, z, finit)\n    or\n  NumRu::Lapack.slaed6  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rb_kniter = argv[0];
  rb_orgati = argv[1];
  rb_rho = argv[2];
  rb_d = argv[3];
  rb_z = argv[4];
  rb_finit = argv[5];

  orgati = (rb_orgati == Qtrue);
  finit = (real)NUM2DBL(rb_finit);
  rho = (real)NUM2DBL(rb_rho);
  if (!NA_IsNArray(rb_z))
    rb_raise(rb_eArgError, "z (5th argument) must be NArray");
  if (NA_RANK(rb_z) != 1)
    rb_raise(rb_eArgError, "rank of z (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_z) != (3))
    rb_raise(rb_eRuntimeError, "shape 0 of z must be %d", 3);
  if (NA_TYPE(rb_z) != NA_SFLOAT)
    rb_z = na_change_type(rb_z, NA_SFLOAT);
  z = NA_PTR_TYPE(rb_z, real*);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (4th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_d) != (3))
    rb_raise(rb_eRuntimeError, "shape 0 of d must be %d", 3);
  if (NA_TYPE(rb_d) != NA_SFLOAT)
    rb_d = na_change_type(rb_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rb_d, real*);
  kniter = NUM2INT(rb_kniter);

  slaed6_(&kniter, &orgati, &rho, d, z, &finit, &tau, &info);

  rb_tau = rb_float_new((double)tau);
  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_tau, rb_info);
}

void
init_lapack_slaed6(VALUE mLapack){
  rb_define_module_function(mLapack, "slaed6", rb_slaed6, -1);
}
