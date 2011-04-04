#include "rb_lapack.h"

extern VOID slasq5_(integer *i0, integer *n0, real *z, integer *pp, real *tau, real *dmin, real *dmin1, real *dmin2, real *dn, real *dnm1, real *dnm2, logical *ieee);

static VALUE
rb_slasq5(int argc, VALUE *argv, VALUE self){
  VALUE rb_i0;
  integer i0; 
  VALUE rb_n0;
  integer n0; 
  VALUE rb_z;
  real *z; 
  VALUE rb_pp;
  integer pp; 
  VALUE rb_tau;
  real tau; 
  VALUE rb_ieee;
  logical ieee; 
  VALUE rb_dmin;
  real dmin; 
  VALUE rb_dmin1;
  real dmin1; 
  VALUE rb_dmin2;
  real dmin2; 
  VALUE rb_dn;
  real dn; 
  VALUE rb_dnm1;
  real dnm1; 
  VALUE rb_dnm2;
  real dnm2; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  dmin, dmin1, dmin2, dn, dnm1, dnm2 = NumRu::Lapack.slasq5( i0, n0, z, pp, tau, ieee)\n    or\n  NumRu::Lapack.slasq5  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rb_i0 = argv[0];
  rb_n0 = argv[1];
  rb_z = argv[2];
  rb_pp = argv[3];
  rb_tau = argv[4];
  rb_ieee = argv[5];

  pp = NUM2INT(rb_pp);
  n0 = NUM2INT(rb_n0);
  tau = (real)NUM2DBL(rb_tau);
  ieee = (rb_ieee == Qtrue);
  i0 = NUM2INT(rb_i0);
  if (!NA_IsNArray(rb_z))
    rb_raise(rb_eArgError, "z (3th argument) must be NArray");
  if (NA_RANK(rb_z) != 1)
    rb_raise(rb_eArgError, "rank of z (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_z) != (4*n0))
    rb_raise(rb_eRuntimeError, "shape 0 of z must be %d", 4*n0);
  if (NA_TYPE(rb_z) != NA_SFLOAT)
    rb_z = na_change_type(rb_z, NA_SFLOAT);
  z = NA_PTR_TYPE(rb_z, real*);

  slasq5_(&i0, &n0, z, &pp, &tau, &dmin, &dmin1, &dmin2, &dn, &dnm1, &dnm2, &ieee);

  rb_dmin = rb_float_new((double)dmin);
  rb_dmin1 = rb_float_new((double)dmin1);
  rb_dmin2 = rb_float_new((double)dmin2);
  rb_dn = rb_float_new((double)dn);
  rb_dnm1 = rb_float_new((double)dnm1);
  rb_dnm2 = rb_float_new((double)dnm2);
  return rb_ary_new3(6, rb_dmin, rb_dmin1, rb_dmin2, rb_dn, rb_dnm1, rb_dnm2);
}

void
init_lapack_slasq5(VALUE mLapack){
  rb_define_module_function(mLapack, "slasq5", rb_slasq5, -1);
}
