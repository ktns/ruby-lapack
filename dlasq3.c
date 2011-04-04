#include "rb_lapack.h"

extern VOID dlasq3_(integer *i0, integer *n0, doublereal *z, integer *pp, doublereal *dmin, doublereal *sigma, doublereal *desig, doublereal *qmax, integer *nfail, integer *iter, integer *ndiv, logical *ieee, integer *ttype, doublereal *dmin1, doublereal *dmin2, doublereal *dn, doublereal *dn1, doublereal *dn2, doublereal *g, doublereal *tau);

static VALUE
rb_dlasq3(int argc, VALUE *argv, VALUE self){
  VALUE rb_i0;
  integer i0; 
  VALUE rb_n0;
  integer n0; 
  VALUE rb_z;
  doublereal *z; 
  VALUE rb_pp;
  integer pp; 
  VALUE rb_desig;
  doublereal desig; 
  VALUE rb_qmax;
  doublereal qmax; 
  VALUE rb_ieee;
  logical ieee; 
  VALUE rb_ttype;
  integer ttype; 
  VALUE rb_dmin1;
  doublereal dmin1; 
  VALUE rb_dmin2;
  doublereal dmin2; 
  VALUE rb_dn;
  doublereal dn; 
  VALUE rb_dn1;
  doublereal dn1; 
  VALUE rb_dn2;
  doublereal dn2; 
  VALUE rb_g;
  doublereal g; 
  VALUE rb_tau;
  doublereal tau; 
  VALUE rb_dmin;
  doublereal dmin; 
  VALUE rb_sigma;
  doublereal sigma; 
  VALUE rb_nfail;
  integer nfail; 
  VALUE rb_iter;
  integer iter; 
  VALUE rb_ndiv;
  integer ndiv; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  dmin, sigma, nfail, iter, ndiv, n0, pp, desig, ttype, dmin1, dmin2, dn, dn1, dn2, g, tau = NumRu::Lapack.dlasq3( i0, n0, z, pp, desig, qmax, ieee, ttype, dmin1, dmin2, dn, dn1, dn2, g, tau)\n    or\n  NumRu::Lapack.dlasq3  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 15)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 15)", argc);
  rb_i0 = argv[0];
  rb_n0 = argv[1];
  rb_z = argv[2];
  rb_pp = argv[3];
  rb_desig = argv[4];
  rb_qmax = argv[5];
  rb_ieee = argv[6];
  rb_ttype = argv[7];
  rb_dmin1 = argv[8];
  rb_dmin2 = argv[9];
  rb_dn = argv[10];
  rb_dn1 = argv[11];
  rb_dn2 = argv[12];
  rb_g = argv[13];
  rb_tau = argv[14];

  pp = NUM2INT(rb_pp);
  n0 = NUM2INT(rb_n0);
  ttype = NUM2INT(rb_ttype);
  qmax = NUM2DBL(rb_qmax);
  dmin1 = NUM2DBL(rb_dmin1);
  desig = NUM2DBL(rb_desig);
  dmin2 = NUM2DBL(rb_dmin2);
  dn = NUM2DBL(rb_dn);
  dn1 = NUM2DBL(rb_dn1);
  i0 = NUM2INT(rb_i0);
  tau = NUM2DBL(rb_tau);
  dn2 = NUM2DBL(rb_dn2);
  ieee = (rb_ieee == Qtrue);
  g = NUM2DBL(rb_g);
  if (!NA_IsNArray(rb_z))
    rb_raise(rb_eArgError, "z (3th argument) must be NArray");
  if (NA_RANK(rb_z) != 1)
    rb_raise(rb_eArgError, "rank of z (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_z) != (4*n0))
    rb_raise(rb_eRuntimeError, "shape 0 of z must be %d", 4*n0);
  if (NA_TYPE(rb_z) != NA_DFLOAT)
    rb_z = na_change_type(rb_z, NA_DFLOAT);
  z = NA_PTR_TYPE(rb_z, doublereal*);

  dlasq3_(&i0, &n0, z, &pp, &dmin, &sigma, &desig, &qmax, &nfail, &iter, &ndiv, &ieee, &ttype, &dmin1, &dmin2, &dn, &dn1, &dn2, &g, &tau);

  rb_dmin = rb_float_new((double)dmin);
  rb_sigma = rb_float_new((double)sigma);
  rb_nfail = INT2NUM(nfail);
  rb_iter = INT2NUM(iter);
  rb_ndiv = INT2NUM(ndiv);
  rb_n0 = INT2NUM(n0);
  rb_pp = INT2NUM(pp);
  rb_desig = rb_float_new((double)desig);
  rb_ttype = INT2NUM(ttype);
  rb_dmin1 = rb_float_new((double)dmin1);
  rb_dmin2 = rb_float_new((double)dmin2);
  rb_dn = rb_float_new((double)dn);
  rb_dn1 = rb_float_new((double)dn1);
  rb_dn2 = rb_float_new((double)dn2);
  rb_g = rb_float_new((double)g);
  rb_tau = rb_float_new((double)tau);
  return rb_ary_new3(16, rb_dmin, rb_sigma, rb_nfail, rb_iter, rb_ndiv, rb_n0, rb_pp, rb_desig, rb_ttype, rb_dmin1, rb_dmin2, rb_dn, rb_dn1, rb_dn2, rb_g, rb_tau);
}

void
init_lapack_dlasq3(VALUE mLapack){
  rb_define_module_function(mLapack, "dlasq3", rb_dlasq3, -1);
}
