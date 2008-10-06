#include "rb_lapack.h"

static VALUE
rb_slasq3(int argc, VALUE *argv, VALUE self){
  VALUE rb_i0;
  integer i0; 
  VALUE rb_n0;
  integer n0; 
  VALUE rb_z;
  real *z; 
  VALUE rb_pp;
  integer pp; 
  VALUE rb_desig;
  real desig; 
  VALUE rb_qmax;
  real qmax; 
  VALUE rb_ieee;
  logical ieee; 
  VALUE rb_dmin;
  real dmin; 
  VALUE rb_sigma;
  real sigma; 
  VALUE rb_nfail;
  integer nfail; 
  VALUE rb_iter;
  integer iter; 
  VALUE rb_ndiv;
  integer ndiv; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  dmin, sigma, nfail, iter, ndiv, desig = NumRu::Lapack.slasq3( i0, n0, z, pp, desig, qmax, ieee)\n    or\n  NumRu::Lapack.slasq3  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLASQ3( I0, N0, Z, PP, DMIN, SIGMA, DESIG, QMAX, NFAIL, ITER, NDIV, IEEE )\n\n*  Purpose\n*  =======\n*\n*  SLASQ3 checks for deflation, computes a shift (TAU) and calls dqds.\n*  In case of failure it changes shifts, and tries again until output\n*  is positive.\n*\n\n*  Arguments\n*  =========\n*\n*  I0     (input) INTEGER\n*         First index.\n*\n*  N0     (input) INTEGER\n*         Last index.\n*\n*  Z      (input) REAL array, dimension ( 4*N )\n*         Z holds the qd array.\n*\n*  PP     (input) INTEGER\n*         PP=0 for ping, PP=1 for pong.\n*\n*  DMIN   (output) REAL\n*         Minimum value of d.\n*\n*  SIGMA  (output) REAL\n*         Sum of shifts used in current segment.\n*\n*  DESIG  (input/output) REAL\n*         Lower order part of SIGMA\n*\n*  QMAX   (input) REAL\n*         Maximum value of q.\n*\n*  NFAIL  (output) INTEGER\n*         Number of times shift was too big.\n*\n*  ITER   (output) INTEGER\n*         Number of iterations.\n*\n*  NDIV   (output) INTEGER\n*         Number of divisions.\n*\n*  TTYPE  (output) INTEGER\n*         Shift type.\n*\n*  IEEE   (input) LOGICAL\n*         Flag for IEEE or non IEEE arithmetic (passed to SLASQ5).\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 7)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 7)", argc);
  rb_i0 = argv[0];
  rb_n0 = argv[1];
  rb_z = argv[2];
  rb_pp = argv[3];
  rb_desig = argv[4];
  rb_qmax = argv[5];
  rb_ieee = argv[6];

  i0 = NUM2INT(rb_i0);
  n0 = NUM2INT(rb_n0);
  pp = NUM2INT(rb_pp);
  desig = (real)NUM2DBL(rb_desig);
  qmax = (real)NUM2DBL(rb_qmax);
  ieee = (rb_ieee == Qtrue);
  if (!NA_IsNArray(rb_z))
    rb_raise(rb_eArgError, "z (3th argument) must be NArray");
  if (NA_RANK(rb_z) != 1)
    rb_raise(rb_eArgError, "rank of z (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_z) != (4*n0))
    rb_raise(rb_eRuntimeError, "shape 0 of z must be %d", 4*n0);
  if (NA_TYPE(rb_z) != NA_SFLOAT)
    rb_z = na_change_type(rb_z, NA_SFLOAT);
  z = NA_PTR_TYPE(rb_z, real*);

  slasq3_(&i0, &n0, z, &pp, &dmin, &sigma, &desig, &qmax, &nfail, &iter, &ndiv, &ieee);

  rb_dmin = rb_float_new((double)dmin);
  rb_sigma = rb_float_new((double)sigma);
  rb_nfail = INT2NUM(nfail);
  rb_iter = INT2NUM(iter);
  rb_ndiv = INT2NUM(ndiv);
  rb_desig = rb_float_new((double)desig);
  return rb_ary_new3(6, rb_dmin, rb_sigma, rb_nfail, rb_iter, rb_ndiv, rb_desig);
}

void
init_lapack_slasq3(VALUE mLapack){
  rb_define_module_function(mLapack, "slasq3", rb_slasq3, -1);
}
