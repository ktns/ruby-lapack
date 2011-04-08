#include "rb_lapack.h"

extern VOID dlasq3_(integer *i0, integer *n0, doublereal *z, integer *pp, doublereal *dmin, doublereal *sigma, doublereal *desig, doublereal *qmax, integer *nfail, integer *iter, integer *ndiv, logical *ieee, integer *ttype, doublereal *dmin1, doublereal *dmin2, doublereal *dn, doublereal *dn1, doublereal *dn2, doublereal *g, doublereal *tau);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dlasq3(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_i0;
  integer i0; 
  VALUE rblapack_n0;
  integer n0; 
  VALUE rblapack_z;
  doublereal *z; 
  VALUE rblapack_pp;
  integer pp; 
  VALUE rblapack_desig;
  doublereal desig; 
  VALUE rblapack_qmax;
  doublereal qmax; 
  VALUE rblapack_ieee;
  logical ieee; 
  VALUE rblapack_ttype;
  integer ttype; 
  VALUE rblapack_dmin1;
  doublereal dmin1; 
  VALUE rblapack_dmin2;
  doublereal dmin2; 
  VALUE rblapack_dn;
  doublereal dn; 
  VALUE rblapack_dn1;
  doublereal dn1; 
  VALUE rblapack_dn2;
  doublereal dn2; 
  VALUE rblapack_g;
  doublereal g; 
  VALUE rblapack_tau;
  doublereal tau; 
  VALUE rblapack_dmin;
  doublereal dmin; 
  VALUE rblapack_sigma;
  doublereal sigma; 
  VALUE rblapack_nfail;
  integer nfail; 
  VALUE rblapack_iter;
  integer iter; 
  VALUE rblapack_ndiv;
  integer ndiv; 


  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  dmin, sigma, nfail, iter, ndiv, n0, pp, desig, ttype, dmin1, dmin2, dn, dn1, dn2, g, tau = NumRu::Lapack.dlasq3( i0, n0, z, pp, desig, qmax, ieee, ttype, dmin1, dmin2, dn, dn1, dn2, g, tau, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLASQ3( I0, N0, Z, PP, DMIN, SIGMA, DESIG, QMAX, NFAIL, ITER, NDIV, IEEE, TTYPE, DMIN1, DMIN2, DN, DN1, DN2, G, TAU )\n\n*  Purpose\n*  =======\n*\n*  DLASQ3 checks for deflation, computes a shift (TAU) and calls dqds.\n*  In case of failure it changes shifts, and tries again until output\n*  is positive.\n*\n\n*  Arguments\n*  =========\n*\n*  I0     (input) INTEGER\n*         First index.\n*\n*  N0     (input/output) INTEGER\n*         Last index.\n*\n*  Z      (input) DOUBLE PRECISION array, dimension ( 4*N )\n*         Z holds the qd array.\n*\n*  PP     (input/output) INTEGER\n*         PP=0 for ping, PP=1 for pong.\n*         PP=2 indicates that flipping was applied to the Z array   \n*         and that the initial tests for deflation should not be \n*         performed.\n*\n*  DMIN   (output) DOUBLE PRECISION\n*         Minimum value of d.\n*\n*  SIGMA  (output) DOUBLE PRECISION\n*         Sum of shifts used in current segment.\n*\n*  DESIG  (input/output) DOUBLE PRECISION\n*         Lower order part of SIGMA\n*\n*  QMAX   (input) DOUBLE PRECISION\n*         Maximum value of q.\n*\n*  NFAIL  (output) INTEGER\n*         Number of times shift was too big.\n*\n*  ITER   (output) INTEGER\n*         Number of iterations.\n*\n*  NDIV   (output) INTEGER\n*         Number of divisions.\n*\n*  IEEE   (input) LOGICAL\n*         Flag for IEEE or non IEEE arithmetic (passed to DLASQ5).\n*\n*  TTYPE  (input/output) INTEGER\n*         Shift type.\n*\n*  DMIN1  (input/output) DOUBLE PRECISION\n*\n*  DMIN2  (input/output) DOUBLE PRECISION\n*\n*  DN     (input/output) DOUBLE PRECISION\n*\n*  DN1    (input/output) DOUBLE PRECISION\n*\n*  DN2    (input/output) DOUBLE PRECISION\n*\n*  G      (input/output) DOUBLE PRECISION\n*\n*  TAU    (input/output) DOUBLE PRECISION\n*\n*         These are passed as arguments in order to save their values\n*         between calls to DLASQ3.\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  dmin, sigma, nfail, iter, ndiv, n0, pp, desig, ttype, dmin1, dmin2, dn, dn1, dn2, g, tau = NumRu::Lapack.dlasq3( i0, n0, z, pp, desig, qmax, ieee, ttype, dmin1, dmin2, dn, dn1, dn2, g, tau, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 15)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 15)", argc);
  rblapack_i0 = argv[0];
  rblapack_n0 = argv[1];
  rblapack_z = argv[2];
  rblapack_pp = argv[3];
  rblapack_desig = argv[4];
  rblapack_qmax = argv[5];
  rblapack_ieee = argv[6];
  rblapack_ttype = argv[7];
  rblapack_dmin1 = argv[8];
  rblapack_dmin2 = argv[9];
  rblapack_dn = argv[10];
  rblapack_dn1 = argv[11];
  rblapack_dn2 = argv[12];
  rblapack_g = argv[13];
  rblapack_tau = argv[14];
  if (rb_options != Qnil) {
  }

  pp = NUM2INT(rblapack_pp);
  n0 = NUM2INT(rblapack_n0);
  ttype = NUM2INT(rblapack_ttype);
  qmax = NUM2DBL(rblapack_qmax);
  dmin1 = NUM2DBL(rblapack_dmin1);
  desig = NUM2DBL(rblapack_desig);
  dmin2 = NUM2DBL(rblapack_dmin2);
  dn = NUM2DBL(rblapack_dn);
  dn1 = NUM2DBL(rblapack_dn1);
  i0 = NUM2INT(rblapack_i0);
  tau = NUM2DBL(rblapack_tau);
  dn2 = NUM2DBL(rblapack_dn2);
  ieee = (rblapack_ieee == Qtrue);
  g = NUM2DBL(rblapack_g);
  if (!NA_IsNArray(rblapack_z))
    rb_raise(rb_eArgError, "z (3th argument) must be NArray");
  if (NA_RANK(rblapack_z) != 1)
    rb_raise(rb_eArgError, "rank of z (3th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_z) != (4*n0))
    rb_raise(rb_eRuntimeError, "shape 0 of z must be %d", 4*n0);
  if (NA_TYPE(rblapack_z) != NA_DFLOAT)
    rblapack_z = na_change_type(rblapack_z, NA_DFLOAT);
  z = NA_PTR_TYPE(rblapack_z, doublereal*);

  dlasq3_(&i0, &n0, z, &pp, &dmin, &sigma, &desig, &qmax, &nfail, &iter, &ndiv, &ieee, &ttype, &dmin1, &dmin2, &dn, &dn1, &dn2, &g, &tau);

  rblapack_dmin = rb_float_new((double)dmin);
  rblapack_sigma = rb_float_new((double)sigma);
  rblapack_nfail = INT2NUM(nfail);
  rblapack_iter = INT2NUM(iter);
  rblapack_ndiv = INT2NUM(ndiv);
  rblapack_n0 = INT2NUM(n0);
  rblapack_pp = INT2NUM(pp);
  rblapack_desig = rb_float_new((double)desig);
  rblapack_ttype = INT2NUM(ttype);
  rblapack_dmin1 = rb_float_new((double)dmin1);
  rblapack_dmin2 = rb_float_new((double)dmin2);
  rblapack_dn = rb_float_new((double)dn);
  rblapack_dn1 = rb_float_new((double)dn1);
  rblapack_dn2 = rb_float_new((double)dn2);
  rblapack_g = rb_float_new((double)g);
  rblapack_tau = rb_float_new((double)tau);
  return rb_ary_new3(16, rblapack_dmin, rblapack_sigma, rblapack_nfail, rblapack_iter, rblapack_ndiv, rblapack_n0, rblapack_pp, rblapack_desig, rblapack_ttype, rblapack_dmin1, rblapack_dmin2, rblapack_dn, rblapack_dn1, rblapack_dn2, rblapack_g, rblapack_tau);
}

void
init_lapack_dlasq3(VALUE mLapack){
  rb_define_module_function(mLapack, "dlasq3", rblapack_dlasq3, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
