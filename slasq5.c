#include "rb_lapack.h"

extern VOID slasq5_(integer *i0, integer *n0, real *z, integer *pp, real *tau, real *dmin, real *dmin1, real *dmin2, real *dn, real *dnm1, real *dnm2, logical *ieee);

static VALUE sHelp, sUsage;

static VALUE
rblapack_slasq5(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_i0;
  integer i0; 
  VALUE rblapack_n0;
  integer n0; 
  VALUE rblapack_z;
  real *z; 
  VALUE rblapack_pp;
  integer pp; 
  VALUE rblapack_tau;
  real tau; 
  VALUE rblapack_ieee;
  logical ieee; 
  VALUE rblapack_dmin;
  real dmin; 
  VALUE rblapack_dmin1;
  real dmin1; 
  VALUE rblapack_dmin2;
  real dmin2; 
  VALUE rblapack_dn;
  real dn; 
  VALUE rblapack_dnm1;
  real dnm1; 
  VALUE rblapack_dnm2;
  real dnm2; 


  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  dmin, dmin1, dmin2, dn, dnm1, dnm2 = NumRu::Lapack.slasq5( i0, n0, z, pp, tau, ieee, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLASQ5( I0, N0, Z, PP, TAU, DMIN, DMIN1, DMIN2, DN, DNM1, DNM2, IEEE )\n\n*  Purpose\n*  =======\n*\n*  SLASQ5 computes one dqds transform in ping-pong form, one\n*  version for IEEE machines another for non IEEE machines.\n*\n\n*  Arguments\n*  =========\n*\n*  I0    (input) INTEGER\n*        First index.\n*\n*  N0    (input) INTEGER\n*        Last index.\n*\n*  Z     (input) REAL array, dimension ( 4*N )\n*        Z holds the qd array. EMIN is stored in Z(4*N0) to avoid\n*        an extra argument.\n*\n*  PP    (input) INTEGER\n*        PP=0 for ping, PP=1 for pong.\n*\n*  TAU   (input) REAL\n*        This is the shift.\n*\n*  DMIN  (output) REAL\n*        Minimum value of d.\n*\n*  DMIN1 (output) REAL\n*        Minimum value of d, excluding D( N0 ).\n*\n*  DMIN2 (output) REAL\n*        Minimum value of d, excluding D( N0 ) and D( N0-1 ).\n*\n*  DN    (output) REAL\n*        d(N0), the last value of d.\n*\n*  DNM1  (output) REAL\n*        d(N0-1).\n*\n*  DNM2  (output) REAL\n*        d(N0-2).\n*\n*  IEEE  (input) LOGICAL\n*        Flag for IEEE or non IEEE arithmetic.\n*\n\n*  =====================================================================\n*\n*     .. Parameter ..\n      REAL               ZERO\n      PARAMETER          ( ZERO = 0.0E0 )\n*     ..\n*     .. Local Scalars ..\n      INTEGER            J4, J4P2\n      REAL               D, EMIN, TEMP\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          MIN\n*     ..\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  dmin, dmin1, dmin2, dn, dnm1, dnm2 = NumRu::Lapack.slasq5( i0, n0, z, pp, tau, ieee, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rblapack_i0 = argv[0];
  rblapack_n0 = argv[1];
  rblapack_z = argv[2];
  rblapack_pp = argv[3];
  rblapack_tau = argv[4];
  rblapack_ieee = argv[5];
  if (rb_options != Qnil) {
  }

  pp = NUM2INT(rblapack_pp);
  n0 = NUM2INT(rblapack_n0);
  tau = (real)NUM2DBL(rblapack_tau);
  ieee = (rblapack_ieee == Qtrue);
  i0 = NUM2INT(rblapack_i0);
  if (!NA_IsNArray(rblapack_z))
    rb_raise(rb_eArgError, "z (3th argument) must be NArray");
  if (NA_RANK(rblapack_z) != 1)
    rb_raise(rb_eArgError, "rank of z (3th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_z) != (4*n0))
    rb_raise(rb_eRuntimeError, "shape 0 of z must be %d", 4*n0);
  if (NA_TYPE(rblapack_z) != NA_SFLOAT)
    rblapack_z = na_change_type(rblapack_z, NA_SFLOAT);
  z = NA_PTR_TYPE(rblapack_z, real*);

  slasq5_(&i0, &n0, z, &pp, &tau, &dmin, &dmin1, &dmin2, &dn, &dnm1, &dnm2, &ieee);

  rblapack_dmin = rb_float_new((double)dmin);
  rblapack_dmin1 = rb_float_new((double)dmin1);
  rblapack_dmin2 = rb_float_new((double)dmin2);
  rblapack_dn = rb_float_new((double)dn);
  rblapack_dnm1 = rb_float_new((double)dnm1);
  rblapack_dnm2 = rb_float_new((double)dnm2);
  return rb_ary_new3(6, rblapack_dmin, rblapack_dmin1, rblapack_dmin2, rblapack_dn, rblapack_dnm1, rblapack_dnm2);
}

void
init_lapack_slasq5(VALUE mLapack){
  rb_define_module_function(mLapack, "slasq5", rblapack_slasq5, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
