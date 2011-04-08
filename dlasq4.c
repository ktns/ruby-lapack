#include "rb_lapack.h"

extern VOID dlasq4_(integer *i0, integer *n0, doublereal *z, integer *pp, integer *n0in, doublereal *dmin, doublereal *dmin1, doublereal *dmin2, doublereal *dn, doublereal *dn1, doublereal *dn2, doublereal *tau, integer *ttype, real *g);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dlasq4(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_i0;
  integer i0; 
  VALUE rblapack_n0;
  integer n0; 
  VALUE rblapack_z;
  doublereal *z; 
  VALUE rblapack_pp;
  integer pp; 
  VALUE rblapack_n0in;
  integer n0in; 
  VALUE rblapack_dmin;
  doublereal dmin; 
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
  real g; 
  VALUE rblapack_tau;
  doublereal tau; 
  VALUE rblapack_ttype;
  integer ttype; 


  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  tau, ttype, g = NumRu::Lapack.dlasq4( i0, n0, z, pp, n0in, dmin, dmin1, dmin2, dn, dn1, dn2, g, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLASQ4( I0, N0, Z, PP, N0IN, DMIN, DMIN1, DMIN2, DN, DN1, DN2, TAU, TTYPE, G )\n\n*  Purpose\n*  =======\n*\n*  DLASQ4 computes an approximation TAU to the smallest eigenvalue\n*  using values of d from the previous transform.\n*\n\n*  I0    (input) INTEGER\n*        First index.\n*\n*  N0    (input) INTEGER\n*        Last index.\n*\n*  Z     (input) DOUBLE PRECISION array, dimension ( 4*N )\n*        Z holds the qd array.\n*\n*  PP    (input) INTEGER\n*        PP=0 for ping, PP=1 for pong.\n*\n*  NOIN  (input) INTEGER\n*        The value of N0 at start of EIGTEST.\n*\n*  DMIN  (input) DOUBLE PRECISION\n*        Minimum value of d.\n*\n*  DMIN1 (input) DOUBLE PRECISION\n*        Minimum value of d, excluding D( N0 ).\n*\n*  DMIN2 (input) DOUBLE PRECISION\n*        Minimum value of d, excluding D( N0 ) and D( N0-1 ).\n*\n*  DN    (input) DOUBLE PRECISION\n*        d(N)\n*\n*  DN1   (input) DOUBLE PRECISION\n*        d(N-1)\n*\n*  DN2   (input) DOUBLE PRECISION\n*        d(N-2)\n*\n*  TAU   (output) DOUBLE PRECISION\n*        This is the shift.\n*\n*  TTYPE (output) INTEGER\n*        Shift type.\n*\n*  G     (input/output) REAL\n*        G is passed as an argument in order to save its value between\n*        calls to DLASQ4.\n*\n\n*  Further Details\n*  ===============\n*  CNST1 = 9/16\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  tau, ttype, g = NumRu::Lapack.dlasq4( i0, n0, z, pp, n0in, dmin, dmin1, dmin2, dn, dn1, dn2, g, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 12)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 12)", argc);
  rblapack_i0 = argv[0];
  rblapack_n0 = argv[1];
  rblapack_z = argv[2];
  rblapack_pp = argv[3];
  rblapack_n0in = argv[4];
  rblapack_dmin = argv[5];
  rblapack_dmin1 = argv[6];
  rblapack_dmin2 = argv[7];
  rblapack_dn = argv[8];
  rblapack_dn1 = argv[9];
  rblapack_dn2 = argv[10];
  rblapack_g = argv[11];
  if (rb_options != Qnil) {
  }

  pp = NUM2INT(rblapack_pp);
  n0 = NUM2INT(rblapack_n0);
  dn = NUM2DBL(rblapack_dn);
  dmin1 = NUM2DBL(rblapack_dmin1);
  dmin = NUM2DBL(rblapack_dmin);
  dmin2 = NUM2DBL(rblapack_dmin2);
  dn2 = NUM2DBL(rblapack_dn2);
  dn1 = NUM2DBL(rblapack_dn1);
  n0in = NUM2INT(rblapack_n0in);
  i0 = NUM2INT(rblapack_i0);
  g = (real)NUM2DBL(rblapack_g);
  if (!NA_IsNArray(rblapack_z))
    rb_raise(rb_eArgError, "z (3th argument) must be NArray");
  if (NA_RANK(rblapack_z) != 1)
    rb_raise(rb_eArgError, "rank of z (3th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_z) != (4*n0))
    rb_raise(rb_eRuntimeError, "shape 0 of z must be %d", 4*n0);
  if (NA_TYPE(rblapack_z) != NA_DFLOAT)
    rblapack_z = na_change_type(rblapack_z, NA_DFLOAT);
  z = NA_PTR_TYPE(rblapack_z, doublereal*);

  dlasq4_(&i0, &n0, z, &pp, &n0in, &dmin, &dmin1, &dmin2, &dn, &dn1, &dn2, &tau, &ttype, &g);

  rblapack_tau = rb_float_new((double)tau);
  rblapack_ttype = INT2NUM(ttype);
  rblapack_g = rb_float_new((double)g);
  return rb_ary_new3(3, rblapack_tau, rblapack_ttype, rblapack_g);
}

void
init_lapack_dlasq4(VALUE mLapack){
  rb_define_module_function(mLapack, "dlasq4", rblapack_dlasq4, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
