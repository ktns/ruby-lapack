#include "rb_lapack.h"

extern VOID slasq4_(integer *i0, integer *n0, real *z, integer *pp, integer *n0in, real *dmin, real *dmin1, real *dmin2, real *dn, real *dn1, real *dn2, real *tau, integer *ttype, real *g);

static VALUE sHelp, sUsage;

static VALUE
rblapack_slasq4(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_i0;
  integer i0; 
  VALUE rblapack_n0;
  integer n0; 
  VALUE rblapack_z;
  real *z; 
  VALUE rblapack_pp;
  integer pp; 
  VALUE rblapack_n0in;
  integer n0in; 
  VALUE rblapack_dmin;
  real dmin; 
  VALUE rblapack_dmin1;
  real dmin1; 
  VALUE rblapack_dmin2;
  real dmin2; 
  VALUE rblapack_dn;
  real dn; 
  VALUE rblapack_dn1;
  real dn1; 
  VALUE rblapack_dn2;
  real dn2; 
  VALUE rblapack_g;
  real g; 
  VALUE rblapack_tau;
  real tau; 
  VALUE rblapack_ttype;
  integer ttype; 


  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  tau, ttype, g = NumRu::Lapack.slasq4( i0, n0, z, pp, n0in, dmin, dmin1, dmin2, dn, dn1, dn2, g, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLASQ4( I0, N0, Z, PP, N0IN, DMIN, DMIN1, DMIN2, DN, DN1, DN2, TAU, TTYPE, G )\n\n*  Purpose\n*  =======\n*\n*  SLASQ4 computes an approximation TAU to the smallest eigenvalue\n*  using values of d from the previous transform.\n*\n\n*  I0    (input) INTEGER\n*        First index.\n*\n*  N0    (input) INTEGER\n*        Last index.\n*\n*  Z     (input) REAL array, dimension ( 4*N )\n*        Z holds the qd array.\n*\n*  PP    (input) INTEGER\n*        PP=0 for ping, PP=1 for pong.\n*\n*  NOIN  (input) INTEGER\n*        The value of N0 at start of EIGTEST.\n*\n*  DMIN  (input) REAL\n*        Minimum value of d.\n*\n*  DMIN1 (input) REAL\n*        Minimum value of d, excluding D( N0 ).\n*\n*  DMIN2 (input) REAL\n*        Minimum value of d, excluding D( N0 ) and D( N0-1 ).\n*\n*  DN    (input) REAL\n*        d(N)\n*\n*  DN1   (input) REAL\n*        d(N-1)\n*\n*  DN2   (input) REAL\n*        d(N-2)\n*\n*  TAU   (output) REAL\n*        This is the shift.\n*\n*  TTYPE (output) INTEGER\n*        Shift type.\n*\n*  G     (input/output) REAL\n*        G is passed as an argument in order to save its value between\n*        calls to SLASQ4.\n*\n\n*  Further Details\n*  ===============\n*  CNST1 = 9/16\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  tau, ttype, g = NumRu::Lapack.slasq4( i0, n0, z, pp, n0in, dmin, dmin1, dmin2, dn, dn1, dn2, g, [:usage => usage, :help => help])\n");
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
  dn = (real)NUM2DBL(rblapack_dn);
  dmin1 = (real)NUM2DBL(rblapack_dmin1);
  dmin = (real)NUM2DBL(rblapack_dmin);
  dmin2 = (real)NUM2DBL(rblapack_dmin2);
  dn2 = (real)NUM2DBL(rblapack_dn2);
  dn1 = (real)NUM2DBL(rblapack_dn1);
  n0in = NUM2INT(rblapack_n0in);
  i0 = NUM2INT(rblapack_i0);
  g = (real)NUM2DBL(rblapack_g);
  if (!NA_IsNArray(rblapack_z))
    rb_raise(rb_eArgError, "z (3th argument) must be NArray");
  if (NA_RANK(rblapack_z) != 1)
    rb_raise(rb_eArgError, "rank of z (3th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_z) != (4*n0))
    rb_raise(rb_eRuntimeError, "shape 0 of z must be %d", 4*n0);
  if (NA_TYPE(rblapack_z) != NA_SFLOAT)
    rblapack_z = na_change_type(rblapack_z, NA_SFLOAT);
  z = NA_PTR_TYPE(rblapack_z, real*);

  slasq4_(&i0, &n0, z, &pp, &n0in, &dmin, &dmin1, &dmin2, &dn, &dn1, &dn2, &tau, &ttype, &g);

  rblapack_tau = rb_float_new((double)tau);
  rblapack_ttype = INT2NUM(ttype);
  rblapack_g = rb_float_new((double)g);
  return rb_ary_new3(3, rblapack_tau, rblapack_ttype, rblapack_g);
}

void
init_lapack_slasq4(VALUE mLapack){
  rb_define_module_function(mLapack, "slasq4", rblapack_slasq4, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
