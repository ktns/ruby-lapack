#include "rb_lapack.h"

extern VOID slasq4_(integer *i0, integer *n0, real *z, integer *pp, integer *n0in, real *dmin, real *dmin1, real *dmin2, real *dn, real *dn1, real *dn2, real *tau, integer *ttype, real *g);

static VALUE
rb_slasq4(int argc, VALUE *argv, VALUE self){
  VALUE rb_i0;
  integer i0; 
  VALUE rb_n0;
  integer n0; 
  VALUE rb_z;
  real *z; 
  VALUE rb_pp;
  integer pp; 
  VALUE rb_n0in;
  integer n0in; 
  VALUE rb_dmin;
  real dmin; 
  VALUE rb_dmin1;
  real dmin1; 
  VALUE rb_dmin2;
  real dmin2; 
  VALUE rb_dn;
  real dn; 
  VALUE rb_dn1;
  real dn1; 
  VALUE rb_dn2;
  real dn2; 
  VALUE rb_g;
  real g; 
  VALUE rb_tau;
  real tau; 
  VALUE rb_ttype;
  integer ttype; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  tau, ttype, g = NumRu::Lapack.slasq4( i0, n0, z, pp, n0in, dmin, dmin1, dmin2, dn, dn1, dn2, g)\n    or\n  NumRu::Lapack.slasq4  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLASQ4( I0, N0, Z, PP, N0IN, DMIN, DMIN1, DMIN2, DN, DN1, DN2, TAU, TTYPE, G )\n\n*  Purpose\n*  =======\n*\n*  SLASQ4 computes an approximation TAU to the smallest eigenvalue\n*  using values of d from the previous transform.\n*\n\n*  I0    (input) INTEGER\n*        First index.\n*\n*  N0    (input) INTEGER\n*        Last index.\n*\n*  Z     (input) REAL array, dimension ( 4*N )\n*        Z holds the qd array.\n*\n*  PP    (input) INTEGER\n*        PP=0 for ping, PP=1 for pong.\n*\n*  NOIN  (input) INTEGER\n*        The value of N0 at start of EIGTEST.\n*\n*  DMIN  (input) REAL\n*        Minimum value of d.\n*\n*  DMIN1 (input) REAL\n*        Minimum value of d, excluding D( N0 ).\n*\n*  DMIN2 (input) REAL\n*        Minimum value of d, excluding D( N0 ) and D( N0-1 ).\n*\n*  DN    (input) REAL\n*        d(N)\n*\n*  DN1   (input) REAL\n*        d(N-1)\n*\n*  DN2   (input) REAL\n*        d(N-2)\n*\n*  TAU   (output) REAL\n*        This is the shift.\n*\n*  TTYPE (output) INTEGER\n*        Shift type.\n*\n*  G     (input/output) REAL\n*        G is passed as an argument in order to save its value between\n*        calls to SLASQ4.\n*\n\n*  Further Details\n*  ===============\n*  CNST1 = 9/16\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 12)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 12)", argc);
  rb_i0 = argv[0];
  rb_n0 = argv[1];
  rb_z = argv[2];
  rb_pp = argv[3];
  rb_n0in = argv[4];
  rb_dmin = argv[5];
  rb_dmin1 = argv[6];
  rb_dmin2 = argv[7];
  rb_dn = argv[8];
  rb_dn1 = argv[9];
  rb_dn2 = argv[10];
  rb_g = argv[11];

  pp = NUM2INT(rb_pp);
  n0 = NUM2INT(rb_n0);
  dn = (real)NUM2DBL(rb_dn);
  dmin1 = (real)NUM2DBL(rb_dmin1);
  dmin = (real)NUM2DBL(rb_dmin);
  dmin2 = (real)NUM2DBL(rb_dmin2);
  dn2 = (real)NUM2DBL(rb_dn2);
  dn1 = (real)NUM2DBL(rb_dn1);
  n0in = NUM2INT(rb_n0in);
  i0 = NUM2INT(rb_i0);
  g = (real)NUM2DBL(rb_g);
  if (!NA_IsNArray(rb_z))
    rb_raise(rb_eArgError, "z (3th argument) must be NArray");
  if (NA_RANK(rb_z) != 1)
    rb_raise(rb_eArgError, "rank of z (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_z) != (4*n0))
    rb_raise(rb_eRuntimeError, "shape 0 of z must be %d", 4*n0);
  if (NA_TYPE(rb_z) != NA_SFLOAT)
    rb_z = na_change_type(rb_z, NA_SFLOAT);
  z = NA_PTR_TYPE(rb_z, real*);

  slasq4_(&i0, &n0, z, &pp, &n0in, &dmin, &dmin1, &dmin2, &dn, &dn1, &dn2, &tau, &ttype, &g);

  rb_tau = rb_float_new((double)tau);
  rb_ttype = INT2NUM(ttype);
  rb_g = rb_float_new((double)g);
  return rb_ary_new3(3, rb_tau, rb_ttype, rb_g);
}

void
init_lapack_slasq4(VALUE mLapack){
  rb_define_module_function(mLapack, "slasq4", rb_slasq4, -1);
}
