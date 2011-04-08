#include "rb_lapack.h"

extern VOID slarrk_(integer *n, integer *iw, real *gl, real *gu, real *d, real *e2, real *pivmin, real *reltol, real *w, real *werr, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_slarrk(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_iw;
  integer iw; 
  VALUE rblapack_gl;
  real gl; 
  VALUE rblapack_gu;
  real gu; 
  VALUE rblapack_d;
  real *d; 
  VALUE rblapack_e2;
  real *e2; 
  VALUE rblapack_pivmin;
  real pivmin; 
  VALUE rblapack_reltol;
  real reltol; 
  VALUE rblapack_w;
  real w; 
  VALUE rblapack_werr;
  real werr; 
  VALUE rblapack_info;
  integer info; 

  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  w, werr, info = NumRu::Lapack.slarrk( iw, gl, gu, d, e2, pivmin, reltol, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLARRK( N, IW, GL, GU, D, E2, PIVMIN, RELTOL, W, WERR, INFO)\n\n*  Purpose\n*  =======\n*\n*  SLARRK computes one eigenvalue of a symmetric tridiagonal\n*  matrix T to suitable accuracy. This is an auxiliary code to be\n*  called from SSTEMR.\n*\n*  To avoid overflow, the matrix must be scaled so that its\n*  largest element is no greater than overflow**(1/2) *\n*  underflow**(1/4) in absolute value, and for greatest\n*  accuracy, it should not be much smaller than that.\n*\n*  See W. Kahan \"Accurate Eigenvalues of a Symmetric Tridiagonal\n*  Matrix\", Report CS41, Computer Science Dept., Stanford\n*  University, July 21, 1966.\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The order of the tridiagonal matrix T.  N >= 0.\n*\n*  IW      (input) INTEGER\n*          The index of the eigenvalues to be returned.\n*\n*  GL      (input) REAL            \n*  GU      (input) REAL            \n*          An upper and a lower bound on the eigenvalue.\n*\n*  D       (input) REAL             array, dimension (N)\n*          The n diagonal elements of the tridiagonal matrix T.\n*\n*  E2      (input) REAL             array, dimension (N-1)\n*          The (n-1) squared off-diagonal elements of the tridiagonal matrix T.\n*\n*  PIVMIN  (input) REAL            \n*          The minimum pivot allowed in the Sturm sequence for T.\n*\n*  RELTOL  (input) REAL            \n*          The minimum relative width of an interval.  When an interval\n*          is narrower than RELTOL times the larger (in\n*          magnitude) endpoint, then it is considered to be\n*          sufficiently small, i.e., converged.  Note: this should\n*          always be at least radix*machine epsilon.\n*\n*  W       (output) REAL            \n*\n*  WERR    (output) REAL            \n*          The error bound on the corresponding eigenvalue approximation\n*          in W.\n*\n*  INFO    (output) INTEGER\n*          = 0:       Eigenvalue converged\n*          = -1:      Eigenvalue did NOT converge\n*\n*  Internal Parameters\n*  ===================\n*\n*  FUDGE   REAL            , default = 2\n*          A \"fudge factor\" to widen the Gershgorin intervals.\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  w, werr, info = NumRu::Lapack.slarrk( iw, gl, gu, d, e2, pivmin, reltol, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 7)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 7)", argc);
  rblapack_iw = argv[0];
  rblapack_gl = argv[1];
  rblapack_gu = argv[2];
  rblapack_d = argv[3];
  rblapack_e2 = argv[4];
  rblapack_pivmin = argv[5];
  rblapack_reltol = argv[6];
  if (rb_options != Qnil) {
  }

  pivmin = (real)NUM2DBL(rblapack_pivmin);
  gu = (real)NUM2DBL(rblapack_gu);
  iw = NUM2INT(rblapack_iw);
  gl = (real)NUM2DBL(rblapack_gl);
  if (!NA_IsNArray(rblapack_d))
    rb_raise(rb_eArgError, "d (4th argument) must be NArray");
  if (NA_RANK(rblapack_d) != 1)
    rb_raise(rb_eArgError, "rank of d (4th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_d);
  if (NA_TYPE(rblapack_d) != NA_SFLOAT)
    rblapack_d = na_change_type(rblapack_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rblapack_d, real*);
  reltol = (real)NUM2DBL(rblapack_reltol);
  if (!NA_IsNArray(rblapack_e2))
    rb_raise(rb_eArgError, "e2 (5th argument) must be NArray");
  if (NA_RANK(rblapack_e2) != 1)
    rb_raise(rb_eArgError, "rank of e2 (5th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_e2) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of e2 must be %d", n-1);
  if (NA_TYPE(rblapack_e2) != NA_SFLOAT)
    rblapack_e2 = na_change_type(rblapack_e2, NA_SFLOAT);
  e2 = NA_PTR_TYPE(rblapack_e2, real*);

  slarrk_(&n, &iw, &gl, &gu, d, e2, &pivmin, &reltol, &w, &werr, &info);

  rblapack_w = rb_float_new((double)w);
  rblapack_werr = rb_float_new((double)werr);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(3, rblapack_w, rblapack_werr, rblapack_info);
}

void
init_lapack_slarrk(VALUE mLapack){
  rb_define_module_function(mLapack, "slarrk", rblapack_slarrk, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
