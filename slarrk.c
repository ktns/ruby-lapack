#include "rb_lapack.h"

static VALUE
rb_slarrk(int argc, VALUE *argv, VALUE self){
  VALUE rb_iw;
  integer iw; 
  VALUE rb_gl;
  real gl; 
  VALUE rb_gu;
  real gu; 
  VALUE rb_d;
  real *d; 
  VALUE rb_e2;
  real *e2; 
  VALUE rb_pivmin;
  real pivmin; 
  VALUE rb_reltol;
  real reltol; 
  VALUE rb_w;
  real w; 
  VALUE rb_werr;
  real werr; 
  VALUE rb_info;
  integer info; 

  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  w, werr, info = NumRu::Lapack.slarrk( iw, gl, gu, d, e2, pivmin, reltol)\n    or\n  NumRu::Lapack.slarrk  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLARRK( N, IW, GL, GU, D, E2, PIVMIN, RELTOL, W, WERR, INFO)\n\n*  Purpose\n*  =======\n*\n*  SLARRK computes one eigenvalue of a symmetric tridiagonal\n*  matrix T to suitable accuracy. This is an auxiliary code to be\n*  called from SSTEMR.\n*\n*  To avoid overflow, the matrix must be scaled so that its\n*  largest element is no greater than overflow**(1/2) *\n*  underflow**(1/4) in absolute value, and for greatest\n*  accuracy, it should not be much smaller than that.\n*\n*  See W. Kahan \"Accurate Eigenvalues of a Symmetric Tridiagonal\n*  Matrix\", Report CS41, Computer Science Dept., Stanford\n*  University, July 21, 1966.\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The order of the tridiagonal matrix T.  N >= 0.\n*\n*  IW      (input) INTEGER\n*          The index of the eigenvalues to be returned.\n*\n*  GL      (input) REAL            \n*  GU      (input) REAL            \n*          An upper and a lower bound on the eigenvalue.\n*\n*  D       (input) REAL             array, dimension (N)\n*          The n diagonal elements of the tridiagonal matrix T.\n*\n*  E2      (input) REAL             array, dimension (N-1)\n*          The (n-1) squared off-diagonal elements of the tridiagonal matrix T.\n*\n*  PIVMIN  (input) REAL            \n*          The minimum pivot allowed in the Sturm sequence for T.\n*\n*  RELTOL  (input) REAL            \n*          The minimum relative width of an interval.  When an interval\n*          is narrower than RELTOL times the larger (in\n*          magnitude) endpoint, then it is considered to be\n*          sufficiently small, i.e., converged.  Note: this should\n*          always be at least radix*machine epsilon.\n*\n*  W       (output) REAL            \n*\n*  WERR    (output) REAL            \n*          The error bound on the corresponding eigenvalue approximation\n*          in W.\n*\n*  INFO    (output) INTEGER\n*          = 0:       Eigenvalue converged\n*          = -1:      Eigenvalue did NOT converge\n*\n*  Internal Parameters\n*  ===================\n*\n*  FUDGE   REAL            , default = 2\n*          A \"fudge factor\" to widen the Gershgorin intervals.\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 7)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 7)", argc);
  rb_iw = argv[0];
  rb_gl = argv[1];
  rb_gu = argv[2];
  rb_d = argv[3];
  rb_e2 = argv[4];
  rb_pivmin = argv[5];
  rb_reltol = argv[6];

  iw = NUM2INT(rb_iw);
  gl = (real)NUM2DBL(rb_gl);
  gu = (real)NUM2DBL(rb_gu);
  pivmin = (real)NUM2DBL(rb_pivmin);
  reltol = (real)NUM2DBL(rb_reltol);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (4th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (4th argument) must be %d", 1);
  n = NA_SHAPE0(rb_d);
  if (NA_TYPE(rb_d) != NA_SFLOAT)
    rb_d = na_change_type(rb_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rb_d, real*);
  if (!NA_IsNArray(rb_e2))
    rb_raise(rb_eArgError, "e2 (5th argument) must be NArray");
  if (NA_RANK(rb_e2) != 1)
    rb_raise(rb_eArgError, "rank of e2 (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_e2) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of e2 must be %d", n-1);
  if (NA_TYPE(rb_e2) != NA_SFLOAT)
    rb_e2 = na_change_type(rb_e2, NA_SFLOAT);
  e2 = NA_PTR_TYPE(rb_e2, real*);

  slarrk_(&n, &iw, &gl, &gu, d, e2, &pivmin, &reltol, &w, &werr, &info);

  rb_w = rb_float_new((double)w);
  rb_werr = rb_float_new((double)werr);
  rb_info = INT2NUM(info);
  return rb_ary_new3(3, rb_w, rb_werr, rb_info);
}

void
init_lapack_slarrk(VALUE mLapack){
  rb_define_module_function(mLapack, "slarrk", rb_slarrk, -1);
}
