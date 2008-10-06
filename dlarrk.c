#include "rb_lapack.h"

static VALUE
rb_dlarrk(int argc, VALUE *argv, VALUE self){
  VALUE rb_iw;
  integer iw; 
  VALUE rb_gl;
  doublereal gl; 
  VALUE rb_gu;
  doublereal gu; 
  VALUE rb_d;
  doublereal *d; 
  VALUE rb_e2;
  doublereal *e2; 
  VALUE rb_pivmin;
  doublereal pivmin; 
  VALUE rb_reltol;
  doublereal reltol; 
  VALUE rb_w;
  doublereal w; 
  VALUE rb_werr;
  doublereal werr; 
  VALUE rb_info;
  integer info; 

  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  w, werr, info = NumRu::Lapack.dlarrk( iw, gl, gu, d, e2, pivmin, reltol)\n    or\n  NumRu::Lapack.dlarrk  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLARRK( N, IW, GL, GU, D, E2, PIVMIN, RELTOL, W, WERR, INFO)\n\n*  Purpose\n*  =======\n*\n*  DLARRK computes one eigenvalue of a symmetric tridiagonal\n*  matrix T to suitable accuracy. This is an auxiliary code to be\n*  called from DSTEMR.\n*\n*  To avoid overflow, the matrix must be scaled so that its\n*  largest element is no greater than overflow**(1/2) *\n*  underflow**(1/4) in absolute value, and for greatest\n*  accuracy, it should not be much smaller than that.\n*\n*  See W. Kahan \"Accurate Eigenvalues of a Symmetric Tridiagonal\n*  Matrix\", Report CS41, Computer Science Dept., Stanford\n*  University, July 21, 1966.\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The order of the tridiagonal matrix T.  N >= 0.\n*\n*  IW      (input) INTEGER\n*          The index of the eigenvalues to be returned.\n*\n*  GL      (input) DOUBLE PRECISION\n*  GU      (input) DOUBLE PRECISION\n*          An upper and a lower bound on the eigenvalue.\n*\n*  D       (input) DOUBLE PRECISION array, dimension (N)\n*          The n diagonal elements of the tridiagonal matrix T.\n*\n*  E2      (input) DOUBLE PRECISION array, dimension (N-1)\n*          The (n-1) squared off-diagonal elements of the tridiagonal matrix T.\n*\n*  PIVMIN  (input) DOUBLE PRECISION\n*          The minimum pivot allowed in the Sturm sequence for T.\n*\n*  RELTOL  (input) DOUBLE PRECISION\n*          The minimum relative width of an interval.  When an interval\n*          is narrower than RELTOL times the larger (in\n*          magnitude) endpoint, then it is considered to be\n*          sufficiently small, i.e., converged.  Note: this should\n*          always be at least radix*machine epsilon.\n*\n*  W       (output) DOUBLE PRECISION\n*\n*  WERR    (output) DOUBLE PRECISION\n*          The error bound on the corresponding eigenvalue approximation\n*          in W.\n*\n*  INFO    (output) INTEGER\n*          = 0:       Eigenvalue converged\n*          = -1:      Eigenvalue did NOT converge\n*\n*  Internal Parameters\n*  ===================\n*\n*  FUDGE   DOUBLE PRECISION, default = 2\n*          A \"fudge factor\" to widen the Gershgorin intervals.\n*\n\n*  =====================================================================\n*\n\n");
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
  gl = NUM2DBL(rb_gl);
  gu = NUM2DBL(rb_gu);
  pivmin = NUM2DBL(rb_pivmin);
  reltol = NUM2DBL(rb_reltol);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (4th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (4th argument) must be %d", 1);
  n = NA_SHAPE0(rb_d);
  if (NA_TYPE(rb_d) != NA_DFLOAT)
    rb_d = na_change_type(rb_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rb_d, doublereal*);
  if (!NA_IsNArray(rb_e2))
    rb_raise(rb_eArgError, "e2 (5th argument) must be NArray");
  if (NA_RANK(rb_e2) != 1)
    rb_raise(rb_eArgError, "rank of e2 (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_e2) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of e2 must be %d", n-1);
  if (NA_TYPE(rb_e2) != NA_DFLOAT)
    rb_e2 = na_change_type(rb_e2, NA_DFLOAT);
  e2 = NA_PTR_TYPE(rb_e2, doublereal*);

  dlarrk_(&n, &iw, &gl, &gu, d, e2, &pivmin, &reltol, &w, &werr, &info);

  rb_w = rb_float_new((double)w);
  rb_werr = rb_float_new((double)werr);
  rb_info = INT2NUM(info);
  return rb_ary_new3(3, rb_w, rb_werr, rb_info);
}

void
init_lapack_dlarrk(VALUE mLapack){
  rb_define_module_function(mLapack, "dlarrk", rb_dlarrk, -1);
}
