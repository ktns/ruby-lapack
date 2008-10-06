#include "rb_lapack.h"

static VALUE
rb_dlarrj(int argc, VALUE *argv, VALUE self){
  VALUE rb_d;
  doublereal *d; 
  VALUE rb_e2;
  doublereal *e2; 
  VALUE rb_ifirst;
  integer ifirst; 
  VALUE rb_ilast;
  integer ilast; 
  VALUE rb_rtol;
  doublereal rtol; 
  VALUE rb_offset;
  integer offset; 
  VALUE rb_w;
  doublereal *w; 
  VALUE rb_werr;
  doublereal *werr; 
  VALUE rb_pivmin;
  doublereal pivmin; 
  VALUE rb_spdiam;
  doublereal spdiam; 
  VALUE rb_info;
  integer info; 
  VALUE rb_w_out__;
  doublereal *w_out__;
  VALUE rb_werr_out__;
  doublereal *werr_out__;
  doublereal *work;
  integer *iwork;

  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, w, werr = NumRu::Lapack.dlarrj( d, e2, ifirst, ilast, rtol, offset, w, werr, pivmin, spdiam)\n    or\n  NumRu::Lapack.dlarrj  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLARRJ( N, D, E2, IFIRST, ILAST, RTOL, OFFSET, W, WERR, WORK, IWORK, PIVMIN, SPDIAM, INFO )\n\n*  Purpose\n*  =======\n*\n*  Given the initial eigenvalue approximations of T, DLARRJ\n*  does  bisection to refine the eigenvalues of T,\n*  W( IFIRST-OFFSET ) through W( ILAST-OFFSET ), to more accuracy. Initial\n*  guesses for these eigenvalues are input in W, the corresponding estimate\n*  of the error in these guesses in WERR. During bisection, intervals\n*  [left, right] are maintained by storing their mid-points and\n*  semi-widths in the arrays W and WERR respectively.\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The order of the matrix.\n*\n*  D       (input) DOUBLE PRECISION array, dimension (N)\n*          The N diagonal elements of T.\n*\n*  E2      (input) DOUBLE PRECISION array, dimension (N-1)\n*          The Squares of the (N-1) subdiagonal elements of T.\n*\n*  IFIRST  (input) INTEGER\n*          The index of the first eigenvalue to be computed.\n*\n*  ILAST   (input) INTEGER\n*          The index of the last eigenvalue to be computed.\n*\n*  RTOL   (input) DOUBLE PRECISION\n*          Tolerance for the convergence of the bisection intervals.\n*          An interval [LEFT,RIGHT] has converged if\n*          RIGHT-LEFT.LT.RTOL*MAX(|LEFT|,|RIGHT|).\n*\n*  OFFSET  (input) INTEGER\n*          Offset for the arrays W and WERR, i.e., the IFIRST-OFFSET\n*          through ILAST-OFFSET elements of these arrays are to be used.\n*\n*  W       (input/output) DOUBLE PRECISION array, dimension (N)\n*          On input, W( IFIRST-OFFSET ) through W( ILAST-OFFSET ) are\n*          estimates of the eigenvalues of L D L^T indexed IFIRST through\n*          ILAST.\n*          On output, these estimates are refined.\n*\n*  WERR    (input/output) DOUBLE PRECISION array, dimension (N)\n*          On input, WERR( IFIRST-OFFSET ) through WERR( ILAST-OFFSET ) are\n*          the errors in the estimates of the corresponding elements in W.\n*          On output, these errors are refined.\n*\n*  WORK    (workspace) DOUBLE PRECISION array, dimension (2*N)\n*          Workspace.\n*\n*  IWORK   (workspace) INTEGER array, dimension (2*N)\n*          Workspace.\n*\n*  PIVMIN  (input) DOUBLE PRECISION\n*          The minimum pivot in the Sturm sequence for T.\n*\n*  SPDIAM  (input) DOUBLE PRECISION\n*          The spectral diameter of T.\n*\n*  INFO    (output) INTEGER\n*          Error flag.\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Beresford Parlett, University of California, Berkeley, USA\n*     Jim Demmel, University of California, Berkeley, USA\n*     Inderjit Dhillon, University of Texas, Austin, USA\n*     Osni Marques, LBNL/NERSC, USA\n*     Christof Voemel, University of California, Berkeley, USA\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 10)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 10)", argc);
  rb_d = argv[0];
  rb_e2 = argv[1];
  rb_ifirst = argv[2];
  rb_ilast = argv[3];
  rb_rtol = argv[4];
  rb_offset = argv[5];
  rb_w = argv[6];
  rb_werr = argv[7];
  rb_pivmin = argv[8];
  rb_spdiam = argv[9];

  ifirst = NUM2INT(rb_ifirst);
  ilast = NUM2INT(rb_ilast);
  rtol = NUM2DBL(rb_rtol);
  offset = NUM2INT(rb_offset);
  pivmin = NUM2DBL(rb_pivmin);
  spdiam = NUM2DBL(rb_spdiam);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (1th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (1th argument) must be %d", 1);
  n = NA_SHAPE0(rb_d);
  if (NA_TYPE(rb_d) != NA_DFLOAT)
    rb_d = na_change_type(rb_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rb_d, doublereal*);
  if (!NA_IsNArray(rb_e2))
    rb_raise(rb_eArgError, "e2 (2th argument) must be NArray");
  if (NA_RANK(rb_e2) != 1)
    rb_raise(rb_eArgError, "rank of e2 (2th argument) must be %d", 1);
  if (NA_SHAPE0(rb_e2) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of e2 must be %d", n-1);
  if (NA_TYPE(rb_e2) != NA_DFLOAT)
    rb_e2 = na_change_type(rb_e2, NA_DFLOAT);
  e2 = NA_PTR_TYPE(rb_e2, doublereal*);
  if (!NA_IsNArray(rb_w))
    rb_raise(rb_eArgError, "w (7th argument) must be NArray");
  if (NA_RANK(rb_w) != 1)
    rb_raise(rb_eArgError, "rank of w (7th argument) must be %d", 1);
  if (NA_SHAPE0(rb_w) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of w must be the same as shape 0 of d");
  if (NA_TYPE(rb_w) != NA_DFLOAT)
    rb_w = na_change_type(rb_w, NA_DFLOAT);
  w = NA_PTR_TYPE(rb_w, doublereal*);
  if (!NA_IsNArray(rb_werr))
    rb_raise(rb_eArgError, "werr (8th argument) must be NArray");
  if (NA_RANK(rb_werr) != 1)
    rb_raise(rb_eArgError, "rank of werr (8th argument) must be %d", 1);
  if (NA_SHAPE0(rb_werr) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of werr must be the same as shape 0 of d");
  if (NA_TYPE(rb_werr) != NA_DFLOAT)
    rb_werr = na_change_type(rb_werr, NA_DFLOAT);
  werr = NA_PTR_TYPE(rb_werr, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_w_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  w_out__ = NA_PTR_TYPE(rb_w_out__, doublereal*);
  MEMCPY(w_out__, w, doublereal, NA_TOTAL(rb_w));
  rb_w = rb_w_out__;
  w = w_out__;
  {
    int shape[1];
    shape[0] = n;
    rb_werr_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  werr_out__ = NA_PTR_TYPE(rb_werr_out__, doublereal*);
  MEMCPY(werr_out__, werr, doublereal, NA_TOTAL(rb_werr));
  rb_werr = rb_werr_out__;
  werr = werr_out__;
  work = ALLOC_N(doublereal, (2*n));
  iwork = ALLOC_N(integer, (2*n));

  dlarrj_(&n, d, e2, &ifirst, &ilast, &rtol, &offset, w, werr, work, iwork, &pivmin, &spdiam, &info);

  free(work);
  free(iwork);
  rb_info = INT2NUM(info);
  return rb_ary_new3(3, rb_info, rb_w, rb_werr);
}

void
init_lapack_dlarrj(VALUE mLapack){
  rb_define_module_function(mLapack, "dlarrj", rb_dlarrj, -1);
}
