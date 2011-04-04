#include "rb_lapack.h"

extern VOID slarrj_(integer *n, real *d, real *e2, integer *ifirst, integer *ilast, real *rtol, integer *offset, real *w, real *werr, real *work, integer *iwork, real *pivmin, real *spdiam, integer *info);

static VALUE
rb_slarrj(int argc, VALUE *argv, VALUE self){
  VALUE rb_d;
  real *d; 
  VALUE rb_e2;
  real *e2; 
  VALUE rb_ifirst;
  integer ifirst; 
  VALUE rb_ilast;
  integer ilast; 
  VALUE rb_rtol;
  real rtol; 
  VALUE rb_offset;
  integer offset; 
  VALUE rb_w;
  real *w; 
  VALUE rb_werr;
  real *werr; 
  VALUE rb_pivmin;
  real pivmin; 
  VALUE rb_spdiam;
  real spdiam; 
  VALUE rb_info;
  integer info; 
  VALUE rb_w_out__;
  real *w_out__;
  VALUE rb_werr_out__;
  real *werr_out__;
  real *work;
  integer *iwork;

  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, w, werr = NumRu::Lapack.slarrj( d, e2, ifirst, ilast, rtol, offset, w, werr, pivmin, spdiam)\n    or\n  NumRu::Lapack.slarrj  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLARRJ( N, D, E2, IFIRST, ILAST, RTOL, OFFSET, W, WERR, WORK, IWORK, PIVMIN, SPDIAM, INFO )\n\n*  Purpose\n*  =======\n*\n*  Given the initial eigenvalue approximations of T, SLARRJ\n*  does  bisection to refine the eigenvalues of T,\n*  W( IFIRST-OFFSET ) through W( ILAST-OFFSET ), to more accuracy. Initial\n*  guesses for these eigenvalues are input in W, the corresponding estimate\n*  of the error in these guesses in WERR. During bisection, intervals\n*  [left, right] are maintained by storing their mid-points and\n*  semi-widths in the arrays W and WERR respectively.\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The order of the matrix.\n*\n*  D       (input) REAL             array, dimension (N)\n*          The N diagonal elements of T.\n*\n*  E2      (input) REAL             array, dimension (N-1)\n*          The Squares of the (N-1) subdiagonal elements of T.\n*\n*  IFIRST  (input) INTEGER\n*          The index of the first eigenvalue to be computed.\n*\n*  ILAST   (input) INTEGER\n*          The index of the last eigenvalue to be computed.\n*\n*  RTOL   (input) REAL            \n*          Tolerance for the convergence of the bisection intervals.\n*          An interval [LEFT,RIGHT] has converged if\n*          RIGHT-LEFT.LT.RTOL*MAX(|LEFT|,|RIGHT|).\n*\n*  OFFSET  (input) INTEGER\n*          Offset for the arrays W and WERR, i.e., the IFIRST-OFFSET\n*          through ILAST-OFFSET elements of these arrays are to be used.\n*\n*  W       (input/output) REAL             array, dimension (N)\n*          On input, W( IFIRST-OFFSET ) through W( ILAST-OFFSET ) are\n*          estimates of the eigenvalues of L D L^T indexed IFIRST through\n*          ILAST.\n*          On output, these estimates are refined.\n*\n*  WERR    (input/output) REAL             array, dimension (N)\n*          On input, WERR( IFIRST-OFFSET ) through WERR( ILAST-OFFSET ) are\n*          the errors in the estimates of the corresponding elements in W.\n*          On output, these errors are refined.\n*\n*  WORK    (workspace) REAL             array, dimension (2*N)\n*          Workspace.\n*\n*  IWORK   (workspace) INTEGER array, dimension (2*N)\n*          Workspace.\n*\n*  PIVMIN  (input) REAL\n*          The minimum pivot in the Sturm sequence for T.\n*\n*  SPDIAM  (input) REAL\n*          The spectral diameter of T.\n*\n*  INFO    (output) INTEGER\n*          Error flag.\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Beresford Parlett, University of California, Berkeley, USA\n*     Jim Demmel, University of California, Berkeley, USA\n*     Inderjit Dhillon, University of Texas, Austin, USA\n*     Osni Marques, LBNL/NERSC, USA\n*     Christof Voemel, University of California, Berkeley, USA\n*\n*  =====================================================================\n*\n\n");
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

  ilast = NUM2INT(rb_ilast);
  if (!NA_IsNArray(rb_w))
    rb_raise(rb_eArgError, "w (7th argument) must be NArray");
  if (NA_RANK(rb_w) != 1)
    rb_raise(rb_eArgError, "rank of w (7th argument) must be %d", 1);
  n = NA_SHAPE0(rb_w);
  if (NA_TYPE(rb_w) != NA_SFLOAT)
    rb_w = na_change_type(rb_w, NA_SFLOAT);
  w = NA_PTR_TYPE(rb_w, real*);
  rtol = (real)NUM2DBL(rb_rtol);
  offset = NUM2INT(rb_offset);
  spdiam = (real)NUM2DBL(rb_spdiam);
  pivmin = (real)NUM2DBL(rb_pivmin);
  if (!NA_IsNArray(rb_werr))
    rb_raise(rb_eArgError, "werr (8th argument) must be NArray");
  if (NA_RANK(rb_werr) != 1)
    rb_raise(rb_eArgError, "rank of werr (8th argument) must be %d", 1);
  if (NA_SHAPE0(rb_werr) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of werr must be the same as shape 0 of w");
  if (NA_TYPE(rb_werr) != NA_SFLOAT)
    rb_werr = na_change_type(rb_werr, NA_SFLOAT);
  werr = NA_PTR_TYPE(rb_werr, real*);
  ifirst = NUM2INT(rb_ifirst);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (1th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (1th argument) must be %d", 1);
  if (NA_SHAPE0(rb_d) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of d must be the same as shape 0 of w");
  if (NA_TYPE(rb_d) != NA_SFLOAT)
    rb_d = na_change_type(rb_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rb_d, real*);
  if (!NA_IsNArray(rb_e2))
    rb_raise(rb_eArgError, "e2 (2th argument) must be NArray");
  if (NA_RANK(rb_e2) != 1)
    rb_raise(rb_eArgError, "rank of e2 (2th argument) must be %d", 1);
  if (NA_SHAPE0(rb_e2) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of e2 must be %d", n-1);
  if (NA_TYPE(rb_e2) != NA_SFLOAT)
    rb_e2 = na_change_type(rb_e2, NA_SFLOAT);
  e2 = NA_PTR_TYPE(rb_e2, real*);
  {
    int shape[1];
    shape[0] = n;
    rb_w_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  w_out__ = NA_PTR_TYPE(rb_w_out__, real*);
  MEMCPY(w_out__, w, real, NA_TOTAL(rb_w));
  rb_w = rb_w_out__;
  w = w_out__;
  {
    int shape[1];
    shape[0] = n;
    rb_werr_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  werr_out__ = NA_PTR_TYPE(rb_werr_out__, real*);
  MEMCPY(werr_out__, werr, real, NA_TOTAL(rb_werr));
  rb_werr = rb_werr_out__;
  werr = werr_out__;
  work = ALLOC_N(real, (2*n));
  iwork = ALLOC_N(integer, (2*n));

  slarrj_(&n, d, e2, &ifirst, &ilast, &rtol, &offset, w, werr, work, iwork, &pivmin, &spdiam, &info);

  free(work);
  free(iwork);
  rb_info = INT2NUM(info);
  return rb_ary_new3(3, rb_info, rb_w, rb_werr);
}

void
init_lapack_slarrj(VALUE mLapack){
  rb_define_module_function(mLapack, "slarrj", rb_slarrj, -1);
}
