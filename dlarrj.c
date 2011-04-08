#include "rb_lapack.h"

extern VOID dlarrj_(integer *n, doublereal *d, doublereal *e2, integer *ifirst, integer *ilast, doublereal *rtol, integer *offset, doublereal *w, doublereal *werr, doublereal *work, integer *iwork, doublereal *pivmin, doublereal *spdiam, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dlarrj(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_d;
  doublereal *d; 
  VALUE rblapack_e2;
  doublereal *e2; 
  VALUE rblapack_ifirst;
  integer ifirst; 
  VALUE rblapack_ilast;
  integer ilast; 
  VALUE rblapack_rtol;
  doublereal rtol; 
  VALUE rblapack_offset;
  integer offset; 
  VALUE rblapack_w;
  doublereal *w; 
  VALUE rblapack_werr;
  doublereal *werr; 
  VALUE rblapack_pivmin;
  doublereal pivmin; 
  VALUE rblapack_spdiam;
  doublereal spdiam; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_w_out__;
  doublereal *w_out__;
  VALUE rblapack_werr_out__;
  doublereal *werr_out__;
  doublereal *work;
  integer *iwork;

  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, w, werr = NumRu::Lapack.dlarrj( d, e2, ifirst, ilast, rtol, offset, w, werr, pivmin, spdiam, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLARRJ( N, D, E2, IFIRST, ILAST, RTOL, OFFSET, W, WERR, WORK, IWORK, PIVMIN, SPDIAM, INFO )\n\n*  Purpose\n*  =======\n*\n*  Given the initial eigenvalue approximations of T, DLARRJ\n*  does  bisection to refine the eigenvalues of T,\n*  W( IFIRST-OFFSET ) through W( ILAST-OFFSET ), to more accuracy. Initial\n*  guesses for these eigenvalues are input in W, the corresponding estimate\n*  of the error in these guesses in WERR. During bisection, intervals\n*  [left, right] are maintained by storing their mid-points and\n*  semi-widths in the arrays W and WERR respectively.\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The order of the matrix.\n*\n*  D       (input) DOUBLE PRECISION array, dimension (N)\n*          The N diagonal elements of T.\n*\n*  E2      (input) DOUBLE PRECISION array, dimension (N-1)\n*          The Squares of the (N-1) subdiagonal elements of T.\n*\n*  IFIRST  (input) INTEGER\n*          The index of the first eigenvalue to be computed.\n*\n*  ILAST   (input) INTEGER\n*          The index of the last eigenvalue to be computed.\n*\n*  RTOL    (input) DOUBLE PRECISION\n*          Tolerance for the convergence of the bisection intervals.\n*          An interval [LEFT,RIGHT] has converged if\n*          RIGHT-LEFT.LT.RTOL*MAX(|LEFT|,|RIGHT|).\n*\n*  OFFSET  (input) INTEGER\n*          Offset for the arrays W and WERR, i.e., the IFIRST-OFFSET\n*          through ILAST-OFFSET elements of these arrays are to be used.\n*\n*  W       (input/output) DOUBLE PRECISION array, dimension (N)\n*          On input, W( IFIRST-OFFSET ) through W( ILAST-OFFSET ) are\n*          estimates of the eigenvalues of L D L^T indexed IFIRST through\n*          ILAST.\n*          On output, these estimates are refined.\n*\n*  WERR    (input/output) DOUBLE PRECISION array, dimension (N)\n*          On input, WERR( IFIRST-OFFSET ) through WERR( ILAST-OFFSET ) are\n*          the errors in the estimates of the corresponding elements in W.\n*          On output, these errors are refined.\n*\n*  WORK    (workspace) DOUBLE PRECISION array, dimension (2*N)\n*          Workspace.\n*\n*  IWORK   (workspace) INTEGER array, dimension (2*N)\n*          Workspace.\n*\n*  PIVMIN  (input) DOUBLE PRECISION\n*          The minimum pivot in the Sturm sequence for T.\n*\n*  SPDIAM  (input) DOUBLE PRECISION\n*          The spectral diameter of T.\n*\n*  INFO    (output) INTEGER\n*          Error flag.\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Beresford Parlett, University of California, Berkeley, USA\n*     Jim Demmel, University of California, Berkeley, USA\n*     Inderjit Dhillon, University of Texas, Austin, USA\n*     Osni Marques, LBNL/NERSC, USA\n*     Christof Voemel, University of California, Berkeley, USA\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, w, werr = NumRu::Lapack.dlarrj( d, e2, ifirst, ilast, rtol, offset, w, werr, pivmin, spdiam, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 10)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 10)", argc);
  rblapack_d = argv[0];
  rblapack_e2 = argv[1];
  rblapack_ifirst = argv[2];
  rblapack_ilast = argv[3];
  rblapack_rtol = argv[4];
  rblapack_offset = argv[5];
  rblapack_w = argv[6];
  rblapack_werr = argv[7];
  rblapack_pivmin = argv[8];
  rblapack_spdiam = argv[9];
  if (rb_options != Qnil) {
  }

  ilast = NUM2INT(rblapack_ilast);
  if (!NA_IsNArray(rblapack_w))
    rb_raise(rb_eArgError, "w (7th argument) must be NArray");
  if (NA_RANK(rblapack_w) != 1)
    rb_raise(rb_eArgError, "rank of w (7th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_w);
  if (NA_TYPE(rblapack_w) != NA_DFLOAT)
    rblapack_w = na_change_type(rblapack_w, NA_DFLOAT);
  w = NA_PTR_TYPE(rblapack_w, doublereal*);
  rtol = NUM2DBL(rblapack_rtol);
  offset = NUM2INT(rblapack_offset);
  spdiam = NUM2DBL(rblapack_spdiam);
  pivmin = NUM2DBL(rblapack_pivmin);
  if (!NA_IsNArray(rblapack_werr))
    rb_raise(rb_eArgError, "werr (8th argument) must be NArray");
  if (NA_RANK(rblapack_werr) != 1)
    rb_raise(rb_eArgError, "rank of werr (8th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_werr) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of werr must be the same as shape 0 of w");
  if (NA_TYPE(rblapack_werr) != NA_DFLOAT)
    rblapack_werr = na_change_type(rblapack_werr, NA_DFLOAT);
  werr = NA_PTR_TYPE(rblapack_werr, doublereal*);
  ifirst = NUM2INT(rblapack_ifirst);
  if (!NA_IsNArray(rblapack_d))
    rb_raise(rb_eArgError, "d (1th argument) must be NArray");
  if (NA_RANK(rblapack_d) != 1)
    rb_raise(rb_eArgError, "rank of d (1th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_d) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of d must be the same as shape 0 of w");
  if (NA_TYPE(rblapack_d) != NA_DFLOAT)
    rblapack_d = na_change_type(rblapack_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rblapack_d, doublereal*);
  if (!NA_IsNArray(rblapack_e2))
    rb_raise(rb_eArgError, "e2 (2th argument) must be NArray");
  if (NA_RANK(rblapack_e2) != 1)
    rb_raise(rb_eArgError, "rank of e2 (2th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_e2) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of e2 must be %d", n-1);
  if (NA_TYPE(rblapack_e2) != NA_DFLOAT)
    rblapack_e2 = na_change_type(rblapack_e2, NA_DFLOAT);
  e2 = NA_PTR_TYPE(rblapack_e2, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_w_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  w_out__ = NA_PTR_TYPE(rblapack_w_out__, doublereal*);
  MEMCPY(w_out__, w, doublereal, NA_TOTAL(rblapack_w));
  rblapack_w = rblapack_w_out__;
  w = w_out__;
  {
    int shape[1];
    shape[0] = n;
    rblapack_werr_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  werr_out__ = NA_PTR_TYPE(rblapack_werr_out__, doublereal*);
  MEMCPY(werr_out__, werr, doublereal, NA_TOTAL(rblapack_werr));
  rblapack_werr = rblapack_werr_out__;
  werr = werr_out__;
  work = ALLOC_N(doublereal, (2*n));
  iwork = ALLOC_N(integer, (2*n));

  dlarrj_(&n, d, e2, &ifirst, &ilast, &rtol, &offset, w, werr, work, iwork, &pivmin, &spdiam, &info);

  free(work);
  free(iwork);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(3, rblapack_info, rblapack_w, rblapack_werr);
}

void
init_lapack_dlarrj(VALUE mLapack){
  rb_define_module_function(mLapack, "dlarrj", rblapack_dlarrj, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
