#include "rb_lapack.h"

extern VOID slarrb_(integer *n, real *d, real *lld, integer *ifirst, integer *ilast, real *rtol1, real *rtol2, integer *offset, real *w, real *wgap, real *werr, real *work, integer *iwork, real *pivmin, real *spdiam, integer *twist, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_slarrb(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_d;
  real *d; 
  VALUE rblapack_lld;
  real *lld; 
  VALUE rblapack_ifirst;
  integer ifirst; 
  VALUE rblapack_ilast;
  integer ilast; 
  VALUE rblapack_rtol1;
  real rtol1; 
  VALUE rblapack_rtol2;
  real rtol2; 
  VALUE rblapack_offset;
  integer offset; 
  VALUE rblapack_w;
  real *w; 
  VALUE rblapack_wgap;
  real *wgap; 
  VALUE rblapack_werr;
  real *werr; 
  VALUE rblapack_pivmin;
  real pivmin; 
  VALUE rblapack_spdiam;
  real spdiam; 
  VALUE rblapack_twist;
  integer twist; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_w_out__;
  real *w_out__;
  VALUE rblapack_wgap_out__;
  real *wgap_out__;
  VALUE rblapack_werr_out__;
  real *werr_out__;
  real *work;
  integer *iwork;

  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, w, wgap, werr = NumRu::Lapack.slarrb( d, lld, ifirst, ilast, rtol1, rtol2, offset, w, wgap, werr, pivmin, spdiam, twist, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLARRB( N, D, LLD, IFIRST, ILAST, RTOL1, RTOL2, OFFSET, W, WGAP, WERR, WORK, IWORK, PIVMIN, SPDIAM, TWIST, INFO )\n\n*  Purpose\n*  =======\n*\n*  Given the relatively robust representation(RRR) L D L^T, SLARRB\n*  does \"limited\" bisection to refine the eigenvalues of L D L^T,\n*  W( IFIRST-OFFSET ) through W( ILAST-OFFSET ), to more accuracy. Initial\n*  guesses for these eigenvalues are input in W, the corresponding estimate\n*  of the error in these guesses and their gaps are input in WERR\n*  and WGAP, respectively. During bisection, intervals\n*  [left, right] are maintained by storing their mid-points and\n*  semi-widths in the arrays W and WERR respectively.\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The order of the matrix.\n*\n*  D       (input) REAL             array, dimension (N)\n*          The N diagonal elements of the diagonal matrix D.\n*\n*  LLD     (input) REAL             array, dimension (N-1)\n*          The (N-1) elements L(i)*L(i)*D(i).\n*\n*  IFIRST  (input) INTEGER\n*          The index of the first eigenvalue to be computed.\n*\n*  ILAST   (input) INTEGER\n*          The index of the last eigenvalue to be computed.\n*\n*  RTOL1   (input) REAL            \n*  RTOL2   (input) REAL            \n*          Tolerance for the convergence of the bisection intervals.\n*          An interval [LEFT,RIGHT] has converged if\n*          RIGHT-LEFT.LT.MAX( RTOL1*GAP, RTOL2*MAX(|LEFT|,|RIGHT|) )\n*          where GAP is the (estimated) distance to the nearest\n*          eigenvalue.\n*\n*  OFFSET  (input) INTEGER\n*          Offset for the arrays W, WGAP and WERR, i.e., the IFIRST-OFFSET\n*          through ILAST-OFFSET elements of these arrays are to be used.\n*\n*  W       (input/output) REAL             array, dimension (N)\n*          On input, W( IFIRST-OFFSET ) through W( ILAST-OFFSET ) are\n*          estimates of the eigenvalues of L D L^T indexed IFIRST throug\n*          ILAST.\n*          On output, these estimates are refined.\n*\n*  WGAP    (input/output) REAL             array, dimension (N-1)\n*          On input, the (estimated) gaps between consecutive\n*          eigenvalues of L D L^T, i.e., WGAP(I-OFFSET) is the gap between\n*          eigenvalues I and I+1. Note that if IFIRST.EQ.ILAST\n*          then WGAP(IFIRST-OFFSET) must be set to ZERO.\n*          On output, these gaps are refined.\n*\n*  WERR    (input/output) REAL             array, dimension (N)\n*          On input, WERR( IFIRST-OFFSET ) through WERR( ILAST-OFFSET ) are\n*          the errors in the estimates of the corresponding elements in W.\n*          On output, these errors are refined.\n*\n*  WORK    (workspace) REAL             array, dimension (2*N)\n*          Workspace.\n*\n*  IWORK   (workspace) INTEGER array, dimension (2*N)\n*          Workspace.\n*\n*  PIVMIN  (input) REAL\n*          The minimum pivot in the Sturm sequence.\n*\n*  SPDIAM  (input) REAL\n*          The spectral diameter of the matrix.\n*\n*  TWIST   (input) INTEGER\n*          The twist index for the twisted factorization that is used\n*          for the negcount.\n*          TWIST = N: Compute negcount from L D L^T - LAMBDA I = L+ D+ L+^T\n*          TWIST = 1: Compute negcount from L D L^T - LAMBDA I = U- D- U-^T\n*          TWIST = R: Compute negcount from L D L^T - LAMBDA I = N(r) D(r) N(r)\n*\n*  INFO    (output) INTEGER\n*          Error flag.\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Beresford Parlett, University of California, Berkeley, USA\n*     Jim Demmel, University of California, Berkeley, USA\n*     Inderjit Dhillon, University of Texas, Austin, USA\n*     Osni Marques, LBNL/NERSC, USA\n*     Christof Voemel, University of California, Berkeley, USA\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, w, wgap, werr = NumRu::Lapack.slarrb( d, lld, ifirst, ilast, rtol1, rtol2, offset, w, wgap, werr, pivmin, spdiam, twist, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 13)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 13)", argc);
  rblapack_d = argv[0];
  rblapack_lld = argv[1];
  rblapack_ifirst = argv[2];
  rblapack_ilast = argv[3];
  rblapack_rtol1 = argv[4];
  rblapack_rtol2 = argv[5];
  rblapack_offset = argv[6];
  rblapack_w = argv[7];
  rblapack_wgap = argv[8];
  rblapack_werr = argv[9];
  rblapack_pivmin = argv[10];
  rblapack_spdiam = argv[11];
  rblapack_twist = argv[12];
  if (rb_options != Qnil) {
  }

  ilast = NUM2INT(rblapack_ilast);
  if (!NA_IsNArray(rblapack_w))
    rb_raise(rb_eArgError, "w (8th argument) must be NArray");
  if (NA_RANK(rblapack_w) != 1)
    rb_raise(rb_eArgError, "rank of w (8th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_w);
  if (NA_TYPE(rblapack_w) != NA_SFLOAT)
    rblapack_w = na_change_type(rblapack_w, NA_SFLOAT);
  w = NA_PTR_TYPE(rblapack_w, real*);
  rtol2 = (real)NUM2DBL(rblapack_rtol2);
  spdiam = (real)NUM2DBL(rblapack_spdiam);
  offset = NUM2INT(rblapack_offset);
  if (!NA_IsNArray(rblapack_d))
    rb_raise(rb_eArgError, "d (1th argument) must be NArray");
  if (NA_RANK(rblapack_d) != 1)
    rb_raise(rb_eArgError, "rank of d (1th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_d) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of d must be the same as shape 0 of w");
  if (NA_TYPE(rblapack_d) != NA_SFLOAT)
    rblapack_d = na_change_type(rblapack_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rblapack_d, real*);
  pivmin = (real)NUM2DBL(rblapack_pivmin);
  twist = NUM2INT(rblapack_twist);
  rtol1 = (real)NUM2DBL(rblapack_rtol1);
  ifirst = NUM2INT(rblapack_ifirst);
  if (!NA_IsNArray(rblapack_werr))
    rb_raise(rb_eArgError, "werr (10th argument) must be NArray");
  if (NA_RANK(rblapack_werr) != 1)
    rb_raise(rb_eArgError, "rank of werr (10th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_werr) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of werr must be the same as shape 0 of w");
  if (NA_TYPE(rblapack_werr) != NA_SFLOAT)
    rblapack_werr = na_change_type(rblapack_werr, NA_SFLOAT);
  werr = NA_PTR_TYPE(rblapack_werr, real*);
  if (!NA_IsNArray(rblapack_lld))
    rb_raise(rb_eArgError, "lld (2th argument) must be NArray");
  if (NA_RANK(rblapack_lld) != 1)
    rb_raise(rb_eArgError, "rank of lld (2th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_lld) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of lld must be %d", n-1);
  if (NA_TYPE(rblapack_lld) != NA_SFLOAT)
    rblapack_lld = na_change_type(rblapack_lld, NA_SFLOAT);
  lld = NA_PTR_TYPE(rblapack_lld, real*);
  if (!NA_IsNArray(rblapack_wgap))
    rb_raise(rb_eArgError, "wgap (9th argument) must be NArray");
  if (NA_RANK(rblapack_wgap) != 1)
    rb_raise(rb_eArgError, "rank of wgap (9th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_wgap) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of wgap must be %d", n-1);
  if (NA_TYPE(rblapack_wgap) != NA_SFLOAT)
    rblapack_wgap = na_change_type(rblapack_wgap, NA_SFLOAT);
  wgap = NA_PTR_TYPE(rblapack_wgap, real*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_w_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  w_out__ = NA_PTR_TYPE(rblapack_w_out__, real*);
  MEMCPY(w_out__, w, real, NA_TOTAL(rblapack_w));
  rblapack_w = rblapack_w_out__;
  w = w_out__;
  {
    int shape[1];
    shape[0] = n-1;
    rblapack_wgap_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  wgap_out__ = NA_PTR_TYPE(rblapack_wgap_out__, real*);
  MEMCPY(wgap_out__, wgap, real, NA_TOTAL(rblapack_wgap));
  rblapack_wgap = rblapack_wgap_out__;
  wgap = wgap_out__;
  {
    int shape[1];
    shape[0] = n;
    rblapack_werr_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  werr_out__ = NA_PTR_TYPE(rblapack_werr_out__, real*);
  MEMCPY(werr_out__, werr, real, NA_TOTAL(rblapack_werr));
  rblapack_werr = rblapack_werr_out__;
  werr = werr_out__;
  work = ALLOC_N(real, (2*n));
  iwork = ALLOC_N(integer, (2*n));

  slarrb_(&n, d, lld, &ifirst, &ilast, &rtol1, &rtol2, &offset, w, wgap, werr, work, iwork, &pivmin, &spdiam, &twist, &info);

  free(work);
  free(iwork);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(4, rblapack_info, rblapack_w, rblapack_wgap, rblapack_werr);
}

void
init_lapack_slarrb(VALUE mLapack){
  rb_define_module_function(mLapack, "slarrb", rblapack_slarrb, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
