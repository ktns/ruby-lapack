#include "rb_lapack.h"

static VALUE
rb_slarrf(int argc, VALUE *argv, VALUE self){
  VALUE rb_d;
  real *d; 
  VALUE rb_l;
  real *l; 
  VALUE rb_ld;
  real *ld; 
  VALUE rb_clstrt;
  integer clstrt; 
  VALUE rb_clend;
  integer clend; 
  VALUE rb_w;
  real *w; 
  VALUE rb_wgap;
  real *wgap; 
  VALUE rb_werr;
  real *werr; 
  VALUE rb_spdiam;
  real spdiam; 
  VALUE rb_clgapl;
  real clgapl; 
  VALUE rb_clgapr;
  real clgapr; 
  VALUE rb_pivmin;
  real pivmin; 
  VALUE rb_sigma;
  real sigma; 
  VALUE rb_dplus;
  real *dplus; 
  VALUE rb_lplus;
  real *lplus; 
  VALUE rb_info;
  integer info; 
  VALUE rb_wgap_out__;
  real *wgap_out__;
  real *work;

  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  sigma, dplus, lplus, info, wgap = NumRu::Lapack.slarrf( d, l, ld, clstrt, clend, w, wgap, werr, spdiam, clgapl, clgapr, pivmin)\n    or\n  NumRu::Lapack.slarrf  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLARRF( N, D, L, LD, CLSTRT, CLEND, W, WGAP, WERR, SPDIAM, CLGAPL, CLGAPR, PIVMIN, SIGMA, DPLUS, LPLUS, WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  Given the initial representation L D L^T and its cluster of close\n*  eigenvalues (in a relative measure), W( CLSTRT ), W( CLSTRT+1 ), ...\n*  W( CLEND ), SLARRF finds a new relatively robust representation\n*  L D L^T - SIGMA I = L(+) D(+) L(+)^T such that at least one of the\n*  eigenvalues of L(+) D(+) L(+)^T is relatively isolated.\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The order of the matrix (subblock, if the matrix splitted).\n*\n*  D       (input) REAL             array, dimension (N)\n*          The N diagonal elements of the diagonal matrix D.\n*\n*  L       (input) REAL             array, dimension (N-1)\n*          The (N-1) subdiagonal elements of the unit bidiagonal\n*          matrix L.\n*\n*  LD      (input) REAL             array, dimension (N-1)\n*          The (N-1) elements L(i)*D(i).\n*\n*  CLSTRT  (input) INTEGER\n*          The index of the first eigenvalue in the cluster.\n*\n*  CLEND   (input) INTEGER\n*          The index of the last eigenvalue in the cluster.\n*\n*  W       (input) REAL             array, dimension >=  (CLEND-CLSTRT+1)\n*          The eigenvalue APPROXIMATIONS of L D L^T in ascending order.\n*          W( CLSTRT ) through W( CLEND ) form the cluster of relatively\n*          close eigenalues.\n*\n*  WGAP    (input/output) REAL             array, dimension >=  (CLEND-CLSTRT+1)\n*          The separation from the right neighbor eigenvalue in W.\n*\n*  WERR    (input) REAL             array, dimension >=  (CLEND-CLSTRT+1)\n*          WERR contain the semiwidth of the uncertainty\n*          interval of the corresponding eigenvalue APPROXIMATION in W\n*\n*  SPDIAM (input) estimate of the spectral diameter obtained from the\n*          Gerschgorin intervals\n*\n*  CLGAPL, CLGAPR (input) absolute gap on each end of the cluster.\n*          Set by the calling routine to protect against shifts too close\n*          to eigenvalues outside the cluster.\n*\n*  PIVMIN  (input) DOUBLE PRECISION\n*          The minimum pivot allowed in the Sturm sequence.\n*\n*  SIGMA   (output) REAL            \n*          The shift used to form L(+) D(+) L(+)^T.\n*\n*  DPLUS   (output) REAL             array, dimension (N)\n*          The N diagonal elements of the diagonal matrix D(+).\n*\n*  LPLUS   (output) REAL             array, dimension (N-1)\n*          The first (N-1) elements of LPLUS contain the subdiagonal\n*          elements of the unit bidiagonal matrix L(+).\n*\n*  WORK    (workspace) REAL             array, dimension (2*N)\n*          Workspace.\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Beresford Parlett, University of California, Berkeley, USA\n*     Jim Demmel, University of California, Berkeley, USA\n*     Inderjit Dhillon, University of Texas, Austin, USA\n*     Osni Marques, LBNL/NERSC, USA\n*     Christof Voemel, University of California, Berkeley, USA\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 12)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 12)", argc);
  rb_d = argv[0];
  rb_l = argv[1];
  rb_ld = argv[2];
  rb_clstrt = argv[3];
  rb_clend = argv[4];
  rb_w = argv[5];
  rb_wgap = argv[6];
  rb_werr = argv[7];
  rb_spdiam = argv[8];
  rb_clgapl = argv[9];
  rb_clgapr = argv[10];
  rb_pivmin = argv[11];

  clstrt = NUM2INT(rb_clstrt);
  clend = NUM2INT(rb_clend);
  spdiam = (real)NUM2DBL(rb_spdiam);
  clgapl = (real)NUM2DBL(rb_clgapl);
  clgapr = (real)NUM2DBL(rb_clgapr);
  pivmin = (real)NUM2DBL(rb_pivmin);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (1th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (1th argument) must be %d", 1);
  n = NA_SHAPE0(rb_d);
  if (NA_TYPE(rb_d) != NA_SFLOAT)
    rb_d = na_change_type(rb_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rb_d, real*);
  if (!NA_IsNArray(rb_l))
    rb_raise(rb_eArgError, "l (2th argument) must be NArray");
  if (NA_RANK(rb_l) != 1)
    rb_raise(rb_eArgError, "rank of l (2th argument) must be %d", 1);
  if (NA_SHAPE0(rb_l) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of l must be %d", n-1);
  if (NA_TYPE(rb_l) != NA_SFLOAT)
    rb_l = na_change_type(rb_l, NA_SFLOAT);
  l = NA_PTR_TYPE(rb_l, real*);
  if (!NA_IsNArray(rb_ld))
    rb_raise(rb_eArgError, "ld (3th argument) must be NArray");
  if (NA_RANK(rb_ld) != 1)
    rb_raise(rb_eArgError, "rank of ld (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_ld) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of ld must be %d", n-1);
  if (NA_TYPE(rb_ld) != NA_SFLOAT)
    rb_ld = na_change_type(rb_ld, NA_SFLOAT);
  ld = NA_PTR_TYPE(rb_ld, real*);
  if (!NA_IsNArray(rb_w))
    rb_raise(rb_eArgError, "w (6th argument) must be NArray");
  if (NA_RANK(rb_w) != 1)
    rb_raise(rb_eArgError, "rank of w (6th argument) must be %d", 1);
  if (NA_SHAPE0(rb_w) != (clend-clstrt+1))
    rb_raise(rb_eRuntimeError, "shape 0 of w must be %d", clend-clstrt+1);
  if (NA_TYPE(rb_w) != NA_SFLOAT)
    rb_w = na_change_type(rb_w, NA_SFLOAT);
  w = NA_PTR_TYPE(rb_w, real*);
  if (!NA_IsNArray(rb_wgap))
    rb_raise(rb_eArgError, "wgap (7th argument) must be NArray");
  if (NA_RANK(rb_wgap) != 1)
    rb_raise(rb_eArgError, "rank of wgap (7th argument) must be %d", 1);
  if (NA_SHAPE0(rb_wgap) != (clend-clstrt+1))
    rb_raise(rb_eRuntimeError, "shape 0 of wgap must be %d", clend-clstrt+1);
  if (NA_TYPE(rb_wgap) != NA_SFLOAT)
    rb_wgap = na_change_type(rb_wgap, NA_SFLOAT);
  wgap = NA_PTR_TYPE(rb_wgap, real*);
  if (!NA_IsNArray(rb_werr))
    rb_raise(rb_eArgError, "werr (8th argument) must be NArray");
  if (NA_RANK(rb_werr) != 1)
    rb_raise(rb_eArgError, "rank of werr (8th argument) must be %d", 1);
  if (NA_SHAPE0(rb_werr) != (clend-clstrt+1))
    rb_raise(rb_eRuntimeError, "shape 0 of werr must be %d", clend-clstrt+1);
  if (NA_TYPE(rb_werr) != NA_SFLOAT)
    rb_werr = na_change_type(rb_werr, NA_SFLOAT);
  werr = NA_PTR_TYPE(rb_werr, real*);
  {
    int shape[1];
    shape[0] = n;
    rb_dplus = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  dplus = NA_PTR_TYPE(rb_dplus, real*);
  {
    int shape[1];
    shape[0] = n-1;
    rb_lplus = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  lplus = NA_PTR_TYPE(rb_lplus, real*);
  {
    int shape[1];
    shape[0] = clend-clstrt+1;
    rb_wgap_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  wgap_out__ = NA_PTR_TYPE(rb_wgap_out__, real*);
  MEMCPY(wgap_out__, wgap, real, NA_TOTAL(rb_wgap));
  rb_wgap = rb_wgap_out__;
  wgap = wgap_out__;
  work = ALLOC_N(real, (2*n));

  slarrf_(&n, d, l, ld, &clstrt, &clend, w, wgap, werr, &spdiam, &clgapl, &clgapr, &pivmin, &sigma, dplus, lplus, work, &info);

  free(work);
  rb_sigma = rb_float_new((double)sigma);
  rb_info = INT2NUM(info);
  return rb_ary_new3(5, rb_sigma, rb_dplus, rb_lplus, rb_info, rb_wgap);
}

void
init_lapack_slarrf(VALUE mLapack){
  rb_define_module_function(mLapack, "slarrf", rb_slarrf, -1);
}
