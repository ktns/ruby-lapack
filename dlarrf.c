#include "rb_lapack.h"

extern VOID dlarrf_(integer *n, doublereal *d, doublereal *l, doublereal *ld, integer *clstrt, integer *clend, doublereal *w, doublereal *wgap, doublereal *werr, doublereal *spdiam, doublereal *clgapl, doublereal *clgapr, doublereal *pivmin, doublereal *sigma, doublereal *dplus, doublereal *lplus, doublereal *work, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dlarrf(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_d;
  doublereal *d; 
  VALUE rblapack_l;
  doublereal *l; 
  VALUE rblapack_ld;
  doublereal *ld; 
  VALUE rblapack_clstrt;
  integer clstrt; 
  VALUE rblapack_clend;
  integer clend; 
  VALUE rblapack_w;
  doublereal *w; 
  VALUE rblapack_wgap;
  doublereal *wgap; 
  VALUE rblapack_werr;
  doublereal *werr; 
  VALUE rblapack_spdiam;
  doublereal spdiam; 
  VALUE rblapack_clgapl;
  doublereal clgapl; 
  VALUE rblapack_clgapr;
  doublereal clgapr; 
  VALUE rblapack_pivmin;
  doublereal pivmin; 
  VALUE rblapack_sigma;
  doublereal sigma; 
  VALUE rblapack_dplus;
  doublereal *dplus; 
  VALUE rblapack_lplus;
  doublereal *lplus; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_wgap_out__;
  doublereal *wgap_out__;
  doublereal *work;

  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  sigma, dplus, lplus, info, wgap = NumRu::Lapack.dlarrf( d, l, ld, clstrt, clend, w, wgap, werr, spdiam, clgapl, clgapr, pivmin, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLARRF( N, D, L, LD, CLSTRT, CLEND, W, WGAP, WERR, SPDIAM, CLGAPL, CLGAPR, PIVMIN, SIGMA, DPLUS, LPLUS, WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  Given the initial representation L D L^T and its cluster of close\n*  eigenvalues (in a relative measure), W( CLSTRT ), W( CLSTRT+1 ), ...\n*  W( CLEND ), DLARRF finds a new relatively robust representation\n*  L D L^T - SIGMA I = L(+) D(+) L(+)^T such that at least one of the\n*  eigenvalues of L(+) D(+) L(+)^T is relatively isolated.\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The order of the matrix (subblock, if the matrix splitted).\n*\n*  D       (input) DOUBLE PRECISION array, dimension (N)\n*          The N diagonal elements of the diagonal matrix D.\n*\n*  L       (input) DOUBLE PRECISION array, dimension (N-1)\n*          The (N-1) subdiagonal elements of the unit bidiagonal\n*          matrix L.\n*\n*  LD      (input) DOUBLE PRECISION array, dimension (N-1)\n*          The (N-1) elements L(i)*D(i).\n*\n*  CLSTRT  (input) INTEGER\n*          The index of the first eigenvalue in the cluster.\n*\n*  CLEND   (input) INTEGER\n*          The index of the last eigenvalue in the cluster.\n*\n*  W       (input) DOUBLE PRECISION array, dimension\n*          dimension is >=  (CLEND-CLSTRT+1)\n*          The eigenvalue APPROXIMATIONS of L D L^T in ascending order.\n*          W( CLSTRT ) through W( CLEND ) form the cluster of relatively\n*          close eigenalues.\n*\n*  WGAP    (input/output) DOUBLE PRECISION array, dimension\n*          dimension is >=  (CLEND-CLSTRT+1)\n*          The separation from the right neighbor eigenvalue in W.\n*\n*  WERR    (input) DOUBLE PRECISION array, dimension\n*          dimension is  >=  (CLEND-CLSTRT+1)\n*          WERR contain the semiwidth of the uncertainty\n*          interval of the corresponding eigenvalue APPROXIMATION in W\n*\n*  SPDIAM  (input) DOUBLE PRECISION\n*          estimate of the spectral diameter obtained from the\n*          Gerschgorin intervals\n*\n*  CLGAPL  (input) DOUBLE PRECISION\n*\n*  CLGAPR  (input) DOUBLE PRECISION\n*          absolute gap on each end of the cluster.\n*          Set by the calling routine to protect against shifts too close\n*          to eigenvalues outside the cluster.\n*\n*  PIVMIN  (input) DOUBLE PRECISION\n*          The minimum pivot allowed in the Sturm sequence.\n*\n*  SIGMA   (output) DOUBLE PRECISION\n*          The shift used to form L(+) D(+) L(+)^T.\n*\n*  DPLUS   (output) DOUBLE PRECISION array, dimension (N)\n*          The N diagonal elements of the diagonal matrix D(+).\n*\n*  LPLUS   (output) DOUBLE PRECISION array, dimension (N-1)\n*          The first (N-1) elements of LPLUS contain the subdiagonal\n*          elements of the unit bidiagonal matrix L(+).\n*\n*  WORK    (workspace) DOUBLE PRECISION array, dimension (2*N)\n*          Workspace.\n*\n*  INFO    (output) INTEGER\n*          Signals processing OK (=0) or failure (=1)\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Beresford Parlett, University of California, Berkeley, USA\n*     Jim Demmel, University of California, Berkeley, USA\n*     Inderjit Dhillon, University of Texas, Austin, USA\n*     Osni Marques, LBNL/NERSC, USA\n*     Christof Voemel, University of California, Berkeley, USA\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  sigma, dplus, lplus, info, wgap = NumRu::Lapack.dlarrf( d, l, ld, clstrt, clend, w, wgap, werr, spdiam, clgapl, clgapr, pivmin, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 12)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 12)", argc);
  rblapack_d = argv[0];
  rblapack_l = argv[1];
  rblapack_ld = argv[2];
  rblapack_clstrt = argv[3];
  rblapack_clend = argv[4];
  rblapack_w = argv[5];
  rblapack_wgap = argv[6];
  rblapack_werr = argv[7];
  rblapack_spdiam = argv[8];
  rblapack_clgapl = argv[9];
  rblapack_clgapr = argv[10];
  rblapack_pivmin = argv[11];
  if (rb_options != Qnil) {
  }

  pivmin = NUM2DBL(rblapack_pivmin);
  clgapl = NUM2DBL(rblapack_clgapl);
  clend = NUM2INT(rblapack_clend);
  clgapr = NUM2DBL(rblapack_clgapr);
  spdiam = NUM2DBL(rblapack_spdiam);
  clstrt = NUM2INT(rblapack_clstrt);
  if (!NA_IsNArray(rblapack_d))
    rb_raise(rb_eArgError, "d (1th argument) must be NArray");
  if (NA_RANK(rblapack_d) != 1)
    rb_raise(rb_eArgError, "rank of d (1th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_d);
  if (NA_TYPE(rblapack_d) != NA_DFLOAT)
    rblapack_d = na_change_type(rblapack_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rblapack_d, doublereal*);
  if (!NA_IsNArray(rblapack_ld))
    rb_raise(rb_eArgError, "ld (3th argument) must be NArray");
  if (NA_RANK(rblapack_ld) != 1)
    rb_raise(rb_eArgError, "rank of ld (3th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_ld) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of ld must be %d", n-1);
  if (NA_TYPE(rblapack_ld) != NA_DFLOAT)
    rblapack_ld = na_change_type(rblapack_ld, NA_DFLOAT);
  ld = NA_PTR_TYPE(rblapack_ld, doublereal*);
  if (!NA_IsNArray(rblapack_werr))
    rb_raise(rb_eArgError, "werr (8th argument) must be NArray");
  if (NA_RANK(rblapack_werr) != 1)
    rb_raise(rb_eArgError, "rank of werr (8th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_werr) != (clend-clstrt+1))
    rb_raise(rb_eRuntimeError, "shape 0 of werr must be %d", clend-clstrt+1);
  if (NA_TYPE(rblapack_werr) != NA_DFLOAT)
    rblapack_werr = na_change_type(rblapack_werr, NA_DFLOAT);
  werr = NA_PTR_TYPE(rblapack_werr, doublereal*);
  if (!NA_IsNArray(rblapack_wgap))
    rb_raise(rb_eArgError, "wgap (7th argument) must be NArray");
  if (NA_RANK(rblapack_wgap) != 1)
    rb_raise(rb_eArgError, "rank of wgap (7th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_wgap) != (clend-clstrt+1))
    rb_raise(rb_eRuntimeError, "shape 0 of wgap must be %d", clend-clstrt+1);
  if (NA_TYPE(rblapack_wgap) != NA_DFLOAT)
    rblapack_wgap = na_change_type(rblapack_wgap, NA_DFLOAT);
  wgap = NA_PTR_TYPE(rblapack_wgap, doublereal*);
  if (!NA_IsNArray(rblapack_l))
    rb_raise(rb_eArgError, "l (2th argument) must be NArray");
  if (NA_RANK(rblapack_l) != 1)
    rb_raise(rb_eArgError, "rank of l (2th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_l) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of l must be %d", n-1);
  if (NA_TYPE(rblapack_l) != NA_DFLOAT)
    rblapack_l = na_change_type(rblapack_l, NA_DFLOAT);
  l = NA_PTR_TYPE(rblapack_l, doublereal*);
  if (!NA_IsNArray(rblapack_w))
    rb_raise(rb_eArgError, "w (6th argument) must be NArray");
  if (NA_RANK(rblapack_w) != 1)
    rb_raise(rb_eArgError, "rank of w (6th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_w) != (clend-clstrt+1))
    rb_raise(rb_eRuntimeError, "shape 0 of w must be %d", clend-clstrt+1);
  if (NA_TYPE(rblapack_w) != NA_DFLOAT)
    rblapack_w = na_change_type(rblapack_w, NA_DFLOAT);
  w = NA_PTR_TYPE(rblapack_w, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_dplus = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  dplus = NA_PTR_TYPE(rblapack_dplus, doublereal*);
  {
    int shape[1];
    shape[0] = n-1;
    rblapack_lplus = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  lplus = NA_PTR_TYPE(rblapack_lplus, doublereal*);
  {
    int shape[1];
    shape[0] = clend-clstrt+1;
    rblapack_wgap_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  wgap_out__ = NA_PTR_TYPE(rblapack_wgap_out__, doublereal*);
  MEMCPY(wgap_out__, wgap, doublereal, NA_TOTAL(rblapack_wgap));
  rblapack_wgap = rblapack_wgap_out__;
  wgap = wgap_out__;
  work = ALLOC_N(doublereal, (2*n));

  dlarrf_(&n, d, l, ld, &clstrt, &clend, w, wgap, werr, &spdiam, &clgapl, &clgapr, &pivmin, &sigma, dplus, lplus, work, &info);

  free(work);
  rblapack_sigma = rb_float_new((double)sigma);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(5, rblapack_sigma, rblapack_dplus, rblapack_lplus, rblapack_info, rblapack_wgap);
}

void
init_lapack_dlarrf(VALUE mLapack){
  rb_define_module_function(mLapack, "dlarrf", rblapack_dlarrf, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
