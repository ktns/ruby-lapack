#include "rb_lapack.h"

extern VOID dlarrv_(integer *n, doublereal *vl, doublereal *vu, doublereal *d, doublereal *l, doublereal *pivmin, integer *isplit, integer *m, integer *dol, integer *dou, doublereal *minrgp, doublereal *rtol1, doublereal *rtol2, doublereal *w, doublereal *werr, doublereal *wgap, integer *iblock, integer *indexw, doublereal *gers, doublereal *z, integer *ldz, integer *isuppz, doublereal *work, integer *iwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dlarrv(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_vl;
  doublereal vl; 
  VALUE rblapack_vu;
  doublereal vu; 
  VALUE rblapack_d;
  doublereal *d; 
  VALUE rblapack_l;
  doublereal *l; 
  VALUE rblapack_pivmin;
  doublereal pivmin; 
  VALUE rblapack_isplit;
  integer *isplit; 
  VALUE rblapack_m;
  integer m; 
  VALUE rblapack_dol;
  integer dol; 
  VALUE rblapack_dou;
  integer dou; 
  VALUE rblapack_minrgp;
  doublereal minrgp; 
  VALUE rblapack_rtol1;
  doublereal rtol1; 
  VALUE rblapack_rtol2;
  doublereal rtol2; 
  VALUE rblapack_w;
  doublereal *w; 
  VALUE rblapack_werr;
  doublereal *werr; 
  VALUE rblapack_wgap;
  doublereal *wgap; 
  VALUE rblapack_iblock;
  integer *iblock; 
  VALUE rblapack_indexw;
  integer *indexw; 
  VALUE rblapack_gers;
  doublereal *gers; 
  VALUE rblapack_z;
  doublereal *z; 
  VALUE rblapack_isuppz;
  integer *isuppz; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_d_out__;
  doublereal *d_out__;
  VALUE rblapack_l_out__;
  doublereal *l_out__;
  VALUE rblapack_w_out__;
  doublereal *w_out__;
  VALUE rblapack_werr_out__;
  doublereal *werr_out__;
  VALUE rblapack_wgap_out__;
  doublereal *wgap_out__;
  doublereal *work;
  integer *iwork;

  integer n;
  integer ldz;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  z, isuppz, info, d, l, w, werr, wgap = NumRu::Lapack.dlarrv( vl, vu, d, l, pivmin, isplit, m, dol, dou, minrgp, rtol1, rtol2, w, werr, wgap, iblock, indexw, gers, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLARRV( N, VL, VU, D, L, PIVMIN, ISPLIT, M, DOL, DOU, MINRGP, RTOL1, RTOL2, W, WERR, WGAP, IBLOCK, INDEXW, GERS, Z, LDZ, ISUPPZ, WORK, IWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  DLARRV computes the eigenvectors of the tridiagonal matrix\n*  T = L D L^T given L, D and APPROXIMATIONS to the eigenvalues of L D L^T.\n*  The input eigenvalues should have been computed by DLARRE.\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The order of the matrix.  N >= 0.\n*\n*  VL      (input) DOUBLE PRECISION\n*  VU      (input) DOUBLE PRECISION\n*          Lower and upper bounds of the interval that contains the desired\n*          eigenvalues. VL < VU. Needed to compute gaps on the left or right\n*          end of the extremal eigenvalues in the desired RANGE.\n*\n*  D       (input/output) DOUBLE PRECISION array, dimension (N)\n*          On entry, the N diagonal elements of the diagonal matrix D.\n*          On exit, D may be overwritten.\n*\n*  L       (input/output) DOUBLE PRECISION array, dimension (N)\n*          On entry, the (N-1) subdiagonal elements of the unit\n*          bidiagonal matrix L are in elements 1 to N-1 of L\n*          (if the matrix is not splitted.) At the end of each block\n*          is stored the corresponding shift as given by DLARRE.\n*          On exit, L is overwritten.\n*\n*  PIVMIN  (input) DOUBLE PRECISION\n*          The minimum pivot allowed in the Sturm sequence.\n*\n*  ISPLIT  (input) INTEGER array, dimension (N)\n*          The splitting points, at which T breaks up into blocks.\n*          The first block consists of rows/columns 1 to\n*          ISPLIT( 1 ), the second of rows/columns ISPLIT( 1 )+1\n*          through ISPLIT( 2 ), etc.\n*\n*  M       (input) INTEGER\n*          The total number of input eigenvalues.  0 <= M <= N.\n*\n*  DOL     (input) INTEGER\n*  DOU     (input) INTEGER\n*          If the user wants to compute only selected eigenvectors from all\n*          the eigenvalues supplied, he can specify an index range DOL:DOU.\n*          Or else the setting DOL=1, DOU=M should be applied.\n*          Note that DOL and DOU refer to the order in which the eigenvalues\n*          are stored in W.\n*          If the user wants to compute only selected eigenpairs, then\n*          the columns DOL-1 to DOU+1 of the eigenvector space Z contain the\n*          computed eigenvectors. All other columns of Z are set to zero.\n*\n*  MINRGP  (input) DOUBLE PRECISION\n*\n*  RTOL1   (input) DOUBLE PRECISION\n*  RTOL2   (input) DOUBLE PRECISION\n*           Parameters for bisection.\n*           An interval [LEFT,RIGHT] has converged if\n*           RIGHT-LEFT.LT.MAX( RTOL1*GAP, RTOL2*MAX(|LEFT|,|RIGHT|) )\n*\n*  W       (input/output) DOUBLE PRECISION array, dimension (N)\n*          The first M elements of W contain the APPROXIMATE eigenvalues for\n*          which eigenvectors are to be computed.  The eigenvalues\n*          should be grouped by split-off block and ordered from\n*          smallest to largest within the block ( The output array\n*          W from DLARRE is expected here ). Furthermore, they are with\n*          respect to the shift of the corresponding root representation\n*          for their block. On exit, W holds the eigenvalues of the\n*          UNshifted matrix.\n*\n*  WERR    (input/output) DOUBLE PRECISION array, dimension (N)\n*          The first M elements contain the semiwidth of the uncertainty\n*          interval of the corresponding eigenvalue in W\n*\n*  WGAP    (input/output) DOUBLE PRECISION array, dimension (N)\n*          The separation from the right neighbor eigenvalue in W.\n*\n*  IBLOCK  (input) INTEGER array, dimension (N)\n*          The indices of the blocks (submatrices) associated with the\n*          corresponding eigenvalues in W; IBLOCK(i)=1 if eigenvalue\n*          W(i) belongs to the first block from the top, =2 if W(i)\n*          belongs to the second block, etc.\n*\n*  INDEXW  (input) INTEGER array, dimension (N)\n*          The indices of the eigenvalues within each block (submatrix);\n*          for example, INDEXW(i)= 10 and IBLOCK(i)=2 imply that the\n*          i-th eigenvalue W(i) is the 10-th eigenvalue in the second block.\n*\n*  GERS    (input) DOUBLE PRECISION array, dimension (2*N)\n*          The N Gerschgorin intervals (the i-th Gerschgorin interval\n*          is (GERS(2*i-1), GERS(2*i)). The Gerschgorin intervals should\n*          be computed from the original UNshifted matrix.\n*\n*  Z       (output) DOUBLE PRECISION array, dimension (LDZ, max(1,M) )\n*          If INFO = 0, the first M columns of Z contain the\n*          orthonormal eigenvectors of the matrix T\n*          corresponding to the input eigenvalues, with the i-th\n*          column of Z holding the eigenvector associated with W(i).\n*          Note: the user must ensure that at least max(1,M) columns are\n*          supplied in the array Z.\n*\n*  LDZ     (input) INTEGER\n*          The leading dimension of the array Z.  LDZ >= 1, and if\n*          JOBZ = 'V', LDZ >= max(1,N).\n*\n*  ISUPPZ  (output) INTEGER array, dimension ( 2*max(1,M) )\n*          The support of the eigenvectors in Z, i.e., the indices\n*          indicating the nonzero elements in Z. The I-th eigenvector\n*          is nonzero only in elements ISUPPZ( 2*I-1 ) through\n*          ISUPPZ( 2*I ).\n*\n*  WORK    (workspace) DOUBLE PRECISION array, dimension (12*N)\n*\n*  IWORK   (workspace) INTEGER array, dimension (7*N)\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*\n*          > 0:  A problem occured in DLARRV.\n*          < 0:  One of the called subroutines signaled an internal problem.\n*                Needs inspection of the corresponding parameter IINFO\n*                for further information.\n*\n*          =-1:  Problem in DLARRB when refining a child's eigenvalues.\n*          =-2:  Problem in DLARRF when computing the RRR of a child.\n*                When a child is inside a tight cluster, it can be difficult\n*                to find an RRR. A partial remedy from the user's point of\n*                view is to make the parameter MINRGP smaller and recompile.\n*                However, as the orthogonality of the computed vectors is\n*                proportional to 1/MINRGP, the user should be aware that\n*                he might be trading in precision when he decreases MINRGP.\n*          =-3:  Problem in DLARRB when refining a single eigenvalue\n*                after the Rayleigh correction was rejected.\n*          = 5:  The Rayleigh Quotient Iteration failed to converge to\n*                full accuracy in MAXITR steps.\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Beresford Parlett, University of California, Berkeley, USA\n*     Jim Demmel, University of California, Berkeley, USA\n*     Inderjit Dhillon, University of Texas, Austin, USA\n*     Osni Marques, LBNL/NERSC, USA\n*     Christof Voemel, University of California, Berkeley, USA\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  z, isuppz, info, d, l, w, werr, wgap = NumRu::Lapack.dlarrv( vl, vu, d, l, pivmin, isplit, m, dol, dou, minrgp, rtol1, rtol2, w, werr, wgap, iblock, indexw, gers, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 18)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 18)", argc);
  rblapack_vl = argv[0];
  rblapack_vu = argv[1];
  rblapack_d = argv[2];
  rblapack_l = argv[3];
  rblapack_pivmin = argv[4];
  rblapack_isplit = argv[5];
  rblapack_m = argv[6];
  rblapack_dol = argv[7];
  rblapack_dou = argv[8];
  rblapack_minrgp = argv[9];
  rblapack_rtol1 = argv[10];
  rblapack_rtol2 = argv[11];
  rblapack_w = argv[12];
  rblapack_werr = argv[13];
  rblapack_wgap = argv[14];
  rblapack_iblock = argv[15];
  rblapack_indexw = argv[16];
  rblapack_gers = argv[17];
  if (rb_options != Qnil) {
  }

  vl = NUM2DBL(rblapack_vl);
  if (!NA_IsNArray(rblapack_w))
    rb_raise(rb_eArgError, "w (13th argument) must be NArray");
  if (NA_RANK(rblapack_w) != 1)
    rb_raise(rb_eArgError, "rank of w (13th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_w);
  if (NA_TYPE(rblapack_w) != NA_DFLOAT)
    rblapack_w = na_change_type(rblapack_w, NA_DFLOAT);
  w = NA_PTR_TYPE(rblapack_w, doublereal*);
  dol = NUM2INT(rblapack_dol);
  if (!NA_IsNArray(rblapack_l))
    rb_raise(rb_eArgError, "l (4th argument) must be NArray");
  if (NA_RANK(rblapack_l) != 1)
    rb_raise(rb_eArgError, "rank of l (4th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_l) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of l must be the same as shape 0 of w");
  if (NA_TYPE(rblapack_l) != NA_DFLOAT)
    rblapack_l = na_change_type(rblapack_l, NA_DFLOAT);
  l = NA_PTR_TYPE(rblapack_l, doublereal*);
  pivmin = NUM2DBL(rblapack_pivmin);
  dou = NUM2INT(rblapack_dou);
  if (!NA_IsNArray(rblapack_wgap))
    rb_raise(rb_eArgError, "wgap (15th argument) must be NArray");
  if (NA_RANK(rblapack_wgap) != 1)
    rb_raise(rb_eArgError, "rank of wgap (15th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_wgap) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of wgap must be the same as shape 0 of w");
  if (NA_TYPE(rblapack_wgap) != NA_DFLOAT)
    rblapack_wgap = na_change_type(rblapack_wgap, NA_DFLOAT);
  wgap = NA_PTR_TYPE(rblapack_wgap, doublereal*);
  m = NUM2INT(rblapack_m);
  minrgp = NUM2DBL(rblapack_minrgp);
  rtol2 = NUM2DBL(rblapack_rtol2);
  if (!NA_IsNArray(rblapack_isplit))
    rb_raise(rb_eArgError, "isplit (6th argument) must be NArray");
  if (NA_RANK(rblapack_isplit) != 1)
    rb_raise(rb_eArgError, "rank of isplit (6th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_isplit) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of isplit must be the same as shape 0 of w");
  if (NA_TYPE(rblapack_isplit) != NA_LINT)
    rblapack_isplit = na_change_type(rblapack_isplit, NA_LINT);
  isplit = NA_PTR_TYPE(rblapack_isplit, integer*);
  if (!NA_IsNArray(rblapack_indexw))
    rb_raise(rb_eArgError, "indexw (17th argument) must be NArray");
  if (NA_RANK(rblapack_indexw) != 1)
    rb_raise(rb_eArgError, "rank of indexw (17th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_indexw) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of indexw must be the same as shape 0 of w");
  if (NA_TYPE(rblapack_indexw) != NA_LINT)
    rblapack_indexw = na_change_type(rblapack_indexw, NA_LINT);
  indexw = NA_PTR_TYPE(rblapack_indexw, integer*);
  if (!NA_IsNArray(rblapack_werr))
    rb_raise(rb_eArgError, "werr (14th argument) must be NArray");
  if (NA_RANK(rblapack_werr) != 1)
    rb_raise(rb_eArgError, "rank of werr (14th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_werr) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of werr must be the same as shape 0 of w");
  if (NA_TYPE(rblapack_werr) != NA_DFLOAT)
    rblapack_werr = na_change_type(rblapack_werr, NA_DFLOAT);
  werr = NA_PTR_TYPE(rblapack_werr, doublereal*);
  if (!NA_IsNArray(rblapack_iblock))
    rb_raise(rb_eArgError, "iblock (16th argument) must be NArray");
  if (NA_RANK(rblapack_iblock) != 1)
    rb_raise(rb_eArgError, "rank of iblock (16th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_iblock) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of iblock must be the same as shape 0 of w");
  if (NA_TYPE(rblapack_iblock) != NA_LINT)
    rblapack_iblock = na_change_type(rblapack_iblock, NA_LINT);
  iblock = NA_PTR_TYPE(rblapack_iblock, integer*);
  rtol1 = NUM2DBL(rblapack_rtol1);
  vu = NUM2DBL(rblapack_vu);
  if (!NA_IsNArray(rblapack_d))
    rb_raise(rb_eArgError, "d (3th argument) must be NArray");
  if (NA_RANK(rblapack_d) != 1)
    rb_raise(rb_eArgError, "rank of d (3th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_d) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of d must be the same as shape 0 of w");
  if (NA_TYPE(rblapack_d) != NA_DFLOAT)
    rblapack_d = na_change_type(rblapack_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rblapack_d, doublereal*);
  if (!NA_IsNArray(rblapack_gers))
    rb_raise(rb_eArgError, "gers (18th argument) must be NArray");
  if (NA_RANK(rblapack_gers) != 1)
    rb_raise(rb_eArgError, "rank of gers (18th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_gers) != (2*n))
    rb_raise(rb_eRuntimeError, "shape 0 of gers must be %d", 2*n);
  if (NA_TYPE(rblapack_gers) != NA_DFLOAT)
    rblapack_gers = na_change_type(rblapack_gers, NA_DFLOAT);
  gers = NA_PTR_TYPE(rblapack_gers, doublereal*);
  ldz = n;
  {
    int shape[2];
    shape[0] = ldz;
    shape[1] = MAX(1,m);
    rblapack_z = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  z = NA_PTR_TYPE(rblapack_z, doublereal*);
  {
    int shape[1];
    shape[0] = 2*MAX(1,m);
    rblapack_isuppz = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  isuppz = NA_PTR_TYPE(rblapack_isuppz, integer*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_d_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  d_out__ = NA_PTR_TYPE(rblapack_d_out__, doublereal*);
  MEMCPY(d_out__, d, doublereal, NA_TOTAL(rblapack_d));
  rblapack_d = rblapack_d_out__;
  d = d_out__;
  {
    int shape[1];
    shape[0] = n;
    rblapack_l_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  l_out__ = NA_PTR_TYPE(rblapack_l_out__, doublereal*);
  MEMCPY(l_out__, l, doublereal, NA_TOTAL(rblapack_l));
  rblapack_l = rblapack_l_out__;
  l = l_out__;
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
  {
    int shape[1];
    shape[0] = n;
    rblapack_wgap_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  wgap_out__ = NA_PTR_TYPE(rblapack_wgap_out__, doublereal*);
  MEMCPY(wgap_out__, wgap, doublereal, NA_TOTAL(rblapack_wgap));
  rblapack_wgap = rblapack_wgap_out__;
  wgap = wgap_out__;
  work = ALLOC_N(doublereal, (12*n));
  iwork = ALLOC_N(integer, (7*n));

  dlarrv_(&n, &vl, &vu, d, l, &pivmin, isplit, &m, &dol, &dou, &minrgp, &rtol1, &rtol2, w, werr, wgap, iblock, indexw, gers, z, &ldz, isuppz, work, iwork, &info);

  free(work);
  free(iwork);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(8, rblapack_z, rblapack_isuppz, rblapack_info, rblapack_d, rblapack_l, rblapack_w, rblapack_werr, rblapack_wgap);
}

void
init_lapack_dlarrv(VALUE mLapack){
  rb_define_module_function(mLapack, "dlarrv", rblapack_dlarrv, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
