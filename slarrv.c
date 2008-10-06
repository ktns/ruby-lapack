#include "rb_lapack.h"

static VALUE
rb_slarrv(int argc, VALUE *argv, VALUE self){
  VALUE rb_vl;
  real vl; 
  VALUE rb_vu;
  real vu; 
  VALUE rb_d;
  real *d; 
  VALUE rb_l;
  real *l; 
  VALUE rb_pivmin;
  real pivmin; 
  VALUE rb_isplit;
  integer *isplit; 
  VALUE rb_m;
  integer m; 
  VALUE rb_dol;
  integer dol; 
  VALUE rb_dou;
  integer dou; 
  VALUE rb_minrgp;
  real minrgp; 
  VALUE rb_rtol1;
  real rtol1; 
  VALUE rb_rtol2;
  real rtol2; 
  VALUE rb_w;
  real *w; 
  VALUE rb_werr;
  real *werr; 
  VALUE rb_wgap;
  real *wgap; 
  VALUE rb_iblock;
  integer *iblock; 
  VALUE rb_indexw;
  integer *indexw; 
  VALUE rb_gers;
  real *gers; 
  VALUE rb_z;
  real *z; 
  VALUE rb_isuppz;
  integer *isuppz; 
  VALUE rb_info;
  integer info; 
  VALUE rb_d_out__;
  real *d_out__;
  VALUE rb_l_out__;
  real *l_out__;
  VALUE rb_w_out__;
  real *w_out__;
  VALUE rb_werr_out__;
  real *werr_out__;
  VALUE rb_wgap_out__;
  real *wgap_out__;
  real *work;
  integer *iwork;

  integer n;
  integer ldz;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  z, isuppz, info, d, l, w, werr, wgap = NumRu::Lapack.slarrv( vl, vu, d, l, pivmin, isplit, m, dol, dou, minrgp, rtol1, rtol2, w, werr, wgap, iblock, indexw, gers)\n    or\n  NumRu::Lapack.slarrv  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLARRV( N, VL, VU, D, L, PIVMIN, ISPLIT, M, DOL, DOU, MINRGP, RTOL1, RTOL2, W, WERR, WGAP, IBLOCK, INDEXW, GERS, Z, LDZ, ISUPPZ, WORK, IWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  SLARRV computes the eigenvectors of the tridiagonal matrix\n*  T = L D L^T given L, D and APPROXIMATIONS to the eigenvalues of L D L^T.\n*  The input eigenvalues should have been computed by SLARRE.\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The order of the matrix.  N >= 0.\n*\n*  VL      (input) REAL            \n*  VU      (input) REAL            \n*          Lower and upper bounds of the interval that contains the desired\n*          eigenvalues. VL < VU. Needed to compute gaps on the left or right\n*          end of the extremal eigenvalues in the desired RANGE.\n*\n*  D       (input/output) REAL             array, dimension (N)\n*          On entry, the N diagonal elements of the diagonal matrix D.\n*          On exit, D may be overwritten.\n*\n*  L       (input/output) REAL             array, dimension (N)\n*          On entry, the (N-1) subdiagonal elements of the unit\n*          bidiagonal matrix L are in elements 1 to N-1 of L\n*          (if the matrix is not splitted.) At the end of each block\n*          is stored the corresponding shift as given by SLARRE.\n*          On exit, L is overwritten.\n*\n*  PIVMIN  (in) DOUBLE PRECISION\n*          The minimum pivot allowed in the Sturm sequence.\n*\n*  ISPLIT  (input) INTEGER array, dimension (N)\n*          The splitting points, at which T breaks up into blocks.\n*          The first block consists of rows/columns 1 to\n*          ISPLIT( 1 ), the second of rows/columns ISPLIT( 1 )+1\n*          through ISPLIT( 2 ), etc.\n*\n*  M       (input) INTEGER\n*          The total number of input eigenvalues.  0 <= M <= N.\n*\n*  DOL     (input) INTEGER\n*  DOU     (input) INTEGER\n*          If the user wants to compute only selected eigenvectors from all\n*          the eigenvalues supplied, he can specify an index range DOL:DOU.\n*          Or else the setting DOL=1, DOU=M should be applied.\n*          Note that DOL and DOU refer to the order in which the eigenvalues\n*          are stored in W.\n*          If the user wants to compute only selected eigenpairs, then\n*          the columns DOL-1 to DOU+1 of the eigenvector space Z contain the\n*          computed eigenvectors. All other columns of Z are set to zero.\n*\n*  MINRGP  (input) REAL            \n*\n*  RTOL1   (input) REAL            \n*  RTOL2   (input) REAL            \n*           Parameters for bisection.\n*           An interval [LEFT,RIGHT] has converged if\n*           RIGHT-LEFT.LT.MAX( RTOL1*GAP, RTOL2*MAX(|LEFT|,|RIGHT|) )\n*\n*  W       (input/output) REAL             array, dimension (N)\n*          The first M elements of W contain the APPROXIMATE eigenvalues for\n*          which eigenvectors are to be computed.  The eigenvalues\n*          should be grouped by split-off block and ordered from\n*          smallest to largest within the block ( The output array\n*          W from SLARRE is expected here ). Furthermore, they are with\n*          respect to the shift of the corresponding root representation\n*          for their block. On exit, W holds the eigenvalues of the\n*          UNshifted matrix.\n*\n*  WERR    (input/output) REAL             array, dimension (N)\n*          The first M elements contain the semiwidth of the uncertainty\n*          interval of the corresponding eigenvalue in W\n*\n*  WGAP    (input/output) REAL             array, dimension (N)\n*          The separation from the right neighbor eigenvalue in W.\n*\n*  IBLOCK  (input) INTEGER array, dimension (N)\n*          The indices of the blocks (submatrices) associated with the\n*          corresponding eigenvalues in W; IBLOCK(i)=1 if eigenvalue\n*          W(i) belongs to the first block from the top, =2 if W(i)\n*          belongs to the second block, etc.\n*\n*  INDEXW  (input) INTEGER array, dimension (N)\n*          The indices of the eigenvalues within each block (submatrix);\n*          for example, INDEXW(i)= 10 and IBLOCK(i)=2 imply that the\n*          i-th eigenvalue W(i) is the 10-th eigenvalue in the second block.\n*\n*  GERS    (input) REAL             array, dimension (2*N)\n*          The N Gerschgorin intervals (the i-th Gerschgorin interval\n*          is (GERS(2*i-1), GERS(2*i)). The Gerschgorin intervals should\n*          be computed from the original UNshifted matrix.\n*\n*  Z       (output) REAL             array, dimension (LDZ, max(1,M) )\n*          If INFO = 0, the first M columns of Z contain the\n*          orthonormal eigenvectors of the matrix T\n*          corresponding to the input eigenvalues, with the i-th\n*          column of Z holding the eigenvector associated with W(i).\n*          Note: the user must ensure that at least max(1,M) columns are\n*          supplied in the array Z.\n*\n*  LDZ     (input) INTEGER\n*          The leading dimension of the array Z.  LDZ >= 1, and if\n*          JOBZ = 'V', LDZ >= max(1,N).\n*\n*  ISUPPZ  (output) INTEGER array, dimension ( 2*max(1,M) )\n*          The support of the eigenvectors in Z, i.e., the indices\n*          indicating the nonzero elements in Z. The I-th eigenvector\n*          is nonzero only in elements ISUPPZ( 2*I-1 ) through\n*          ISUPPZ( 2*I ).\n*\n*  WORK    (workspace) REAL             array, dimension (12*N)\n*\n*  IWORK   (workspace) INTEGER array, dimension (7*N)\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*\n*          > 0:  A problem occured in SLARRV.\n*          < 0:  One of the called subroutines signaled an internal problem.\n*                Needs inspection of the corresponding parameter IINFO\n*                for further information.\n*\n*          =-1:  Problem in SLARRB when refining a child's eigenvalues.\n*          =-2:  Problem in SLARRF when computing the RRR of a child.\n*                When a child is inside a tight cluster, it can be difficult\n*                to find an RRR. A partial remedy from the user's point of\n*                view is to make the parameter MINRGP smaller and recompile.\n*                However, as the orthogonality of the computed vectors is\n*                proportional to 1/MINRGP, the user should be aware that\n*                he might be trading in precision when he decreases MINRGP.\n*          =-3:  Problem in SLARRB when refining a single eigenvalue\n*                after the Rayleigh correction was rejected.\n*          = 5:  The Rayleigh Quotient Iteration failed to converge to\n*                full accuracy in MAXITR steps.\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Beresford Parlett, University of California, Berkeley, USA\n*     Jim Demmel, University of California, Berkeley, USA\n*     Inderjit Dhillon, University of Texas, Austin, USA\n*     Osni Marques, LBNL/NERSC, USA\n*     Christof Voemel, University of California, Berkeley, USA\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 18)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 18)", argc);
  rb_vl = argv[0];
  rb_vu = argv[1];
  rb_d = argv[2];
  rb_l = argv[3];
  rb_pivmin = argv[4];
  rb_isplit = argv[5];
  rb_m = argv[6];
  rb_dol = argv[7];
  rb_dou = argv[8];
  rb_minrgp = argv[9];
  rb_rtol1 = argv[10];
  rb_rtol2 = argv[11];
  rb_w = argv[12];
  rb_werr = argv[13];
  rb_wgap = argv[14];
  rb_iblock = argv[15];
  rb_indexw = argv[16];
  rb_gers = argv[17];

  vl = (real)NUM2DBL(rb_vl);
  vu = (real)NUM2DBL(rb_vu);
  pivmin = (real)NUM2DBL(rb_pivmin);
  m = NUM2INT(rb_m);
  dol = NUM2INT(rb_dol);
  dou = NUM2INT(rb_dou);
  minrgp = (real)NUM2DBL(rb_minrgp);
  rtol1 = (real)NUM2DBL(rb_rtol1);
  rtol2 = (real)NUM2DBL(rb_rtol2);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (3th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (3th argument) must be %d", 1);
  n = NA_SHAPE0(rb_d);
  if (NA_TYPE(rb_d) != NA_SFLOAT)
    rb_d = na_change_type(rb_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rb_d, real*);
  if (!NA_IsNArray(rb_l))
    rb_raise(rb_eArgError, "l (4th argument) must be NArray");
  if (NA_RANK(rb_l) != 1)
    rb_raise(rb_eArgError, "rank of l (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_l) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of l must be the same as shape 0 of d");
  if (NA_TYPE(rb_l) != NA_SFLOAT)
    rb_l = na_change_type(rb_l, NA_SFLOAT);
  l = NA_PTR_TYPE(rb_l, real*);
  if (!NA_IsNArray(rb_isplit))
    rb_raise(rb_eArgError, "isplit (6th argument) must be NArray");
  if (NA_RANK(rb_isplit) != 1)
    rb_raise(rb_eArgError, "rank of isplit (6th argument) must be %d", 1);
  if (NA_SHAPE0(rb_isplit) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of isplit must be the same as shape 0 of d");
  if (NA_TYPE(rb_isplit) != NA_LINT)
    rb_isplit = na_change_type(rb_isplit, NA_LINT);
  isplit = NA_PTR_TYPE(rb_isplit, integer*);
  if (!NA_IsNArray(rb_w))
    rb_raise(rb_eArgError, "w (13th argument) must be NArray");
  if (NA_RANK(rb_w) != 1)
    rb_raise(rb_eArgError, "rank of w (13th argument) must be %d", 1);
  if (NA_SHAPE0(rb_w) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of w must be the same as shape 0 of d");
  if (NA_TYPE(rb_w) != NA_SFLOAT)
    rb_w = na_change_type(rb_w, NA_SFLOAT);
  w = NA_PTR_TYPE(rb_w, real*);
  if (!NA_IsNArray(rb_werr))
    rb_raise(rb_eArgError, "werr (14th argument) must be NArray");
  if (NA_RANK(rb_werr) != 1)
    rb_raise(rb_eArgError, "rank of werr (14th argument) must be %d", 1);
  if (NA_SHAPE0(rb_werr) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of werr must be the same as shape 0 of d");
  if (NA_TYPE(rb_werr) != NA_SFLOAT)
    rb_werr = na_change_type(rb_werr, NA_SFLOAT);
  werr = NA_PTR_TYPE(rb_werr, real*);
  if (!NA_IsNArray(rb_wgap))
    rb_raise(rb_eArgError, "wgap (15th argument) must be NArray");
  if (NA_RANK(rb_wgap) != 1)
    rb_raise(rb_eArgError, "rank of wgap (15th argument) must be %d", 1);
  if (NA_SHAPE0(rb_wgap) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of wgap must be the same as shape 0 of d");
  if (NA_TYPE(rb_wgap) != NA_SFLOAT)
    rb_wgap = na_change_type(rb_wgap, NA_SFLOAT);
  wgap = NA_PTR_TYPE(rb_wgap, real*);
  if (!NA_IsNArray(rb_iblock))
    rb_raise(rb_eArgError, "iblock (16th argument) must be NArray");
  if (NA_RANK(rb_iblock) != 1)
    rb_raise(rb_eArgError, "rank of iblock (16th argument) must be %d", 1);
  if (NA_SHAPE0(rb_iblock) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of iblock must be the same as shape 0 of d");
  if (NA_TYPE(rb_iblock) != NA_LINT)
    rb_iblock = na_change_type(rb_iblock, NA_LINT);
  iblock = NA_PTR_TYPE(rb_iblock, integer*);
  if (!NA_IsNArray(rb_indexw))
    rb_raise(rb_eArgError, "indexw (17th argument) must be NArray");
  if (NA_RANK(rb_indexw) != 1)
    rb_raise(rb_eArgError, "rank of indexw (17th argument) must be %d", 1);
  if (NA_SHAPE0(rb_indexw) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of indexw must be the same as shape 0 of d");
  if (NA_TYPE(rb_indexw) != NA_LINT)
    rb_indexw = na_change_type(rb_indexw, NA_LINT);
  indexw = NA_PTR_TYPE(rb_indexw, integer*);
  if (!NA_IsNArray(rb_gers))
    rb_raise(rb_eArgError, "gers (18th argument) must be NArray");
  if (NA_RANK(rb_gers) != 1)
    rb_raise(rb_eArgError, "rank of gers (18th argument) must be %d", 1);
  if (NA_SHAPE0(rb_gers) != (2*n))
    rb_raise(rb_eRuntimeError, "shape 0 of gers must be %d", 2*n);
  if (NA_TYPE(rb_gers) != NA_SFLOAT)
    rb_gers = na_change_type(rb_gers, NA_SFLOAT);
  gers = NA_PTR_TYPE(rb_gers, real*);
  ldz = n;
  {
    int shape[2];
    shape[0] = ldz;
    shape[1] = MAX(1,m);
    rb_z = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  z = NA_PTR_TYPE(rb_z, real*);
  {
    int shape[1];
    shape[0] = 2*MAX(1,m);
    rb_isuppz = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  isuppz = NA_PTR_TYPE(rb_isuppz, integer*);
  {
    int shape[1];
    shape[0] = n;
    rb_d_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  d_out__ = NA_PTR_TYPE(rb_d_out__, real*);
  MEMCPY(d_out__, d, real, NA_TOTAL(rb_d));
  rb_d = rb_d_out__;
  d = d_out__;
  {
    int shape[1];
    shape[0] = n;
    rb_l_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  l_out__ = NA_PTR_TYPE(rb_l_out__, real*);
  MEMCPY(l_out__, l, real, NA_TOTAL(rb_l));
  rb_l = rb_l_out__;
  l = l_out__;
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
  {
    int shape[1];
    shape[0] = n;
    rb_wgap_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  wgap_out__ = NA_PTR_TYPE(rb_wgap_out__, real*);
  MEMCPY(wgap_out__, wgap, real, NA_TOTAL(rb_wgap));
  rb_wgap = rb_wgap_out__;
  wgap = wgap_out__;
  work = ALLOC_N(real, (12*n));
  iwork = ALLOC_N(integer, (7*n));

  slarrv_(&n, &vl, &vu, d, l, &pivmin, isplit, &m, &dol, &dou, &minrgp, &rtol1, &rtol2, w, werr, wgap, iblock, indexw, gers, z, &ldz, isuppz, work, iwork, &info);

  free(work);
  free(iwork);
  rb_info = INT2NUM(info);
  return rb_ary_new3(8, rb_z, rb_isuppz, rb_info, rb_d, rb_l, rb_w, rb_werr, rb_wgap);
}

void
init_lapack_slarrv(VALUE mLapack){
  rb_define_module_function(mLapack, "slarrv", rb_slarrv, -1);
}
