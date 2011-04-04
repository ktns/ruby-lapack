#include "rb_lapack.h"

extern VOID dlarre_(char *range, integer *n, doublereal *vl, doublereal *vu, integer *il, integer *iu, doublereal *d, doublereal *e, doublereal *e2, doublereal *rtol1, doublereal *rtol2, doublereal *spltol, integer *nsplit, integer *isplit, integer *m, doublereal *w, doublereal *werr, doublereal *wgap, integer *iblock, integer *indexw, doublereal *gers, doublereal *pivmin, doublereal *work, integer *iwork, integer *info);

static VALUE
rb_dlarre(int argc, VALUE *argv, VALUE self){
  VALUE rb_range;
  char range; 
  VALUE rb_vl;
  doublereal vl; 
  VALUE rb_vu;
  doublereal vu; 
  VALUE rb_il;
  integer il; 
  VALUE rb_iu;
  integer iu; 
  VALUE rb_d;
  doublereal *d; 
  VALUE rb_e;
  doublereal *e; 
  VALUE rb_e2;
  doublereal *e2; 
  VALUE rb_rtol1;
  doublereal rtol1; 
  VALUE rb_rtol2;
  doublereal rtol2; 
  VALUE rb_spltol;
  doublereal spltol; 
  VALUE rb_nsplit;
  integer nsplit; 
  VALUE rb_isplit;
  integer *isplit; 
  VALUE rb_m;
  integer m; 
  VALUE rb_w;
  doublereal *w; 
  VALUE rb_werr;
  doublereal *werr; 
  VALUE rb_wgap;
  doublereal *wgap; 
  VALUE rb_iblock;
  integer *iblock; 
  VALUE rb_indexw;
  integer *indexw; 
  VALUE rb_gers;
  doublereal *gers; 
  VALUE rb_pivmin;
  doublereal pivmin; 
  VALUE rb_info;
  integer info; 
  VALUE rb_d_out__;
  doublereal *d_out__;
  VALUE rb_e_out__;
  doublereal *e_out__;
  VALUE rb_e2_out__;
  doublereal *e2_out__;
  doublereal *work;
  integer *iwork;

  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  nsplit, isplit, m, w, werr, wgap, iblock, indexw, gers, pivmin, info, vl, vu, d, e, e2 = NumRu::Lapack.dlarre( range, vl, vu, il, iu, d, e, e2, rtol1, rtol2, spltol)\n    or\n  NumRu::Lapack.dlarre  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLARRE( RANGE, N, VL, VU, IL, IU, D, E, E2, RTOL1, RTOL2, SPLTOL, NSPLIT, ISPLIT, M, W, WERR, WGAP, IBLOCK, INDEXW, GERS, PIVMIN, WORK, IWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  To find the desired eigenvalues of a given real symmetric\n*  tridiagonal matrix T, DLARRE sets any \"small\" off-diagonal\n*  elements to zero, and for each unreduced block T_i, it finds\n*  (a) a suitable shift at one end of the block's spectrum,\n*  (b) the base representation, T_i - sigma_i I = L_i D_i L_i^T, and\n*  (c) eigenvalues of each L_i D_i L_i^T.\n*  The representations and eigenvalues found are then used by\n*  DSTEMR to compute the eigenvectors of T.\n*  The accuracy varies depending on whether bisection is used to\n*  find a few eigenvalues or the dqds algorithm (subroutine DLASQ2) to\n*  conpute all and then discard any unwanted one.\n*  As an added benefit, DLARRE also outputs the n\n*  Gerschgorin intervals for the matrices L_i D_i L_i^T.\n*\n\n*  Arguments\n*  =========\n*\n*  RANGE   (input) CHARACTER*1\n*          = 'A': (\"All\")   all eigenvalues will be found.\n*          = 'V': (\"Value\") all eigenvalues in the half-open interval\n*                           (VL, VU] will be found.\n*          = 'I': (\"Index\") the IL-th through IU-th eigenvalues (of the\n*                           entire matrix) will be found.\n*\n*  N       (input) INTEGER\n*          The order of the matrix. N > 0.\n*\n*  VL      (input/output) DOUBLE PRECISION\n*  VU      (input/output) DOUBLE PRECISION\n*          If RANGE='V', the lower and upper bounds for the eigenvalues.\n*          Eigenvalues less than or equal to VL, or greater than VU,\n*          will not be returned.  VL < VU.\n*          If RANGE='I' or ='A', DLARRE computes bounds on the desired\n*          part of the spectrum.\n*\n*  IL      (input) INTEGER\n*  IU      (input) INTEGER\n*          If RANGE='I', the indices (in ascending order) of the\n*          smallest and largest eigenvalues to be returned.\n*          1 <= IL <= IU <= N.\n*\n*  D       (input/output) DOUBLE PRECISION array, dimension (N)\n*          On entry, the N diagonal elements of the tridiagonal\n*          matrix T.\n*          On exit, the N diagonal elements of the diagonal\n*          matrices D_i.\n*\n*  E       (input/output) DOUBLE PRECISION array, dimension (N)\n*          On entry, the first (N-1) entries contain the subdiagonal\n*          elements of the tridiagonal matrix T; E(N) need not be set.\n*          On exit, E contains the subdiagonal elements of the unit\n*          bidiagonal matrices L_i. The entries E( ISPLIT( I ) ),\n*          1 <= I <= NSPLIT, contain the base points sigma_i on output.\n*\n*  E2      (input/output) DOUBLE PRECISION array, dimension (N)\n*          On entry, the first (N-1) entries contain the SQUARES of the\n*          subdiagonal elements of the tridiagonal matrix T;\n*          E2(N) need not be set.\n*          On exit, the entries E2( ISPLIT( I ) ),\n*          1 <= I <= NSPLIT, have been set to zero\n*\n*  RTOL1   (input) DOUBLE PRECISION\n*  RTOL2   (input) DOUBLE PRECISION\n*           Parameters for bisection.\n*           An interval [LEFT,RIGHT] has converged if\n*           RIGHT-LEFT.LT.MAX( RTOL1*GAP, RTOL2*MAX(|LEFT|,|RIGHT|) )\n*\n*  SPLTOL  (input) DOUBLE PRECISION\n*          The threshold for splitting.\n*\n*  NSPLIT  (output) INTEGER\n*          The number of blocks T splits into. 1 <= NSPLIT <= N.\n*\n*  ISPLIT  (output) INTEGER array, dimension (N)\n*          The splitting points, at which T breaks up into blocks.\n*          The first block consists of rows/columns 1 to ISPLIT(1),\n*          the second of rows/columns ISPLIT(1)+1 through ISPLIT(2),\n*          etc., and the NSPLIT-th consists of rows/columns\n*          ISPLIT(NSPLIT-1)+1 through ISPLIT(NSPLIT)=N.\n*\n*  M       (output) INTEGER\n*          The total number of eigenvalues (of all L_i D_i L_i^T)\n*          found.\n*\n*  W       (output) DOUBLE PRECISION array, dimension (N)\n*          The first M elements contain the eigenvalues. The\n*          eigenvalues of each of the blocks, L_i D_i L_i^T, are\n*          sorted in ascending order ( DLARRE may use the\n*          remaining N-M elements as workspace).\n*\n*  WERR    (output) DOUBLE PRECISION array, dimension (N)\n*          The error bound on the corresponding eigenvalue in W.\n*\n*  WGAP    (output) DOUBLE PRECISION array, dimension (N)\n*          The separation from the right neighbor eigenvalue in W.\n*          The gap is only with respect to the eigenvalues of the same block\n*          as each block has its own representation tree.\n*          Exception: at the right end of a block we store the left gap\n*\n*  IBLOCK  (output) INTEGER array, dimension (N)\n*          The indices of the blocks (submatrices) associated with the\n*          corresponding eigenvalues in W; IBLOCK(i)=1 if eigenvalue\n*          W(i) belongs to the first block from the top, =2 if W(i)\n*          belongs to the second block, etc.\n*\n*  INDEXW  (output) INTEGER array, dimension (N)\n*          The indices of the eigenvalues within each block (submatrix);\n*          for example, INDEXW(i)= 10 and IBLOCK(i)=2 imply that the\n*          i-th eigenvalue W(i) is the 10-th eigenvalue in block 2\n*\n*  GERS    (output) DOUBLE PRECISION array, dimension (2*N)\n*          The N Gerschgorin intervals (the i-th Gerschgorin interval\n*          is (GERS(2*i-1), GERS(2*i)).\n*\n*  PIVMIN  (output) DOUBLE PRECISION\n*          The minimum pivot in the Sturm sequence for T.\n*\n*  WORK    (workspace) DOUBLE PRECISION array, dimension (6*N)\n*          Workspace.\n*\n*  IWORK   (workspace) INTEGER array, dimension (5*N)\n*          Workspace.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          > 0:  A problem occured in DLARRE.\n*          < 0:  One of the called subroutines signaled an internal problem.\n*                Needs inspection of the corresponding parameter IINFO\n*                for further information.\n*\n*          =-1:  Problem in DLARRD.\n*          = 2:  No base representation could be found in MAXTRY iterations.\n*                Increasing MAXTRY and recompilation might be a remedy.\n*          =-3:  Problem in DLARRB when computing the refined root\n*                representation for DLASQ2.\n*          =-4:  Problem in DLARRB when preforming bisection on the\n*                desired part of the spectrum.\n*          =-5:  Problem in DLASQ2.\n*          =-6:  Problem in DLASQ2.\n*\n\n*  Further Details\n*  The base representations are required to suffer very little\n*  element growth and consequently define all their eigenvalues to\n*  high relative accuracy.\n*  ===============\n*\n*  Based on contributions by\n*     Beresford Parlett, University of California, Berkeley, USA\n*     Jim Demmel, University of California, Berkeley, USA\n*     Inderjit Dhillon, University of Texas, Austin, USA\n*     Osni Marques, LBNL/NERSC, USA\n*     Christof Voemel, University of California, Berkeley, USA\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 11)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 11)", argc);
  rb_range = argv[0];
  rb_vl = argv[1];
  rb_vu = argv[2];
  rb_il = argv[3];
  rb_iu = argv[4];
  rb_d = argv[5];
  rb_e = argv[6];
  rb_e2 = argv[7];
  rb_rtol1 = argv[8];
  rb_rtol2 = argv[9];
  rb_spltol = argv[10];

  vl = NUM2DBL(rb_vl);
  iu = NUM2INT(rb_iu);
  if (!NA_IsNArray(rb_e2))
    rb_raise(rb_eArgError, "e2 (8th argument) must be NArray");
  if (NA_RANK(rb_e2) != 1)
    rb_raise(rb_eArgError, "rank of e2 (8th argument) must be %d", 1);
  n = NA_SHAPE0(rb_e2);
  if (NA_TYPE(rb_e2) != NA_DFLOAT)
    rb_e2 = na_change_type(rb_e2, NA_DFLOAT);
  e2 = NA_PTR_TYPE(rb_e2, doublereal*);
  il = NUM2INT(rb_il);
  spltol = NUM2DBL(rb_spltol);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (6th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (6th argument) must be %d", 1);
  if (NA_SHAPE0(rb_d) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of d must be the same as shape 0 of e2");
  if (NA_TYPE(rb_d) != NA_DFLOAT)
    rb_d = na_change_type(rb_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rb_d, doublereal*);
  if (!NA_IsNArray(rb_e))
    rb_raise(rb_eArgError, "e (7th argument) must be NArray");
  if (NA_RANK(rb_e) != 1)
    rb_raise(rb_eArgError, "rank of e (7th argument) must be %d", 1);
  if (NA_SHAPE0(rb_e) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of e must be the same as shape 0 of e2");
  if (NA_TYPE(rb_e) != NA_DFLOAT)
    rb_e = na_change_type(rb_e, NA_DFLOAT);
  e = NA_PTR_TYPE(rb_e, doublereal*);
  range = StringValueCStr(rb_range)[0];
  rtol1 = NUM2DBL(rb_rtol1);
  vu = NUM2DBL(rb_vu);
  rtol2 = NUM2DBL(rb_rtol2);
  {
    int shape[1];
    shape[0] = n;
    rb_isplit = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  isplit = NA_PTR_TYPE(rb_isplit, integer*);
  {
    int shape[1];
    shape[0] = n;
    rb_w = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  w = NA_PTR_TYPE(rb_w, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_werr = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  werr = NA_PTR_TYPE(rb_werr, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_wgap = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  wgap = NA_PTR_TYPE(rb_wgap, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_iblock = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  iblock = NA_PTR_TYPE(rb_iblock, integer*);
  {
    int shape[1];
    shape[0] = n;
    rb_indexw = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  indexw = NA_PTR_TYPE(rb_indexw, integer*);
  {
    int shape[1];
    shape[0] = 2*n;
    rb_gers = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  gers = NA_PTR_TYPE(rb_gers, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_d_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  d_out__ = NA_PTR_TYPE(rb_d_out__, doublereal*);
  MEMCPY(d_out__, d, doublereal, NA_TOTAL(rb_d));
  rb_d = rb_d_out__;
  d = d_out__;
  {
    int shape[1];
    shape[0] = n;
    rb_e_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  e_out__ = NA_PTR_TYPE(rb_e_out__, doublereal*);
  MEMCPY(e_out__, e, doublereal, NA_TOTAL(rb_e));
  rb_e = rb_e_out__;
  e = e_out__;
  {
    int shape[1];
    shape[0] = n;
    rb_e2_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  e2_out__ = NA_PTR_TYPE(rb_e2_out__, doublereal*);
  MEMCPY(e2_out__, e2, doublereal, NA_TOTAL(rb_e2));
  rb_e2 = rb_e2_out__;
  e2 = e2_out__;
  work = ALLOC_N(doublereal, (6*n));
  iwork = ALLOC_N(integer, (5*n));

  dlarre_(&range, &n, &vl, &vu, &il, &iu, d, e, e2, &rtol1, &rtol2, &spltol, &nsplit, isplit, &m, w, werr, wgap, iblock, indexw, gers, &pivmin, work, iwork, &info);

  free(work);
  free(iwork);
  rb_nsplit = INT2NUM(nsplit);
  rb_m = INT2NUM(m);
  rb_pivmin = rb_float_new((double)pivmin);
  rb_info = INT2NUM(info);
  rb_vl = rb_float_new((double)vl);
  rb_vu = rb_float_new((double)vu);
  return rb_ary_new3(16, rb_nsplit, rb_isplit, rb_m, rb_w, rb_werr, rb_wgap, rb_iblock, rb_indexw, rb_gers, rb_pivmin, rb_info, rb_vl, rb_vu, rb_d, rb_e, rb_e2);
}

void
init_lapack_dlarre(VALUE mLapack){
  rb_define_module_function(mLapack, "dlarre", rb_dlarre, -1);
}
