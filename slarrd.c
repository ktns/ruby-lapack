#include "rb_lapack.h"

static VALUE
rb_slarrd(int argc, VALUE *argv, VALUE self){
  VALUE rb_range;
  char range; 
  VALUE rb_order;
  char order; 
  VALUE rb_vl;
  real vl; 
  VALUE rb_vu;
  real vu; 
  VALUE rb_il;
  integer il; 
  VALUE rb_iu;
  integer iu; 
  VALUE rb_gers;
  real *gers; 
  VALUE rb_reltol;
  real reltol; 
  VALUE rb_d;
  real *d; 
  VALUE rb_e;
  real *e; 
  VALUE rb_e2;
  real *e2; 
  VALUE rb_pivmin;
  real pivmin; 
  VALUE rb_nsplit;
  integer nsplit; 
  VALUE rb_isplit;
  integer *isplit; 
  VALUE rb_m;
  integer m; 
  VALUE rb_w;
  real *w; 
  VALUE rb_werr;
  real *werr; 
  VALUE rb_wl;
  real wl; 
  VALUE rb_wu;
  real wu; 
  VALUE rb_iblock;
  integer *iblock; 
  VALUE rb_indexw;
  integer *indexw; 
  VALUE rb_info;
  integer info; 
  real *work;
  integer *iwork;

  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  m, w, werr, wl, wu, iblock, indexw, info = NumRu::Lapack.slarrd( range, order, vl, vu, il, iu, gers, reltol, d, e, e2, pivmin, nsplit, isplit)\n    or\n  NumRu::Lapack.slarrd  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLARRD( RANGE, ORDER, N, VL, VU, IL, IU, GERS, RELTOL, D, E, E2, PIVMIN, NSPLIT, ISPLIT, M, W, WERR, WL, WU, IBLOCK, INDEXW, WORK, IWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  SLARRD computes the eigenvalues of a symmetric tridiagonal\n*  matrix T to suitable accuracy. This is an auxiliary code to be\n*  called from SSTEMR.\n*  The user may ask for all eigenvalues, all eigenvalues\n*  in the half-open interval (VL, VU], or the IL-th through IU-th\n*  eigenvalues.\n*\n*  To avoid overflow, the matrix must be scaled so that its\n*  largest element is no greater than overflow**(1/2) *\n*  underflow**(1/4) in absolute value, and for greatest\n*  accuracy, it should not be much smaller than that.\n*\n*  See W. Kahan \"Accurate Eigenvalues of a Symmetric Tridiagonal\n*  Matrix\", Report CS41, Computer Science Dept., Stanford\n*  University, July 21, 1966.\n*\n\n*  Arguments\n*  =========\n*\n*  RANGE   (input) CHARACTER\n*          = 'A': (\"All\")   all eigenvalues will be found.\n*          = 'V': (\"Value\") all eigenvalues in the half-open interval\n*                           (VL, VU] will be found.\n*          = 'I': (\"Index\") the IL-th through IU-th eigenvalues (of the\n*                           entire matrix) will be found.\n*\n*  ORDER   (input) CHARACTER\n*          = 'B': (\"By Block\") the eigenvalues will be grouped by\n*                              split-off block (see IBLOCK, ISPLIT) and\n*                              ordered from smallest to largest within\n*                              the block.\n*          = 'E': (\"Entire matrix\")\n*                              the eigenvalues for the entire matrix\n*                              will be ordered from smallest to\n*                              largest.\n*\n*  N       (input) INTEGER\n*          The order of the tridiagonal matrix T.  N >= 0.\n*\n*  VL      (input) REAL            \n*  VU      (input) REAL            \n*          If RANGE='V', the lower and upper bounds of the interval to\n*          be searched for eigenvalues.  Eigenvalues less than or equal\n*          to VL, or greater than VU, will not be returned.  VL < VU.\n*          Not referenced if RANGE = 'A' or 'I'.\n*\n*  IL      (input) INTEGER\n*  IU      (input) INTEGER\n*          If RANGE='I', the indices (in ascending order) of the\n*          smallest and largest eigenvalues to be returned.\n*          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.\n*          Not referenced if RANGE = 'A' or 'V'.\n*\n*  GERS    (input) REAL             array, dimension (2*N)\n*          The N Gerschgorin intervals (the i-th Gerschgorin interval\n*          is (GERS(2*i-1), GERS(2*i)).\n*\n*  RELTOL  (input) REAL            \n*          The minimum relative width of an interval.  When an interval\n*          is narrower than RELTOL times the larger (in\n*          magnitude) endpoint, then it is considered to be\n*          sufficiently small, i.e., converged.  Note: this should\n*          always be at least radix*machine epsilon.\n*\n*  D       (input) REAL             array, dimension (N)\n*          The n diagonal elements of the tridiagonal matrix T.\n*\n*  E       (input) REAL             array, dimension (N-1)\n*          The (n-1) off-diagonal elements of the tridiagonal matrix T.\n*\n*  E2      (input) REAL             array, dimension (N-1)\n*          The (n-1) squared off-diagonal elements of the tridiagonal matrix T.\n*\n*  PIVMIN  (input) REAL            \n*          The minimum pivot allowed in the Sturm sequence for T.\n*\n*  NSPLIT  (input) INTEGER\n*          The number of diagonal blocks in the matrix T.\n*          1 <= NSPLIT <= N.\n*\n*  ISPLIT  (input) INTEGER array, dimension (N)\n*          The splitting points, at which T breaks up into submatrices.\n*          The first submatrix consists of rows/columns 1 to ISPLIT(1),\n*          the second of rows/columns ISPLIT(1)+1 through ISPLIT(2),\n*          etc., and the NSPLIT-th consists of rows/columns\n*          ISPLIT(NSPLIT-1)+1 through ISPLIT(NSPLIT)=N.\n*          (Only the first NSPLIT elements will actually be used, but\n*          since the user cannot know a priori what value NSPLIT will\n*          have, N words must be reserved for ISPLIT.)\n*\n*  M       (output) INTEGER\n*          The actual number of eigenvalues found. 0 <= M <= N.\n*          (See also the description of INFO=2,3.)\n*\n*  W       (output) REAL             array, dimension (N)\n*          On exit, the first M elements of W will contain the\n*          eigenvalue approximations. SLARRD computes an interval\n*          I_j = (a_j, b_j] that includes eigenvalue j. The eigenvalue\n*          approximation is given as the interval midpoint\n*          W(j)= ( a_j + b_j)/2. The corresponding error is bounded by\n*          WERR(j) = abs( a_j - b_j)/2\n*\n*  WERR    (output) REAL             array, dimension (N)\n*          The error bound on the corresponding eigenvalue approximation\n*          in W.\n*\n*  WL      (output) REAL            \n*  WU      (output) REAL            \n*          The interval (WL, WU] contains all the wanted eigenvalues.\n*          If RANGE='V', then WL=VL and WU=VU.\n*          If RANGE='A', then WL and WU are the global Gerschgorin bounds\n*                        on the spectrum.\n*          If RANGE='I', then WL and WU are computed by SLAEBZ from the\n*                        index range specified.\n*\n*  IBLOCK  (output) INTEGER array, dimension (N)\n*          At each row/column j where E(j) is zero or small, the\n*          matrix T is considered to split into a block diagonal\n*          matrix.  On exit, if INFO = 0, IBLOCK(i) specifies to which\n*          block (from 1 to the number of blocks) the eigenvalue W(i)\n*          belongs.  (SLARRD may use the remaining N-M elements as\n*          workspace.)\n*\n*  INDEXW  (output) INTEGER array, dimension (N)\n*          The indices of the eigenvalues within each block (submatrix);\n*          for example, INDEXW(i)= j and IBLOCK(i)=k imply that the\n*          i-th eigenvalue W(i) is the j-th eigenvalue in block k.\n*\n*  WORK    (workspace) REAL             array, dimension (4*N)\n*\n*  IWORK   (workspace) INTEGER array, dimension (3*N)\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*          > 0:  some or all of the eigenvalues failed to converge or\n*                were not computed:\n*                =1 or 3: Bisection failed to converge for some\n*                        eigenvalues; these eigenvalues are flagged by a\n*                        negative block number.  The effect is that the\n*                        eigenvalues may not be as accurate as the\n*                        absolute and relative tolerances.  This is\n*                        generally caused by unexpectedly inaccurate\n*                        arithmetic.\n*                =2 or 3: RANGE='I' only: Not all of the eigenvalues\n*                        IL:IU were found.\n*                        Effect: M < IU+1-IL\n*                        Cause:  non-monotonic arithmetic, causing the\n*                                Sturm sequence to be non-monotonic.\n*                        Cure:   recalculate, using RANGE='A', and pick\n*                                out eigenvalues IL:IU.  In some cases,\n*                                increasing the PARAMETER \"FUDGE\" may\n*                                make things work.\n*                = 4:    RANGE='I', and the Gershgorin interval\n*                        initially used was too small.  No eigenvalues\n*                        were computed.\n*                        Probable cause: your machine has sloppy\n*                                        floating-point arithmetic.\n*                        Cure: Increase the PARAMETER \"FUDGE\",\n*                              recompile, and try again.\n*\n*  Internal Parameters\n*  ===================\n*\n*  FUDGE   REAL            , default = 2\n*          A \"fudge factor\" to widen the Gershgorin intervals.  Ideally,\n*          a value of 1 should work, but on machines with sloppy\n*          arithmetic, this needs to be larger.  The default for\n*          publicly released versions should be large enough to handle\n*          the worst machine around.  Note that this has no effect\n*          on accuracy of the solution.\n*\n*  Based on contributions by\n*     W. Kahan, University of California, Berkeley, USA\n*     Beresford Parlett, University of California, Berkeley, USA\n*     Jim Demmel, University of California, Berkeley, USA\n*     Inderjit Dhillon, University of Texas, Austin, USA\n*     Osni Marques, LBNL/NERSC, USA\n*     Christof Voemel, University of California, Berkeley, USA\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 14)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 14)", argc);
  rb_range = argv[0];
  rb_order = argv[1];
  rb_vl = argv[2];
  rb_vu = argv[3];
  rb_il = argv[4];
  rb_iu = argv[5];
  rb_gers = argv[6];
  rb_reltol = argv[7];
  rb_d = argv[8];
  rb_e = argv[9];
  rb_e2 = argv[10];
  rb_pivmin = argv[11];
  rb_nsplit = argv[12];
  rb_isplit = argv[13];

  reltol = (real)NUM2DBL(rb_reltol);
  range = StringValueCStr(rb_range)[0];
  order = StringValueCStr(rb_order)[0];
  vu = (real)NUM2DBL(rb_vu);
  iu = NUM2INT(rb_iu);
  vl = (real)NUM2DBL(rb_vl);
  pivmin = (real)NUM2DBL(rb_pivmin);
  nsplit = NUM2INT(rb_nsplit);
  il = NUM2INT(rb_il);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (1th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (1th argument) must be %d", 1);
  n = NA_SHAPE0(rb_d);
  if (NA_TYPE(rb_d) != NA_SFLOAT)
    rb_d = na_change_type(rb_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rb_d, real*);
  if (!NA_IsNArray(rb_isplit))
    rb_raise(rb_eArgError, "isplit (2th argument) must be NArray");
  if (NA_RANK(rb_isplit) != 1)
    rb_raise(rb_eArgError, "rank of isplit (2th argument) must be %d", 1);
  if (NA_SHAPE0(rb_isplit) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of isplit must be the same as shape 0 of d");
  if (NA_TYPE(rb_isplit) != NA_LINT)
    rb_isplit = na_change_type(rb_isplit, NA_LINT);
  isplit = NA_PTR_TYPE(rb_isplit, integer*);
  if (!NA_IsNArray(rb_e))
    rb_raise(rb_eArgError, "e (7th argument) must be NArray");
  if (NA_RANK(rb_e) != 1)
    rb_raise(rb_eArgError, "rank of e (7th argument) must be %d", 1);
  if (NA_SHAPE0(rb_e) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of e must be %d", n-1);
  if (NA_TYPE(rb_e) != NA_SFLOAT)
    rb_e = na_change_type(rb_e, NA_SFLOAT);
  e = NA_PTR_TYPE(rb_e, real*);
  if (!NA_IsNArray(rb_e2))
    rb_raise(rb_eArgError, "e2 (9th argument) must be NArray");
  if (NA_RANK(rb_e2) != 1)
    rb_raise(rb_eArgError, "rank of e2 (9th argument) must be %d", 1);
  if (NA_SHAPE0(rb_e2) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of e2 must be %d", n-1);
  if (NA_TYPE(rb_e2) != NA_SFLOAT)
    rb_e2 = na_change_type(rb_e2, NA_SFLOAT);
  e2 = NA_PTR_TYPE(rb_e2, real*);
  if (!NA_IsNArray(rb_gers))
    rb_raise(rb_eArgError, "gers (12th argument) must be NArray");
  if (NA_RANK(rb_gers) != 1)
    rb_raise(rb_eArgError, "rank of gers (12th argument) must be %d", 1);
  if (NA_SHAPE0(rb_gers) != (2*n))
    rb_raise(rb_eRuntimeError, "shape 0 of gers must be %d", 2*n);
  if (NA_TYPE(rb_gers) != NA_SFLOAT)
    rb_gers = na_change_type(rb_gers, NA_SFLOAT);
  gers = NA_PTR_TYPE(rb_gers, real*);
  {
    int shape[1];
    shape[0] = n;
    rb_w = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  w = NA_PTR_TYPE(rb_w, real*);
  {
    int shape[1];
    shape[0] = n;
    rb_werr = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  werr = NA_PTR_TYPE(rb_werr, real*);
  {
    int shape[1];
    shape[0] = DIM_LEN(n);
    rb_iblock = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  iblock = NA_PTR_TYPE(rb_iblock, integer*);
  {
    int shape[1];
    shape[0] = DIM_LEN(n);
    rb_indexw = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  indexw = NA_PTR_TYPE(rb_indexw, integer*);
  work = ALLOC_N(real, (4*n));
  iwork = ALLOC_N(integer, (3*n));

  slarrd_(&range, &order, &n, &vl, &vu, &il, &iu, gers, &reltol, d, e, e2, &pivmin, &nsplit, isplit, &m, w, werr, &wl, &wu, iblock, indexw, work, iwork, &info);

  free(work);
  free(iwork);
  rb_m = INT2NUM(m);
  rb_wl = rb_float_new((double)wl);
  rb_wu = rb_float_new((double)wu);
  rb_info = INT2NUM(info);
  return rb_ary_new3(8, rb_m, rb_w, rb_werr, rb_wl, rb_wu, rb_iblock, rb_indexw, rb_info);
}

void
init_lapack_slarrd(VALUE mLapack){
  rb_define_module_function(mLapack, "slarrd", rb_slarrd, -1);
}
