#include "rb_lapack.h"

extern VOID sstebz_(char *range, char *order, integer *n, real *vl, real *vu, integer *il, integer *iu, real *abstol, real *d, real *e, integer *m, integer *nsplit, real *w, integer *iblock, integer *isplit, real *work, integer *iwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_sstebz(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_range;
  char range; 
  VALUE rblapack_order;
  char order; 
  VALUE rblapack_vl;
  real vl; 
  VALUE rblapack_vu;
  real vu; 
  VALUE rblapack_il;
  integer il; 
  VALUE rblapack_iu;
  integer iu; 
  VALUE rblapack_abstol;
  real abstol; 
  VALUE rblapack_d;
  real *d; 
  VALUE rblapack_e;
  real *e; 
  VALUE rblapack_m;
  integer m; 
  VALUE rblapack_nsplit;
  integer nsplit; 
  VALUE rblapack_w;
  real *w; 
  VALUE rblapack_iblock;
  integer *iblock; 
  VALUE rblapack_isplit;
  integer *isplit; 
  VALUE rblapack_info;
  integer info; 
  real *work;
  integer *iwork;

  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  m, nsplit, w, iblock, isplit, info = NumRu::Lapack.sstebz( range, order, vl, vu, il, iu, abstol, d, e, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE SSTEBZ( RANGE, ORDER, N, VL, VU, IL, IU, ABSTOL, D, E, M, NSPLIT, W, IBLOCK, ISPLIT, WORK, IWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  SSTEBZ computes the eigenvalues of a symmetric tridiagonal\n*  matrix T.  The user may ask for all eigenvalues, all eigenvalues\n*  in the half-open interval (VL, VU], or the IL-th through IU-th\n*  eigenvalues.\n*\n*  To avoid overflow, the matrix must be scaled so that its\n*  largest element is no greater than overflow**(1/2) *\n*  underflow**(1/4) in absolute value, and for greatest\n*  accuracy, it should not be much smaller than that.\n*\n*  See W. Kahan \"Accurate Eigenvalues of a Symmetric Tridiagonal\n*  Matrix\", Report CS41, Computer Science Dept., Stanford\n*  University, July 21, 1966.\n*\n\n*  Arguments\n*  =========\n*\n*  RANGE   (input) CHARACTER*1\n*          = 'A': (\"All\")   all eigenvalues will be found.\n*          = 'V': (\"Value\") all eigenvalues in the half-open interval\n*                           (VL, VU] will be found.\n*          = 'I': (\"Index\") the IL-th through IU-th eigenvalues (of the\n*                           entire matrix) will be found.\n*\n*  ORDER   (input) CHARACTER*1\n*          = 'B': (\"By Block\") the eigenvalues will be grouped by\n*                              split-off block (see IBLOCK, ISPLIT) and\n*                              ordered from smallest to largest within\n*                              the block.\n*          = 'E': (\"Entire matrix\")\n*                              the eigenvalues for the entire matrix\n*                              will be ordered from smallest to\n*                              largest.\n*\n*  N       (input) INTEGER\n*          The order of the tridiagonal matrix T.  N >= 0.\n*\n*  VL      (input) REAL\n*  VU      (input) REAL\n*          If RANGE='V', the lower and upper bounds of the interval to\n*          be searched for eigenvalues.  Eigenvalues less than or equal\n*          to VL, or greater than VU, will not be returned.  VL < VU.\n*          Not referenced if RANGE = 'A' or 'I'.\n*\n*  IL      (input) INTEGER\n*  IU      (input) INTEGER\n*          If RANGE='I', the indices (in ascending order) of the\n*          smallest and largest eigenvalues to be returned.\n*          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.\n*          Not referenced if RANGE = 'A' or 'V'.\n*\n*  ABSTOL  (input) REAL\n*          The absolute tolerance for the eigenvalues.  An eigenvalue\n*          (or cluster) is considered to be located if it has been\n*          determined to lie in an interval whose width is ABSTOL or\n*          less.  If ABSTOL is less than or equal to zero, then ULP*|T|\n*          will be used, where |T| means the 1-norm of T.\n*\n*          Eigenvalues will be computed most accurately when ABSTOL is\n*          set to twice the underflow threshold 2*SLAMCH('S'), not zero.\n*\n*  D       (input) REAL array, dimension (N)\n*          The n diagonal elements of the tridiagonal matrix T.\n*\n*  E       (input) REAL array, dimension (N-1)\n*          The (n-1) off-diagonal elements of the tridiagonal matrix T.\n*\n*  M       (output) INTEGER\n*          The actual number of eigenvalues found. 0 <= M <= N.\n*          (See also the description of INFO=2,3.)\n*\n*  NSPLIT  (output) INTEGER\n*          The number of diagonal blocks in the matrix T.\n*          1 <= NSPLIT <= N.\n*\n*  W       (output) REAL array, dimension (N)\n*          On exit, the first M elements of W will contain the\n*          eigenvalues.  (SSTEBZ may use the remaining N-M elements as\n*          workspace.)\n*\n*  IBLOCK  (output) INTEGER array, dimension (N)\n*          At each row/column j where E(j) is zero or small, the\n*          matrix T is considered to split into a block diagonal\n*          matrix.  On exit, if INFO = 0, IBLOCK(i) specifies to which\n*          block (from 1 to the number of blocks) the eigenvalue W(i)\n*          belongs.  (SSTEBZ may use the remaining N-M elements as\n*          workspace.)\n*\n*  ISPLIT  (output) INTEGER array, dimension (N)\n*          The splitting points, at which T breaks up into submatrices.\n*          The first submatrix consists of rows/columns 1 to ISPLIT(1),\n*          the second of rows/columns ISPLIT(1)+1 through ISPLIT(2),\n*          etc., and the NSPLIT-th consists of rows/columns\n*          ISPLIT(NSPLIT-1)+1 through ISPLIT(NSPLIT)=N.\n*          (Only the first NSPLIT elements will actually be used, but\n*          since the user cannot know a priori what value NSPLIT will\n*          have, N words must be reserved for ISPLIT.)\n*\n*  WORK    (workspace) REAL array, dimension (4*N)\n*\n*  IWORK   (workspace) INTEGER array, dimension (3*N)\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*          > 0:  some or all of the eigenvalues failed to converge or\n*                were not computed:\n*                =1 or 3: Bisection failed to converge for some\n*                        eigenvalues; these eigenvalues are flagged by a\n*                        negative block number.  The effect is that the\n*                        eigenvalues may not be as accurate as the\n*                        absolute and relative tolerances.  This is\n*                        generally caused by unexpectedly inaccurate\n*                        arithmetic.\n*                =2 or 3: RANGE='I' only: Not all of the eigenvalues\n*                        IL:IU were found.\n*                        Effect: M < IU+1-IL\n*                        Cause:  non-monotonic arithmetic, causing the\n*                                Sturm sequence to be non-monotonic.\n*                        Cure:   recalculate, using RANGE='A', and pick\n*                                out eigenvalues IL:IU.  In some cases,\n*                                increasing the PARAMETER \"FUDGE\" may\n*                                make things work.\n*                = 4:    RANGE='I', and the Gershgorin interval\n*                        initially used was too small.  No eigenvalues\n*                        were computed.\n*                        Probable cause: your machine has sloppy\n*                                        floating-point arithmetic.\n*                        Cure: Increase the PARAMETER \"FUDGE\",\n*                              recompile, and try again.\n*\n*  Internal Parameters\n*  ===================\n*\n*  RELFAC  REAL, default = 2.0e0\n*          The relative tolerance.  An interval (a,b] lies within\n*          \"relative tolerance\" if  b-a < RELFAC*ulp*max(|a|,|b|),\n*          where \"ulp\" is the machine precision (distance from 1 to\n*          the next larger floating point number.)\n*\n*  FUDGE   REAL, default = 2\n*          A \"fudge factor\" to widen the Gershgorin intervals.  Ideally,\n*          a value of 1 should work, but on machines with sloppy\n*          arithmetic, this needs to be larger.  The default for\n*          publicly released versions should be large enough to handle\n*          the worst machine around.  Note that this has no effect\n*          on accuracy of the solution.\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  m, nsplit, w, iblock, isplit, info = NumRu::Lapack.sstebz( range, order, vl, vu, il, iu, abstol, d, e, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 9)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 9)", argc);
  rblapack_range = argv[0];
  rblapack_order = argv[1];
  rblapack_vl = argv[2];
  rblapack_vu = argv[3];
  rblapack_il = argv[4];
  rblapack_iu = argv[5];
  rblapack_abstol = argv[6];
  rblapack_d = argv[7];
  rblapack_e = argv[8];
  if (rb_options != Qnil) {
  }

  abstol = (real)NUM2DBL(rblapack_abstol);
  vl = (real)NUM2DBL(rblapack_vl);
  iu = NUM2INT(rblapack_iu);
  il = NUM2INT(rblapack_il);
  range = StringValueCStr(rblapack_range)[0];
  vu = (real)NUM2DBL(rblapack_vu);
  order = StringValueCStr(rblapack_order)[0];
  if (!NA_IsNArray(rblapack_d))
    rb_raise(rb_eArgError, "d (8th argument) must be NArray");
  if (NA_RANK(rblapack_d) != 1)
    rb_raise(rb_eArgError, "rank of d (8th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_d);
  if (NA_TYPE(rblapack_d) != NA_SFLOAT)
    rblapack_d = na_change_type(rblapack_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rblapack_d, real*);
  if (!NA_IsNArray(rblapack_e))
    rb_raise(rb_eArgError, "e (9th argument) must be NArray");
  if (NA_RANK(rblapack_e) != 1)
    rb_raise(rb_eArgError, "rank of e (9th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_e) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of e must be %d", n-1);
  if (NA_TYPE(rblapack_e) != NA_SFLOAT)
    rblapack_e = na_change_type(rblapack_e, NA_SFLOAT);
  e = NA_PTR_TYPE(rblapack_e, real*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_w = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  w = NA_PTR_TYPE(rblapack_w, real*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_iblock = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  iblock = NA_PTR_TYPE(rblapack_iblock, integer*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_isplit = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  isplit = NA_PTR_TYPE(rblapack_isplit, integer*);
  work = ALLOC_N(real, (4*n));
  iwork = ALLOC_N(integer, (3*n));

  sstebz_(&range, &order, &n, &vl, &vu, &il, &iu, &abstol, d, e, &m, &nsplit, w, iblock, isplit, work, iwork, &info);

  free(work);
  free(iwork);
  rblapack_m = INT2NUM(m);
  rblapack_nsplit = INT2NUM(nsplit);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(6, rblapack_m, rblapack_nsplit, rblapack_w, rblapack_iblock, rblapack_isplit, rblapack_info);
}

void
init_lapack_sstebz(VALUE mLapack){
  rb_define_module_function(mLapack, "sstebz", rblapack_sstebz, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
