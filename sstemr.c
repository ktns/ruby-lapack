#include "rb_lapack.h"

static VALUE
rb_sstemr(int argc, VALUE *argv, VALUE self){
  VALUE rb_jobz;
  char jobz; 
  VALUE rb_range;
  char range; 
  VALUE rb_d;
  real *d; 
  VALUE rb_e;
  real *e; 
  VALUE rb_vl;
  real vl; 
  VALUE rb_vu;
  real vu; 
  VALUE rb_il;
  integer il; 
  VALUE rb_iu;
  integer iu; 
  VALUE rb_nzc;
  integer nzc; 
  VALUE rb_tryrac;
  logical tryrac; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_liwork;
  integer liwork; 
  VALUE rb_m;
  integer m; 
  VALUE rb_w;
  real *w; 
  VALUE rb_z;
  real *z; 
  VALUE rb_isuppz;
  integer *isuppz; 
  VALUE rb_work;
  real *work; 
  VALUE rb_iwork;
  integer *iwork; 
  VALUE rb_info;
  integer info; 
  VALUE rb_d_out__;
  real *d_out__;
  VALUE rb_e_out__;
  real *e_out__;

  integer n;
  integer ldz;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  m, w, z, isuppz, work, iwork, info, d, e, tryrac = NumRu::Lapack.sstemr( jobz, range, d, e, vl, vu, il, iu, nzc, tryrac, lwork, liwork)\n    or\n  NumRu::Lapack.sstemr  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SSTEMR( JOBZ, RANGE, N, D, E, VL, VU, IL, IU, M, W, Z, LDZ, NZC, ISUPPZ, TRYRAC, WORK, LWORK, IWORK, LIWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  SSTEMR computes selected eigenvalues and, optionally, eigenvectors\n*  of a real symmetric tridiagonal matrix T. Any such unreduced matrix has\n*  a well defined set of pairwise different real eigenvalues, the corresponding\n*  real eigenvectors are pairwise orthogonal.\n*\n*  The spectrum may be computed either completely or partially by specifying\n*  either an interval (VL,VU] or a range of indices IL:IU for the desired\n*  eigenvalues.\n*\n*  Depending on the number of desired eigenvalues, these are computed either\n*  by bisection or the dqds algorithm. Numerically orthogonal eigenvectors are\n*  computed by the use of various suitable L D L^T factorizations near clusters\n*  of close eigenvalues (referred to as RRRs, Relatively Robust\n*  Representations). An informal sketch of the algorithm follows.\n*\n*  For each unreduced block (submatrix) of T,\n*     (a) Compute T - sigma I  = L D L^T, so that L and D\n*         define all the wanted eigenvalues to high relative accuracy.\n*         This means that small relative changes in the entries of D and L\n*         cause only small relative changes in the eigenvalues and\n*         eigenvectors. The standard (unfactored) representation of the\n*         tridiagonal matrix T does not have this property in general.\n*     (b) Compute the eigenvalues to suitable accuracy.\n*         If the eigenvectors are desired, the algorithm attains full\n*         accuracy of the computed eigenvalues only right before\n*         the corresponding vectors have to be computed, see steps c) and d).\n*     (c) For each cluster of close eigenvalues, select a new\n*         shift close to the cluster, find a new factorization, and refine\n*         the shifted eigenvalues to suitable accuracy.\n*     (d) For each eigenvalue with a large enough relative separation compute\n*         the corresponding eigenvector by forming a rank revealing twisted\n*         factorization. Go back to (c) for any clusters that remain.\n*\n*  For more details, see:\n*  - Inderjit S. Dhillon and Beresford N. Parlett: \"Multiple representations\n*    to compute orthogonal eigenvectors of symmetric tridiagonal matrices,\"\n*    Linear Algebra and its Applications, 387(1), pp. 1-28, August 2004.\n*  - Inderjit Dhillon and Beresford Parlett: \"Orthogonal Eigenvectors and\n*    Relative Gaps,\" SIAM Journal on Matrix Analysis and Applications, Vol. 25,\n*    2004.  Also LAPACK Working Note 154.\n*  - Inderjit Dhillon: \"A new O(n^2) algorithm for the symmetric\n*    tridiagonal eigenvalue/eigenvector problem\",\n*    Computer Science Division Technical Report No. UCB/CSD-97-971,\n*    UC Berkeley, May 1997.\n*\n*  Further Details\n*  1.SSTEMR works only on machines which follow IEEE-754\n*  floating-point standard in their handling of infinities and NaNs.\n*  This permits the use of efficient inner loops avoiding a check for\n*  zero divisors.\n*\n\n*  Arguments\n*  =========\n*\n*  JOBZ    (input) CHARACTER*1\n*          = 'N':  Compute eigenvalues only;\n*          = 'V':  Compute eigenvalues and eigenvectors.\n*\n*  RANGE   (input) CHARACTER*1\n*          = 'A': all eigenvalues will be found.\n*          = 'V': all eigenvalues in the half-open interval (VL,VU]\n*                 will be found.\n*          = 'I': the IL-th through IU-th eigenvalues will be found.\n*\n*  N       (input) INTEGER\n*          The order of the matrix.  N >= 0.\n*\n*  D       (input/output) REAL array, dimension (N)\n*          On entry, the N diagonal elements of the tridiagonal matrix\n*          T. On exit, D is overwritten.\n*\n*  E       (input/output) REAL array, dimension (N)\n*          On entry, the (N-1) subdiagonal elements of the tridiagonal\n*          matrix T in elements 1 to N-1 of E. E(N) need not be set on\n*          input, but is used internally as workspace.\n*          On exit, E is overwritten.\n*\n*  VL      (input) REAL\n*  VU      (input) REAL\n*          If RANGE='V', the lower and upper bounds of the interval to\n*          be searched for eigenvalues. VL < VU.\n*          Not referenced if RANGE = 'A' or 'I'.\n*\n*  IL      (input) INTEGER\n*  IU      (input) INTEGER\n*          If RANGE='I', the indices (in ascending order) of the\n*          smallest and largest eigenvalues to be returned.\n*          1 <= IL <= IU <= N, if N > 0.\n*          Not referenced if RANGE = 'A' or 'V'.\n*\n*  M       (output) INTEGER\n*          The total number of eigenvalues found.  0 <= M <= N.\n*          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1.\n*\n*  W       (output) REAL array, dimension (N)\n*          The first M elements contain the selected eigenvalues in\n*          ascending order.\n*\n*  Z       (output) REAL array, dimension (LDZ, max(1,M) )\n*          If JOBZ = 'V', and if INFO = 0, then the first M columns of Z\n*          contain the orthonormal eigenvectors of the matrix T\n*          corresponding to the selected eigenvalues, with the i-th\n*          column of Z holding the eigenvector associated with W(i).\n*          If JOBZ = 'N', then Z is not referenced.\n*          Note: the user must ensure that at least max(1,M) columns are\n*          supplied in the array Z; if RANGE = 'V', the exact value of M\n*          is not known in advance and can be computed with a workspace\n*          query by setting NZC = -1, see below.\n*\n*  LDZ     (input) INTEGER\n*          The leading dimension of the array Z.  LDZ >= 1, and if\n*          JOBZ = 'V', then LDZ >= max(1,N).\n*\n*  NZC     (input) INTEGER\n*          The number of eigenvectors to be held in the array Z.\n*          If RANGE = 'A', then NZC >= max(1,N).\n*          If RANGE = 'V', then NZC >= the number of eigenvalues in (VL,VU].\n*          If RANGE = 'I', then NZC >= IU-IL+1.\n*          If NZC = -1, then a workspace query is assumed; the\n*          routine calculates the number of columns of the array Z that\n*          are needed to hold the eigenvectors.\n*          This value is returned as the first entry of the Z array, and\n*          no error message related to NZC is issued by XERBLA.\n*\n*  ISUPPZ  (output) INTEGER ARRAY, dimension ( 2*max(1,M) )\n*          The support of the eigenvectors in Z, i.e., the indices\n*          indicating the nonzero elements in Z. The i-th computed eigenvector\n*          is nonzero only in elements ISUPPZ( 2*i-1 ) through\n*          ISUPPZ( 2*i ). This is relevant in the case when the matrix\n*          is split. ISUPPZ is only accessed when JOBZ is 'V' and N > 0.\n*\n*  TRYRAC  (input/output) LOGICAL\n*          If TRYRAC.EQ..TRUE., indicates that the code should check whether\n*          the tridiagonal matrix defines its eigenvalues to high relative\n*          accuracy.  If so, the code uses relative-accuracy preserving\n*          algorithms that might be (a bit) slower depending on the matrix.\n*          If the matrix does not define its eigenvalues to high relative\n*          accuracy, the code can uses possibly faster algorithms.\n*          If TRYRAC.EQ..FALSE., the code is not required to guarantee\n*          relatively accurate eigenvalues and can use the fastest possible\n*          techniques.\n*          On exit, a .TRUE. TRYRAC will be set to .FALSE. if the matrix\n*          does not define its eigenvalues to high relative accuracy.\n*\n*  WORK    (workspace/output) REAL array, dimension (LWORK)\n*          On exit, if INFO = 0, WORK(1) returns the optimal\n*          (and minimal) LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK. LWORK >= max(1,18*N)\n*          if JOBZ = 'V', and LWORK >= max(1,12*N) if JOBZ = 'N'.\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  IWORK   (workspace/output) INTEGER array, dimension (LIWORK)\n*          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.\n*\n*  LIWORK  (input) INTEGER\n*          The dimension of the array IWORK.  LIWORK >= max(1,10*N)\n*          if the eigenvectors are desired, and LIWORK >= max(1,8*N)\n*          if only the eigenvalues are to be computed.\n*          If LIWORK = -1, then a workspace query is assumed; the\n*          routine only calculates the optimal size of the IWORK array,\n*          returns this value as the first entry of the IWORK array, and\n*          no error message related to LIWORK is issued by XERBLA.\n*\n*  INFO    (output) INTEGER\n*          On exit, INFO\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*          > 0:  if INFO = 1X, internal error in SLARRE,\n*                if INFO = 2X, internal error in SLARRV.\n*                Here, the digit X = ABS( IINFO ) < 10, where IINFO is\n*                the nonzero error code returned by SLARRE or\n*                SLARRV, respectively.\n*\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Beresford Parlett, University of California, Berkeley, USA\n*     Jim Demmel, University of California, Berkeley, USA\n*     Inderjit Dhillon, University of Texas, Austin, USA\n*     Osni Marques, LBNL/NERSC, USA\n*     Christof Voemel, University of California, Berkeley, USA\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 12)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 12)", argc);
  rb_jobz = argv[0];
  rb_range = argv[1];
  rb_d = argv[2];
  rb_e = argv[3];
  rb_vl = argv[4];
  rb_vu = argv[5];
  rb_il = argv[6];
  rb_iu = argv[7];
  rb_nzc = argv[8];
  rb_tryrac = argv[9];
  rb_lwork = argv[10];
  rb_liwork = argv[11];

  jobz = StringValueCStr(rb_jobz)[0];
  range = StringValueCStr(rb_range)[0];
  vl = (real)NUM2DBL(rb_vl);
  vu = (real)NUM2DBL(rb_vu);
  il = NUM2INT(rb_il);
  iu = NUM2INT(rb_iu);
  nzc = NUM2INT(rb_nzc);
  tryrac = (rb_tryrac == Qtrue);
  lwork = NUM2INT(rb_lwork);
  liwork = NUM2INT(rb_liwork);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (3th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (3th argument) must be %d", 1);
  n = NA_SHAPE0(rb_d);
  if (NA_TYPE(rb_d) != NA_SFLOAT)
    rb_d = na_change_type(rb_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rb_d, real*);
  if (!NA_IsNArray(rb_e))
    rb_raise(rb_eArgError, "e (4th argument) must be NArray");
  if (NA_RANK(rb_e) != 1)
    rb_raise(rb_eArgError, "rank of e (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_e) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of e must be the same as shape 0 of d");
  if (NA_TYPE(rb_e) != NA_SFLOAT)
    rb_e = na_change_type(rb_e, NA_SFLOAT);
  e = NA_PTR_TYPE(rb_e, real*);
  {
    int shape[1];
    shape[0] = n;
    rb_w = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  w = NA_PTR_TYPE(rb_w, real*);
  ldz = lsame_(&jobz,"V") ? MAX(1,n) : 1;
  m = lsame_(&range,"A") ? n : lsame_(&range,"I") ? iu-il+1 : 0;
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
    shape[0] = MAX(1,lwork);
    rb_work = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rb_work, real*);
  {
    int shape[1];
    shape[0] = MAX(1,liwork);
    rb_iwork = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  iwork = NA_PTR_TYPE(rb_iwork, integer*);
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
    rb_e_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  e_out__ = NA_PTR_TYPE(rb_e_out__, real*);
  MEMCPY(e_out__, e, real, NA_TOTAL(rb_e));
  rb_e = rb_e_out__;
  e = e_out__;

  sstemr_(&jobz, &range, &n, d, e, &vl, &vu, &il, &iu, &m, w, z, &ldz, &nzc, isuppz, &tryrac, work, &lwork, iwork, &liwork, &info);

  rb_m = INT2NUM(m);
  rb_info = INT2NUM(info);
  rb_tryrac = tryrac ? Qtrue : Qfalse;
  return rb_ary_new3(10, rb_m, rb_w, rb_z, rb_isuppz, rb_work, rb_iwork, rb_info, rb_d, rb_e, rb_tryrac);
}

void
init_lapack_sstemr(VALUE mLapack){
  rb_define_module_function(mLapack, "sstemr", rb_sstemr, -1);
}
