#include "rb_lapack.h"

extern VOID dstemr_(char *jobz, char *range, integer *n, doublereal *d, doublereal *e, doublereal *vl, doublereal *vu, integer *il, integer *iu, integer *m, doublereal *w, doublereal *z, integer *ldz, integer *nzc, integer *isuppz, logical *tryrac, doublereal *work, integer *lwork, integer *iwork, integer *liwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dstemr(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_jobz;
  char jobz; 
  VALUE rblapack_range;
  char range; 
  VALUE rblapack_d;
  doublereal *d; 
  VALUE rblapack_e;
  doublereal *e; 
  VALUE rblapack_vl;
  doublereal vl; 
  VALUE rblapack_vu;
  doublereal vu; 
  VALUE rblapack_il;
  integer il; 
  VALUE rblapack_iu;
  integer iu; 
  VALUE rblapack_nzc;
  integer nzc; 
  VALUE rblapack_tryrac;
  logical tryrac; 
  VALUE rblapack_lwork;
  integer lwork; 
  VALUE rblapack_liwork;
  integer liwork; 
  VALUE rblapack_m;
  integer m; 
  VALUE rblapack_w;
  doublereal *w; 
  VALUE rblapack_z;
  doublereal *z; 
  VALUE rblapack_isuppz;
  integer *isuppz; 
  VALUE rblapack_work;
  doublereal *work; 
  VALUE rblapack_iwork;
  integer *iwork; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_d_out__;
  doublereal *d_out__;
  VALUE rblapack_e_out__;
  doublereal *e_out__;

  integer n;
  integer ldz;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  m, w, z, isuppz, work, iwork, info, d, e, tryrac = NumRu::Lapack.dstemr( jobz, range, d, e, vl, vu, il, iu, nzc, tryrac, lwork, liwork, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DSTEMR( JOBZ, RANGE, N, D, E, VL, VU, IL, IU, M, W, Z, LDZ, NZC, ISUPPZ, TRYRAC, WORK, LWORK, IWORK, LIWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  DSTEMR computes selected eigenvalues and, optionally, eigenvectors\n*  of a real symmetric tridiagonal matrix T. Any such unreduced matrix has\n*  a well defined set of pairwise different real eigenvalues, the corresponding\n*  real eigenvectors are pairwise orthogonal.\n*\n*  The spectrum may be computed either completely or partially by specifying\n*  either an interval (VL,VU] or a range of indices IL:IU for the desired\n*  eigenvalues.\n*\n*  Depending on the number of desired eigenvalues, these are computed either\n*  by bisection or the dqds algorithm. Numerically orthogonal eigenvectors are\n*  computed by the use of various suitable L D L^T factorizations near clusters\n*  of close eigenvalues (referred to as RRRs, Relatively Robust\n*  Representations). An informal sketch of the algorithm follows.\n*\n*  For each unreduced block (submatrix) of T,\n*     (a) Compute T - sigma I  = L D L^T, so that L and D\n*         define all the wanted eigenvalues to high relative accuracy.\n*         This means that small relative changes in the entries of D and L\n*         cause only small relative changes in the eigenvalues and\n*         eigenvectors. The standard (unfactored) representation of the\n*         tridiagonal matrix T does not have this property in general.\n*     (b) Compute the eigenvalues to suitable accuracy.\n*         If the eigenvectors are desired, the algorithm attains full\n*         accuracy of the computed eigenvalues only right before\n*         the corresponding vectors have to be computed, see steps c) and d).\n*     (c) For each cluster of close eigenvalues, select a new\n*         shift close to the cluster, find a new factorization, and refine\n*         the shifted eigenvalues to suitable accuracy.\n*     (d) For each eigenvalue with a large enough relative separation compute\n*         the corresponding eigenvector by forming a rank revealing twisted\n*         factorization. Go back to (c) for any clusters that remain.\n*\n*  For more details, see:\n*  - Inderjit S. Dhillon and Beresford N. Parlett: \"Multiple representations\n*    to compute orthogonal eigenvectors of symmetric tridiagonal matrices,\"\n*    Linear Algebra and its Applications, 387(1), pp. 1-28, August 2004.\n*  - Inderjit Dhillon and Beresford Parlett: \"Orthogonal Eigenvectors and\n*    Relative Gaps,\" SIAM Journal on Matrix Analysis and Applications, Vol. 25,\n*    2004.  Also LAPACK Working Note 154.\n*  - Inderjit Dhillon: \"A new O(n^2) algorithm for the symmetric\n*    tridiagonal eigenvalue/eigenvector problem\",\n*    Computer Science Division Technical Report No. UCB/CSD-97-971,\n*    UC Berkeley, May 1997.\n*\n*  Further Details\n*  1.DSTEMR works only on machines which follow IEEE-754\n*  floating-point standard in their handling of infinities and NaNs.\n*  This permits the use of efficient inner loops avoiding a check for\n*  zero divisors.\n*\n\n*  Arguments\n*  =========\n*\n*  JOBZ    (input) CHARACTER*1\n*          = 'N':  Compute eigenvalues only;\n*          = 'V':  Compute eigenvalues and eigenvectors.\n*\n*  RANGE   (input) CHARACTER*1\n*          = 'A': all eigenvalues will be found.\n*          = 'V': all eigenvalues in the half-open interval (VL,VU]\n*                 will be found.\n*          = 'I': the IL-th through IU-th eigenvalues will be found.\n*\n*  N       (input) INTEGER\n*          The order of the matrix.  N >= 0.\n*\n*  D       (input/output) DOUBLE PRECISION array, dimension (N)\n*          On entry, the N diagonal elements of the tridiagonal matrix\n*          T. On exit, D is overwritten.\n*\n*  E       (input/output) DOUBLE PRECISION array, dimension (N)\n*          On entry, the (N-1) subdiagonal elements of the tridiagonal\n*          matrix T in elements 1 to N-1 of E. E(N) need not be set on\n*          input, but is used internally as workspace.\n*          On exit, E is overwritten.\n*\n*  VL      (input) DOUBLE PRECISION\n*  VU      (input) DOUBLE PRECISION\n*          If RANGE='V', the lower and upper bounds of the interval to\n*          be searched for eigenvalues. VL < VU.\n*          Not referenced if RANGE = 'A' or 'I'.\n*\n*  IL      (input) INTEGER\n*  IU      (input) INTEGER\n*          If RANGE='I', the indices (in ascending order) of the\n*          smallest and largest eigenvalues to be returned.\n*          1 <= IL <= IU <= N, if N > 0.\n*          Not referenced if RANGE = 'A' or 'V'.\n*\n*  M       (output) INTEGER\n*          The total number of eigenvalues found.  0 <= M <= N.\n*          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1.\n*\n*  W       (output) DOUBLE PRECISION array, dimension (N)\n*          The first M elements contain the selected eigenvalues in\n*          ascending order.\n*\n*  Z       (output) DOUBLE PRECISION array, dimension (LDZ, max(1,M) )\n*          If JOBZ = 'V', and if INFO = 0, then the first M columns of Z\n*          contain the orthonormal eigenvectors of the matrix T\n*          corresponding to the selected eigenvalues, with the i-th\n*          column of Z holding the eigenvector associated with W(i).\n*          If JOBZ = 'N', then Z is not referenced.\n*          Note: the user must ensure that at least max(1,M) columns are\n*          supplied in the array Z; if RANGE = 'V', the exact value of M\n*          is not known in advance and can be computed with a workspace\n*          query by setting NZC = -1, see below.\n*\n*  LDZ     (input) INTEGER\n*          The leading dimension of the array Z.  LDZ >= 1, and if\n*          JOBZ = 'V', then LDZ >= max(1,N).\n*\n*  NZC     (input) INTEGER\n*          The number of eigenvectors to be held in the array Z.\n*          If RANGE = 'A', then NZC >= max(1,N).\n*          If RANGE = 'V', then NZC >= the number of eigenvalues in (VL,VU].\n*          If RANGE = 'I', then NZC >= IU-IL+1.\n*          If NZC = -1, then a workspace query is assumed; the\n*          routine calculates the number of columns of the array Z that\n*          are needed to hold the eigenvectors.\n*          This value is returned as the first entry of the Z array, and\n*          no error message related to NZC is issued by XERBLA.\n*\n*  ISUPPZ  (output) INTEGER ARRAY, dimension ( 2*max(1,M) )\n*          The support of the eigenvectors in Z, i.e., the indices\n*          indicating the nonzero elements in Z. The i-th computed eigenvector\n*          is nonzero only in elements ISUPPZ( 2*i-1 ) through\n*          ISUPPZ( 2*i ). This is relevant in the case when the matrix\n*          is split. ISUPPZ is only accessed when JOBZ is 'V' and N > 0.\n*\n*  TRYRAC  (input/output) LOGICAL\n*          If TRYRAC.EQ..TRUE., indicates that the code should check whether\n*          the tridiagonal matrix defines its eigenvalues to high relative\n*          accuracy.  If so, the code uses relative-accuracy preserving\n*          algorithms that might be (a bit) slower depending on the matrix.\n*          If the matrix does not define its eigenvalues to high relative\n*          accuracy, the code can uses possibly faster algorithms.\n*          If TRYRAC.EQ..FALSE., the code is not required to guarantee\n*          relatively accurate eigenvalues and can use the fastest possible\n*          techniques.\n*          On exit, a .TRUE. TRYRAC will be set to .FALSE. if the matrix\n*          does not define its eigenvalues to high relative accuracy.\n*\n*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)\n*          On exit, if INFO = 0, WORK(1) returns the optimal\n*          (and minimal) LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK. LWORK >= max(1,18*N)\n*          if JOBZ = 'V', and LWORK >= max(1,12*N) if JOBZ = 'N'.\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  IWORK   (workspace/output) INTEGER array, dimension (LIWORK)\n*          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.\n*\n*  LIWORK  (input) INTEGER\n*          The dimension of the array IWORK.  LIWORK >= max(1,10*N)\n*          if the eigenvectors are desired, and LIWORK >= max(1,8*N)\n*          if only the eigenvalues are to be computed.\n*          If LIWORK = -1, then a workspace query is assumed; the\n*          routine only calculates the optimal size of the IWORK array,\n*          returns this value as the first entry of the IWORK array, and\n*          no error message related to LIWORK is issued by XERBLA.\n*\n*  INFO    (output) INTEGER\n*          On exit, INFO\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*          > 0:  if INFO = 1X, internal error in DLARRE,\n*                if INFO = 2X, internal error in DLARRV.\n*                Here, the digit X = ABS( IINFO ) < 10, where IINFO is\n*                the nonzero error code returned by DLARRE or\n*                DLARRV, respectively.\n*\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Beresford Parlett, University of California, Berkeley, USA\n*     Jim Demmel, University of California, Berkeley, USA\n*     Inderjit Dhillon, University of Texas, Austin, USA\n*     Osni Marques, LBNL/NERSC, USA\n*     Christof Voemel, University of California, Berkeley, USA\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  m, w, z, isuppz, work, iwork, info, d, e, tryrac = NumRu::Lapack.dstemr( jobz, range, d, e, vl, vu, il, iu, nzc, tryrac, lwork, liwork, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 12)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 12)", argc);
  rblapack_jobz = argv[0];
  rblapack_range = argv[1];
  rblapack_d = argv[2];
  rblapack_e = argv[3];
  rblapack_vl = argv[4];
  rblapack_vu = argv[5];
  rblapack_il = argv[6];
  rblapack_iu = argv[7];
  rblapack_nzc = argv[8];
  rblapack_tryrac = argv[9];
  rblapack_lwork = argv[10];
  rblapack_liwork = argv[11];
  if (rb_options != Qnil) {
  }

  vl = NUM2DBL(rblapack_vl);
  nzc = NUM2INT(rblapack_nzc);
  iu = NUM2INT(rblapack_iu);
  jobz = StringValueCStr(rblapack_jobz)[0];
  vu = NUM2DBL(rblapack_vu);
  liwork = NUM2INT(rblapack_liwork);
  il = NUM2INT(rblapack_il);
  tryrac = (rblapack_tryrac == Qtrue);
  lwork = NUM2INT(rblapack_lwork);
  if (!NA_IsNArray(rblapack_e))
    rb_raise(rb_eArgError, "e (4th argument) must be NArray");
  if (NA_RANK(rblapack_e) != 1)
    rb_raise(rb_eArgError, "rank of e (4th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_e);
  if (NA_TYPE(rblapack_e) != NA_DFLOAT)
    rblapack_e = na_change_type(rblapack_e, NA_DFLOAT);
  e = NA_PTR_TYPE(rblapack_e, doublereal*);
  range = StringValueCStr(rblapack_range)[0];
  if (!NA_IsNArray(rblapack_d))
    rb_raise(rb_eArgError, "d (3th argument) must be NArray");
  if (NA_RANK(rblapack_d) != 1)
    rb_raise(rb_eArgError, "rank of d (3th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_d) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of d must be the same as shape 0 of e");
  if (NA_TYPE(rblapack_d) != NA_DFLOAT)
    rblapack_d = na_change_type(rblapack_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rblapack_d, doublereal*);
  m = lsame_(&range,"A") ? n : lsame_(&range,"I") ? iu-il+1 : 0;
  ldz = lsame_(&jobz,"V") ? MAX(1,n) : 1;
  {
    int shape[1];
    shape[0] = n;
    rblapack_w = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  w = NA_PTR_TYPE(rblapack_w, doublereal*);
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
    shape[0] = MAX(1,lwork);
    rblapack_work = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rblapack_work, doublereal*);
  {
    int shape[1];
    shape[0] = MAX(1,liwork);
    rblapack_iwork = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  iwork = NA_PTR_TYPE(rblapack_iwork, integer*);
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
    rblapack_e_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  e_out__ = NA_PTR_TYPE(rblapack_e_out__, doublereal*);
  MEMCPY(e_out__, e, doublereal, NA_TOTAL(rblapack_e));
  rblapack_e = rblapack_e_out__;
  e = e_out__;

  dstemr_(&jobz, &range, &n, d, e, &vl, &vu, &il, &iu, &m, w, z, &ldz, &nzc, isuppz, &tryrac, work, &lwork, iwork, &liwork, &info);

  rblapack_m = INT2NUM(m);
  rblapack_info = INT2NUM(info);
  rblapack_tryrac = tryrac ? Qtrue : Qfalse;
  return rb_ary_new3(10, rblapack_m, rblapack_w, rblapack_z, rblapack_isuppz, rblapack_work, rblapack_iwork, rblapack_info, rblapack_d, rblapack_e, rblapack_tryrac);
}

void
init_lapack_dstemr(VALUE mLapack){
  rb_define_module_function(mLapack, "dstemr", rblapack_dstemr, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
