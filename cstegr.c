#include "rb_lapack.h"

extern VOID cstegr_(char *jobz, char *range, integer *n, real *d, real *e, real *vl, real *vu, integer *il, integer *iu, real *abstol, integer *m, real *w, complex *z, integer *ldz, integer *isuppz, real *work, integer *lwork, integer *iwork, integer *liwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_cstegr(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_jobz;
  char jobz; 
  VALUE rblapack_range;
  char range; 
  VALUE rblapack_d;
  real *d; 
  VALUE rblapack_e;
  real *e; 
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
  VALUE rblapack_lwork;
  integer lwork; 
  VALUE rblapack_liwork;
  integer liwork; 
  VALUE rblapack_m;
  integer m; 
  VALUE rblapack_w;
  real *w; 
  VALUE rblapack_z;
  complex *z; 
  VALUE rblapack_isuppz;
  integer *isuppz; 
  VALUE rblapack_work;
  real *work; 
  VALUE rblapack_iwork;
  integer *iwork; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_d_out__;
  real *d_out__;
  VALUE rblapack_e_out__;
  real *e_out__;

  integer n;
  integer ldz;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  m, w, z, isuppz, work, iwork, info, d, e = NumRu::Lapack.cstegr( jobz, range, d, e, vl, vu, il, iu, abstol, lwork, liwork, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE CSTEGR( JOBZ, RANGE, N, D, E, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK, IWORK, LIWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  CSTEGR computes selected eigenvalues and, optionally, eigenvectors\n*  of a real symmetric tridiagonal matrix T. Any such unreduced matrix has\n*  a well defined set of pairwise different real eigenvalues, the corresponding\n*  real eigenvectors are pairwise orthogonal.\n*\n*  The spectrum may be computed either completely or partially by specifying\n*  either an interval (VL,VU] or a range of indices IL:IU for the desired\n*  eigenvalues.\n*\n*  CSTEGR is a compatability wrapper around the improved CSTEMR routine.\n*  See SSTEMR for further details.\n*\n*  One important change is that the ABSTOL parameter no longer provides any\n*  benefit and hence is no longer used.\n*\n*  Note : CSTEGR and CSTEMR work only on machines which follow\n*  IEEE-754 floating-point standard in their handling of infinities and\n*  NaNs.  Normal execution may create these exceptiona values and hence\n*  may abort due to a floating point exception in environments which\n*  do not conform to the IEEE-754 standard.\n*\n\n*  Arguments\n*  =========\n*\n*  JOBZ    (input) CHARACTER*1\n*          = 'N':  Compute eigenvalues only;\n*          = 'V':  Compute eigenvalues and eigenvectors.\n*\n*  RANGE   (input) CHARACTER*1\n*          = 'A': all eigenvalues will be found.\n*          = 'V': all eigenvalues in the half-open interval (VL,VU]\n*                 will be found.\n*          = 'I': the IL-th through IU-th eigenvalues will be found.\n*\n*  N       (input) INTEGER\n*          The order of the matrix.  N >= 0.\n*\n*  D       (input/output) REAL array, dimension (N)\n*          On entry, the N diagonal elements of the tridiagonal matrix\n*          T. On exit, D is overwritten.\n*\n*  E       (input/output) REAL array, dimension (N)\n*          On entry, the (N-1) subdiagonal elements of the tridiagonal\n*          matrix T in elements 1 to N-1 of E. E(N) need not be set on\n*          input, but is used internally as workspace.\n*          On exit, E is overwritten.\n*\n*  VL      (input) REAL\n*  VU      (input) REAL\n*          If RANGE='V', the lower and upper bounds of the interval to\n*          be searched for eigenvalues. VL < VU.\n*          Not referenced if RANGE = 'A' or 'I'.\n*\n*  IL      (input) INTEGER\n*  IU      (input) INTEGER\n*          If RANGE='I', the indices (in ascending order) of the\n*          smallest and largest eigenvalues to be returned.\n*          1 <= IL <= IU <= N, if N > 0.\n*          Not referenced if RANGE = 'A' or 'V'.\n*\n*  ABSTOL  (input) REAL\n*          Unused.  Was the absolute error tolerance for the\n*          eigenvalues/eigenvectors in previous versions.\n*\n*  M       (output) INTEGER\n*          The total number of eigenvalues found.  0 <= M <= N.\n*          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1.\n*\n*  W       (output) REAL array, dimension (N)\n*          The first M elements contain the selected eigenvalues in\n*          ascending order.\n*\n*  Z       (output) COMPLEX array, dimension (LDZ, max(1,M) )\n*          If JOBZ = 'V', and if INFO = 0, then the first M columns of Z\n*          contain the orthonormal eigenvectors of the matrix T\n*          corresponding to the selected eigenvalues, with the i-th\n*          column of Z holding the eigenvector associated with W(i).\n*          If JOBZ = 'N', then Z is not referenced.\n*          Note: the user must ensure that at least max(1,M) columns are\n*          supplied in the array Z; if RANGE = 'V', the exact value of M\n*          is not known in advance and an upper bound must be used.\n*          Supplying N columns is always safe.\n*\n*  LDZ     (input) INTEGER\n*          The leading dimension of the array Z.  LDZ >= 1, and if\n*          JOBZ = 'V', then LDZ >= max(1,N).\n*\n*  ISUPPZ  (output) INTEGER ARRAY, dimension ( 2*max(1,M) )\n*          The support of the eigenvectors in Z, i.e., the indices\n*          indicating the nonzero elements in Z. The i-th computed eigenvector\n*          is nonzero only in elements ISUPPZ( 2*i-1 ) through\n*          ISUPPZ( 2*i ). This is relevant in the case when the matrix\n*          is split. ISUPPZ is only accessed when JOBZ is 'V' and N > 0.\n*\n*  WORK    (workspace/output) REAL array, dimension (LWORK)\n*          On exit, if INFO = 0, WORK(1) returns the optimal\n*          (and minimal) LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK. LWORK >= max(1,18*N)\n*          if JOBZ = 'V', and LWORK >= max(1,12*N) if JOBZ = 'N'.\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  IWORK   (workspace/output) INTEGER array, dimension (LIWORK)\n*          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.\n*\n*  LIWORK  (input) INTEGER\n*          The dimension of the array IWORK.  LIWORK >= max(1,10*N)\n*          if the eigenvectors are desired, and LIWORK >= max(1,8*N)\n*          if only the eigenvalues are to be computed.\n*          If LIWORK = -1, then a workspace query is assumed; the\n*          routine only calculates the optimal size of the IWORK array,\n*          returns this value as the first entry of the IWORK array, and\n*          no error message related to LIWORK is issued by XERBLA.\n*\n*  INFO    (output) INTEGER\n*          On exit, INFO\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*          > 0:  if INFO = 1X, internal error in SLARRE,\n*                if INFO = 2X, internal error in CLARRV.\n*                Here, the digit X = ABS( IINFO ) < 10, where IINFO is\n*                the nonzero error code returned by SLARRE or\n*                CLARRV, respectively.\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Inderjit Dhillon, IBM Almaden, USA\n*     Osni Marques, LBNL/NERSC, USA\n*     Christof Voemel, LBNL/NERSC, USA\n*\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      LOGICAL TRYRAC\n*     ..\n*     .. External Subroutines ..\n      EXTERNAL CSTEMR\n*     ..\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  m, w, z, isuppz, work, iwork, info, d, e = NumRu::Lapack.cstegr( jobz, range, d, e, vl, vu, il, iu, abstol, lwork, liwork, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 11)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 11)", argc);
  rblapack_jobz = argv[0];
  rblapack_range = argv[1];
  rblapack_d = argv[2];
  rblapack_e = argv[3];
  rblapack_vl = argv[4];
  rblapack_vu = argv[5];
  rblapack_il = argv[6];
  rblapack_iu = argv[7];
  rblapack_abstol = argv[8];
  rblapack_lwork = argv[9];
  rblapack_liwork = argv[10];
  if (rb_options != Qnil) {
  }

  abstol = (real)NUM2DBL(rblapack_abstol);
  vl = (real)NUM2DBL(rblapack_vl);
  iu = NUM2INT(rblapack_iu);
  jobz = StringValueCStr(rblapack_jobz)[0];
  vu = (real)NUM2DBL(rblapack_vu);
  liwork = NUM2INT(rblapack_liwork);
  range = StringValueCStr(rblapack_range)[0];
  if (!NA_IsNArray(rblapack_d))
    rb_raise(rb_eArgError, "d (3th argument) must be NArray");
  if (NA_RANK(rblapack_d) != 1)
    rb_raise(rb_eArgError, "rank of d (3th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_d);
  if (NA_TYPE(rblapack_d) != NA_SFLOAT)
    rblapack_d = na_change_type(rblapack_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rblapack_d, real*);
  lwork = NUM2INT(rblapack_lwork);
  if (!NA_IsNArray(rblapack_e))
    rb_raise(rb_eArgError, "e (4th argument) must be NArray");
  if (NA_RANK(rblapack_e) != 1)
    rb_raise(rb_eArgError, "rank of e (4th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_e) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of e must be the same as shape 0 of d");
  if (NA_TYPE(rblapack_e) != NA_SFLOAT)
    rblapack_e = na_change_type(rblapack_e, NA_SFLOAT);
  e = NA_PTR_TYPE(rblapack_e, real*);
  il = NUM2INT(rblapack_il);
  m = lsame_(&range,"A") ? n : lsame_(&range,"I") ? iu-il+1 : 0;
  ldz = lsame_(&jobz,"V") ? MAX(1,n) : 1;
  {
    int shape[1];
    shape[0] = n;
    rblapack_w = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  w = NA_PTR_TYPE(rblapack_w, real*);
  {
    int shape[2];
    shape[0] = ldz;
    shape[1] = MAX(1,m);
    rblapack_z = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  z = NA_PTR_TYPE(rblapack_z, complex*);
  {
    int shape[1];
    shape[0] = 2*MAX(1,m);
    rblapack_isuppz = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  isuppz = NA_PTR_TYPE(rblapack_isuppz, integer*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rblapack_work = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rblapack_work, real*);
  {
    int shape[1];
    shape[0] = MAX(1,liwork);
    rblapack_iwork = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  iwork = NA_PTR_TYPE(rblapack_iwork, integer*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_d_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  d_out__ = NA_PTR_TYPE(rblapack_d_out__, real*);
  MEMCPY(d_out__, d, real, NA_TOTAL(rblapack_d));
  rblapack_d = rblapack_d_out__;
  d = d_out__;
  {
    int shape[1];
    shape[0] = n;
    rblapack_e_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  e_out__ = NA_PTR_TYPE(rblapack_e_out__, real*);
  MEMCPY(e_out__, e, real, NA_TOTAL(rblapack_e));
  rblapack_e = rblapack_e_out__;
  e = e_out__;

  cstegr_(&jobz, &range, &n, d, e, &vl, &vu, &il, &iu, &abstol, &m, w, z, &ldz, isuppz, work, &lwork, iwork, &liwork, &info);

  rblapack_m = INT2NUM(m);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(9, rblapack_m, rblapack_w, rblapack_z, rblapack_isuppz, rblapack_work, rblapack_iwork, rblapack_info, rblapack_d, rblapack_e);
}

void
init_lapack_cstegr(VALUE mLapack){
  rb_define_module_function(mLapack, "cstegr", rblapack_cstegr, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
