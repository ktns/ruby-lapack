#include "rb_lapack.h"

extern VOID dstegr_(char *jobz, char *range, integer *n, doublereal *d, doublereal *e, doublereal *vl, doublereal *vu, integer *il, integer *iu, doublereal *abstol, integer *m, doublereal *w, doublereal *z, integer *ldz, integer *isuppz, doublereal *work, integer *lwork, integer *iwork, integer *liwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dstegr(int argc, VALUE *argv, VALUE self){
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
  VALUE rblapack_abstol;
  doublereal abstol; 
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
      printf("%s\n", "USAGE:\n  m, w, z, isuppz, work, iwork, info, d, e = NumRu::Lapack.dstegr( jobz, range, d, e, vl, vu, il, iu, abstol, lwork, liwork, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DSTEGR( JOBZ, RANGE, N, D, E, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK, IWORK, LIWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  DSTEGR computes selected eigenvalues and, optionally, eigenvectors\n*  of a real symmetric tridiagonal matrix T. Any such unreduced matrix has\n*  a well defined set of pairwise different real eigenvalues, the corresponding\n*  real eigenvectors are pairwise orthogonal.\n*\n*  The spectrum may be computed either completely or partially by specifying\n*  either an interval (VL,VU] or a range of indices IL:IU for the desired\n*  eigenvalues.\n*\n*  DSTEGR is a compatability wrapper around the improved DSTEMR routine.\n*  See DSTEMR for further details.\n*\n*  One important change is that the ABSTOL parameter no longer provides any\n*  benefit and hence is no longer used.\n*\n*  Note : DSTEGR and DSTEMR work only on machines which follow\n*  IEEE-754 floating-point standard in their handling of infinities and\n*  NaNs.  Normal execution may create these exceptiona values and hence\n*  may abort due to a floating point exception in environments which\n*  do not conform to the IEEE-754 standard.\n*\n\n*  Arguments\n*  =========\n*\n*  JOBZ    (input) CHARACTER*1\n*          = 'N':  Compute eigenvalues only;\n*          = 'V':  Compute eigenvalues and eigenvectors.\n*\n*  RANGE   (input) CHARACTER*1\n*          = 'A': all eigenvalues will be found.\n*          = 'V': all eigenvalues in the half-open interval (VL,VU]\n*                 will be found.\n*          = 'I': the IL-th through IU-th eigenvalues will be found.\n*\n*  N       (input) INTEGER\n*          The order of the matrix.  N >= 0.\n*\n*  D       (input/output) DOUBLE PRECISION array, dimension (N)\n*          On entry, the N diagonal elements of the tridiagonal matrix\n*          T. On exit, D is overwritten.\n*\n*  E       (input/output) DOUBLE PRECISION array, dimension (N)\n*          On entry, the (N-1) subdiagonal elements of the tridiagonal\n*          matrix T in elements 1 to N-1 of E. E(N) need not be set on\n*          input, but is used internally as workspace.\n*          On exit, E is overwritten.\n*\n*  VL      (input) DOUBLE PRECISION\n*  VU      (input) DOUBLE PRECISION\n*          If RANGE='V', the lower and upper bounds of the interval to\n*          be searched for eigenvalues. VL < VU.\n*          Not referenced if RANGE = 'A' or 'I'.\n*\n*  IL      (input) INTEGER\n*  IU      (input) INTEGER\n*          If RANGE='I', the indices (in ascending order) of the\n*          smallest and largest eigenvalues to be returned.\n*          1 <= IL <= IU <= N, if N > 0.\n*          Not referenced if RANGE = 'A' or 'V'.\n*\n*  ABSTOL  (input) DOUBLE PRECISION\n*          Unused.  Was the absolute error tolerance for the\n*          eigenvalues/eigenvectors in previous versions.\n*\n*  M       (output) INTEGER\n*          The total number of eigenvalues found.  0 <= M <= N.\n*          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1.\n*\n*  W       (output) DOUBLE PRECISION array, dimension (N)\n*          The first M elements contain the selected eigenvalues in\n*          ascending order.\n*\n*  Z       (output) DOUBLE PRECISION array, dimension (LDZ, max(1,M) )\n*          If JOBZ = 'V', and if INFO = 0, then the first M columns of Z\n*          contain the orthonormal eigenvectors of the matrix T\n*          corresponding to the selected eigenvalues, with the i-th\n*          column of Z holding the eigenvector associated with W(i).\n*          If JOBZ = 'N', then Z is not referenced.\n*          Note: the user must ensure that at least max(1,M) columns are\n*          supplied in the array Z; if RANGE = 'V', the exact value of M\n*          is not known in advance and an upper bound must be used.\n*          Supplying N columns is always safe.\n*\n*  LDZ     (input) INTEGER\n*          The leading dimension of the array Z.  LDZ >= 1, and if\n*          JOBZ = 'V', then LDZ >= max(1,N).\n*\n*  ISUPPZ  (output) INTEGER ARRAY, dimension ( 2*max(1,M) )\n*          The support of the eigenvectors in Z, i.e., the indices\n*          indicating the nonzero elements in Z. The i-th computed eigenvector\n*          is nonzero only in elements ISUPPZ( 2*i-1 ) through\n*          ISUPPZ( 2*i ). This is relevant in the case when the matrix\n*          is split. ISUPPZ is only accessed when JOBZ is 'V' and N > 0.\n*\n*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)\n*          On exit, if INFO = 0, WORK(1) returns the optimal\n*          (and minimal) LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK. LWORK >= max(1,18*N)\n*          if JOBZ = 'V', and LWORK >= max(1,12*N) if JOBZ = 'N'.\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  IWORK   (workspace/output) INTEGER array, dimension (LIWORK)\n*          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.\n*\n*  LIWORK  (input) INTEGER\n*          The dimension of the array IWORK.  LIWORK >= max(1,10*N)\n*          if the eigenvectors are desired, and LIWORK >= max(1,8*N)\n*          if only the eigenvalues are to be computed.\n*          If LIWORK = -1, then a workspace query is assumed; the\n*          routine only calculates the optimal size of the IWORK array,\n*          returns this value as the first entry of the IWORK array, and\n*          no error message related to LIWORK is issued by XERBLA.\n*\n*  INFO    (output) INTEGER\n*          On exit, INFO\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*          > 0:  if INFO = 1X, internal error in DLARRE,\n*                if INFO = 2X, internal error in DLARRV.\n*                Here, the digit X = ABS( IINFO ) < 10, where IINFO is\n*                the nonzero error code returned by DLARRE or\n*                DLARRV, respectively.\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Inderjit Dhillon, IBM Almaden, USA\n*     Osni Marques, LBNL/NERSC, USA\n*     Christof Voemel, LBNL/NERSC, USA\n*\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      LOGICAL TRYRAC\n*     ..\n*     .. External Subroutines ..\n      EXTERNAL DSTEMR\n*     ..\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  m, w, z, isuppz, work, iwork, info, d, e = NumRu::Lapack.dstegr( jobz, range, d, e, vl, vu, il, iu, abstol, lwork, liwork, [:usage => usage, :help => help])\n");
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

  abstol = NUM2DBL(rblapack_abstol);
  vl = NUM2DBL(rblapack_vl);
  iu = NUM2INT(rblapack_iu);
  jobz = StringValueCStr(rblapack_jobz)[0];
  vu = NUM2DBL(rblapack_vu);
  liwork = NUM2INT(rblapack_liwork);
  range = StringValueCStr(rblapack_range)[0];
  if (!NA_IsNArray(rblapack_d))
    rb_raise(rb_eArgError, "d (3th argument) must be NArray");
  if (NA_RANK(rblapack_d) != 1)
    rb_raise(rb_eArgError, "rank of d (3th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_d);
  if (NA_TYPE(rblapack_d) != NA_DFLOAT)
    rblapack_d = na_change_type(rblapack_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rblapack_d, doublereal*);
  lwork = NUM2INT(rblapack_lwork);
  if (!NA_IsNArray(rblapack_e))
    rb_raise(rb_eArgError, "e (4th argument) must be NArray");
  if (NA_RANK(rblapack_e) != 1)
    rb_raise(rb_eArgError, "rank of e (4th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_e) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of e must be the same as shape 0 of d");
  if (NA_TYPE(rblapack_e) != NA_DFLOAT)
    rblapack_e = na_change_type(rblapack_e, NA_DFLOAT);
  e = NA_PTR_TYPE(rblapack_e, doublereal*);
  il = NUM2INT(rblapack_il);
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

  dstegr_(&jobz, &range, &n, d, e, &vl, &vu, &il, &iu, &abstol, &m, w, z, &ldz, isuppz, work, &lwork, iwork, &liwork, &info);

  rblapack_m = INT2NUM(m);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(9, rblapack_m, rblapack_w, rblapack_z, rblapack_isuppz, rblapack_work, rblapack_iwork, rblapack_info, rblapack_d, rblapack_e);
}

void
init_lapack_dstegr(VALUE mLapack){
  rb_define_module_function(mLapack, "dstegr", rblapack_dstegr, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
