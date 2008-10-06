#include "rb_lapack.h"

static VALUE
rb_dstegr(int argc, VALUE *argv, VALUE self){
  VALUE rb_jobz;
  char jobz; 
  VALUE rb_range;
  char range; 
  VALUE rb_d;
  doublereal *d; 
  VALUE rb_e;
  doublereal *e; 
  VALUE rb_vl;
  doublereal vl; 
  VALUE rb_vu;
  doublereal vu; 
  VALUE rb_il;
  integer il; 
  VALUE rb_iu;
  integer iu; 
  VALUE rb_abstol;
  doublereal abstol; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_liwork;
  integer liwork; 
  VALUE rb_m;
  integer m; 
  VALUE rb_w;
  doublereal *w; 
  VALUE rb_z;
  doublereal *z; 
  VALUE rb_isuppz;
  integer *isuppz; 
  VALUE rb_work;
  doublereal *work; 
  VALUE rb_iwork;
  integer *iwork; 
  VALUE rb_info;
  integer info; 
  VALUE rb_d_out__;
  doublereal *d_out__;
  VALUE rb_e_out__;
  doublereal *e_out__;

  integer n;
  integer ldz;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  m, w, z, isuppz, work, iwork, info, d, e = NumRu::Lapack.dstegr( jobz, range, d, e, vl, vu, il, iu, abstol, lwork, liwork)\n    or\n  NumRu::Lapack.dstegr  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DSTEGR( JOBZ, RANGE, N, D, E, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK, IWORK, LIWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  DSTEGR computes selected eigenvalues and, optionally, eigenvectors\n*  of a real symmetric tridiagonal matrix T. Any such unreduced matrix has\n*  a well defined set of pairwise different real eigenvalues, the corresponding\n*  real eigenvectors are pairwise orthogonal.\n*\n*  The spectrum may be computed either completely or partially by specifying\n*  either an interval (VL,VU] or a range of indices IL:IU for the desired\n*  eigenvalues.\n*\n*  DSTEGR is a compatability wrapper around the improved DSTEMR routine.\n*  See DSTEMR for further details.\n*\n*  One important change is that the ABSTOL parameter no longer provides any\n*  benefit and hence is no longer used.\n*\n*  Note : DSTEGR and DSTEMR work only on machines which follow\n*  IEEE-754 floating-point standard in their handling of infinities and\n*  NaNs.  Normal execution may create these exceptiona values and hence\n*  may abort due to a floating point exception in environments which\n*  do not conform to the IEEE-754 standard.\n*\n\n*  Arguments\n*  =========\n*\n*  JOBZ    (input) CHARACTER*1\n*          = 'N':  Compute eigenvalues only;\n*          = 'V':  Compute eigenvalues and eigenvectors.\n*\n*  RANGE   (input) CHARACTER*1\n*          = 'A': all eigenvalues will be found.\n*          = 'V': all eigenvalues in the half-open interval (VL,VU]\n*                 will be found.\n*          = 'I': the IL-th through IU-th eigenvalues will be found.\n*\n*  N       (input) INTEGER\n*          The order of the matrix.  N >= 0.\n*\n*  D       (input/output) DOUBLE PRECISION array, dimension (N)\n*          On entry, the N diagonal elements of the tridiagonal matrix\n*          T. On exit, D is overwritten.\n*\n*  E       (input/output) DOUBLE PRECISION array, dimension (N)\n*          On entry, the (N-1) subdiagonal elements of the tridiagonal\n*          matrix T in elements 1 to N-1 of E. E(N) need not be set on\n*          input, but is used internally as workspace.\n*          On exit, E is overwritten.\n*\n*  VL      (input) DOUBLE PRECISION\n*  VU      (input) DOUBLE PRECISION\n*          If RANGE='V', the lower and upper bounds of the interval to\n*          be searched for eigenvalues. VL < VU.\n*          Not referenced if RANGE = 'A' or 'I'.\n*\n*  IL      (input) INTEGER\n*  IU      (input) INTEGER\n*          If RANGE='I', the indices (in ascending order) of the\n*          smallest and largest eigenvalues to be returned.\n*          1 <= IL <= IU <= N, if N > 0.\n*          Not referenced if RANGE = 'A' or 'V'.\n*\n*  ABSTOL  (input) DOUBLE PRECISION\n*          Unused.  Was the absolute error tolerance for the\n*          eigenvalues/eigenvectors in previous versions.\n*\n*  M       (output) INTEGER\n*          The total number of eigenvalues found.  0 <= M <= N.\n*          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1.\n*\n*  W       (output) DOUBLE PRECISION array, dimension (N)\n*          The first M elements contain the selected eigenvalues in\n*          ascending order.\n*\n*  Z       (output) DOUBLE PRECISION array, dimension (LDZ, max(1,M) )\n*          If JOBZ = 'V', and if INFO = 0, then the first M columns of Z\n*          contain the orthonormal eigenvectors of the matrix T\n*          corresponding to the selected eigenvalues, with the i-th\n*          column of Z holding the eigenvector associated with W(i).\n*          If JOBZ = 'N', then Z is not referenced.\n*          Note: the user must ensure that at least max(1,M) columns are\n*          supplied in the array Z; if RANGE = 'V', the exact value of M\n*          is not known in advance and an upper bound must be used.\n*          Supplying N columns is always safe.\n*\n*  LDZ     (input) INTEGER\n*          The leading dimension of the array Z.  LDZ >= 1, and if\n*          JOBZ = 'V', then LDZ >= max(1,N).\n*\n*  ISUPPZ  (output) INTEGER ARRAY, dimension ( 2*max(1,M) )\n*          The support of the eigenvectors in Z, i.e., the indices\n*          indicating the nonzero elements in Z. The i-th computed eigenvector\n*          is nonzero only in elements ISUPPZ( 2*i-1 ) through\n*          ISUPPZ( 2*i ). This is relevant in the case when the matrix\n*          is split. ISUPPZ is only accessed when JOBZ is 'V' and N > 0.\n*\n*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)\n*          On exit, if INFO = 0, WORK(1) returns the optimal\n*          (and minimal) LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK. LWORK >= max(1,18*N)\n*          if JOBZ = 'V', and LWORK >= max(1,12*N) if JOBZ = 'N'.\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  IWORK   (workspace/output) INTEGER array, dimension (LIWORK)\n*          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.\n*\n*  LIWORK  (input) INTEGER\n*          The dimension of the array IWORK.  LIWORK >= max(1,10*N)\n*          if the eigenvectors are desired, and LIWORK >= max(1,8*N)\n*          if only the eigenvalues are to be computed.\n*          If LIWORK = -1, then a workspace query is assumed; the\n*          routine only calculates the optimal size of the IWORK array,\n*          returns this value as the first entry of the IWORK array, and\n*          no error message related to LIWORK is issued by XERBLA.\n*\n*  INFO    (output) INTEGER\n*          On exit, INFO\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*          > 0:  if INFO = 1X, internal error in DLARRE,\n*                if INFO = 2X, internal error in DLARRV.\n*                Here, the digit X = ABS( IINFO ) < 10, where IINFO is\n*                the nonzero error code returned by DLARRE or\n*                DLARRV, respectively.\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Inderjit Dhillon, IBM Almaden, USA\n*     Osni Marques, LBNL/NERSC, USA\n*     Christof Voemel, LBNL/NERSC, USA\n*\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      LOGICAL TRYRAC\n*     ..\n*     .. External Subroutines ..\n      EXTERNAL DSTEMR\n*     ..\n\n");
    return Qnil;
  }
  if (argc != 11)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 11)", argc);
  rb_jobz = argv[0];
  rb_range = argv[1];
  rb_d = argv[2];
  rb_e = argv[3];
  rb_vl = argv[4];
  rb_vu = argv[5];
  rb_il = argv[6];
  rb_iu = argv[7];
  rb_abstol = argv[8];
  rb_lwork = argv[9];
  rb_liwork = argv[10];

  jobz = StringValueCStr(rb_jobz)[0];
  range = StringValueCStr(rb_range)[0];
  vl = NUM2DBL(rb_vl);
  vu = NUM2DBL(rb_vu);
  il = NUM2INT(rb_il);
  iu = NUM2INT(rb_iu);
  abstol = NUM2DBL(rb_abstol);
  lwork = NUM2INT(rb_lwork);
  liwork = NUM2INT(rb_liwork);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (3th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (3th argument) must be %d", 1);
  n = NA_SHAPE0(rb_d);
  if (NA_TYPE(rb_d) != NA_DFLOAT)
    rb_d = na_change_type(rb_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rb_d, doublereal*);
  if (!NA_IsNArray(rb_e))
    rb_raise(rb_eArgError, "e (4th argument) must be NArray");
  if (NA_RANK(rb_e) != 1)
    rb_raise(rb_eArgError, "rank of e (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_e) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of e must be the same as shape 0 of d");
  if (NA_TYPE(rb_e) != NA_DFLOAT)
    rb_e = na_change_type(rb_e, NA_DFLOAT);
  e = NA_PTR_TYPE(rb_e, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_w = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  w = NA_PTR_TYPE(rb_w, doublereal*);
  ldz = lsame_(&jobz,"V") ? MAX(1,n) : 1;
  m = lsame_(&range,"A") ? n : lsame_(&range,"I") ? iu-il+1 : 0;
  {
    int shape[2];
    shape[0] = ldz;
    shape[1] = MAX(1,m);
    rb_z = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  z = NA_PTR_TYPE(rb_z, doublereal*);
  {
    int shape[1];
    shape[0] = 2*MAX(1,m);
    rb_isuppz = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  isuppz = NA_PTR_TYPE(rb_isuppz, integer*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rb_work = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rb_work, doublereal*);
  {
    int shape[1];
    shape[0] = MAX(1,liwork);
    rb_iwork = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  iwork = NA_PTR_TYPE(rb_iwork, integer*);
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

  dstegr_(&jobz, &range, &n, d, e, &vl, &vu, &il, &iu, &abstol, &m, w, z, &ldz, isuppz, work, &lwork, iwork, &liwork, &info);

  rb_m = INT2NUM(m);
  rb_info = INT2NUM(info);
  return rb_ary_new3(9, rb_m, rb_w, rb_z, rb_isuppz, rb_work, rb_iwork, rb_info, rb_d, rb_e);
}

void
init_lapack_dstegr(VALUE mLapack){
  rb_define_module_function(mLapack, "dstegr", rb_dstegr, -1);
}
