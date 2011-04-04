#include "rb_lapack.h"

extern VOID dstevd_(char *jobz, integer *n, doublereal *d, doublereal *e, doublereal *z, integer *ldz, doublereal *work, integer *lwork, integer *iwork, integer *liwork, integer *info);

static VALUE
rb_dstevd(int argc, VALUE *argv, VALUE self){
  VALUE rb_jobz;
  char jobz; 
  VALUE rb_d;
  doublereal *d; 
  VALUE rb_e;
  doublereal *e; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_liwork;
  integer liwork; 
  VALUE rb_z;
  doublereal *z; 
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
    printf("%s\n", "USAGE:\n  z, work, iwork, info, d, e = NumRu::Lapack.dstevd( jobz, d, e, lwork, liwork)\n    or\n  NumRu::Lapack.dstevd  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DSTEVD( JOBZ, N, D, E, Z, LDZ, WORK, LWORK, IWORK, LIWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  DSTEVD computes all eigenvalues and, optionally, eigenvectors of a\n*  real symmetric tridiagonal matrix. If eigenvectors are desired, it\n*  uses a divide and conquer algorithm.\n*\n*  The divide and conquer algorithm makes very mild assumptions about\n*  floating point arithmetic. It will work on machines with a guard\n*  digit in add/subtract, or on those binary machines without guard\n*  digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or\n*  Cray-2. It could conceivably fail on hexadecimal or decimal machines\n*  without guard digits, but we know of none.\n*\n\n*  Arguments\n*  =========\n*\n*  JOBZ    (input) CHARACTER*1\n*          = 'N':  Compute eigenvalues only;\n*          = 'V':  Compute eigenvalues and eigenvectors.\n*\n*  N       (input) INTEGER\n*          The order of the matrix.  N >= 0.\n*\n*  D       (input/output) DOUBLE PRECISION array, dimension (N)\n*          On entry, the n diagonal elements of the tridiagonal matrix\n*          A.\n*          On exit, if INFO = 0, the eigenvalues in ascending order.\n*\n*  E       (input/output) DOUBLE PRECISION array, dimension (N-1)\n*          On entry, the (n-1) subdiagonal elements of the tridiagonal\n*          matrix A, stored in elements 1 to N-1 of E.\n*          On exit, the contents of E are destroyed.\n*\n*  Z       (output) DOUBLE PRECISION array, dimension (LDZ, N)\n*          If JOBZ = 'V', then if INFO = 0, Z contains the orthonormal\n*          eigenvectors of the matrix A, with the i-th column of Z\n*          holding the eigenvector associated with D(i).\n*          If JOBZ = 'N', then Z is not referenced.\n*\n*  LDZ     (input) INTEGER\n*          The leading dimension of the array Z.  LDZ >= 1, and if\n*          JOBZ = 'V', LDZ >= max(1,N).\n*\n*  WORK    (workspace/output) DOUBLE PRECISION array,\n*                                         dimension (LWORK)\n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK.\n*          If JOBZ  = 'N' or N <= 1 then LWORK must be at least 1.\n*          If JOBZ  = 'V' and N > 1 then LWORK must be at least\n*                         ( 1 + 4*N + N**2 ).\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal sizes of the WORK and IWORK\n*          arrays, returns these values as the first entries of the WORK\n*          and IWORK arrays, and no error message related to LWORK or\n*          LIWORK is issued by XERBLA.\n*\n*  IWORK   (workspace/output) INTEGER array, dimension (MAX(1,LIWORK))\n*          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.\n*\n*  LIWORK  (input) INTEGER\n*          The dimension of the array IWORK.\n*          If JOBZ  = 'N' or N <= 1 then LIWORK must be at least 1.\n*          If JOBZ  = 'V' and N > 1 then LIWORK must be at least 3+5*N.\n*\n*          If LIWORK = -1, then a workspace query is assumed; the\n*          routine only calculates the optimal sizes of the WORK and\n*          IWORK arrays, returns these values as the first entries of\n*          the WORK and IWORK arrays, and no error message related to\n*          LWORK or LIWORK is issued by XERBLA.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*          > 0:  if INFO = i, the algorithm failed to converge; i\n*                off-diagonal elements of E did not converge to zero.\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_jobz = argv[0];
  rb_d = argv[1];
  rb_e = argv[2];
  rb_lwork = argv[3];
  rb_liwork = argv[4];

  jobz = StringValueCStr(rb_jobz)[0];
  liwork = NUM2INT(rb_liwork);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (2th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (2th argument) must be %d", 1);
  n = NA_SHAPE0(rb_d);
  if (NA_TYPE(rb_d) != NA_DFLOAT)
    rb_d = na_change_type(rb_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rb_d, doublereal*);
  lwork = NUM2INT(rb_lwork);
  if (!NA_IsNArray(rb_e))
    rb_raise(rb_eArgError, "e (3th argument) must be NArray");
  if (NA_RANK(rb_e) != 1)
    rb_raise(rb_eArgError, "rank of e (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_e) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of e must be %d", n-1);
  if (NA_TYPE(rb_e) != NA_DFLOAT)
    rb_e = na_change_type(rb_e, NA_DFLOAT);
  e = NA_PTR_TYPE(rb_e, doublereal*);
  ldz = lsame_(&jobz,"V") ? MAX(1,n) : 1;
  {
    int shape[2];
    shape[0] = ldz;
    shape[1] = n;
    rb_z = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  z = NA_PTR_TYPE(rb_z, doublereal*);
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
    shape[0] = n-1;
    rb_e_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  e_out__ = NA_PTR_TYPE(rb_e_out__, doublereal*);
  MEMCPY(e_out__, e, doublereal, NA_TOTAL(rb_e));
  rb_e = rb_e_out__;
  e = e_out__;

  dstevd_(&jobz, &n, d, e, z, &ldz, work, &lwork, iwork, &liwork, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(6, rb_z, rb_work, rb_iwork, rb_info, rb_d, rb_e);
}

void
init_lapack_dstevd(VALUE mLapack){
  rb_define_module_function(mLapack, "dstevd", rb_dstevd, -1);
}
