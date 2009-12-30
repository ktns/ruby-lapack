#include "rb_lapack.h"

static VALUE
rb_cstedc(int argc, VALUE *argv, VALUE self){
  VALUE rb_compz;
  char compz; 
  VALUE rb_d;
  real *d; 
  VALUE rb_e;
  real *e; 
  VALUE rb_z;
  complex *z; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_lrwork;
  integer lrwork; 
  VALUE rb_liwork;
  integer liwork; 
  VALUE rb_work;
  complex *work; 
  VALUE rb_rwork;
  real *rwork; 
  VALUE rb_iwork;
  integer *iwork; 
  VALUE rb_info;
  integer info; 
  VALUE rb_d_out__;
  real *d_out__;
  VALUE rb_e_out__;
  real *e_out__;
  VALUE rb_z_out__;
  complex *z_out__;

  integer n;
  integer ldz;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  work, rwork, iwork, info, d, e, z = NumRu::Lapack.cstedc( compz, d, e, z, lwork, lrwork, liwork)\n    or\n  NumRu::Lapack.cstedc  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE CSTEDC( COMPZ, N, D, E, Z, LDZ, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  CSTEDC computes all eigenvalues and, optionally, eigenvectors of a\n*  symmetric tridiagonal matrix using the divide and conquer method.\n*  The eigenvectors of a full or band complex Hermitian matrix can also\n*  be found if CHETRD or CHPTRD or CHBTRD has been used to reduce this\n*  matrix to tridiagonal form.\n*\n*  This code makes very mild assumptions about floating point\n*  arithmetic. It will work on machines with a guard digit in\n*  add/subtract, or on those binary machines without guard digits\n*  which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2.\n*  It could conceivably fail on hexadecimal or decimal machines\n*  without guard digits, but we know of none.  See SLAED3 for details.\n*\n\n*  Arguments\n*  =========\n*\n*  COMPZ   (input) CHARACTER*1\n*          = 'N':  Compute eigenvalues only.\n*          = 'I':  Compute eigenvectors of tridiagonal matrix also.\n*          = 'V':  Compute eigenvectors of original Hermitian matrix\n*                  also.  On entry, Z contains the unitary matrix used\n*                  to reduce the original matrix to tridiagonal form.\n*\n*  N       (input) INTEGER\n*          The dimension of the symmetric tridiagonal matrix.  N >= 0.\n*\n*  D       (input/output) REAL array, dimension (N)\n*          On entry, the diagonal elements of the tridiagonal matrix.\n*          On exit, if INFO = 0, the eigenvalues in ascending order.\n*\n*  E       (input/output) REAL array, dimension (N-1)\n*          On entry, the subdiagonal elements of the tridiagonal matrix.\n*          On exit, E has been destroyed.\n*\n*  Z       (input/output) COMPLEX array, dimension (LDZ,N)\n*          On entry, if COMPZ = 'V', then Z contains the unitary\n*          matrix used in the reduction to tridiagonal form.\n*          On exit, if INFO = 0, then if COMPZ = 'V', Z contains the\n*          orthonormal eigenvectors of the original Hermitian matrix,\n*          and if COMPZ = 'I', Z contains the orthonormal eigenvectors\n*          of the symmetric tridiagonal matrix.\n*          If  COMPZ = 'N', then Z is not referenced.\n*\n*  LDZ     (input) INTEGER\n*          The leading dimension of the array Z.  LDZ >= 1.\n*          If eigenvectors are desired, then LDZ >= max(1,N).\n*\n*  WORK    (workspace/output) COMPLEX    array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK.\n*          If COMPZ = 'N' or 'I', or N <= 1, LWORK must be at least 1.\n*          If COMPZ = 'V' and N > 1, LWORK must be at least N*N.\n*          Note that for COMPZ = 'V', then if N is less than or\n*          equal to the minimum divide size, usually 25, then LWORK need\n*          only be 1.\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal sizes of the WORK, RWORK and\n*          IWORK arrays, returns these values as the first entries of\n*          the WORK, RWORK and IWORK arrays, and no error message\n*          related to LWORK or LRWORK or LIWORK is issued by XERBLA.\n*\n*  RWORK   (workspace/output) REAL array, dimension (MAX(1,LRWORK))\n*          On exit, if INFO = 0, RWORK(1) returns the optimal LRWORK.\n*\n*  LRWORK  (input) INTEGER\n*          The dimension of the array RWORK.\n*          If COMPZ = 'N' or N <= 1, LRWORK must be at least 1.\n*          If COMPZ = 'V' and N > 1, LRWORK must be at least\n*                         1 + 3*N + 2*N*lg N + 3*N**2 ,\n*                         where lg( N ) = smallest integer k such\n*                         that 2**k >= N.\n*          If COMPZ = 'I' and N > 1, LRWORK must be at least\n*                         1 + 4*N + 2*N**2 .\n*          Note that for COMPZ = 'I' or 'V', then if N is less than or\n*          equal to the minimum divide size, usually 25, then LRWORK\n*          need only be max(1,2*(N-1)).\n*\n*          If LRWORK = -1, then a workspace query is assumed; the\n*          routine only calculates the optimal sizes of the WORK, RWORK\n*          and IWORK arrays, returns these values as the first entries\n*          of the WORK, RWORK and IWORK arrays, and no error message\n*          related to LWORK or LRWORK or LIWORK is issued by XERBLA.\n*\n*  IWORK   (workspace/output) INTEGER array, dimension (MAX(1,LIWORK))\n*          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.\n*\n*  LIWORK  (input) INTEGER\n*          The dimension of the array IWORK.\n*          If COMPZ = 'N' or N <= 1, LIWORK must be at least 1.\n*          If COMPZ = 'V' or N > 1,  LIWORK must be at least\n*                                    6 + 6*N + 5*N*lg N.\n*          If COMPZ = 'I' or N > 1,  LIWORK must be at least\n*                                    3 + 5*N .\n*          Note that for COMPZ = 'I' or 'V', then if N is less than or\n*          equal to the minimum divide size, usually 25, then LIWORK\n*          need only be 1.\n*\n*          If LIWORK = -1, then a workspace query is assumed; the\n*          routine only calculates the optimal sizes of the WORK, RWORK\n*          and IWORK arrays, returns these values as the first entries\n*          of the WORK, RWORK and IWORK arrays, and no error message\n*          related to LWORK or LRWORK or LIWORK is issued by XERBLA.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*          > 0:  The algorithm failed to compute an eigenvalue while\n*                working on the submatrix lying in rows and columns\n*                INFO/(N+1) through mod(INFO,N+1).\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Jeff Rutter, Computer Science Division, University of California\n*     at Berkeley, USA\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 7)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 7)", argc);
  rb_compz = argv[0];
  rb_d = argv[1];
  rb_e = argv[2];
  rb_z = argv[3];
  rb_lwork = argv[4];
  rb_lrwork = argv[5];
  rb_liwork = argv[6];

  compz = StringValueCStr(rb_compz)[0];
  lwork = NUM2INT(rb_lwork);
  lrwork = NUM2INT(rb_lrwork);
  liwork = NUM2INT(rb_liwork);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (2th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (2th argument) must be %d", 1);
  n = NA_SHAPE0(rb_d);
  if (NA_TYPE(rb_d) != NA_SFLOAT)
    rb_d = na_change_type(rb_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rb_d, real*);
  if (!NA_IsNArray(rb_e))
    rb_raise(rb_eArgError, "e (3th argument) must be NArray");
  if (NA_RANK(rb_e) != 1)
    rb_raise(rb_eArgError, "rank of e (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_e) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of e must be %d", n-1);
  if (NA_TYPE(rb_e) != NA_SFLOAT)
    rb_e = na_change_type(rb_e, NA_SFLOAT);
  e = NA_PTR_TYPE(rb_e, real*);
  if (!NA_IsNArray(rb_z))
    rb_raise(rb_eArgError, "z (4th argument) must be NArray");
  if (NA_RANK(rb_z) != 2)
    rb_raise(rb_eArgError, "rank of z (4th argument) must be %d", 2);
  ldz = NA_SHAPE0(rb_z);
  if (NA_SHAPE1(rb_z) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of z must be the same as shape 0 of d");
  if (NA_TYPE(rb_z) != NA_SCOMPLEX)
    rb_z = na_change_type(rb_z, NA_SCOMPLEX);
  z = NA_PTR_TYPE(rb_z, complex*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rb_work = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rb_work, complex*);
  {
    int shape[1];
    shape[0] = MAX(1,lrwork);
    rb_rwork = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  rwork = NA_PTR_TYPE(rb_rwork, real*);
  {
    int shape[1];
    shape[0] = DIM_LEN(MAX(1,liwork));
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
    shape[0] = n-1;
    rb_e_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  e_out__ = NA_PTR_TYPE(rb_e_out__, real*);
  MEMCPY(e_out__, e, real, NA_TOTAL(rb_e));
  rb_e = rb_e_out__;
  e = e_out__;
  {
    int shape[2];
    shape[0] = ldz;
    shape[1] = n;
    rb_z_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  z_out__ = NA_PTR_TYPE(rb_z_out__, complex*);
  MEMCPY(z_out__, z, complex, NA_TOTAL(rb_z));
  rb_z = rb_z_out__;
  z = z_out__;

  cstedc_(&compz, &n, d, e, z, &ldz, work, &lwork, rwork, &lrwork, iwork, &liwork, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(7, rb_work, rb_rwork, rb_iwork, rb_info, rb_d, rb_e, rb_z);
}

void
init_lapack_cstedc(VALUE mLapack){
  rb_define_module_function(mLapack, "cstedc", rb_cstedc, -1);
}
