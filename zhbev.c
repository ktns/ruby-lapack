#include "rb_lapack.h"

static VALUE
rb_zhbev(int argc, VALUE *argv, VALUE self){
  VALUE rb_jobz;
  char jobz; 
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_kd;
  integer kd; 
  VALUE rb_ab;
  doublecomplex *ab; 
  VALUE rb_w;
  doublereal *w; 
  VALUE rb_z;
  doublecomplex *z; 
  VALUE rb_info;
  integer info; 
  VALUE rb_ab_out__;
  doublecomplex *ab_out__;
  doublecomplex *work;
  doublereal *rwork;

  integer ldab;
  integer n;
  integer ldz;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  w, z, info, ab = NumRu::Lapack.zhbev( jobz, uplo, kd, ab)\n    or\n  NumRu::Lapack.zhbev  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZHBEV( JOBZ, UPLO, N, KD, AB, LDAB, W, Z, LDZ, WORK, RWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  ZHBEV computes all the eigenvalues and, optionally, eigenvectors of\n*  a complex Hermitian band matrix A.\n*\n\n*  Arguments\n*  =========\n*\n*  JOBZ    (input) CHARACTER*1\n*          = 'N':  Compute eigenvalues only;\n*          = 'V':  Compute eigenvalues and eigenvectors.\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U':  Upper triangle of A is stored;\n*          = 'L':  Lower triangle of A is stored.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  KD      (input) INTEGER\n*          The number of superdiagonals of the matrix A if UPLO = 'U',\n*          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.\n*\n*  AB      (input/output) COMPLEX*16 array, dimension (LDAB, N)\n*          On entry, the upper or lower triangle of the Hermitian band\n*          matrix A, stored in the first KD+1 rows of the array.  The\n*          j-th column of A is stored in the j-th column of the array AB\n*          as follows:\n*          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;\n*          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).\n*\n*          On exit, AB is overwritten by values generated during the\n*          reduction to tridiagonal form.  If UPLO = 'U', the first\n*          superdiagonal and the diagonal of the tridiagonal matrix T\n*          are returned in rows KD and KD+1 of AB, and if UPLO = 'L',\n*          the diagonal and first subdiagonal of T are returned in the\n*          first two rows of AB.\n*\n*  LDAB    (input) INTEGER\n*          The leading dimension of the array AB.  LDAB >= KD + 1.\n*\n*  W       (output) DOUBLE PRECISION array, dimension (N)\n*          If INFO = 0, the eigenvalues in ascending order.\n*\n*  Z       (output) COMPLEX*16 array, dimension (LDZ, N)\n*          If JOBZ = 'V', then if INFO = 0, Z contains the orthonormal\n*          eigenvectors of the matrix A, with the i-th column of Z\n*          holding the eigenvector associated with W(i).\n*          If JOBZ = 'N', then Z is not referenced.\n*\n*  LDZ     (input) INTEGER\n*          The leading dimension of the array Z.  LDZ >= 1, and if\n*          JOBZ = 'V', LDZ >= max(1,N).\n*\n*  WORK    (workspace) COMPLEX*16 array, dimension (N)\n*\n*  RWORK   (workspace) DOUBLE PRECISION array, dimension (max(1,3*N-2))\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*          > 0:  if INFO = i, the algorithm failed to converge; i\n*                off-diagonal elements of an intermediate tridiagonal\n*                form did not converge to zero.\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_jobz = argv[0];
  rb_uplo = argv[1];
  rb_kd = argv[2];
  rb_ab = argv[3];

  jobz = StringValueCStr(rb_jobz)[0];
  uplo = StringValueCStr(rb_uplo)[0];
  kd = NUM2INT(rb_kd);
  if (!NA_IsNArray(rb_ab))
    rb_raise(rb_eArgError, "ab (4th argument) must be NArray");
  if (NA_RANK(rb_ab) != 2)
    rb_raise(rb_eArgError, "rank of ab (4th argument) must be %d", 2);
  ldab = NA_SHAPE0(rb_ab);
  n = NA_SHAPE1(rb_ab);
  if (NA_TYPE(rb_ab) != NA_DCOMPLEX)
    rb_ab = na_change_type(rb_ab, NA_DCOMPLEX);
  ab = NA_PTR_TYPE(rb_ab, doublecomplex*);
  {
    int shape[1];
    shape[0] = n;
    rb_w = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  w = NA_PTR_TYPE(rb_w, doublereal*);
  ldz = lsame_(&jobz,"V") ? MAX(1,n) : 1;
  {
    int shape[2];
    shape[0] = ldz;
    shape[1] = n;
    rb_z = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  z = NA_PTR_TYPE(rb_z, doublecomplex*);
  {
    int shape[2];
    shape[0] = ldab;
    shape[1] = n;
    rb_ab_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  ab_out__ = NA_PTR_TYPE(rb_ab_out__, doublecomplex*);
  MEMCPY(ab_out__, ab, doublecomplex, NA_TOTAL(rb_ab));
  rb_ab = rb_ab_out__;
  ab = ab_out__;
  work = ALLOC_N(doublecomplex, (n));
  rwork = ALLOC_N(doublereal, (MAX(1,3*n-2)));

  zhbev_(&jobz, &uplo, &n, &kd, ab, &ldab, w, z, &ldz, work, rwork, &info);

  free(work);
  free(rwork);
  rb_info = INT2NUM(info);
  return rb_ary_new3(4, rb_w, rb_z, rb_info, rb_ab);
}

void
init_lapack_zhbev(VALUE mLapack){
  rb_define_module_function(mLapack, "zhbev", rb_zhbev, -1);
}
