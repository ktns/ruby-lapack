#include "rb_lapack.h"

static VALUE
rb_sspev(int argc, VALUE *argv, VALUE self){
  VALUE rb_jobz;
  char jobz; 
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_ap;
  real *ap; 
  VALUE rb_w;
  real *w; 
  VALUE rb_z;
  real *z; 
  VALUE rb_info;
  integer info; 
  VALUE rb_ap_out__;
  real *ap_out__;
  real *work;

  integer ldap;
  integer n;
  integer ldz;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  w, z, info, ap = NumRu::Lapack.sspev( jobz, uplo, ap)\n    or\n  NumRu::Lapack.sspev  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SSPEV( JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  SSPEV computes all the eigenvalues and, optionally, eigenvectors of a\n*  real symmetric matrix A in packed storage.\n*\n\n*  Arguments\n*  =========\n*\n*  JOBZ    (input) CHARACTER*1\n*          = 'N':  Compute eigenvalues only;\n*          = 'V':  Compute eigenvalues and eigenvectors.\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U':  Upper triangle of A is stored;\n*          = 'L':  Lower triangle of A is stored.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  AP      (input/output) REAL array, dimension (N*(N+1)/2)\n*          On entry, the upper or lower triangle of the symmetric matrix\n*          A, packed columnwise in a linear array.  The j-th column of A\n*          is stored in the array AP as follows:\n*          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;\n*          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.\n*\n*          On exit, AP is overwritten by values generated during the\n*          reduction to tridiagonal form.  If UPLO = 'U', the diagonal\n*          and first superdiagonal of the tridiagonal matrix T overwrite\n*          the corresponding elements of A, and if UPLO = 'L', the\n*          diagonal and first subdiagonal of T overwrite the\n*          corresponding elements of A.\n*\n*  W       (output) REAL array, dimension (N)\n*          If INFO = 0, the eigenvalues in ascending order.\n*\n*  Z       (output) REAL array, dimension (LDZ, N)\n*          If JOBZ = 'V', then if INFO = 0, Z contains the orthonormal\n*          eigenvectors of the matrix A, with the i-th column of Z\n*          holding the eigenvector associated with W(i).\n*          If JOBZ = 'N', then Z is not referenced.\n*\n*  LDZ     (input) INTEGER\n*          The leading dimension of the array Z.  LDZ >= 1, and if\n*          JOBZ = 'V', LDZ >= max(1,N).\n*\n*  WORK    (workspace) REAL array, dimension (3*N)\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*          > 0:  if INFO = i, the algorithm failed to converge; i\n*                off-diagonal elements of an intermediate tridiagonal\n*                form did not converge to zero.\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_jobz = argv[0];
  rb_uplo = argv[1];
  rb_ap = argv[2];

  jobz = StringValueCStr(rb_jobz)[0];
  uplo = StringValueCStr(rb_uplo)[0];
  if (!NA_IsNArray(rb_ap))
    rb_raise(rb_eArgError, "ap (3th argument) must be NArray");
  if (NA_RANK(rb_ap) != 1)
    rb_raise(rb_eArgError, "rank of ap (3th argument) must be %d", 1);
  ldap = NA_SHAPE0(rb_ap);
  if (NA_TYPE(rb_ap) != NA_SFLOAT)
    rb_ap = na_change_type(rb_ap, NA_SFLOAT);
  ap = NA_PTR_TYPE(rb_ap, real*);
  n = (int)(sqrt((double)8*ldap+1)-1)/2;
  {
    int shape[1];
    shape[0] = n;
    rb_w = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  w = NA_PTR_TYPE(rb_w, real*);
  ldz = lsame_(&jobz,"V") ? MAX(1,n) : 1;
  {
    int shape[2];
    shape[0] = ldz;
    shape[1] = n;
    rb_z = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  z = NA_PTR_TYPE(rb_z, real*);
  {
    int shape[1];
    shape[0] = ldap;
    rb_ap_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  ap_out__ = NA_PTR_TYPE(rb_ap_out__, real*);
  MEMCPY(ap_out__, ap, real, NA_TOTAL(rb_ap));
  rb_ap = rb_ap_out__;
  ap = ap_out__;
  work = ALLOC_N(real, (3*n));

  sspev_(&jobz, &uplo, &n, ap, w, z, &ldz, work, &info);

  free(work);
  rb_info = INT2NUM(info);
  return rb_ary_new3(4, rb_w, rb_z, rb_info, rb_ap);
}

void
init_lapack_sspev(VALUE mLapack){
  rb_define_module_function(mLapack, "sspev", rb_sspev, -1);
}
