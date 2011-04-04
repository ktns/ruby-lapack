#include "rb_lapack.h"

extern VOID sstev_(char *jobz, integer *n, real *d, real *e, real *z, integer *ldz, real *work, integer *info);

static VALUE
rb_sstev(int argc, VALUE *argv, VALUE self){
  VALUE rb_jobz;
  char jobz; 
  VALUE rb_d;
  real *d; 
  VALUE rb_e;
  real *e; 
  VALUE rb_z;
  real *z; 
  VALUE rb_info;
  integer info; 
  VALUE rb_d_out__;
  real *d_out__;
  VALUE rb_e_out__;
  real *e_out__;
  real *work;

  integer n;
  integer ldz;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  z, info, d, e = NumRu::Lapack.sstev( jobz, d, e)\n    or\n  NumRu::Lapack.sstev  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SSTEV( JOBZ, N, D, E, Z, LDZ, WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  SSTEV computes all eigenvalues and, optionally, eigenvectors of a\n*  real symmetric tridiagonal matrix A.\n*\n\n*  Arguments\n*  =========\n*\n*  JOBZ    (input) CHARACTER*1\n*          = 'N':  Compute eigenvalues only;\n*          = 'V':  Compute eigenvalues and eigenvectors.\n*\n*  N       (input) INTEGER\n*          The order of the matrix.  N >= 0.\n*\n*  D       (input/output) REAL array, dimension (N)\n*          On entry, the n diagonal elements of the tridiagonal matrix\n*          A.\n*          On exit, if INFO = 0, the eigenvalues in ascending order.\n*\n*  E       (input/output) REAL array, dimension (N-1)\n*          On entry, the (n-1) subdiagonal elements of the tridiagonal\n*          matrix A, stored in elements 1 to N-1 of E.\n*          On exit, the contents of E are destroyed.\n*\n*  Z       (output) REAL array, dimension (LDZ, N)\n*          If JOBZ = 'V', then if INFO = 0, Z contains the orthonormal\n*          eigenvectors of the matrix A, with the i-th column of Z\n*          holding the eigenvector associated with D(i).\n*          If JOBZ = 'N', then Z is not referenced.\n*\n*  LDZ     (input) INTEGER\n*          The leading dimension of the array Z.  LDZ >= 1, and if\n*          JOBZ = 'V', LDZ >= max(1,N).\n*\n*  WORK    (workspace) REAL array, dimension (max(1,2*N-2))\n*          If JOBZ = 'N', WORK is not referenced.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*          > 0:  if INFO = i, the algorithm failed to converge; i\n*                off-diagonal elements of E did not converge to zero.\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_jobz = argv[0];
  rb_d = argv[1];
  rb_e = argv[2];

  jobz = StringValueCStr(rb_jobz)[0];
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
  work = ALLOC_N(real, (lsame_(&jobz,"N") ? 0 : MAX(1,2*n-2)));

  sstev_(&jobz, &n, d, e, z, &ldz, work, &info);

  free(work);
  rb_info = INT2NUM(info);
  return rb_ary_new3(4, rb_z, rb_info, rb_d, rb_e);
}

void
init_lapack_sstev(VALUE mLapack){
  rb_define_module_function(mLapack, "sstev", rb_sstev, -1);
}
