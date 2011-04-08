#include "rb_lapack.h"

extern VOID sstev_(char *jobz, integer *n, real *d, real *e, real *z, integer *ldz, real *work, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_sstev(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_jobz;
  char jobz; 
  VALUE rblapack_d;
  real *d; 
  VALUE rblapack_e;
  real *e; 
  VALUE rblapack_z;
  real *z; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_d_out__;
  real *d_out__;
  VALUE rblapack_e_out__;
  real *e_out__;
  real *work;

  integer n;
  integer ldz;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  z, info, d, e = NumRu::Lapack.sstev( jobz, d, e, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE SSTEV( JOBZ, N, D, E, Z, LDZ, WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  SSTEV computes all eigenvalues and, optionally, eigenvectors of a\n*  real symmetric tridiagonal matrix A.\n*\n\n*  Arguments\n*  =========\n*\n*  JOBZ    (input) CHARACTER*1\n*          = 'N':  Compute eigenvalues only;\n*          = 'V':  Compute eigenvalues and eigenvectors.\n*\n*  N       (input) INTEGER\n*          The order of the matrix.  N >= 0.\n*\n*  D       (input/output) REAL array, dimension (N)\n*          On entry, the n diagonal elements of the tridiagonal matrix\n*          A.\n*          On exit, if INFO = 0, the eigenvalues in ascending order.\n*\n*  E       (input/output) REAL array, dimension (N-1)\n*          On entry, the (n-1) subdiagonal elements of the tridiagonal\n*          matrix A, stored in elements 1 to N-1 of E.\n*          On exit, the contents of E are destroyed.\n*\n*  Z       (output) REAL array, dimension (LDZ, N)\n*          If JOBZ = 'V', then if INFO = 0, Z contains the orthonormal\n*          eigenvectors of the matrix A, with the i-th column of Z\n*          holding the eigenvector associated with D(i).\n*          If JOBZ = 'N', then Z is not referenced.\n*\n*  LDZ     (input) INTEGER\n*          The leading dimension of the array Z.  LDZ >= 1, and if\n*          JOBZ = 'V', LDZ >= max(1,N).\n*\n*  WORK    (workspace) REAL array, dimension (max(1,2*N-2))\n*          If JOBZ = 'N', WORK is not referenced.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*          > 0:  if INFO = i, the algorithm failed to converge; i\n*                off-diagonal elements of E did not converge to zero.\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  z, info, d, e = NumRu::Lapack.sstev( jobz, d, e, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rblapack_jobz = argv[0];
  rblapack_d = argv[1];
  rblapack_e = argv[2];
  if (rb_options != Qnil) {
  }

  jobz = StringValueCStr(rblapack_jobz)[0];
  if (!NA_IsNArray(rblapack_d))
    rb_raise(rb_eArgError, "d (2th argument) must be NArray");
  if (NA_RANK(rblapack_d) != 1)
    rb_raise(rb_eArgError, "rank of d (2th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_d);
  if (NA_TYPE(rblapack_d) != NA_SFLOAT)
    rblapack_d = na_change_type(rblapack_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rblapack_d, real*);
  if (!NA_IsNArray(rblapack_e))
    rb_raise(rb_eArgError, "e (3th argument) must be NArray");
  if (NA_RANK(rblapack_e) != 1)
    rb_raise(rb_eArgError, "rank of e (3th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_e) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of e must be %d", n-1);
  if (NA_TYPE(rblapack_e) != NA_SFLOAT)
    rblapack_e = na_change_type(rblapack_e, NA_SFLOAT);
  e = NA_PTR_TYPE(rblapack_e, real*);
  ldz = lsame_(&jobz,"V") ? MAX(1,n) : 1;
  {
    int shape[2];
    shape[0] = ldz;
    shape[1] = n;
    rblapack_z = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  z = NA_PTR_TYPE(rblapack_z, real*);
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
    shape[0] = n-1;
    rblapack_e_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  e_out__ = NA_PTR_TYPE(rblapack_e_out__, real*);
  MEMCPY(e_out__, e, real, NA_TOTAL(rblapack_e));
  rblapack_e = rblapack_e_out__;
  e = e_out__;
  work = ALLOC_N(real, (lsame_(&jobz,"N") ? 0 : MAX(1,2*n-2)));

  sstev_(&jobz, &n, d, e, z, &ldz, work, &info);

  free(work);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(4, rblapack_z, rblapack_info, rblapack_d, rblapack_e);
}

void
init_lapack_sstev(VALUE mLapack){
  rb_define_module_function(mLapack, "sstev", rblapack_sstev, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
