#include "rb_lapack.h"

extern VOID spteqr_(char *compz, integer *n, real *d, real *e, real *z, integer *ldz, real *work, integer *info);

static VALUE
rb_spteqr(int argc, VALUE *argv, VALUE self){
  VALUE rb_compz;
  char compz; 
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
  VALUE rb_z_out__;
  real *z_out__;
  real *work;

  integer n;
  integer ldz;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, d, e, z = NumRu::Lapack.spteqr( compz, d, e, z)\n    or\n  NumRu::Lapack.spteqr  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SPTEQR( COMPZ, N, D, E, Z, LDZ, WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  SPTEQR computes all eigenvalues and, optionally, eigenvectors of a\n*  symmetric positive definite tridiagonal matrix by first factoring the\n*  matrix using SPTTRF, and then calling SBDSQR to compute the singular\n*  values of the bidiagonal factor.\n*\n*  This routine computes the eigenvalues of the positive definite\n*  tridiagonal matrix to high relative accuracy.  This means that if the\n*  eigenvalues range over many orders of magnitude in size, then the\n*  small eigenvalues and corresponding eigenvectors will be computed\n*  more accurately than, for example, with the standard QR method.\n*\n*  The eigenvectors of a full or band symmetric positive definite matrix\n*  can also be found if SSYTRD, SSPTRD, or SSBTRD has been used to\n*  reduce this matrix to tridiagonal form. (The reduction to tridiagonal\n*  form, however, may preclude the possibility of obtaining high\n*  relative accuracy in the small eigenvalues of the original matrix, if\n*  these eigenvalues range over many orders of magnitude.)\n*\n\n*  Arguments\n*  =========\n*\n*  COMPZ   (input) CHARACTER*1\n*          = 'N':  Compute eigenvalues only.\n*          = 'V':  Compute eigenvectors of original symmetric\n*                  matrix also.  Array Z contains the orthogonal\n*                  matrix used to reduce the original matrix to\n*                  tridiagonal form.\n*          = 'I':  Compute eigenvectors of tridiagonal matrix also.\n*\n*  N       (input) INTEGER\n*          The order of the matrix.  N >= 0.\n*\n*  D       (input/output) REAL array, dimension (N)\n*          On entry, the n diagonal elements of the tridiagonal\n*          matrix.\n*          On normal exit, D contains the eigenvalues, in descending\n*          order.\n*\n*  E       (input/output) REAL array, dimension (N-1)\n*          On entry, the (n-1) subdiagonal elements of the tridiagonal\n*          matrix.\n*          On exit, E has been destroyed.\n*\n*  Z       (input/output) REAL array, dimension (LDZ, N)\n*          On entry, if COMPZ = 'V', the orthogonal matrix used in the\n*          reduction to tridiagonal form.\n*          On exit, if COMPZ = 'V', the orthonormal eigenvectors of the\n*          original symmetric matrix;\n*          if COMPZ = 'I', the orthonormal eigenvectors of the\n*          tridiagonal matrix.\n*          If INFO > 0 on exit, Z contains the eigenvectors associated\n*          with only the stored eigenvalues.\n*          If  COMPZ = 'N', then Z is not referenced.\n*\n*  LDZ     (input) INTEGER\n*          The leading dimension of the array Z.  LDZ >= 1, and if\n*          COMPZ = 'V' or 'I', LDZ >= max(1,N).\n*\n*  WORK    (workspace) REAL array, dimension (4*N)\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*          > 0:  if INFO = i, and i is:\n*                <= N  the Cholesky factorization of the matrix could\n*                      not be performed because the i-th principal minor\n*                      was not positive definite.\n*                > N   the SVD algorithm failed to converge;\n*                      if INFO = N+i, i off-diagonal elements of the\n*                      bidiagonal factor did not converge to zero.\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_compz = argv[0];
  rb_d = argv[1];
  rb_e = argv[2];
  rb_z = argv[3];

  compz = StringValueCStr(rb_compz)[0];
  if (!NA_IsNArray(rb_z))
    rb_raise(rb_eArgError, "z (4th argument) must be NArray");
  if (NA_RANK(rb_z) != 2)
    rb_raise(rb_eArgError, "rank of z (4th argument) must be %d", 2);
  n = NA_SHAPE1(rb_z);
  ldz = NA_SHAPE0(rb_z);
  if (NA_TYPE(rb_z) != NA_SFLOAT)
    rb_z = na_change_type(rb_z, NA_SFLOAT);
  z = NA_PTR_TYPE(rb_z, real*);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (2th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (2th argument) must be %d", 1);
  if (NA_SHAPE0(rb_d) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of d must be the same as shape 1 of z");
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
    rb_z_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  z_out__ = NA_PTR_TYPE(rb_z_out__, real*);
  MEMCPY(z_out__, z, real, NA_TOTAL(rb_z));
  rb_z = rb_z_out__;
  z = z_out__;
  work = ALLOC_N(real, (4*n));

  spteqr_(&compz, &n, d, e, z, &ldz, work, &info);

  free(work);
  rb_info = INT2NUM(info);
  return rb_ary_new3(4, rb_info, rb_d, rb_e, rb_z);
}

void
init_lapack_spteqr(VALUE mLapack){
  rb_define_module_function(mLapack, "spteqr", rb_spteqr, -1);
}
