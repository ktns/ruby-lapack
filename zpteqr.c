#include "rb_lapack.h"

extern VOID zpteqr_(char *compz, integer *n, doublereal *d, doublereal *e, doublecomplex *z, integer *ldz, doublereal *work, integer *info);

static VALUE
rb_zpteqr(int argc, VALUE *argv, VALUE self){
  VALUE rb_compz;
  char compz; 
  VALUE rb_d;
  doublereal *d; 
  VALUE rb_e;
  doublereal *e; 
  VALUE rb_z;
  doublecomplex *z; 
  VALUE rb_info;
  integer info; 
  VALUE rb_d_out__;
  doublereal *d_out__;
  VALUE rb_e_out__;
  doublereal *e_out__;
  VALUE rb_z_out__;
  doublecomplex *z_out__;
  doublereal *work;

  integer n;
  integer ldz;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, d, e, z = NumRu::Lapack.zpteqr( compz, d, e, z)\n    or\n  NumRu::Lapack.zpteqr  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZPTEQR( COMPZ, N, D, E, Z, LDZ, WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  ZPTEQR computes all eigenvalues and, optionally, eigenvectors of a\n*  symmetric positive definite tridiagonal matrix by first factoring the\n*  matrix using DPTTRF and then calling ZBDSQR to compute the singular\n*  values of the bidiagonal factor.\n*\n*  This routine computes the eigenvalues of the positive definite\n*  tridiagonal matrix to high relative accuracy.  This means that if the\n*  eigenvalues range over many orders of magnitude in size, then the\n*  small eigenvalues and corresponding eigenvectors will be computed\n*  more accurately than, for example, with the standard QR method.\n*\n*  The eigenvectors of a full or band positive definite Hermitian matrix\n*  can also be found if ZHETRD, ZHPTRD, or ZHBTRD has been used to\n*  reduce this matrix to tridiagonal form.  (The reduction to\n*  tridiagonal form, however, may preclude the possibility of obtaining\n*  high relative accuracy in the small eigenvalues of the original\n*  matrix, if these eigenvalues range over many orders of magnitude.)\n*\n\n*  Arguments\n*  =========\n*\n*  COMPZ   (input) CHARACTER*1\n*          = 'N':  Compute eigenvalues only.\n*          = 'V':  Compute eigenvectors of original Hermitian\n*                  matrix also.  Array Z contains the unitary matrix\n*                  used to reduce the original matrix to tridiagonal\n*                  form.\n*          = 'I':  Compute eigenvectors of tridiagonal matrix also.\n*\n*  N       (input) INTEGER\n*          The order of the matrix.  N >= 0.\n*\n*  D       (input/output) DOUBLE PRECISION array, dimension (N)\n*          On entry, the n diagonal elements of the tridiagonal matrix.\n*          On normal exit, D contains the eigenvalues, in descending\n*          order.\n*\n*  E       (input/output) DOUBLE PRECISION array, dimension (N-1)\n*          On entry, the (n-1) subdiagonal elements of the tridiagonal\n*          matrix.\n*          On exit, E has been destroyed.\n*\n*  Z       (input/output) COMPLEX*16 array, dimension (LDZ, N)\n*          On entry, if COMPZ = 'V', the unitary matrix used in the\n*          reduction to tridiagonal form.\n*          On exit, if COMPZ = 'V', the orthonormal eigenvectors of the\n*          original Hermitian matrix;\n*          if COMPZ = 'I', the orthonormal eigenvectors of the\n*          tridiagonal matrix.\n*          If INFO > 0 on exit, Z contains the eigenvectors associated\n*          with only the stored eigenvalues.\n*          If  COMPZ = 'N', then Z is not referenced.\n*\n*  LDZ     (input) INTEGER\n*          The leading dimension of the array Z.  LDZ >= 1, and if\n*          COMPZ = 'V' or 'I', LDZ >= max(1,N).\n*\n*  WORK    (workspace) DOUBLE PRECISION array, dimension (4*N)\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*          > 0:  if INFO = i, and i is:\n*                <= N  the Cholesky factorization of the matrix could\n*                      not be performed because the i-th principal minor\n*                      was not positive definite.\n*                > N   the SVD algorithm failed to converge;\n*                      if INFO = N+i, i off-diagonal elements of the\n*                      bidiagonal factor did not converge to zero.\n*\n\n*  ====================================================================\n*\n\n");
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
  if (NA_TYPE(rb_z) != NA_DCOMPLEX)
    rb_z = na_change_type(rb_z, NA_DCOMPLEX);
  z = NA_PTR_TYPE(rb_z, doublecomplex*);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (2th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (2th argument) must be %d", 1);
  if (NA_SHAPE0(rb_d) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of d must be the same as shape 1 of z");
  if (NA_TYPE(rb_d) != NA_DFLOAT)
    rb_d = na_change_type(rb_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rb_d, doublereal*);
  if (!NA_IsNArray(rb_e))
    rb_raise(rb_eArgError, "e (3th argument) must be NArray");
  if (NA_RANK(rb_e) != 1)
    rb_raise(rb_eArgError, "rank of e (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_e) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of e must be %d", n-1);
  if (NA_TYPE(rb_e) != NA_DFLOAT)
    rb_e = na_change_type(rb_e, NA_DFLOAT);
  e = NA_PTR_TYPE(rb_e, doublereal*);
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
  {
    int shape[2];
    shape[0] = ldz;
    shape[1] = n;
    rb_z_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  z_out__ = NA_PTR_TYPE(rb_z_out__, doublecomplex*);
  MEMCPY(z_out__, z, doublecomplex, NA_TOTAL(rb_z));
  rb_z = rb_z_out__;
  z = z_out__;
  work = ALLOC_N(doublereal, (4*n));

  zpteqr_(&compz, &n, d, e, z, &ldz, work, &info);

  free(work);
  rb_info = INT2NUM(info);
  return rb_ary_new3(4, rb_info, rb_d, rb_e, rb_z);
}

void
init_lapack_zpteqr(VALUE mLapack){
  rb_define_module_function(mLapack, "zpteqr", rb_zpteqr, -1);
}
