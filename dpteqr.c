#include "rb_lapack.h"

extern VOID dpteqr_(char *compz, integer *n, doublereal *d, doublereal *e, doublereal *z, integer *ldz, doublereal *work, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dpteqr(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_compz;
  char compz; 
  VALUE rblapack_d;
  doublereal *d; 
  VALUE rblapack_e;
  doublereal *e; 
  VALUE rblapack_z;
  doublereal *z; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_d_out__;
  doublereal *d_out__;
  VALUE rblapack_e_out__;
  doublereal *e_out__;
  VALUE rblapack_z_out__;
  doublereal *z_out__;
  doublereal *work;

  integer n;
  integer ldz;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, d, e, z = NumRu::Lapack.dpteqr( compz, d, e, z, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DPTEQR( COMPZ, N, D, E, Z, LDZ, WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  DPTEQR computes all eigenvalues and, optionally, eigenvectors of a\n*  symmetric positive definite tridiagonal matrix by first factoring the\n*  matrix using DPTTRF, and then calling DBDSQR to compute the singular\n*  values of the bidiagonal factor.\n*\n*  This routine computes the eigenvalues of the positive definite\n*  tridiagonal matrix to high relative accuracy.  This means that if the\n*  eigenvalues range over many orders of magnitude in size, then the\n*  small eigenvalues and corresponding eigenvectors will be computed\n*  more accurately than, for example, with the standard QR method.\n*\n*  The eigenvectors of a full or band symmetric positive definite matrix\n*  can also be found if DSYTRD, DSPTRD, or DSBTRD has been used to\n*  reduce this matrix to tridiagonal form. (The reduction to tridiagonal\n*  form, however, may preclude the possibility of obtaining high\n*  relative accuracy in the small eigenvalues of the original matrix, if\n*  these eigenvalues range over many orders of magnitude.)\n*\n\n*  Arguments\n*  =========\n*\n*  COMPZ   (input) CHARACTER*1\n*          = 'N':  Compute eigenvalues only.\n*          = 'V':  Compute eigenvectors of original symmetric\n*                  matrix also.  Array Z contains the orthogonal\n*                  matrix used to reduce the original matrix to\n*                  tridiagonal form.\n*          = 'I':  Compute eigenvectors of tridiagonal matrix also.\n*\n*  N       (input) INTEGER\n*          The order of the matrix.  N >= 0.\n*\n*  D       (input/output) DOUBLE PRECISION array, dimension (N)\n*          On entry, the n diagonal elements of the tridiagonal\n*          matrix.\n*          On normal exit, D contains the eigenvalues, in descending\n*          order.\n*\n*  E       (input/output) DOUBLE PRECISION array, dimension (N-1)\n*          On entry, the (n-1) subdiagonal elements of the tridiagonal\n*          matrix.\n*          On exit, E has been destroyed.\n*\n*  Z       (input/output) DOUBLE PRECISION array, dimension (LDZ, N)\n*          On entry, if COMPZ = 'V', the orthogonal matrix used in the\n*          reduction to tridiagonal form.\n*          On exit, if COMPZ = 'V', the orthonormal eigenvectors of the\n*          original symmetric matrix;\n*          if COMPZ = 'I', the orthonormal eigenvectors of the\n*          tridiagonal matrix.\n*          If INFO > 0 on exit, Z contains the eigenvectors associated\n*          with only the stored eigenvalues.\n*          If  COMPZ = 'N', then Z is not referenced.\n*\n*  LDZ     (input) INTEGER\n*          The leading dimension of the array Z.  LDZ >= 1, and if\n*          COMPZ = 'V' or 'I', LDZ >= max(1,N).\n*\n*  WORK    (workspace) DOUBLE PRECISION array, dimension (4*N)\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*          > 0:  if INFO = i, and i is:\n*                <= N  the Cholesky factorization of the matrix could\n*                      not be performed because the i-th principal minor\n*                      was not positive definite.\n*                > N   the SVD algorithm failed to converge;\n*                      if INFO = N+i, i off-diagonal elements of the\n*                      bidiagonal factor did not converge to zero.\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, d, e, z = NumRu::Lapack.dpteqr( compz, d, e, z, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rblapack_compz = argv[0];
  rblapack_d = argv[1];
  rblapack_e = argv[2];
  rblapack_z = argv[3];
  if (rb_options != Qnil) {
  }

  compz = StringValueCStr(rblapack_compz)[0];
  if (!NA_IsNArray(rblapack_z))
    rb_raise(rb_eArgError, "z (4th argument) must be NArray");
  if (NA_RANK(rblapack_z) != 2)
    rb_raise(rb_eArgError, "rank of z (4th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_z);
  ldz = NA_SHAPE0(rblapack_z);
  if (NA_TYPE(rblapack_z) != NA_DFLOAT)
    rblapack_z = na_change_type(rblapack_z, NA_DFLOAT);
  z = NA_PTR_TYPE(rblapack_z, doublereal*);
  if (!NA_IsNArray(rblapack_d))
    rb_raise(rb_eArgError, "d (2th argument) must be NArray");
  if (NA_RANK(rblapack_d) != 1)
    rb_raise(rb_eArgError, "rank of d (2th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_d) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of d must be the same as shape 1 of z");
  if (NA_TYPE(rblapack_d) != NA_DFLOAT)
    rblapack_d = na_change_type(rblapack_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rblapack_d, doublereal*);
  if (!NA_IsNArray(rblapack_e))
    rb_raise(rb_eArgError, "e (3th argument) must be NArray");
  if (NA_RANK(rblapack_e) != 1)
    rb_raise(rb_eArgError, "rank of e (3th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_e) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of e must be %d", n-1);
  if (NA_TYPE(rblapack_e) != NA_DFLOAT)
    rblapack_e = na_change_type(rblapack_e, NA_DFLOAT);
  e = NA_PTR_TYPE(rblapack_e, doublereal*);
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
    shape[0] = n-1;
    rblapack_e_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  e_out__ = NA_PTR_TYPE(rblapack_e_out__, doublereal*);
  MEMCPY(e_out__, e, doublereal, NA_TOTAL(rblapack_e));
  rblapack_e = rblapack_e_out__;
  e = e_out__;
  {
    int shape[2];
    shape[0] = ldz;
    shape[1] = n;
    rblapack_z_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  z_out__ = NA_PTR_TYPE(rblapack_z_out__, doublereal*);
  MEMCPY(z_out__, z, doublereal, NA_TOTAL(rblapack_z));
  rblapack_z = rblapack_z_out__;
  z = z_out__;
  work = ALLOC_N(doublereal, (4*n));

  dpteqr_(&compz, &n, d, e, z, &ldz, work, &info);

  free(work);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(4, rblapack_info, rblapack_d, rblapack_e, rblapack_z);
}

void
init_lapack_dpteqr(VALUE mLapack){
  rb_define_module_function(mLapack, "dpteqr", rblapack_dpteqr, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
