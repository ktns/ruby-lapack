#include "rb_lapack.h"

extern VOID dsteqr_(char *compz, integer *n, doublereal *d, doublereal *e, doublereal *z, integer *ldz, doublereal *work, integer *info);

static VALUE
rb_dsteqr(int argc, VALUE *argv, VALUE self){
  VALUE rb_compz;
  char compz; 
  VALUE rb_d;
  doublereal *d; 
  VALUE rb_e;
  doublereal *e; 
  VALUE rb_z;
  doublereal *z; 
  VALUE rb_info;
  integer info; 
  VALUE rb_d_out__;
  doublereal *d_out__;
  VALUE rb_e_out__;
  doublereal *e_out__;
  VALUE rb_z_out__;
  doublereal *z_out__;
  doublereal *work;

  integer n;
  integer ldz;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, d, e, z = NumRu::Lapack.dsteqr( compz, d, e, z)\n    or\n  NumRu::Lapack.dsteqr  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DSTEQR( COMPZ, N, D, E, Z, LDZ, WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  DSTEQR computes all eigenvalues and, optionally, eigenvectors of a\n*  symmetric tridiagonal matrix using the implicit QL or QR method.\n*  The eigenvectors of a full or band symmetric matrix can also be found\n*  if DSYTRD or DSPTRD or DSBTRD has been used to reduce this matrix to\n*  tridiagonal form.\n*\n\n*  Arguments\n*  =========\n*\n*  COMPZ   (input) CHARACTER*1\n*          = 'N':  Compute eigenvalues only.\n*          = 'V':  Compute eigenvalues and eigenvectors of the original\n*                  symmetric matrix.  On entry, Z must contain the\n*                  orthogonal matrix used to reduce the original matrix\n*                  to tridiagonal form.\n*          = 'I':  Compute eigenvalues and eigenvectors of the\n*                  tridiagonal matrix.  Z is initialized to the identity\n*                  matrix.\n*\n*  N       (input) INTEGER\n*          The order of the matrix.  N >= 0.\n*\n*  D       (input/output) DOUBLE PRECISION array, dimension (N)\n*          On entry, the diagonal elements of the tridiagonal matrix.\n*          On exit, if INFO = 0, the eigenvalues in ascending order.\n*\n*  E       (input/output) DOUBLE PRECISION array, dimension (N-1)\n*          On entry, the (n-1) subdiagonal elements of the tridiagonal\n*          matrix.\n*          On exit, E has been destroyed.\n*\n*  Z       (input/output) DOUBLE PRECISION array, dimension (LDZ, N)\n*          On entry, if  COMPZ = 'V', then Z contains the orthogonal\n*          matrix used in the reduction to tridiagonal form.\n*          On exit, if INFO = 0, then if  COMPZ = 'V', Z contains the\n*          orthonormal eigenvectors of the original symmetric matrix,\n*          and if COMPZ = 'I', Z contains the orthonormal eigenvectors\n*          of the symmetric tridiagonal matrix.\n*          If COMPZ = 'N', then Z is not referenced.\n*\n*  LDZ     (input) INTEGER\n*          The leading dimension of the array Z.  LDZ >= 1, and if\n*          eigenvectors are desired, then  LDZ >= max(1,N).\n*\n*  WORK    (workspace) DOUBLE PRECISION array, dimension (max(1,2*N-2))\n*          If COMPZ = 'N', then WORK is not referenced.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*          > 0:  the algorithm has failed to find all the eigenvalues in\n*                a total of 30*N iterations; if INFO = i, then i\n*                elements of E have not converged to zero; on exit, D\n*                and E contain the elements of a symmetric tridiagonal\n*                matrix which is orthogonally similar to the original\n*                matrix.\n*\n\n*  =====================================================================\n*\n\n");
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
  if (NA_TYPE(rb_z) != NA_DFLOAT)
    rb_z = na_change_type(rb_z, NA_DFLOAT);
  z = NA_PTR_TYPE(rb_z, doublereal*);
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
    rb_z_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  z_out__ = NA_PTR_TYPE(rb_z_out__, doublereal*);
  MEMCPY(z_out__, z, doublereal, NA_TOTAL(rb_z));
  rb_z = rb_z_out__;
  z = z_out__;
  work = ALLOC_N(doublereal, (lsame_(&compz,"N") ? 0 : MAX(1,2*n-2)));

  dsteqr_(&compz, &n, d, e, z, &ldz, work, &info);

  free(work);
  rb_info = INT2NUM(info);
  return rb_ary_new3(4, rb_info, rb_d, rb_e, rb_z);
}

void
init_lapack_dsteqr(VALUE mLapack){
  rb_define_module_function(mLapack, "dsteqr", rb_dsteqr, -1);
}
