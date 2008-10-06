#include "rb_lapack.h"

static VALUE
rb_zptcon(int argc, VALUE *argv, VALUE self){
  VALUE rb_d;
  doublereal *d; 
  VALUE rb_e;
  doublecomplex *e; 
  VALUE rb_anorm;
  doublereal anorm; 
  VALUE rb_rcond;
  doublereal rcond; 
  VALUE rb_info;
  integer info; 
  doublereal *rwork;

  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  rcond, info = NumRu::Lapack.zptcon( d, e, anorm)\n    or\n  NumRu::Lapack.zptcon  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZPTCON( N, D, E, ANORM, RCOND, RWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  ZPTCON computes the reciprocal of the condition number (in the\n*  1-norm) of a complex Hermitian positive definite tridiagonal matrix\n*  using the factorization A = L*D*L**H or A = U**H*D*U computed by\n*  ZPTTRF.\n*\n*  Norm(inv(A)) is computed by a direct method, and the reciprocal of\n*  the condition number is computed as\n*                   RCOND = 1 / (ANORM * norm(inv(A))).\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  D       (input) DOUBLE PRECISION array, dimension (N)\n*          The n diagonal elements of the diagonal matrix D from the\n*          factorization of A, as computed by ZPTTRF.\n*\n*  E       (input) COMPLEX*16 array, dimension (N-1)\n*          The (n-1) off-diagonal elements of the unit bidiagonal factor\n*          U or L from the factorization of A, as computed by ZPTTRF.\n*\n*  ANORM   (input) DOUBLE PRECISION\n*          The 1-norm of the original matrix A.\n*\n*  RCOND   (output) DOUBLE PRECISION\n*          The reciprocal of the condition number of the matrix A,\n*          computed as RCOND = 1/(ANORM * AINVNM), where AINVNM is the\n*          1-norm of inv(A) computed in this routine.\n*\n*  RWORK   (workspace) DOUBLE PRECISION array, dimension (N)\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*\n\n*  Further Details\n*  ===============\n*\n*  The method used is described in Nicholas J. Higham, \"Efficient\n*  Algorithms for Computing the Condition Number of a Tridiagonal\n*  Matrix\", SIAM J. Sci. Stat. Comput., Vol. 7, No. 1, January 1986.\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_d = argv[0];
  rb_e = argv[1];
  rb_anorm = argv[2];

  anorm = NUM2DBL(rb_anorm);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (1th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (1th argument) must be %d", 1);
  n = NA_SHAPE0(rb_d);
  if (NA_TYPE(rb_d) != NA_DFLOAT)
    rb_d = na_change_type(rb_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rb_d, doublereal*);
  if (!NA_IsNArray(rb_e))
    rb_raise(rb_eArgError, "e (2th argument) must be NArray");
  if (NA_RANK(rb_e) != 1)
    rb_raise(rb_eArgError, "rank of e (2th argument) must be %d", 1);
  if (NA_SHAPE0(rb_e) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of e must be %d", n-1);
  if (NA_TYPE(rb_e) != NA_DCOMPLEX)
    rb_e = na_change_type(rb_e, NA_DCOMPLEX);
  e = NA_PTR_TYPE(rb_e, doublecomplex*);
  rwork = ALLOC_N(doublereal, (n));

  zptcon_(&n, d, e, &anorm, &rcond, rwork, &info);

  free(rwork);
  rb_rcond = rb_float_new((double)rcond);
  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_rcond, rb_info);
}

void
init_lapack_zptcon(VALUE mLapack){
  rb_define_module_function(mLapack, "zptcon", rb_zptcon, -1);
}
