#include "rb_lapack.h"

extern VOID sptcon_(integer *n, real *d, real *e, real *anorm, real *rcond, real *work, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_sptcon(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_d;
  real *d; 
  VALUE rblapack_e;
  real *e; 
  VALUE rblapack_anorm;
  real anorm; 
  VALUE rblapack_rcond;
  real rcond; 
  VALUE rblapack_info;
  integer info; 
  real *work;

  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  rcond, info = NumRu::Lapack.sptcon( d, e, anorm, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE SPTCON( N, D, E, ANORM, RCOND, WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  SPTCON computes the reciprocal of the condition number (in the\n*  1-norm) of a real symmetric positive definite tridiagonal matrix\n*  using the factorization A = L*D*L**T or A = U**T*D*U computed by\n*  SPTTRF.\n*\n*  Norm(inv(A)) is computed by a direct method, and the reciprocal of\n*  the condition number is computed as\n*               RCOND = 1 / (ANORM * norm(inv(A))).\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  D       (input) REAL array, dimension (N)\n*          The n diagonal elements of the diagonal matrix D from the\n*          factorization of A, as computed by SPTTRF.\n*\n*  E       (input) REAL array, dimension (N-1)\n*          The (n-1) off-diagonal elements of the unit bidiagonal factor\n*          U or L from the factorization of A,  as computed by SPTTRF.\n*\n*  ANORM   (input) REAL\n*          The 1-norm of the original matrix A.\n*\n*  RCOND   (output) REAL\n*          The reciprocal of the condition number of the matrix A,\n*          computed as RCOND = 1/(ANORM * AINVNM), where AINVNM is the\n*          1-norm of inv(A) computed in this routine.\n*\n*  WORK    (workspace) REAL array, dimension (N)\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*\n\n*  Further Details\n*  ===============\n*\n*  The method used is described in Nicholas J. Higham, \"Efficient\n*  Algorithms for Computing the Condition Number of a Tridiagonal\n*  Matrix\", SIAM J. Sci. Stat. Comput., Vol. 7, No. 1, January 1986.\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  rcond, info = NumRu::Lapack.sptcon( d, e, anorm, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rblapack_d = argv[0];
  rblapack_e = argv[1];
  rblapack_anorm = argv[2];
  if (rb_options != Qnil) {
  }

  anorm = (real)NUM2DBL(rblapack_anorm);
  if (!NA_IsNArray(rblapack_d))
    rb_raise(rb_eArgError, "d (1th argument) must be NArray");
  if (NA_RANK(rblapack_d) != 1)
    rb_raise(rb_eArgError, "rank of d (1th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_d);
  if (NA_TYPE(rblapack_d) != NA_SFLOAT)
    rblapack_d = na_change_type(rblapack_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rblapack_d, real*);
  if (!NA_IsNArray(rblapack_e))
    rb_raise(rb_eArgError, "e (2th argument) must be NArray");
  if (NA_RANK(rblapack_e) != 1)
    rb_raise(rb_eArgError, "rank of e (2th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_e) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of e must be %d", n-1);
  if (NA_TYPE(rblapack_e) != NA_SFLOAT)
    rblapack_e = na_change_type(rblapack_e, NA_SFLOAT);
  e = NA_PTR_TYPE(rblapack_e, real*);
  work = ALLOC_N(real, (n));

  sptcon_(&n, d, e, &anorm, &rcond, work, &info);

  free(work);
  rblapack_rcond = rb_float_new((double)rcond);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(2, rblapack_rcond, rblapack_info);
}

void
init_lapack_sptcon(VALUE mLapack){
  rb_define_module_function(mLapack, "sptcon", rblapack_sptcon, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
