#include "rb_lapack.h"

extern VOID cgtcon_(char *norm, integer *n, complex *dl, complex *d, complex *du, complex *du2, integer *ipiv, real *anorm, real *rcond, complex *work, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_cgtcon(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_norm;
  char norm; 
  VALUE rblapack_dl;
  complex *dl; 
  VALUE rblapack_d;
  complex *d; 
  VALUE rblapack_du;
  complex *du; 
  VALUE rblapack_du2;
  complex *du2; 
  VALUE rblapack_ipiv;
  integer *ipiv; 
  VALUE rblapack_anorm;
  real anorm; 
  VALUE rblapack_rcond;
  real rcond; 
  VALUE rblapack_info;
  integer info; 
  complex *work;

  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  rcond, info = NumRu::Lapack.cgtcon( norm, dl, d, du, du2, ipiv, anorm, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE CGTCON( NORM, N, DL, D, DU, DU2, IPIV, ANORM, RCOND, WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  CGTCON estimates the reciprocal of the condition number of a complex\n*  tridiagonal matrix A using the LU factorization as computed by\n*  CGTTRF.\n*\n*  An estimate is obtained for norm(inv(A)), and the reciprocal of the\n*  condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))).\n*\n\n*  Arguments\n*  =========\n*\n*  NORM    (input) CHARACTER*1\n*          Specifies whether the 1-norm condition number or the\n*          infinity-norm condition number is required:\n*          = '1' or 'O':  1-norm;\n*          = 'I':         Infinity-norm.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  DL      (input) COMPLEX array, dimension (N-1)\n*          The (n-1) multipliers that define the matrix L from the\n*          LU factorization of A as computed by CGTTRF.\n*\n*  D       (input) COMPLEX array, dimension (N)\n*          The n diagonal elements of the upper triangular matrix U from\n*          the LU factorization of A.\n*\n*  DU      (input) COMPLEX array, dimension (N-1)\n*          The (n-1) elements of the first superdiagonal of U.\n*\n*  DU2     (input) COMPLEX array, dimension (N-2)\n*          The (n-2) elements of the second superdiagonal of U.\n*\n*  IPIV    (input) INTEGER array, dimension (N)\n*          The pivot indices; for 1 <= i <= n, row i of the matrix was\n*          interchanged with row IPIV(i).  IPIV(i) will always be either\n*          i or i+1; IPIV(i) = i indicates a row interchange was not\n*          required.\n*\n*  ANORM   (input) REAL\n*          If NORM = '1' or 'O', the 1-norm of the original matrix A.\n*          If NORM = 'I', the infinity-norm of the original matrix A.\n*\n*  RCOND   (output) REAL\n*          The reciprocal of the condition number of the matrix A,\n*          computed as RCOND = 1/(ANORM * AINVNM), where AINVNM is an\n*          estimate of the 1-norm of inv(A) computed in this routine.\n*\n*  WORK    (workspace) COMPLEX array, dimension (2*N)\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  rcond, info = NumRu::Lapack.cgtcon( norm, dl, d, du, du2, ipiv, anorm, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 7)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 7)", argc);
  rblapack_norm = argv[0];
  rblapack_dl = argv[1];
  rblapack_d = argv[2];
  rblapack_du = argv[3];
  rblapack_du2 = argv[4];
  rblapack_ipiv = argv[5];
  rblapack_anorm = argv[6];
  if (rb_options != Qnil) {
  }

  anorm = (real)NUM2DBL(rblapack_anorm);
  if (!NA_IsNArray(rblapack_ipiv))
    rb_raise(rb_eArgError, "ipiv (6th argument) must be NArray");
  if (NA_RANK(rblapack_ipiv) != 1)
    rb_raise(rb_eArgError, "rank of ipiv (6th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_ipiv);
  if (NA_TYPE(rblapack_ipiv) != NA_LINT)
    rblapack_ipiv = na_change_type(rblapack_ipiv, NA_LINT);
  ipiv = NA_PTR_TYPE(rblapack_ipiv, integer*);
  norm = StringValueCStr(rblapack_norm)[0];
  if (!NA_IsNArray(rblapack_d))
    rb_raise(rb_eArgError, "d (3th argument) must be NArray");
  if (NA_RANK(rblapack_d) != 1)
    rb_raise(rb_eArgError, "rank of d (3th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_d) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of d must be the same as shape 0 of ipiv");
  if (NA_TYPE(rblapack_d) != NA_SCOMPLEX)
    rblapack_d = na_change_type(rblapack_d, NA_SCOMPLEX);
  d = NA_PTR_TYPE(rblapack_d, complex*);
  if (!NA_IsNArray(rblapack_du))
    rb_raise(rb_eArgError, "du (4th argument) must be NArray");
  if (NA_RANK(rblapack_du) != 1)
    rb_raise(rb_eArgError, "rank of du (4th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_du) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of du must be %d", n-1);
  if (NA_TYPE(rblapack_du) != NA_SCOMPLEX)
    rblapack_du = na_change_type(rblapack_du, NA_SCOMPLEX);
  du = NA_PTR_TYPE(rblapack_du, complex*);
  if (!NA_IsNArray(rblapack_du2))
    rb_raise(rb_eArgError, "du2 (5th argument) must be NArray");
  if (NA_RANK(rblapack_du2) != 1)
    rb_raise(rb_eArgError, "rank of du2 (5th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_du2) != (n-2))
    rb_raise(rb_eRuntimeError, "shape 0 of du2 must be %d", n-2);
  if (NA_TYPE(rblapack_du2) != NA_SCOMPLEX)
    rblapack_du2 = na_change_type(rblapack_du2, NA_SCOMPLEX);
  du2 = NA_PTR_TYPE(rblapack_du2, complex*);
  if (!NA_IsNArray(rblapack_dl))
    rb_raise(rb_eArgError, "dl (2th argument) must be NArray");
  if (NA_RANK(rblapack_dl) != 1)
    rb_raise(rb_eArgError, "rank of dl (2th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_dl) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of dl must be %d", n-1);
  if (NA_TYPE(rblapack_dl) != NA_SCOMPLEX)
    rblapack_dl = na_change_type(rblapack_dl, NA_SCOMPLEX);
  dl = NA_PTR_TYPE(rblapack_dl, complex*);
  work = ALLOC_N(complex, (2*n));

  cgtcon_(&norm, &n, dl, d, du, du2, ipiv, &anorm, &rcond, work, &info);

  free(work);
  rblapack_rcond = rb_float_new((double)rcond);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(2, rblapack_rcond, rblapack_info);
}

void
init_lapack_cgtcon(VALUE mLapack){
  rb_define_module_function(mLapack, "cgtcon", rblapack_cgtcon, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
