#include "rb_lapack.h"

extern VOID sgtcon_(char *norm, integer *n, real *dl, real *d, real *du, real *du2, integer *ipiv, real *anorm, real *rcond, real *work, integer *iwork, integer *info);

static VALUE
rb_sgtcon(int argc, VALUE *argv, VALUE self){
  VALUE rb_norm;
  char norm; 
  VALUE rb_dl;
  real *dl; 
  VALUE rb_d;
  real *d; 
  VALUE rb_du;
  real *du; 
  VALUE rb_du2;
  real *du2; 
  VALUE rb_ipiv;
  integer *ipiv; 
  VALUE rb_anorm;
  real anorm; 
  VALUE rb_rcond;
  real rcond; 
  VALUE rb_info;
  integer info; 
  real *work;
  integer *iwork;

  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  rcond, info = NumRu::Lapack.sgtcon( norm, dl, d, du, du2, ipiv, anorm)\n    or\n  NumRu::Lapack.sgtcon  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SGTCON( NORM, N, DL, D, DU, DU2, IPIV, ANORM, RCOND, WORK, IWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  SGTCON estimates the reciprocal of the condition number of a real\n*  tridiagonal matrix A using the LU factorization as computed by\n*  SGTTRF.\n*\n*  An estimate is obtained for norm(inv(A)), and the reciprocal of the\n*  condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))).\n*\n\n*  Arguments\n*  =========\n*\n*  NORM    (input) CHARACTER*1\n*          Specifies whether the 1-norm condition number or the\n*          infinity-norm condition number is required:\n*          = '1' or 'O':  1-norm;\n*          = 'I':         Infinity-norm.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  DL      (input) REAL array, dimension (N-1)\n*          The (n-1) multipliers that define the matrix L from the\n*          LU factorization of A as computed by SGTTRF.\n*\n*  D       (input) REAL array, dimension (N)\n*          The n diagonal elements of the upper triangular matrix U from\n*          the LU factorization of A.\n*\n*  DU      (input) REAL array, dimension (N-1)\n*          The (n-1) elements of the first superdiagonal of U.\n*\n*  DU2     (input) REAL array, dimension (N-2)\n*          The (n-2) elements of the second superdiagonal of U.\n*\n*  IPIV    (input) INTEGER array, dimension (N)\n*          The pivot indices; for 1 <= i <= n, row i of the matrix was\n*          interchanged with row IPIV(i).  IPIV(i) will always be either\n*          i or i+1; IPIV(i) = i indicates a row interchange was not\n*          required.\n*\n*  ANORM   (input) REAL\n*          If NORM = '1' or 'O', the 1-norm of the original matrix A.\n*          If NORM = 'I', the infinity-norm of the original matrix A.\n*\n*  RCOND   (output) REAL\n*          The reciprocal of the condition number of the matrix A,\n*          computed as RCOND = 1/(ANORM * AINVNM), where AINVNM is an\n*          estimate of the 1-norm of inv(A) computed in this routine.\n*\n*  WORK    (workspace) REAL array, dimension (2*N)\n*\n*  IWORK   (workspace) INTEGER array, dimension (N)\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 7)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 7)", argc);
  rb_norm = argv[0];
  rb_dl = argv[1];
  rb_d = argv[2];
  rb_du = argv[3];
  rb_du2 = argv[4];
  rb_ipiv = argv[5];
  rb_anorm = argv[6];

  anorm = (real)NUM2DBL(rb_anorm);
  if (!NA_IsNArray(rb_ipiv))
    rb_raise(rb_eArgError, "ipiv (6th argument) must be NArray");
  if (NA_RANK(rb_ipiv) != 1)
    rb_raise(rb_eArgError, "rank of ipiv (6th argument) must be %d", 1);
  n = NA_SHAPE0(rb_ipiv);
  if (NA_TYPE(rb_ipiv) != NA_LINT)
    rb_ipiv = na_change_type(rb_ipiv, NA_LINT);
  ipiv = NA_PTR_TYPE(rb_ipiv, integer*);
  norm = StringValueCStr(rb_norm)[0];
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (3th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_d) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of d must be the same as shape 0 of ipiv");
  if (NA_TYPE(rb_d) != NA_SFLOAT)
    rb_d = na_change_type(rb_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rb_d, real*);
  if (!NA_IsNArray(rb_du))
    rb_raise(rb_eArgError, "du (4th argument) must be NArray");
  if (NA_RANK(rb_du) != 1)
    rb_raise(rb_eArgError, "rank of du (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_du) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of du must be %d", n-1);
  if (NA_TYPE(rb_du) != NA_SFLOAT)
    rb_du = na_change_type(rb_du, NA_SFLOAT);
  du = NA_PTR_TYPE(rb_du, real*);
  if (!NA_IsNArray(rb_du2))
    rb_raise(rb_eArgError, "du2 (5th argument) must be NArray");
  if (NA_RANK(rb_du2) != 1)
    rb_raise(rb_eArgError, "rank of du2 (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_du2) != (n-2))
    rb_raise(rb_eRuntimeError, "shape 0 of du2 must be %d", n-2);
  if (NA_TYPE(rb_du2) != NA_SFLOAT)
    rb_du2 = na_change_type(rb_du2, NA_SFLOAT);
  du2 = NA_PTR_TYPE(rb_du2, real*);
  if (!NA_IsNArray(rb_dl))
    rb_raise(rb_eArgError, "dl (2th argument) must be NArray");
  if (NA_RANK(rb_dl) != 1)
    rb_raise(rb_eArgError, "rank of dl (2th argument) must be %d", 1);
  if (NA_SHAPE0(rb_dl) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of dl must be %d", n-1);
  if (NA_TYPE(rb_dl) != NA_SFLOAT)
    rb_dl = na_change_type(rb_dl, NA_SFLOAT);
  dl = NA_PTR_TYPE(rb_dl, real*);
  work = ALLOC_N(real, (2*n));
  iwork = ALLOC_N(integer, (n));

  sgtcon_(&norm, &n, dl, d, du, du2, ipiv, &anorm, &rcond, work, iwork, &info);

  free(work);
  free(iwork);
  rb_rcond = rb_float_new((double)rcond);
  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_rcond, rb_info);
}

void
init_lapack_sgtcon(VALUE mLapack){
  rb_define_module_function(mLapack, "sgtcon", rb_sgtcon, -1);
}
