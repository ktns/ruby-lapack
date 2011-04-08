#include "rb_lapack.h"

extern VOID stpcon_(char *norm, char *uplo, char *diag, integer *n, real *ap, real *rcond, real *work, integer *iwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_stpcon(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_norm;
  char norm; 
  VALUE rblapack_uplo;
  char uplo; 
  VALUE rblapack_diag;
  char diag; 
  VALUE rblapack_ap;
  real *ap; 
  VALUE rblapack_rcond;
  real rcond; 
  VALUE rblapack_info;
  integer info; 
  real *work;
  integer *iwork;

  integer ldap;
  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  rcond, info = NumRu::Lapack.stpcon( norm, uplo, diag, ap, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE STPCON( NORM, UPLO, DIAG, N, AP, RCOND, WORK, IWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  STPCON estimates the reciprocal of the condition number of a packed\n*  triangular matrix A, in either the 1-norm or the infinity-norm.\n*\n*  The norm of A is computed and an estimate is obtained for\n*  norm(inv(A)), then the reciprocal of the condition number is\n*  computed as\n*     RCOND = 1 / ( norm(A) * norm(inv(A)) ).\n*\n\n*  Arguments\n*  =========\n*\n*  NORM    (input) CHARACTER*1\n*          Specifies whether the 1-norm condition number or the\n*          infinity-norm condition number is required:\n*          = '1' or 'O':  1-norm;\n*          = 'I':         Infinity-norm.\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U':  A is upper triangular;\n*          = 'L':  A is lower triangular.\n*\n*  DIAG    (input) CHARACTER*1\n*          = 'N':  A is non-unit triangular;\n*          = 'U':  A is unit triangular.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  AP      (input) REAL array, dimension (N*(N+1)/2)\n*          The upper or lower triangular matrix A, packed columnwise in\n*          a linear array.  The j-th column of A is stored in the array\n*          AP as follows:\n*          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;\n*          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.\n*          If DIAG = 'U', the diagonal elements of A are not referenced\n*          and are assumed to be 1.\n*\n*  RCOND   (output) REAL\n*          The reciprocal of the condition number of the matrix A,\n*          computed as RCOND = 1/(norm(A) * norm(inv(A))).\n*\n*  WORK    (workspace) REAL array, dimension (3*N)\n*\n*  IWORK   (workspace) INTEGER array, dimension (N)\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  rcond, info = NumRu::Lapack.stpcon( norm, uplo, diag, ap, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rblapack_norm = argv[0];
  rblapack_uplo = argv[1];
  rblapack_diag = argv[2];
  rblapack_ap = argv[3];
  if (rb_options != Qnil) {
  }

  diag = StringValueCStr(rblapack_diag)[0];
  uplo = StringValueCStr(rblapack_uplo)[0];
  norm = StringValueCStr(rblapack_norm)[0];
  if (!NA_IsNArray(rblapack_ap))
    rb_raise(rb_eArgError, "ap (4th argument) must be NArray");
  if (NA_RANK(rblapack_ap) != 1)
    rb_raise(rb_eArgError, "rank of ap (4th argument) must be %d", 1);
  ldap = NA_SHAPE0(rblapack_ap);
  if (NA_TYPE(rblapack_ap) != NA_SFLOAT)
    rblapack_ap = na_change_type(rblapack_ap, NA_SFLOAT);
  ap = NA_PTR_TYPE(rblapack_ap, real*);
  n = ((int)sqrtf(ldap*8+1.0f)-1)/2;
  work = ALLOC_N(real, (3*n));
  iwork = ALLOC_N(integer, (n));

  stpcon_(&norm, &uplo, &diag, &n, ap, &rcond, work, iwork, &info);

  free(work);
  free(iwork);
  rblapack_rcond = rb_float_new((double)rcond);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(2, rblapack_rcond, rblapack_info);
}

void
init_lapack_stpcon(VALUE mLapack){
  rb_define_module_function(mLapack, "stpcon", rblapack_stpcon, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
