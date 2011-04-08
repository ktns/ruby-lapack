#include "rb_lapack.h"

extern VOID cspcon_(char *uplo, integer *n, complex *ap, integer *ipiv, real *anorm, real *rcond, complex *work, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_cspcon(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_uplo;
  char uplo; 
  VALUE rblapack_ap;
  complex *ap; 
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
      printf("%s\n", "USAGE:\n  rcond, info = NumRu::Lapack.cspcon( uplo, ap, ipiv, anorm, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE CSPCON( UPLO, N, AP, IPIV, ANORM, RCOND, WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  CSPCON estimates the reciprocal of the condition number (in the\n*  1-norm) of a complex symmetric packed matrix A using the\n*  factorization A = U*D*U**T or A = L*D*L**T computed by CSPTRF.\n*\n*  An estimate is obtained for norm(inv(A)), and the reciprocal of the\n*  condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))).\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO    (input) CHARACTER*1\n*          Specifies whether the details of the factorization are stored\n*          as an upper or lower triangular matrix.\n*          = 'U':  Upper triangular, form is A = U*D*U**T;\n*          = 'L':  Lower triangular, form is A = L*D*L**T.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  AP      (input) COMPLEX array, dimension (N*(N+1)/2)\n*          The block diagonal matrix D and the multipliers used to\n*          obtain the factor U or L as computed by CSPTRF, stored as a\n*          packed triangular matrix.\n*\n*  IPIV    (input) INTEGER array, dimension (N)\n*          Details of the interchanges and the block structure of D\n*          as determined by CSPTRF.\n*\n*  ANORM   (input) REAL\n*          The 1-norm of the original matrix A.\n*\n*  RCOND   (output) REAL\n*          The reciprocal of the condition number of the matrix A,\n*          computed as RCOND = 1/(ANORM * AINVNM), where AINVNM is an\n*          estimate of the 1-norm of inv(A) computed in this routine.\n*\n*  WORK    (workspace) COMPLEX array, dimension (2*N)\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  rcond, info = NumRu::Lapack.cspcon( uplo, ap, ipiv, anorm, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rblapack_uplo = argv[0];
  rblapack_ap = argv[1];
  rblapack_ipiv = argv[2];
  rblapack_anorm = argv[3];
  if (rb_options != Qnil) {
  }

  anorm = (real)NUM2DBL(rblapack_anorm);
  if (!NA_IsNArray(rblapack_ipiv))
    rb_raise(rb_eArgError, "ipiv (3th argument) must be NArray");
  if (NA_RANK(rblapack_ipiv) != 1)
    rb_raise(rb_eArgError, "rank of ipiv (3th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_ipiv);
  if (NA_TYPE(rblapack_ipiv) != NA_LINT)
    rblapack_ipiv = na_change_type(rblapack_ipiv, NA_LINT);
  ipiv = NA_PTR_TYPE(rblapack_ipiv, integer*);
  uplo = StringValueCStr(rblapack_uplo)[0];
  if (!NA_IsNArray(rblapack_ap))
    rb_raise(rb_eArgError, "ap (2th argument) must be NArray");
  if (NA_RANK(rblapack_ap) != 1)
    rb_raise(rb_eArgError, "rank of ap (2th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_ap) != (n*(n+1)/2))
    rb_raise(rb_eRuntimeError, "shape 0 of ap must be %d", n*(n+1)/2);
  if (NA_TYPE(rblapack_ap) != NA_SCOMPLEX)
    rblapack_ap = na_change_type(rblapack_ap, NA_SCOMPLEX);
  ap = NA_PTR_TYPE(rblapack_ap, complex*);
  work = ALLOC_N(complex, (2*n));

  cspcon_(&uplo, &n, ap, ipiv, &anorm, &rcond, work, &info);

  free(work);
  rblapack_rcond = rb_float_new((double)rcond);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(2, rblapack_rcond, rblapack_info);
}

void
init_lapack_cspcon(VALUE mLapack){
  rb_define_module_function(mLapack, "cspcon", rblapack_cspcon, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
