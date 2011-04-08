#include "rb_lapack.h"

extern VOID ssytrs2_(char *uplo, integer *n, integer *nrhs, real *a, integer *lda, integer *ipiv, real *b, integer *ldb, real *work, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_ssytrs2(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_uplo;
  char uplo; 
  VALUE rblapack_a;
  real *a; 
  VALUE rblapack_ipiv;
  integer *ipiv; 
  VALUE rblapack_b;
  real *b; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_b_out__;
  real *b_out__;
  real *work;

  integer lda;
  integer n;
  integer ldb;
  integer nrhs;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, b = NumRu::Lapack.ssytrs2( uplo, a, ipiv, b, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE SSYTRS2( UPLO, N, NRHS, A, LDA, IPIV, B, LDB,  WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  SSYTRS2 solves a system of linear equations A*X = B with a real\n*  symmetric matrix A using the factorization A = U*D*U**T or\n*  A = L*D*L**T computed by SSYTRF and converted by SSYCONV.\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO    (input) CHARACTER*1\n*          Specifies whether the details of the factorization are stored\n*          as an upper or lower triangular matrix.\n*          = 'U':  Upper triangular, form is A = U*D*U**T;\n*          = 'L':  Lower triangular, form is A = L*D*L**T.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  NRHS    (input) INTEGER\n*          The number of right hand sides, i.e., the number of columns\n*          of the matrix B.  NRHS >= 0.\n*\n*  A       (input) REAL array, dimension (LDA,N)\n*          The block diagonal matrix D and the multipliers used to\n*          obtain the factor U or L as computed by SSYTRF.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,N).\n*\n*  IPIV    (input) INTEGER array, dimension (N)\n*          Details of the interchanges and the block structure of D\n*          as determined by SSYTRF.\n*\n*  B       (input/output) REAL array, dimension (LDB,NRHS)\n*          On entry, the right hand side matrix B.\n*          On exit, the solution matrix X.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B.  LDB >= max(1,N).\n*\n*  WORK    (workspace) REAL array, dimension (N)\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, b = NumRu::Lapack.ssytrs2( uplo, a, ipiv, b, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rblapack_uplo = argv[0];
  rblapack_a = argv[1];
  rblapack_ipiv = argv[2];
  rblapack_b = argv[3];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_ipiv))
    rb_raise(rb_eArgError, "ipiv (3th argument) must be NArray");
  if (NA_RANK(rblapack_ipiv) != 1)
    rb_raise(rb_eArgError, "rank of ipiv (3th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_ipiv);
  if (NA_TYPE(rblapack_ipiv) != NA_LINT)
    rblapack_ipiv = na_change_type(rblapack_ipiv, NA_LINT);
  ipiv = NA_PTR_TYPE(rblapack_ipiv, integer*);
  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_a) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of a must be the same as shape 0 of ipiv");
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_SFLOAT)
    rblapack_a = na_change_type(rblapack_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rblapack_a, real*);
  if (!NA_IsNArray(rblapack_b))
    rb_raise(rb_eArgError, "b (4th argument) must be NArray");
  if (NA_RANK(rblapack_b) != 2)
    rb_raise(rb_eArgError, "rank of b (4th argument) must be %d", 2);
  nrhs = NA_SHAPE1(rblapack_b);
  ldb = NA_SHAPE0(rblapack_b);
  if (NA_TYPE(rblapack_b) != NA_SFLOAT)
    rblapack_b = na_change_type(rblapack_b, NA_SFLOAT);
  b = NA_PTR_TYPE(rblapack_b, real*);
  uplo = StringValueCStr(rblapack_uplo)[0];
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = nrhs;
    rblapack_b_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rblapack_b_out__, real*);
  MEMCPY(b_out__, b, real, NA_TOTAL(rblapack_b));
  rblapack_b = rblapack_b_out__;
  b = b_out__;
  work = ALLOC_N(real, (n));

  ssytrs2_(&uplo, &n, &nrhs, a, &lda, ipiv, b, &ldb, work, &info);

  free(work);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(2, rblapack_info, rblapack_b);
}

void
init_lapack_ssytrs2(VALUE mLapack){
  rb_define_module_function(mLapack, "ssytrs2", rblapack_ssytrs2, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
