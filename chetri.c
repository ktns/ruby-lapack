#include "rb_lapack.h"

extern VOID chetri_(char *uplo, integer *n, complex *a, integer *lda, integer *ipiv, complex *work, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_chetri(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_uplo;
  char uplo; 
  VALUE rblapack_a;
  complex *a; 
  VALUE rblapack_ipiv;
  integer *ipiv; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_a_out__;
  complex *a_out__;
  complex *work;

  integer lda;
  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, a = NumRu::Lapack.chetri( uplo, a, ipiv, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE CHETRI( UPLO, N, A, LDA, IPIV, WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  CHETRI computes the inverse of a complex Hermitian indefinite matrix\n*  A using the factorization A = U*D*U**H or A = L*D*L**H computed by\n*  CHETRF.\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO    (input) CHARACTER*1\n*          Specifies whether the details of the factorization are stored\n*          as an upper or lower triangular matrix.\n*          = 'U':  Upper triangular, form is A = U*D*U**H;\n*          = 'L':  Lower triangular, form is A = L*D*L**H.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  A       (input/output) COMPLEX array, dimension (LDA,N)\n*          On entry, the block diagonal matrix D and the multipliers\n*          used to obtain the factor U or L as computed by CHETRF.\n*\n*          On exit, if INFO = 0, the (Hermitian) inverse of the original\n*          matrix.  If UPLO = 'U', the upper triangular part of the\n*          inverse is formed and the part of A below the diagonal is not\n*          referenced; if UPLO = 'L' the lower triangular part of the\n*          inverse is formed and the part of A above the diagonal is\n*          not referenced.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,N).\n*\n*  IPIV    (input) INTEGER array, dimension (N)\n*          Details of the interchanges and the block structure of D\n*          as determined by CHETRF.\n*\n*  WORK    (workspace) COMPLEX array, dimension (N)\n*\n*  INFO    (output) INTEGER\n*          = 0: successful exit\n*          < 0: if INFO = -i, the i-th argument had an illegal value\n*          > 0: if INFO = i, D(i,i) = 0; the matrix is singular and its\n*               inverse could not be computed.\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, a = NumRu::Lapack.chetri( uplo, a, ipiv, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rblapack_uplo = argv[0];
  rblapack_a = argv[1];
  rblapack_ipiv = argv[2];
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
  if (NA_TYPE(rblapack_a) != NA_SCOMPLEX)
    rblapack_a = na_change_type(rblapack_a, NA_SCOMPLEX);
  a = NA_PTR_TYPE(rblapack_a, complex*);
  uplo = StringValueCStr(rblapack_uplo)[0];
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rblapack_a_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rblapack_a_out__, complex*);
  MEMCPY(a_out__, a, complex, NA_TOTAL(rblapack_a));
  rblapack_a = rblapack_a_out__;
  a = a_out__;
  work = ALLOC_N(complex, (n));

  chetri_(&uplo, &n, a, &lda, ipiv, work, &info);

  free(work);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(2, rblapack_info, rblapack_a);
}

void
init_lapack_chetri(VALUE mLapack){
  rb_define_module_function(mLapack, "chetri", rblapack_chetri, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
