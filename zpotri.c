#include "rb_lapack.h"

extern VOID zpotri_(char *uplo, integer *n, doublecomplex *a, integer *lda, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_zpotri(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_uplo;
  char uplo; 
  VALUE rblapack_a;
  doublecomplex *a; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_a_out__;
  doublecomplex *a_out__;

  integer lda;
  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, a = NumRu::Lapack.zpotri( uplo, a, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZPOTRI( UPLO, N, A, LDA, INFO )\n\n*  Purpose\n*  =======\n*\n*  ZPOTRI computes the inverse of a complex Hermitian positive definite\n*  matrix A using the Cholesky factorization A = U**H*U or A = L*L**H\n*  computed by ZPOTRF.\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U':  Upper triangle of A is stored;\n*          = 'L':  Lower triangle of A is stored.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)\n*          On entry, the triangular factor U or L from the Cholesky\n*          factorization A = U**H*U or A = L*L**H, as computed by\n*          ZPOTRF.\n*          On exit, the upper or lower triangle of the (Hermitian)\n*          inverse of A, overwriting the input factor U or L.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,N).\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*          > 0:  if INFO = i, the (i,i) element of the factor U or L is\n*                zero, and the inverse could not be computed.\n*\n\n*  =====================================================================\n*\n*     .. External Functions ..\n      LOGICAL            LSAME\n      EXTERNAL           LSAME\n*     ..\n*     .. External Subroutines ..\n      EXTERNAL           XERBLA, ZLAUUM, ZTRTRI\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          MAX\n*     ..\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, a = NumRu::Lapack.zpotri( uplo, a, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rblapack_uplo = argv[0];
  rblapack_a = argv[1];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_a);
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_DCOMPLEX)
    rblapack_a = na_change_type(rblapack_a, NA_DCOMPLEX);
  a = NA_PTR_TYPE(rblapack_a, doublecomplex*);
  uplo = StringValueCStr(rblapack_uplo)[0];
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rblapack_a_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rblapack_a_out__, doublecomplex*);
  MEMCPY(a_out__, a, doublecomplex, NA_TOTAL(rblapack_a));
  rblapack_a = rblapack_a_out__;
  a = a_out__;

  zpotri_(&uplo, &n, a, &lda, &info);

  rblapack_info = INT2NUM(info);
  return rb_ary_new3(2, rblapack_info, rblapack_a);
}

void
init_lapack_zpotri(VALUE mLapack){
  rb_define_module_function(mLapack, "zpotri", rblapack_zpotri, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
