#include "rb_lapack.h"

extern VOID zsyr_(char *uplo, integer *n, doublecomplex *alpha, doublecomplex *x, integer *incx, doublecomplex *a, integer *lda);

static VALUE sHelp, sUsage;

static VALUE
rblapack_zsyr(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_uplo;
  char uplo; 
  VALUE rblapack_alpha;
  doublecomplex alpha; 
  VALUE rblapack_x;
  doublecomplex *x; 
  VALUE rblapack_incx;
  integer incx; 
  VALUE rblapack_a;
  doublecomplex *a; 
  VALUE rblapack_a_out__;
  doublecomplex *a_out__;

  integer lda;
  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  a = NumRu::Lapack.zsyr( uplo, alpha, x, incx, a, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZSYR( UPLO, N, ALPHA, X, INCX, A, LDA )\n\n*  Purpose\n*  =======\n*\n*  ZSYR   performs the symmetric rank 1 operation\n*\n*     A := alpha*x*( x' ) + A,\n*\n*  where alpha is a complex scalar, x is an n element vector and A is an\n*  n by n symmetric matrix.\n*\n\n*  Arguments\n*  ==========\n*\n*  UPLO     (input) CHARACTER*1\n*           On entry, UPLO specifies whether the upper or lower\n*           triangular part of the array A is to be referenced as\n*           follows:\n*\n*              UPLO = 'U' or 'u'   Only the upper triangular part of A\n*                                  is to be referenced.\n*\n*              UPLO = 'L' or 'l'   Only the lower triangular part of A\n*                                  is to be referenced.\n*\n*           Unchanged on exit.\n*\n*  N        (input) INTEGER\n*           On entry, N specifies the order of the matrix A.\n*           N must be at least zero.\n*           Unchanged on exit.\n*\n*  ALPHA    (input) COMPLEX*16\n*           On entry, ALPHA specifies the scalar alpha.\n*           Unchanged on exit.\n*\n*  X        (input) COMPLEX*16 array, dimension at least\n*           ( 1 + ( N - 1 )*abs( INCX ) ).\n*           Before entry, the incremented array X must contain the N-\n*           element vector x.\n*           Unchanged on exit.\n*\n*  INCX     (input) INTEGER\n*           On entry, INCX specifies the increment for the elements of\n*           X. INCX must not be zero.\n*           Unchanged on exit.\n*\n*  A        (input/output) COMPLEX*16 array, dimension ( LDA, N )\n*           Before entry, with  UPLO = 'U' or 'u', the leading n by n\n*           upper triangular part of the array A must contain the upper\n*           triangular part of the symmetric matrix and the strictly\n*           lower triangular part of A is not referenced. On exit, the\n*           upper triangular part of the array A is overwritten by the\n*           upper triangular part of the updated matrix.\n*           Before entry, with UPLO = 'L' or 'l', the leading n by n\n*           lower triangular part of the array A must contain the lower\n*           triangular part of the symmetric matrix and the strictly\n*           upper triangular part of A is not referenced. On exit, the\n*           lower triangular part of the array A is overwritten by the\n*           lower triangular part of the updated matrix.\n*\n*  LDA      (input) INTEGER\n*           On entry, LDA specifies the first dimension of A as declared\n*           in the calling (sub) program. LDA must be at least\n*           max( 1, N ).\n*           Unchanged on exit.\n*\n\n* =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  a = NumRu::Lapack.zsyr( uplo, alpha, x, incx, a, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rblapack_uplo = argv[0];
  rblapack_alpha = argv[1];
  rblapack_x = argv[2];
  rblapack_incx = argv[3];
  rblapack_a = argv[4];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (5th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (5th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_a);
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_DCOMPLEX)
    rblapack_a = na_change_type(rblapack_a, NA_DCOMPLEX);
  a = NA_PTR_TYPE(rblapack_a, doublecomplex*);
  uplo = StringValueCStr(rblapack_uplo)[0];
  alpha.r = NUM2DBL(rb_funcall(rblapack_alpha, rb_intern("real"), 0));
  alpha.i = NUM2DBL(rb_funcall(rblapack_alpha, rb_intern("imag"), 0));
  incx = NUM2INT(rblapack_incx);
  if (!NA_IsNArray(rblapack_x))
    rb_raise(rb_eArgError, "x (3th argument) must be NArray");
  if (NA_RANK(rblapack_x) != 1)
    rb_raise(rb_eArgError, "rank of x (3th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_x) != (1 + ( n - 1 )*abs( incx )))
    rb_raise(rb_eRuntimeError, "shape 0 of x must be %d", 1 + ( n - 1 )*abs( incx ));
  if (NA_TYPE(rblapack_x) != NA_DCOMPLEX)
    rblapack_x = na_change_type(rblapack_x, NA_DCOMPLEX);
  x = NA_PTR_TYPE(rblapack_x, doublecomplex*);
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

  zsyr_(&uplo, &n, &alpha, x, &incx, a, &lda);

  return rblapack_a;
}

void
init_lapack_zsyr(VALUE mLapack){
  rb_define_module_function(mLapack, "zsyr", rblapack_zsyr, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
