#include "rb_lapack.h"

extern VOID zsyr_(char *uplo, integer *n, doublecomplex *alpha, doublecomplex *x, integer *incx, doublecomplex *a, integer *lda);

static VALUE
rb_zsyr(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_alpha;
  doublecomplex alpha; 
  VALUE rb_x;
  doublecomplex *x; 
  VALUE rb_incx;
  integer incx; 
  VALUE rb_a;
  doublecomplex *a; 
  VALUE rb_a_out__;
  doublecomplex *a_out__;

  integer lda;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  a = NumRu::Lapack.zsyr( uplo, alpha, x, incx, a)\n    or\n  NumRu::Lapack.zsyr  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZSYR( UPLO, N, ALPHA, X, INCX, A, LDA )\n\n*  Purpose\n*  =======\n*\n*  ZSYR   performs the symmetric rank 1 operation\n*\n*     A := alpha*x*( x' ) + A,\n*\n*  where alpha is a complex scalar, x is an n element vector and A is an\n*  n by n symmetric matrix.\n*\n\n*  Arguments\n*  ==========\n*\n*  UPLO     (input) CHARACTER*1\n*           On entry, UPLO specifies whether the upper or lower\n*           triangular part of the array A is to be referenced as\n*           follows:\n*\n*              UPLO = 'U' or 'u'   Only the upper triangular part of A\n*                                  is to be referenced.\n*\n*              UPLO = 'L' or 'l'   Only the lower triangular part of A\n*                                  is to be referenced.\n*\n*           Unchanged on exit.\n*\n*  N        (input) INTEGER\n*           On entry, N specifies the order of the matrix A.\n*           N must be at least zero.\n*           Unchanged on exit.\n*\n*  ALPHA    (input) COMPLEX*16\n*           On entry, ALPHA specifies the scalar alpha.\n*           Unchanged on exit.\n*\n*  X        (input) COMPLEX*16 array, dimension at least\n*           ( 1 + ( N - 1 )*abs( INCX ) ).\n*           Before entry, the incremented array X must contain the N-\n*           element vector x.\n*           Unchanged on exit.\n*\n*  INCX     (input) INTEGER\n*           On entry, INCX specifies the increment for the elements of\n*           X. INCX must not be zero.\n*           Unchanged on exit.\n*\n*  A        (input/output) COMPLEX*16 array, dimension ( LDA, N )\n*           Before entry, with  UPLO = 'U' or 'u', the leading n by n\n*           upper triangular part of the array A must contain the upper\n*           triangular part of the symmetric matrix and the strictly\n*           lower triangular part of A is not referenced. On exit, the\n*           upper triangular part of the array A is overwritten by the\n*           upper triangular part of the updated matrix.\n*           Before entry, with UPLO = 'L' or 'l', the leading n by n\n*           lower triangular part of the array A must contain the lower\n*           triangular part of the symmetric matrix and the strictly\n*           upper triangular part of A is not referenced. On exit, the\n*           lower triangular part of the array A is overwritten by the\n*           lower triangular part of the updated matrix.\n*\n*  LDA      (input) INTEGER\n*           On entry, LDA specifies the first dimension of A as declared\n*           in the calling (sub) program. LDA must be at least\n*           max( 1, N ).\n*           Unchanged on exit.\n*\n\n* =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_uplo = argv[0];
  rb_alpha = argv[1];
  rb_x = argv[2];
  rb_incx = argv[3];
  rb_a = argv[4];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (5th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (5th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_DCOMPLEX)
    rb_a = na_change_type(rb_a, NA_DCOMPLEX);
  a = NA_PTR_TYPE(rb_a, doublecomplex*);
  uplo = StringValueCStr(rb_uplo)[0];
  alpha.r = NUM2DBL(rb_funcall(rb_alpha, rb_intern("real"), 0));
  alpha.i = NUM2DBL(rb_funcall(rb_alpha, rb_intern("imag"), 0));
  incx = NUM2INT(rb_incx);
  if (!NA_IsNArray(rb_x))
    rb_raise(rb_eArgError, "x (3th argument) must be NArray");
  if (NA_RANK(rb_x) != 1)
    rb_raise(rb_eArgError, "rank of x (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_x) != (1 + ( n - 1 )*abs( incx )))
    rb_raise(rb_eRuntimeError, "shape 0 of x must be %d", 1 + ( n - 1 )*abs( incx ));
  if (NA_TYPE(rb_x) != NA_DCOMPLEX)
    rb_x = na_change_type(rb_x, NA_DCOMPLEX);
  x = NA_PTR_TYPE(rb_x, doublecomplex*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rb_a_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rb_a_out__, doublecomplex*);
  MEMCPY(a_out__, a, doublecomplex, NA_TOTAL(rb_a));
  rb_a = rb_a_out__;
  a = a_out__;

  zsyr_(&uplo, &n, &alpha, x, &incx, a, &lda);

  return rb_a;
}

void
init_lapack_zsyr(VALUE mLapack){
  rb_define_module_function(mLapack, "zsyr", rb_zsyr, -1);
}
