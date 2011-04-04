#include "rb_lapack.h"

extern VOID csymv_(char *uplo, integer *n, complex *alpha, complex *a, integer *lda, complex *x, integer *incx, complex *beta, complex *y, integer *incy);

static VALUE
rb_csymv(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_alpha;
  complex alpha; 
  VALUE rb_a;
  complex *a; 
  VALUE rb_x;
  complex *x; 
  VALUE rb_incx;
  integer incx; 
  VALUE rb_beta;
  complex beta; 
  VALUE rb_y;
  complex *y; 
  VALUE rb_incy;
  integer incy; 
  VALUE rb_y_out__;
  complex *y_out__;

  integer lda;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  y = NumRu::Lapack.csymv( uplo, alpha, a, x, incx, beta, y, incy)\n    or\n  NumRu::Lapack.csymv  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE CSYMV( UPLO, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY )\n\n*  Purpose\n*  =======\n*\n*  CSYMV  performs the matrix-vector  operation\n*\n*     y := alpha*A*x + beta*y,\n*\n*  where alpha and beta are scalars, x and y are n element vectors and\n*  A is an n by n symmetric matrix.\n*\n\n*  Arguments\n*  ==========\n*\n*  UPLO     (input) CHARACTER*1\n*           On entry, UPLO specifies whether the upper or lower\n*           triangular part of the array A is to be referenced as\n*           follows:\n*\n*              UPLO = 'U' or 'u'   Only the upper triangular part of A\n*                                  is to be referenced.\n*\n*              UPLO = 'L' or 'l'   Only the lower triangular part of A\n*                                  is to be referenced.\n*\n*           Unchanged on exit.\n*\n*  N        (input) INTEGER\n*           On entry, N specifies the order of the matrix A.\n*           N must be at least zero.\n*           Unchanged on exit.\n*\n*  ALPHA    (input) COMPLEX\n*           On entry, ALPHA specifies the scalar alpha.\n*           Unchanged on exit.\n*\n*  A        (input) COMPLEX array, dimension ( LDA, N )\n*           Before entry, with  UPLO = 'U' or 'u', the leading n by n\n*           upper triangular part of the array A must contain the upper\n*           triangular part of the symmetric matrix and the strictly\n*           lower triangular part of A is not referenced.\n*           Before entry, with UPLO = 'L' or 'l', the leading n by n\n*           lower triangular part of the array A must contain the lower\n*           triangular part of the symmetric matrix and the strictly\n*           upper triangular part of A is not referenced.\n*           Unchanged on exit.\n*\n*  LDA      (input) INTEGER\n*           On entry, LDA specifies the first dimension of A as declared\n*           in the calling (sub) program. LDA must be at least\n*           max( 1, N ).\n*           Unchanged on exit.\n*\n*  X        (input) COMPLEX array, dimension at least\n*           ( 1 + ( N - 1 )*abs( INCX ) ).\n*           Before entry, the incremented array X must contain the N-\n*           element vector x.\n*           Unchanged on exit.\n*\n*  INCX     (input) INTEGER\n*           On entry, INCX specifies the increment for the elements of\n*           X. INCX must not be zero.\n*           Unchanged on exit.\n*\n*  BETA     (input) COMPLEX\n*           On entry, BETA specifies the scalar beta. When BETA is\n*           supplied as zero then Y need not be set on input.\n*           Unchanged on exit.\n*\n*  Y        (input/output) COMPLEX array, dimension at least\n*           ( 1 + ( N - 1 )*abs( INCY ) ).\n*           Before entry, the incremented array Y must contain the n\n*           element vector y. On exit, Y is overwritten by the updated\n*           vector y.\n*\n*  INCY     (input) INTEGER\n*           On entry, INCY specifies the increment for the elements of\n*           Y. INCY must not be zero.\n*           Unchanged on exit.\n*\n\n* =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 8)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 8)", argc);
  rb_uplo = argv[0];
  rb_alpha = argv[1];
  rb_a = argv[2];
  rb_x = argv[3];
  rb_incx = argv[4];
  rb_beta = argv[5];
  rb_y = argv[6];
  rb_incy = argv[7];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SCOMPLEX)
    rb_a = na_change_type(rb_a, NA_SCOMPLEX);
  a = NA_PTR_TYPE(rb_a, complex*);
  uplo = StringValueCStr(rb_uplo)[0];
  alpha.r = (real)NUM2DBL(rb_funcall(rb_alpha, rb_intern("real"), 0));
  alpha.i = (real)NUM2DBL(rb_funcall(rb_alpha, rb_intern("imag"), 0));
  beta.r = (real)NUM2DBL(rb_funcall(rb_beta, rb_intern("real"), 0));
  beta.i = (real)NUM2DBL(rb_funcall(rb_beta, rb_intern("imag"), 0));
  incy = NUM2INT(rb_incy);
  incx = NUM2INT(rb_incx);
  if (!NA_IsNArray(rb_y))
    rb_raise(rb_eArgError, "y (7th argument) must be NArray");
  if (NA_RANK(rb_y) != 1)
    rb_raise(rb_eArgError, "rank of y (7th argument) must be %d", 1);
  if (NA_SHAPE0(rb_y) != (1 + ( n - 1 )*abs( incy )))
    rb_raise(rb_eRuntimeError, "shape 0 of y must be %d", 1 + ( n - 1 )*abs( incy ));
  if (NA_TYPE(rb_y) != NA_SCOMPLEX)
    rb_y = na_change_type(rb_y, NA_SCOMPLEX);
  y = NA_PTR_TYPE(rb_y, complex*);
  if (!NA_IsNArray(rb_x))
    rb_raise(rb_eArgError, "x (4th argument) must be NArray");
  if (NA_RANK(rb_x) != 1)
    rb_raise(rb_eArgError, "rank of x (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_x) != (1 + ( n - 1 )*abs( incx )))
    rb_raise(rb_eRuntimeError, "shape 0 of x must be %d", 1 + ( n - 1 )*abs( incx ));
  if (NA_TYPE(rb_x) != NA_SCOMPLEX)
    rb_x = na_change_type(rb_x, NA_SCOMPLEX);
  x = NA_PTR_TYPE(rb_x, complex*);
  {
    int shape[1];
    shape[0] = 1 + ( n - 1 )*abs( incy );
    rb_y_out__ = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  y_out__ = NA_PTR_TYPE(rb_y_out__, complex*);
  MEMCPY(y_out__, y, complex, NA_TOTAL(rb_y));
  rb_y = rb_y_out__;
  y = y_out__;

  csymv_(&uplo, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);

  return rb_y;
}

void
init_lapack_csymv(VALUE mLapack){
  rb_define_module_function(mLapack, "csymv", rb_csymv, -1);
}
