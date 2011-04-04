#include "rb_lapack.h"

extern VOID sla_syamv_(integer *uplo, integer *n, real *alpha, real *a, integer *lda, real *x, integer *incx, real *beta, real *y, integer *incy);

static VALUE
rb_sla_syamv(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  integer uplo; 
  VALUE rb_alpha;
  real alpha; 
  VALUE rb_a;
  real *a; 
  VALUE rb_x;
  real *x; 
  VALUE rb_incx;
  integer incx; 
  VALUE rb_beta;
  real beta; 
  VALUE rb_y;
  real *y; 
  VALUE rb_incy;
  integer incy; 
  VALUE rb_y_out__;
  real *y_out__;

  integer lda;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  y = NumRu::Lapack.sla_syamv( uplo, alpha, a, x, incx, beta, y, incy)\n    or\n  NumRu::Lapack.sla_syamv  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLA_SYAMV( UPLO, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY )\n\n*  Purpose\n*  =======\n*\n*  SLA_SYAMV  performs the matrix-vector operation\n*\n*          y := alpha*abs(A)*abs(x) + beta*abs(y),\n*\n*  where alpha and beta are scalars, x and y are vectors and A is an\n*  n by n symmetric matrix.\n*\n*  This function is primarily used in calculating error bounds.\n*  To protect against underflow during evaluation, components in\n*  the resulting vector are perturbed away from zero by (N+1)\n*  times the underflow threshold.  To prevent unnecessarily large\n*  errors for block-structure embedded in general matrices,\n*  \"symbolically\" zero components are not perturbed.  A zero\n*  entry is considered \"symbolic\" if all multiplications involved\n*  in computing that entry have at least one zero multiplicand.\n*\n\n*  Arguments\n*  ==========\n*\n*  UPLO    (input) INTEGER\n*           On entry, UPLO specifies whether the upper or lower\n*           triangular part of the array A is to be referenced as\n*           follows:\n*\n*              UPLO = BLAS_UPPER   Only the upper triangular part of A\n*                                  is to be referenced.\n*\n*              UPLO = BLAS_LOWER   Only the lower triangular part of A\n*                                  is to be referenced.\n*\n*           Unchanged on exit.\n*\n*  N       (input) INTEGER\n*           On entry, N specifies the number of columns of the matrix A.\n*           N must be at least zero.\n*           Unchanged on exit.\n*\n*  ALPHA   (input) REAL            .\n*           On entry, ALPHA specifies the scalar alpha.\n*           Unchanged on exit.\n*\n*  A      - REAL             array of DIMENSION ( LDA, n ).\n*           Before entry, the leading m by n part of the array A must\n*           contain the matrix of coefficients.\n*           Unchanged on exit.\n*\n*  LDA     (input) INTEGER\n*           On entry, LDA specifies the first dimension of A as declared\n*           in the calling (sub) program. LDA must be at least\n*           max( 1, n ).\n*           Unchanged on exit.\n*\n*  X       (input) REAL array, dimension\n*           ( 1 + ( n - 1 )*abs( INCX ) )\n*           Before entry, the incremented array X must contain the\n*           vector x.\n*           Unchanged on exit.\n*\n*  INCX    (input) INTEGER\n*           On entry, INCX specifies the increment for the elements of\n*           X. INCX must not be zero.\n*           Unchanged on exit.\n*\n*  BETA    (input) REAL            .\n*           On entry, BETA specifies the scalar beta. When BETA is\n*           supplied as zero then Y need not be set on input.\n*           Unchanged on exit.\n*\n*  Y       (input/output) REAL array, dimension\n*           ( 1 + ( n - 1 )*abs( INCY ) )\n*           Before entry with BETA non-zero, the incremented array Y\n*           must contain the vector y. On exit, Y is overwritten by the\n*           updated vector y.\n*\n*  INCY    (input) INTEGER\n*           On entry, INCY specifies the increment for the elements of\n*           Y. INCY must not be zero.\n*           Unchanged on exit.\n*\n\n*  Further Details\n*  ===============\n*\n*  Level 2 Blas routine.\n*\n*  -- Written on 22-October-1986.\n*     Jack Dongarra, Argonne National Lab.\n*     Jeremy Du Croz, Nag Central Office.\n*     Sven Hammarling, Nag Central Office.\n*     Richard Hanson, Sandia National Labs.\n*  -- Modified for the absolute-value product, April 2006\n*     Jason Riedy, UC Berkeley\n*\n*  =====================================================================\n*\n\n");
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
  if (lda != (MAX(1, n)))
    rb_raise(rb_eRuntimeError, "shape 0 of a must be %d", MAX(1, n));
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);
  uplo = NUM2INT(rb_uplo);
  alpha = (real)NUM2DBL(rb_alpha);
  beta = (real)NUM2DBL(rb_beta);
  incy = NUM2INT(rb_incy);
  incx = NUM2INT(rb_incx);
  lda = MAX(1, n);
  if (!NA_IsNArray(rb_y))
    rb_raise(rb_eArgError, "y (7th argument) must be NArray");
  if (NA_RANK(rb_y) != 1)
    rb_raise(rb_eArgError, "rank of y (7th argument) must be %d", 1);
  if (NA_SHAPE0(rb_y) != (1 + ( n - 1 )*abs( incy )))
    rb_raise(rb_eRuntimeError, "shape 0 of y must be %d", 1 + ( n - 1 )*abs( incy ));
  if (NA_TYPE(rb_y) != NA_SFLOAT)
    rb_y = na_change_type(rb_y, NA_SFLOAT);
  y = NA_PTR_TYPE(rb_y, real*);
  if (!NA_IsNArray(rb_x))
    rb_raise(rb_eArgError, "x (4th argument) must be NArray");
  if (NA_RANK(rb_x) != 1)
    rb_raise(rb_eArgError, "rank of x (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_x) != (1 + ( n - 1 )*abs( incx )))
    rb_raise(rb_eRuntimeError, "shape 0 of x must be %d", 1 + ( n - 1 )*abs( incx ));
  if (NA_TYPE(rb_x) != NA_SFLOAT)
    rb_x = na_change_type(rb_x, NA_SFLOAT);
  x = NA_PTR_TYPE(rb_x, real*);
  {
    int shape[1];
    shape[0] = 1 + ( n - 1 )*abs( incy );
    rb_y_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  y_out__ = NA_PTR_TYPE(rb_y_out__, real*);
  MEMCPY(y_out__, y, real, NA_TOTAL(rb_y));
  rb_y = rb_y_out__;
  y = y_out__;

  sla_syamv_(&uplo, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);

  return rb_y;
}

void
init_lapack_sla_syamv(VALUE mLapack){
  rb_define_module_function(mLapack, "sla_syamv", rb_sla_syamv, -1);
}
