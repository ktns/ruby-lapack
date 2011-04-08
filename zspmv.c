#include "rb_lapack.h"

extern VOID zspmv_(char *uplo, integer *n, doublecomplex *alpha, doublecomplex *ap, doublecomplex *x, integer *incx, doublecomplex *beta, doublecomplex *y, integer *incy);

static VALUE sHelp, sUsage;

static VALUE
rblapack_zspmv(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_uplo;
  char uplo; 
  VALUE rblapack_n;
  integer n; 
  VALUE rblapack_alpha;
  doublecomplex alpha; 
  VALUE rblapack_ap;
  doublecomplex *ap; 
  VALUE rblapack_x;
  doublecomplex *x; 
  VALUE rblapack_incx;
  integer incx; 
  VALUE rblapack_beta;
  doublecomplex beta; 
  VALUE rblapack_y;
  doublecomplex *y; 
  VALUE rblapack_incy;
  integer incy; 
  VALUE rblapack_y_out__;
  doublecomplex *y_out__;


  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  y = NumRu::Lapack.zspmv( uplo, n, alpha, ap, x, incx, beta, y, incy, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZSPMV( UPLO, N, ALPHA, AP, X, INCX, BETA, Y, INCY )\n\n*  Purpose\n*  =======\n*\n*  ZSPMV  performs the matrix-vector operation\n*\n*     y := alpha*A*x + beta*y,\n*\n*  where alpha and beta are scalars, x and y are n element vectors and\n*  A is an n by n symmetric matrix, supplied in packed form.\n*\n\n*  Arguments\n*  ==========\n*\n*  UPLO     (input) CHARACTER*1\n*           On entry, UPLO specifies whether the upper or lower\n*           triangular part of the matrix A is supplied in the packed\n*           array AP as follows:\n*\n*              UPLO = 'U' or 'u'   The upper triangular part of A is\n*                                  supplied in AP.\n*\n*              UPLO = 'L' or 'l'   The lower triangular part of A is\n*                                  supplied in AP.\n*\n*           Unchanged on exit.\n*\n*  N        (input) INTEGER\n*           On entry, N specifies the order of the matrix A.\n*           N must be at least zero.\n*           Unchanged on exit.\n*\n*  ALPHA    (input) COMPLEX*16\n*           On entry, ALPHA specifies the scalar alpha.\n*           Unchanged on exit.\n*\n*  AP       (input) COMPLEX*16 array, dimension at least\n*           ( ( N*( N + 1 ) )/2 ).\n*           Before entry, with UPLO = 'U' or 'u', the array AP must\n*           contain the upper triangular part of the symmetric matrix\n*           packed sequentially, column by column, so that AP( 1 )\n*           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 )\n*           and a( 2, 2 ) respectively, and so on.\n*           Before entry, with UPLO = 'L' or 'l', the array AP must\n*           contain the lower triangular part of the symmetric matrix\n*           packed sequentially, column by column, so that AP( 1 )\n*           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 )\n*           and a( 3, 1 ) respectively, and so on.\n*           Unchanged on exit.\n*\n*  X        (input) COMPLEX*16 array, dimension at least\n*           ( 1 + ( N - 1 )*abs( INCX ) ).\n*           Before entry, the incremented array X must contain the N-\n*           element vector x.\n*           Unchanged on exit.\n*\n*  INCX     (input) INTEGER\n*           On entry, INCX specifies the increment for the elements of\n*           X. INCX must not be zero.\n*           Unchanged on exit.\n*\n*  BETA     (input) COMPLEX*16\n*           On entry, BETA specifies the scalar beta. When BETA is\n*           supplied as zero then Y need not be set on input.\n*           Unchanged on exit.\n*\n*  Y        (input/output) COMPLEX*16 array, dimension at least\n*           ( 1 + ( N - 1 )*abs( INCY ) ).\n*           Before entry, the incremented array Y must contain the n\n*           element vector y. On exit, Y is overwritten by the updated\n*           vector y.\n*\n*  INCY     (input) INTEGER\n*           On entry, INCY specifies the increment for the elements of\n*           Y. INCY must not be zero.\n*           Unchanged on exit.\n*\n\n* =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  y = NumRu::Lapack.zspmv( uplo, n, alpha, ap, x, incx, beta, y, incy, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 9)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 9)", argc);
  rblapack_uplo = argv[0];
  rblapack_n = argv[1];
  rblapack_alpha = argv[2];
  rblapack_ap = argv[3];
  rblapack_x = argv[4];
  rblapack_incx = argv[5];
  rblapack_beta = argv[6];
  rblapack_y = argv[7];
  rblapack_incy = argv[8];
  if (rb_options != Qnil) {
  }

  uplo = StringValueCStr(rblapack_uplo)[0];
  alpha.r = NUM2DBL(rb_funcall(rblapack_alpha, rb_intern("real"), 0));
  alpha.i = NUM2DBL(rb_funcall(rblapack_alpha, rb_intern("imag"), 0));
  n = NUM2INT(rblapack_n);
  incy = NUM2INT(rblapack_incy);
  incx = NUM2INT(rblapack_incx);
  beta.r = NUM2DBL(rb_funcall(rblapack_beta, rb_intern("real"), 0));
  beta.i = NUM2DBL(rb_funcall(rblapack_beta, rb_intern("imag"), 0));
  if (!NA_IsNArray(rblapack_ap))
    rb_raise(rb_eArgError, "ap (4th argument) must be NArray");
  if (NA_RANK(rblapack_ap) != 1)
    rb_raise(rb_eArgError, "rank of ap (4th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_ap) != (( n*( n + 1 ) )/2))
    rb_raise(rb_eRuntimeError, "shape 0 of ap must be %d", ( n*( n + 1 ) )/2);
  if (NA_TYPE(rblapack_ap) != NA_DCOMPLEX)
    rblapack_ap = na_change_type(rblapack_ap, NA_DCOMPLEX);
  ap = NA_PTR_TYPE(rblapack_ap, doublecomplex*);
  if (!NA_IsNArray(rblapack_y))
    rb_raise(rb_eArgError, "y (8th argument) must be NArray");
  if (NA_RANK(rblapack_y) != 1)
    rb_raise(rb_eArgError, "rank of y (8th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_y) != (1 + ( n - 1 )*abs( incy )))
    rb_raise(rb_eRuntimeError, "shape 0 of y must be %d", 1 + ( n - 1 )*abs( incy ));
  if (NA_TYPE(rblapack_y) != NA_DCOMPLEX)
    rblapack_y = na_change_type(rblapack_y, NA_DCOMPLEX);
  y = NA_PTR_TYPE(rblapack_y, doublecomplex*);
  if (!NA_IsNArray(rblapack_x))
    rb_raise(rb_eArgError, "x (5th argument) must be NArray");
  if (NA_RANK(rblapack_x) != 1)
    rb_raise(rb_eArgError, "rank of x (5th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_x) != (1 + ( n - 1 )*abs( incx )))
    rb_raise(rb_eRuntimeError, "shape 0 of x must be %d", 1 + ( n - 1 )*abs( incx ));
  if (NA_TYPE(rblapack_x) != NA_DCOMPLEX)
    rblapack_x = na_change_type(rblapack_x, NA_DCOMPLEX);
  x = NA_PTR_TYPE(rblapack_x, doublecomplex*);
  {
    int shape[1];
    shape[0] = 1 + ( n - 1 )*abs( incy );
    rblapack_y_out__ = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  y_out__ = NA_PTR_TYPE(rblapack_y_out__, doublecomplex*);
  MEMCPY(y_out__, y, doublecomplex, NA_TOTAL(rblapack_y));
  rblapack_y = rblapack_y_out__;
  y = y_out__;

  zspmv_(&uplo, &n, &alpha, ap, x, &incx, &beta, y, &incy);

  return rblapack_y;
}

void
init_lapack_zspmv(VALUE mLapack){
  rb_define_module_function(mLapack, "zspmv", rblapack_zspmv, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
