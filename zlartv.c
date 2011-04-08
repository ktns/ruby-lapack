#include "rb_lapack.h"

extern VOID zlartv_(integer *n, doublecomplex *x, integer *incx, doublecomplex *y, integer *incy, doublereal *c, doublecomplex *s, integer *incc);

static VALUE sHelp, sUsage;

static VALUE
rblapack_zlartv(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_n;
  integer n; 
  VALUE rblapack_x;
  doublecomplex *x; 
  VALUE rblapack_incx;
  integer incx; 
  VALUE rblapack_y;
  doublecomplex *y; 
  VALUE rblapack_incy;
  integer incy; 
  VALUE rblapack_c;
  doublereal *c; 
  VALUE rblapack_s;
  doublecomplex *s; 
  VALUE rblapack_incc;
  integer incc; 
  VALUE rblapack_x_out__;
  doublecomplex *x_out__;
  VALUE rblapack_y_out__;
  doublecomplex *y_out__;


  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  x, y = NumRu::Lapack.zlartv( n, x, incx, y, incy, c, s, incc, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZLARTV( N, X, INCX, Y, INCY, C, S, INCC )\n\n*  Purpose\n*  =======\n*\n*  ZLARTV applies a vector of complex plane rotations with real cosines\n*  to elements of the complex vectors x and y. For i = 1,2,...,n\n*\n*     ( x(i) ) := (        c(i)   s(i) ) ( x(i) )\n*     ( y(i) )    ( -conjg(s(i))  c(i) ) ( y(i) )\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The number of plane rotations to be applied.\n*\n*  X       (input/output) COMPLEX*16 array, dimension (1+(N-1)*INCX)\n*          The vector x.\n*\n*  INCX    (input) INTEGER\n*          The increment between elements of X. INCX > 0.\n*\n*  Y       (input/output) COMPLEX*16 array, dimension (1+(N-1)*INCY)\n*          The vector y.\n*\n*  INCY    (input) INTEGER\n*          The increment between elements of Y. INCY > 0.\n*\n*  C       (input) DOUBLE PRECISION array, dimension (1+(N-1)*INCC)\n*          The cosines of the plane rotations.\n*\n*  S       (input) COMPLEX*16 array, dimension (1+(N-1)*INCC)\n*          The sines of the plane rotations.\n*\n*  INCC    (input) INTEGER\n*          The increment between elements of C and S. INCC > 0.\n*\n\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      INTEGER            I, IC, IX, IY\n      COMPLEX*16         XI, YI\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          DCONJG\n*     ..\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  x, y = NumRu::Lapack.zlartv( n, x, incx, y, incy, c, s, incc, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 8)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 8)", argc);
  rblapack_n = argv[0];
  rblapack_x = argv[1];
  rblapack_incx = argv[2];
  rblapack_y = argv[3];
  rblapack_incy = argv[4];
  rblapack_c = argv[5];
  rblapack_s = argv[6];
  rblapack_incc = argv[7];
  if (rb_options != Qnil) {
  }

  incx = NUM2INT(rblapack_incx);
  incy = NUM2INT(rblapack_incy);
  incc = NUM2INT(rblapack_incc);
  n = NUM2INT(rblapack_n);
  if (!NA_IsNArray(rblapack_x))
    rb_raise(rb_eArgError, "x (2th argument) must be NArray");
  if (NA_RANK(rblapack_x) != 1)
    rb_raise(rb_eArgError, "rank of x (2th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_x) != (1+(n-1)*incx))
    rb_raise(rb_eRuntimeError, "shape 0 of x must be %d", 1+(n-1)*incx);
  if (NA_TYPE(rblapack_x) != NA_DCOMPLEX)
    rblapack_x = na_change_type(rblapack_x, NA_DCOMPLEX);
  x = NA_PTR_TYPE(rblapack_x, doublecomplex*);
  if (!NA_IsNArray(rblapack_y))
    rb_raise(rb_eArgError, "y (4th argument) must be NArray");
  if (NA_RANK(rblapack_y) != 1)
    rb_raise(rb_eArgError, "rank of y (4th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_y) != (1+(n-1)*incy))
    rb_raise(rb_eRuntimeError, "shape 0 of y must be %d", 1+(n-1)*incy);
  if (NA_TYPE(rblapack_y) != NA_DCOMPLEX)
    rblapack_y = na_change_type(rblapack_y, NA_DCOMPLEX);
  y = NA_PTR_TYPE(rblapack_y, doublecomplex*);
  if (!NA_IsNArray(rblapack_c))
    rb_raise(rb_eArgError, "c (6th argument) must be NArray");
  if (NA_RANK(rblapack_c) != 1)
    rb_raise(rb_eArgError, "rank of c (6th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_c) != (1+(n-1)*incc))
    rb_raise(rb_eRuntimeError, "shape 0 of c must be %d", 1+(n-1)*incc);
  if (NA_TYPE(rblapack_c) != NA_DFLOAT)
    rblapack_c = na_change_type(rblapack_c, NA_DFLOAT);
  c = NA_PTR_TYPE(rblapack_c, doublereal*);
  if (!NA_IsNArray(rblapack_s))
    rb_raise(rb_eArgError, "s (7th argument) must be NArray");
  if (NA_RANK(rblapack_s) != 1)
    rb_raise(rb_eArgError, "rank of s (7th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_s) != (1+(n-1)*incc))
    rb_raise(rb_eRuntimeError, "shape 0 of s must be %d", 1+(n-1)*incc);
  if (NA_TYPE(rblapack_s) != NA_DCOMPLEX)
    rblapack_s = na_change_type(rblapack_s, NA_DCOMPLEX);
  s = NA_PTR_TYPE(rblapack_s, doublecomplex*);
  {
    int shape[1];
    shape[0] = 1+(n-1)*incx;
    rblapack_x_out__ = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  x_out__ = NA_PTR_TYPE(rblapack_x_out__, doublecomplex*);
  MEMCPY(x_out__, x, doublecomplex, NA_TOTAL(rblapack_x));
  rblapack_x = rblapack_x_out__;
  x = x_out__;
  {
    int shape[1];
    shape[0] = 1+(n-1)*incy;
    rblapack_y_out__ = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  y_out__ = NA_PTR_TYPE(rblapack_y_out__, doublecomplex*);
  MEMCPY(y_out__, y, doublecomplex, NA_TOTAL(rblapack_y));
  rblapack_y = rblapack_y_out__;
  y = y_out__;

  zlartv_(&n, x, &incx, y, &incy, c, s, &incc);

  return rb_ary_new3(2, rblapack_x, rblapack_y);
}

void
init_lapack_zlartv(VALUE mLapack){
  rb_define_module_function(mLapack, "zlartv", rblapack_zlartv, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
