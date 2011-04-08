#include "rb_lapack.h"

extern VOID dlargv_(integer *n, doublereal *x, integer *incx, doublereal *y, integer *incy, doublereal *c, integer *incc);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dlargv(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_n;
  integer n; 
  VALUE rblapack_x;
  doublereal *x; 
  VALUE rblapack_incx;
  integer incx; 
  VALUE rblapack_y;
  doublereal *y; 
  VALUE rblapack_incy;
  integer incy; 
  VALUE rblapack_incc;
  integer incc; 
  VALUE rblapack_c;
  doublereal *c; 
  VALUE rblapack_x_out__;
  doublereal *x_out__;
  VALUE rblapack_y_out__;
  doublereal *y_out__;


  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  c, x, y = NumRu::Lapack.dlargv( n, x, incx, y, incy, incc, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLARGV( N, X, INCX, Y, INCY, C, INCC )\n\n*  Purpose\n*  =======\n*\n*  DLARGV generates a vector of real plane rotations, determined by\n*  elements of the real vectors x and y. For i = 1,2,...,n\n*\n*     (  c(i)  s(i) ) ( x(i) ) = ( a(i) )\n*     ( -s(i)  c(i) ) ( y(i) ) = (   0  )\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The number of plane rotations to be generated.\n*\n*  X       (input/output) DOUBLE PRECISION array,\n*                         dimension (1+(N-1)*INCX)\n*          On entry, the vector x.\n*          On exit, x(i) is overwritten by a(i), for i = 1,...,n.\n*\n*  INCX    (input) INTEGER\n*          The increment between elements of X. INCX > 0.\n*\n*  Y       (input/output) DOUBLE PRECISION array,\n*                         dimension (1+(N-1)*INCY)\n*          On entry, the vector y.\n*          On exit, the sines of the plane rotations.\n*\n*  INCY    (input) INTEGER\n*          The increment between elements of Y. INCY > 0.\n*\n*  C       (output) DOUBLE PRECISION array, dimension (1+(N-1)*INCC)\n*          The cosines of the plane rotations.\n*\n*  INCC    (input) INTEGER\n*          The increment between elements of C. INCC > 0.\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  c, x, y = NumRu::Lapack.dlargv( n, x, incx, y, incy, incc, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rblapack_n = argv[0];
  rblapack_x = argv[1];
  rblapack_incx = argv[2];
  rblapack_y = argv[3];
  rblapack_incy = argv[4];
  rblapack_incc = argv[5];
  if (rb_options != Qnil) {
  }

  incy = NUM2INT(rblapack_incy);
  incc = NUM2INT(rblapack_incc);
  n = NUM2INT(rblapack_n);
  incx = NUM2INT(rblapack_incx);
  if (!NA_IsNArray(rblapack_y))
    rb_raise(rb_eArgError, "y (4th argument) must be NArray");
  if (NA_RANK(rblapack_y) != 1)
    rb_raise(rb_eArgError, "rank of y (4th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_y) != (1+(n-1)*incy))
    rb_raise(rb_eRuntimeError, "shape 0 of y must be %d", 1+(n-1)*incy);
  if (NA_TYPE(rblapack_y) != NA_DFLOAT)
    rblapack_y = na_change_type(rblapack_y, NA_DFLOAT);
  y = NA_PTR_TYPE(rblapack_y, doublereal*);
  if (!NA_IsNArray(rblapack_x))
    rb_raise(rb_eArgError, "x (2th argument) must be NArray");
  if (NA_RANK(rblapack_x) != 1)
    rb_raise(rb_eArgError, "rank of x (2th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_x) != (1+(n-1)*incx))
    rb_raise(rb_eRuntimeError, "shape 0 of x must be %d", 1+(n-1)*incx);
  if (NA_TYPE(rblapack_x) != NA_DFLOAT)
    rblapack_x = na_change_type(rblapack_x, NA_DFLOAT);
  x = NA_PTR_TYPE(rblapack_x, doublereal*);
  {
    int shape[1];
    shape[0] = 1+(n-1)*incc;
    rblapack_c = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  c = NA_PTR_TYPE(rblapack_c, doublereal*);
  {
    int shape[1];
    shape[0] = 1+(n-1)*incx;
    rblapack_x_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  x_out__ = NA_PTR_TYPE(rblapack_x_out__, doublereal*);
  MEMCPY(x_out__, x, doublereal, NA_TOTAL(rblapack_x));
  rblapack_x = rblapack_x_out__;
  x = x_out__;
  {
    int shape[1];
    shape[0] = 1+(n-1)*incy;
    rblapack_y_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  y_out__ = NA_PTR_TYPE(rblapack_y_out__, doublereal*);
  MEMCPY(y_out__, y, doublereal, NA_TOTAL(rblapack_y));
  rblapack_y = rblapack_y_out__;
  y = y_out__;

  dlargv_(&n, x, &incx, y, &incy, c, &incc);

  return rb_ary_new3(3, rblapack_c, rblapack_x, rblapack_y);
}

void
init_lapack_dlargv(VALUE mLapack){
  rb_define_module_function(mLapack, "dlargv", rblapack_dlargv, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
