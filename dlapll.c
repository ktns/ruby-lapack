#include "rb_lapack.h"

extern VOID dlapll_(integer *n, doublereal *x, integer *incx, doublereal *y, integer *incy, doublereal *ssmin);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dlapll(int argc, VALUE *argv, VALUE self){
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
  VALUE rblapack_ssmin;
  doublereal ssmin; 
  VALUE rblapack_x_out__;
  doublereal *x_out__;
  VALUE rblapack_y_out__;
  doublereal *y_out__;


  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  ssmin, x, y = NumRu::Lapack.dlapll( n, x, incx, y, incy, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLAPLL( N, X, INCX, Y, INCY, SSMIN )\n\n*  Purpose\n*  =======\n*\n*  Given two column vectors X and Y, let\n*\n*                       A = ( X Y ).\n*\n*  The subroutine first computes the QR factorization of A = Q*R,\n*  and then computes the SVD of the 2-by-2 upper triangular matrix R.\n*  The smaller singular value of R is returned in SSMIN, which is used\n*  as the measurement of the linear dependency of the vectors X and Y.\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The length of the vectors X and Y.\n*\n*  X       (input/output) DOUBLE PRECISION array,\n*                         dimension (1+(N-1)*INCX)\n*          On entry, X contains the N-vector X.\n*          On exit, X is overwritten.\n*\n*  INCX    (input) INTEGER\n*          The increment between successive elements of X. INCX > 0.\n*\n*  Y       (input/output) DOUBLE PRECISION array,\n*                         dimension (1+(N-1)*INCY)\n*          On entry, Y contains the N-vector Y.\n*          On exit, Y is overwritten.\n*\n*  INCY    (input) INTEGER\n*          The increment between successive elements of Y. INCY > 0.\n*\n*  SSMIN   (output) DOUBLE PRECISION\n*          The smallest singular value of the N-by-2 matrix A = ( X Y ).\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  ssmin, x, y = NumRu::Lapack.dlapll( n, x, incx, y, incy, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rblapack_n = argv[0];
  rblapack_x = argv[1];
  rblapack_incx = argv[2];
  rblapack_y = argv[3];
  rblapack_incy = argv[4];
  if (rb_options != Qnil) {
  }

  incy = NUM2INT(rblapack_incy);
  incx = NUM2INT(rblapack_incx);
  n = NUM2INT(rblapack_n);
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

  dlapll_(&n, x, &incx, y, &incy, &ssmin);

  rblapack_ssmin = rb_float_new((double)ssmin);
  return rb_ary_new3(3, rblapack_ssmin, rblapack_x, rblapack_y);
}

void
init_lapack_dlapll(VALUE mLapack){
  rb_define_module_function(mLapack, "dlapll", rblapack_dlapll, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
