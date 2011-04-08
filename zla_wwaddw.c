#include "rb_lapack.h"

extern VOID zla_wwaddw_(integer *n, doublecomplex *x, doublecomplex *y, doublecomplex *w);

static VALUE sHelp, sUsage;

static VALUE
rblapack_zla_wwaddw(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_x;
  doublecomplex *x; 
  VALUE rblapack_y;
  doublecomplex *y; 
  VALUE rblapack_w;
  doublecomplex *w; 
  VALUE rblapack_x_out__;
  doublecomplex *x_out__;
  VALUE rblapack_y_out__;
  doublecomplex *y_out__;

  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  x, y = NumRu::Lapack.zla_wwaddw( x, y, w, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZLA_WWADDW( N, X, Y, W )\n\n*     Purpose\n*     =======\n*\n*     ZLA_WWADDW adds a vector W into a doubled-single vector (X, Y).\n*\n*     This works for all extant IBM's hex and binary floating point\n*     arithmetics, but not for decimal.\n*\n\n*     Arguments\n*     =========\n*\n*     N      (input) INTEGER\n*            The length of vectors X, Y, and W.\n*\n*     X      (input/output) COMPLEX*16 array, dimension (N)\n*            The first part of the doubled-single accumulation vector.\n*\n*     Y      (input/output) COMPLEX*16 array, dimension (N)\n*            The second part of the doubled-single accumulation vector.\n*\n*     W      (input) COMPLEX*16 array, dimension (N)\n*            The vector to be added.\n*\n\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      COMPLEX*16         S\n      INTEGER            I\n*     ..\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  x, y = NumRu::Lapack.zla_wwaddw( x, y, w, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rblapack_x = argv[0];
  rblapack_y = argv[1];
  rblapack_w = argv[2];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_w))
    rb_raise(rb_eArgError, "w (3th argument) must be NArray");
  if (NA_RANK(rblapack_w) != 1)
    rb_raise(rb_eArgError, "rank of w (3th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_w);
  if (NA_TYPE(rblapack_w) != NA_DCOMPLEX)
    rblapack_w = na_change_type(rblapack_w, NA_DCOMPLEX);
  w = NA_PTR_TYPE(rblapack_w, doublecomplex*);
  if (!NA_IsNArray(rblapack_x))
    rb_raise(rb_eArgError, "x (1th argument) must be NArray");
  if (NA_RANK(rblapack_x) != 1)
    rb_raise(rb_eArgError, "rank of x (1th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_x) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of x must be the same as shape 0 of w");
  if (NA_TYPE(rblapack_x) != NA_DCOMPLEX)
    rblapack_x = na_change_type(rblapack_x, NA_DCOMPLEX);
  x = NA_PTR_TYPE(rblapack_x, doublecomplex*);
  if (!NA_IsNArray(rblapack_y))
    rb_raise(rb_eArgError, "y (2th argument) must be NArray");
  if (NA_RANK(rblapack_y) != 1)
    rb_raise(rb_eArgError, "rank of y (2th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_y) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of y must be the same as shape 0 of w");
  if (NA_TYPE(rblapack_y) != NA_DCOMPLEX)
    rblapack_y = na_change_type(rblapack_y, NA_DCOMPLEX);
  y = NA_PTR_TYPE(rblapack_y, doublecomplex*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_x_out__ = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  x_out__ = NA_PTR_TYPE(rblapack_x_out__, doublecomplex*);
  MEMCPY(x_out__, x, doublecomplex, NA_TOTAL(rblapack_x));
  rblapack_x = rblapack_x_out__;
  x = x_out__;
  {
    int shape[1];
    shape[0] = n;
    rblapack_y_out__ = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  y_out__ = NA_PTR_TYPE(rblapack_y_out__, doublecomplex*);
  MEMCPY(y_out__, y, doublecomplex, NA_TOTAL(rblapack_y));
  rblapack_y = rblapack_y_out__;
  y = y_out__;

  zla_wwaddw_(&n, x, y, w);

  return rb_ary_new3(2, rblapack_x, rblapack_y);
}

void
init_lapack_zla_wwaddw(VALUE mLapack){
  rb_define_module_function(mLapack, "zla_wwaddw", rblapack_zla_wwaddw, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
