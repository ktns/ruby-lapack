#include "rb_lapack.h"

static VALUE
rb_srscl(int argc, VALUE *argv, VALUE self){
  VALUE rb_n;
  integer n; 
  VALUE rb_sa;
  real sa; 
  VALUE rb_sx;
  real *sx; 
  VALUE rb_incx;
  integer incx; 
  VALUE rb_sx_out__;
  real *sx_out__;


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  sx = NumRu::Lapack.srscl( n, sa, sx, incx)\n    or\n  NumRu::Lapack.srscl  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SRSCL( N, SA, SX, INCX )\n\n*  Purpose\n*  =======\n*\n*  SRSCL multiplies an n-element real vector x by the real scalar 1/a.\n*  This is done without overflow or underflow as long as\n*  the final result x/a does not overflow or underflow.\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The number of components of the vector x.\n*\n*  SA      (input) REAL\n*          The scalar a which is used to divide each component of x.\n*          SA must be >= 0, or the subroutine will divide by zero.\n*\n*  SX      (input/output) REAL array, dimension\n*                         (1+(N-1)*abs(INCX))\n*          The n-element vector x.\n*\n*  INCX    (input) INTEGER\n*          The increment between successive values of the vector SX.\n*          > 0:  SX(1) = X(1) and SX(1+(i-1)*INCX) = x(i),     1< i<= n\n*\n\n* =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_n = argv[0];
  rb_sa = argv[1];
  rb_sx = argv[2];
  rb_incx = argv[3];

  n = NUM2INT(rb_n);
  sa = (real)NUM2DBL(rb_sa);
  incx = NUM2INT(rb_incx);
  if (!NA_IsNArray(rb_sx))
    rb_raise(rb_eArgError, "sx (3th argument) must be NArray");
  if (NA_RANK(rb_sx) != 1)
    rb_raise(rb_eArgError, "rank of sx (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_sx) != (1+(n-1)*abs(incx)))
    rb_raise(rb_eRuntimeError, "shape 0 of sx must be %d", 1+(n-1)*abs(incx));
  if (NA_TYPE(rb_sx) != NA_SFLOAT)
    rb_sx = na_change_type(rb_sx, NA_SFLOAT);
  sx = NA_PTR_TYPE(rb_sx, real*);
  {
    int shape[1];
    shape[0] = 1+(n-1)*abs(incx);
    rb_sx_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  sx_out__ = NA_PTR_TYPE(rb_sx_out__, real*);
  MEMCPY(sx_out__, sx, real, NA_TOTAL(rb_sx));
  rb_sx = rb_sx_out__;
  sx = sx_out__;

  srscl_(&n, &sa, sx, &incx);

  return rb_sx;
}

void
init_lapack_srscl(VALUE mLapack){
  rb_define_module_function(mLapack, "srscl", rb_srscl, -1);
}
