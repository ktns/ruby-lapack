#include "rb_lapack.h"

extern VOID drscl_(integer *n, doublereal *sa, doublereal *sx, integer *incx);

static VALUE
rb_drscl(int argc, VALUE *argv, VALUE self){
  VALUE rb_n;
  integer n; 
  VALUE rb_sa;
  doublereal sa; 
  VALUE rb_sx;
  doublereal *sx; 
  VALUE rb_incx;
  integer incx; 
  VALUE rb_sx_out__;
  doublereal *sx_out__;


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  sx = NumRu::Lapack.drscl( n, sa, sx, incx)\n    or\n  NumRu::Lapack.drscl  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DRSCL( N, SA, SX, INCX )\n\n*  Purpose\n*  =======\n*\n*  DRSCL multiplies an n-element real vector x by the real scalar 1/a.\n*  This is done without overflow or underflow as long as\n*  the final result x/a does not overflow or underflow.\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The number of components of the vector x.\n*\n*  SA      (input) DOUBLE PRECISION\n*          The scalar a which is used to divide each component of x.\n*          SA must be >= 0, or the subroutine will divide by zero.\n*\n*  SX      (input/output) DOUBLE PRECISION array, dimension\n*                         (1+(N-1)*abs(INCX))\n*          The n-element vector x.\n*\n*  INCX    (input) INTEGER\n*          The increment between successive values of the vector SX.\n*          > 0:  SX(1) = X(1) and SX(1+(i-1)*INCX) = x(i),     1< i<= n\n*\n\n* =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_n = argv[0];
  rb_sa = argv[1];
  rb_sx = argv[2];
  rb_incx = argv[3];

  sa = NUM2DBL(rb_sa);
  n = NUM2INT(rb_n);
  incx = NUM2INT(rb_incx);
  if (!NA_IsNArray(rb_sx))
    rb_raise(rb_eArgError, "sx (3th argument) must be NArray");
  if (NA_RANK(rb_sx) != 1)
    rb_raise(rb_eArgError, "rank of sx (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_sx) != (1+(n-1)*abs(incx)))
    rb_raise(rb_eRuntimeError, "shape 0 of sx must be %d", 1+(n-1)*abs(incx));
  if (NA_TYPE(rb_sx) != NA_DFLOAT)
    rb_sx = na_change_type(rb_sx, NA_DFLOAT);
  sx = NA_PTR_TYPE(rb_sx, doublereal*);
  {
    int shape[1];
    shape[0] = 1+(n-1)*abs(incx);
    rb_sx_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  sx_out__ = NA_PTR_TYPE(rb_sx_out__, doublereal*);
  MEMCPY(sx_out__, sx, doublereal, NA_TOTAL(rb_sx));
  rb_sx = rb_sx_out__;
  sx = sx_out__;

  drscl_(&n, &sa, sx, &incx);

  return rb_sx;
}

void
init_lapack_drscl(VALUE mLapack){
  rb_define_module_function(mLapack, "drscl", rb_drscl, -1);
}
