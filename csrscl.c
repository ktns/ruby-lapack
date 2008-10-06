#include "rb_lapack.h"

static VALUE
rb_csrscl(int argc, VALUE *argv, VALUE self){
  VALUE rb_n;
  integer n; 
  VALUE rb_sa;
  real sa; 
  VALUE rb_sx;
  complex *sx; 
  VALUE rb_incx;
  integer incx; 
  VALUE rb_sx_out__;
  complex *sx_out__;


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  sx = NumRu::Lapack.csrscl( n, sa, sx, incx)\n    or\n  NumRu::Lapack.csrscl  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE CSRSCL( N, SA, SX, INCX )\n\n*  Purpose\n*  =======\n*\n*  CSRSCL multiplies an n-element complex vector x by the real scalar\n*  1/a.  This is done without overflow or underflow as long as\n*  the final result x/a does not overflow or underflow.\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The number of components of the vector x.\n*\n*  SA      (input) REAL\n*          The scalar a which is used to divide each component of x.\n*          SA must be >= 0, or the subroutine will divide by zero.\n*\n*  SX      (input/output) COMPLEX array, dimension\n*                         (1+(N-1)*abs(INCX))\n*          The n-element vector x.\n*\n*  INCX    (input) INTEGER\n*          The increment between successive values of the vector SX.\n*          > 0:  SX(1) = X(1) and SX(1+(i-1)*INCX) = x(i),     1< i<= n\n*\n\n* =====================================================================\n*\n\n");
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
  if (NA_TYPE(rb_sx) != NA_SCOMPLEX)
    rb_sx = na_change_type(rb_sx, NA_SCOMPLEX);
  sx = NA_PTR_TYPE(rb_sx, complex*);
  {
    int shape[1];
    shape[0] = 1+(n-1)*abs(incx);
    rb_sx_out__ = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  sx_out__ = NA_PTR_TYPE(rb_sx_out__, complex*);
  MEMCPY(sx_out__, sx, complex, NA_TOTAL(rb_sx));
  rb_sx = rb_sx_out__;
  sx = sx_out__;

  csrscl_(&n, &sa, sx, &incx);

  return rb_sx;
}

void
init_lapack_csrscl(VALUE mLapack){
  rb_define_module_function(mLapack, "csrscl", rb_csrscl, -1);
}
