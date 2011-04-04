#include "rb_lapack.h"

extern doublereal dzsum1_(integer *n, doublecomplex *cx, integer *incx);

static VALUE
rb_dzsum1(int argc, VALUE *argv, VALUE self){
  VALUE rb_cx;
  doublecomplex *cx; 
  VALUE rb_incx;
  integer incx; 
  VALUE rb___out__;
  doublereal __out__; 

  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.dzsum1( cx, incx)\n    or\n  NumRu::Lapack.dzsum1  # print help\n\n\nFORTRAN MANUAL\n      DOUBLE PRECISION FUNCTION DZSUM1( N, CX, INCX )\n\n*  Purpose\n*  =======\n*\n*  DZSUM1 takes the sum of the absolute values of a complex\n*  vector and returns a double precision result.\n*\n*  Based on DZASUM from the Level 1 BLAS.\n*  The change is to use the 'genuine' absolute value.\n*\n*  Contributed by Nick Higham for use with ZLACON.\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The number of elements in the vector CX.\n*\n*  CX      (input) COMPLEX*16 array, dimension (N)\n*          The vector whose elements will be summed.\n*\n*  INCX    (input) INTEGER\n*          The spacing between successive values of CX.  INCX > 0.\n*\n\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      INTEGER            I, NINCX\n      DOUBLE PRECISION   STEMP\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          ABS\n*     ..\n\n");
    return Qnil;
  }
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rb_cx = argv[0];
  rb_incx = argv[1];

  incx = NUM2INT(rb_incx);
  if (!NA_IsNArray(rb_cx))
    rb_raise(rb_eArgError, "cx (1th argument) must be NArray");
  if (NA_RANK(rb_cx) != 1)
    rb_raise(rb_eArgError, "rank of cx (1th argument) must be %d", 1);
  n = NA_SHAPE0(rb_cx);
  if (NA_TYPE(rb_cx) != NA_DCOMPLEX)
    rb_cx = na_change_type(rb_cx, NA_DCOMPLEX);
  cx = NA_PTR_TYPE(rb_cx, doublecomplex*);

  __out__ = dzsum1_(&n, cx, &incx);

  rb___out__ = rb_float_new((double)__out__);
  return rb___out__;
}

void
init_lapack_dzsum1(VALUE mLapack){
  rb_define_module_function(mLapack, "dzsum1", rb_dzsum1, -1);
}
