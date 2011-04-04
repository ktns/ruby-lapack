#include "rb_lapack.h"

extern integer icmax1_(integer *n, complex *cx, integer *incx);

static VALUE
rb_icmax1(int argc, VALUE *argv, VALUE self){
  VALUE rb_cx;
  complex *cx; 
  VALUE rb_incx;
  integer incx; 
  VALUE rb___out__;
  integer __out__; 

  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.icmax1( cx, incx)\n    or\n  NumRu::Lapack.icmax1  # print help\n\n\nFORTRAN MANUAL\n      INTEGER          FUNCTION ICMAX1( N, CX, INCX )\n\n*  Purpose\n*  =======\n*\n*  ICMAX1 finds the index of the element whose real part has maximum\n*  absolute value.\n*\n*  Based on ICAMAX from Level 1 BLAS.\n*  The change is to use the 'genuine' absolute value.\n*\n*  Contributed by Nick Higham for use with CLACON.\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The number of elements in the vector CX.\n*\n*  CX      (input) COMPLEX array, dimension (N)\n*          The vector whose elements will be summed.\n*\n*  INCX    (input) INTEGER\n*          The spacing between successive values of CX.  INCX >= 1.\n*\n\n* =====================================================================\n*\n*     .. Local Scalars ..\n      INTEGER            I, IX\n      REAL               SMAX\n      COMPLEX            ZDUM\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          ABS\n*     ..\n*     .. Statement Functions ..\n      REAL               CABS1\n*     ..\n*     .. Statement Function definitions ..\n*\n*     NEXT LINE IS THE ONLY MODIFICATION.\n      CABS1( ZDUM ) = ABS( ZDUM )\n*     ..\n\n");
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
  if (NA_TYPE(rb_cx) != NA_SCOMPLEX)
    rb_cx = na_change_type(rb_cx, NA_SCOMPLEX);
  cx = NA_PTR_TYPE(rb_cx, complex*);

  __out__ = icmax1_(&n, cx, &incx);

  rb___out__ = INT2NUM(__out__);
  return rb___out__;
}

void
init_lapack_icmax1(VALUE mLapack){
  rb_define_module_function(mLapack, "icmax1", rb_icmax1, -1);
}
