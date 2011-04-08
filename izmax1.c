#include "rb_lapack.h"

extern integer izmax1_(integer *n, doublecomplex *cx, integer *incx);

static VALUE sHelp, sUsage;

static VALUE
rblapack_izmax1(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_cx;
  doublecomplex *cx; 
  VALUE rblapack_incx;
  integer incx; 
  VALUE rblapack___out__;
  integer __out__; 

  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.izmax1( cx, incx, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      INTEGER          FUNCTION IZMAX1( N, CX, INCX )\n\n*  Purpose\n*  =======\n*\n*  IZMAX1 finds the index of the element whose real part has maximum\n*  absolute value.\n*\n*  Based on IZAMAX from Level 1 BLAS.\n*  The change is to use the 'genuine' absolute value.\n*\n*  Contributed by Nick Higham for use with ZLACON.\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The number of elements in the vector CX.\n*\n*  CX      (input) COMPLEX*16 array, dimension (N)\n*          The vector whose elements will be summed.\n*\n*  INCX    (input) INTEGER\n*          The spacing between successive values of CX.  INCX >= 1.\n*\n\n* =====================================================================\n*\n*     .. Local Scalars ..\n      INTEGER            I, IX\n      DOUBLE PRECISION   SMAX\n      COMPLEX*16         ZDUM\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          ABS\n*     ..\n*     .. Statement Functions ..\n      DOUBLE PRECISION   CABS1\n*     ..\n*     .. Statement Function definitions ..\n*\n*     NEXT LINE IS THE ONLY MODIFICATION.\n      CABS1( ZDUM ) = ABS( ZDUM )\n*     ..\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.izmax1( cx, incx, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rblapack_cx = argv[0];
  rblapack_incx = argv[1];
  if (rb_options != Qnil) {
  }

  incx = NUM2INT(rblapack_incx);
  if (!NA_IsNArray(rblapack_cx))
    rb_raise(rb_eArgError, "cx (1th argument) must be NArray");
  if (NA_RANK(rblapack_cx) != 1)
    rb_raise(rb_eArgError, "rank of cx (1th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_cx);
  if (NA_TYPE(rblapack_cx) != NA_DCOMPLEX)
    rblapack_cx = na_change_type(rblapack_cx, NA_DCOMPLEX);
  cx = NA_PTR_TYPE(rblapack_cx, doublecomplex*);

  __out__ = izmax1_(&n, cx, &incx);

  rblapack___out__ = INT2NUM(__out__);
  return rblapack___out__;
}

void
init_lapack_izmax1(VALUE mLapack){
  rb_define_module_function(mLapack, "izmax1", rblapack_izmax1, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
