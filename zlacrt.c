#include "rb_lapack.h"

extern VOID zlacrt_(integer *n, doublecomplex *cx, integer *incx, doublecomplex *cy, integer *incy, doublecomplex *c, doublecomplex *s);

static VALUE sHelp, sUsage;

static VALUE
rblapack_zlacrt(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_cx;
  doublecomplex *cx; 
  VALUE rblapack_incx;
  integer incx; 
  VALUE rblapack_cy;
  doublecomplex *cy; 
  VALUE rblapack_incy;
  integer incy; 
  VALUE rblapack_c;
  doublecomplex c; 
  VALUE rblapack_s;
  doublecomplex s; 
  VALUE rblapack_cx_out__;
  doublecomplex *cx_out__;
  VALUE rblapack_cy_out__;
  doublecomplex *cy_out__;

  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  cx, cy = NumRu::Lapack.zlacrt( cx, incx, cy, incy, c, s, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZLACRT( N, CX, INCX, CY, INCY, C, S )\n\n*  Purpose\n*  =======\n*\n*  ZLACRT performs the operation\n*\n*     (  c  s )( x )  ==> ( x )\n*     ( -s  c )( y )      ( y )\n*\n*  where c and s are complex and the vectors x and y are complex.\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The number of elements in the vectors CX and CY.\n*\n*  CX      (input/output) COMPLEX*16 array, dimension (N)\n*          On input, the vector x.\n*          On output, CX is overwritten with c*x + s*y.\n*\n*  INCX    (input) INTEGER\n*          The increment between successive values of CX.  INCX <> 0.\n*\n*  CY      (input/output) COMPLEX*16 array, dimension (N)\n*          On input, the vector y.\n*          On output, CY is overwritten with -s*x + c*y.\n*\n*  INCY    (input) INTEGER\n*          The increment between successive values of CY.  INCY <> 0.\n*\n*  C       (input) COMPLEX*16\n*  S       (input) COMPLEX*16\n*          C and S define the matrix\n*             [  C   S  ].\n*             [ -S   C  ]\n*\n\n* =====================================================================\n*\n*     .. Local Scalars ..\n      INTEGER            I, IX, IY\n      COMPLEX*16         CTEMP\n*     ..\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  cx, cy = NumRu::Lapack.zlacrt( cx, incx, cy, incy, c, s, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rblapack_cx = argv[0];
  rblapack_incx = argv[1];
  rblapack_cy = argv[2];
  rblapack_incy = argv[3];
  rblapack_c = argv[4];
  rblapack_s = argv[5];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_cy))
    rb_raise(rb_eArgError, "cy (3th argument) must be NArray");
  if (NA_RANK(rblapack_cy) != 1)
    rb_raise(rb_eArgError, "rank of cy (3th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_cy);
  if (NA_TYPE(rblapack_cy) != NA_DCOMPLEX)
    rblapack_cy = na_change_type(rblapack_cy, NA_DCOMPLEX);
  cy = NA_PTR_TYPE(rblapack_cy, doublecomplex*);
  c.r = NUM2DBL(rb_funcall(rblapack_c, rb_intern("real"), 0));
  c.i = NUM2DBL(rb_funcall(rblapack_c, rb_intern("imag"), 0));
  incx = NUM2INT(rblapack_incx);
  incy = NUM2INT(rblapack_incy);
  s.r = NUM2DBL(rb_funcall(rblapack_s, rb_intern("real"), 0));
  s.i = NUM2DBL(rb_funcall(rblapack_s, rb_intern("imag"), 0));
  if (!NA_IsNArray(rblapack_cx))
    rb_raise(rb_eArgError, "cx (1th argument) must be NArray");
  if (NA_RANK(rblapack_cx) != 1)
    rb_raise(rb_eArgError, "rank of cx (1th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_cx) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of cx must be the same as shape 0 of cy");
  if (NA_TYPE(rblapack_cx) != NA_DCOMPLEX)
    rblapack_cx = na_change_type(rblapack_cx, NA_DCOMPLEX);
  cx = NA_PTR_TYPE(rblapack_cx, doublecomplex*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_cx_out__ = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  cx_out__ = NA_PTR_TYPE(rblapack_cx_out__, doublecomplex*);
  MEMCPY(cx_out__, cx, doublecomplex, NA_TOTAL(rblapack_cx));
  rblapack_cx = rblapack_cx_out__;
  cx = cx_out__;
  {
    int shape[1];
    shape[0] = n;
    rblapack_cy_out__ = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  cy_out__ = NA_PTR_TYPE(rblapack_cy_out__, doublecomplex*);
  MEMCPY(cy_out__, cy, doublecomplex, NA_TOTAL(rblapack_cy));
  rblapack_cy = rblapack_cy_out__;
  cy = cy_out__;

  zlacrt_(&n, cx, &incx, cy, &incy, &c, &s);

  return rb_ary_new3(2, rblapack_cx, rblapack_cy);
}

void
init_lapack_zlacrt(VALUE mLapack){
  rb_define_module_function(mLapack, "zlacrt", rblapack_zlacrt, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
