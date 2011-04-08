#include "rb_lapack.h"

extern VOID crot_(integer *n, complex *cx, integer *incx, complex *cy, integer *incy, real *c, complex *s);

static VALUE sHelp, sUsage;

static VALUE
rblapack_crot(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_cx;
  complex *cx; 
  VALUE rblapack_incx;
  integer incx; 
  VALUE rblapack_cy;
  complex *cy; 
  VALUE rblapack_incy;
  integer incy; 
  VALUE rblapack_c;
  real c; 
  VALUE rblapack_s;
  complex s; 
  VALUE rblapack_cx_out__;
  complex *cx_out__;
  VALUE rblapack_cy_out__;
  complex *cy_out__;

  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  cx, cy = NumRu::Lapack.crot( cx, incx, cy, incy, c, s, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE CROT( N, CX, INCX, CY, INCY, C, S )\n\n*  Purpose\n*  =======\n*\n*  CROT   applies a plane rotation, where the cos (C) is real and the\n*  sin (S) is complex, and the vectors CX and CY are complex.\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The number of elements in the vectors CX and CY.\n*\n*  CX      (input/output) COMPLEX array, dimension (N)\n*          On input, the vector X.\n*          On output, CX is overwritten with C*X + S*Y.\n*\n*  INCX    (input) INTEGER\n*          The increment between successive values of CY.  INCX <> 0.\n*\n*  CY      (input/output) COMPLEX array, dimension (N)\n*          On input, the vector Y.\n*          On output, CY is overwritten with -CONJG(S)*X + C*Y.\n*\n*  INCY    (input) INTEGER\n*          The increment between successive values of CY.  INCX <> 0.\n*\n*  C       (input) REAL\n*  S       (input) COMPLEX\n*          C and S define a rotation\n*             [  C          S  ]\n*             [ -conjg(S)   C  ]\n*          where C*C + S*CONJG(S) = 1.0.\n*\n\n* =====================================================================\n*\n*     .. Local Scalars ..\n      INTEGER            I, IX, IY\n      COMPLEX            STEMP\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          CONJG\n*     ..\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  cx, cy = NumRu::Lapack.crot( cx, incx, cy, incy, c, s, [:usage => usage, :help => help])\n");
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
  if (NA_TYPE(rblapack_cy) != NA_SCOMPLEX)
    rblapack_cy = na_change_type(rblapack_cy, NA_SCOMPLEX);
  cy = NA_PTR_TYPE(rblapack_cy, complex*);
  c = (real)NUM2DBL(rblapack_c);
  incx = NUM2INT(rblapack_incx);
  incy = NUM2INT(rblapack_incy);
  s.r = (real)NUM2DBL(rb_funcall(rblapack_s, rb_intern("real"), 0));
  s.i = (real)NUM2DBL(rb_funcall(rblapack_s, rb_intern("imag"), 0));
  if (!NA_IsNArray(rblapack_cx))
    rb_raise(rb_eArgError, "cx (1th argument) must be NArray");
  if (NA_RANK(rblapack_cx) != 1)
    rb_raise(rb_eArgError, "rank of cx (1th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_cx) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of cx must be the same as shape 0 of cy");
  if (NA_TYPE(rblapack_cx) != NA_SCOMPLEX)
    rblapack_cx = na_change_type(rblapack_cx, NA_SCOMPLEX);
  cx = NA_PTR_TYPE(rblapack_cx, complex*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_cx_out__ = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  cx_out__ = NA_PTR_TYPE(rblapack_cx_out__, complex*);
  MEMCPY(cx_out__, cx, complex, NA_TOTAL(rblapack_cx));
  rblapack_cx = rblapack_cx_out__;
  cx = cx_out__;
  {
    int shape[1];
    shape[0] = n;
    rblapack_cy_out__ = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  cy_out__ = NA_PTR_TYPE(rblapack_cy_out__, complex*);
  MEMCPY(cy_out__, cy, complex, NA_TOTAL(rblapack_cy));
  rblapack_cy = rblapack_cy_out__;
  cy = cy_out__;

  crot_(&n, cx, &incx, cy, &incy, &c, &s);

  return rb_ary_new3(2, rblapack_cx, rblapack_cy);
}

void
init_lapack_crot(VALUE mLapack){
  rb_define_module_function(mLapack, "crot", rblapack_crot, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
