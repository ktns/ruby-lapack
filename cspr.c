#include "rb_lapack.h"

extern VOID cspr_(char *uplo, integer *n, complex *alpha, complex *x, integer *incx, complex *ap);

static VALUE
rb_cspr(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_n;
  integer n; 
  VALUE rb_alpha;
  complex alpha; 
  VALUE rb_x;
  complex *x; 
  VALUE rb_incx;
  integer incx; 
  VALUE rb_ap;
  complex *ap; 
  VALUE rb_ap_out__;
  complex *ap_out__;


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  ap = NumRu::Lapack.cspr( uplo, n, alpha, x, incx, ap)\n    or\n  NumRu::Lapack.cspr  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rb_uplo = argv[0];
  rb_n = argv[1];
  rb_alpha = argv[2];
  rb_x = argv[3];
  rb_incx = argv[4];
  rb_ap = argv[5];

  uplo = StringValueCStr(rb_uplo)[0];
  n = NUM2INT(rb_n);
  alpha.r = (real)NUM2DBL(rb_funcall(rb_alpha, rb_intern("real"), 0));
  alpha.i = (real)NUM2DBL(rb_funcall(rb_alpha, rb_intern("imag"), 0));
  incx = NUM2INT(rb_incx);
  if (!NA_IsNArray(rb_ap))
    rb_raise(rb_eArgError, "ap (6th argument) must be NArray");
  if (NA_RANK(rb_ap) != 1)
    rb_raise(rb_eArgError, "rank of ap (6th argument) must be %d", 1);
  if (NA_SHAPE0(rb_ap) != (( n*( n + 1 ) )/2))
    rb_raise(rb_eRuntimeError, "shape 0 of ap must be %d", ( n*( n + 1 ) )/2);
  if (NA_TYPE(rb_ap) != NA_SCOMPLEX)
    rb_ap = na_change_type(rb_ap, NA_SCOMPLEX);
  ap = NA_PTR_TYPE(rb_ap, complex*);
  if (!NA_IsNArray(rb_x))
    rb_raise(rb_eArgError, "x (4th argument) must be NArray");
  if (NA_RANK(rb_x) != 1)
    rb_raise(rb_eArgError, "rank of x (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_x) != (1 + ( n - 1 )*abs( incx )))
    rb_raise(rb_eRuntimeError, "shape 0 of x must be %d", 1 + ( n - 1 )*abs( incx ));
  if (NA_TYPE(rb_x) != NA_SCOMPLEX)
    rb_x = na_change_type(rb_x, NA_SCOMPLEX);
  x = NA_PTR_TYPE(rb_x, complex*);
  {
    int shape[1];
    shape[0] = ( n*( n + 1 ) )/2;
    rb_ap_out__ = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  ap_out__ = NA_PTR_TYPE(rb_ap_out__, complex*);
  MEMCPY(ap_out__, ap, complex, NA_TOTAL(rb_ap));
  rb_ap = rb_ap_out__;
  ap = ap_out__;

  cspr_(&uplo, &n, &alpha, x, &incx, ap);

  return rb_ap;
}

void
init_lapack_cspr(VALUE mLapack){
  rb_define_module_function(mLapack, "cspr", rb_cspr, -1);
}
