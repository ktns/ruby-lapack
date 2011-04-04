#include "rb_lapack.h"

extern VOID cgttrf_(integer *n, complex *dl, complex *d, complex *du, complex *du2, integer *ipiv, integer *info);

static VALUE
rb_cgttrf(int argc, VALUE *argv, VALUE self){
  VALUE rb_dl;
  complex *dl; 
  VALUE rb_d;
  complex *d; 
  VALUE rb_du;
  complex *du; 
  VALUE rb_du2;
  complex *du2; 
  VALUE rb_ipiv;
  integer *ipiv; 
  VALUE rb_info;
  integer info; 
  VALUE rb_dl_out__;
  complex *dl_out__;
  VALUE rb_d_out__;
  complex *d_out__;
  VALUE rb_du_out__;
  complex *du_out__;

  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  du2, ipiv, info, dl, d, du = NumRu::Lapack.cgttrf( dl, d, du)\n    or\n  NumRu::Lapack.cgttrf  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_dl = argv[0];
  rb_d = argv[1];
  rb_du = argv[2];

  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (2th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (2th argument) must be %d", 1);
  n = NA_SHAPE0(rb_d);
  if (NA_TYPE(rb_d) != NA_SCOMPLEX)
    rb_d = na_change_type(rb_d, NA_SCOMPLEX);
  d = NA_PTR_TYPE(rb_d, complex*);
  if (!NA_IsNArray(rb_du))
    rb_raise(rb_eArgError, "du (3th argument) must be NArray");
  if (NA_RANK(rb_du) != 1)
    rb_raise(rb_eArgError, "rank of du (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_du) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of du must be %d", n-1);
  if (NA_TYPE(rb_du) != NA_SCOMPLEX)
    rb_du = na_change_type(rb_du, NA_SCOMPLEX);
  du = NA_PTR_TYPE(rb_du, complex*);
  if (!NA_IsNArray(rb_dl))
    rb_raise(rb_eArgError, "dl (1th argument) must be NArray");
  if (NA_RANK(rb_dl) != 1)
    rb_raise(rb_eArgError, "rank of dl (1th argument) must be %d", 1);
  if (NA_SHAPE0(rb_dl) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of dl must be %d", n-1);
  if (NA_TYPE(rb_dl) != NA_SCOMPLEX)
    rb_dl = na_change_type(rb_dl, NA_SCOMPLEX);
  dl = NA_PTR_TYPE(rb_dl, complex*);
  {
    int shape[1];
    shape[0] = n-2;
    rb_du2 = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  du2 = NA_PTR_TYPE(rb_du2, complex*);
  {
    int shape[1];
    shape[0] = n;
    rb_ipiv = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  ipiv = NA_PTR_TYPE(rb_ipiv, integer*);
  {
    int shape[1];
    shape[0] = n-1;
    rb_dl_out__ = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  dl_out__ = NA_PTR_TYPE(rb_dl_out__, complex*);
  MEMCPY(dl_out__, dl, complex, NA_TOTAL(rb_dl));
  rb_dl = rb_dl_out__;
  dl = dl_out__;
  {
    int shape[1];
    shape[0] = n;
    rb_d_out__ = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  d_out__ = NA_PTR_TYPE(rb_d_out__, complex*);
  MEMCPY(d_out__, d, complex, NA_TOTAL(rb_d));
  rb_d = rb_d_out__;
  d = d_out__;
  {
    int shape[1];
    shape[0] = n-1;
    rb_du_out__ = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  du_out__ = NA_PTR_TYPE(rb_du_out__, complex*);
  MEMCPY(du_out__, du, complex, NA_TOTAL(rb_du));
  rb_du = rb_du_out__;
  du = du_out__;

  cgttrf_(&n, dl, d, du, du2, ipiv, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(6, rb_du2, rb_ipiv, rb_info, rb_dl, rb_d, rb_du);
}

void
init_lapack_cgttrf(VALUE mLapack){
  rb_define_module_function(mLapack, "cgttrf", rb_cgttrf, -1);
}
