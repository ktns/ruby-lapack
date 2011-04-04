#include "rb_lapack.h"

extern VOID stftri_(char *transr, char *uplo, char *diag, integer *n, real *a, integer *info);

static VALUE
rb_stftri(int argc, VALUE *argv, VALUE self){
  VALUE rb_transr;
  char transr; 
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_diag;
  char diag; 
  VALUE rb_n;
  integer n; 
  VALUE rb_a;
  real *a; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  real *a_out__;


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, a = NumRu::Lapack.stftri( transr, uplo, diag, n, a)\n    or\n  NumRu::Lapack.stftri  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_transr = argv[0];
  rb_uplo = argv[1];
  rb_diag = argv[2];
  rb_n = argv[3];
  rb_a = argv[4];

  transr = StringValueCStr(rb_transr)[0];
  diag = StringValueCStr(rb_diag)[0];
  uplo = StringValueCStr(rb_uplo)[0];
  n = NUM2INT(rb_n);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (5th argument) must be NArray");
  if (NA_RANK(rb_a) != 1)
    rb_raise(rb_eArgError, "rank of a (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_a) != (n*(n+1)/2))
    rb_raise(rb_eRuntimeError, "shape 0 of a must be %d", n*(n+1)/2);
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);
  {
    int shape[1];
    shape[0] = n*(n+1)/2;
    rb_a_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rb_a_out__, real*);
  MEMCPY(a_out__, a, real, NA_TOTAL(rb_a));
  rb_a = rb_a_out__;
  a = a_out__;

  stftri_(&transr, &uplo, &diag, &n, a, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_info, rb_a);
}

void
init_lapack_stftri(VALUE mLapack){
  rb_define_module_function(mLapack, "stftri", rb_stftri, -1);
}
