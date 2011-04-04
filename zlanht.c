#include "rb_lapack.h"

extern doublereal zlanht_(char *norm, integer *n, doublereal *d, doublecomplex *e);

static VALUE
rb_zlanht(int argc, VALUE *argv, VALUE self){
  VALUE rb_norm;
  char norm; 
  VALUE rb_d;
  doublereal *d; 
  VALUE rb_e;
  doublecomplex *e; 
  VALUE rb___out__;
  doublereal __out__; 

  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.zlanht( norm, d, e)\n    or\n  NumRu::Lapack.zlanht  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_norm = argv[0];
  rb_d = argv[1];
  rb_e = argv[2];

  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (2th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (2th argument) must be %d", 1);
  n = NA_SHAPE0(rb_d);
  if (NA_TYPE(rb_d) != NA_DFLOAT)
    rb_d = na_change_type(rb_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rb_d, doublereal*);
  norm = StringValueCStr(rb_norm)[0];
  if (!NA_IsNArray(rb_e))
    rb_raise(rb_eArgError, "e (3th argument) must be NArray");
  if (NA_RANK(rb_e) != 1)
    rb_raise(rb_eArgError, "rank of e (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_e) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of e must be %d", n-1);
  if (NA_TYPE(rb_e) != NA_DCOMPLEX)
    rb_e = na_change_type(rb_e, NA_DCOMPLEX);
  e = NA_PTR_TYPE(rb_e, doublecomplex*);

  __out__ = zlanht_(&norm, &n, d, e);

  rb___out__ = rb_float_new((double)__out__);
  return rb___out__;
}

void
init_lapack_zlanht(VALUE mLapack){
  rb_define_module_function(mLapack, "zlanht", rb_zlanht, -1);
}
