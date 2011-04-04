#include "rb_lapack.h"

extern VOID dlasrt_(char *id, integer *n, doublereal *d, integer *info);

static VALUE
rb_dlasrt(int argc, VALUE *argv, VALUE self){
  VALUE rb_id;
  char id; 
  VALUE rb_d;
  doublereal *d; 
  VALUE rb_info;
  integer info; 
  VALUE rb_d_out__;
  doublereal *d_out__;

  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, d = NumRu::Lapack.dlasrt( id, d)\n    or\n  NumRu::Lapack.dlasrt  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rb_id = argv[0];
  rb_d = argv[1];

  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (2th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (2th argument) must be %d", 1);
  n = NA_SHAPE0(rb_d);
  if (NA_TYPE(rb_d) != NA_DFLOAT)
    rb_d = na_change_type(rb_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rb_d, doublereal*);
  id = StringValueCStr(rb_id)[0];
  {
    int shape[1];
    shape[0] = n;
    rb_d_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  d_out__ = NA_PTR_TYPE(rb_d_out__, doublereal*);
  MEMCPY(d_out__, d, doublereal, NA_TOTAL(rb_d));
  rb_d = rb_d_out__;
  d = d_out__;

  dlasrt_(&id, &n, d, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_info, rb_d);
}

void
init_lapack_dlasrt(VALUE mLapack){
  rb_define_module_function(mLapack, "dlasrt", rb_dlasrt, -1);
}
