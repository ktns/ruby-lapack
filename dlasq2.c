#include "rb_lapack.h"

extern VOID dlasq2_(integer *n, doublereal *z, integer *info);

static VALUE
rb_dlasq2(int argc, VALUE *argv, VALUE self){
  VALUE rb_n;
  integer n; 
  VALUE rb_z;
  doublereal *z; 
  VALUE rb_info;
  integer info; 
  VALUE rb_z_out__;
  doublereal *z_out__;


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, z = NumRu::Lapack.dlasq2( n, z)\n    or\n  NumRu::Lapack.dlasq2  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rb_n = argv[0];
  rb_z = argv[1];

  n = NUM2INT(rb_n);
  if (!NA_IsNArray(rb_z))
    rb_raise(rb_eArgError, "z (2th argument) must be NArray");
  if (NA_RANK(rb_z) != 1)
    rb_raise(rb_eArgError, "rank of z (2th argument) must be %d", 1);
  if (NA_SHAPE0(rb_z) != (4*n))
    rb_raise(rb_eRuntimeError, "shape 0 of z must be %d", 4*n);
  if (NA_TYPE(rb_z) != NA_DFLOAT)
    rb_z = na_change_type(rb_z, NA_DFLOAT);
  z = NA_PTR_TYPE(rb_z, doublereal*);
  {
    int shape[1];
    shape[0] = 4*n;
    rb_z_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  z_out__ = NA_PTR_TYPE(rb_z_out__, doublereal*);
  MEMCPY(z_out__, z, doublereal, NA_TOTAL(rb_z));
  rb_z = rb_z_out__;
  z = z_out__;

  dlasq2_(&n, z, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_info, rb_z);
}

void
init_lapack_dlasq2(VALUE mLapack){
  rb_define_module_function(mLapack, "dlasq2", rb_dlasq2, -1);
}
