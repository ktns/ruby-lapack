#include "rb_lapack.h"

extern VOID slasq2_(integer *n, real *z, integer *info);

static VALUE
rb_slasq2(int argc, VALUE *argv, VALUE self){
  VALUE rb_n;
  integer n; 
  VALUE rb_z;
  real *z; 
  VALUE rb_info;
  integer info; 
  VALUE rb_z_out__;
  real *z_out__;


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, z = NumRu::Lapack.slasq2( n, z)\n    or\n  NumRu::Lapack.slasq2  # print help\n\n\nFORTRAN MANUAL\n\n");
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
  if (NA_TYPE(rb_z) != NA_SFLOAT)
    rb_z = na_change_type(rb_z, NA_SFLOAT);
  z = NA_PTR_TYPE(rb_z, real*);
  {
    int shape[1];
    shape[0] = 4*n;
    rb_z_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  z_out__ = NA_PTR_TYPE(rb_z_out__, real*);
  MEMCPY(z_out__, z, real, NA_TOTAL(rb_z));
  rb_z = rb_z_out__;
  z = z_out__;

  slasq2_(&n, z, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_info, rb_z);
}

void
init_lapack_slasq2(VALUE mLapack){
  rb_define_module_function(mLapack, "slasq2", rb_slasq2, -1);
}
