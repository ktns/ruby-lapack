#include "rb_lapack.h"

extern VOID slamrg_(integer *n1, integer *n2, real *a, integer *strd1, integer *strd2, integer *index);

static VALUE
rb_slamrg(int argc, VALUE *argv, VALUE self){
  VALUE rb_n1;
  integer n1; 
  VALUE rb_n2;
  integer n2; 
  VALUE rb_a;
  real *a; 
  VALUE rb_strd1;
  integer strd1; 
  VALUE rb_strd2;
  integer strd2; 
  VALUE rb_index;
  integer *index; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  index = NumRu::Lapack.slamrg( n1, n2, a, strd1, strd2)\n    or\n  NumRu::Lapack.slamrg  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_n1 = argv[0];
  rb_n2 = argv[1];
  rb_a = argv[2];
  rb_strd1 = argv[3];
  rb_strd2 = argv[4];

  strd1 = NUM2INT(rb_strd1);
  n2 = NUM2INT(rb_n2);
  n1 = NUM2INT(rb_n1);
  strd2 = NUM2INT(rb_strd2);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rb_a) != 1)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_a) != (n1+n2))
    rb_raise(rb_eRuntimeError, "shape 0 of a must be %d", n1+n2);
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);
  {
    int shape[1];
    shape[0] = n1+n2;
    rb_index = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  index = NA_PTR_TYPE(rb_index, integer*);

  slamrg_(&n1, &n2, a, &strd1, &strd2, index);

  return rb_index;
}

void
init_lapack_slamrg(VALUE mLapack){
  rb_define_module_function(mLapack, "slamrg", rb_slamrg, -1);
}
