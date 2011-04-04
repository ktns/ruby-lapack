#include "rb_lapack.h"

extern VOID slaqr1_(integer *n, real *h, integer *ldh, real *sr1, real *si1, real *sr2, real *si2, real *v);

static VALUE
rb_slaqr1(int argc, VALUE *argv, VALUE self){
  VALUE rb_h;
  real *h; 
  VALUE rb_sr1;
  real sr1; 
  VALUE rb_si1;
  real si1; 
  VALUE rb_sr2;
  real sr2; 
  VALUE rb_si2;
  real si2; 
  VALUE rb_v;
  real *v; 

  integer ldh;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  v = NumRu::Lapack.slaqr1( h, sr1, si1, sr2, si2)\n    or\n  NumRu::Lapack.slaqr1  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_h = argv[0];
  rb_sr1 = argv[1];
  rb_si1 = argv[2];
  rb_sr2 = argv[3];
  rb_si2 = argv[4];

  si1 = (real)NUM2DBL(rb_si1);
  si2 = (real)NUM2DBL(rb_si2);
  sr1 = (real)NUM2DBL(rb_sr1);
  sr2 = (real)NUM2DBL(rb_sr2);
  if (!NA_IsNArray(rb_h))
    rb_raise(rb_eArgError, "h (1th argument) must be NArray");
  if (NA_RANK(rb_h) != 2)
    rb_raise(rb_eArgError, "rank of h (1th argument) must be %d", 2);
  n = NA_SHAPE1(rb_h);
  ldh = NA_SHAPE0(rb_h);
  if (NA_TYPE(rb_h) != NA_SFLOAT)
    rb_h = na_change_type(rb_h, NA_SFLOAT);
  h = NA_PTR_TYPE(rb_h, real*);
  {
    int shape[1];
    shape[0] = n;
    rb_v = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  v = NA_PTR_TYPE(rb_v, real*);

  slaqr1_(&n, h, &ldh, &sr1, &si1, &sr2, &si2, v);

  return rb_v;
}

void
init_lapack_slaqr1(VALUE mLapack){
  rb_define_module_function(mLapack, "slaqr1", rb_slaqr1, -1);
}
