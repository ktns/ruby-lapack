#include "rb_lapack.h"

extern VOID dlaqr1_(integer *n, doublereal *h, integer *ldh, doublereal *sr1, doublereal *si1, doublereal *sr2, doublereal *si2, doublereal *v);

static VALUE
rb_dlaqr1(int argc, VALUE *argv, VALUE self){
  VALUE rb_h;
  doublereal *h; 
  VALUE rb_sr1;
  doublereal sr1; 
  VALUE rb_si1;
  doublereal si1; 
  VALUE rb_sr2;
  doublereal sr2; 
  VALUE rb_si2;
  doublereal si2; 
  VALUE rb_v;
  doublereal *v; 

  integer ldh;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  v = NumRu::Lapack.dlaqr1( h, sr1, si1, sr2, si2)\n    or\n  NumRu::Lapack.dlaqr1  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_h = argv[0];
  rb_sr1 = argv[1];
  rb_si1 = argv[2];
  rb_sr2 = argv[3];
  rb_si2 = argv[4];

  si1 = NUM2DBL(rb_si1);
  si2 = NUM2DBL(rb_si2);
  sr1 = NUM2DBL(rb_sr1);
  sr2 = NUM2DBL(rb_sr2);
  if (!NA_IsNArray(rb_h))
    rb_raise(rb_eArgError, "h (1th argument) must be NArray");
  if (NA_RANK(rb_h) != 2)
    rb_raise(rb_eArgError, "rank of h (1th argument) must be %d", 2);
  n = NA_SHAPE1(rb_h);
  ldh = NA_SHAPE0(rb_h);
  if (NA_TYPE(rb_h) != NA_DFLOAT)
    rb_h = na_change_type(rb_h, NA_DFLOAT);
  h = NA_PTR_TYPE(rb_h, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_v = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  v = NA_PTR_TYPE(rb_v, doublereal*);

  dlaqr1_(&n, h, &ldh, &sr1, &si1, &sr2, &si2, v);

  return rb_v;
}

void
init_lapack_dlaqr1(VALUE mLapack){
  rb_define_module_function(mLapack, "dlaqr1", rb_dlaqr1, -1);
}
