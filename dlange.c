#include "rb_lapack.h"

extern doublereal dlange_(char *norm, integer *m, integer *n, doublereal *a, integer *lda, doublereal *work);

static VALUE
rb_dlange(int argc, VALUE *argv, VALUE self){
  VALUE rb_norm;
  char norm; 
  VALUE rb_m;
  integer m; 
  VALUE rb_a;
  doublereal *a; 
  VALUE rb___out__;
  doublereal __out__; 
  doublereal *work;

  integer lda;
  integer n;
  integer lwork;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.dlange( norm, m, a)\n    or\n  NumRu::Lapack.dlange  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_norm = argv[0];
  rb_m = argv[1];
  rb_a = argv[2];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_DFLOAT)
    rb_a = na_change_type(rb_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rb_a, doublereal*);
  m = NUM2INT(rb_m);
  norm = StringValueCStr(rb_norm)[0];
  lwork = lsame_(&norm,"I") ? m : 0;
  work = ALLOC_N(doublereal, (MAX(1,lwork)));

  __out__ = dlange_(&norm, &m, &n, a, &lda, work);

  free(work);
  rb___out__ = rb_float_new((double)__out__);
  return rb___out__;
}

void
init_lapack_dlange(VALUE mLapack){
  rb_define_module_function(mLapack, "dlange", rb_dlange, -1);
}
