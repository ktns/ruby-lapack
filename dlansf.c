#include "rb_lapack.h"

extern doublereal dlansf_(char *norm, char *transr, char *uplo, integer *n, doublereal *a, doublereal *work);

static VALUE
rb_dlansf(int argc, VALUE *argv, VALUE self){
  VALUE rb_norm;
  char norm; 
  VALUE rb_transr;
  char transr; 
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_n;
  integer n; 
  VALUE rb_a;
  doublereal *a; 
  VALUE rb___out__;
  doublereal __out__; 
  doublereal *work;

  integer lwork;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.dlansf( norm, transr, uplo, n, a)\n    or\n  NumRu::Lapack.dlansf  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_norm = argv[0];
  rb_transr = argv[1];
  rb_uplo = argv[2];
  rb_n = argv[3];
  rb_a = argv[4];

  transr = StringValueCStr(rb_transr)[0];
  n = NUM2INT(rb_n);
  uplo = StringValueCStr(rb_uplo)[0];
  norm = StringValueCStr(rb_norm)[0];
  lwork = ((lsame_(&norm,"I")) || ((('1') || ('o')))) ? n : 0;
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (5th argument) must be NArray");
  if (NA_RANK(rb_a) != 1)
    rb_raise(rb_eArgError, "rank of a (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_a) != (n*(n+1)/2))
    rb_raise(rb_eRuntimeError, "shape 0 of a must be %d", n*(n+1)/2);
  if (NA_TYPE(rb_a) != NA_DFLOAT)
    rb_a = na_change_type(rb_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rb_a, doublereal*);
  work = ALLOC_N(doublereal, (MAX(1,lwork)));

  __out__ = dlansf_(&norm, &transr, &uplo, &n, a, work);

  free(work);
  rb___out__ = rb_float_new((double)__out__);
  return rb___out__;
}

void
init_lapack_dlansf(VALUE mLapack){
  rb_define_module_function(mLapack, "dlansf", rb_dlansf, -1);
}
