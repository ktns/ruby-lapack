#include "rb_lapack.h"

extern real clanhs_(char *norm, integer *n, complex *a, integer *lda, real *work);

static VALUE
rb_clanhs(int argc, VALUE *argv, VALUE self){
  VALUE rb_norm;
  char norm; 
  VALUE rb_a;
  complex *a; 
  VALUE rb___out__;
  real __out__; 
  real *work;

  integer lda;
  integer n;
  integer lwork;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.clanhs( norm, a)\n    or\n  NumRu::Lapack.clanhs  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rb_norm = argv[0];
  rb_a = argv[1];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SCOMPLEX)
    rb_a = na_change_type(rb_a, NA_SCOMPLEX);
  a = NA_PTR_TYPE(rb_a, complex*);
  norm = StringValueCStr(rb_norm)[0];
  lwork = lsame_(&norm,"I") ? n : 0;
  work = ALLOC_N(real, (MAX(1,lwork)));

  __out__ = clanhs_(&norm, &n, a, &lda, work);

  free(work);
  rb___out__ = rb_float_new((double)__out__);
  return rb___out__;
}

void
init_lapack_clanhs(VALUE mLapack){
  rb_define_module_function(mLapack, "clanhs", rb_clanhs, -1);
}
