#include "rb_lapack.h"

extern doublereal dlantr_(char *norm, char *uplo, char *diag, integer *m, integer *n, doublereal *a, integer *lda, doublereal *work);

static VALUE
rb_dlantr(int argc, VALUE *argv, VALUE self){
  VALUE rb_norm;
  char norm; 
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_diag;
  char diag; 
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
    printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.dlantr( norm, uplo, diag, m, a)\n    or\n  NumRu::Lapack.dlantr  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_norm = argv[0];
  rb_uplo = argv[1];
  rb_diag = argv[2];
  rb_m = argv[3];
  rb_a = argv[4];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (5th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (5th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_DFLOAT)
    rb_a = na_change_type(rb_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rb_a, doublereal*);
  m = NUM2INT(rb_m);
  diag = StringValueCStr(rb_diag)[0];
  norm = StringValueCStr(rb_norm)[0];
  uplo = StringValueCStr(rb_uplo)[0];
  lwork = lsame_(&norm,"I") ? m : 0;
  work = ALLOC_N(doublereal, (MAX(1,lwork)));

  __out__ = dlantr_(&norm, &uplo, &diag, &m, &n, a, &lda, work);

  free(work);
  rb___out__ = rb_float_new((double)__out__);
  return rb___out__;
}

void
init_lapack_dlantr(VALUE mLapack){
  rb_define_module_function(mLapack, "dlantr", rb_dlantr, -1);
}
