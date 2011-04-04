#include "rb_lapack.h"

extern VOID dsyequb_(char *uplo, integer *n, doublereal *a, integer *lda, doublereal *s, doublereal *scond, doublereal *amax, doublereal *work, integer *info);

static VALUE
rb_dsyequb(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_a;
  doublereal *a; 
  VALUE rb_s;
  doublereal *s; 
  VALUE rb_scond;
  doublereal scond; 
  VALUE rb_amax;
  doublereal amax; 
  VALUE rb_info;
  integer info; 
  doublereal *work;

  integer lda;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  s, scond, amax, info = NumRu::Lapack.dsyequb( uplo, a)\n    or\n  NumRu::Lapack.dsyequb  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rb_uplo = argv[0];
  rb_a = argv[1];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_DFLOAT)
    rb_a = na_change_type(rb_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rb_a, doublereal*);
  uplo = StringValueCStr(rb_uplo)[0];
  {
    int shape[1];
    shape[0] = n;
    rb_s = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  s = NA_PTR_TYPE(rb_s, doublereal*);
  work = ALLOC_N(doublereal, (3*n));

  dsyequb_(&uplo, &n, a, &lda, s, &scond, &amax, work, &info);

  free(work);
  rb_scond = rb_float_new((double)scond);
  rb_amax = rb_float_new((double)amax);
  rb_info = INT2NUM(info);
  return rb_ary_new3(4, rb_s, rb_scond, rb_amax, rb_info);
}

void
init_lapack_dsyequb(VALUE mLapack){
  rb_define_module_function(mLapack, "dsyequb", rb_dsyequb, -1);
}
