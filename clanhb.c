#include "rb_lapack.h"

extern real clanhb_(char *norm, char *uplo, integer *n, integer *k, complex *ab, integer *ldab, real *work);

static VALUE
rb_clanhb(int argc, VALUE *argv, VALUE self){
  VALUE rb_norm;
  char norm; 
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_k;
  integer k; 
  VALUE rb_ab;
  complex *ab; 
  VALUE rb___out__;
  real __out__; 
  real *work;

  integer ldab;
  integer n;
  integer lwork;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.clanhb( norm, uplo, k, ab)\n    or\n  NumRu::Lapack.clanhb  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_norm = argv[0];
  rb_uplo = argv[1];
  rb_k = argv[2];
  rb_ab = argv[3];

  k = NUM2INT(rb_k);
  if (!NA_IsNArray(rb_ab))
    rb_raise(rb_eArgError, "ab (4th argument) must be NArray");
  if (NA_RANK(rb_ab) != 2)
    rb_raise(rb_eArgError, "rank of ab (4th argument) must be %d", 2);
  n = NA_SHAPE1(rb_ab);
  ldab = NA_SHAPE0(rb_ab);
  if (NA_TYPE(rb_ab) != NA_SCOMPLEX)
    rb_ab = na_change_type(rb_ab, NA_SCOMPLEX);
  ab = NA_PTR_TYPE(rb_ab, complex*);
  norm = StringValueCStr(rb_norm)[0];
  uplo = StringValueCStr(rb_uplo)[0];
  lwork = ((lsame_(&norm,"I")) || ((('1') || ('o')))) ? n : 0;
  work = ALLOC_N(real, (MAX(1,lwork)));

  __out__ = clanhb_(&norm, &uplo, &n, &k, ab, &ldab, work);

  free(work);
  rb___out__ = rb_float_new((double)__out__);
  return rb___out__;
}

void
init_lapack_clanhb(VALUE mLapack){
  rb_define_module_function(mLapack, "clanhb", rb_clanhb, -1);
}
