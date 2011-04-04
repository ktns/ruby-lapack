#include "rb_lapack.h"

extern doublereal zlantb_(char *norm, char *uplo, char *diag, integer *n, integer *k, doublecomplex *ab, integer *ldab, doublereal *work);

static VALUE
rb_zlantb(int argc, VALUE *argv, VALUE self){
  VALUE rb_norm;
  char norm; 
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_diag;
  char diag; 
  VALUE rb_k;
  integer k; 
  VALUE rb_ab;
  doublecomplex *ab; 
  VALUE rb___out__;
  doublereal __out__; 
  doublereal *work;

  integer ldab;
  integer n;
  integer lwork;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.zlantb( norm, uplo, diag, k, ab)\n    or\n  NumRu::Lapack.zlantb  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_norm = argv[0];
  rb_uplo = argv[1];
  rb_diag = argv[2];
  rb_k = argv[3];
  rb_ab = argv[4];

  k = NUM2INT(rb_k);
  if (!NA_IsNArray(rb_ab))
    rb_raise(rb_eArgError, "ab (5th argument) must be NArray");
  if (NA_RANK(rb_ab) != 2)
    rb_raise(rb_eArgError, "rank of ab (5th argument) must be %d", 2);
  n = NA_SHAPE1(rb_ab);
  ldab = NA_SHAPE0(rb_ab);
  if (NA_TYPE(rb_ab) != NA_DCOMPLEX)
    rb_ab = na_change_type(rb_ab, NA_DCOMPLEX);
  ab = NA_PTR_TYPE(rb_ab, doublecomplex*);
  diag = StringValueCStr(rb_diag)[0];
  norm = StringValueCStr(rb_norm)[0];
  uplo = StringValueCStr(rb_uplo)[0];
  lwork = lsame_(&norm,"I") ? n : 0;
  work = ALLOC_N(doublereal, (MAX(1,lwork)));

  __out__ = zlantb_(&norm, &uplo, &diag, &n, &k, ab, &ldab, work);

  free(work);
  rb___out__ = rb_float_new((double)__out__);
  return rb___out__;
}

void
init_lapack_zlantb(VALUE mLapack){
  rb_define_module_function(mLapack, "zlantb", rb_zlantb, -1);
}
