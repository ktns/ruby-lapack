#include "rb_lapack.h"

extern VOID zlat2c_(char *uplo, integer *n, doublecomplex *a, integer *lda, complex *sa, integer *ldsa, integer *info);

static VALUE
rb_zlat2c(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_a;
  doublecomplex *a; 
  VALUE rb_sa;
  complex *sa; 
  VALUE rb_info;
  integer info; 

  integer lda;
  integer n;
  integer ldsa;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  sa, info = NumRu::Lapack.zlat2c( uplo, a)\n    or\n  NumRu::Lapack.zlat2c  # print help\n\n\nFORTRAN MANUAL\n\n");
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
  if (NA_TYPE(rb_a) != NA_DCOMPLEX)
    rb_a = na_change_type(rb_a, NA_DCOMPLEX);
  a = NA_PTR_TYPE(rb_a, doublecomplex*);
  uplo = StringValueCStr(rb_uplo)[0];
  ldsa = MAX(1,n);
  {
    int shape[2];
    shape[0] = ldsa;
    shape[1] = n;
    rb_sa = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  sa = NA_PTR_TYPE(rb_sa, complex*);

  zlat2c_(&uplo, &n, a, &lda, sa, &ldsa, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_sa, rb_info);
}

void
init_lapack_zlat2c(VALUE mLapack){
  rb_define_module_function(mLapack, "zlat2c", rb_zlat2c, -1);
}
