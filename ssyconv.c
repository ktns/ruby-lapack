#include "rb_lapack.h"

extern VOID ssyconv_(char *uplo, char *way, integer *n, real *a, integer *lda, integer *ipiv, real *work, integer *info);

static VALUE
rb_ssyconv(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_way;
  char way; 
  VALUE rb_a;
  real *a; 
  VALUE rb_ipiv;
  integer *ipiv; 
  VALUE rb_info;
  integer info; 
  real *work;

  integer lda;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info = NumRu::Lapack.ssyconv( uplo, way, a, ipiv)\n    or\n  NumRu::Lapack.ssyconv  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_uplo = argv[0];
  rb_way = argv[1];
  rb_a = argv[2];
  rb_ipiv = argv[3];

  if (!NA_IsNArray(rb_ipiv))
    rb_raise(rb_eArgError, "ipiv (4th argument) must be NArray");
  if (NA_RANK(rb_ipiv) != 1)
    rb_raise(rb_eArgError, "rank of ipiv (4th argument) must be %d", 1);
  n = NA_SHAPE0(rb_ipiv);
  if (NA_TYPE(rb_ipiv) != NA_LINT)
    rb_ipiv = na_change_type(rb_ipiv, NA_LINT);
  ipiv = NA_PTR_TYPE(rb_ipiv, integer*);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  if (NA_SHAPE1(rb_a) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of a must be the same as shape 0 of ipiv");
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);
  way = StringValueCStr(rb_way)[0];
  uplo = StringValueCStr(rb_uplo)[0];
  work = ALLOC_N(real, (MAX(1,n)));

  ssyconv_(&uplo, &way, &n, a, &lda, ipiv, work, &info);

  free(work);
  rb_info = INT2NUM(info);
  return rb_info;
}

void
init_lapack_ssyconv(VALUE mLapack){
  rb_define_module_function(mLapack, "ssyconv", rb_ssyconv, -1);
}
