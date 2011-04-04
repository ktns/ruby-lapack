#include "rb_lapack.h"

extern VOID stpttr_(char *uplo, integer *n, real *ap, real *a, integer *lda, integer *info);

static VALUE
rb_stpttr(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_ap;
  real *ap; 
  VALUE rb_a;
  real *a; 
  VALUE rb_info;
  integer info; 

  integer ldap;
  integer lda;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  a, info = NumRu::Lapack.stpttr( uplo, ap)\n    or\n  NumRu::Lapack.stpttr  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rb_uplo = argv[0];
  rb_ap = argv[1];

  uplo = StringValueCStr(rb_uplo)[0];
  if (!NA_IsNArray(rb_ap))
    rb_raise(rb_eArgError, "ap (2th argument) must be NArray");
  if (NA_RANK(rb_ap) != 1)
    rb_raise(rb_eArgError, "rank of ap (2th argument) must be %d", 1);
  ldap = NA_SHAPE0(rb_ap);
  if (NA_TYPE(rb_ap) != NA_SFLOAT)
    rb_ap = na_change_type(rb_ap, NA_SFLOAT);
  ap = NA_PTR_TYPE(rb_ap, real*);
  n = ((int)sqrtf(ldap*8+1.0f)-1)/2;
  lda = MAX(1,n);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rb_a = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  a = NA_PTR_TYPE(rb_a, real*);

  stpttr_(&uplo, &n, ap, a, &lda, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_a, rb_info);
}

void
init_lapack_stpttr(VALUE mLapack){
  rb_define_module_function(mLapack, "stpttr", rb_stpttr, -1);
}
