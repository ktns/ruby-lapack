#include "rb_lapack.h"

extern VOID ztpttf_(char *transr, char *uplo, integer *n, doublecomplex *ap, doublecomplex *arf, integer *info);

static VALUE
rb_ztpttf(int argc, VALUE *argv, VALUE self){
  VALUE rb_transr;
  char transr; 
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_n;
  integer n; 
  VALUE rb_ap;
  doublecomplex *ap; 
  VALUE rb_arf;
  doublecomplex *arf; 
  VALUE rb_info;
  integer info; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  arf, info = NumRu::Lapack.ztpttf( transr, uplo, n, ap)\n    or\n  NumRu::Lapack.ztpttf  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_transr = argv[0];
  rb_uplo = argv[1];
  rb_n = argv[2];
  rb_ap = argv[3];

  n = NUM2INT(rb_n);
  transr = StringValueCStr(rb_transr)[0];
  uplo = StringValueCStr(rb_uplo)[0];
  if (!NA_IsNArray(rb_ap))
    rb_raise(rb_eArgError, "ap (4th argument) must be NArray");
  if (NA_RANK(rb_ap) != 1)
    rb_raise(rb_eArgError, "rank of ap (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_ap) != (( n*(n+1)/2 )))
    rb_raise(rb_eRuntimeError, "shape 0 of ap must be %d", ( n*(n+1)/2 ));
  if (NA_TYPE(rb_ap) != NA_DCOMPLEX)
    rb_ap = na_change_type(rb_ap, NA_DCOMPLEX);
  ap = NA_PTR_TYPE(rb_ap, doublecomplex*);
  {
    int shape[1];
    shape[0] = ( n*(n+1)/2 );
    rb_arf = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  arf = NA_PTR_TYPE(rb_arf, doublecomplex*);

  ztpttf_(&transr, &uplo, &n, ap, arf, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_arf, rb_info);
}

void
init_lapack_ztpttf(VALUE mLapack){
  rb_define_module_function(mLapack, "ztpttf", rb_ztpttf, -1);
}