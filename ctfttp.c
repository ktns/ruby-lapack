#include "rb_lapack.h"

extern VOID ctfttp_(char *transr, char *uplo, integer *n, complex *arf, complex *ap, integer *info);

static VALUE
rb_ctfttp(int argc, VALUE *argv, VALUE self){
  VALUE rb_transr;
  char transr; 
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_n;
  integer n; 
  VALUE rb_arf;
  complex *arf; 
  VALUE rb_ap;
  complex *ap; 
  VALUE rb_info;
  integer info; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  ap, info = NumRu::Lapack.ctfttp( transr, uplo, n, arf)\n    or\n  NumRu::Lapack.ctfttp  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_transr = argv[0];
  rb_uplo = argv[1];
  rb_n = argv[2];
  rb_arf = argv[3];

  n = NUM2INT(rb_n);
  transr = StringValueCStr(rb_transr)[0];
  uplo = StringValueCStr(rb_uplo)[0];
  if (!NA_IsNArray(rb_arf))
    rb_raise(rb_eArgError, "arf (4th argument) must be NArray");
  if (NA_RANK(rb_arf) != 1)
    rb_raise(rb_eArgError, "rank of arf (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_arf) != (( n*(n+1)/2 )))
    rb_raise(rb_eRuntimeError, "shape 0 of arf must be %d", ( n*(n+1)/2 ));
  if (NA_TYPE(rb_arf) != NA_SCOMPLEX)
    rb_arf = na_change_type(rb_arf, NA_SCOMPLEX);
  arf = NA_PTR_TYPE(rb_arf, complex*);
  {
    int shape[1];
    shape[0] = ( n*(n+1)/2 );
    rb_ap = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  ap = NA_PTR_TYPE(rb_ap, complex*);

  ctfttp_(&transr, &uplo, &n, arf, ap, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_ap, rb_info);
}

void
init_lapack_ctfttp(VALUE mLapack){
  rb_define_module_function(mLapack, "ctfttp", rb_ctfttp, -1);
}
