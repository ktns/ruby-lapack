#include "rb_lapack.h"

extern VOID ctfttr_(char *transr, char *uplo, integer *n, complex *arf, complex *a, integer *lda, integer *info);

static VALUE
rb_ctfttr(int argc, VALUE *argv, VALUE self){
  VALUE rb_transr;
  char transr; 
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_arf;
  complex *arf; 
  VALUE rb_a;
  complex *a; 
  VALUE rb_info;
  integer info; 

  integer ldarf;
  integer lda;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  a, info = NumRu::Lapack.ctfttr( transr, uplo, arf)\n    or\n  NumRu::Lapack.ctfttr  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_transr = argv[0];
  rb_uplo = argv[1];
  rb_arf = argv[2];

  transr = StringValueCStr(rb_transr)[0];
  uplo = StringValueCStr(rb_uplo)[0];
  if (!NA_IsNArray(rb_arf))
    rb_raise(rb_eArgError, "arf (3th argument) must be NArray");
  if (NA_RANK(rb_arf) != 1)
    rb_raise(rb_eArgError, "rank of arf (3th argument) must be %d", 1);
  ldarf = NA_SHAPE0(rb_arf);
  if (NA_TYPE(rb_arf) != NA_SCOMPLEX)
    rb_arf = na_change_type(rb_arf, NA_SCOMPLEX);
  arf = NA_PTR_TYPE(rb_arf, complex*);
  n = ((int)sqrtf(ldarf*8+1.0f)-1)/2;
  lda = MAX(1,n);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rb_a = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  a = NA_PTR_TYPE(rb_a, complex*);

  ctfttr_(&transr, &uplo, &n, arf, a, &lda, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_a, rb_info);
}

void
init_lapack_ctfttr(VALUE mLapack){
  rb_define_module_function(mLapack, "ctfttr", rb_ctfttr, -1);
}
