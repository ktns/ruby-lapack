#include "rb_lapack.h"

extern VOID dtfttr_(char *transr, char *uplo, integer *n, doublereal *arf, doublereal *a, integer *lda, integer *info);

static VALUE
rb_dtfttr(int argc, VALUE *argv, VALUE self){
  VALUE rb_transr;
  char transr; 
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_arf;
  doublereal *arf; 
  VALUE rb_a;
  doublereal *a; 
  VALUE rb_info;
  integer info; 
  integer ldarf;

  integer lda;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  a, info = NumRu::Lapack.dtfttr( transr, uplo, arf)\n    or\n  NumRu::Lapack.dtfttr  # print help\n\n\nFORTRAN MANUAL\n\n");
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
  if (NA_TYPE(rb_arf) != NA_DFLOAT)
    rb_arf = na_change_type(rb_arf, NA_DFLOAT);
  arf = NA_PTR_TYPE(rb_arf, doublereal*);
  n = ((int)sqrtf(8*ldarf+1.0f)-1)/2;
  lda = MAX(1,n);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rb_a = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  a = NA_PTR_TYPE(rb_a, doublereal*);

  dtfttr_(&transr, &uplo, &n, arf, a, &lda, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_a, rb_info);
}

void
init_lapack_dtfttr(VALUE mLapack){
  rb_define_module_function(mLapack, "dtfttr", rb_dtfttr, -1);
}
