#include "rb_lapack.h"

extern VOID sgebal_(char *job, integer *n, real *a, integer *lda, integer *ilo, integer *ihi, real *scale, integer *info);

static VALUE
rb_sgebal(int argc, VALUE *argv, VALUE self){
  VALUE rb_job;
  char job; 
  VALUE rb_a;
  real *a; 
  VALUE rb_ilo;
  integer ilo; 
  VALUE rb_ihi;
  integer ihi; 
  VALUE rb_scale;
  real *scale; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  real *a_out__;

  integer lda;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  ilo, ihi, scale, info, a = NumRu::Lapack.sgebal( job, a)\n    or\n  NumRu::Lapack.sgebal  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rb_job = argv[0];
  rb_a = argv[1];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);
  job = StringValueCStr(rb_job)[0];
  {
    int shape[1];
    shape[0] = n;
    rb_scale = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  scale = NA_PTR_TYPE(rb_scale, real*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rb_a_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rb_a_out__, real*);
  MEMCPY(a_out__, a, real, NA_TOTAL(rb_a));
  rb_a = rb_a_out__;
  a = a_out__;

  sgebal_(&job, &n, a, &lda, &ilo, &ihi, scale, &info);

  rb_ilo = INT2NUM(ilo);
  rb_ihi = INT2NUM(ihi);
  rb_info = INT2NUM(info);
  return rb_ary_new3(5, rb_ilo, rb_ihi, rb_scale, rb_info, rb_a);
}

void
init_lapack_sgebal(VALUE mLapack){
  rb_define_module_function(mLapack, "sgebal", rb_sgebal, -1);
}
