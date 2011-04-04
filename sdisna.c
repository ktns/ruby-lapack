#include "rb_lapack.h"

extern VOID sdisna_(char *job, integer *m, integer *n, real *d, real *sep, integer *info);

static VALUE
rb_sdisna(int argc, VALUE *argv, VALUE self){
  VALUE rb_job;
  char job; 
  VALUE rb_n;
  integer n; 
  VALUE rb_d;
  real *d; 
  VALUE rb_sep;
  real *sep; 
  VALUE rb_info;
  integer info; 

  integer m;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  sep, info = NumRu::Lapack.sdisna( job, n, d)\n    or\n  NumRu::Lapack.sdisna  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_job = argv[0];
  rb_n = argv[1];
  rb_d = argv[2];

  n = NUM2INT(rb_n);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (3th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (3th argument) must be %d", 1);
  m = NA_SHAPE0(rb_d);
  if (NA_TYPE(rb_d) != NA_SFLOAT)
    rb_d = na_change_type(rb_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rb_d, real*);
  job = StringValueCStr(rb_job)[0];
  {
    int shape[1];
    shape[0] = lsame_(&job,"E") ? m : ((lsame_(&job,"L")) || (lsame_(&job,"R"))) ? MIN(m,n) : 0;
    rb_sep = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  sep = NA_PTR_TYPE(rb_sep, real*);

  sdisna_(&job, &m, &n, d, sep, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_sep, rb_info);
}

void
init_lapack_sdisna(VALUE mLapack){
  rb_define_module_function(mLapack, "sdisna", rb_sdisna, -1);
}
