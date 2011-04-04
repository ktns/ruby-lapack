#include "rb_lapack.h"

extern VOID clarnv_(integer *idist, integer *iseed, integer *n, complex *x);

static VALUE
rb_clarnv(int argc, VALUE *argv, VALUE self){
  VALUE rb_idist;
  integer idist; 
  VALUE rb_iseed;
  integer *iseed; 
  VALUE rb_n;
  integer n; 
  VALUE rb_x;
  complex *x; 
  VALUE rb_iseed_out__;
  integer *iseed_out__;


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  x, iseed = NumRu::Lapack.clarnv( idist, iseed, n)\n    or\n  NumRu::Lapack.clarnv  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_idist = argv[0];
  rb_iseed = argv[1];
  rb_n = argv[2];

  n = NUM2INT(rb_n);
  idist = NUM2INT(rb_idist);
  if (!NA_IsNArray(rb_iseed))
    rb_raise(rb_eArgError, "iseed (2th argument) must be NArray");
  if (NA_RANK(rb_iseed) != 1)
    rb_raise(rb_eArgError, "rank of iseed (2th argument) must be %d", 1);
  if (NA_SHAPE0(rb_iseed) != (4))
    rb_raise(rb_eRuntimeError, "shape 0 of iseed must be %d", 4);
  if (NA_TYPE(rb_iseed) != NA_LINT)
    rb_iseed = na_change_type(rb_iseed, NA_LINT);
  iseed = NA_PTR_TYPE(rb_iseed, integer*);
  {
    int shape[1];
    shape[0] = MAX(1,n);
    rb_x = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  x = NA_PTR_TYPE(rb_x, complex*);
  {
    int shape[1];
    shape[0] = 4;
    rb_iseed_out__ = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  iseed_out__ = NA_PTR_TYPE(rb_iseed_out__, integer*);
  MEMCPY(iseed_out__, iseed, integer, NA_TOTAL(rb_iseed));
  rb_iseed = rb_iseed_out__;
  iseed = iseed_out__;

  clarnv_(&idist, iseed, &n, x);

  return rb_ary_new3(2, rb_x, rb_iseed);
}

void
init_lapack_clarnv(VALUE mLapack){
  rb_define_module_function(mLapack, "clarnv", rb_clarnv, -1);
}
