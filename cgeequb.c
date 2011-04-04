#include "rb_lapack.h"

extern VOID cgeequb_(integer *m, integer *n, complex *a, integer *lda, real *r, real *c, real *rowcnd, real *colcnd, real *amax, integer *info);

static VALUE
rb_cgeequb(int argc, VALUE *argv, VALUE self){
  VALUE rb_a;
  complex *a; 
  VALUE rb_r;
  real *r; 
  VALUE rb_c;
  real *c; 
  VALUE rb_rowcnd;
  real rowcnd; 
  VALUE rb_colcnd;
  real colcnd; 
  VALUE rb_amax;
  real amax; 
  VALUE rb_info;
  integer info; 

  integer lda;
  integer n;
  integer m;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  r, c, rowcnd, colcnd, amax, info = NumRu::Lapack.cgeequb( a)\n    or\n  NumRu::Lapack.cgeequb  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 1)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 1)", argc);
  rb_a = argv[0];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (1th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (1th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (lda != (m))
    rb_raise(rb_eRuntimeError, "shape 0 of a must be %d", m);
  m = lda;
  if (NA_TYPE(rb_a) != NA_SCOMPLEX)
    rb_a = na_change_type(rb_a, NA_SCOMPLEX);
  a = NA_PTR_TYPE(rb_a, complex*);
  lda = m;
  {
    int shape[1];
    shape[0] = m;
    rb_r = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  r = NA_PTR_TYPE(rb_r, real*);
  {
    int shape[1];
    shape[0] = n;
    rb_c = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  c = NA_PTR_TYPE(rb_c, real*);

  cgeequb_(&m, &n, a, &lda, r, c, &rowcnd, &colcnd, &amax, &info);

  rb_rowcnd = rb_float_new((double)rowcnd);
  rb_colcnd = rb_float_new((double)colcnd);
  rb_amax = rb_float_new((double)amax);
  rb_info = INT2NUM(info);
  return rb_ary_new3(6, rb_r, rb_c, rb_rowcnd, rb_colcnd, rb_amax, rb_info);
}

void
init_lapack_cgeequb(VALUE mLapack){
  rb_define_module_function(mLapack, "cgeequb", rb_cgeequb, -1);
}
