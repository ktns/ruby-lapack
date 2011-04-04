#include "rb_lapack.h"

extern VOID dlaqtr_(logical *ltran, logical *lreal, integer *n, doublereal *t, integer *ldt, doublereal *b, doublereal *w, doublereal *scale, doublereal *x, doublereal *work, integer *info);

static VALUE
rb_dlaqtr(int argc, VALUE *argv, VALUE self){
  VALUE rb_ltran;
  logical ltran; 
  VALUE rb_lreal;
  logical lreal; 
  VALUE rb_t;
  doublereal *t; 
  VALUE rb_b;
  doublereal *b; 
  VALUE rb_w;
  doublereal w; 
  VALUE rb_x;
  doublereal *x; 
  VALUE rb_scale;
  doublereal scale; 
  VALUE rb_info;
  integer info; 
  VALUE rb_x_out__;
  doublereal *x_out__;
  doublereal *work;

  integer ldt;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  scale, info, x = NumRu::Lapack.dlaqtr( ltran, lreal, t, b, w, x)\n    or\n  NumRu::Lapack.dlaqtr  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rb_ltran = argv[0];
  rb_lreal = argv[1];
  rb_t = argv[2];
  rb_b = argv[3];
  rb_w = argv[4];
  rb_x = argv[5];

  w = NUM2DBL(rb_w);
  lreal = (rb_lreal == Qtrue);
  if (!NA_IsNArray(rb_t))
    rb_raise(rb_eArgError, "t (3th argument) must be NArray");
  if (NA_RANK(rb_t) != 2)
    rb_raise(rb_eArgError, "rank of t (3th argument) must be %d", 2);
  n = NA_SHAPE1(rb_t);
  ldt = NA_SHAPE0(rb_t);
  if (NA_TYPE(rb_t) != NA_DFLOAT)
    rb_t = na_change_type(rb_t, NA_DFLOAT);
  t = NA_PTR_TYPE(rb_t, doublereal*);
  ltran = (rb_ltran == Qtrue);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (4th argument) must be NArray");
  if (NA_RANK(rb_b) != 1)
    rb_raise(rb_eArgError, "rank of b (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_b) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of b must be the same as shape 1 of t");
  if (NA_TYPE(rb_b) != NA_DFLOAT)
    rb_b = na_change_type(rb_b, NA_DFLOAT);
  b = NA_PTR_TYPE(rb_b, doublereal*);
  if (!NA_IsNArray(rb_x))
    rb_raise(rb_eArgError, "x (6th argument) must be NArray");
  if (NA_RANK(rb_x) != 1)
    rb_raise(rb_eArgError, "rank of x (6th argument) must be %d", 1);
  if (NA_SHAPE0(rb_x) != (2*n))
    rb_raise(rb_eRuntimeError, "shape 0 of x must be %d", 2*n);
  if (NA_TYPE(rb_x) != NA_DFLOAT)
    rb_x = na_change_type(rb_x, NA_DFLOAT);
  x = NA_PTR_TYPE(rb_x, doublereal*);
  {
    int shape[1];
    shape[0] = 2*n;
    rb_x_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  x_out__ = NA_PTR_TYPE(rb_x_out__, doublereal*);
  MEMCPY(x_out__, x, doublereal, NA_TOTAL(rb_x));
  rb_x = rb_x_out__;
  x = x_out__;
  work = ALLOC_N(doublereal, (n));

  dlaqtr_(&ltran, &lreal, &n, t, &ldt, b, &w, &scale, x, work, &info);

  free(work);
  rb_scale = rb_float_new((double)scale);
  rb_info = INT2NUM(info);
  return rb_ary_new3(3, rb_scale, rb_info, rb_x);
}

void
init_lapack_dlaqtr(VALUE mLapack){
  rb_define_module_function(mLapack, "dlaqtr", rb_dlaqtr, -1);
}
