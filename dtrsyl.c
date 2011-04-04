#include "rb_lapack.h"

extern VOID dtrsyl_(char *trana, char *tranb, integer *isgn, integer *m, integer *n, doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *c, integer *ldc, doublereal *scale, integer *info);

static VALUE
rb_dtrsyl(int argc, VALUE *argv, VALUE self){
  VALUE rb_trana;
  char trana; 
  VALUE rb_tranb;
  char tranb; 
  VALUE rb_isgn;
  integer isgn; 
  VALUE rb_a;
  doublereal *a; 
  VALUE rb_b;
  doublereal *b; 
  VALUE rb_c;
  doublereal *c; 
  VALUE rb_scale;
  doublereal scale; 
  VALUE rb_info;
  integer info; 
  VALUE rb_c_out__;
  doublereal *c_out__;

  integer lda;
  integer m;
  integer ldb;
  integer n;
  integer ldc;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  scale, info, c = NumRu::Lapack.dtrsyl( trana, tranb, isgn, a, b, c)\n    or\n  NumRu::Lapack.dtrsyl  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rb_trana = argv[0];
  rb_tranb = argv[1];
  rb_isgn = argv[2];
  rb_a = argv[3];
  rb_b = argv[4];
  rb_c = argv[5];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (4th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (4th argument) must be %d", 2);
  m = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_DFLOAT)
    rb_a = na_change_type(rb_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rb_a, doublereal*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (5th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (5th argument) must be %d", 2);
  n = NA_SHAPE1(rb_b);
  ldb = NA_SHAPE0(rb_b);
  if (NA_TYPE(rb_b) != NA_DFLOAT)
    rb_b = na_change_type(rb_b, NA_DFLOAT);
  b = NA_PTR_TYPE(rb_b, doublereal*);
  trana = StringValueCStr(rb_trana)[0];
  if (!NA_IsNArray(rb_c))
    rb_raise(rb_eArgError, "c (6th argument) must be NArray");
  if (NA_RANK(rb_c) != 2)
    rb_raise(rb_eArgError, "rank of c (6th argument) must be %d", 2);
  if (NA_SHAPE1(rb_c) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of c must be the same as shape 1 of b");
  ldc = NA_SHAPE0(rb_c);
  if (NA_TYPE(rb_c) != NA_DFLOAT)
    rb_c = na_change_type(rb_c, NA_DFLOAT);
  c = NA_PTR_TYPE(rb_c, doublereal*);
  tranb = StringValueCStr(rb_tranb)[0];
  isgn = NUM2INT(rb_isgn);
  {
    int shape[2];
    shape[0] = ldc;
    shape[1] = n;
    rb_c_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  c_out__ = NA_PTR_TYPE(rb_c_out__, doublereal*);
  MEMCPY(c_out__, c, doublereal, NA_TOTAL(rb_c));
  rb_c = rb_c_out__;
  c = c_out__;

  dtrsyl_(&trana, &tranb, &isgn, &m, &n, a, &lda, b, &ldb, c, &ldc, &scale, &info);

  rb_scale = rb_float_new((double)scale);
  rb_info = INT2NUM(info);
  return rb_ary_new3(3, rb_scale, rb_info, rb_c);
}

void
init_lapack_dtrsyl(VALUE mLapack){
  rb_define_module_function(mLapack, "dtrsyl", rb_dtrsyl, -1);
}
