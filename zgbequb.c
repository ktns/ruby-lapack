#include "rb_lapack.h"

extern VOID zgbequb_(integer *m, integer *n, integer *kl, integer *ku, doublereal *ab, integer *ldab, doublereal *r, doublereal *c, doublereal *rowcnd, doublereal *colcnd, doublereal *amax, integer *info);

static VALUE
rb_zgbequb(int argc, VALUE *argv, VALUE self){
  VALUE rb_kl;
  integer kl; 
  VALUE rb_ku;
  integer ku; 
  VALUE rb_ab;
  doublereal *ab; 
  VALUE rb_r;
  doublereal *r; 
  VALUE rb_c;
  doublereal *c; 
  VALUE rb_rowcnd;
  doublereal rowcnd; 
  VALUE rb_colcnd;
  doublereal colcnd; 
  VALUE rb_amax;
  doublereal amax; 
  VALUE rb_info;
  integer info; 

  integer ldab;
  integer n;
  integer m;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  r, c, rowcnd, colcnd, amax, info = NumRu::Lapack.zgbequb( kl, ku, ab)\n    or\n  NumRu::Lapack.zgbequb  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_kl = argv[0];
  rb_ku = argv[1];
  rb_ab = argv[2];

  if (!NA_IsNArray(rb_ab))
    rb_raise(rb_eArgError, "ab (3th argument) must be NArray");
  if (NA_RANK(rb_ab) != 2)
    rb_raise(rb_eArgError, "rank of ab (3th argument) must be %d", 2);
  n = NA_SHAPE1(rb_ab);
  ldab = NA_SHAPE0(rb_ab);
  if (ldab != (m))
    rb_raise(rb_eRuntimeError, "shape 0 of ab must be %d", m);
  m = ldab;
  if (NA_TYPE(rb_ab) != NA_DFLOAT)
    rb_ab = na_change_type(rb_ab, NA_DFLOAT);
  ab = NA_PTR_TYPE(rb_ab, doublereal*);
  kl = NUM2INT(rb_kl);
  ku = NUM2INT(rb_ku);
  ldab = m;
  {
    int shape[1];
    shape[0] = m;
    rb_r = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  r = NA_PTR_TYPE(rb_r, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_c = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  c = NA_PTR_TYPE(rb_c, doublereal*);

  zgbequb_(&m, &n, &kl, &ku, ab, &ldab, r, c, &rowcnd, &colcnd, &amax, &info);

  rb_rowcnd = rb_float_new((double)rowcnd);
  rb_colcnd = rb_float_new((double)colcnd);
  rb_amax = rb_float_new((double)amax);
  rb_info = INT2NUM(info);
  return rb_ary_new3(6, rb_r, rb_c, rb_rowcnd, rb_colcnd, rb_amax, rb_info);
}

void
init_lapack_zgbequb(VALUE mLapack){
  rb_define_module_function(mLapack, "zgbequb", rb_zgbequb, -1);
}
