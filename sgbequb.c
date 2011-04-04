#include "rb_lapack.h"

extern VOID sgbequb_(integer *m, integer *n, integer *kl, integer *ku, doublereal *ab, integer *ldab, real *r, real *c, real *rowcnd, real *colcnd, real *amax, integer *info);

static VALUE
rb_sgbequb(int argc, VALUE *argv, VALUE self){
  VALUE rb_kl;
  integer kl; 
  VALUE rb_ku;
  integer ku; 
  VALUE rb_ab;
  doublereal *ab; 
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

  integer ldab;
  integer n;
  integer m;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  r, c, rowcnd, colcnd, amax, info = NumRu::Lapack.sgbequb( kl, ku, ab)\n    or\n  NumRu::Lapack.sgbequb  # print help\n\n\nFORTRAN MANUAL\n\n");
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
    rb_r = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  r = NA_PTR_TYPE(rb_r, real*);
  {
    int shape[1];
    shape[0] = n;
    rb_c = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  c = NA_PTR_TYPE(rb_c, real*);

  sgbequb_(&m, &n, &kl, &ku, ab, &ldab, r, c, &rowcnd, &colcnd, &amax, &info);

  rb_rowcnd = rb_float_new((double)rowcnd);
  rb_colcnd = rb_float_new((double)colcnd);
  rb_amax = rb_float_new((double)amax);
  rb_info = INT2NUM(info);
  return rb_ary_new3(6, rb_r, rb_c, rb_rowcnd, rb_colcnd, rb_amax, rb_info);
}

void
init_lapack_sgbequb(VALUE mLapack){
  rb_define_module_function(mLapack, "sgbequb", rb_sgbequb, -1);
}
