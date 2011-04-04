#include "rb_lapack.h"

extern VOID zlaqgb_(integer *m, integer *n, integer *kl, integer *ku, doublecomplex *ab, integer *ldab, doublereal *r, doublereal *c, doublereal *rowcnd, doublereal *colcnd, doublereal *amax, char *equed);

static VALUE
rb_zlaqgb(int argc, VALUE *argv, VALUE self){
  VALUE rb_kl;
  integer kl; 
  VALUE rb_ku;
  integer ku; 
  VALUE rb_ab;
  doublecomplex *ab; 
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
  VALUE rb_equed;
  char equed; 
  VALUE rb_ab_out__;
  doublecomplex *ab_out__;

  integer ldab;
  integer n;
  integer m;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  equed, ab = NumRu::Lapack.zlaqgb( kl, ku, ab, r, c, rowcnd, colcnd, amax)\n    or\n  NumRu::Lapack.zlaqgb  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 8)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 8)", argc);
  rb_kl = argv[0];
  rb_ku = argv[1];
  rb_ab = argv[2];
  rb_r = argv[3];
  rb_c = argv[4];
  rb_rowcnd = argv[5];
  rb_colcnd = argv[6];
  rb_amax = argv[7];

  if (!NA_IsNArray(rb_ab))
    rb_raise(rb_eArgError, "ab (3th argument) must be NArray");
  if (NA_RANK(rb_ab) != 2)
    rb_raise(rb_eArgError, "rank of ab (3th argument) must be %d", 2);
  n = NA_SHAPE1(rb_ab);
  ldab = NA_SHAPE0(rb_ab);
  if (NA_TYPE(rb_ab) != NA_DCOMPLEX)
    rb_ab = na_change_type(rb_ab, NA_DCOMPLEX);
  ab = NA_PTR_TYPE(rb_ab, doublecomplex*);
  kl = NUM2INT(rb_kl);
  if (!NA_IsNArray(rb_c))
    rb_raise(rb_eArgError, "c (5th argument) must be NArray");
  if (NA_RANK(rb_c) != 1)
    rb_raise(rb_eArgError, "rank of c (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_c) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of c must be the same as shape 1 of ab");
  if (NA_TYPE(rb_c) != NA_DFLOAT)
    rb_c = na_change_type(rb_c, NA_DFLOAT);
  c = NA_PTR_TYPE(rb_c, doublereal*);
  amax = NUM2DBL(rb_amax);
  colcnd = NUM2DBL(rb_colcnd);
  if (!NA_IsNArray(rb_r))
    rb_raise(rb_eArgError, "r (4th argument) must be NArray");
  if (NA_RANK(rb_r) != 1)
    rb_raise(rb_eArgError, "rank of r (4th argument) must be %d", 1);
  m = NA_SHAPE0(rb_r);
  if (NA_TYPE(rb_r) != NA_DFLOAT)
    rb_r = na_change_type(rb_r, NA_DFLOAT);
  r = NA_PTR_TYPE(rb_r, doublereal*);
  rowcnd = NUM2DBL(rb_rowcnd);
  ku = NUM2INT(rb_ku);
  {
    int shape[2];
    shape[0] = ldab;
    shape[1] = n;
    rb_ab_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  ab_out__ = NA_PTR_TYPE(rb_ab_out__, doublecomplex*);
  MEMCPY(ab_out__, ab, doublecomplex, NA_TOTAL(rb_ab));
  rb_ab = rb_ab_out__;
  ab = ab_out__;

  zlaqgb_(&m, &n, &kl, &ku, ab, &ldab, r, c, &rowcnd, &colcnd, &amax, &equed);

  rb_equed = rb_str_new(&equed,1);
  return rb_ary_new3(2, rb_equed, rb_ab);
}

void
init_lapack_zlaqgb(VALUE mLapack){
  rb_define_module_function(mLapack, "zlaqgb", rb_zlaqgb, -1);
}
