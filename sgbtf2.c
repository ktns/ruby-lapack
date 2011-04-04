#include "rb_lapack.h"

extern VOID sgbtf2_(integer *m, integer *n, integer *kl, integer *ku, real *ab, integer *ldab, integer *ipiv, integer *info);

static VALUE
rb_sgbtf2(int argc, VALUE *argv, VALUE self){
  VALUE rb_m;
  integer m; 
  VALUE rb_kl;
  integer kl; 
  VALUE rb_ku;
  integer ku; 
  VALUE rb_ab;
  real *ab; 
  VALUE rb_ipiv;
  integer *ipiv; 
  VALUE rb_info;
  integer info; 
  VALUE rb_ab_out__;
  real *ab_out__;

  integer ldab;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  ipiv, info, ab = NumRu::Lapack.sgbtf2( m, kl, ku, ab)\n    or\n  NumRu::Lapack.sgbtf2  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_m = argv[0];
  rb_kl = argv[1];
  rb_ku = argv[2];
  rb_ab = argv[3];

  if (!NA_IsNArray(rb_ab))
    rb_raise(rb_eArgError, "ab (4th argument) must be NArray");
  if (NA_RANK(rb_ab) != 2)
    rb_raise(rb_eArgError, "rank of ab (4th argument) must be %d", 2);
  n = NA_SHAPE1(rb_ab);
  ldab = NA_SHAPE0(rb_ab);
  if (NA_TYPE(rb_ab) != NA_SFLOAT)
    rb_ab = na_change_type(rb_ab, NA_SFLOAT);
  ab = NA_PTR_TYPE(rb_ab, real*);
  kl = NUM2INT(rb_kl);
  m = NUM2INT(rb_m);
  ku = NUM2INT(rb_ku);
  {
    int shape[1];
    shape[0] = MIN(m,n);
    rb_ipiv = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  ipiv = NA_PTR_TYPE(rb_ipiv, integer*);
  {
    int shape[2];
    shape[0] = ldab;
    shape[1] = n;
    rb_ab_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  ab_out__ = NA_PTR_TYPE(rb_ab_out__, real*);
  MEMCPY(ab_out__, ab, real, NA_TOTAL(rb_ab));
  rb_ab = rb_ab_out__;
  ab = ab_out__;

  sgbtf2_(&m, &n, &kl, &ku, ab, &ldab, ipiv, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(3, rb_ipiv, rb_info, rb_ab);
}

void
init_lapack_sgbtf2(VALUE mLapack){
  rb_define_module_function(mLapack, "sgbtf2", rb_sgbtf2, -1);
}
