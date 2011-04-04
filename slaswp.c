#include "rb_lapack.h"

extern VOID slaswp_(integer *n, real *a, integer *lda, integer *k1, integer *k2, integer *ipiv, integer *incx);

static VALUE
rb_slaswp(int argc, VALUE *argv, VALUE self){
  VALUE rb_a;
  real *a; 
  VALUE rb_k1;
  integer k1; 
  VALUE rb_k2;
  integer k2; 
  VALUE rb_ipiv;
  integer *ipiv; 
  VALUE rb_incx;
  integer incx; 
  VALUE rb_a_out__;
  real *a_out__;

  integer lda;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  a = NumRu::Lapack.slaswp( a, k1, k2, ipiv, incx)\n    or\n  NumRu::Lapack.slaswp  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_a = argv[0];
  rb_k1 = argv[1];
  rb_k2 = argv[2];
  rb_ipiv = argv[3];
  rb_incx = argv[4];

  k2 = NUM2INT(rb_k2);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (1th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (1th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);
  k1 = NUM2INT(rb_k1);
  incx = NUM2INT(rb_incx);
  if (!NA_IsNArray(rb_ipiv))
    rb_raise(rb_eArgError, "ipiv (4th argument) must be NArray");
  if (NA_RANK(rb_ipiv) != 1)
    rb_raise(rb_eArgError, "rank of ipiv (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_ipiv) != (k2*abs(incx)))
    rb_raise(rb_eRuntimeError, "shape 0 of ipiv must be %d", k2*abs(incx));
  if (NA_TYPE(rb_ipiv) != NA_LINT)
    rb_ipiv = na_change_type(rb_ipiv, NA_LINT);
  ipiv = NA_PTR_TYPE(rb_ipiv, integer*);
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

  slaswp_(&n, a, &lda, &k1, &k2, ipiv, &incx);

  return rb_a;
}

void
init_lapack_slaswp(VALUE mLapack){
  rb_define_module_function(mLapack, "slaswp", rb_slaswp, -1);
}
