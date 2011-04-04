#include "rb_lapack.h"

extern VOID dgesc2_(integer *n, doublereal *a, integer *lda, doublereal *rhs, integer *ipiv, integer *jpiv, doublereal *scale);

static VALUE
rb_dgesc2(int argc, VALUE *argv, VALUE self){
  VALUE rb_a;
  doublereal *a; 
  VALUE rb_rhs;
  doublereal *rhs; 
  VALUE rb_ipiv;
  integer *ipiv; 
  VALUE rb_jpiv;
  integer *jpiv; 
  VALUE rb_scale;
  doublereal scale; 
  VALUE rb_rhs_out__;
  doublereal *rhs_out__;

  integer lda;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  scale, rhs = NumRu::Lapack.dgesc2( a, rhs, ipiv, jpiv)\n    or\n  NumRu::Lapack.dgesc2  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_a = argv[0];
  rb_rhs = argv[1];
  rb_ipiv = argv[2];
  rb_jpiv = argv[3];

  if (!NA_IsNArray(rb_ipiv))
    rb_raise(rb_eArgError, "ipiv (3th argument) must be NArray");
  if (NA_RANK(rb_ipiv) != 1)
    rb_raise(rb_eArgError, "rank of ipiv (3th argument) must be %d", 1);
  n = NA_SHAPE0(rb_ipiv);
  if (NA_TYPE(rb_ipiv) != NA_LINT)
    rb_ipiv = na_change_type(rb_ipiv, NA_LINT);
  ipiv = NA_PTR_TYPE(rb_ipiv, integer*);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (1th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (1th argument) must be %d", 2);
  if (NA_SHAPE1(rb_a) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of a must be the same as shape 0 of ipiv");
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_DFLOAT)
    rb_a = na_change_type(rb_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rb_a, doublereal*);
  if (!NA_IsNArray(rb_rhs))
    rb_raise(rb_eArgError, "rhs (2th argument) must be NArray");
  if (NA_RANK(rb_rhs) != 1)
    rb_raise(rb_eArgError, "rank of rhs (2th argument) must be %d", 1);
  if (NA_SHAPE0(rb_rhs) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of rhs must be the same as shape 0 of ipiv");
  if (NA_TYPE(rb_rhs) != NA_DFLOAT)
    rb_rhs = na_change_type(rb_rhs, NA_DFLOAT);
  rhs = NA_PTR_TYPE(rb_rhs, doublereal*);
  if (!NA_IsNArray(rb_jpiv))
    rb_raise(rb_eArgError, "jpiv (4th argument) must be NArray");
  if (NA_RANK(rb_jpiv) != 1)
    rb_raise(rb_eArgError, "rank of jpiv (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_jpiv) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of jpiv must be the same as shape 0 of ipiv");
  if (NA_TYPE(rb_jpiv) != NA_LINT)
    rb_jpiv = na_change_type(rb_jpiv, NA_LINT);
  jpiv = NA_PTR_TYPE(rb_jpiv, integer*);
  {
    int shape[1];
    shape[0] = n;
    rb_rhs_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  rhs_out__ = NA_PTR_TYPE(rb_rhs_out__, doublereal*);
  MEMCPY(rhs_out__, rhs, doublereal, NA_TOTAL(rb_rhs));
  rb_rhs = rb_rhs_out__;
  rhs = rhs_out__;

  dgesc2_(&n, a, &lda, rhs, ipiv, jpiv, &scale);

  rb_scale = rb_float_new((double)scale);
  return rb_ary_new3(2, rb_scale, rb_rhs);
}

void
init_lapack_dgesc2(VALUE mLapack){
  rb_define_module_function(mLapack, "dgesc2", rb_dgesc2, -1);
}
