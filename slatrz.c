#include "rb_lapack.h"

extern VOID slatrz_(integer *m, integer *n, integer *l, real *a, integer *lda, real *tau, real *work);

static VALUE
rb_slatrz(int argc, VALUE *argv, VALUE self){
  VALUE rb_l;
  integer l; 
  VALUE rb_a;
  real *a; 
  VALUE rb_tau;
  real *tau; 
  VALUE rb_a_out__;
  real *a_out__;
  real *work;

  integer lda;
  integer n;
  integer m;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  tau, a = NumRu::Lapack.slatrz( l, a)\n    or\n  NumRu::Lapack.slatrz  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rb_l = argv[0];
  rb_a = argv[1];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);
  l = NUM2INT(rb_l);
  m = lda;
  {
    int shape[1];
    shape[0] = m;
    rb_tau = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  tau = NA_PTR_TYPE(rb_tau, real*);
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
  work = ALLOC_N(real, (m));

  slatrz_(&m, &n, &l, a, &lda, tau, work);

  free(work);
  return rb_ary_new3(2, rb_tau, rb_a);
}

void
init_lapack_slatrz(VALUE mLapack){
  rb_define_module_function(mLapack, "slatrz", rb_slatrz, -1);
}
