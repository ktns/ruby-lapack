#include "rb_lapack.h"

extern VOID slahrd_(integer *n, integer *k, integer *nb, real *a, integer *lda, real *tau, real *t, integer *ldt, real *y, integer *ldy);

static VALUE
rb_slahrd(int argc, VALUE *argv, VALUE self){
  VALUE rb_n;
  integer n; 
  VALUE rb_k;
  integer k; 
  VALUE rb_nb;
  integer nb; 
  VALUE rb_a;
  real *a; 
  VALUE rb_tau;
  real *tau; 
  VALUE rb_t;
  real *t; 
  VALUE rb_y;
  real *y; 
  VALUE rb_a_out__;
  real *a_out__;

  integer lda;
  integer ldt;
  integer ldy;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  tau, t, y, a = NumRu::Lapack.slahrd( n, k, nb, a)\n    or\n  NumRu::Lapack.slahrd  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_n = argv[0];
  rb_k = argv[1];
  rb_nb = argv[2];
  rb_a = argv[3];

  k = NUM2INT(rb_k);
  nb = NUM2INT(rb_nb);
  n = NUM2INT(rb_n);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (4th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (4th argument) must be %d", 2);
  if (NA_SHAPE1(rb_a) != (n-k+1))
    rb_raise(rb_eRuntimeError, "shape 1 of a must be %d", n-k+1);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);
  ldt = nb;
  ldy = n;
  {
    int shape[1];
    shape[0] = MAX(1,nb);
    rb_tau = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  tau = NA_PTR_TYPE(rb_tau, real*);
  {
    int shape[2];
    shape[0] = ldt;
    shape[1] = MAX(1,nb);
    rb_t = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  t = NA_PTR_TYPE(rb_t, real*);
  {
    int shape[2];
    shape[0] = ldy;
    shape[1] = MAX(1,nb);
    rb_y = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  y = NA_PTR_TYPE(rb_y, real*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n-k+1;
    rb_a_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rb_a_out__, real*);
  MEMCPY(a_out__, a, real, NA_TOTAL(rb_a));
  rb_a = rb_a_out__;
  a = a_out__;

  slahrd_(&n, &k, &nb, a, &lda, tau, t, &ldt, y, &ldy);

  return rb_ary_new3(4, rb_tau, rb_t, rb_y, rb_a);
}

void
init_lapack_slahrd(VALUE mLapack){
  rb_define_module_function(mLapack, "slahrd", rb_slahrd, -1);
}
