#include "rb_lapack.h"

extern VOID sorgbr_(char *vect, integer *m, integer *n, integer *k, real *a, integer *lda, real *tau, real *work, integer *lwork, integer *info);

static VALUE
rb_sorgbr(int argc, VALUE *argv, VALUE self){
  VALUE rb_vect;
  char vect; 
  VALUE rb_m;
  integer m; 
  VALUE rb_k;
  integer k; 
  VALUE rb_a;
  real *a; 
  VALUE rb_tau;
  real *tau; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_work;
  real *work; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  real *a_out__;

  integer lda;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  work, info, a = NumRu::Lapack.sorgbr( vect, m, k, a, tau, lwork)\n    or\n  NumRu::Lapack.sorgbr  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rb_vect = argv[0];
  rb_m = argv[1];
  rb_k = argv[2];
  rb_a = argv[3];
  rb_tau = argv[4];
  rb_lwork = argv[5];

  k = NUM2INT(rb_k);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (4th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (4th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);
  m = NUM2INT(rb_m);
  lwork = NUM2INT(rb_lwork);
  vect = StringValueCStr(rb_vect)[0];
  if (!NA_IsNArray(rb_tau))
    rb_raise(rb_eArgError, "tau (5th argument) must be NArray");
  if (NA_RANK(rb_tau) != 1)
    rb_raise(rb_eArgError, "rank of tau (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_tau) != (MIN(m,k)))
    rb_raise(rb_eRuntimeError, "shape 0 of tau must be %d", MIN(m,k));
  if (NA_TYPE(rb_tau) != NA_SFLOAT)
    rb_tau = na_change_type(rb_tau, NA_SFLOAT);
  tau = NA_PTR_TYPE(rb_tau, real*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rb_work = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rb_work, real*);
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

  sorgbr_(&vect, &m, &n, &k, a, &lda, tau, work, &lwork, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(3, rb_work, rb_info, rb_a);
}

void
init_lapack_sorgbr(VALUE mLapack){
  rb_define_module_function(mLapack, "sorgbr", rb_sorgbr, -1);
}
