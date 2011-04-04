#include "rb_lapack.h"

extern VOID slaed3_(integer *k, integer *n, integer *n1, real *d, real *q, integer *ldq, real *rho, real *dlamda, real *q2, integer *indx, integer *ctot, real *w, real *s, integer *info);

static VALUE
rb_slaed3(int argc, VALUE *argv, VALUE self){
  VALUE rb_n1;
  integer n1; 
  VALUE rb_rho;
  real rho; 
  VALUE rb_dlamda;
  real *dlamda; 
  VALUE rb_q2;
  real *q2; 
  VALUE rb_indx;
  integer *indx; 
  VALUE rb_ctot;
  integer *ctot; 
  VALUE rb_w;
  real *w; 
  VALUE rb_d;
  real *d; 
  VALUE rb_q;
  real *q; 
  VALUE rb_info;
  integer info; 
  VALUE rb_dlamda_out__;
  real *dlamda_out__;
  VALUE rb_w_out__;
  real *w_out__;
  real *s;

  integer k;
  integer n;
  integer ldq;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  d, q, info, dlamda, w = NumRu::Lapack.slaed3( n1, rho, dlamda, q2, indx, ctot, w)\n    or\n  NumRu::Lapack.slaed3  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 7)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 7)", argc);
  rb_n1 = argv[0];
  rb_rho = argv[1];
  rb_dlamda = argv[2];
  rb_q2 = argv[3];
  rb_indx = argv[4];
  rb_ctot = argv[5];
  rb_w = argv[6];

  if (!NA_IsNArray(rb_ctot))
    rb_raise(rb_eArgError, "ctot (6th argument) must be NArray");
  if (NA_RANK(rb_ctot) != 1)
    rb_raise(rb_eArgError, "rank of ctot (6th argument) must be %d", 1);
  if (NA_SHAPE0(rb_ctot) != (4))
    rb_raise(rb_eRuntimeError, "shape 0 of ctot must be %d", 4);
  if (NA_TYPE(rb_ctot) != NA_LINT)
    rb_ctot = na_change_type(rb_ctot, NA_LINT);
  ctot = NA_PTR_TYPE(rb_ctot, integer*);
  if (!NA_IsNArray(rb_w))
    rb_raise(rb_eArgError, "w (7th argument) must be NArray");
  if (NA_RANK(rb_w) != 1)
    rb_raise(rb_eArgError, "rank of w (7th argument) must be %d", 1);
  k = NA_SHAPE0(rb_w);
  if (NA_TYPE(rb_w) != NA_SFLOAT)
    rb_w = na_change_type(rb_w, NA_SFLOAT);
  w = NA_PTR_TYPE(rb_w, real*);
  if (!NA_IsNArray(rb_q2))
    rb_raise(rb_eArgError, "q2 (4th argument) must be NArray");
  if (NA_RANK(rb_q2) != 2)
    rb_raise(rb_eArgError, "rank of q2 (4th argument) must be %d", 2);
  n = NA_SHAPE1(rb_q2);
  if (NA_SHAPE0(rb_q2) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of q2 must be the same as shape 1 of q2");
  if (NA_TYPE(rb_q2) != NA_SFLOAT)
    rb_q2 = na_change_type(rb_q2, NA_SFLOAT);
  q2 = NA_PTR_TYPE(rb_q2, real*);
  if (!NA_IsNArray(rb_dlamda))
    rb_raise(rb_eArgError, "dlamda (3th argument) must be NArray");
  if (NA_RANK(rb_dlamda) != 1)
    rb_raise(rb_eArgError, "rank of dlamda (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_dlamda) != k)
    rb_raise(rb_eRuntimeError, "shape 0 of dlamda must be the same as shape 0 of w");
  if (NA_TYPE(rb_dlamda) != NA_SFLOAT)
    rb_dlamda = na_change_type(rb_dlamda, NA_SFLOAT);
  dlamda = NA_PTR_TYPE(rb_dlamda, real*);
  if (!NA_IsNArray(rb_indx))
    rb_raise(rb_eArgError, "indx (5th argument) must be NArray");
  if (NA_RANK(rb_indx) != 1)
    rb_raise(rb_eArgError, "rank of indx (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_indx) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of indx must be the same as shape 1 of q2");
  if (NA_TYPE(rb_indx) != NA_LINT)
    rb_indx = na_change_type(rb_indx, NA_LINT);
  indx = NA_PTR_TYPE(rb_indx, integer*);
  n1 = NUM2INT(rb_n1);
  rho = (real)NUM2DBL(rb_rho);
  ldq = MAX(1,n);
  {
    int shape[1];
    shape[0] = n;
    rb_d = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  d = NA_PTR_TYPE(rb_d, real*);
  {
    int shape[2];
    shape[0] = ldq;
    shape[1] = n;
    rb_q = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  q = NA_PTR_TYPE(rb_q, real*);
  {
    int shape[1];
    shape[0] = k;
    rb_dlamda_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  dlamda_out__ = NA_PTR_TYPE(rb_dlamda_out__, real*);
  MEMCPY(dlamda_out__, dlamda, real, NA_TOTAL(rb_dlamda));
  rb_dlamda = rb_dlamda_out__;
  dlamda = dlamda_out__;
  {
    int shape[1];
    shape[0] = k;
    rb_w_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  w_out__ = NA_PTR_TYPE(rb_w_out__, real*);
  MEMCPY(w_out__, w, real, NA_TOTAL(rb_w));
  rb_w = rb_w_out__;
  w = w_out__;
  s = ALLOC_N(real, (MAX(1,k))*(n1 + 1));

  slaed3_(&k, &n, &n1, d, q, &ldq, &rho, dlamda, q2, indx, ctot, w, s, &info);

  free(s);
  rb_info = INT2NUM(info);
  return rb_ary_new3(5, rb_d, rb_q, rb_info, rb_dlamda, rb_w);
}

void
init_lapack_slaed3(VALUE mLapack){
  rb_define_module_function(mLapack, "slaed3", rb_slaed3, -1);
}
