#include "rb_lapack.h"

extern VOID dlaed3_(integer *k, integer *n, integer *n1, doublereal *d, doublereal *q, integer *ldq, doublereal *rho, doublereal *dlamda, doublereal *q2, integer *indx, integer *ctot, doublereal *w, doublereal *s, integer *info);

static VALUE
rb_dlaed3(int argc, VALUE *argv, VALUE self){
  VALUE rb_n1;
  integer n1; 
  VALUE rb_rho;
  doublereal rho; 
  VALUE rb_dlamda;
  doublereal *dlamda; 
  VALUE rb_q2;
  doublereal *q2; 
  VALUE rb_indx;
  integer *indx; 
  VALUE rb_ctot;
  integer *ctot; 
  VALUE rb_w;
  doublereal *w; 
  VALUE rb_d;
  doublereal *d; 
  VALUE rb_q;
  doublereal *q; 
  VALUE rb_info;
  integer info; 
  VALUE rb_dlamda_out__;
  doublereal *dlamda_out__;
  VALUE rb_w_out__;
  doublereal *w_out__;
  doublereal *s;

  integer k;
  integer n;
  integer ldq;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  d, q, info, dlamda, w = NumRu::Lapack.dlaed3( n1, rho, dlamda, q2, indx, ctot, w)\n    or\n  NumRu::Lapack.dlaed3  # print help\n\n\nFORTRAN MANUAL\n\n");
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
  if (NA_TYPE(rb_w) != NA_DFLOAT)
    rb_w = na_change_type(rb_w, NA_DFLOAT);
  w = NA_PTR_TYPE(rb_w, doublereal*);
  if (!NA_IsNArray(rb_q2))
    rb_raise(rb_eArgError, "q2 (4th argument) must be NArray");
  if (NA_RANK(rb_q2) != 2)
    rb_raise(rb_eArgError, "rank of q2 (4th argument) must be %d", 2);
  n = NA_SHAPE1(rb_q2);
  if (NA_SHAPE0(rb_q2) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of q2 must be the same as shape 1 of q2");
  if (NA_TYPE(rb_q2) != NA_DFLOAT)
    rb_q2 = na_change_type(rb_q2, NA_DFLOAT);
  q2 = NA_PTR_TYPE(rb_q2, doublereal*);
  if (!NA_IsNArray(rb_dlamda))
    rb_raise(rb_eArgError, "dlamda (3th argument) must be NArray");
  if (NA_RANK(rb_dlamda) != 1)
    rb_raise(rb_eArgError, "rank of dlamda (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_dlamda) != k)
    rb_raise(rb_eRuntimeError, "shape 0 of dlamda must be the same as shape 0 of w");
  if (NA_TYPE(rb_dlamda) != NA_DFLOAT)
    rb_dlamda = na_change_type(rb_dlamda, NA_DFLOAT);
  dlamda = NA_PTR_TYPE(rb_dlamda, doublereal*);
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
  rho = NUM2DBL(rb_rho);
  ldq = MAX(1,n);
  {
    int shape[1];
    shape[0] = n;
    rb_d = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  d = NA_PTR_TYPE(rb_d, doublereal*);
  {
    int shape[2];
    shape[0] = ldq;
    shape[1] = n;
    rb_q = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  q = NA_PTR_TYPE(rb_q, doublereal*);
  {
    int shape[1];
    shape[0] = k;
    rb_dlamda_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  dlamda_out__ = NA_PTR_TYPE(rb_dlamda_out__, doublereal*);
  MEMCPY(dlamda_out__, dlamda, doublereal, NA_TOTAL(rb_dlamda));
  rb_dlamda = rb_dlamda_out__;
  dlamda = dlamda_out__;
  {
    int shape[1];
    shape[0] = k;
    rb_w_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  w_out__ = NA_PTR_TYPE(rb_w_out__, doublereal*);
  MEMCPY(w_out__, w, doublereal, NA_TOTAL(rb_w));
  rb_w = rb_w_out__;
  w = w_out__;
  s = ALLOC_N(doublereal, (MAX(1,k))*((n1 + 1)));

  dlaed3_(&k, &n, &n1, d, q, &ldq, &rho, dlamda, q2, indx, ctot, w, s, &info);

  free(s);
  rb_info = INT2NUM(info);
  return rb_ary_new3(5, rb_d, rb_q, rb_info, rb_dlamda, rb_w);
}

void
init_lapack_dlaed3(VALUE mLapack){
  rb_define_module_function(mLapack, "dlaed3", rb_dlaed3, -1);
}
