#include "rb_lapack.h"

extern VOID dlaed2_(integer *k, integer *n, integer *n1, doublereal *d, doublereal *q, integer *ldq, integer *indxq, doublereal *rho, doublereal *z, doublereal *dlamda, doublereal *w, doublereal *q2, integer *indx, integer *indxc, integer *indxp, integer *coltyp, integer *info);

static VALUE
rb_dlaed2(int argc, VALUE *argv, VALUE self){
  VALUE rb_n1;
  integer n1; 
  VALUE rb_d;
  doublereal *d; 
  VALUE rb_q;
  doublereal *q; 
  VALUE rb_indxq;
  integer *indxq; 
  VALUE rb_rho;
  doublereal rho; 
  VALUE rb_z;
  doublereal *z; 
  VALUE rb_k;
  integer k; 
  VALUE rb_dlamda;
  doublereal *dlamda; 
  VALUE rb_w;
  doublereal *w; 
  VALUE rb_q2;
  doublereal *q2; 
  VALUE rb_indxc;
  integer *indxc; 
  VALUE rb_coltyp;
  integer *coltyp; 
  VALUE rb_info;
  integer info; 
  VALUE rb_d_out__;
  doublereal *d_out__;
  VALUE rb_q_out__;
  doublereal *q_out__;
  VALUE rb_indxq_out__;
  integer *indxq_out__;
  integer *indx;
  integer *indxp;

  integer n;
  integer ldq;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  k, dlamda, w, q2, indxc, coltyp, info, d, q, indxq, rho = NumRu::Lapack.dlaed2( n1, d, q, indxq, rho, z)\n    or\n  NumRu::Lapack.dlaed2  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rb_n1 = argv[0];
  rb_d = argv[1];
  rb_q = argv[2];
  rb_indxq = argv[3];
  rb_rho = argv[4];
  rb_z = argv[5];

  if (!NA_IsNArray(rb_indxq))
    rb_raise(rb_eArgError, "indxq (4th argument) must be NArray");
  if (NA_RANK(rb_indxq) != 1)
    rb_raise(rb_eArgError, "rank of indxq (4th argument) must be %d", 1);
  n = NA_SHAPE0(rb_indxq);
  if (NA_TYPE(rb_indxq) != NA_LINT)
    rb_indxq = na_change_type(rb_indxq, NA_LINT);
  indxq = NA_PTR_TYPE(rb_indxq, integer*);
  rho = NUM2DBL(rb_rho);
  n1 = NUM2INT(rb_n1);
  if (!NA_IsNArray(rb_z))
    rb_raise(rb_eArgError, "z (6th argument) must be NArray");
  if (NA_RANK(rb_z) != 1)
    rb_raise(rb_eArgError, "rank of z (6th argument) must be %d", 1);
  if (NA_SHAPE0(rb_z) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of z must be the same as shape 0 of indxq");
  if (NA_TYPE(rb_z) != NA_DFLOAT)
    rb_z = na_change_type(rb_z, NA_DFLOAT);
  z = NA_PTR_TYPE(rb_z, doublereal*);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (2th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (2th argument) must be %d", 1);
  if (NA_SHAPE0(rb_d) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of d must be the same as shape 0 of indxq");
  if (NA_TYPE(rb_d) != NA_DFLOAT)
    rb_d = na_change_type(rb_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rb_d, doublereal*);
  if (!NA_IsNArray(rb_q))
    rb_raise(rb_eArgError, "q (3th argument) must be NArray");
  if (NA_RANK(rb_q) != 2)
    rb_raise(rb_eArgError, "rank of q (3th argument) must be %d", 2);
  if (NA_SHAPE1(rb_q) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of q must be the same as shape 0 of indxq");
  ldq = NA_SHAPE0(rb_q);
  if (NA_TYPE(rb_q) != NA_DFLOAT)
    rb_q = na_change_type(rb_q, NA_DFLOAT);
  q = NA_PTR_TYPE(rb_q, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_dlamda = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  dlamda = NA_PTR_TYPE(rb_dlamda, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_w = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  w = NA_PTR_TYPE(rb_w, doublereal*);
  {
    int shape[1];
    shape[0] = pow(n1,2)+pow(n-n1,2);
    rb_q2 = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  q2 = NA_PTR_TYPE(rb_q2, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_indxc = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  indxc = NA_PTR_TYPE(rb_indxc, integer*);
  {
    int shape[1];
    shape[0] = n;
    rb_coltyp = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  coltyp = NA_PTR_TYPE(rb_coltyp, integer*);
  {
    int shape[1];
    shape[0] = n;
    rb_d_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  d_out__ = NA_PTR_TYPE(rb_d_out__, doublereal*);
  MEMCPY(d_out__, d, doublereal, NA_TOTAL(rb_d));
  rb_d = rb_d_out__;
  d = d_out__;
  {
    int shape[2];
    shape[0] = ldq;
    shape[1] = n;
    rb_q_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  q_out__ = NA_PTR_TYPE(rb_q_out__, doublereal*);
  MEMCPY(q_out__, q, doublereal, NA_TOTAL(rb_q));
  rb_q = rb_q_out__;
  q = q_out__;
  {
    int shape[1];
    shape[0] = n;
    rb_indxq_out__ = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  indxq_out__ = NA_PTR_TYPE(rb_indxq_out__, integer*);
  MEMCPY(indxq_out__, indxq, integer, NA_TOTAL(rb_indxq));
  rb_indxq = rb_indxq_out__;
  indxq = indxq_out__;
  indx = ALLOC_N(integer, (n));
  indxp = ALLOC_N(integer, (n));

  dlaed2_(&k, &n, &n1, d, q, &ldq, indxq, &rho, z, dlamda, w, q2, indx, indxc, indxp, coltyp, &info);

  free(indx);
  free(indxp);
  rb_k = INT2NUM(k);
  rb_info = INT2NUM(info);
  rb_rho = rb_float_new((double)rho);
  return rb_ary_new3(11, rb_k, rb_dlamda, rb_w, rb_q2, rb_indxc, rb_coltyp, rb_info, rb_d, rb_q, rb_indxq, rb_rho);
}

void
init_lapack_dlaed2(VALUE mLapack){
  rb_define_module_function(mLapack, "dlaed2", rb_dlaed2, -1);
}
