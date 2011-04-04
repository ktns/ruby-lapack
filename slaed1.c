#include "rb_lapack.h"

extern VOID slaed1_(integer *n, real *d, real *q, integer *ldq, integer *indxq, real *rho, integer *cutpnt, real *work, integer *iwork, integer *info);

static VALUE
rb_slaed1(int argc, VALUE *argv, VALUE self){
  VALUE rb_d;
  real *d; 
  VALUE rb_q;
  real *q; 
  VALUE rb_indxq;
  integer *indxq; 
  VALUE rb_rho;
  real rho; 
  VALUE rb_cutpnt;
  integer cutpnt; 
  VALUE rb_info;
  integer info; 
  VALUE rb_d_out__;
  real *d_out__;
  VALUE rb_q_out__;
  real *q_out__;
  VALUE rb_indxq_out__;
  integer *indxq_out__;
  real *work;
  integer *iwork;

  integer n;
  integer ldq;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, d, q, indxq = NumRu::Lapack.slaed1( d, q, indxq, rho, cutpnt)\n    or\n  NumRu::Lapack.slaed1  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_d = argv[0];
  rb_q = argv[1];
  rb_indxq = argv[2];
  rb_rho = argv[3];
  rb_cutpnt = argv[4];

  cutpnt = NUM2INT(rb_cutpnt);
  if (!NA_IsNArray(rb_indxq))
    rb_raise(rb_eArgError, "indxq (3th argument) must be NArray");
  if (NA_RANK(rb_indxq) != 1)
    rb_raise(rb_eArgError, "rank of indxq (3th argument) must be %d", 1);
  n = NA_SHAPE0(rb_indxq);
  if (NA_TYPE(rb_indxq) != NA_LINT)
    rb_indxq = na_change_type(rb_indxq, NA_LINT);
  indxq = NA_PTR_TYPE(rb_indxq, integer*);
  rho = (real)NUM2DBL(rb_rho);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (1th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (1th argument) must be %d", 1);
  if (NA_SHAPE0(rb_d) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of d must be the same as shape 0 of indxq");
  if (NA_TYPE(rb_d) != NA_SFLOAT)
    rb_d = na_change_type(rb_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rb_d, real*);
  if (!NA_IsNArray(rb_q))
    rb_raise(rb_eArgError, "q (2th argument) must be NArray");
  if (NA_RANK(rb_q) != 2)
    rb_raise(rb_eArgError, "rank of q (2th argument) must be %d", 2);
  if (NA_SHAPE1(rb_q) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of q must be the same as shape 0 of indxq");
  ldq = NA_SHAPE0(rb_q);
  if (NA_TYPE(rb_q) != NA_SFLOAT)
    rb_q = na_change_type(rb_q, NA_SFLOAT);
  q = NA_PTR_TYPE(rb_q, real*);
  {
    int shape[1];
    shape[0] = n;
    rb_d_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  d_out__ = NA_PTR_TYPE(rb_d_out__, real*);
  MEMCPY(d_out__, d, real, NA_TOTAL(rb_d));
  rb_d = rb_d_out__;
  d = d_out__;
  {
    int shape[2];
    shape[0] = ldq;
    shape[1] = n;
    rb_q_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  q_out__ = NA_PTR_TYPE(rb_q_out__, real*);
  MEMCPY(q_out__, q, real, NA_TOTAL(rb_q));
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
  work = ALLOC_N(real, (4*n + pow(n,2)));
  iwork = ALLOC_N(integer, (4*n));

  slaed1_(&n, d, q, &ldq, indxq, &rho, &cutpnt, work, iwork, &info);

  free(work);
  free(iwork);
  rb_info = INT2NUM(info);
  return rb_ary_new3(4, rb_info, rb_d, rb_q, rb_indxq);
}

void
init_lapack_slaed1(VALUE mLapack){
  rb_define_module_function(mLapack, "slaed1", rb_slaed1, -1);
}
