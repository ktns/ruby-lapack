#include "rb_lapack.h"

extern VOID slaed0_(integer *icompq, integer *qsiz, integer *n, real *d, real *e, real *q, integer *ldq, real *qstore, integer *ldqs, real *work, integer *iwork, integer *info);

static VALUE
rb_slaed0(int argc, VALUE *argv, VALUE self){
  VALUE rb_icompq;
  integer icompq; 
  VALUE rb_qsiz;
  integer qsiz; 
  VALUE rb_d;
  real *d; 
  VALUE rb_e;
  real *e; 
  VALUE rb_q;
  real *q; 
  VALUE rb_info;
  integer info; 
  VALUE rb_d_out__;
  real *d_out__;
  VALUE rb_q_out__;
  real *q_out__;
  real *qstore;
  real *work;
  integer *iwork;

  integer n;
  integer ldq;
  integer ldqs;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, d, q = NumRu::Lapack.slaed0( icompq, qsiz, d, e, q)\n    or\n  NumRu::Lapack.slaed0  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_icompq = argv[0];
  rb_qsiz = argv[1];
  rb_d = argv[2];
  rb_e = argv[3];
  rb_q = argv[4];

  qsiz = NUM2INT(rb_qsiz);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (3th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (3th argument) must be %d", 1);
  n = NA_SHAPE0(rb_d);
  if (NA_TYPE(rb_d) != NA_SFLOAT)
    rb_d = na_change_type(rb_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rb_d, real*);
  if (!NA_IsNArray(rb_q))
    rb_raise(rb_eArgError, "q (5th argument) must be NArray");
  if (NA_RANK(rb_q) != 2)
    rb_raise(rb_eArgError, "rank of q (5th argument) must be %d", 2);
  if (NA_SHAPE1(rb_q) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of q must be the same as shape 0 of d");
  ldq = NA_SHAPE0(rb_q);
  if (NA_TYPE(rb_q) != NA_SFLOAT)
    rb_q = na_change_type(rb_q, NA_SFLOAT);
  q = NA_PTR_TYPE(rb_q, real*);
  icompq = NUM2INT(rb_icompq);
  if (!NA_IsNArray(rb_e))
    rb_raise(rb_eArgError, "e (4th argument) must be NArray");
  if (NA_RANK(rb_e) != 1)
    rb_raise(rb_eArgError, "rank of e (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_e) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of e must be %d", n-1);
  if (NA_TYPE(rb_e) != NA_SFLOAT)
    rb_e = na_change_type(rb_e, NA_SFLOAT);
  e = NA_PTR_TYPE(rb_e, real*);
  ldqs = icompq == 1 ? MAX(1,n) : 1;
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
  qstore = ALLOC_N(real, (ldqs)*(n));
  work = ALLOC_N(real, (((icompq == 0) || (icompq == 1)) ? 1 + 3*n + 2*n*LG(n) + 2*pow(n,2) : icompq == 2 ? 4*n + pow(n,2) : 0));
  iwork = ALLOC_N(integer, (((icompq == 0) || (icompq == 1)) ? 6 + 6*n + 5*n*LG(n) : icompq == 2 ? 3 + 5*n : 0));

  slaed0_(&icompq, &qsiz, &n, d, e, q, &ldq, qstore, &ldqs, work, iwork, &info);

  free(qstore);
  free(work);
  free(iwork);
  rb_info = INT2NUM(info);
  return rb_ary_new3(3, rb_info, rb_d, rb_q);
}

void
init_lapack_slaed0(VALUE mLapack){
  rb_define_module_function(mLapack, "slaed0", rb_slaed0, -1);
}
