#include "rb_lapack.h"

extern VOID claed0_(integer *qsiz, integer *n, real *d, real *e, complex *q, integer *ldq, complex *qstore, integer *ldqs, real *rwork, integer *iwork, integer *info);

static VALUE
rb_claed0(int argc, VALUE *argv, VALUE self){
  VALUE rb_qsiz;
  integer qsiz; 
  VALUE rb_d;
  real *d; 
  VALUE rb_e;
  real *e; 
  VALUE rb_q;
  complex *q; 
  VALUE rb_info;
  integer info; 
  VALUE rb_d_out__;
  real *d_out__;
  VALUE rb_e_out__;
  real *e_out__;
  VALUE rb_q_out__;
  complex *q_out__;
  complex *qstore;
  real *rwork;
  integer *iwork;

  integer n;
  integer ldq;
  integer ldqs;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, d, e, q = NumRu::Lapack.claed0( qsiz, d, e, q)\n    or\n  NumRu::Lapack.claed0  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_qsiz = argv[0];
  rb_d = argv[1];
  rb_e = argv[2];
  rb_q = argv[3];

  qsiz = NUM2INT(rb_qsiz);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (2th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (2th argument) must be %d", 1);
  n = NA_SHAPE0(rb_d);
  if (NA_TYPE(rb_d) != NA_SFLOAT)
    rb_d = na_change_type(rb_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rb_d, real*);
  if (!NA_IsNArray(rb_q))
    rb_raise(rb_eArgError, "q (4th argument) must be NArray");
  if (NA_RANK(rb_q) != 2)
    rb_raise(rb_eArgError, "rank of q (4th argument) must be %d", 2);
  if (NA_SHAPE1(rb_q) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of q must be the same as shape 0 of d");
  ldq = NA_SHAPE0(rb_q);
  if (NA_TYPE(rb_q) != NA_SCOMPLEX)
    rb_q = na_change_type(rb_q, NA_SCOMPLEX);
  q = NA_PTR_TYPE(rb_q, complex*);
  if (!NA_IsNArray(rb_e))
    rb_raise(rb_eArgError, "e (3th argument) must be NArray");
  if (NA_RANK(rb_e) != 1)
    rb_raise(rb_eArgError, "rank of e (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_e) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of e must be %d", n-1);
  if (NA_TYPE(rb_e) != NA_SFLOAT)
    rb_e = na_change_type(rb_e, NA_SFLOAT);
  e = NA_PTR_TYPE(rb_e, real*);
  ldqs = MAX(1,n);
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
    int shape[1];
    shape[0] = n-1;
    rb_e_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  e_out__ = NA_PTR_TYPE(rb_e_out__, real*);
  MEMCPY(e_out__, e, real, NA_TOTAL(rb_e));
  rb_e = rb_e_out__;
  e = e_out__;
  {
    int shape[2];
    shape[0] = ldq;
    shape[1] = n;
    rb_q_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  q_out__ = NA_PTR_TYPE(rb_q_out__, complex*);
  MEMCPY(q_out__, q, complex, NA_TOTAL(rb_q));
  rb_q = rb_q_out__;
  q = q_out__;
  qstore = ALLOC_N(complex, (ldqs)*(n));
  rwork = ALLOC_N(real, (1 + 3*n + 2*n*LG(n) + 3*pow(n,2)));
  iwork = ALLOC_N(integer, (6 + 6*n + 5*n*LG(n)));

  claed0_(&qsiz, &n, d, e, q, &ldq, qstore, &ldqs, rwork, iwork, &info);

  free(qstore);
  free(rwork);
  free(iwork);
  rb_info = INT2NUM(info);
  return rb_ary_new3(4, rb_info, rb_d, rb_e, rb_q);
}

void
init_lapack_claed0(VALUE mLapack){
  rb_define_module_function(mLapack, "claed0", rb_claed0, -1);
}
