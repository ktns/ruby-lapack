#include "rb_lapack.h"

extern VOID sggsvd_(char *jobu, char *jobv, char *jobq, integer *m, integer *n, integer *p, integer *k, integer *l, real *a, integer *lda, real *b, integer *ldb, real *alpha, real *beta, real *u, integer *ldu, real *v, integer *ldv, real *q, integer *ldq, real *work, integer *iwork, integer *info);

static VALUE
rb_sggsvd(int argc, VALUE *argv, VALUE self){
  VALUE rb_jobu;
  char jobu; 
  VALUE rb_jobv;
  char jobv; 
  VALUE rb_jobq;
  char jobq; 
  VALUE rb_a;
  real *a; 
  VALUE rb_b;
  real *b; 
  VALUE rb_k;
  integer k; 
  VALUE rb_l;
  integer l; 
  VALUE rb_alpha;
  real *alpha; 
  VALUE rb_beta;
  real *beta; 
  VALUE rb_u;
  real *u; 
  VALUE rb_v;
  real *v; 
  VALUE rb_q;
  real *q; 
  VALUE rb_iwork;
  integer *iwork; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  real *a_out__;
  VALUE rb_b_out__;
  real *b_out__;
  real *work;

  integer lda;
  integer n;
  integer ldb;
  integer ldu;
  integer m;
  integer ldv;
  integer p;
  integer ldq;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  k, l, alpha, beta, u, v, q, iwork, info, a, b = NumRu::Lapack.sggsvd( jobu, jobv, jobq, a, b)\n    or\n  NumRu::Lapack.sggsvd  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_jobu = argv[0];
  rb_jobv = argv[1];
  rb_jobq = argv[2];
  rb_a = argv[3];
  rb_b = argv[4];

  jobq = StringValueCStr(rb_jobq)[0];
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (4th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (4th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (lda != (m))
    rb_raise(rb_eRuntimeError, "shape 0 of a must be %d", m);
  m = lda;
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (5th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (5th argument) must be %d", 2);
  if (NA_SHAPE1(rb_b) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of b must be the same as shape 1 of a");
  ldb = NA_SHAPE0(rb_b);
  if (ldb != (p))
    rb_raise(rb_eRuntimeError, "shape 0 of b must be %d", p);
  p = ldb;
  if (NA_TYPE(rb_b) != NA_SFLOAT)
    rb_b = na_change_type(rb_b, NA_SFLOAT);
  b = NA_PTR_TYPE(rb_b, real*);
  jobu = StringValueCStr(rb_jobu)[0];
  jobv = StringValueCStr(rb_jobv)[0];
  ldq = lsame_(&jobq,"Q") ? MAX(1,n) : 1;
  ldu = lsame_(&jobu,"U") ? MAX(1,m) : 1;
  ldv = lsame_(&jobv,"V") ? MAX(1,p) : 1;
  lda = m;
  ldb = p;
  {
    int shape[1];
    shape[0] = n;
    rb_alpha = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  alpha = NA_PTR_TYPE(rb_alpha, real*);
  {
    int shape[1];
    shape[0] = n;
    rb_beta = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  beta = NA_PTR_TYPE(rb_beta, real*);
  {
    int shape[2];
    shape[0] = ldu;
    shape[1] = m;
    rb_u = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  u = NA_PTR_TYPE(rb_u, real*);
  {
    int shape[2];
    shape[0] = ldv;
    shape[1] = p;
    rb_v = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  v = NA_PTR_TYPE(rb_v, real*);
  {
    int shape[2];
    shape[0] = ldq;
    shape[1] = n;
    rb_q = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  q = NA_PTR_TYPE(rb_q, real*);
  {
    int shape[1];
    shape[0] = n;
    rb_iwork = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  iwork = NA_PTR_TYPE(rb_iwork, integer*);
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
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = n;
    rb_b_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rb_b_out__, real*);
  MEMCPY(b_out__, b, real, NA_TOTAL(rb_b));
  rb_b = rb_b_out__;
  b = b_out__;
  work = ALLOC_N(real, (MAX(3*n,m)*(p)+n));

  sggsvd_(&jobu, &jobv, &jobq, &m, &n, &p, &k, &l, a, &lda, b, &ldb, alpha, beta, u, &ldu, v, &ldv, q, &ldq, work, iwork, &info);

  free(work);
  rb_k = INT2NUM(k);
  rb_l = INT2NUM(l);
  rb_info = INT2NUM(info);
  return rb_ary_new3(11, rb_k, rb_l, rb_alpha, rb_beta, rb_u, rb_v, rb_q, rb_iwork, rb_info, rb_a, rb_b);
}

void
init_lapack_sggsvd(VALUE mLapack){
  rb_define_module_function(mLapack, "sggsvd", rb_sggsvd, -1);
}
