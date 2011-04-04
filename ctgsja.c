#include "rb_lapack.h"

extern VOID ctgsja_(char *jobu, char *jobv, char *jobq, integer *m, integer *p, integer *n, integer *k, integer *l, complex *a, integer *lda, complex *b, integer *ldb, real *tola, real *tolb, real *alpha, real *beta, complex *u, integer *ldu, complex *v, integer *ldv, complex *q, integer *ldq, complex *work, integer *ncycle, integer *info);

static VALUE
rb_ctgsja(int argc, VALUE *argv, VALUE self){
  VALUE rb_jobu;
  char jobu; 
  VALUE rb_jobv;
  char jobv; 
  VALUE rb_jobq;
  char jobq; 
  VALUE rb_k;
  integer k; 
  VALUE rb_l;
  integer l; 
  VALUE rb_a;
  complex *a; 
  VALUE rb_b;
  complex *b; 
  VALUE rb_tola;
  real tola; 
  VALUE rb_tolb;
  real tolb; 
  VALUE rb_u;
  complex *u; 
  VALUE rb_v;
  complex *v; 
  VALUE rb_q;
  complex *q; 
  VALUE rb_alpha;
  real *alpha; 
  VALUE rb_beta;
  real *beta; 
  VALUE rb_ncycle;
  integer ncycle; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  complex *a_out__;
  VALUE rb_b_out__;
  complex *b_out__;
  VALUE rb_u_out__;
  complex *u_out__;
  VALUE rb_v_out__;
  complex *v_out__;
  VALUE rb_q_out__;
  complex *q_out__;
  complex *work;

  integer lda;
  integer n;
  integer ldb;
  integer ldu;
  integer m;
  integer ldv;
  integer p;
  integer ldq;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  alpha, beta, ncycle, info, a, b, u, v, q = NumRu::Lapack.ctgsja( jobu, jobv, jobq, k, l, a, b, tola, tolb, u, v, q)\n    or\n  NumRu::Lapack.ctgsja  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 12)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 12)", argc);
  rb_jobu = argv[0];
  rb_jobv = argv[1];
  rb_jobq = argv[2];
  rb_k = argv[3];
  rb_l = argv[4];
  rb_a = argv[5];
  rb_b = argv[6];
  rb_tola = argv[7];
  rb_tolb = argv[8];
  rb_u = argv[9];
  rb_v = argv[10];
  rb_q = argv[11];

  if (!NA_IsNArray(rb_v))
    rb_raise(rb_eArgError, "v (11th argument) must be NArray");
  if (NA_RANK(rb_v) != 2)
    rb_raise(rb_eArgError, "rank of v (11th argument) must be %d", 2);
  p = NA_SHAPE1(rb_v);
  ldv = NA_SHAPE0(rb_v);
  if (NA_TYPE(rb_v) != NA_SCOMPLEX)
    rb_v = na_change_type(rb_v, NA_SCOMPLEX);
  v = NA_PTR_TYPE(rb_v, complex*);
  k = NUM2INT(rb_k);
  jobq = StringValueCStr(rb_jobq)[0];
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (6th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (6th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SCOMPLEX)
    rb_a = na_change_type(rb_a, NA_SCOMPLEX);
  a = NA_PTR_TYPE(rb_a, complex*);
  l = NUM2INT(rb_l);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (7th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (7th argument) must be %d", 2);
  if (NA_SHAPE1(rb_b) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of b must be the same as shape 1 of a");
  ldb = NA_SHAPE0(rb_b);
  if (NA_TYPE(rb_b) != NA_SCOMPLEX)
    rb_b = na_change_type(rb_b, NA_SCOMPLEX);
  b = NA_PTR_TYPE(rb_b, complex*);
  jobu = StringValueCStr(rb_jobu)[0];
  jobv = StringValueCStr(rb_jobv)[0];
  if (!NA_IsNArray(rb_q))
    rb_raise(rb_eArgError, "q (12th argument) must be NArray");
  if (NA_RANK(rb_q) != 2)
    rb_raise(rb_eArgError, "rank of q (12th argument) must be %d", 2);
  if (NA_SHAPE1(rb_q) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of q must be the same as shape 1 of a");
  ldq = NA_SHAPE0(rb_q);
  if (NA_TYPE(rb_q) != NA_SCOMPLEX)
    rb_q = na_change_type(rb_q, NA_SCOMPLEX);
  q = NA_PTR_TYPE(rb_q, complex*);
  tola = (real)NUM2DBL(rb_tola);
  tolb = (real)NUM2DBL(rb_tolb);
  if (!NA_IsNArray(rb_u))
    rb_raise(rb_eArgError, "u (10th argument) must be NArray");
  if (NA_RANK(rb_u) != 2)
    rb_raise(rb_eArgError, "rank of u (10th argument) must be %d", 2);
  m = NA_SHAPE1(rb_u);
  ldu = NA_SHAPE0(rb_u);
  if (NA_TYPE(rb_u) != NA_SCOMPLEX)
    rb_u = na_change_type(rb_u, NA_SCOMPLEX);
  u = NA_PTR_TYPE(rb_u, complex*);
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
    shape[0] = lda;
    shape[1] = n;
    rb_a_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rb_a_out__, complex*);
  MEMCPY(a_out__, a, complex, NA_TOTAL(rb_a));
  rb_a = rb_a_out__;
  a = a_out__;
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = n;
    rb_b_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rb_b_out__, complex*);
  MEMCPY(b_out__, b, complex, NA_TOTAL(rb_b));
  rb_b = rb_b_out__;
  b = b_out__;
  {
    int shape[2];
    shape[0] = ldu;
    shape[1] = m;
    rb_u_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  u_out__ = NA_PTR_TYPE(rb_u_out__, complex*);
  MEMCPY(u_out__, u, complex, NA_TOTAL(rb_u));
  rb_u = rb_u_out__;
  u = u_out__;
  {
    int shape[2];
    shape[0] = ldv;
    shape[1] = p;
    rb_v_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  v_out__ = NA_PTR_TYPE(rb_v_out__, complex*);
  MEMCPY(v_out__, v, complex, NA_TOTAL(rb_v));
  rb_v = rb_v_out__;
  v = v_out__;
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
  work = ALLOC_N(complex, (2*n));

  ctgsja_(&jobu, &jobv, &jobq, &m, &p, &n, &k, &l, a, &lda, b, &ldb, &tola, &tolb, alpha, beta, u, &ldu, v, &ldv, q, &ldq, work, &ncycle, &info);

  free(work);
  rb_ncycle = INT2NUM(ncycle);
  rb_info = INT2NUM(info);
  return rb_ary_new3(9, rb_alpha, rb_beta, rb_ncycle, rb_info, rb_a, rb_b, rb_u, rb_v, rb_q);
}

void
init_lapack_ctgsja(VALUE mLapack){
  rb_define_module_function(mLapack, "ctgsja", rb_ctgsja, -1);
}
