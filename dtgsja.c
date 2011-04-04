#include "rb_lapack.h"

extern VOID dtgsja_(char *jobu, char *jobv, char *jobq, integer *m, integer *p, integer *n, integer *k, integer *l, doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *tola, doublereal *tolb, doublereal *alpha, doublereal *beta, doublereal *u, integer *ldu, doublereal *v, integer *ldv, doublereal *q, integer *ldq, doublereal *work, integer *ncycle, integer *info);

static VALUE
rb_dtgsja(int argc, VALUE *argv, VALUE self){
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
  doublereal *a; 
  VALUE rb_b;
  doublereal *b; 
  VALUE rb_tola;
  doublereal tola; 
  VALUE rb_tolb;
  doublereal tolb; 
  VALUE rb_u;
  doublereal *u; 
  VALUE rb_v;
  doublereal *v; 
  VALUE rb_q;
  doublereal *q; 
  VALUE rb_alpha;
  doublereal *alpha; 
  VALUE rb_beta;
  doublereal *beta; 
  VALUE rb_ncycle;
  integer ncycle; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  doublereal *a_out__;
  VALUE rb_b_out__;
  doublereal *b_out__;
  VALUE rb_u_out__;
  doublereal *u_out__;
  VALUE rb_v_out__;
  doublereal *v_out__;
  VALUE rb_q_out__;
  doublereal *q_out__;
  doublereal *work;

  integer lda;
  integer n;
  integer ldb;
  integer ldu;
  integer m;
  integer ldv;
  integer p;
  integer ldq;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  alpha, beta, ncycle, info, a, b, u, v, q = NumRu::Lapack.dtgsja( jobu, jobv, jobq, k, l, a, b, tola, tolb, u, v, q)\n    or\n  NumRu::Lapack.dtgsja  # print help\n\n\nFORTRAN MANUAL\n\n");
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
  if (NA_TYPE(rb_v) != NA_DFLOAT)
    rb_v = na_change_type(rb_v, NA_DFLOAT);
  v = NA_PTR_TYPE(rb_v, doublereal*);
  k = NUM2INT(rb_k);
  jobq = StringValueCStr(rb_jobq)[0];
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (6th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (6th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_DFLOAT)
    rb_a = na_change_type(rb_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rb_a, doublereal*);
  l = NUM2INT(rb_l);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (7th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (7th argument) must be %d", 2);
  if (NA_SHAPE1(rb_b) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of b must be the same as shape 1 of a");
  ldb = NA_SHAPE0(rb_b);
  if (NA_TYPE(rb_b) != NA_DFLOAT)
    rb_b = na_change_type(rb_b, NA_DFLOAT);
  b = NA_PTR_TYPE(rb_b, doublereal*);
  jobu = StringValueCStr(rb_jobu)[0];
  jobv = StringValueCStr(rb_jobv)[0];
  if (!NA_IsNArray(rb_q))
    rb_raise(rb_eArgError, "q (12th argument) must be NArray");
  if (NA_RANK(rb_q) != 2)
    rb_raise(rb_eArgError, "rank of q (12th argument) must be %d", 2);
  if (NA_SHAPE1(rb_q) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of q must be the same as shape 1 of a");
  ldq = NA_SHAPE0(rb_q);
  if (NA_TYPE(rb_q) != NA_DFLOAT)
    rb_q = na_change_type(rb_q, NA_DFLOAT);
  q = NA_PTR_TYPE(rb_q, doublereal*);
  tola = NUM2DBL(rb_tola);
  tolb = NUM2DBL(rb_tolb);
  if (!NA_IsNArray(rb_u))
    rb_raise(rb_eArgError, "u (10th argument) must be NArray");
  if (NA_RANK(rb_u) != 2)
    rb_raise(rb_eArgError, "rank of u (10th argument) must be %d", 2);
  m = NA_SHAPE1(rb_u);
  ldu = NA_SHAPE0(rb_u);
  if (NA_TYPE(rb_u) != NA_DFLOAT)
    rb_u = na_change_type(rb_u, NA_DFLOAT);
  u = NA_PTR_TYPE(rb_u, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_alpha = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  alpha = NA_PTR_TYPE(rb_alpha, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_beta = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  beta = NA_PTR_TYPE(rb_beta, doublereal*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rb_a_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rb_a_out__, doublereal*);
  MEMCPY(a_out__, a, doublereal, NA_TOTAL(rb_a));
  rb_a = rb_a_out__;
  a = a_out__;
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = n;
    rb_b_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rb_b_out__, doublereal*);
  MEMCPY(b_out__, b, doublereal, NA_TOTAL(rb_b));
  rb_b = rb_b_out__;
  b = b_out__;
  {
    int shape[2];
    shape[0] = ldu;
    shape[1] = m;
    rb_u_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  u_out__ = NA_PTR_TYPE(rb_u_out__, doublereal*);
  MEMCPY(u_out__, u, doublereal, NA_TOTAL(rb_u));
  rb_u = rb_u_out__;
  u = u_out__;
  {
    int shape[2];
    shape[0] = ldv;
    shape[1] = p;
    rb_v_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  v_out__ = NA_PTR_TYPE(rb_v_out__, doublereal*);
  MEMCPY(v_out__, v, doublereal, NA_TOTAL(rb_v));
  rb_v = rb_v_out__;
  v = v_out__;
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
  work = ALLOC_N(doublereal, (2*n));

  dtgsja_(&jobu, &jobv, &jobq, &m, &p, &n, &k, &l, a, &lda, b, &ldb, &tola, &tolb, alpha, beta, u, &ldu, v, &ldv, q, &ldq, work, &ncycle, &info);

  free(work);
  rb_ncycle = INT2NUM(ncycle);
  rb_info = INT2NUM(info);
  return rb_ary_new3(9, rb_alpha, rb_beta, rb_ncycle, rb_info, rb_a, rb_b, rb_u, rb_v, rb_q);
}

void
init_lapack_dtgsja(VALUE mLapack){
  rb_define_module_function(mLapack, "dtgsja", rb_dtgsja, -1);
}
