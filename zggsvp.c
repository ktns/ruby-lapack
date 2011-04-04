#include "rb_lapack.h"

extern VOID zggsvp_(char *jobu, char *jobv, char *jobq, integer *m, integer *p, integer *n, doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, doublereal *tola, doublereal *tolb, integer *k, integer *l, doublecomplex *u, integer *ldu, doublecomplex *v, integer *ldv, doublecomplex *q, integer *ldq, integer *iwork, doublereal *rwork, doublecomplex *tau, doublecomplex *work, integer *info);

static VALUE
rb_zggsvp(int argc, VALUE *argv, VALUE self){
  VALUE rb_jobu;
  char jobu; 
  VALUE rb_jobv;
  char jobv; 
  VALUE rb_jobq;
  char jobq; 
  VALUE rb_a;
  doublecomplex *a; 
  VALUE rb_b;
  doublecomplex *b; 
  VALUE rb_tola;
  doublereal tola; 
  VALUE rb_tolb;
  doublereal tolb; 
  VALUE rb_k;
  integer k; 
  VALUE rb_l;
  integer l; 
  VALUE rb_u;
  doublecomplex *u; 
  VALUE rb_v;
  doublecomplex *v; 
  VALUE rb_q;
  doublecomplex *q; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  doublecomplex *a_out__;
  VALUE rb_b_out__;
  doublecomplex *b_out__;
  integer *iwork;
  doublereal *rwork;
  doublecomplex *tau;
  doublecomplex *work;

  integer lda;
  integer n;
  integer ldb;
  integer ldu;
  integer m;
  integer ldv;
  integer p;
  integer ldq;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  k, l, u, v, q, info, a, b = NumRu::Lapack.zggsvp( jobu, jobv, jobq, a, b, tola, tolb)\n    or\n  NumRu::Lapack.zggsvp  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 7)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 7)", argc);
  rb_jobu = argv[0];
  rb_jobv = argv[1];
  rb_jobq = argv[2];
  rb_a = argv[3];
  rb_b = argv[4];
  rb_tola = argv[5];
  rb_tolb = argv[6];

  jobq = StringValueCStr(rb_jobq)[0];
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (4th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (4th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_DCOMPLEX)
    rb_a = na_change_type(rb_a, NA_DCOMPLEX);
  a = NA_PTR_TYPE(rb_a, doublecomplex*);
  tolb = NUM2DBL(rb_tolb);
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
  if (NA_TYPE(rb_b) != NA_DCOMPLEX)
    rb_b = na_change_type(rb_b, NA_DCOMPLEX);
  b = NA_PTR_TYPE(rb_b, doublecomplex*);
  tola = NUM2DBL(rb_tola);
  jobu = StringValueCStr(rb_jobu)[0];
  jobv = StringValueCStr(rb_jobv)[0];
  ldu = lsame_(&jobu,"U") ? MAX(1,m) : 1;
  ldq = lsame_(&jobq,"Q") ? MAX(1,n) : 1;
  ldv = lsame_(&jobv,"V") ? MAX(1,p) : 1;
  m = lda;
  ldb = p;
  {
    int shape[2];
    shape[0] = ldu;
    shape[1] = m;
    rb_u = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  u = NA_PTR_TYPE(rb_u, doublecomplex*);
  {
    int shape[2];
    shape[0] = ldv;
    shape[1] = p;
    rb_v = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  v = NA_PTR_TYPE(rb_v, doublecomplex*);
  {
    int shape[2];
    shape[0] = ldq;
    shape[1] = n;
    rb_q = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  q = NA_PTR_TYPE(rb_q, doublecomplex*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rb_a_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rb_a_out__, doublecomplex*);
  MEMCPY(a_out__, a, doublecomplex, NA_TOTAL(rb_a));
  rb_a = rb_a_out__;
  a = a_out__;
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = n;
    rb_b_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rb_b_out__, doublecomplex*);
  MEMCPY(b_out__, b, doublecomplex, NA_TOTAL(rb_b));
  rb_b = rb_b_out__;
  b = b_out__;
  iwork = ALLOC_N(integer, (n));
  rwork = ALLOC_N(doublereal, (2*n));
  tau = ALLOC_N(doublecomplex, (n));
  work = ALLOC_N(doublecomplex, (MAX(3*n,m)*(p)));

  zggsvp_(&jobu, &jobv, &jobq, &m, &p, &n, a, &lda, b, &ldb, &tola, &tolb, &k, &l, u, &ldu, v, &ldv, q, &ldq, iwork, rwork, tau, work, &info);

  free(iwork);
  free(rwork);
  free(tau);
  free(work);
  rb_k = INT2NUM(k);
  rb_l = INT2NUM(l);
  rb_info = INT2NUM(info);
  return rb_ary_new3(8, rb_k, rb_l, rb_u, rb_v, rb_q, rb_info, rb_a, rb_b);
}

void
init_lapack_zggsvp(VALUE mLapack){
  rb_define_module_function(mLapack, "zggsvp", rb_zggsvp, -1);
}
