#include "rb_lapack.h"

extern VOID cggsvp_(char *jobu, char *jobv, char *jobq, integer *m, integer *p, integer *n, complex *a, integer *lda, complex *b, integer *ldb, real *tola, real *tolb, integer *k, integer *l, complex *u, integer *ldu, complex *v, integer *ldv, complex *q, integer *ldq, integer *iwork, real *rwork, complex *tau, complex *work, integer *info);

static VALUE
rb_cggsvp(int argc, VALUE *argv, VALUE self){
  VALUE rb_jobu;
  char jobu; 
  VALUE rb_jobv;
  char jobv; 
  VALUE rb_jobq;
  char jobq; 
  VALUE rb_a;
  complex *a; 
  VALUE rb_b;
  complex *b; 
  VALUE rb_tola;
  real tola; 
  VALUE rb_tolb;
  real tolb; 
  VALUE rb_k;
  integer k; 
  VALUE rb_l;
  integer l; 
  VALUE rb_u;
  complex *u; 
  VALUE rb_v;
  complex *v; 
  VALUE rb_q;
  complex *q; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  complex *a_out__;
  VALUE rb_b_out__;
  complex *b_out__;
  integer *iwork;
  real *rwork;
  complex *tau;
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
    printf("%s\n", "USAGE:\n  k, l, u, v, q, info, a, b = NumRu::Lapack.cggsvp( jobu, jobv, jobq, a, b, tola, tolb)\n    or\n  NumRu::Lapack.cggsvp  # print help\n\n\nFORTRAN MANUAL\n\n");
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
  if (lda != (m))
    rb_raise(rb_eRuntimeError, "shape 0 of a must be %d", m);
  m = lda;
  if (NA_TYPE(rb_a) != NA_SCOMPLEX)
    rb_a = na_change_type(rb_a, NA_SCOMPLEX);
  a = NA_PTR_TYPE(rb_a, complex*);
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
  if (NA_TYPE(rb_b) != NA_SCOMPLEX)
    rb_b = na_change_type(rb_b, NA_SCOMPLEX);
  b = NA_PTR_TYPE(rb_b, complex*);
  tolb = (real)NUM2DBL(rb_tolb);
  jobu = StringValueCStr(rb_jobu)[0];
  jobv = StringValueCStr(rb_jobv)[0];
  tola = (real)NUM2DBL(rb_tola);
  ldu = lsame_(&jobu,"U") ? MAX(1,m) : 1;
  ldv = lsame_(&jobv,"V") ? MAX(1,p) : 1;
  ldq = lsame_(&jobq,"Q") ? MAX(1,n) : 1;
  lda = m;
  ldb = p;
  {
    int shape[2];
    shape[0] = ldu;
    shape[1] = m;
    rb_u = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  u = NA_PTR_TYPE(rb_u, complex*);
  {
    int shape[2];
    shape[0] = ldv;
    shape[1] = p;
    rb_v = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  v = NA_PTR_TYPE(rb_v, complex*);
  {
    int shape[2];
    shape[0] = ldq;
    shape[1] = n;
    rb_q = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  q = NA_PTR_TYPE(rb_q, complex*);
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
  iwork = ALLOC_N(integer, (n));
  rwork = ALLOC_N(real, (2*n));
  tau = ALLOC_N(complex, (n));
  work = ALLOC_N(complex, (MAX(3*n,m)*(p)));

  cggsvp_(&jobu, &jobv, &jobq, &m, &p, &n, a, &lda, b, &ldb, &tola, &tolb, &k, &l, u, &ldu, v, &ldv, q, &ldq, iwork, rwork, tau, work, &info);

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
init_lapack_cggsvp(VALUE mLapack){
  rb_define_module_function(mLapack, "cggsvp", rb_cggsvp, -1);
}
