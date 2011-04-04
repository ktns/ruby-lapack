#include "rb_lapack.h"

extern VOID dgsvj1_(char *jobv, integer *m, integer *n, integer *n1, doublereal *a, integer *lda, doublereal *d, doublereal *sva, integer *mv, doublereal *v, integer *ldv, doublereal *eps, doublereal *sfmin, doublereal *tol, integer *nsweep, doublereal *work, integer *lwork, integer *info);

static VALUE
rb_dgsvj1(int argc, VALUE *argv, VALUE self){
  VALUE rb_jobv;
  char jobv; 
  VALUE rb_m;
  integer m; 
  VALUE rb_n1;
  integer n1; 
  VALUE rb_a;
  doublereal *a; 
  VALUE rb_d;
  doublereal *d; 
  VALUE rb_sva;
  doublereal *sva; 
  VALUE rb_mv;
  integer mv; 
  VALUE rb_v;
  doublereal *v; 
  VALUE rb_eps;
  doublereal eps; 
  VALUE rb_sfmin;
  doublereal sfmin; 
  VALUE rb_tol;
  doublereal tol; 
  VALUE rb_nsweep;
  integer nsweep; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  doublereal *a_out__;
  VALUE rb_d_out__;
  doublereal *d_out__;
  VALUE rb_sva_out__;
  doublereal *sva_out__;
  VALUE rb_v_out__;
  doublereal *v_out__;
  doublereal *work;

  integer lda;
  integer n;
  integer ldv;
  integer lwork;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, a, d, sva, v = NumRu::Lapack.dgsvj1( jobv, m, n1, a, d, sva, mv, v, eps, sfmin, tol, nsweep)\n    or\n  NumRu::Lapack.dgsvj1  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 12)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 12)", argc);
  rb_jobv = argv[0];
  rb_m = argv[1];
  rb_n1 = argv[2];
  rb_a = argv[3];
  rb_d = argv[4];
  rb_sva = argv[5];
  rb_mv = argv[6];
  rb_v = argv[7];
  rb_eps = argv[8];
  rb_sfmin = argv[9];
  rb_tol = argv[10];
  rb_nsweep = argv[11];

  sfmin = NUM2DBL(rb_sfmin);
  if (!NA_IsNArray(rb_v))
    rb_raise(rb_eArgError, "v (8th argument) must be NArray");
  if (NA_RANK(rb_v) != 2)
    rb_raise(rb_eArgError, "rank of v (8th argument) must be %d", 2);
  n = NA_SHAPE1(rb_v);
  ldv = NA_SHAPE0(rb_v);
  if (NA_TYPE(rb_v) != NA_DFLOAT)
    rb_v = na_change_type(rb_v, NA_DFLOAT);
  v = NA_PTR_TYPE(rb_v, doublereal*);
  nsweep = NUM2INT(rb_nsweep);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (4th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (4th argument) must be %d", 2);
  if (NA_SHAPE1(rb_a) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of a must be the same as shape 1 of v");
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_DFLOAT)
    rb_a = na_change_type(rb_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rb_a, doublereal*);
  m = NUM2INT(rb_m);
  n1 = NUM2INT(rb_n1);
  if (!NA_IsNArray(rb_sva))
    rb_raise(rb_eArgError, "sva (6th argument) must be NArray");
  if (NA_RANK(rb_sva) != 1)
    rb_raise(rb_eArgError, "rank of sva (6th argument) must be %d", 1);
  if (NA_SHAPE0(rb_sva) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of sva must be the same as shape 1 of v");
  if (NA_TYPE(rb_sva) != NA_DFLOAT)
    rb_sva = na_change_type(rb_sva, NA_DFLOAT);
  sva = NA_PTR_TYPE(rb_sva, doublereal*);
  mv = NUM2INT(rb_mv);
  jobv = StringValueCStr(rb_jobv)[0];
  tol = NUM2DBL(rb_tol);
  eps = NUM2DBL(rb_eps);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (5th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_d) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of d must be the same as shape 1 of v");
  if (NA_TYPE(rb_d) != NA_DFLOAT)
    rb_d = na_change_type(rb_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rb_d, doublereal*);
  lwork = m;
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
    int shape[1];
    shape[0] = n;
    rb_d_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  d_out__ = NA_PTR_TYPE(rb_d_out__, doublereal*);
  MEMCPY(d_out__, d, doublereal, NA_TOTAL(rb_d));
  rb_d = rb_d_out__;
  d = d_out__;
  {
    int shape[1];
    shape[0] = n;
    rb_sva_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  sva_out__ = NA_PTR_TYPE(rb_sva_out__, doublereal*);
  MEMCPY(sva_out__, sva, doublereal, NA_TOTAL(rb_sva));
  rb_sva = rb_sva_out__;
  sva = sva_out__;
  {
    int shape[2];
    shape[0] = ldv;
    shape[1] = n;
    rb_v_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  v_out__ = NA_PTR_TYPE(rb_v_out__, doublereal*);
  MEMCPY(v_out__, v, doublereal, NA_TOTAL(rb_v));
  rb_v = rb_v_out__;
  v = v_out__;
  work = ALLOC_N(doublereal, (lwork));

  dgsvj1_(&jobv, &m, &n, &n1, a, &lda, d, sva, &mv, v, &ldv, &eps, &sfmin, &tol, &nsweep, work, &lwork, &info);

  free(work);
  rb_info = INT2NUM(info);
  return rb_ary_new3(5, rb_info, rb_a, rb_d, rb_sva, rb_v);
}

void
init_lapack_dgsvj1(VALUE mLapack){
  rb_define_module_function(mLapack, "dgsvj1", rb_dgsvj1, -1);
}
