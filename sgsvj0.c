#include "rb_lapack.h"

extern VOID sgsvj0_(char *jobv, integer *m, integer *n, real *a, integer *lda, real *d, real *sva, integer *mv, real *v, integer *ldv, integer *eps, integer *sfmin, real *tol, integer *nsweep, real *work, integer *lwork, integer *info);

static VALUE
rb_sgsvj0(int argc, VALUE *argv, VALUE self){
  VALUE rb_jobv;
  char jobv; 
  VALUE rb_m;
  integer m; 
  VALUE rb_a;
  real *a; 
  VALUE rb_d;
  real *d; 
  VALUE rb_sva;
  real *sva; 
  VALUE rb_mv;
  integer mv; 
  VALUE rb_v;
  real *v; 
  VALUE rb_eps;
  integer eps; 
  VALUE rb_sfmin;
  integer sfmin; 
  VALUE rb_tol;
  real tol; 
  VALUE rb_nsweep;
  integer nsweep; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  real *a_out__;
  VALUE rb_d_out__;
  real *d_out__;
  VALUE rb_sva_out__;
  real *sva_out__;
  VALUE rb_v_out__;
  real *v_out__;
  real *work;

  integer lda;
  integer n;
  integer ldv;
  integer lwork;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, a, d, sva, v = NumRu::Lapack.sgsvj0( jobv, m, a, d, sva, mv, v, eps, sfmin, tol, nsweep)\n    or\n  NumRu::Lapack.sgsvj0  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 11)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 11)", argc);
  rb_jobv = argv[0];
  rb_m = argv[1];
  rb_a = argv[2];
  rb_d = argv[3];
  rb_sva = argv[4];
  rb_mv = argv[5];
  rb_v = argv[6];
  rb_eps = argv[7];
  rb_sfmin = argv[8];
  rb_tol = argv[9];
  rb_nsweep = argv[10];

  sfmin = NUM2INT(rb_sfmin);
  if (!NA_IsNArray(rb_v))
    rb_raise(rb_eArgError, "v (7th argument) must be NArray");
  if (NA_RANK(rb_v) != 2)
    rb_raise(rb_eArgError, "rank of v (7th argument) must be %d", 2);
  n = NA_SHAPE1(rb_v);
  ldv = NA_SHAPE0(rb_v);
  if (NA_TYPE(rb_v) != NA_SFLOAT)
    rb_v = na_change_type(rb_v, NA_SFLOAT);
  v = NA_PTR_TYPE(rb_v, real*);
  nsweep = NUM2INT(rb_nsweep);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  if (NA_SHAPE1(rb_a) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of a must be the same as shape 1 of v");
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);
  m = NUM2INT(rb_m);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (4th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_d) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of d must be the same as shape 1 of v");
  if (NA_TYPE(rb_d) != NA_SFLOAT)
    rb_d = na_change_type(rb_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rb_d, real*);
  if (!NA_IsNArray(rb_sva))
    rb_raise(rb_eArgError, "sva (5th argument) must be NArray");
  if (NA_RANK(rb_sva) != 1)
    rb_raise(rb_eArgError, "rank of sva (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_sva) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of sva must be the same as shape 1 of v");
  if (NA_TYPE(rb_sva) != NA_SFLOAT)
    rb_sva = na_change_type(rb_sva, NA_SFLOAT);
  sva = NA_PTR_TYPE(rb_sva, real*);
  jobv = StringValueCStr(rb_jobv)[0];
  tol = (real)NUM2DBL(rb_tol);
  eps = NUM2INT(rb_eps);
  mv = NUM2INT(rb_mv);
  lwork = m;
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
    shape[0] = n;
    rb_sva_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  sva_out__ = NA_PTR_TYPE(rb_sva_out__, real*);
  MEMCPY(sva_out__, sva, real, NA_TOTAL(rb_sva));
  rb_sva = rb_sva_out__;
  sva = sva_out__;
  {
    int shape[2];
    shape[0] = ldv;
    shape[1] = n;
    rb_v_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  v_out__ = NA_PTR_TYPE(rb_v_out__, real*);
  MEMCPY(v_out__, v, real, NA_TOTAL(rb_v));
  rb_v = rb_v_out__;
  v = v_out__;
  work = ALLOC_N(real, (lwork));

  sgsvj0_(&jobv, &m, &n, a, &lda, d, sva, &mv, v, &ldv, &eps, &sfmin, &tol, &nsweep, work, &lwork, &info);

  free(work);
  rb_info = INT2NUM(info);
  return rb_ary_new3(5, rb_info, rb_a, rb_d, rb_sva, rb_v);
}

void
init_lapack_sgsvj0(VALUE mLapack){
  rb_define_module_function(mLapack, "sgsvj0", rb_sgsvj0, -1);
}
