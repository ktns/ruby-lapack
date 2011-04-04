#include "rb_lapack.h"

extern VOID dtrsen_(char *job, char *compq, logical *select, integer *n, doublereal *t, integer *ldt, doublereal *q, integer *ldq, doublereal *wr, doublereal *wi, integer *m, doublereal *s, doublereal *sep, doublereal *work, integer *lwork, integer *iwork, integer *liwork, integer *info);

static VALUE
rb_dtrsen(int argc, VALUE *argv, VALUE self){
  VALUE rb_job;
  char job; 
  VALUE rb_compq;
  char compq; 
  VALUE rb_select;
  logical *select; 
  VALUE rb_t;
  doublereal *t; 
  VALUE rb_q;
  doublereal *q; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_liwork;
  integer liwork; 
  VALUE rb_wr;
  doublereal *wr; 
  VALUE rb_wi;
  doublereal *wi; 
  VALUE rb_m;
  integer m; 
  VALUE rb_s;
  doublereal s; 
  VALUE rb_sep;
  doublereal sep; 
  VALUE rb_work;
  doublereal *work; 
  VALUE rb_info;
  integer info; 
  VALUE rb_t_out__;
  doublereal *t_out__;
  VALUE rb_q_out__;
  doublereal *q_out__;
  integer *iwork;

  integer n;
  integer ldt;
  integer ldq;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  wr, wi, m, s, sep, work, info, t, q = NumRu::Lapack.dtrsen( job, compq, select, t, q, lwork, liwork)\n    or\n  NumRu::Lapack.dtrsen  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 7)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 7)", argc);
  rb_job = argv[0];
  rb_compq = argv[1];
  rb_select = argv[2];
  rb_t = argv[3];
  rb_q = argv[4];
  rb_lwork = argv[5];
  rb_liwork = argv[6];

  liwork = NUM2INT(rb_liwork);
  compq = StringValueCStr(rb_compq)[0];
  lwork = NUM2INT(rb_lwork);
  if (!NA_IsNArray(rb_q))
    rb_raise(rb_eArgError, "q (5th argument) must be NArray");
  if (NA_RANK(rb_q) != 2)
    rb_raise(rb_eArgError, "rank of q (5th argument) must be %d", 2);
  n = NA_SHAPE1(rb_q);
  ldq = NA_SHAPE0(rb_q);
  if (NA_TYPE(rb_q) != NA_DFLOAT)
    rb_q = na_change_type(rb_q, NA_DFLOAT);
  q = NA_PTR_TYPE(rb_q, doublereal*);
  job = StringValueCStr(rb_job)[0];
  if (!NA_IsNArray(rb_select))
    rb_raise(rb_eArgError, "select (3th argument) must be NArray");
  if (NA_RANK(rb_select) != 1)
    rb_raise(rb_eArgError, "rank of select (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_select) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of select must be the same as shape 1 of q");
  if (NA_TYPE(rb_select) != NA_LINT)
    rb_select = na_change_type(rb_select, NA_LINT);
  select = NA_PTR_TYPE(rb_select, logical*);
  if (!NA_IsNArray(rb_t))
    rb_raise(rb_eArgError, "t (4th argument) must be NArray");
  if (NA_RANK(rb_t) != 2)
    rb_raise(rb_eArgError, "rank of t (4th argument) must be %d", 2);
  if (NA_SHAPE1(rb_t) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of t must be the same as shape 1 of q");
  ldt = NA_SHAPE0(rb_t);
  if (NA_TYPE(rb_t) != NA_DFLOAT)
    rb_t = na_change_type(rb_t, NA_DFLOAT);
  t = NA_PTR_TYPE(rb_t, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_wr = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  wr = NA_PTR_TYPE(rb_wr, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_wi = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  wi = NA_PTR_TYPE(rb_wi, doublereal*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rb_work = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rb_work, doublereal*);
  {
    int shape[2];
    shape[0] = ldt;
    shape[1] = n;
    rb_t_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  t_out__ = NA_PTR_TYPE(rb_t_out__, doublereal*);
  MEMCPY(t_out__, t, doublereal, NA_TOTAL(rb_t));
  rb_t = rb_t_out__;
  t = t_out__;
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
  iwork = ALLOC_N(integer, (MAX(1,liwork)));

  dtrsen_(&job, &compq, select, &n, t, &ldt, q, &ldq, wr, wi, &m, &s, &sep, work, &lwork, iwork, &liwork, &info);

  free(iwork);
  rb_m = INT2NUM(m);
  rb_s = rb_float_new((double)s);
  rb_sep = rb_float_new((double)sep);
  rb_info = INT2NUM(info);
  return rb_ary_new3(9, rb_wr, rb_wi, rb_m, rb_s, rb_sep, rb_work, rb_info, rb_t, rb_q);
}

void
init_lapack_dtrsen(VALUE mLapack){
  rb_define_module_function(mLapack, "dtrsen", rb_dtrsen, -1);
}
