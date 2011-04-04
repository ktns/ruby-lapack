#include "rb_lapack.h"

extern VOID dgesvj_(char *joba, char *jobu, char *jobv, integer *m, integer *n, doublereal *a, integer *lda, doublereal *sva, integer *mv, doublereal *v, integer *ldv, doublereal *work, integer *lwork, integer *info);

static VALUE
rb_dgesvj(int argc, VALUE *argv, VALUE self){
  VALUE rb_joba;
  char joba; 
  VALUE rb_jobu;
  char jobu; 
  VALUE rb_jobv;
  char jobv; 
  VALUE rb_m;
  integer m; 
  VALUE rb_a;
  doublereal *a; 
  VALUE rb_mv;
  integer mv; 
  VALUE rb_v;
  doublereal *v; 
  VALUE rb_work;
  doublereal *work; 
  VALUE rb_sva;
  doublereal *sva; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  doublereal *a_out__;
  VALUE rb_v_out__;
  doublereal *v_out__;
  VALUE rb_work_out__;
  doublereal *work_out__;

  integer lda;
  integer n;
  integer ldv;
  integer lwork;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  sva, info, a, v, work = NumRu::Lapack.dgesvj( joba, jobu, jobv, m, a, mv, v, work)\n    or\n  NumRu::Lapack.dgesvj  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 8)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 8)", argc);
  rb_joba = argv[0];
  rb_jobu = argv[1];
  rb_jobv = argv[2];
  rb_m = argv[3];
  rb_a = argv[4];
  rb_mv = argv[5];
  rb_v = argv[6];
  rb_work = argv[7];

  if (!NA_IsNArray(rb_v))
    rb_raise(rb_eArgError, "v (7th argument) must be NArray");
  if (NA_RANK(rb_v) != 2)
    rb_raise(rb_eArgError, "rank of v (7th argument) must be %d", 2);
  n = NA_SHAPE1(rb_v);
  ldv = NA_SHAPE0(rb_v);
  if (NA_TYPE(rb_v) != NA_DFLOAT)
    rb_v = na_change_type(rb_v, NA_DFLOAT);
  v = NA_PTR_TYPE(rb_v, doublereal*);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (5th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (5th argument) must be %d", 2);
  if (NA_SHAPE1(rb_a) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of a must be the same as shape 1 of v");
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_DFLOAT)
    rb_a = na_change_type(rb_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rb_a, doublereal*);
  m = NUM2INT(rb_m);
  jobu = StringValueCStr(rb_jobu)[0];
  mv = NUM2INT(rb_mv);
  jobv = StringValueCStr(rb_jobv)[0];
  joba = StringValueCStr(rb_joba)[0];
  if (!NA_IsNArray(rb_work))
    rb_raise(rb_eArgError, "work (8th argument) must be NArray");
  if (NA_RANK(rb_work) != 1)
    rb_raise(rb_eArgError, "rank of work (8th argument) must be %d", 1);
  lwork = NA_SHAPE0(rb_work);
  if (lwork != (MAX(4,m+n)))
    rb_raise(rb_eRuntimeError, "shape 0 of work must be %d", MAX(4,m+n));
  if (NA_TYPE(rb_work) != NA_DFLOAT)
    rb_work = na_change_type(rb_work, NA_DFLOAT);
  work = NA_PTR_TYPE(rb_work, doublereal*);
  lwork = MAX(4,m+n);
  {
    int shape[1];
    shape[0] = n;
    rb_sva = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  sva = NA_PTR_TYPE(rb_sva, doublereal*);
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
    shape[0] = ldv;
    shape[1] = n;
    rb_v_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  v_out__ = NA_PTR_TYPE(rb_v_out__, doublereal*);
  MEMCPY(v_out__, v, doublereal, NA_TOTAL(rb_v));
  rb_v = rb_v_out__;
  v = v_out__;
  {
    int shape[1];
    shape[0] = lwork;
    rb_work_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  work_out__ = NA_PTR_TYPE(rb_work_out__, doublereal*);
  MEMCPY(work_out__, work, doublereal, NA_TOTAL(rb_work));
  rb_work = rb_work_out__;
  work = work_out__;

  dgesvj_(&joba, &jobu, &jobv, &m, &n, a, &lda, sva, &mv, v, &ldv, work, &lwork, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(5, rb_sva, rb_info, rb_a, rb_v, rb_work);
}

void
init_lapack_dgesvj(VALUE mLapack){
  rb_define_module_function(mLapack, "dgesvj", rb_dgesvj, -1);
}
