#include "rb_lapack.h"

extern VOID sgejsv_(char *joba, char *jobu, char *jobv, char *jobr, char *jobt, char *jobp, integer *m, integer *n, real *a, integer *lda, real *sva, real *u, integer *ldu, real *v, integer *ldv, real *work, integer *lwork, integer *iwork, integer *info);

static VALUE
rb_sgejsv(int argc, VALUE *argv, VALUE self){
  VALUE rb_joba;
  char joba; 
  VALUE rb_jobu;
  char jobu; 
  VALUE rb_jobv;
  char jobv; 
  VALUE rb_jobr;
  char jobr; 
  VALUE rb_jobt;
  char jobt; 
  VALUE rb_jobp;
  char jobp; 
  VALUE rb_m;
  integer m; 
  VALUE rb_a;
  real *a; 
  VALUE rb_work;
  real *work; 
  VALUE rb_sva;
  real *sva; 
  VALUE rb_u;
  real *u; 
  VALUE rb_v;
  real *v; 
  VALUE rb_iwork;
  integer *iwork; 
  VALUE rb_info;
  integer info; 
  VALUE rb_work_out__;
  real *work_out__;

  integer lda;
  integer n;
  integer lwork;
  integer ldu;
  integer ldv;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  sva, u, v, iwork, info, work = NumRu::Lapack.sgejsv( joba, jobu, jobv, jobr, jobt, jobp, m, a, work)\n    or\n  NumRu::Lapack.sgejsv  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 9)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 9)", argc);
  rb_joba = argv[0];
  rb_jobu = argv[1];
  rb_jobv = argv[2];
  rb_jobr = argv[3];
  rb_jobt = argv[4];
  rb_jobp = argv[5];
  rb_m = argv[6];
  rb_a = argv[7];
  rb_work = argv[8];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (8th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (8th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);
  jobr = StringValueCStr(rb_jobr)[0];
  m = NUM2INT(rb_m);
  jobt = StringValueCStr(rb_jobt)[0];
  jobu = StringValueCStr(rb_jobu)[0];
  if (!NA_IsNArray(rb_work))
    rb_raise(rb_eArgError, "work (9th argument) must be NArray");
  if (NA_RANK(rb_work) != 1)
    rb_raise(rb_eArgError, "rank of work (9th argument) must be %d", 1);
  lwork = NA_SHAPE0(rb_work);
  if (NA_TYPE(rb_work) != NA_SFLOAT)
    rb_work = na_change_type(rb_work, NA_SFLOAT);
  work = NA_PTR_TYPE(rb_work, real*);
  joba = StringValueCStr(rb_joba)[0];
  jobp = StringValueCStr(rb_jobp)[0];
  jobv = StringValueCStr(rb_jobv)[0];
  ldv = (lsame_(&jobu,"U")||lsame_(&jobu,"F")||lsame_(&jobu,"W")) ? n : 1;
  ldu = (lsame_(&jobu,"U")||lsame_(&jobu,"F")||lsame_(&jobu,"W")) ? m : 1;
  {
    int shape[1];
    shape[0] = n;
    rb_sva = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  sva = NA_PTR_TYPE(rb_sva, real*);
  {
    int shape[2];
    shape[0] = ldu;
    shape[1] = n;
    rb_u = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  u = NA_PTR_TYPE(rb_u, real*);
  {
    int shape[2];
    shape[0] = ldv;
    shape[1] = n;
    rb_v = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  v = NA_PTR_TYPE(rb_v, real*);
  {
    int shape[1];
    shape[0] = m+3*n;
    rb_iwork = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  iwork = NA_PTR_TYPE(rb_iwork, integer*);
  {
    int shape[1];
    shape[0] = lwork;
    rb_work_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  work_out__ = NA_PTR_TYPE(rb_work_out__, real*);
  MEMCPY(work_out__, work, real, NA_TOTAL(rb_work));
  rb_work = rb_work_out__;
  work = work_out__;

  sgejsv_(&joba, &jobu, &jobv, &jobr, &jobt, &jobp, &m, &n, a, &lda, sva, u, &ldu, v, &ldv, work, &lwork, iwork, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(6, rb_sva, rb_u, rb_v, rb_iwork, rb_info, rb_work);
}

void
init_lapack_sgejsv(VALUE mLapack){
  rb_define_module_function(mLapack, "sgejsv", rb_sgejsv, -1);
}
