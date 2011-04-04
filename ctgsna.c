#include "rb_lapack.h"

extern VOID ctgsna_(char *job, char *howmny, logical *select, integer *n, complex *a, integer *lda, complex *b, integer *ldb, complex *vl, integer *ldvl, complex *vr, integer *ldvr, real *s, real *dif, integer *mm, integer *m, complex *work, integer *lwork, integer *iwork, integer *info);

static VALUE
rb_ctgsna(int argc, VALUE *argv, VALUE self){
  VALUE rb_job;
  char job; 
  VALUE rb_howmny;
  char howmny; 
  VALUE rb_select;
  logical *select; 
  VALUE rb_a;
  complex *a; 
  VALUE rb_b;
  complex *b; 
  VALUE rb_vl;
  complex *vl; 
  VALUE rb_vr;
  complex *vr; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_s;
  real *s; 
  VALUE rb_dif;
  real *dif; 
  VALUE rb_m;
  integer m; 
  VALUE rb_work;
  complex *work; 
  VALUE rb_info;
  integer info; 
  integer *iwork;

  integer n;
  integer lda;
  integer ldb;
  integer ldvl;
  integer ldvr;
  integer mm;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  s, dif, m, work, info = NumRu::Lapack.ctgsna( job, howmny, select, a, b, vl, vr, lwork)\n    or\n  NumRu::Lapack.ctgsna  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 8)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 8)", argc);
  rb_job = argv[0];
  rb_howmny = argv[1];
  rb_select = argv[2];
  rb_a = argv[3];
  rb_b = argv[4];
  rb_vl = argv[5];
  rb_vr = argv[6];
  rb_lwork = argv[7];

  howmny = StringValueCStr(rb_howmny)[0];
  if (!NA_IsNArray(rb_vl))
    rb_raise(rb_eArgError, "vl (6th argument) must be NArray");
  if (NA_RANK(rb_vl) != 2)
    rb_raise(rb_eArgError, "rank of vl (6th argument) must be %d", 2);
  m = NA_SHAPE1(rb_vl);
  ldvl = NA_SHAPE0(rb_vl);
  if (NA_TYPE(rb_vl) != NA_SCOMPLEX)
    rb_vl = na_change_type(rb_vl, NA_SCOMPLEX);
  vl = NA_PTR_TYPE(rb_vl, complex*);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (4th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (4th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
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
  if (NA_TYPE(rb_b) != NA_SCOMPLEX)
    rb_b = na_change_type(rb_b, NA_SCOMPLEX);
  b = NA_PTR_TYPE(rb_b, complex*);
  if (!NA_IsNArray(rb_select))
    rb_raise(rb_eArgError, "select (3th argument) must be NArray");
  if (NA_RANK(rb_select) != 1)
    rb_raise(rb_eArgError, "rank of select (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_select) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of select must be the same as shape 1 of a");
  if (NA_TYPE(rb_select) != NA_LINT)
    rb_select = na_change_type(rb_select, NA_LINT);
  select = NA_PTR_TYPE(rb_select, logical*);
  if (!NA_IsNArray(rb_vr))
    rb_raise(rb_eArgError, "vr (7th argument) must be NArray");
  if (NA_RANK(rb_vr) != 2)
    rb_raise(rb_eArgError, "rank of vr (7th argument) must be %d", 2);
  if (NA_SHAPE1(rb_vr) != m)
    rb_raise(rb_eRuntimeError, "shape 1 of vr must be the same as shape 1 of vl");
  ldvr = NA_SHAPE0(rb_vr);
  if (NA_TYPE(rb_vr) != NA_SCOMPLEX)
    rb_vr = na_change_type(rb_vr, NA_SCOMPLEX);
  vr = NA_PTR_TYPE(rb_vr, complex*);
  job = StringValueCStr(rb_job)[0];
  lwork = NUM2INT(rb_lwork);
  mm = m;
  {
    int shape[1];
    shape[0] = mm;
    rb_s = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  s = NA_PTR_TYPE(rb_s, real*);
  {
    int shape[1];
    shape[0] = mm;
    rb_dif = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  dif = NA_PTR_TYPE(rb_dif, real*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rb_work = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rb_work, complex*);
  iwork = ALLOC_N(integer, (lsame_(&job,"E") ? 0 : n+2));

  ctgsna_(&job, &howmny, select, &n, a, &lda, b, &ldb, vl, &ldvl, vr, &ldvr, s, dif, &mm, &m, work, &lwork, iwork, &info);

  free(iwork);
  rb_m = INT2NUM(m);
  rb_info = INT2NUM(info);
  return rb_ary_new3(5, rb_s, rb_dif, rb_m, rb_work, rb_info);
}

void
init_lapack_ctgsna(VALUE mLapack){
  rb_define_module_function(mLapack, "ctgsna", rb_ctgsna, -1);
}
