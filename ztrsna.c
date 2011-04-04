#include "rb_lapack.h"

extern VOID ztrsna_(char *job, char *howmny, logical *select, integer *n, doublecomplex *t, integer *ldt, doublecomplex *vl, integer *ldvl, doublecomplex *vr, integer *ldvr, doublereal *s, doublereal *sep, integer *mm, integer *m, doublecomplex *work, integer *ldwork, doublereal *rwork, integer *info);

static VALUE
rb_ztrsna(int argc, VALUE *argv, VALUE self){
  VALUE rb_job;
  char job; 
  VALUE rb_howmny;
  char howmny; 
  VALUE rb_select;
  logical *select; 
  VALUE rb_t;
  doublecomplex *t; 
  VALUE rb_vl;
  doublecomplex *vl; 
  VALUE rb_vr;
  doublecomplex *vr; 
  VALUE rb_ldwork;
  integer ldwork; 
  VALUE rb_s;
  doublereal *s; 
  VALUE rb_sep;
  doublereal *sep; 
  VALUE rb_m;
  integer m; 
  VALUE rb_info;
  integer info; 
  doublecomplex *work;
  doublereal *rwork;

  integer n;
  integer ldt;
  integer ldvl;
  integer ldvr;
  integer mm;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  s, sep, m, info = NumRu::Lapack.ztrsna( job, howmny, select, t, vl, vr, ldwork)\n    or\n  NumRu::Lapack.ztrsna  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 7)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 7)", argc);
  rb_job = argv[0];
  rb_howmny = argv[1];
  rb_select = argv[2];
  rb_t = argv[3];
  rb_vl = argv[4];
  rb_vr = argv[5];
  rb_ldwork = argv[6];

  howmny = StringValueCStr(rb_howmny)[0];
  if (!NA_IsNArray(rb_t))
    rb_raise(rb_eArgError, "t (4th argument) must be NArray");
  if (NA_RANK(rb_t) != 2)
    rb_raise(rb_eArgError, "rank of t (4th argument) must be %d", 2);
  n = NA_SHAPE1(rb_t);
  ldt = NA_SHAPE0(rb_t);
  if (NA_TYPE(rb_t) != NA_DCOMPLEX)
    rb_t = na_change_type(rb_t, NA_DCOMPLEX);
  t = NA_PTR_TYPE(rb_t, doublecomplex*);
  if (!NA_IsNArray(rb_vl))
    rb_raise(rb_eArgError, "vl (5th argument) must be NArray");
  if (NA_RANK(rb_vl) != 2)
    rb_raise(rb_eArgError, "rank of vl (5th argument) must be %d", 2);
  m = NA_SHAPE1(rb_vl);
  ldvl = NA_SHAPE0(rb_vl);
  if (NA_TYPE(rb_vl) != NA_DCOMPLEX)
    rb_vl = na_change_type(rb_vl, NA_DCOMPLEX);
  vl = NA_PTR_TYPE(rb_vl, doublecomplex*);
  if (!NA_IsNArray(rb_vr))
    rb_raise(rb_eArgError, "vr (6th argument) must be NArray");
  if (NA_RANK(rb_vr) != 2)
    rb_raise(rb_eArgError, "rank of vr (6th argument) must be %d", 2);
  if (NA_SHAPE1(rb_vr) != m)
    rb_raise(rb_eRuntimeError, "shape 1 of vr must be the same as shape 1 of vl");
  ldvr = NA_SHAPE0(rb_vr);
  if (NA_TYPE(rb_vr) != NA_DCOMPLEX)
    rb_vr = na_change_type(rb_vr, NA_DCOMPLEX);
  vr = NA_PTR_TYPE(rb_vr, doublecomplex*);
  if (!NA_IsNArray(rb_select))
    rb_raise(rb_eArgError, "select (3th argument) must be NArray");
  if (NA_RANK(rb_select) != 1)
    rb_raise(rb_eArgError, "rank of select (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_select) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of select must be the same as shape 1 of t");
  if (NA_TYPE(rb_select) != NA_LINT)
    rb_select = na_change_type(rb_select, NA_LINT);
  select = NA_PTR_TYPE(rb_select, logical*);
  job = StringValueCStr(rb_job)[0];
  ldwork = ((lsame_(&job,"V")) || (lsame_(&job,"B"))) ? n : 1;
  mm = m;
  {
    int shape[1];
    shape[0] = mm;
    rb_s = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  s = NA_PTR_TYPE(rb_s, doublereal*);
  {
    int shape[1];
    shape[0] = mm;
    rb_sep = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  sep = NA_PTR_TYPE(rb_sep, doublereal*);
  work = ALLOC_N(doublecomplex, (lsame_(&job,"E") ? 0 : ldwork)*(lsame_(&job,"E") ? 0 : n+6));
  rwork = ALLOC_N(doublereal, (lsame_(&job,"E") ? 0 : n));

  ztrsna_(&job, &howmny, select, &n, t, &ldt, vl, &ldvl, vr, &ldvr, s, sep, &mm, &m, work, &ldwork, rwork, &info);

  free(work);
  free(rwork);
  rb_m = INT2NUM(m);
  rb_info = INT2NUM(info);
  return rb_ary_new3(4, rb_s, rb_sep, rb_m, rb_info);
}

void
init_lapack_ztrsna(VALUE mLapack){
  rb_define_module_function(mLapack, "ztrsna", rb_ztrsna, -1);
}
