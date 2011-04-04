#include "rb_lapack.h"

extern VOID sgeev_(char *jobvl, char *jobvr, integer *n, real *a, integer *lda, real *wr, real *wi, real *vl, integer *ldvl, real *vr, integer *ldvr, real *work, integer *lwork, integer *info);

static VALUE
rb_sgeev(int argc, VALUE *argv, VALUE self){
  VALUE rb_jobvl;
  char jobvl; 
  VALUE rb_jobvr;
  char jobvr; 
  VALUE rb_a;
  real *a; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_wr;
  real *wr; 
  VALUE rb_wi;
  real *wi; 
  VALUE rb_vl;
  real *vl; 
  VALUE rb_vr;
  real *vr; 
  VALUE rb_work;
  real *work; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  real *a_out__;

  integer lda;
  integer n;
  integer ldvl;
  integer ldvr;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  wr, wi, vl, vr, work, info, a = NumRu::Lapack.sgeev( jobvl, jobvr, a, lwork)\n    or\n  NumRu::Lapack.sgeev  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_jobvl = argv[0];
  rb_jobvr = argv[1];
  rb_a = argv[2];
  rb_lwork = argv[3];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);
  jobvr = StringValueCStr(rb_jobvr)[0];
  jobvl = StringValueCStr(rb_jobvl)[0];
  lwork = NUM2INT(rb_lwork);
  ldvr = lsame_(&jobvr,"V") ? n : 1;
  ldvl = lsame_(&jobvl,"V") ? n : 1;
  {
    int shape[1];
    shape[0] = n;
    rb_wr = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  wr = NA_PTR_TYPE(rb_wr, real*);
  {
    int shape[1];
    shape[0] = n;
    rb_wi = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  wi = NA_PTR_TYPE(rb_wi, real*);
  {
    int shape[2];
    shape[0] = ldvl;
    shape[1] = n;
    rb_vl = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  vl = NA_PTR_TYPE(rb_vl, real*);
  {
    int shape[2];
    shape[0] = ldvr;
    shape[1] = n;
    rb_vr = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  vr = NA_PTR_TYPE(rb_vr, real*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rb_work = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rb_work, real*);
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

  sgeev_(&jobvl, &jobvr, &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr, work, &lwork, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(7, rb_wr, rb_wi, rb_vl, rb_vr, rb_work, rb_info, rb_a);
}

void
init_lapack_sgeev(VALUE mLapack){
  rb_define_module_function(mLapack, "sgeev", rb_sgeev, -1);
}
