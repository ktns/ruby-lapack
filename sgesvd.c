#include "rb_lapack.h"

extern VOID sgesvd_(char *jobu, char *jobvt, integer *m, integer *n, real *a, integer *lda, real *s, real *u, integer *ldu, real *vt, integer *ldvt, real *work, integer *lwork, integer *info);

static VALUE
rb_sgesvd(int argc, VALUE *argv, VALUE self){
  VALUE rb_jobu;
  char jobu; 
  VALUE rb_jobvt;
  char jobvt; 
  VALUE rb_m;
  integer m; 
  VALUE rb_a;
  real *a; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_s;
  real *s; 
  VALUE rb_u;
  real *u; 
  VALUE rb_vt;
  real *vt; 
  VALUE rb_work;
  real *work; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  real *a_out__;

  integer lda;
  integer n;
  integer ldu;
  integer ldvt;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  s, u, vt, work, info, a = NumRu::Lapack.sgesvd( jobu, jobvt, m, a, lwork)\n    or\n  NumRu::Lapack.sgesvd  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_jobu = argv[0];
  rb_jobvt = argv[1];
  rb_m = argv[2];
  rb_a = argv[3];
  rb_lwork = argv[4];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (4th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (4th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);
  m = NUM2INT(rb_m);
  jobvt = StringValueCStr(rb_jobvt)[0];
  lwork = NUM2INT(rb_lwork);
  jobu = StringValueCStr(rb_jobu)[0];
  ldvt = lsame_(&jobvt,"A") ? n : lsame_(&jobvt,"S") ? MIN(m,n) : 1;
  ldu = ((lsame_(&jobu,"S")) || (lsame_(&jobu,"A"))) ? m : 1;
  {
    int shape[1];
    shape[0] = MIN(m,n);
    rb_s = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  s = NA_PTR_TYPE(rb_s, real*);
  {
    int shape[2];
    shape[0] = ldu;
    shape[1] = lsame_(&jobu,"A") ? m : lsame_(&jobu,"S") ? MIN(m,n) : 0;
    rb_u = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  u = NA_PTR_TYPE(rb_u, real*);
  {
    int shape[2];
    shape[0] = ldvt;
    shape[1] = n;
    rb_vt = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  vt = NA_PTR_TYPE(rb_vt, real*);
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

  sgesvd_(&jobu, &jobvt, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(6, rb_s, rb_u, rb_vt, rb_work, rb_info, rb_a);
}

void
init_lapack_sgesvd(VALUE mLapack){
  rb_define_module_function(mLapack, "sgesvd", rb_sgesvd, -1);
}
