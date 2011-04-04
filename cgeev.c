#include "rb_lapack.h"

extern VOID cgeev_(char *jobvl, char *jobvr, integer *n, complex *a, integer *lda, complex *w, complex *vl, integer *ldvl, complex *vr, integer *ldvr, complex *work, integer *lwork, real *rwork, integer *info);

static VALUE
rb_cgeev(int argc, VALUE *argv, VALUE self){
  VALUE rb_jobvl;
  char jobvl; 
  VALUE rb_jobvr;
  char jobvr; 
  VALUE rb_a;
  complex *a; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_w;
  complex *w; 
  VALUE rb_vl;
  complex *vl; 
  VALUE rb_vr;
  complex *vr; 
  VALUE rb_work;
  complex *work; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  complex *a_out__;
  real *rwork;

  integer lda;
  integer n;
  integer ldvl;
  integer ldvr;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  w, vl, vr, work, info, a = NumRu::Lapack.cgeev( jobvl, jobvr, a, lwork)\n    or\n  NumRu::Lapack.cgeev  # print help\n\n\nFORTRAN MANUAL\n\n");
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
  if (NA_TYPE(rb_a) != NA_SCOMPLEX)
    rb_a = na_change_type(rb_a, NA_SCOMPLEX);
  a = NA_PTR_TYPE(rb_a, complex*);
  jobvr = StringValueCStr(rb_jobvr)[0];
  jobvl = StringValueCStr(rb_jobvl)[0];
  lwork = NUM2INT(rb_lwork);
  ldvr = lsame_(&jobvr,"V") ? n : 1;
  ldvl = lsame_(&jobvl,"V") ? n : 1;
  {
    int shape[1];
    shape[0] = n;
    rb_w = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  w = NA_PTR_TYPE(rb_w, complex*);
  {
    int shape[2];
    shape[0] = ldvl;
    shape[1] = n;
    rb_vl = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  vl = NA_PTR_TYPE(rb_vl, complex*);
  {
    int shape[2];
    shape[0] = ldvr;
    shape[1] = n;
    rb_vr = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  vr = NA_PTR_TYPE(rb_vr, complex*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rb_work = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rb_work, complex*);
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
  rwork = ALLOC_N(real, (2*n));

  cgeev_(&jobvl, &jobvr, &n, a, &lda, w, vl, &ldvl, vr, &ldvr, work, &lwork, rwork, &info);

  free(rwork);
  rb_info = INT2NUM(info);
  return rb_ary_new3(6, rb_w, rb_vl, rb_vr, rb_work, rb_info, rb_a);
}

void
init_lapack_cgeev(VALUE mLapack){
  rb_define_module_function(mLapack, "cgeev", rb_cgeev, -1);
}
