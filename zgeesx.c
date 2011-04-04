#include "rb_lapack.h"

static logical
rb_select(doublecomplex *arg0){
  VALUE rb_arg0;

  VALUE rb_ret;
  logical ret;

  rb_arg0 = rb_funcall(rb_gv_get("Complex"), rb_intern("new"), 2, rb_float_new((double)(arg0->r)), rb_float_new((double)(arg0->i)));

  rb_ret = rb_yield_values(1, rb_arg0);

  ret = (rb_ret == Qtrue);
  return ret;
}

extern VOID zgeesx_(char *jobvs, char *sort, L_fp *select, char *sense, integer *n, doublecomplex *a, integer *lda, integer *sdim, doublecomplex *w, doublecomplex *vs, integer *ldvs, doublereal *rconde, doublereal *rcondv, doublecomplex *work, integer *lwork, doublereal *rwork, logical *bwork, integer *info);

static VALUE
rb_zgeesx(int argc, VALUE *argv, VALUE self){
  VALUE rb_jobvs;
  char jobvs; 
  VALUE rb_sort;
  char sort; 
  VALUE rb_sense;
  char sense; 
  VALUE rb_a;
  doublecomplex *a; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_sdim;
  integer sdim; 
  VALUE rb_w;
  doublecomplex *w; 
  VALUE rb_vs;
  doublecomplex *vs; 
  VALUE rb_rconde;
  doublereal rconde; 
  VALUE rb_rcondv;
  doublereal rcondv; 
  VALUE rb_work;
  doublecomplex *work; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  doublecomplex *a_out__;
  doublereal *rwork;
  logical *bwork;

  integer lda;
  integer n;
  integer ldvs;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  sdim, w, vs, rconde, rcondv, work, info, a = NumRu::Lapack.zgeesx( jobvs, sort, sense, a, lwork){|a| ... }\n    or\n  NumRu::Lapack.zgeesx  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_jobvs = argv[0];
  rb_sort = argv[1];
  rb_sense = argv[2];
  rb_a = argv[3];
  rb_lwork = argv[4];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (4th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (4th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_DCOMPLEX)
    rb_a = na_change_type(rb_a, NA_DCOMPLEX);
  a = NA_PTR_TYPE(rb_a, doublecomplex*);
  jobvs = StringValueCStr(rb_jobvs)[0];
  sort = StringValueCStr(rb_sort)[0];
  lwork = NUM2INT(rb_lwork);
  sense = StringValueCStr(rb_sense)[0];
  ldvs = lsame_(&jobvs,"V") ? n : 1;
  {
    int shape[1];
    shape[0] = n;
    rb_w = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  w = NA_PTR_TYPE(rb_w, doublecomplex*);
  {
    int shape[2];
    shape[0] = ldvs;
    shape[1] = n;
    rb_vs = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  vs = NA_PTR_TYPE(rb_vs, doublecomplex*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rb_work = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rb_work, doublecomplex*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rb_a_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rb_a_out__, doublecomplex*);
  MEMCPY(a_out__, a, doublecomplex, NA_TOTAL(rb_a));
  rb_a = rb_a_out__;
  a = a_out__;
  rwork = ALLOC_N(doublereal, (n));
  bwork = ALLOC_N(logical, (lsame_(&sort,"N") ? 0 : n));

  zgeesx_(&jobvs, &sort, rb_select, &sense, &n, a, &lda, &sdim, w, vs, &ldvs, &rconde, &rcondv, work, &lwork, rwork, bwork, &info);

  free(rwork);
  free(bwork);
  rb_sdim = INT2NUM(sdim);
  rb_rconde = rb_float_new((double)rconde);
  rb_rcondv = rb_float_new((double)rcondv);
  rb_info = INT2NUM(info);
  return rb_ary_new3(8, rb_sdim, rb_w, rb_vs, rb_rconde, rb_rcondv, rb_work, rb_info, rb_a);
}

void
init_lapack_zgeesx(VALUE mLapack){
  rb_define_module_function(mLapack, "zgeesx", rb_zgeesx, -1);
}
