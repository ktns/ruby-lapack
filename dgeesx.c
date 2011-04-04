#include "rb_lapack.h"

static logical
rb_select(doublereal *arg0, doublereal *arg1){
  VALUE rb_arg0, rb_arg1;

  VALUE rb_ret;
  logical ret;

  rb_arg0 = rb_float_new((double)(*arg0));
  rb_arg1 = rb_float_new((double)(*arg1));

  rb_ret = rb_yield_values(2, rb_arg0, rb_arg1);

  ret = (rb_ret == Qtrue);
  return ret;
}

extern VOID dgeesx_(char *jobvs, char *sort, L_fp *select, char *sense, integer *n, doublereal *a, integer *lda, integer *sdim, doublereal *wr, doublereal *wi, doublereal *vs, integer *ldvs, doublereal *rconde, doublereal *rcondv, doublereal *work, integer *lwork, integer *iwork, integer *liwork, logical *bwork, integer *info);

static VALUE
rb_dgeesx(int argc, VALUE *argv, VALUE self){
  VALUE rb_jobvs;
  char jobvs; 
  VALUE rb_sort;
  char sort; 
  VALUE rb_sense;
  char sense; 
  VALUE rb_a;
  doublereal *a; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_liwork;
  integer liwork; 
  VALUE rb_sdim;
  integer sdim; 
  VALUE rb_wr;
  doublereal *wr; 
  VALUE rb_wi;
  doublereal *wi; 
  VALUE rb_vs;
  doublereal *vs; 
  VALUE rb_rconde;
  doublereal rconde; 
  VALUE rb_rcondv;
  doublereal rcondv; 
  VALUE rb_work;
  doublereal *work; 
  VALUE rb_iwork;
  integer *iwork; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  doublereal *a_out__;
  logical *bwork;

  integer lda;
  integer n;
  integer ldvs;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  sdim, wr, wi, vs, rconde, rcondv, work, iwork, info, a = NumRu::Lapack.dgeesx( jobvs, sort, sense, a, lwork, liwork){|a,b| ... }\n    or\n  NumRu::Lapack.dgeesx  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rb_jobvs = argv[0];
  rb_sort = argv[1];
  rb_sense = argv[2];
  rb_a = argv[3];
  rb_lwork = argv[4];
  rb_liwork = argv[5];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (4th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (4th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_DFLOAT)
    rb_a = na_change_type(rb_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rb_a, doublereal*);
  jobvs = StringValueCStr(rb_jobvs)[0];
  sort = StringValueCStr(rb_sort)[0];
  lwork = NUM2INT(rb_lwork);
  sense = StringValueCStr(rb_sense)[0];
  liwork = NUM2INT(rb_liwork);
  ldvs = lsame_(&jobvs,"V") ? n : 1;
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
    int shape[2];
    shape[0] = ldvs;
    shape[1] = n;
    rb_vs = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  vs = NA_PTR_TYPE(rb_vs, doublereal*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rb_work = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rb_work, doublereal*);
  {
    int shape[1];
    shape[0] = MAX(1,liwork);
    rb_iwork = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  iwork = NA_PTR_TYPE(rb_iwork, integer*);
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
  bwork = ALLOC_N(logical, (lsame_(&sort,"N") ? 0 : n));

  dgeesx_(&jobvs, &sort, rb_select, &sense, &n, a, &lda, &sdim, wr, wi, vs, &ldvs, &rconde, &rcondv, work, &lwork, iwork, &liwork, bwork, &info);

  free(bwork);
  rb_sdim = INT2NUM(sdim);
  rb_rconde = rb_float_new((double)rconde);
  rb_rcondv = rb_float_new((double)rcondv);
  rb_info = INT2NUM(info);
  return rb_ary_new3(10, rb_sdim, rb_wr, rb_wi, rb_vs, rb_rconde, rb_rcondv, rb_work, rb_iwork, rb_info, rb_a);
}

void
init_lapack_dgeesx(VALUE mLapack){
  rb_define_module_function(mLapack, "dgeesx", rb_dgeesx, -1);
}
