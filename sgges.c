#include "rb_lapack.h"

static logical
rb_selctg(real *arg0, real *arg1, real *arg2){
  VALUE rb_arg0, rb_arg1, rb_arg2;

  VALUE rb_ret;
  logical ret;

  rb_arg0 = rb_float_new((double)(*arg0));
  rb_arg1 = rb_float_new((double)(*arg1));
  rb_arg2 = rb_float_new((double)(*arg2));

  rb_ret = rb_yield_values(3, rb_arg0, rb_arg1, rb_arg2);

  ret = (rb_ret == Qtrue);
  return ret;
}

extern VOID sgges_(char *jobvsl, char *jobvsr, char *sort, L_fp *selctg, integer *n, real *a, integer *lda, real *b, integer *ldb, integer *sdim, real *alphar, real *alphai, real *beta, real *vsl, integer *ldvsl, real *vsr, integer *ldvsr, real *work, integer *lwork, logical *bwork, integer *info);

static VALUE
rb_sgges(int argc, VALUE *argv, VALUE self){
  VALUE rb_jobvsl;
  char jobvsl; 
  VALUE rb_jobvsr;
  char jobvsr; 
  VALUE rb_sort;
  char sort; 
  VALUE rb_a;
  real *a; 
  VALUE rb_b;
  real *b; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_sdim;
  integer sdim; 
  VALUE rb_alphar;
  real *alphar; 
  VALUE rb_alphai;
  real *alphai; 
  VALUE rb_beta;
  real *beta; 
  VALUE rb_vsl;
  real *vsl; 
  VALUE rb_vsr;
  real *vsr; 
  VALUE rb_work;
  real *work; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  real *a_out__;
  VALUE rb_b_out__;
  real *b_out__;
  logical *bwork;

  integer lda;
  integer n;
  integer ldb;
  integer ldvsl;
  integer ldvsr;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  sdim, alphar, alphai, beta, vsl, vsr, work, info, a, b = NumRu::Lapack.sgges( jobvsl, jobvsr, sort, a, b, lwork){|a,b,c| ... }\n    or\n  NumRu::Lapack.sgges  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rb_jobvsl = argv[0];
  rb_jobvsr = argv[1];
  rb_sort = argv[2];
  rb_a = argv[3];
  rb_b = argv[4];
  rb_lwork = argv[5];

  jobvsl = StringValueCStr(rb_jobvsl)[0];
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (4th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (4th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (5th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (5th argument) must be %d", 2);
  if (NA_SHAPE1(rb_b) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of b must be the same as shape 1 of a");
  ldb = NA_SHAPE0(rb_b);
  if (NA_TYPE(rb_b) != NA_SFLOAT)
    rb_b = na_change_type(rb_b, NA_SFLOAT);
  b = NA_PTR_TYPE(rb_b, real*);
  jobvsr = StringValueCStr(rb_jobvsr)[0];
  sort = StringValueCStr(rb_sort)[0];
  lwork = NUM2INT(rb_lwork);
  ldvsr = lsame_(&jobvsr,"V") ? n : 1;
  ldvsl = lsame_(&jobvsl,"V") ? n : 1;
  {
    int shape[1];
    shape[0] = n;
    rb_alphar = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  alphar = NA_PTR_TYPE(rb_alphar, real*);
  {
    int shape[1];
    shape[0] = n;
    rb_alphai = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  alphai = NA_PTR_TYPE(rb_alphai, real*);
  {
    int shape[1];
    shape[0] = n;
    rb_beta = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  beta = NA_PTR_TYPE(rb_beta, real*);
  {
    int shape[2];
    shape[0] = ldvsl;
    shape[1] = n;
    rb_vsl = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  vsl = NA_PTR_TYPE(rb_vsl, real*);
  {
    int shape[2];
    shape[0] = ldvsr;
    shape[1] = n;
    rb_vsr = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  vsr = NA_PTR_TYPE(rb_vsr, real*);
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
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = n;
    rb_b_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rb_b_out__, real*);
  MEMCPY(b_out__, b, real, NA_TOTAL(rb_b));
  rb_b = rb_b_out__;
  b = b_out__;
  bwork = ALLOC_N(logical, (lsame_(&sort,"N") ? 0 : n));

  sgges_(&jobvsl, &jobvsr, &sort, rb_selctg, &n, a, &lda, b, &ldb, &sdim, alphar, alphai, beta, vsl, &ldvsl, vsr, &ldvsr, work, &lwork, bwork, &info);

  free(bwork);
  rb_sdim = INT2NUM(sdim);
  rb_info = INT2NUM(info);
  return rb_ary_new3(10, rb_sdim, rb_alphar, rb_alphai, rb_beta, rb_vsl, rb_vsr, rb_work, rb_info, rb_a, rb_b);
}

void
init_lapack_sgges(VALUE mLapack){
  rb_define_module_function(mLapack, "sgges", rb_sgges, -1);
}
