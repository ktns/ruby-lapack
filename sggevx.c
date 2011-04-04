#include "rb_lapack.h"

extern VOID sggevx_(char *balanc, char *jobvl, char *jobvr, char *sense, integer *n, real *a, integer *lda, real *b, integer *ldb, real *alphar, real *alphai, real *beta, real *vl, integer *ldvl, real *vr, integer *ldvr, integer *ilo, integer *ihi, real *lscale, real *rscale, real *abnrm, real *bbnrm, real *rconde, real *rcondv, real *work, integer *lwork, integer *iwork, logical *bwork, integer *info);

static VALUE
rb_sggevx(int argc, VALUE *argv, VALUE self){
  VALUE rb_balanc;
  char balanc; 
  VALUE rb_jobvl;
  char jobvl; 
  VALUE rb_jobvr;
  char jobvr; 
  VALUE rb_sense;
  char sense; 
  VALUE rb_a;
  real *a; 
  VALUE rb_b;
  real *b; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_alphar;
  real *alphar; 
  VALUE rb_alphai;
  real *alphai; 
  VALUE rb_beta;
  real *beta; 
  VALUE rb_vl;
  real *vl; 
  VALUE rb_vr;
  real *vr; 
  VALUE rb_ilo;
  integer ilo; 
  VALUE rb_ihi;
  integer ihi; 
  VALUE rb_lscale;
  real *lscale; 
  VALUE rb_rscale;
  real *rscale; 
  VALUE rb_abnrm;
  real abnrm; 
  VALUE rb_bbnrm;
  real bbnrm; 
  VALUE rb_rconde;
  real *rconde; 
  VALUE rb_rcondv;
  real *rcondv; 
  VALUE rb_work;
  real *work; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  real *a_out__;
  VALUE rb_b_out__;
  real *b_out__;
  integer *iwork;
  logical *bwork;

  integer lda;
  integer n;
  integer ldb;
  integer ldvl;
  integer ldvr;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  alphar, alphai, beta, vl, vr, ilo, ihi, lscale, rscale, abnrm, bbnrm, rconde, rcondv, work, info, a, b = NumRu::Lapack.sggevx( balanc, jobvl, jobvr, sense, a, b, lwork)\n    or\n  NumRu::Lapack.sggevx  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 7)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 7)", argc);
  rb_balanc = argv[0];
  rb_jobvl = argv[1];
  rb_jobvr = argv[2];
  rb_sense = argv[3];
  rb_a = argv[4];
  rb_b = argv[5];
  rb_lwork = argv[6];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (5th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (5th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);
  jobvr = StringValueCStr(rb_jobvr)[0];
  balanc = StringValueCStr(rb_balanc)[0];
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (6th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (6th argument) must be %d", 2);
  if (NA_SHAPE1(rb_b) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of b must be the same as shape 1 of a");
  ldb = NA_SHAPE0(rb_b);
  if (NA_TYPE(rb_b) != NA_SFLOAT)
    rb_b = na_change_type(rb_b, NA_SFLOAT);
  b = NA_PTR_TYPE(rb_b, real*);
  jobvl = StringValueCStr(rb_jobvl)[0];
  lwork = NUM2INT(rb_lwork);
  sense = StringValueCStr(rb_sense)[0];
  ldvl = lsame_(&jobvl,"V") ? n : 1;
  ldvr = lsame_(&jobvr,"V") ? n : 1;
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
    shape[0] = n;
    rb_lscale = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  lscale = NA_PTR_TYPE(rb_lscale, real*);
  {
    int shape[1];
    shape[0] = n;
    rb_rscale = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  rscale = NA_PTR_TYPE(rb_rscale, real*);
  {
    int shape[1];
    shape[0] = n;
    rb_rconde = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  rconde = NA_PTR_TYPE(rb_rconde, real*);
  {
    int shape[1];
    shape[0] = n;
    rb_rcondv = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  rcondv = NA_PTR_TYPE(rb_rcondv, real*);
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
  iwork = ALLOC_N(integer, (lsame_(&sense,"E") ? 0 : n+6));
  bwork = ALLOC_N(logical, (lsame_(&sense,"N") ? 0 : n));

  sggevx_(&balanc, &jobvl, &jobvr, &sense, &n, a, &lda, b, &ldb, alphar, alphai, beta, vl, &ldvl, vr, &ldvr, &ilo, &ihi, lscale, rscale, &abnrm, &bbnrm, rconde, rcondv, work, &lwork, iwork, bwork, &info);

  free(iwork);
  free(bwork);
  rb_ilo = INT2NUM(ilo);
  rb_ihi = INT2NUM(ihi);
  rb_abnrm = rb_float_new((double)abnrm);
  rb_bbnrm = rb_float_new((double)bbnrm);
  rb_info = INT2NUM(info);
  return rb_ary_new3(17, rb_alphar, rb_alphai, rb_beta, rb_vl, rb_vr, rb_ilo, rb_ihi, rb_lscale, rb_rscale, rb_abnrm, rb_bbnrm, rb_rconde, rb_rcondv, rb_work, rb_info, rb_a, rb_b);
}

void
init_lapack_sggevx(VALUE mLapack){
  rb_define_module_function(mLapack, "sggevx", rb_sggevx, -1);
}
