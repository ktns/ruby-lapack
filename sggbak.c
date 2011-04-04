#include "rb_lapack.h"

extern VOID sggbak_(char *job, char *side, integer *n, integer *ilo, integer *ihi, real *lscale, real *rscale, integer *m, real *v, integer *ldv, integer *info);

static VALUE
rb_sggbak(int argc, VALUE *argv, VALUE self){
  VALUE rb_job;
  char job; 
  VALUE rb_side;
  char side; 
  VALUE rb_ilo;
  integer ilo; 
  VALUE rb_ihi;
  integer ihi; 
  VALUE rb_lscale;
  real *lscale; 
  VALUE rb_rscale;
  real *rscale; 
  VALUE rb_v;
  real *v; 
  VALUE rb_info;
  integer info; 
  VALUE rb_v_out__;
  real *v_out__;

  integer n;
  integer ldv;
  integer m;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, v = NumRu::Lapack.sggbak( job, side, ilo, ihi, lscale, rscale, v)\n    or\n  NumRu::Lapack.sggbak  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 7)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 7)", argc);
  rb_job = argv[0];
  rb_side = argv[1];
  rb_ilo = argv[2];
  rb_ihi = argv[3];
  rb_lscale = argv[4];
  rb_rscale = argv[5];
  rb_v = argv[6];

  if (!NA_IsNArray(rb_v))
    rb_raise(rb_eArgError, "v (7th argument) must be NArray");
  if (NA_RANK(rb_v) != 2)
    rb_raise(rb_eArgError, "rank of v (7th argument) must be %d", 2);
  m = NA_SHAPE1(rb_v);
  ldv = NA_SHAPE0(rb_v);
  if (NA_TYPE(rb_v) != NA_SFLOAT)
    rb_v = na_change_type(rb_v, NA_SFLOAT);
  v = NA_PTR_TYPE(rb_v, real*);
  ilo = NUM2INT(rb_ilo);
  side = StringValueCStr(rb_side)[0];
  if (!NA_IsNArray(rb_rscale))
    rb_raise(rb_eArgError, "rscale (6th argument) must be NArray");
  if (NA_RANK(rb_rscale) != 1)
    rb_raise(rb_eArgError, "rank of rscale (6th argument) must be %d", 1);
  n = NA_SHAPE0(rb_rscale);
  if (NA_TYPE(rb_rscale) != NA_SFLOAT)
    rb_rscale = na_change_type(rb_rscale, NA_SFLOAT);
  rscale = NA_PTR_TYPE(rb_rscale, real*);
  job = StringValueCStr(rb_job)[0];
  if (!NA_IsNArray(rb_lscale))
    rb_raise(rb_eArgError, "lscale (5th argument) must be NArray");
  if (NA_RANK(rb_lscale) != 1)
    rb_raise(rb_eArgError, "rank of lscale (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_lscale) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of lscale must be the same as shape 0 of rscale");
  if (NA_TYPE(rb_lscale) != NA_SFLOAT)
    rb_lscale = na_change_type(rb_lscale, NA_SFLOAT);
  lscale = NA_PTR_TYPE(rb_lscale, real*);
  ihi = NUM2INT(rb_ihi);
  {
    int shape[2];
    shape[0] = ldv;
    shape[1] = m;
    rb_v_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  v_out__ = NA_PTR_TYPE(rb_v_out__, real*);
  MEMCPY(v_out__, v, real, NA_TOTAL(rb_v));
  rb_v = rb_v_out__;
  v = v_out__;

  sggbak_(&job, &side, &n, &ilo, &ihi, lscale, rscale, &m, v, &ldv, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_info, rb_v);
}

void
init_lapack_sggbak(VALUE mLapack){
  rb_define_module_function(mLapack, "sggbak", rb_sggbak, -1);
}
