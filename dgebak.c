#include "rb_lapack.h"

extern VOID dgebak_(char *job, char *side, integer *n, integer *ilo, integer *ihi, doublereal *scale, integer *m, doublereal *v, integer *ldv, integer *info);

static VALUE
rb_dgebak(int argc, VALUE *argv, VALUE self){
  VALUE rb_job;
  char job; 
  VALUE rb_side;
  char side; 
  VALUE rb_ilo;
  integer ilo; 
  VALUE rb_ihi;
  integer ihi; 
  VALUE rb_scale;
  doublereal *scale; 
  VALUE rb_v;
  doublereal *v; 
  VALUE rb_info;
  integer info; 
  VALUE rb_v_out__;
  doublereal *v_out__;

  integer n;
  integer ldv;
  integer m;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, v = NumRu::Lapack.dgebak( job, side, ilo, ihi, scale, v)\n    or\n  NumRu::Lapack.dgebak  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rb_job = argv[0];
  rb_side = argv[1];
  rb_ilo = argv[2];
  rb_ihi = argv[3];
  rb_scale = argv[4];
  rb_v = argv[5];

  if (!NA_IsNArray(rb_v))
    rb_raise(rb_eArgError, "v (6th argument) must be NArray");
  if (NA_RANK(rb_v) != 2)
    rb_raise(rb_eArgError, "rank of v (6th argument) must be %d", 2);
  m = NA_SHAPE1(rb_v);
  ldv = NA_SHAPE0(rb_v);
  if (NA_TYPE(rb_v) != NA_DFLOAT)
    rb_v = na_change_type(rb_v, NA_DFLOAT);
  v = NA_PTR_TYPE(rb_v, doublereal*);
  if (!NA_IsNArray(rb_scale))
    rb_raise(rb_eArgError, "scale (5th argument) must be NArray");
  if (NA_RANK(rb_scale) != 1)
    rb_raise(rb_eArgError, "rank of scale (5th argument) must be %d", 1);
  n = NA_SHAPE0(rb_scale);
  if (NA_TYPE(rb_scale) != NA_DFLOAT)
    rb_scale = na_change_type(rb_scale, NA_DFLOAT);
  scale = NA_PTR_TYPE(rb_scale, doublereal*);
  ilo = NUM2INT(rb_ilo);
  side = StringValueCStr(rb_side)[0];
  job = StringValueCStr(rb_job)[0];
  ihi = NUM2INT(rb_ihi);
  {
    int shape[2];
    shape[0] = ldv;
    shape[1] = m;
    rb_v_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  v_out__ = NA_PTR_TYPE(rb_v_out__, doublereal*);
  MEMCPY(v_out__, v, doublereal, NA_TOTAL(rb_v));
  rb_v = rb_v_out__;
  v = v_out__;

  dgebak_(&job, &side, &n, &ilo, &ihi, scale, &m, v, &ldv, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_info, rb_v);
}

void
init_lapack_dgebak(VALUE mLapack){
  rb_define_module_function(mLapack, "dgebak", rb_dgebak, -1);
}
