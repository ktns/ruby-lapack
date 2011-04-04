#include "rb_lapack.h"

extern real sla_gercond_(char *trans, integer *n, real *a, integer *lda, real *af, integer *ldaf, integer *ipiv, integer *cmode, real *c, integer *info, real *work, integer *iwork);

static VALUE
rb_sla_gercond(int argc, VALUE *argv, VALUE self){
  VALUE rb_trans;
  char trans; 
  VALUE rb_a;
  real *a; 
  VALUE rb_af;
  real *af; 
  VALUE rb_ipiv;
  integer *ipiv; 
  VALUE rb_cmode;
  integer cmode; 
  VALUE rb_c;
  real *c; 
  VALUE rb_work;
  real *work; 
  VALUE rb_iwork;
  integer *iwork; 
  VALUE rb_info;
  integer info; 
  VALUE rb___out__;
  real __out__; 

  integer lda;
  integer n;
  integer ldaf;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, __out__ = NumRu::Lapack.sla_gercond( trans, a, af, ipiv, cmode, c, work, iwork)\n    or\n  NumRu::Lapack.sla_gercond  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 8)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 8)", argc);
  rb_trans = argv[0];
  rb_a = argv[1];
  rb_af = argv[2];
  rb_ipiv = argv[3];
  rb_cmode = argv[4];
  rb_c = argv[5];
  rb_work = argv[6];
  rb_iwork = argv[7];

  if (!NA_IsNArray(rb_ipiv))
    rb_raise(rb_eArgError, "ipiv (4th argument) must be NArray");
  if (NA_RANK(rb_ipiv) != 1)
    rb_raise(rb_eArgError, "rank of ipiv (4th argument) must be %d", 1);
  n = NA_SHAPE0(rb_ipiv);
  if (NA_TYPE(rb_ipiv) != NA_LINT)
    rb_ipiv = na_change_type(rb_ipiv, NA_LINT);
  ipiv = NA_PTR_TYPE(rb_ipiv, integer*);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  if (NA_SHAPE1(rb_a) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of a must be the same as shape 0 of ipiv");
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);
  if (!NA_IsNArray(rb_c))
    rb_raise(rb_eArgError, "c (6th argument) must be NArray");
  if (NA_RANK(rb_c) != 1)
    rb_raise(rb_eArgError, "rank of c (6th argument) must be %d", 1);
  if (NA_SHAPE0(rb_c) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of c must be the same as shape 0 of ipiv");
  if (NA_TYPE(rb_c) != NA_SFLOAT)
    rb_c = na_change_type(rb_c, NA_SFLOAT);
  c = NA_PTR_TYPE(rb_c, real*);
  if (!NA_IsNArray(rb_af))
    rb_raise(rb_eArgError, "af (3th argument) must be NArray");
  if (NA_RANK(rb_af) != 2)
    rb_raise(rb_eArgError, "rank of af (3th argument) must be %d", 2);
  if (NA_SHAPE1(rb_af) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of af must be the same as shape 0 of ipiv");
  ldaf = NA_SHAPE0(rb_af);
  if (NA_TYPE(rb_af) != NA_SFLOAT)
    rb_af = na_change_type(rb_af, NA_SFLOAT);
  af = NA_PTR_TYPE(rb_af, real*);
  trans = StringValueCStr(rb_trans)[0];
  cmode = NUM2INT(rb_cmode);
  if (!NA_IsNArray(rb_iwork))
    rb_raise(rb_eArgError, "iwork (8th argument) must be NArray");
  if (NA_RANK(rb_iwork) != 1)
    rb_raise(rb_eArgError, "rank of iwork (8th argument) must be %d", 1);
  if (NA_SHAPE0(rb_iwork) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of iwork must be the same as shape 0 of ipiv");
  if (NA_TYPE(rb_iwork) != NA_LINT)
    rb_iwork = na_change_type(rb_iwork, NA_LINT);
  iwork = NA_PTR_TYPE(rb_iwork, integer*);
  if (!NA_IsNArray(rb_work))
    rb_raise(rb_eArgError, "work (7th argument) must be NArray");
  if (NA_RANK(rb_work) != 1)
    rb_raise(rb_eArgError, "rank of work (7th argument) must be %d", 1);
  if (NA_SHAPE0(rb_work) != (3*n))
    rb_raise(rb_eRuntimeError, "shape 0 of work must be %d", 3*n);
  if (NA_TYPE(rb_work) != NA_SFLOAT)
    rb_work = na_change_type(rb_work, NA_SFLOAT);
  work = NA_PTR_TYPE(rb_work, real*);

  __out__ = sla_gercond_(&trans, &n, a, &lda, af, &ldaf, ipiv, &cmode, c, &info, work, iwork);

  rb_info = INT2NUM(info);
  rb___out__ = rb_float_new((double)__out__);
  return rb_ary_new3(2, rb_info, rb___out__);
}

void
init_lapack_sla_gercond(VALUE mLapack){
  rb_define_module_function(mLapack, "sla_gercond", rb_sla_gercond, -1);
}
