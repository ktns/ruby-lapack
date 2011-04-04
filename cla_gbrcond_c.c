#include "rb_lapack.h"

extern real cla_gbrcond_c_(char *trans, integer *n, integer *kl, integer *ku, complex *ab, integer *ldab, complex *afb, integer *ldafb, integer *ipiv, real *c, logical *capply, integer *info, complex *work, real *rwork);

static VALUE
rb_cla_gbrcond_c(int argc, VALUE *argv, VALUE self){
  VALUE rb_trans;
  char trans; 
  VALUE rb_kl;
  integer kl; 
  VALUE rb_ku;
  integer ku; 
  VALUE rb_ab;
  complex *ab; 
  VALUE rb_afb;
  complex *afb; 
  VALUE rb_ipiv;
  integer *ipiv; 
  VALUE rb_c;
  real *c; 
  VALUE rb_capply;
  logical capply; 
  VALUE rb_work;
  complex *work; 
  VALUE rb_rwork;
  real *rwork; 
  VALUE rb_info;
  integer info; 
  VALUE rb___out__;
  real __out__; 

  integer ldab;
  integer n;
  integer ldafb;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, __out__ = NumRu::Lapack.cla_gbrcond_c( trans, kl, ku, ab, afb, ipiv, c, capply, work, rwork)\n    or\n  NumRu::Lapack.cla_gbrcond_c  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 10)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 10)", argc);
  rb_trans = argv[0];
  rb_kl = argv[1];
  rb_ku = argv[2];
  rb_ab = argv[3];
  rb_afb = argv[4];
  rb_ipiv = argv[5];
  rb_c = argv[6];
  rb_capply = argv[7];
  rb_work = argv[8];
  rb_rwork = argv[9];

  if (!NA_IsNArray(rb_ipiv))
    rb_raise(rb_eArgError, "ipiv (6th argument) must be NArray");
  if (NA_RANK(rb_ipiv) != 1)
    rb_raise(rb_eArgError, "rank of ipiv (6th argument) must be %d", 1);
  n = NA_SHAPE0(rb_ipiv);
  if (NA_TYPE(rb_ipiv) != NA_LINT)
    rb_ipiv = na_change_type(rb_ipiv, NA_LINT);
  ipiv = NA_PTR_TYPE(rb_ipiv, integer*);
  if (!NA_IsNArray(rb_ab))
    rb_raise(rb_eArgError, "ab (4th argument) must be NArray");
  if (NA_RANK(rb_ab) != 2)
    rb_raise(rb_eArgError, "rank of ab (4th argument) must be %d", 2);
  if (NA_SHAPE1(rb_ab) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of ab must be the same as shape 0 of ipiv");
  ldab = NA_SHAPE0(rb_ab);
  if (NA_TYPE(rb_ab) != NA_SCOMPLEX)
    rb_ab = na_change_type(rb_ab, NA_SCOMPLEX);
  ab = NA_PTR_TYPE(rb_ab, complex*);
  kl = NUM2INT(rb_kl);
  if (!NA_IsNArray(rb_rwork))
    rb_raise(rb_eArgError, "rwork (10th argument) must be NArray");
  if (NA_RANK(rb_rwork) != 1)
    rb_raise(rb_eArgError, "rank of rwork (10th argument) must be %d", 1);
  if (NA_SHAPE0(rb_rwork) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of rwork must be the same as shape 0 of ipiv");
  if (NA_TYPE(rb_rwork) != NA_SFLOAT)
    rb_rwork = na_change_type(rb_rwork, NA_SFLOAT);
  rwork = NA_PTR_TYPE(rb_rwork, real*);
  if (!NA_IsNArray(rb_c))
    rb_raise(rb_eArgError, "c (7th argument) must be NArray");
  if (NA_RANK(rb_c) != 1)
    rb_raise(rb_eArgError, "rank of c (7th argument) must be %d", 1);
  if (NA_SHAPE0(rb_c) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of c must be the same as shape 0 of ipiv");
  if (NA_TYPE(rb_c) != NA_SFLOAT)
    rb_c = na_change_type(rb_c, NA_SFLOAT);
  c = NA_PTR_TYPE(rb_c, real*);
  ku = NUM2INT(rb_ku);
  if (!NA_IsNArray(rb_afb))
    rb_raise(rb_eArgError, "afb (5th argument) must be NArray");
  if (NA_RANK(rb_afb) != 2)
    rb_raise(rb_eArgError, "rank of afb (5th argument) must be %d", 2);
  if (NA_SHAPE1(rb_afb) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of afb must be the same as shape 0 of ipiv");
  ldafb = NA_SHAPE0(rb_afb);
  if (NA_TYPE(rb_afb) != NA_SCOMPLEX)
    rb_afb = na_change_type(rb_afb, NA_SCOMPLEX);
  afb = NA_PTR_TYPE(rb_afb, complex*);
  trans = StringValueCStr(rb_trans)[0];
  capply = (rb_capply == Qtrue);
  if (!NA_IsNArray(rb_work))
    rb_raise(rb_eArgError, "work (9th argument) must be NArray");
  if (NA_RANK(rb_work) != 1)
    rb_raise(rb_eArgError, "rank of work (9th argument) must be %d", 1);
  if (NA_SHAPE0(rb_work) != (2*n))
    rb_raise(rb_eRuntimeError, "shape 0 of work must be %d", 2*n);
  if (NA_TYPE(rb_work) != NA_SCOMPLEX)
    rb_work = na_change_type(rb_work, NA_SCOMPLEX);
  work = NA_PTR_TYPE(rb_work, complex*);

  __out__ = cla_gbrcond_c_(&trans, &n, &kl, &ku, ab, &ldab, afb, &ldafb, ipiv, c, &capply, &info, work, rwork);

  rb_info = INT2NUM(info);
  rb___out__ = rb_float_new((double)__out__);
  return rb_ary_new3(2, rb_info, rb___out__);
}

void
init_lapack_cla_gbrcond_c(VALUE mLapack){
  rb_define_module_function(mLapack, "cla_gbrcond_c", rb_cla_gbrcond_c, -1);
}
