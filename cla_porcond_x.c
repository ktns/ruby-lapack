#include "rb_lapack.h"

extern real cla_porcond_x_(char *uplo, integer *n, complex *a, integer *lda, complex *af, integer *ldaf, complex *x, integer *info, complex *work, real *rwork);

static VALUE
rb_cla_porcond_x(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_a;
  complex *a; 
  VALUE rb_af;
  complex *af; 
  VALUE rb_x;
  complex *x; 
  VALUE rb_work;
  complex *work; 
  VALUE rb_rwork;
  real *rwork; 
  VALUE rb_info;
  integer info; 
  VALUE rb___out__;
  real __out__; 

  integer lda;
  integer n;
  integer ldaf;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, __out__ = NumRu::Lapack.cla_porcond_x( uplo, a, af, x, work, rwork)\n    or\n  NumRu::Lapack.cla_porcond_x  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rb_uplo = argv[0];
  rb_a = argv[1];
  rb_af = argv[2];
  rb_x = argv[3];
  rb_work = argv[4];
  rb_rwork = argv[5];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SCOMPLEX)
    rb_a = na_change_type(rb_a, NA_SCOMPLEX);
  a = NA_PTR_TYPE(rb_a, complex*);
  if (!NA_IsNArray(rb_x))
    rb_raise(rb_eArgError, "x (4th argument) must be NArray");
  if (NA_RANK(rb_x) != 1)
    rb_raise(rb_eArgError, "rank of x (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_x) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of x must be the same as shape 1 of a");
  if (NA_TYPE(rb_x) != NA_SCOMPLEX)
    rb_x = na_change_type(rb_x, NA_SCOMPLEX);
  x = NA_PTR_TYPE(rb_x, complex*);
  if (!NA_IsNArray(rb_rwork))
    rb_raise(rb_eArgError, "rwork (6th argument) must be NArray");
  if (NA_RANK(rb_rwork) != 1)
    rb_raise(rb_eArgError, "rank of rwork (6th argument) must be %d", 1);
  if (NA_SHAPE0(rb_rwork) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of rwork must be the same as shape 1 of a");
  if (NA_TYPE(rb_rwork) != NA_SFLOAT)
    rb_rwork = na_change_type(rb_rwork, NA_SFLOAT);
  rwork = NA_PTR_TYPE(rb_rwork, real*);
  uplo = StringValueCStr(rb_uplo)[0];
  if (!NA_IsNArray(rb_af))
    rb_raise(rb_eArgError, "af (3th argument) must be NArray");
  if (NA_RANK(rb_af) != 2)
    rb_raise(rb_eArgError, "rank of af (3th argument) must be %d", 2);
  if (NA_SHAPE1(rb_af) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of af must be the same as shape 1 of a");
  ldaf = NA_SHAPE0(rb_af);
  if (NA_TYPE(rb_af) != NA_SCOMPLEX)
    rb_af = na_change_type(rb_af, NA_SCOMPLEX);
  af = NA_PTR_TYPE(rb_af, complex*);
  if (!NA_IsNArray(rb_work))
    rb_raise(rb_eArgError, "work (5th argument) must be NArray");
  if (NA_RANK(rb_work) != 1)
    rb_raise(rb_eArgError, "rank of work (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_work) != (2*n))
    rb_raise(rb_eRuntimeError, "shape 0 of work must be %d", 2*n);
  if (NA_TYPE(rb_work) != NA_SCOMPLEX)
    rb_work = na_change_type(rb_work, NA_SCOMPLEX);
  work = NA_PTR_TYPE(rb_work, complex*);

  __out__ = cla_porcond_x_(&uplo, &n, a, &lda, af, &ldaf, x, &info, work, rwork);

  rb_info = INT2NUM(info);
  rb___out__ = rb_float_new((double)__out__);
  return rb_ary_new3(2, rb_info, rb___out__);
}

void
init_lapack_cla_porcond_x(VALUE mLapack){
  rb_define_module_function(mLapack, "cla_porcond_x", rb_cla_porcond_x, -1);
}
