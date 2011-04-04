#include "rb_lapack.h"

extern real cla_porpvgrw_(char *uplo, integer *ncols, complex *a, integer *lda, complex *af, integer *ldaf, complex *work);

static VALUE
rb_cla_porpvgrw(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_ncols;
  integer ncols; 
  VALUE rb_a;
  complex *a; 
  VALUE rb_af;
  complex *af; 
  VALUE rb_work;
  complex *work; 
  VALUE rb___out__;
  real __out__; 

  integer lda;
  integer n;
  integer ldaf;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.cla_porpvgrw( uplo, ncols, a, af, work)\n    or\n  NumRu::Lapack.cla_porpvgrw  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_uplo = argv[0];
  rb_ncols = argv[1];
  rb_a = argv[2];
  rb_af = argv[3];
  rb_work = argv[4];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SCOMPLEX)
    rb_a = na_change_type(rb_a, NA_SCOMPLEX);
  a = NA_PTR_TYPE(rb_a, complex*);
  ncols = NUM2INT(rb_ncols);
  uplo = StringValueCStr(rb_uplo)[0];
  if (!NA_IsNArray(rb_af))
    rb_raise(rb_eArgError, "af (4th argument) must be NArray");
  if (NA_RANK(rb_af) != 2)
    rb_raise(rb_eArgError, "rank of af (4th argument) must be %d", 2);
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

  __out__ = cla_porpvgrw_(&uplo, &ncols, a, &lda, af, &ldaf, work);

  rb___out__ = rb_float_new((double)__out__);
  return rb___out__;
}

void
init_lapack_cla_porpvgrw(VALUE mLapack){
  rb_define_module_function(mLapack, "cla_porpvgrw", rb_cla_porpvgrw, -1);
}