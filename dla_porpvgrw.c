#include "rb_lapack.h"

extern doublereal dla_porpvgrw_(char *uplo, integer *ncols, doublereal *a, integer *lda, doublereal *af, integer *ldaf, doublereal *work);

static VALUE
rb_dla_porpvgrw(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_ncols;
  integer ncols; 
  VALUE rb_a;
  doublereal *a; 
  VALUE rb_af;
  doublereal *af; 
  VALUE rb_work;
  doublereal *work; 
  VALUE rb___out__;
  doublereal __out__; 

  integer lda;
  integer n;
  integer ldaf;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.dla_porpvgrw( uplo, ncols, a, af, work)\n    or\n  NumRu::Lapack.dla_porpvgrw  # print help\n\n\nFORTRAN MANUAL\n\n");
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
  if (NA_TYPE(rb_a) != NA_DFLOAT)
    rb_a = na_change_type(rb_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rb_a, doublereal*);
  ncols = NUM2INT(rb_ncols);
  uplo = StringValueCStr(rb_uplo)[0];
  if (!NA_IsNArray(rb_af))
    rb_raise(rb_eArgError, "af (4th argument) must be NArray");
  if (NA_RANK(rb_af) != 2)
    rb_raise(rb_eArgError, "rank of af (4th argument) must be %d", 2);
  if (NA_SHAPE1(rb_af) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of af must be the same as shape 1 of a");
  ldaf = NA_SHAPE0(rb_af);
  if (NA_TYPE(rb_af) != NA_DFLOAT)
    rb_af = na_change_type(rb_af, NA_DFLOAT);
  af = NA_PTR_TYPE(rb_af, doublereal*);
  if (!NA_IsNArray(rb_work))
    rb_raise(rb_eArgError, "work (5th argument) must be NArray");
  if (NA_RANK(rb_work) != 1)
    rb_raise(rb_eArgError, "rank of work (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_work) != (2*n))
    rb_raise(rb_eRuntimeError, "shape 0 of work must be %d", 2*n);
  if (NA_TYPE(rb_work) != NA_DFLOAT)
    rb_work = na_change_type(rb_work, NA_DFLOAT);
  work = NA_PTR_TYPE(rb_work, doublereal*);

  __out__ = dla_porpvgrw_(&uplo, &ncols, a, &lda, af, &ldaf, work);

  rb___out__ = rb_float_new((double)__out__);
  return rb___out__;
}

void
init_lapack_dla_porpvgrw(VALUE mLapack){
  rb_define_module_function(mLapack, "dla_porpvgrw", rb_dla_porpvgrw, -1);
}
