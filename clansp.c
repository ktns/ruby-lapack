#include "rb_lapack.h"

extern real clansp_(char *norm, char *uplo, integer *n, complex *ap, real *work);

static VALUE
rb_clansp(int argc, VALUE *argv, VALUE self){
  VALUE rb_norm;
  char norm; 
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_n;
  integer n; 
  VALUE rb_ap;
  complex *ap; 
  VALUE rb___out__;
  real __out__; 
  real *work;


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.clansp( norm, uplo, n, ap)\n    or\n  NumRu::Lapack.clansp  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_norm = argv[0];
  rb_uplo = argv[1];
  rb_n = argv[2];
  rb_ap = argv[3];

  n = NUM2INT(rb_n);
  uplo = StringValueCStr(rb_uplo)[0];
  norm = StringValueCStr(rb_norm)[0];
  if (!NA_IsNArray(rb_ap))
    rb_raise(rb_eArgError, "ap (4th argument) must be NArray");
  if (NA_RANK(rb_ap) != 1)
    rb_raise(rb_eArgError, "rank of ap (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_ap) != (n*(n+1)/2))
    rb_raise(rb_eRuntimeError, "shape 0 of ap must be %d", n*(n+1)/2);
  if (NA_TYPE(rb_ap) != NA_SCOMPLEX)
    rb_ap = na_change_type(rb_ap, NA_SCOMPLEX);
  ap = NA_PTR_TYPE(rb_ap, complex*);
  work = ALLOC_N(real, (MAX(1,(lsame_(&norm,"I")||lsame_(&norm,"1")||lsame_(&norm,"O")) ? n : 0)));

  __out__ = clansp_(&norm, &uplo, &n, ap, work);

  free(work);
  rb___out__ = rb_float_new((double)__out__);
  return rb___out__;
}

void
init_lapack_clansp(VALUE mLapack){
  rb_define_module_function(mLapack, "clansp", rb_clansp, -1);
}
