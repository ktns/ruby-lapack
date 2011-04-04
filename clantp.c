#include "rb_lapack.h"

extern real clantp_(char *norm, char *uplo, char *diag, integer *n, complex *ap, real *work);

static VALUE
rb_clantp(int argc, VALUE *argv, VALUE self){
  VALUE rb_norm;
  char norm; 
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_diag;
  char diag; 
  VALUE rb_n;
  integer n; 
  VALUE rb_ap;
  complex *ap; 
  VALUE rb___out__;
  real __out__; 
  real *work;

  integer lwork;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.clantp( norm, uplo, diag, n, ap)\n    or\n  NumRu::Lapack.clantp  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_norm = argv[0];
  rb_uplo = argv[1];
  rb_diag = argv[2];
  rb_n = argv[3];
  rb_ap = argv[4];

  diag = StringValueCStr(rb_diag)[0];
  n = NUM2INT(rb_n);
  norm = StringValueCStr(rb_norm)[0];
  uplo = StringValueCStr(rb_uplo)[0];
  if (!NA_IsNArray(rb_ap))
    rb_raise(rb_eArgError, "ap (5th argument) must be NArray");
  if (NA_RANK(rb_ap) != 1)
    rb_raise(rb_eArgError, "rank of ap (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_ap) != (n*(n+1)/2))
    rb_raise(rb_eRuntimeError, "shape 0 of ap must be %d", n*(n+1)/2);
  if (NA_TYPE(rb_ap) != NA_SCOMPLEX)
    rb_ap = na_change_type(rb_ap, NA_SCOMPLEX);
  ap = NA_PTR_TYPE(rb_ap, complex*);
  lwork = lsame_(&norm,"I") ? n : 0;
  work = ALLOC_N(real, (MAX(1,lwork)));

  __out__ = clantp_(&norm, &uplo, &diag, &n, ap, work);

  free(work);
  rb___out__ = rb_float_new((double)__out__);
  return rb___out__;
}

void
init_lapack_clantp(VALUE mLapack){
  rb_define_module_function(mLapack, "clantp", rb_clantp, -1);
}
