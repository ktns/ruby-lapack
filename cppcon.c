#include "rb_lapack.h"

extern VOID cppcon_(char *uplo, integer *n, complex *ap, real *anorm, real *rcond, complex *work, real *rwork, integer *info);

static VALUE
rb_cppcon(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_ap;
  complex *ap; 
  VALUE rb_anorm;
  real anorm; 
  VALUE rb_rcond;
  real rcond; 
  VALUE rb_info;
  integer info; 
  complex *work;
  real *rwork;

  integer ldap;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  rcond, info = NumRu::Lapack.cppcon( uplo, ap, anorm)\n    or\n  NumRu::Lapack.cppcon  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_uplo = argv[0];
  rb_ap = argv[1];
  rb_anorm = argv[2];

  anorm = (real)NUM2DBL(rb_anorm);
  uplo = StringValueCStr(rb_uplo)[0];
  if (!NA_IsNArray(rb_ap))
    rb_raise(rb_eArgError, "ap (2th argument) must be NArray");
  if (NA_RANK(rb_ap) != 1)
    rb_raise(rb_eArgError, "rank of ap (2th argument) must be %d", 1);
  ldap = NA_SHAPE0(rb_ap);
  if (NA_TYPE(rb_ap) != NA_SCOMPLEX)
    rb_ap = na_change_type(rb_ap, NA_SCOMPLEX);
  ap = NA_PTR_TYPE(rb_ap, complex*);
  n = ((int)sqrtf(ldap*8+1.0f)-1)/2;
  work = ALLOC_N(complex, (2*n));
  rwork = ALLOC_N(real, (n));

  cppcon_(&uplo, &n, ap, &anorm, &rcond, work, rwork, &info);

  free(work);
  free(rwork);
  rb_rcond = rb_float_new((double)rcond);
  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_rcond, rb_info);
}

void
init_lapack_cppcon(VALUE mLapack){
  rb_define_module_function(mLapack, "cppcon", rb_cppcon, -1);
}
