#include "rb_lapack.h"

extern VOID csptrf_(char *uplo, integer *n, complex *ap, integer *ipiv, integer *info);

static VALUE
rb_csptrf(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_ap;
  complex *ap; 
  VALUE rb_ipiv;
  integer *ipiv; 
  VALUE rb_info;
  integer info; 
  VALUE rb_ap_out__;
  complex *ap_out__;

  integer ldap;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  ipiv, info, ap = NumRu::Lapack.csptrf( uplo, ap)\n    or\n  NumRu::Lapack.csptrf  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rb_uplo = argv[0];
  rb_ap = argv[1];

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
  {
    int shape[1];
    shape[0] = n;
    rb_ipiv = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  ipiv = NA_PTR_TYPE(rb_ipiv, integer*);
  {
    int shape[1];
    shape[0] = ldap;
    rb_ap_out__ = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  ap_out__ = NA_PTR_TYPE(rb_ap_out__, complex*);
  MEMCPY(ap_out__, ap, complex, NA_TOTAL(rb_ap));
  rb_ap = rb_ap_out__;
  ap = ap_out__;

  csptrf_(&uplo, &n, ap, ipiv, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(3, rb_ipiv, rb_info, rb_ap);
}

void
init_lapack_csptrf(VALUE mLapack){
  rb_define_module_function(mLapack, "csptrf", rb_csptrf, -1);
}
