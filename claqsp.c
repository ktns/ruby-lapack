#include "rb_lapack.h"

extern VOID claqsp_(char *uplo, integer *n, complex *ap, real *s, real *scond, real *amax, char *equed);

static VALUE
rb_claqsp(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_ap;
  complex *ap; 
  VALUE rb_s;
  real *s; 
  VALUE rb_scond;
  real scond; 
  VALUE rb_amax;
  real amax; 
  VALUE rb_equed;
  char equed; 
  VALUE rb_ap_out__;
  complex *ap_out__;

  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  equed, ap = NumRu::Lapack.claqsp( uplo, ap, s, scond, amax)\n    or\n  NumRu::Lapack.claqsp  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_uplo = argv[0];
  rb_ap = argv[1];
  rb_s = argv[2];
  rb_scond = argv[3];
  rb_amax = argv[4];

  scond = (real)NUM2DBL(rb_scond);
  uplo = StringValueCStr(rb_uplo)[0];
  if (!NA_IsNArray(rb_s))
    rb_raise(rb_eArgError, "s (3th argument) must be NArray");
  if (NA_RANK(rb_s) != 1)
    rb_raise(rb_eArgError, "rank of s (3th argument) must be %d", 1);
  n = NA_SHAPE0(rb_s);
  if (NA_TYPE(rb_s) != NA_SFLOAT)
    rb_s = na_change_type(rb_s, NA_SFLOAT);
  s = NA_PTR_TYPE(rb_s, real*);
  amax = (real)NUM2DBL(rb_amax);
  if (!NA_IsNArray(rb_ap))
    rb_raise(rb_eArgError, "ap (2th argument) must be NArray");
  if (NA_RANK(rb_ap) != 1)
    rb_raise(rb_eArgError, "rank of ap (2th argument) must be %d", 1);
  if (NA_SHAPE0(rb_ap) != (n*(n+1)/2))
    rb_raise(rb_eRuntimeError, "shape 0 of ap must be %d", n*(n+1)/2);
  if (NA_TYPE(rb_ap) != NA_SCOMPLEX)
    rb_ap = na_change_type(rb_ap, NA_SCOMPLEX);
  ap = NA_PTR_TYPE(rb_ap, complex*);
  {
    int shape[1];
    shape[0] = n*(n+1)/2;
    rb_ap_out__ = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  ap_out__ = NA_PTR_TYPE(rb_ap_out__, complex*);
  MEMCPY(ap_out__, ap, complex, NA_TOTAL(rb_ap));
  rb_ap = rb_ap_out__;
  ap = ap_out__;

  claqsp_(&uplo, &n, ap, s, &scond, &amax, &equed);

  rb_equed = rb_str_new(&equed,1);
  return rb_ary_new3(2, rb_equed, rb_ap);
}

void
init_lapack_claqsp(VALUE mLapack){
  rb_define_module_function(mLapack, "claqsp", rb_claqsp, -1);
}
