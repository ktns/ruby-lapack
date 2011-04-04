#include "rb_lapack.h"

extern VOID dlaqsp_(char *uplo, integer *n, doublereal *ap, doublereal *s, doublereal *scond, doublereal *amax, char *equed);

static VALUE
rb_dlaqsp(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_ap;
  doublereal *ap; 
  VALUE rb_s;
  doublereal *s; 
  VALUE rb_scond;
  doublereal scond; 
  VALUE rb_amax;
  doublereal amax; 
  VALUE rb_equed;
  char equed; 
  VALUE rb_ap_out__;
  doublereal *ap_out__;

  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  equed, ap = NumRu::Lapack.dlaqsp( uplo, ap, s, scond, amax)\n    or\n  NumRu::Lapack.dlaqsp  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_uplo = argv[0];
  rb_ap = argv[1];
  rb_s = argv[2];
  rb_scond = argv[3];
  rb_amax = argv[4];

  scond = NUM2DBL(rb_scond);
  uplo = StringValueCStr(rb_uplo)[0];
  if (!NA_IsNArray(rb_s))
    rb_raise(rb_eArgError, "s (3th argument) must be NArray");
  if (NA_RANK(rb_s) != 1)
    rb_raise(rb_eArgError, "rank of s (3th argument) must be %d", 1);
  n = NA_SHAPE0(rb_s);
  if (NA_TYPE(rb_s) != NA_DFLOAT)
    rb_s = na_change_type(rb_s, NA_DFLOAT);
  s = NA_PTR_TYPE(rb_s, doublereal*);
  amax = NUM2DBL(rb_amax);
  if (!NA_IsNArray(rb_ap))
    rb_raise(rb_eArgError, "ap (2th argument) must be NArray");
  if (NA_RANK(rb_ap) != 1)
    rb_raise(rb_eArgError, "rank of ap (2th argument) must be %d", 1);
  if (NA_SHAPE0(rb_ap) != (n*(n+1)/2))
    rb_raise(rb_eRuntimeError, "shape 0 of ap must be %d", n*(n+1)/2);
  if (NA_TYPE(rb_ap) != NA_DFLOAT)
    rb_ap = na_change_type(rb_ap, NA_DFLOAT);
  ap = NA_PTR_TYPE(rb_ap, doublereal*);
  {
    int shape[1];
    shape[0] = n*(n+1)/2;
    rb_ap_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  ap_out__ = NA_PTR_TYPE(rb_ap_out__, doublereal*);
  MEMCPY(ap_out__, ap, doublereal, NA_TOTAL(rb_ap));
  rb_ap = rb_ap_out__;
  ap = ap_out__;

  dlaqsp_(&uplo, &n, ap, s, &scond, &amax, &equed);

  rb_equed = rb_str_new(&equed,1);
  return rb_ary_new3(2, rb_equed, rb_ap);
}

void
init_lapack_dlaqsp(VALUE mLapack){
  rb_define_module_function(mLapack, "dlaqsp", rb_dlaqsp, -1);
}
