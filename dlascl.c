#include "rb_lapack.h"

extern VOID dlascl_(char *type, integer *kl, integer *ku, doublereal *cfrom, doublereal *cto, integer *m, integer *n, doublereal *a, integer *lda, integer *info);

static VALUE
rb_dlascl(int argc, VALUE *argv, VALUE self){
  VALUE rb_type;
  char type; 
  VALUE rb_kl;
  integer kl; 
  VALUE rb_ku;
  integer ku; 
  VALUE rb_cfrom;
  doublereal cfrom; 
  VALUE rb_cto;
  doublereal cto; 
  VALUE rb_m;
  integer m; 
  VALUE rb_a;
  doublereal *a; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  doublereal *a_out__;

  integer lda;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, a = NumRu::Lapack.dlascl( type, kl, ku, cfrom, cto, m, a)\n    or\n  NumRu::Lapack.dlascl  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 7)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 7)", argc);
  rb_type = argv[0];
  rb_kl = argv[1];
  rb_ku = argv[2];
  rb_cfrom = argv[3];
  rb_cto = argv[4];
  rb_m = argv[5];
  rb_a = argv[6];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (7th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (7th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_DFLOAT)
    rb_a = na_change_type(rb_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rb_a, doublereal*);
  kl = NUM2INT(rb_kl);
  m = NUM2INT(rb_m);
  cfrom = NUM2DBL(rb_cfrom);
  type = StringValueCStr(rb_type)[0];
  cto = NUM2DBL(rb_cto);
  ku = NUM2INT(rb_ku);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rb_a_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rb_a_out__, doublereal*);
  MEMCPY(a_out__, a, doublereal, NA_TOTAL(rb_a));
  rb_a = rb_a_out__;
  a = a_out__;

  dlascl_(&type, &kl, &ku, &cfrom, &cto, &m, &n, a, &lda, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_info, rb_a);
}

void
init_lapack_dlascl(VALUE mLapack){
  rb_define_module_function(mLapack, "dlascl", rb_dlascl, -1);
}
