#include "rb_lapack.h"

extern VOID zlaset_(char *uplo, integer *m, integer *n, doublecomplex *alpha, doublecomplex *beta, doublecomplex *a, integer *lda);

static VALUE
rb_zlaset(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_m;
  integer m; 
  VALUE rb_alpha;
  doublecomplex alpha; 
  VALUE rb_beta;
  doublecomplex beta; 
  VALUE rb_a;
  doublecomplex *a; 
  VALUE rb_a_out__;
  doublecomplex *a_out__;

  integer lda;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  a = NumRu::Lapack.zlaset( uplo, m, alpha, beta, a)\n    or\n  NumRu::Lapack.zlaset  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_uplo = argv[0];
  rb_m = argv[1];
  rb_alpha = argv[2];
  rb_beta = argv[3];
  rb_a = argv[4];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (5th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (5th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_DCOMPLEX)
    rb_a = na_change_type(rb_a, NA_DCOMPLEX);
  a = NA_PTR_TYPE(rb_a, doublecomplex*);
  m = NUM2INT(rb_m);
  beta.r = NUM2DBL(rb_funcall(rb_beta, rb_intern("real"), 0));
  beta.i = NUM2DBL(rb_funcall(rb_beta, rb_intern("imag"), 0));
  alpha.r = NUM2DBL(rb_funcall(rb_alpha, rb_intern("real"), 0));
  alpha.i = NUM2DBL(rb_funcall(rb_alpha, rb_intern("imag"), 0));
  uplo = StringValueCStr(rb_uplo)[0];
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rb_a_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rb_a_out__, doublecomplex*);
  MEMCPY(a_out__, a, doublecomplex, NA_TOTAL(rb_a));
  rb_a = rb_a_out__;
  a = a_out__;

  zlaset_(&uplo, &m, &n, &alpha, &beta, a, &lda);

  return rb_a;
}

void
init_lapack_zlaset(VALUE mLapack){
  rb_define_module_function(mLapack, "zlaset", rb_zlaset, -1);
}
