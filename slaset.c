#include "rb_lapack.h"

extern VOID slaset_(char *uplo, integer *m, integer *n, real *alpha, real *beta, real *a, integer *lda);

static VALUE
rb_slaset(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_m;
  integer m; 
  VALUE rb_alpha;
  real alpha; 
  VALUE rb_beta;
  real beta; 
  VALUE rb_a;
  real *a; 
  VALUE rb_a_out__;
  real *a_out__;

  integer lda;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  a = NumRu::Lapack.slaset( uplo, m, alpha, beta, a)\n    or\n  NumRu::Lapack.slaset  # print help\n\n\nFORTRAN MANUAL\n\n");
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
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);
  m = NUM2INT(rb_m);
  beta = (real)NUM2DBL(rb_beta);
  alpha = (real)NUM2DBL(rb_alpha);
  uplo = StringValueCStr(rb_uplo)[0];
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rb_a_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rb_a_out__, real*);
  MEMCPY(a_out__, a, real, NA_TOTAL(rb_a));
  rb_a = rb_a_out__;
  a = a_out__;

  slaset_(&uplo, &m, &n, &alpha, &beta, a, &lda);

  return rb_a;
}

void
init_lapack_slaset(VALUE mLapack){
  rb_define_module_function(mLapack, "slaset", rb_slaset, -1);
}
