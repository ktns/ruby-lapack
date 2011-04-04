#include "rb_lapack.h"

extern VOID chfrk_(char *transr, char *uplo, char *trans, integer *n, integer *k, real *alpha, complex *a, integer *lda, real *beta, complex *c);

static VALUE
rb_chfrk(int argc, VALUE *argv, VALUE self){
  VALUE rb_transr;
  char transr; 
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_trans;
  char trans; 
  VALUE rb_n;
  integer n; 
  VALUE rb_k;
  integer k; 
  VALUE rb_alpha;
  real alpha; 
  VALUE rb_a;
  complex *a; 
  VALUE rb_beta;
  real beta; 
  VALUE rb_c;
  complex *c; 
  VALUE rb_c_out__;
  complex *c_out__;

  integer lda;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  c = NumRu::Lapack.chfrk( transr, uplo, trans, n, k, alpha, a, beta, c)\n    or\n  NumRu::Lapack.chfrk  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 9)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 9)", argc);
  rb_transr = argv[0];
  rb_uplo = argv[1];
  rb_trans = argv[2];
  rb_n = argv[3];
  rb_k = argv[4];
  rb_alpha = argv[5];
  rb_a = argv[6];
  rb_beta = argv[7];
  rb_c = argv[8];

  k = NUM2INT(rb_k);
  uplo = StringValueCStr(rb_uplo)[0];
  trans = StringValueCStr(rb_trans)[0];
  n = NUM2INT(rb_n);
  beta = (real)NUM2DBL(rb_beta);
  alpha = (real)NUM2DBL(rb_alpha);
  transr = StringValueCStr(rb_transr)[0];
  if (!NA_IsNArray(rb_c))
    rb_raise(rb_eArgError, "c (9th argument) must be NArray");
  if (NA_RANK(rb_c) != 1)
    rb_raise(rb_eArgError, "rank of c (9th argument) must be %d", 1);
  if (NA_SHAPE0(rb_c) != (n*(n+1)/2))
    rb_raise(rb_eRuntimeError, "shape 0 of c must be %d", n*(n+1)/2);
  if (NA_TYPE(rb_c) != NA_SCOMPLEX)
    rb_c = na_change_type(rb_c, NA_SCOMPLEX);
  c = NA_PTR_TYPE(rb_c, complex*);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (7th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (7th argument) must be %d", 2);
  if (NA_SHAPE1(rb_a) != (lsame_(&trans,"N") ? k : n))
    rb_raise(rb_eRuntimeError, "shape 1 of a must be %d", lsame_(&trans,"N") ? k : n);
  lda = NA_SHAPE0(rb_a);
  if (lda != (lsame_(&trans,"N") ? MAX(1,n) : MAX(1,k)))
    rb_raise(rb_eRuntimeError, "shape 0 of a must be %d", lsame_(&trans,"N") ? MAX(1,n) : MAX(1,k));
  if (NA_TYPE(rb_a) != NA_SCOMPLEX)
    rb_a = na_change_type(rb_a, NA_SCOMPLEX);
  a = NA_PTR_TYPE(rb_a, complex*);
  lda = lsame_(&trans,"N") ? MAX(1,n) : MAX(1,k);
  {
    int shape[1];
    shape[0] = n*(n+1)/2;
    rb_c_out__ = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  c_out__ = NA_PTR_TYPE(rb_c_out__, complex*);
  MEMCPY(c_out__, c, complex, NA_TOTAL(rb_c));
  rb_c = rb_c_out__;
  c = c_out__;

  chfrk_(&transr, &uplo, &trans, &n, &k, &alpha, a, &lda, &beta, c);

  return rb_c;
}

void
init_lapack_chfrk(VALUE mLapack){
  rb_define_module_function(mLapack, "chfrk", rb_chfrk, -1);
}
