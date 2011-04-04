#include "rb_lapack.h"

extern VOID ztfsm_(char *transr, char *side, char *uplo, char *trans, char *diag, integer *m, integer *n, doublecomplex *alpha, doublecomplex *a, doublecomplex *b, integer *ldb);

static VALUE
rb_ztfsm(int argc, VALUE *argv, VALUE self){
  VALUE rb_transr;
  char transr; 
  VALUE rb_side;
  char side; 
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_trans;
  char trans; 
  VALUE rb_diag;
  char diag; 
  VALUE rb_m;
  integer m; 
  VALUE rb_alpha;
  doublecomplex alpha; 
  VALUE rb_a;
  doublecomplex *a; 
  VALUE rb_b;
  doublecomplex *b; 
  VALUE rb_b_out__;
  doublecomplex *b_out__;

  integer ldb;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  b = NumRu::Lapack.ztfsm( transr, side, uplo, trans, diag, m, alpha, a, b)\n    or\n  NumRu::Lapack.ztfsm  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 9)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 9)", argc);
  rb_transr = argv[0];
  rb_side = argv[1];
  rb_uplo = argv[2];
  rb_trans = argv[3];
  rb_diag = argv[4];
  rb_m = argv[5];
  rb_alpha = argv[6];
  rb_a = argv[7];
  rb_b = argv[8];

  transr = StringValueCStr(rb_transr)[0];
  side = StringValueCStr(rb_side)[0];
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (9th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (9th argument) must be %d", 2);
  n = NA_SHAPE1(rb_b);
  ldb = NA_SHAPE0(rb_b);
  if (NA_TYPE(rb_b) != NA_DCOMPLEX)
    rb_b = na_change_type(rb_b, NA_DCOMPLEX);
  b = NA_PTR_TYPE(rb_b, doublecomplex*);
  m = NUM2INT(rb_m);
  uplo = StringValueCStr(rb_uplo)[0];
  alpha.r = NUM2DBL(rb_funcall(rb_alpha, rb_intern("real"), 0));
  alpha.i = NUM2DBL(rb_funcall(rb_alpha, rb_intern("imag"), 0));
  trans = StringValueCStr(rb_trans)[0];
  diag = StringValueCStr(rb_diag)[0];
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (8th argument) must be NArray");
  if (NA_RANK(rb_a) != 1)
    rb_raise(rb_eArgError, "rank of a (8th argument) must be %d", 1);
  if (NA_SHAPE0(rb_a) != (n*(n+1)/2))
    rb_raise(rb_eRuntimeError, "shape 0 of a must be %d", n*(n+1)/2);
  if (NA_TYPE(rb_a) != NA_DCOMPLEX)
    rb_a = na_change_type(rb_a, NA_DCOMPLEX);
  a = NA_PTR_TYPE(rb_a, doublecomplex*);
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = n;
    rb_b_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rb_b_out__, doublecomplex*);
  MEMCPY(b_out__, b, doublecomplex, NA_TOTAL(rb_b));
  rb_b = rb_b_out__;
  b = b_out__;

  ztfsm_(&transr, &side, &uplo, &trans, &diag, &m, &n, &alpha, a, b, &ldb);

  return rb_b;
}

void
init_lapack_ztfsm(VALUE mLapack){
  rb_define_module_function(mLapack, "ztfsm", rb_ztfsm, -1);
}
