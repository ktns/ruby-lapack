#include "rb_lapack.h"

extern VOID stfsm_(char *transr, char *side, char *uplo, char *trans, char *diag, integer *m, integer *n, real *alpha, real *a, real *b, integer *ldb);

static VALUE
rb_stfsm(int argc, VALUE *argv, VALUE self){
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
  real alpha; 
  VALUE rb_a;
  real *a; 
  VALUE rb_b;
  real *b; 
  VALUE rb_b_out__;
  real *b_out__;

  integer nt;
  integer ldb;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  b = NumRu::Lapack.stfsm( transr, side, uplo, trans, diag, m, alpha, a, b)\n    or\n  NumRu::Lapack.stfsm  # print help\n\n\nFORTRAN MANUAL\n\n");
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

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (8th argument) must be NArray");
  if (NA_RANK(rb_a) != 1)
    rb_raise(rb_eArgError, "rank of a (8th argument) must be %d", 1);
  nt = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);
  side = StringValueCStr(rb_side)[0];
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (9th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (9th argument) must be %d", 2);
  n = NA_SHAPE1(rb_b);
  ldb = NA_SHAPE0(rb_b);
  if (ldb != (MAX(1,m)))
    rb_raise(rb_eRuntimeError, "shape 0 of b must be %d", MAX(1,m));
  if (NA_TYPE(rb_b) != NA_SFLOAT)
    rb_b = na_change_type(rb_b, NA_SFLOAT);
  b = NA_PTR_TYPE(rb_b, real*);
  m = NUM2INT(rb_m);
  diag = StringValueCStr(rb_diag)[0];
  uplo = StringValueCStr(rb_uplo)[0];
  alpha = (real)NUM2DBL(rb_alpha);
  trans = StringValueCStr(rb_trans)[0];
  transr = StringValueCStr(rb_transr)[0];
  ldb = MAX(1,m);
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = n;
    rb_b_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rb_b_out__, real*);
  MEMCPY(b_out__, b, real, NA_TOTAL(rb_b));
  rb_b = rb_b_out__;
  b = b_out__;

  stfsm_(&transr, &side, &uplo, &trans, &diag, &m, &n, &alpha, a, b, &ldb);

  return rb_b;
}

void
init_lapack_stfsm(VALUE mLapack){
  rb_define_module_function(mLapack, "stfsm", rb_stfsm, -1);
}
