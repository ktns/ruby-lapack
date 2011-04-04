#include "rb_lapack.h"

extern VOID clacrm_(integer *m, integer *n, complex *a, integer *lda, real *b, integer *ldb, complex *c, integer *ldc, real *rwork);

static VALUE
rb_clacrm(int argc, VALUE *argv, VALUE self){
  VALUE rb_m;
  integer m; 
  VALUE rb_a;
  complex *a; 
  VALUE rb_b;
  real *b; 
  VALUE rb_c;
  complex *c; 
  real *rwork;

  integer lda;
  integer n;
  integer ldb;
  integer ldc;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  c = NumRu::Lapack.clacrm( m, a, b)\n    or\n  NumRu::Lapack.clacrm  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_m = argv[0];
  rb_a = argv[1];
  rb_b = argv[2];

  m = NUM2INT(rb_m);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SCOMPLEX)
    rb_a = na_change_type(rb_a, NA_SCOMPLEX);
  a = NA_PTR_TYPE(rb_a, complex*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (3th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (3th argument) must be %d", 2);
  if (NA_SHAPE1(rb_b) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of b must be the same as shape 1 of a");
  ldb = NA_SHAPE0(rb_b);
  if (NA_TYPE(rb_b) != NA_SFLOAT)
    rb_b = na_change_type(rb_b, NA_SFLOAT);
  b = NA_PTR_TYPE(rb_b, real*);
  ldc = MAX(1,n);
  {
    int shape[2];
    shape[0] = ldc;
    shape[1] = n;
    rb_c = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  c = NA_PTR_TYPE(rb_c, complex*);
  rwork = ALLOC_N(real, (2*m*n));

  clacrm_(&m, &n, a, &lda, b, &ldb, c, &ldc, rwork);

  free(rwork);
  return rb_c;
}

void
init_lapack_clacrm(VALUE mLapack){
  rb_define_module_function(mLapack, "clacrm", rb_clacrm, -1);
}
