#include "rb_lapack.h"

extern VOID zcgesv_(integer *n, integer *nrhs, doublecomplex *a, integer *lda, integer *ipiv, doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx, doublecomplex *work, complex *swork, doublereal *rwork, integer *iter, integer *info);

static VALUE
rb_zcgesv(int argc, VALUE *argv, VALUE self){
  VALUE rb_a;
  doublecomplex *a; 
  VALUE rb_b;
  doublecomplex *b; 
  VALUE rb_ipiv;
  integer *ipiv; 
  VALUE rb_x;
  doublecomplex *x; 
  VALUE rb_iter;
  integer iter; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  doublecomplex *a_out__;
  doublecomplex *work;
  complex *swork;
  doublereal *rwork;

  integer lda;
  integer n;
  integer ldb;
  integer nrhs;
  integer ldx;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  ipiv, x, iter, info, a = NumRu::Lapack.zcgesv( a, b)\n    or\n  NumRu::Lapack.zcgesv  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rb_a = argv[0];
  rb_b = argv[1];

  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (2th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (2th argument) must be %d", 2);
  nrhs = NA_SHAPE1(rb_b);
  ldb = NA_SHAPE0(rb_b);
  if (NA_TYPE(rb_b) != NA_DCOMPLEX)
    rb_b = na_change_type(rb_b, NA_DCOMPLEX);
  b = NA_PTR_TYPE(rb_b, doublecomplex*);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (1th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (1th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_DCOMPLEX)
    rb_a = na_change_type(rb_a, NA_DCOMPLEX);
  a = NA_PTR_TYPE(rb_a, doublecomplex*);
  ldx = MAX(1,n);
  {
    int shape[1];
    shape[0] = n;
    rb_ipiv = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  ipiv = NA_PTR_TYPE(rb_ipiv, integer*);
  {
    int shape[2];
    shape[0] = ldx;
    shape[1] = nrhs;
    rb_x = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  x = NA_PTR_TYPE(rb_x, doublecomplex*);
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
  work = ALLOC_N(doublecomplex, (n*nrhs));
  swork = ALLOC_N(complex, (n*(n+nrhs)));
  rwork = ALLOC_N(doublereal, (n));

  zcgesv_(&n, &nrhs, a, &lda, ipiv, b, &ldb, x, &ldx, work, swork, rwork, &iter, &info);

  free(work);
  free(swork);
  free(rwork);
  rb_iter = INT2NUM(iter);
  rb_info = INT2NUM(info);
  return rb_ary_new3(5, rb_ipiv, rb_x, rb_iter, rb_info, rb_a);
}

void
init_lapack_zcgesv(VALUE mLapack){
  rb_define_module_function(mLapack, "zcgesv", rb_zcgesv, -1);
}
