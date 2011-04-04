#include "rb_lapack.h"

extern VOID spstrf_(char *uplo, integer *n, real *a, integer *lda, integer *piv, integer *rank, real *tol, real *work, integer *info);

static VALUE
rb_spstrf(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_a;
  real *a; 
  VALUE rb_tol;
  real tol; 
  VALUE rb_piv;
  integer *piv; 
  VALUE rb_rank;
  integer rank; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  real *a_out__;
  real *work;

  integer lda;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  piv, rank, info, a = NumRu::Lapack.spstrf( uplo, a, tol)\n    or\n  NumRu::Lapack.spstrf  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_uplo = argv[0];
  rb_a = argv[1];
  rb_tol = argv[2];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);
  tol = (real)NUM2DBL(rb_tol);
  uplo = StringValueCStr(rb_uplo)[0];
  {
    int shape[1];
    shape[0] = n;
    rb_piv = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  piv = NA_PTR_TYPE(rb_piv, integer*);
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
  work = ALLOC_N(real, (2*n));

  spstrf_(&uplo, &n, a, &lda, piv, &rank, &tol, work, &info);

  free(work);
  rb_rank = INT2NUM(rank);
  rb_info = INT2NUM(info);
  return rb_ary_new3(4, rb_piv, rb_rank, rb_info, rb_a);
}

void
init_lapack_spstrf(VALUE mLapack){
  rb_define_module_function(mLapack, "spstrf", rb_spstrf, -1);
}
