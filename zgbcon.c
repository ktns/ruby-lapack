#include "rb_lapack.h"

extern VOID zgbcon_(char *norm, integer *n, integer *kl, integer *ku, doublecomplex *ab, integer *ldab, integer *ipiv, doublereal *anorm, doublereal *rcond, doublecomplex *work, doublereal *rwork, integer *info);

static VALUE
rb_zgbcon(int argc, VALUE *argv, VALUE self){
  VALUE rb_norm;
  char norm; 
  VALUE rb_kl;
  integer kl; 
  VALUE rb_ku;
  integer ku; 
  VALUE rb_ab;
  doublecomplex *ab; 
  VALUE rb_ipiv;
  integer *ipiv; 
  VALUE rb_anorm;
  doublereal anorm; 
  VALUE rb_rcond;
  doublereal rcond; 
  VALUE rb_info;
  integer info; 
  doublecomplex *work;
  doublereal *rwork;

  integer ldab;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  rcond, info = NumRu::Lapack.zgbcon( norm, kl, ku, ab, ipiv, anorm)\n    or\n  NumRu::Lapack.zgbcon  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rb_norm = argv[0];
  rb_kl = argv[1];
  rb_ku = argv[2];
  rb_ab = argv[3];
  rb_ipiv = argv[4];
  rb_anorm = argv[5];

  anorm = NUM2DBL(rb_anorm);
  if (!NA_IsNArray(rb_ipiv))
    rb_raise(rb_eArgError, "ipiv (5th argument) must be NArray");
  if (NA_RANK(rb_ipiv) != 1)
    rb_raise(rb_eArgError, "rank of ipiv (5th argument) must be %d", 1);
  n = NA_SHAPE0(rb_ipiv);
  if (NA_TYPE(rb_ipiv) != NA_LINT)
    rb_ipiv = na_change_type(rb_ipiv, NA_LINT);
  ipiv = NA_PTR_TYPE(rb_ipiv, integer*);
  if (!NA_IsNArray(rb_ab))
    rb_raise(rb_eArgError, "ab (4th argument) must be NArray");
  if (NA_RANK(rb_ab) != 2)
    rb_raise(rb_eArgError, "rank of ab (4th argument) must be %d", 2);
  if (NA_SHAPE1(rb_ab) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of ab must be the same as shape 0 of ipiv");
  ldab = NA_SHAPE0(rb_ab);
  if (NA_TYPE(rb_ab) != NA_DCOMPLEX)
    rb_ab = na_change_type(rb_ab, NA_DCOMPLEX);
  ab = NA_PTR_TYPE(rb_ab, doublecomplex*);
  kl = NUM2INT(rb_kl);
  norm = StringValueCStr(rb_norm)[0];
  ku = NUM2INT(rb_ku);
  work = ALLOC_N(doublecomplex, (2*n));
  rwork = ALLOC_N(doublereal, (n));

  zgbcon_(&norm, &n, &kl, &ku, ab, &ldab, ipiv, &anorm, &rcond, work, rwork, &info);

  free(work);
  free(rwork);
  rb_rcond = rb_float_new((double)rcond);
  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_rcond, rb_info);
}

void
init_lapack_zgbcon(VALUE mLapack){
  rb_define_module_function(mLapack, "zgbcon", rb_zgbcon, -1);
}
