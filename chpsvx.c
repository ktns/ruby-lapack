#include "rb_lapack.h"

extern VOID chpsvx_(char *fact, char *uplo, integer *n, integer *nrhs, complex *ap, complex *afp, integer *ipiv, complex *b, integer *ldb, complex *x, integer *ldx, real *rcond, real *ferr, real *berr, complex *work, real *rwork, integer *info);

static VALUE
rb_chpsvx(int argc, VALUE *argv, VALUE self){
  VALUE rb_fact;
  char fact; 
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_ap;
  complex *ap; 
  VALUE rb_afp;
  complex *afp; 
  VALUE rb_ipiv;
  integer *ipiv; 
  VALUE rb_b;
  complex *b; 
  VALUE rb_x;
  complex *x; 
  VALUE rb_rcond;
  real rcond; 
  VALUE rb_ferr;
  real *ferr; 
  VALUE rb_berr;
  real *berr; 
  VALUE rb_info;
  integer info; 
  VALUE rb_afp_out__;
  complex *afp_out__;
  VALUE rb_ipiv_out__;
  integer *ipiv_out__;
  complex *work;
  real *rwork;

  integer n;
  integer ldb;
  integer nrhs;
  integer ldx;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  x, rcond, ferr, berr, info, afp, ipiv = NumRu::Lapack.chpsvx( fact, uplo, ap, afp, ipiv, b)\n    or\n  NumRu::Lapack.chpsvx  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rb_fact = argv[0];
  rb_uplo = argv[1];
  rb_ap = argv[2];
  rb_afp = argv[3];
  rb_ipiv = argv[4];
  rb_b = argv[5];

  uplo = StringValueCStr(rb_uplo)[0];
  if (!NA_IsNArray(rb_ipiv))
    rb_raise(rb_eArgError, "ipiv (5th argument) must be NArray");
  if (NA_RANK(rb_ipiv) != 1)
    rb_raise(rb_eArgError, "rank of ipiv (5th argument) must be %d", 1);
  n = NA_SHAPE0(rb_ipiv);
  if (NA_TYPE(rb_ipiv) != NA_LINT)
    rb_ipiv = na_change_type(rb_ipiv, NA_LINT);
  ipiv = NA_PTR_TYPE(rb_ipiv, integer*);
  fact = StringValueCStr(rb_fact)[0];
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (6th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (6th argument) must be %d", 2);
  nrhs = NA_SHAPE1(rb_b);
  ldb = NA_SHAPE0(rb_b);
  if (NA_TYPE(rb_b) != NA_SCOMPLEX)
    rb_b = na_change_type(rb_b, NA_SCOMPLEX);
  b = NA_PTR_TYPE(rb_b, complex*);
  if (!NA_IsNArray(rb_afp))
    rb_raise(rb_eArgError, "afp (4th argument) must be NArray");
  if (NA_RANK(rb_afp) != 1)
    rb_raise(rb_eArgError, "rank of afp (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_afp) != (n*(n+1)/2))
    rb_raise(rb_eRuntimeError, "shape 0 of afp must be %d", n*(n+1)/2);
  if (NA_TYPE(rb_afp) != NA_SCOMPLEX)
    rb_afp = na_change_type(rb_afp, NA_SCOMPLEX);
  afp = NA_PTR_TYPE(rb_afp, complex*);
  if (!NA_IsNArray(rb_ap))
    rb_raise(rb_eArgError, "ap (3th argument) must be NArray");
  if (NA_RANK(rb_ap) != 1)
    rb_raise(rb_eArgError, "rank of ap (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_ap) != (n*(n+1)/2))
    rb_raise(rb_eRuntimeError, "shape 0 of ap must be %d", n*(n+1)/2);
  if (NA_TYPE(rb_ap) != NA_SCOMPLEX)
    rb_ap = na_change_type(rb_ap, NA_SCOMPLEX);
  ap = NA_PTR_TYPE(rb_ap, complex*);
  ldx = MAX(1,n);
  {
    int shape[2];
    shape[0] = ldx;
    shape[1] = nrhs;
    rb_x = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  x = NA_PTR_TYPE(rb_x, complex*);
  {
    int shape[1];
    shape[0] = nrhs;
    rb_ferr = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  ferr = NA_PTR_TYPE(rb_ferr, real*);
  {
    int shape[1];
    shape[0] = nrhs;
    rb_berr = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  berr = NA_PTR_TYPE(rb_berr, real*);
  {
    int shape[1];
    shape[0] = n*(n+1)/2;
    rb_afp_out__ = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  afp_out__ = NA_PTR_TYPE(rb_afp_out__, complex*);
  MEMCPY(afp_out__, afp, complex, NA_TOTAL(rb_afp));
  rb_afp = rb_afp_out__;
  afp = afp_out__;
  {
    int shape[1];
    shape[0] = n;
    rb_ipiv_out__ = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  ipiv_out__ = NA_PTR_TYPE(rb_ipiv_out__, integer*);
  MEMCPY(ipiv_out__, ipiv, integer, NA_TOTAL(rb_ipiv));
  rb_ipiv = rb_ipiv_out__;
  ipiv = ipiv_out__;
  work = ALLOC_N(complex, (2*n));
  rwork = ALLOC_N(real, (n));

  chpsvx_(&fact, &uplo, &n, &nrhs, ap, afp, ipiv, b, &ldb, x, &ldx, &rcond, ferr, berr, work, rwork, &info);

  free(work);
  free(rwork);
  rb_rcond = rb_float_new((double)rcond);
  rb_info = INT2NUM(info);
  return rb_ary_new3(7, rb_x, rb_rcond, rb_ferr, rb_berr, rb_info, rb_afp, rb_ipiv);
}

void
init_lapack_chpsvx(VALUE mLapack){
  rb_define_module_function(mLapack, "chpsvx", rb_chpsvx, -1);
}
