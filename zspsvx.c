#include "rb_lapack.h"

extern VOID zspsvx_(char *fact, char *uplo, integer *n, integer *nrhs, doublecomplex *ap, doublecomplex *afp, integer *ipiv, doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx, doublereal *rcond, doublereal *ferr, doublereal *berr, doublecomplex *work, doublereal *rwork, integer *info);

static VALUE
rb_zspsvx(int argc, VALUE *argv, VALUE self){
  VALUE rb_fact;
  char fact; 
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_ap;
  doublecomplex *ap; 
  VALUE rb_afp;
  doublecomplex *afp; 
  VALUE rb_ipiv;
  integer *ipiv; 
  VALUE rb_b;
  doublecomplex *b; 
  VALUE rb_x;
  doublecomplex *x; 
  VALUE rb_rcond;
  doublereal rcond; 
  VALUE rb_ferr;
  doublereal *ferr; 
  VALUE rb_berr;
  doublereal *berr; 
  VALUE rb_info;
  integer info; 
  VALUE rb_afp_out__;
  doublecomplex *afp_out__;
  VALUE rb_ipiv_out__;
  integer *ipiv_out__;
  doublecomplex *work;
  doublereal *rwork;

  integer n;
  integer ldb;
  integer nrhs;
  integer ldx;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  x, rcond, ferr, berr, info, afp, ipiv = NumRu::Lapack.zspsvx( fact, uplo, ap, afp, ipiv, b)\n    or\n  NumRu::Lapack.zspsvx  # print help\n\n\nFORTRAN MANUAL\n\n");
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
  if (NA_TYPE(rb_b) != NA_DCOMPLEX)
    rb_b = na_change_type(rb_b, NA_DCOMPLEX);
  b = NA_PTR_TYPE(rb_b, doublecomplex*);
  if (!NA_IsNArray(rb_afp))
    rb_raise(rb_eArgError, "afp (4th argument) must be NArray");
  if (NA_RANK(rb_afp) != 1)
    rb_raise(rb_eArgError, "rank of afp (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_afp) != (n*(n+1)/2))
    rb_raise(rb_eRuntimeError, "shape 0 of afp must be %d", n*(n+1)/2);
  if (NA_TYPE(rb_afp) != NA_DCOMPLEX)
    rb_afp = na_change_type(rb_afp, NA_DCOMPLEX);
  afp = NA_PTR_TYPE(rb_afp, doublecomplex*);
  if (!NA_IsNArray(rb_ap))
    rb_raise(rb_eArgError, "ap (3th argument) must be NArray");
  if (NA_RANK(rb_ap) != 1)
    rb_raise(rb_eArgError, "rank of ap (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_ap) != (n*(n+1)/2))
    rb_raise(rb_eRuntimeError, "shape 0 of ap must be %d", n*(n+1)/2);
  if (NA_TYPE(rb_ap) != NA_DCOMPLEX)
    rb_ap = na_change_type(rb_ap, NA_DCOMPLEX);
  ap = NA_PTR_TYPE(rb_ap, doublecomplex*);
  ldx = MAX(1,n);
  {
    int shape[2];
    shape[0] = ldx;
    shape[1] = nrhs;
    rb_x = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  x = NA_PTR_TYPE(rb_x, doublecomplex*);
  {
    int shape[1];
    shape[0] = nrhs;
    rb_ferr = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  ferr = NA_PTR_TYPE(rb_ferr, doublereal*);
  {
    int shape[1];
    shape[0] = nrhs;
    rb_berr = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  berr = NA_PTR_TYPE(rb_berr, doublereal*);
  {
    int shape[1];
    shape[0] = n*(n+1)/2;
    rb_afp_out__ = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  afp_out__ = NA_PTR_TYPE(rb_afp_out__, doublecomplex*);
  MEMCPY(afp_out__, afp, doublecomplex, NA_TOTAL(rb_afp));
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
  work = ALLOC_N(doublecomplex, (2*n));
  rwork = ALLOC_N(doublereal, (n));

  zspsvx_(&fact, &uplo, &n, &nrhs, ap, afp, ipiv, b, &ldb, x, &ldx, &rcond, ferr, berr, work, rwork, &info);

  free(work);
  free(rwork);
  rb_rcond = rb_float_new((double)rcond);
  rb_info = INT2NUM(info);
  return rb_ary_new3(7, rb_x, rb_rcond, rb_ferr, rb_berr, rb_info, rb_afp, rb_ipiv);
}

void
init_lapack_zspsvx(VALUE mLapack){
  rb_define_module_function(mLapack, "zspsvx", rb_zspsvx, -1);
}
