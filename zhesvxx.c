#include "rb_lapack.h"

extern VOID zhesvxx_(char *fact, char *uplo, integer *n, integer *nrhs, doublecomplex *a, integer *lda, doublecomplex *af, integer *ldaf, integer *ipiv, char *equed, doublereal *s, doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx, doublereal *rcond, doublereal *rpvgrw, doublereal *berr, integer *n_err_bnds, doublereal *err_bnds_norm, doublereal *err_bnds_comp, integer *nparams, doublereal *params, doublecomplex *work, doublereal *rwork, integer *info);

static VALUE
rb_zhesvxx(int argc, VALUE *argv, VALUE self){
  VALUE rb_fact;
  char fact; 
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_a;
  doublecomplex *a; 
  VALUE rb_af;
  doublecomplex *af; 
  VALUE rb_ipiv;
  integer *ipiv; 
  VALUE rb_equed;
  char equed; 
  VALUE rb_s;
  doublereal *s; 
  VALUE rb_b;
  doublecomplex *b; 
  VALUE rb_params;
  doublereal *params; 
  VALUE rb_x;
  doublecomplex *x; 
  VALUE rb_rcond;
  doublereal rcond; 
  VALUE rb_rpvgrw;
  doublereal rpvgrw; 
  VALUE rb_berr;
  doublereal *berr; 
  VALUE rb_err_bnds_norm;
  doublereal *err_bnds_norm; 
  VALUE rb_err_bnds_comp;
  doublereal *err_bnds_comp; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  doublecomplex *a_out__;
  VALUE rb_af_out__;
  doublecomplex *af_out__;
  VALUE rb_ipiv_out__;
  integer *ipiv_out__;
  VALUE rb_s_out__;
  doublereal *s_out__;
  VALUE rb_b_out__;
  doublecomplex *b_out__;
  VALUE rb_params_out__;
  doublereal *params_out__;
  doublecomplex *work;
  doublereal *rwork;

  integer lda;
  integer n;
  integer ldaf;
  integer ldb;
  integer nrhs;
  integer nparams;
  integer ldx;
  integer n_err_bnds;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  x, rcond, rpvgrw, berr, err_bnds_norm, err_bnds_comp, info, a, af, ipiv, equed, s, b, params = NumRu::Lapack.zhesvxx( fact, uplo, a, af, ipiv, equed, s, b, params)\n    or\n  NumRu::Lapack.zhesvxx  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 9)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 9)", argc);
  rb_fact = argv[0];
  rb_uplo = argv[1];
  rb_a = argv[2];
  rb_af = argv[3];
  rb_ipiv = argv[4];
  rb_equed = argv[5];
  rb_s = argv[6];
  rb_b = argv[7];
  rb_params = argv[8];

  uplo = StringValueCStr(rb_uplo)[0];
  if (!NA_IsNArray(rb_ipiv))
    rb_raise(rb_eArgError, "ipiv (5th argument) must be NArray");
  if (NA_RANK(rb_ipiv) != 1)
    rb_raise(rb_eArgError, "rank of ipiv (5th argument) must be %d", 1);
  n = NA_SHAPE0(rb_ipiv);
  if (NA_TYPE(rb_ipiv) != NA_LINT)
    rb_ipiv = na_change_type(rb_ipiv, NA_LINT);
  ipiv = NA_PTR_TYPE(rb_ipiv, integer*);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  if (NA_SHAPE1(rb_a) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of a must be the same as shape 0 of ipiv");
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_DCOMPLEX)
    rb_a = na_change_type(rb_a, NA_DCOMPLEX);
  a = NA_PTR_TYPE(rb_a, doublecomplex*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (8th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (8th argument) must be %d", 2);
  nrhs = NA_SHAPE1(rb_b);
  ldb = NA_SHAPE0(rb_b);
  if (NA_TYPE(rb_b) != NA_DCOMPLEX)
    rb_b = na_change_type(rb_b, NA_DCOMPLEX);
  b = NA_PTR_TYPE(rb_b, doublecomplex*);
  fact = StringValueCStr(rb_fact)[0];
  equed = StringValueCStr(rb_equed)[0];
  n_err_bnds = 3;
  if (!NA_IsNArray(rb_params))
    rb_raise(rb_eArgError, "params (9th argument) must be NArray");
  if (NA_RANK(rb_params) != 1)
    rb_raise(rb_eArgError, "rank of params (9th argument) must be %d", 1);
  nparams = NA_SHAPE0(rb_params);
  if (NA_TYPE(rb_params) != NA_DFLOAT)
    rb_params = na_change_type(rb_params, NA_DFLOAT);
  params = NA_PTR_TYPE(rb_params, doublereal*);
  if (!NA_IsNArray(rb_s))
    rb_raise(rb_eArgError, "s (7th argument) must be NArray");
  if (NA_RANK(rb_s) != 1)
    rb_raise(rb_eArgError, "rank of s (7th argument) must be %d", 1);
  if (NA_SHAPE0(rb_s) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of s must be the same as shape 0 of ipiv");
  if (NA_TYPE(rb_s) != NA_DFLOAT)
    rb_s = na_change_type(rb_s, NA_DFLOAT);
  s = NA_PTR_TYPE(rb_s, doublereal*);
  if (!NA_IsNArray(rb_af))
    rb_raise(rb_eArgError, "af (4th argument) must be NArray");
  if (NA_RANK(rb_af) != 2)
    rb_raise(rb_eArgError, "rank of af (4th argument) must be %d", 2);
  if (NA_SHAPE1(rb_af) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of af must be the same as shape 0 of ipiv");
  ldaf = NA_SHAPE0(rb_af);
  if (NA_TYPE(rb_af) != NA_DCOMPLEX)
    rb_af = na_change_type(rb_af, NA_DCOMPLEX);
  af = NA_PTR_TYPE(rb_af, doublecomplex*);
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
    rb_berr = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  berr = NA_PTR_TYPE(rb_berr, doublereal*);
  {
    int shape[2];
    shape[0] = nrhs;
    shape[1] = n_err_bnds;
    rb_err_bnds_norm = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  err_bnds_norm = NA_PTR_TYPE(rb_err_bnds_norm, doublereal*);
  {
    int shape[2];
    shape[0] = nrhs;
    shape[1] = n_err_bnds;
    rb_err_bnds_comp = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  err_bnds_comp = NA_PTR_TYPE(rb_err_bnds_comp, doublereal*);
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
  {
    int shape[2];
    shape[0] = ldaf;
    shape[1] = n;
    rb_af_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  af_out__ = NA_PTR_TYPE(rb_af_out__, doublecomplex*);
  MEMCPY(af_out__, af, doublecomplex, NA_TOTAL(rb_af));
  rb_af = rb_af_out__;
  af = af_out__;
  {
    int shape[1];
    shape[0] = n;
    rb_ipiv_out__ = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  ipiv_out__ = NA_PTR_TYPE(rb_ipiv_out__, integer*);
  MEMCPY(ipiv_out__, ipiv, integer, NA_TOTAL(rb_ipiv));
  rb_ipiv = rb_ipiv_out__;
  ipiv = ipiv_out__;
  {
    int shape[1];
    shape[0] = n;
    rb_s_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  s_out__ = NA_PTR_TYPE(rb_s_out__, doublereal*);
  MEMCPY(s_out__, s, doublereal, NA_TOTAL(rb_s));
  rb_s = rb_s_out__;
  s = s_out__;
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = nrhs;
    rb_b_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rb_b_out__, doublecomplex*);
  MEMCPY(b_out__, b, doublecomplex, NA_TOTAL(rb_b));
  rb_b = rb_b_out__;
  b = b_out__;
  {
    int shape[1];
    shape[0] = nparams;
    rb_params_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  params_out__ = NA_PTR_TYPE(rb_params_out__, doublereal*);
  MEMCPY(params_out__, params, doublereal, NA_TOTAL(rb_params));
  rb_params = rb_params_out__;
  params = params_out__;
  work = ALLOC_N(doublecomplex, (2*n));
  rwork = ALLOC_N(doublereal, (2*n));

  zhesvxx_(&fact, &uplo, &n, &nrhs, a, &lda, af, &ldaf, ipiv, &equed, s, b, &ldb, x, &ldx, &rcond, &rpvgrw, berr, &n_err_bnds, err_bnds_norm, err_bnds_comp, &nparams, params, work, rwork, &info);

  free(work);
  free(rwork);
  rb_rcond = rb_float_new((double)rcond);
  rb_rpvgrw = rb_float_new((double)rpvgrw);
  rb_info = INT2NUM(info);
  rb_equed = rb_str_new(&equed,1);
  return rb_ary_new3(14, rb_x, rb_rcond, rb_rpvgrw, rb_berr, rb_err_bnds_norm, rb_err_bnds_comp, rb_info, rb_a, rb_af, rb_ipiv, rb_equed, rb_s, rb_b, rb_params);
}

void
init_lapack_zhesvxx(VALUE mLapack){
  rb_define_module_function(mLapack, "zhesvxx", rb_zhesvxx, -1);
}
