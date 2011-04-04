#include "rb_lapack.h"

extern VOID cgesvxx_(char *fact, char *trans, integer *n, integer *nrhs, complex *a, integer *lda, complex *af, integer *ldaf, integer *ipiv, char *equed, real *r, real *c, complex *b, integer *ldb, complex *x, integer *ldx, real *rcond, real *rpvgrw, real *berr, integer *n_err_bnds, real *err_bnds_norm, real *err_bnds_comp, integer *nparams, real *params, complex *work, real *rwork, integer *info);

static VALUE
rb_cgesvxx(int argc, VALUE *argv, VALUE self){
  VALUE rb_fact;
  char fact; 
  VALUE rb_trans;
  char trans; 
  VALUE rb_a;
  complex *a; 
  VALUE rb_af;
  complex *af; 
  VALUE rb_ipiv;
  integer *ipiv; 
  VALUE rb_equed;
  char equed; 
  VALUE rb_r;
  real *r; 
  VALUE rb_c;
  real *c; 
  VALUE rb_b;
  complex *b; 
  VALUE rb_params;
  real *params; 
  VALUE rb_x;
  complex *x; 
  VALUE rb_rcond;
  real rcond; 
  VALUE rb_rpvgrw;
  real rpvgrw; 
  VALUE rb_berr;
  real *berr; 
  VALUE rb_err_bnds_norm;
  real *err_bnds_norm; 
  VALUE rb_err_bnds_comp;
  real *err_bnds_comp; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  complex *a_out__;
  VALUE rb_af_out__;
  complex *af_out__;
  VALUE rb_ipiv_out__;
  integer *ipiv_out__;
  VALUE rb_r_out__;
  real *r_out__;
  VALUE rb_c_out__;
  real *c_out__;
  VALUE rb_b_out__;
  complex *b_out__;
  VALUE rb_params_out__;
  real *params_out__;
  complex *work;
  real *rwork;

  integer lda;
  integer n;
  integer ldaf;
  integer ldb;
  integer nrhs;
  integer nparams;
  integer ldx;
  integer n_err_bnds;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  x, rcond, rpvgrw, berr, err_bnds_norm, err_bnds_comp, info, a, af, ipiv, equed, r, c, b, params = NumRu::Lapack.cgesvxx( fact, trans, a, af, ipiv, equed, r, c, b, params)\n    or\n  NumRu::Lapack.cgesvxx  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 10)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 10)", argc);
  rb_fact = argv[0];
  rb_trans = argv[1];
  rb_a = argv[2];
  rb_af = argv[3];
  rb_ipiv = argv[4];
  rb_equed = argv[5];
  rb_r = argv[6];
  rb_c = argv[7];
  rb_b = argv[8];
  rb_params = argv[9];

  trans = StringValueCStr(rb_trans)[0];
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
  if (NA_TYPE(rb_a) != NA_SCOMPLEX)
    rb_a = na_change_type(rb_a, NA_SCOMPLEX);
  a = NA_PTR_TYPE(rb_a, complex*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (9th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (9th argument) must be %d", 2);
  nrhs = NA_SHAPE1(rb_b);
  ldb = NA_SHAPE0(rb_b);
  if (NA_TYPE(rb_b) != NA_SCOMPLEX)
    rb_b = na_change_type(rb_b, NA_SCOMPLEX);
  b = NA_PTR_TYPE(rb_b, complex*);
  if (!NA_IsNArray(rb_c))
    rb_raise(rb_eArgError, "c (8th argument) must be NArray");
  if (NA_RANK(rb_c) != 1)
    rb_raise(rb_eArgError, "rank of c (8th argument) must be %d", 1);
  if (NA_SHAPE0(rb_c) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of c must be the same as shape 0 of ipiv");
  if (NA_TYPE(rb_c) != NA_SFLOAT)
    rb_c = na_change_type(rb_c, NA_SFLOAT);
  c = NA_PTR_TYPE(rb_c, real*);
  fact = StringValueCStr(rb_fact)[0];
  equed = StringValueCStr(rb_equed)[0];
  if (!NA_IsNArray(rb_r))
    rb_raise(rb_eArgError, "r (7th argument) must be NArray");
  if (NA_RANK(rb_r) != 1)
    rb_raise(rb_eArgError, "rank of r (7th argument) must be %d", 1);
  if (NA_SHAPE0(rb_r) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of r must be the same as shape 0 of ipiv");
  if (NA_TYPE(rb_r) != NA_SFLOAT)
    rb_r = na_change_type(rb_r, NA_SFLOAT);
  r = NA_PTR_TYPE(rb_r, real*);
  n_err_bnds = 3;
  if (!NA_IsNArray(rb_params))
    rb_raise(rb_eArgError, "params (10th argument) must be NArray");
  if (NA_RANK(rb_params) != 1)
    rb_raise(rb_eArgError, "rank of params (10th argument) must be %d", 1);
  nparams = NA_SHAPE0(rb_params);
  if (NA_TYPE(rb_params) != NA_SFLOAT)
    rb_params = na_change_type(rb_params, NA_SFLOAT);
  params = NA_PTR_TYPE(rb_params, real*);
  if (!NA_IsNArray(rb_af))
    rb_raise(rb_eArgError, "af (4th argument) must be NArray");
  if (NA_RANK(rb_af) != 2)
    rb_raise(rb_eArgError, "rank of af (4th argument) must be %d", 2);
  if (NA_SHAPE1(rb_af) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of af must be the same as shape 0 of ipiv");
  ldaf = NA_SHAPE0(rb_af);
  if (NA_TYPE(rb_af) != NA_SCOMPLEX)
    rb_af = na_change_type(rb_af, NA_SCOMPLEX);
  af = NA_PTR_TYPE(rb_af, complex*);
  ldx = n;
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
    rb_berr = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  berr = NA_PTR_TYPE(rb_berr, real*);
  {
    int shape[2];
    shape[0] = nrhs;
    shape[1] = n_err_bnds;
    rb_err_bnds_norm = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  err_bnds_norm = NA_PTR_TYPE(rb_err_bnds_norm, real*);
  {
    int shape[2];
    shape[0] = nrhs;
    shape[1] = n_err_bnds;
    rb_err_bnds_comp = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  err_bnds_comp = NA_PTR_TYPE(rb_err_bnds_comp, real*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rb_a_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rb_a_out__, complex*);
  MEMCPY(a_out__, a, complex, NA_TOTAL(rb_a));
  rb_a = rb_a_out__;
  a = a_out__;
  {
    int shape[2];
    shape[0] = ldaf;
    shape[1] = n;
    rb_af_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  af_out__ = NA_PTR_TYPE(rb_af_out__, complex*);
  MEMCPY(af_out__, af, complex, NA_TOTAL(rb_af));
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
    rb_r_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  r_out__ = NA_PTR_TYPE(rb_r_out__, real*);
  MEMCPY(r_out__, r, real, NA_TOTAL(rb_r));
  rb_r = rb_r_out__;
  r = r_out__;
  {
    int shape[1];
    shape[0] = n;
    rb_c_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  c_out__ = NA_PTR_TYPE(rb_c_out__, real*);
  MEMCPY(c_out__, c, real, NA_TOTAL(rb_c));
  rb_c = rb_c_out__;
  c = c_out__;
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = nrhs;
    rb_b_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rb_b_out__, complex*);
  MEMCPY(b_out__, b, complex, NA_TOTAL(rb_b));
  rb_b = rb_b_out__;
  b = b_out__;
  {
    int shape[1];
    shape[0] = nparams;
    rb_params_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  params_out__ = NA_PTR_TYPE(rb_params_out__, real*);
  MEMCPY(params_out__, params, real, NA_TOTAL(rb_params));
  rb_params = rb_params_out__;
  params = params_out__;
  work = ALLOC_N(complex, (2*n));
  rwork = ALLOC_N(real, (2*n));

  cgesvxx_(&fact, &trans, &n, &nrhs, a, &lda, af, &ldaf, ipiv, &equed, r, c, b, &ldb, x, &ldx, &rcond, &rpvgrw, berr, &n_err_bnds, err_bnds_norm, err_bnds_comp, &nparams, params, work, rwork, &info);

  free(work);
  free(rwork);
  rb_rcond = rb_float_new((double)rcond);
  rb_rpvgrw = rb_float_new((double)rpvgrw);
  rb_info = INT2NUM(info);
  rb_equed = rb_str_new(&equed,1);
  return rb_ary_new3(15, rb_x, rb_rcond, rb_rpvgrw, rb_berr, rb_err_bnds_norm, rb_err_bnds_comp, rb_info, rb_a, rb_af, rb_ipiv, rb_equed, rb_r, rb_c, rb_b, rb_params);
}

void
init_lapack_cgesvxx(VALUE mLapack){
  rb_define_module_function(mLapack, "cgesvxx", rb_cgesvxx, -1);
}
