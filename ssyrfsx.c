#include "rb_lapack.h"

extern VOID ssyrfsx_(char *uplo, char *equed, integer *n, integer *nrhs, real *a, integer *lda, real *af, integer *ldaf, integer *ipiv, real *s, real *b, integer *ldb, real *x, integer *ldx, real *rcond, real *berr, integer *n_err_bnds, real *err_bnds_norm, real *err_bnds_comp, integer *nparams, real *params, real *work, integer *iwork, integer *info);

static VALUE
rb_ssyrfsx(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_equed;
  char equed; 
  VALUE rb_a;
  real *a; 
  VALUE rb_af;
  real *af; 
  VALUE rb_ipiv;
  integer *ipiv; 
  VALUE rb_s;
  real *s; 
  VALUE rb_b;
  real *b; 
  VALUE rb_x;
  real *x; 
  VALUE rb_params;
  real *params; 
  VALUE rb_rcond;
  real rcond; 
  VALUE rb_berr;
  real *berr; 
  VALUE rb_err_bnds_norm;
  real *err_bnds_norm; 
  VALUE rb_err_bnds_comp;
  real *err_bnds_comp; 
  VALUE rb_info;
  integer info; 
  VALUE rb_s_out__;
  real *s_out__;
  VALUE rb_x_out__;
  real *x_out__;
  VALUE rb_params_out__;
  real *params_out__;
  real *work;
  integer *iwork;

  integer lda;
  integer n;
  integer ldaf;
  integer ldb;
  integer nrhs;
  integer ldx;
  integer nparams;
  integer n_err_bnds;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  rcond, berr, err_bnds_norm, err_bnds_comp, info, s, x, params = NumRu::Lapack.ssyrfsx( uplo, equed, a, af, ipiv, s, b, x, params)\n    or\n  NumRu::Lapack.ssyrfsx  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 9)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 9)", argc);
  rb_uplo = argv[0];
  rb_equed = argv[1];
  rb_a = argv[2];
  rb_af = argv[3];
  rb_ipiv = argv[4];
  rb_s = argv[5];
  rb_b = argv[6];
  rb_x = argv[7];
  rb_params = argv[8];

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
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);
  if (!NA_IsNArray(rb_x))
    rb_raise(rb_eArgError, "x (8th argument) must be NArray");
  if (NA_RANK(rb_x) != 2)
    rb_raise(rb_eArgError, "rank of x (8th argument) must be %d", 2);
  nrhs = NA_SHAPE1(rb_x);
  ldx = NA_SHAPE0(rb_x);
  if (NA_TYPE(rb_x) != NA_SFLOAT)
    rb_x = na_change_type(rb_x, NA_SFLOAT);
  x = NA_PTR_TYPE(rb_x, real*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (7th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (7th argument) must be %d", 2);
  if (NA_SHAPE1(rb_b) != nrhs)
    rb_raise(rb_eRuntimeError, "shape 1 of b must be the same as shape 1 of x");
  ldb = NA_SHAPE0(rb_b);
  if (NA_TYPE(rb_b) != NA_SFLOAT)
    rb_b = na_change_type(rb_b, NA_SFLOAT);
  b = NA_PTR_TYPE(rb_b, real*);
  if (!NA_IsNArray(rb_af))
    rb_raise(rb_eArgError, "af (4th argument) must be NArray");
  if (NA_RANK(rb_af) != 2)
    rb_raise(rb_eArgError, "rank of af (4th argument) must be %d", 2);
  if (NA_SHAPE1(rb_af) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of af must be the same as shape 0 of ipiv");
  ldaf = NA_SHAPE0(rb_af);
  if (NA_TYPE(rb_af) != NA_SFLOAT)
    rb_af = na_change_type(rb_af, NA_SFLOAT);
  af = NA_PTR_TYPE(rb_af, real*);
  n_err_bnds = 3;
  if (!NA_IsNArray(rb_params))
    rb_raise(rb_eArgError, "params (9th argument) must be NArray");
  if (NA_RANK(rb_params) != 1)
    rb_raise(rb_eArgError, "rank of params (9th argument) must be %d", 1);
  nparams = NA_SHAPE0(rb_params);
  if (NA_TYPE(rb_params) != NA_SFLOAT)
    rb_params = na_change_type(rb_params, NA_SFLOAT);
  params = NA_PTR_TYPE(rb_params, real*);
  if (!NA_IsNArray(rb_s))
    rb_raise(rb_eArgError, "s (6th argument) must be NArray");
  if (NA_RANK(rb_s) != 1)
    rb_raise(rb_eArgError, "rank of s (6th argument) must be %d", 1);
  if (NA_SHAPE0(rb_s) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of s must be the same as shape 0 of ipiv");
  if (NA_TYPE(rb_s) != NA_SFLOAT)
    rb_s = na_change_type(rb_s, NA_SFLOAT);
  s = NA_PTR_TYPE(rb_s, real*);
  equed = StringValueCStr(rb_equed)[0];
  uplo = StringValueCStr(rb_uplo)[0];
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
    int shape[1];
    shape[0] = n;
    rb_s_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  s_out__ = NA_PTR_TYPE(rb_s_out__, real*);
  MEMCPY(s_out__, s, real, NA_TOTAL(rb_s));
  rb_s = rb_s_out__;
  s = s_out__;
  {
    int shape[2];
    shape[0] = ldx;
    shape[1] = nrhs;
    rb_x_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  x_out__ = NA_PTR_TYPE(rb_x_out__, real*);
  MEMCPY(x_out__, x, real, NA_TOTAL(rb_x));
  rb_x = rb_x_out__;
  x = x_out__;
  {
    int shape[1];
    shape[0] = nparams;
    rb_params_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  params_out__ = NA_PTR_TYPE(rb_params_out__, real*);
  MEMCPY(params_out__, params, real, NA_TOTAL(rb_params));
  rb_params = rb_params_out__;
  params = params_out__;
  work = ALLOC_N(real, (4*n));
  iwork = ALLOC_N(integer, (n));

  ssyrfsx_(&uplo, &equed, &n, &nrhs, a, &lda, af, &ldaf, ipiv, s, b, &ldb, x, &ldx, &rcond, berr, &n_err_bnds, err_bnds_norm, err_bnds_comp, &nparams, params, work, iwork, &info);

  free(work);
  free(iwork);
  rb_rcond = rb_float_new((double)rcond);
  rb_info = INT2NUM(info);
  return rb_ary_new3(8, rb_rcond, rb_berr, rb_err_bnds_norm, rb_err_bnds_comp, rb_info, rb_s, rb_x, rb_params);
}

void
init_lapack_ssyrfsx(VALUE mLapack){
  rb_define_module_function(mLapack, "ssyrfsx", rb_ssyrfsx, -1);
}
