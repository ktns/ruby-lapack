#include "rb_lapack.h"

extern VOID dgerfsx_(char *trans, char *equed, integer *n, integer *nrhs, doublereal *a, integer *lda, doublereal *af, integer *ldaf, integer *ipiv, doublereal *r, doublereal *c, doublereal *b, integer *ldb, doublereal *x, integer *ldx, doublereal *rcond, doublereal *berr, integer *n_err_bnds, doublereal *err_bnds_norm, doublereal *err_bnds_comp, integer *nparams, doublereal *params, doublereal *work, integer *iwork, integer *info);

static VALUE
rb_dgerfsx(int argc, VALUE *argv, VALUE self){
  VALUE rb_trans;
  char trans; 
  VALUE rb_equed;
  char equed; 
  VALUE rb_a;
  doublereal *a; 
  VALUE rb_af;
  doublereal *af; 
  VALUE rb_ipiv;
  integer *ipiv; 
  VALUE rb_r;
  doublereal *r; 
  VALUE rb_c;
  doublereal *c; 
  VALUE rb_b;
  doublereal *b; 
  VALUE rb_x;
  doublereal *x; 
  VALUE rb_params;
  doublereal *params; 
  VALUE rb_rcond;
  doublereal rcond; 
  VALUE rb_berr;
  doublereal *berr; 
  VALUE rb_err_bnds_norm;
  doublereal *err_bnds_norm; 
  VALUE rb_err_bnds_comp;
  doublereal *err_bnds_comp; 
  VALUE rb_info;
  integer info; 
  VALUE rb_x_out__;
  doublereal *x_out__;
  VALUE rb_params_out__;
  doublereal *params_out__;
  doublereal *work;
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
    printf("%s\n", "USAGE:\n  rcond, berr, err_bnds_norm, err_bnds_comp, info, x, params = NumRu::Lapack.dgerfsx( trans, equed, a, af, ipiv, r, c, b, x, params)\n    or\n  NumRu::Lapack.dgerfsx  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 10)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 10)", argc);
  rb_trans = argv[0];
  rb_equed = argv[1];
  rb_a = argv[2];
  rb_af = argv[3];
  rb_ipiv = argv[4];
  rb_r = argv[5];
  rb_c = argv[6];
  rb_b = argv[7];
  rb_x = argv[8];
  rb_params = argv[9];

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
  if (NA_TYPE(rb_a) != NA_DFLOAT)
    rb_a = na_change_type(rb_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rb_a, doublereal*);
  if (!NA_IsNArray(rb_x))
    rb_raise(rb_eArgError, "x (9th argument) must be NArray");
  if (NA_RANK(rb_x) != 2)
    rb_raise(rb_eArgError, "rank of x (9th argument) must be %d", 2);
  nrhs = NA_SHAPE1(rb_x);
  ldx = NA_SHAPE0(rb_x);
  if (NA_TYPE(rb_x) != NA_DFLOAT)
    rb_x = na_change_type(rb_x, NA_DFLOAT);
  x = NA_PTR_TYPE(rb_x, doublereal*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (8th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (8th argument) must be %d", 2);
  if (NA_SHAPE1(rb_b) != nrhs)
    rb_raise(rb_eRuntimeError, "shape 1 of b must be the same as shape 1 of x");
  ldb = NA_SHAPE0(rb_b);
  if (NA_TYPE(rb_b) != NA_DFLOAT)
    rb_b = na_change_type(rb_b, NA_DFLOAT);
  b = NA_PTR_TYPE(rb_b, doublereal*);
  if (!NA_IsNArray(rb_c))
    rb_raise(rb_eArgError, "c (7th argument) must be NArray");
  if (NA_RANK(rb_c) != 1)
    rb_raise(rb_eArgError, "rank of c (7th argument) must be %d", 1);
  if (NA_SHAPE0(rb_c) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of c must be the same as shape 0 of ipiv");
  if (NA_TYPE(rb_c) != NA_DFLOAT)
    rb_c = na_change_type(rb_c, NA_DFLOAT);
  c = NA_PTR_TYPE(rb_c, doublereal*);
  if (!NA_IsNArray(rb_af))
    rb_raise(rb_eArgError, "af (4th argument) must be NArray");
  if (NA_RANK(rb_af) != 2)
    rb_raise(rb_eArgError, "rank of af (4th argument) must be %d", 2);
  if (NA_SHAPE1(rb_af) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of af must be the same as shape 0 of ipiv");
  ldaf = NA_SHAPE0(rb_af);
  if (NA_TYPE(rb_af) != NA_DFLOAT)
    rb_af = na_change_type(rb_af, NA_DFLOAT);
  af = NA_PTR_TYPE(rb_af, doublereal*);
  if (!NA_IsNArray(rb_r))
    rb_raise(rb_eArgError, "r (6th argument) must be NArray");
  if (NA_RANK(rb_r) != 1)
    rb_raise(rb_eArgError, "rank of r (6th argument) must be %d", 1);
  if (NA_SHAPE0(rb_r) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of r must be the same as shape 0 of ipiv");
  if (NA_TYPE(rb_r) != NA_DFLOAT)
    rb_r = na_change_type(rb_r, NA_DFLOAT);
  r = NA_PTR_TYPE(rb_r, doublereal*);
  n_err_bnds = 3;
  if (!NA_IsNArray(rb_params))
    rb_raise(rb_eArgError, "params (10th argument) must be NArray");
  if (NA_RANK(rb_params) != 1)
    rb_raise(rb_eArgError, "rank of params (10th argument) must be %d", 1);
  nparams = NA_SHAPE0(rb_params);
  if (NA_TYPE(rb_params) != NA_DFLOAT)
    rb_params = na_change_type(rb_params, NA_DFLOAT);
  params = NA_PTR_TYPE(rb_params, doublereal*);
  equed = StringValueCStr(rb_equed)[0];
  trans = StringValueCStr(rb_trans)[0];
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
    shape[0] = ldx;
    shape[1] = nrhs;
    rb_x_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  x_out__ = NA_PTR_TYPE(rb_x_out__, doublereal*);
  MEMCPY(x_out__, x, doublereal, NA_TOTAL(rb_x));
  rb_x = rb_x_out__;
  x = x_out__;
  {
    int shape[1];
    shape[0] = nparams;
    rb_params_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  params_out__ = NA_PTR_TYPE(rb_params_out__, doublereal*);
  MEMCPY(params_out__, params, doublereal, NA_TOTAL(rb_params));
  rb_params = rb_params_out__;
  params = params_out__;
  work = ALLOC_N(doublereal, (4*n));
  iwork = ALLOC_N(integer, (n));

  dgerfsx_(&trans, &equed, &n, &nrhs, a, &lda, af, &ldaf, ipiv, r, c, b, &ldb, x, &ldx, &rcond, berr, &n_err_bnds, err_bnds_norm, err_bnds_comp, &nparams, params, work, iwork, &info);

  free(work);
  free(iwork);
  rb_rcond = rb_float_new((double)rcond);
  rb_info = INT2NUM(info);
  return rb_ary_new3(7, rb_rcond, rb_berr, rb_err_bnds_norm, rb_err_bnds_comp, rb_info, rb_x, rb_params);
}

void
init_lapack_dgerfsx(VALUE mLapack){
  rb_define_module_function(mLapack, "dgerfsx", rb_dgerfsx, -1);
}
