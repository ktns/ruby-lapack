#include "rb_lapack.h"

extern VOID dla_porfsx_extended_(integer *prec_type, char *uplo, integer *n, integer *nrhs, doublereal *a, integer *lda, doublereal *af, integer *ldaf, logical *colequ, doublereal *c, doublereal *b, integer *ldb, doublereal *y, integer *ldy, doublereal *berr_out, integer *n_norms, doublereal *err_bnds_norm, doublereal *err_bnds_comp, doublereal *res, doublereal *ayb, doublereal *dy, doublereal *y_tail, doublereal *rcond, integer *ithresh, doublereal *rthresh, doublereal *dz_ub, logical *ignore_cwise, integer *info);

static VALUE
rb_dla_porfsx_extended(int argc, VALUE *argv, VALUE self){
  VALUE rb_prec_type;
  integer prec_type; 
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_a;
  doublereal *a; 
  VALUE rb_af;
  doublereal *af; 
  VALUE rb_colequ;
  logical colequ; 
  VALUE rb_c;
  doublereal *c; 
  VALUE rb_b;
  doublereal *b; 
  VALUE rb_y;
  doublereal *y; 
  VALUE rb_n_norms;
  integer n_norms; 
  VALUE rb_err_bnds_norm;
  doublereal *err_bnds_norm; 
  VALUE rb_err_bnds_comp;
  doublereal *err_bnds_comp; 
  VALUE rb_res;
  doublereal *res; 
  VALUE rb_ayb;
  doublereal *ayb; 
  VALUE rb_dy;
  doublereal *dy; 
  VALUE rb_y_tail;
  doublereal *y_tail; 
  VALUE rb_rcond;
  doublereal rcond; 
  VALUE rb_ithresh;
  integer ithresh; 
  VALUE rb_rthresh;
  doublereal rthresh; 
  VALUE rb_dz_ub;
  doublereal dz_ub; 
  VALUE rb_ignore_cwise;
  logical ignore_cwise; 
  VALUE rb_berr_out;
  doublereal *berr_out; 
  VALUE rb_info;
  integer info; 
  VALUE rb_y_out__;
  doublereal *y_out__;
  VALUE rb_err_bnds_norm_out__;
  doublereal *err_bnds_norm_out__;
  VALUE rb_err_bnds_comp_out__;
  doublereal *err_bnds_comp_out__;

  integer lda;
  integer n;
  integer ldaf;
  integer ldb;
  integer nrhs;
  integer ldy;
  integer n_err_bnds;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  berr_out, info, y, err_bnds_norm, err_bnds_comp = NumRu::Lapack.dla_porfsx_extended( prec_type, uplo, a, af, colequ, c, b, y, n_norms, err_bnds_norm, err_bnds_comp, res, ayb, dy, y_tail, rcond, ithresh, rthresh, dz_ub, ignore_cwise)\n    or\n  NumRu::Lapack.dla_porfsx_extended  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 20)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 20)", argc);
  rb_prec_type = argv[0];
  rb_uplo = argv[1];
  rb_a = argv[2];
  rb_af = argv[3];
  rb_colequ = argv[4];
  rb_c = argv[5];
  rb_b = argv[6];
  rb_y = argv[7];
  rb_n_norms = argv[8];
  rb_err_bnds_norm = argv[9];
  rb_err_bnds_comp = argv[10];
  rb_res = argv[11];
  rb_ayb = argv[12];
  rb_dy = argv[13];
  rb_y_tail = argv[14];
  rb_rcond = argv[15];
  rb_ithresh = argv[16];
  rb_rthresh = argv[17];
  rb_dz_ub = argv[18];
  rb_ignore_cwise = argv[19];

  rcond = NUM2DBL(rb_rcond);
  if (!NA_IsNArray(rb_res))
    rb_raise(rb_eArgError, "res (12th argument) must be NArray");
  if (NA_RANK(rb_res) != 1)
    rb_raise(rb_eArgError, "rank of res (12th argument) must be %d", 1);
  n = NA_SHAPE0(rb_res);
  if (NA_TYPE(rb_res) != NA_DFLOAT)
    rb_res = na_change_type(rb_res, NA_DFLOAT);
  res = NA_PTR_TYPE(rb_res, doublereal*);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  if (NA_SHAPE1(rb_a) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of a must be the same as shape 0 of res");
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_DFLOAT)
    rb_a = na_change_type(rb_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rb_a, doublereal*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (7th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (7th argument) must be %d", 2);
  nrhs = NA_SHAPE1(rb_b);
  ldb = NA_SHAPE0(rb_b);
  if (NA_TYPE(rb_b) != NA_DFLOAT)
    rb_b = na_change_type(rb_b, NA_DFLOAT);
  b = NA_PTR_TYPE(rb_b, doublereal*);
  colequ = (rb_colequ == Qtrue);
  if (!NA_IsNArray(rb_y))
    rb_raise(rb_eArgError, "y (8th argument) must be NArray");
  if (NA_RANK(rb_y) != 2)
    rb_raise(rb_eArgError, "rank of y (8th argument) must be %d", 2);
  if (NA_SHAPE1(rb_y) != nrhs)
    rb_raise(rb_eRuntimeError, "shape 1 of y must be the same as shape 1 of b");
  ldy = NA_SHAPE0(rb_y);
  if (NA_TYPE(rb_y) != NA_DFLOAT)
    rb_y = na_change_type(rb_y, NA_DFLOAT);
  y = NA_PTR_TYPE(rb_y, doublereal*);
  if (!NA_IsNArray(rb_c))
    rb_raise(rb_eArgError, "c (6th argument) must be NArray");
  if (NA_RANK(rb_c) != 1)
    rb_raise(rb_eArgError, "rank of c (6th argument) must be %d", 1);
  if (NA_SHAPE0(rb_c) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of c must be the same as shape 0 of res");
  if (NA_TYPE(rb_c) != NA_DFLOAT)
    rb_c = na_change_type(rb_c, NA_DFLOAT);
  c = NA_PTR_TYPE(rb_c, doublereal*);
  dz_ub = NUM2DBL(rb_dz_ub);
  if (!NA_IsNArray(rb_y_tail))
    rb_raise(rb_eArgError, "y_tail (15th argument) must be NArray");
  if (NA_RANK(rb_y_tail) != 1)
    rb_raise(rb_eArgError, "rank of y_tail (15th argument) must be %d", 1);
  if (NA_SHAPE0(rb_y_tail) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of y_tail must be the same as shape 0 of res");
  if (NA_TYPE(rb_y_tail) != NA_DFLOAT)
    rb_y_tail = na_change_type(rb_y_tail, NA_DFLOAT);
  y_tail = NA_PTR_TYPE(rb_y_tail, doublereal*);
  if (!NA_IsNArray(rb_af))
    rb_raise(rb_eArgError, "af (4th argument) must be NArray");
  if (NA_RANK(rb_af) != 2)
    rb_raise(rb_eArgError, "rank of af (4th argument) must be %d", 2);
  if (NA_SHAPE1(rb_af) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of af must be the same as shape 0 of res");
  ldaf = NA_SHAPE0(rb_af);
  if (NA_TYPE(rb_af) != NA_DFLOAT)
    rb_af = na_change_type(rb_af, NA_DFLOAT);
  af = NA_PTR_TYPE(rb_af, doublereal*);
  prec_type = NUM2INT(rb_prec_type);
  rthresh = NUM2DBL(rb_rthresh);
  ithresh = NUM2INT(rb_ithresh);
  if (!NA_IsNArray(rb_err_bnds_comp))
    rb_raise(rb_eArgError, "err_bnds_comp (11th argument) must be NArray");
  if (NA_RANK(rb_err_bnds_comp) != 2)
    rb_raise(rb_eArgError, "rank of err_bnds_comp (11th argument) must be %d", 2);
  n_err_bnds = NA_SHAPE1(rb_err_bnds_comp);
  if (NA_SHAPE0(rb_err_bnds_comp) != nrhs)
    rb_raise(rb_eRuntimeError, "shape 0 of err_bnds_comp must be the same as shape 1 of b");
  if (NA_TYPE(rb_err_bnds_comp) != NA_DFLOAT)
    rb_err_bnds_comp = na_change_type(rb_err_bnds_comp, NA_DFLOAT);
  err_bnds_comp = NA_PTR_TYPE(rb_err_bnds_comp, doublereal*);
  ignore_cwise = (rb_ignore_cwise == Qtrue);
  if (!NA_IsNArray(rb_dy))
    rb_raise(rb_eArgError, "dy (14th argument) must be NArray");
  if (NA_RANK(rb_dy) != 1)
    rb_raise(rb_eArgError, "rank of dy (14th argument) must be %d", 1);
  if (NA_SHAPE0(rb_dy) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of dy must be the same as shape 0 of res");
  if (NA_TYPE(rb_dy) != NA_DFLOAT)
    rb_dy = na_change_type(rb_dy, NA_DFLOAT);
  dy = NA_PTR_TYPE(rb_dy, doublereal*);
  if (!NA_IsNArray(rb_ayb))
    rb_raise(rb_eArgError, "ayb (13th argument) must be NArray");
  if (NA_RANK(rb_ayb) != 1)
    rb_raise(rb_eArgError, "rank of ayb (13th argument) must be %d", 1);
  if (NA_SHAPE0(rb_ayb) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of ayb must be the same as shape 0 of res");
  if (NA_TYPE(rb_ayb) != NA_DFLOAT)
    rb_ayb = na_change_type(rb_ayb, NA_DFLOAT);
  ayb = NA_PTR_TYPE(rb_ayb, doublereal*);
  if (!NA_IsNArray(rb_err_bnds_norm))
    rb_raise(rb_eArgError, "err_bnds_norm (10th argument) must be NArray");
  if (NA_RANK(rb_err_bnds_norm) != 2)
    rb_raise(rb_eArgError, "rank of err_bnds_norm (10th argument) must be %d", 2);
  if (NA_SHAPE1(rb_err_bnds_norm) != n_err_bnds)
    rb_raise(rb_eRuntimeError, "shape 1 of err_bnds_norm must be the same as shape 1 of err_bnds_comp");
  if (NA_SHAPE0(rb_err_bnds_norm) != nrhs)
    rb_raise(rb_eRuntimeError, "shape 0 of err_bnds_norm must be the same as shape 1 of b");
  if (NA_TYPE(rb_err_bnds_norm) != NA_DFLOAT)
    rb_err_bnds_norm = na_change_type(rb_err_bnds_norm, NA_DFLOAT);
  err_bnds_norm = NA_PTR_TYPE(rb_err_bnds_norm, doublereal*);
  n_norms = NUM2INT(rb_n_norms);
  uplo = StringValueCStr(rb_uplo)[0];
  {
    int shape[1];
    shape[0] = nrhs;
    rb_berr_out = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  berr_out = NA_PTR_TYPE(rb_berr_out, doublereal*);
  {
    int shape[2];
    shape[0] = ldy;
    shape[1] = nrhs;
    rb_y_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  y_out__ = NA_PTR_TYPE(rb_y_out__, doublereal*);
  MEMCPY(y_out__, y, doublereal, NA_TOTAL(rb_y));
  rb_y = rb_y_out__;
  y = y_out__;
  {
    int shape[2];
    shape[0] = nrhs;
    shape[1] = n_err_bnds;
    rb_err_bnds_norm_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  err_bnds_norm_out__ = NA_PTR_TYPE(rb_err_bnds_norm_out__, doublereal*);
  MEMCPY(err_bnds_norm_out__, err_bnds_norm, doublereal, NA_TOTAL(rb_err_bnds_norm));
  rb_err_bnds_norm = rb_err_bnds_norm_out__;
  err_bnds_norm = err_bnds_norm_out__;
  {
    int shape[2];
    shape[0] = nrhs;
    shape[1] = n_err_bnds;
    rb_err_bnds_comp_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  err_bnds_comp_out__ = NA_PTR_TYPE(rb_err_bnds_comp_out__, doublereal*);
  MEMCPY(err_bnds_comp_out__, err_bnds_comp, doublereal, NA_TOTAL(rb_err_bnds_comp));
  rb_err_bnds_comp = rb_err_bnds_comp_out__;
  err_bnds_comp = err_bnds_comp_out__;

  dla_porfsx_extended_(&prec_type, &uplo, &n, &nrhs, a, &lda, af, &ldaf, &colequ, c, b, &ldb, y, &ldy, berr_out, &n_norms, err_bnds_norm, err_bnds_comp, res, ayb, dy, y_tail, &rcond, &ithresh, &rthresh, &dz_ub, &ignore_cwise, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(5, rb_berr_out, rb_info, rb_y, rb_err_bnds_norm, rb_err_bnds_comp);
}

void
init_lapack_dla_porfsx_extended(VALUE mLapack){
  rb_define_module_function(mLapack, "dla_porfsx_extended", rb_dla_porfsx_extended, -1);
}
