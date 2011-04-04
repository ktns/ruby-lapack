#include "rb_lapack.h"

extern VOID cgtrfs_(char *trans, integer *n, integer *nrhs, complex *dl, complex *d, complex *du, complex *dlf, complex *df, complex *duf, complex *du2, integer *ipiv, complex *b, integer *ldb, complex *x, integer *ldx, real *ferr, real *berr, complex *work, real *rwork, integer *info);

static VALUE
rb_cgtrfs(int argc, VALUE *argv, VALUE self){
  VALUE rb_trans;
  char trans; 
  VALUE rb_dl;
  complex *dl; 
  VALUE rb_d;
  complex *d; 
  VALUE rb_du;
  complex *du; 
  VALUE rb_dlf;
  complex *dlf; 
  VALUE rb_df;
  complex *df; 
  VALUE rb_duf;
  complex *duf; 
  VALUE rb_du2;
  complex *du2; 
  VALUE rb_ipiv;
  integer *ipiv; 
  VALUE rb_b;
  complex *b; 
  VALUE rb_x;
  complex *x; 
  VALUE rb_ferr;
  real *ferr; 
  VALUE rb_berr;
  real *berr; 
  VALUE rb_info;
  integer info; 
  VALUE rb_x_out__;
  complex *x_out__;
  complex *work;
  real *rwork;

  integer n;
  integer ldb;
  integer nrhs;
  integer ldx;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  ferr, berr, info, x = NumRu::Lapack.cgtrfs( trans, dl, d, du, dlf, df, duf, du2, ipiv, b, x)\n    or\n  NumRu::Lapack.cgtrfs  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 11)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 11)", argc);
  rb_trans = argv[0];
  rb_dl = argv[1];
  rb_d = argv[2];
  rb_du = argv[3];
  rb_dlf = argv[4];
  rb_df = argv[5];
  rb_duf = argv[6];
  rb_du2 = argv[7];
  rb_ipiv = argv[8];
  rb_b = argv[9];
  rb_x = argv[10];

  if (!NA_IsNArray(rb_ipiv))
    rb_raise(rb_eArgError, "ipiv (9th argument) must be NArray");
  if (NA_RANK(rb_ipiv) != 1)
    rb_raise(rb_eArgError, "rank of ipiv (9th argument) must be %d", 1);
  n = NA_SHAPE0(rb_ipiv);
  if (NA_TYPE(rb_ipiv) != NA_LINT)
    rb_ipiv = na_change_type(rb_ipiv, NA_LINT);
  ipiv = NA_PTR_TYPE(rb_ipiv, integer*);
  trans = StringValueCStr(rb_trans)[0];
  if (!NA_IsNArray(rb_x))
    rb_raise(rb_eArgError, "x (11th argument) must be NArray");
  if (NA_RANK(rb_x) != 2)
    rb_raise(rb_eArgError, "rank of x (11th argument) must be %d", 2);
  nrhs = NA_SHAPE1(rb_x);
  ldx = NA_SHAPE0(rb_x);
  if (NA_TYPE(rb_x) != NA_SCOMPLEX)
    rb_x = na_change_type(rb_x, NA_SCOMPLEX);
  x = NA_PTR_TYPE(rb_x, complex*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (10th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (10th argument) must be %d", 2);
  if (NA_SHAPE1(rb_b) != nrhs)
    rb_raise(rb_eRuntimeError, "shape 1 of b must be the same as shape 1 of x");
  ldb = NA_SHAPE0(rb_b);
  if (NA_TYPE(rb_b) != NA_SCOMPLEX)
    rb_b = na_change_type(rb_b, NA_SCOMPLEX);
  b = NA_PTR_TYPE(rb_b, complex*);
  if (!NA_IsNArray(rb_df))
    rb_raise(rb_eArgError, "df (6th argument) must be NArray");
  if (NA_RANK(rb_df) != 1)
    rb_raise(rb_eArgError, "rank of df (6th argument) must be %d", 1);
  if (NA_SHAPE0(rb_df) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of df must be the same as shape 0 of ipiv");
  if (NA_TYPE(rb_df) != NA_SCOMPLEX)
    rb_df = na_change_type(rb_df, NA_SCOMPLEX);
  df = NA_PTR_TYPE(rb_df, complex*);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (3th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_d) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of d must be the same as shape 0 of ipiv");
  if (NA_TYPE(rb_d) != NA_SCOMPLEX)
    rb_d = na_change_type(rb_d, NA_SCOMPLEX);
  d = NA_PTR_TYPE(rb_d, complex*);
  if (!NA_IsNArray(rb_du))
    rb_raise(rb_eArgError, "du (4th argument) must be NArray");
  if (NA_RANK(rb_du) != 1)
    rb_raise(rb_eArgError, "rank of du (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_du) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of du must be %d", n-1);
  if (NA_TYPE(rb_du) != NA_SCOMPLEX)
    rb_du = na_change_type(rb_du, NA_SCOMPLEX);
  du = NA_PTR_TYPE(rb_du, complex*);
  if (!NA_IsNArray(rb_dlf))
    rb_raise(rb_eArgError, "dlf (5th argument) must be NArray");
  if (NA_RANK(rb_dlf) != 1)
    rb_raise(rb_eArgError, "rank of dlf (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_dlf) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of dlf must be %d", n-1);
  if (NA_TYPE(rb_dlf) != NA_SCOMPLEX)
    rb_dlf = na_change_type(rb_dlf, NA_SCOMPLEX);
  dlf = NA_PTR_TYPE(rb_dlf, complex*);
  if (!NA_IsNArray(rb_dl))
    rb_raise(rb_eArgError, "dl (2th argument) must be NArray");
  if (NA_RANK(rb_dl) != 1)
    rb_raise(rb_eArgError, "rank of dl (2th argument) must be %d", 1);
  if (NA_SHAPE0(rb_dl) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of dl must be %d", n-1);
  if (NA_TYPE(rb_dl) != NA_SCOMPLEX)
    rb_dl = na_change_type(rb_dl, NA_SCOMPLEX);
  dl = NA_PTR_TYPE(rb_dl, complex*);
  if (!NA_IsNArray(rb_duf))
    rb_raise(rb_eArgError, "duf (7th argument) must be NArray");
  if (NA_RANK(rb_duf) != 1)
    rb_raise(rb_eArgError, "rank of duf (7th argument) must be %d", 1);
  if (NA_SHAPE0(rb_duf) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of duf must be %d", n-1);
  if (NA_TYPE(rb_duf) != NA_SCOMPLEX)
    rb_duf = na_change_type(rb_duf, NA_SCOMPLEX);
  duf = NA_PTR_TYPE(rb_duf, complex*);
  if (!NA_IsNArray(rb_du2))
    rb_raise(rb_eArgError, "du2 (8th argument) must be NArray");
  if (NA_RANK(rb_du2) != 1)
    rb_raise(rb_eArgError, "rank of du2 (8th argument) must be %d", 1);
  if (NA_SHAPE0(rb_du2) != (n-2))
    rb_raise(rb_eRuntimeError, "shape 0 of du2 must be %d", n-2);
  if (NA_TYPE(rb_du2) != NA_SCOMPLEX)
    rb_du2 = na_change_type(rb_du2, NA_SCOMPLEX);
  du2 = NA_PTR_TYPE(rb_du2, complex*);
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
    int shape[2];
    shape[0] = ldx;
    shape[1] = nrhs;
    rb_x_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  x_out__ = NA_PTR_TYPE(rb_x_out__, complex*);
  MEMCPY(x_out__, x, complex, NA_TOTAL(rb_x));
  rb_x = rb_x_out__;
  x = x_out__;
  work = ALLOC_N(complex, (2*n));
  rwork = ALLOC_N(real, (n));

  cgtrfs_(&trans, &n, &nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, &ldb, x, &ldx, ferr, berr, work, rwork, &info);

  free(work);
  free(rwork);
  rb_info = INT2NUM(info);
  return rb_ary_new3(4, rb_ferr, rb_berr, rb_info, rb_x);
}

void
init_lapack_cgtrfs(VALUE mLapack){
  rb_define_module_function(mLapack, "cgtrfs", rb_cgtrfs, -1);
}
