#include "rb_lapack.h"

extern VOID cptrfs_(char *uplo, integer *n, integer *nrhs, real *d, complex *e, real *df, complex *ef, complex *b, integer *ldb, complex *x, integer *ldx, real *ferr, real *berr, complex *work, real *rwork, integer *info);

static VALUE
rb_cptrfs(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_d;
  real *d; 
  VALUE rb_e;
  complex *e; 
  VALUE rb_df;
  real *df; 
  VALUE rb_ef;
  complex *ef; 
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
    printf("%s\n", "USAGE:\n  ferr, berr, info, x = NumRu::Lapack.cptrfs( uplo, d, e, df, ef, b, x)\n    or\n  NumRu::Lapack.cptrfs  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 7)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 7)", argc);
  rb_uplo = argv[0];
  rb_d = argv[1];
  rb_e = argv[2];
  rb_df = argv[3];
  rb_ef = argv[4];
  rb_b = argv[5];
  rb_x = argv[6];

  uplo = StringValueCStr(rb_uplo)[0];
  if (!NA_IsNArray(rb_x))
    rb_raise(rb_eArgError, "x (7th argument) must be NArray");
  if (NA_RANK(rb_x) != 2)
    rb_raise(rb_eArgError, "rank of x (7th argument) must be %d", 2);
  nrhs = NA_SHAPE1(rb_x);
  ldx = NA_SHAPE0(rb_x);
  if (NA_TYPE(rb_x) != NA_SCOMPLEX)
    rb_x = na_change_type(rb_x, NA_SCOMPLEX);
  x = NA_PTR_TYPE(rb_x, complex*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (6th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (6th argument) must be %d", 2);
  if (NA_SHAPE1(rb_b) != nrhs)
    rb_raise(rb_eRuntimeError, "shape 1 of b must be the same as shape 1 of x");
  ldb = NA_SHAPE0(rb_b);
  if (NA_TYPE(rb_b) != NA_SCOMPLEX)
    rb_b = na_change_type(rb_b, NA_SCOMPLEX);
  b = NA_PTR_TYPE(rb_b, complex*);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (2th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (2th argument) must be %d", 1);
  n = NA_SHAPE0(rb_d);
  if (NA_TYPE(rb_d) != NA_SFLOAT)
    rb_d = na_change_type(rb_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rb_d, real*);
  if (!NA_IsNArray(rb_df))
    rb_raise(rb_eArgError, "df (4th argument) must be NArray");
  if (NA_RANK(rb_df) != 1)
    rb_raise(rb_eArgError, "rank of df (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_df) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of df must be the same as shape 0 of d");
  if (NA_TYPE(rb_df) != NA_SFLOAT)
    rb_df = na_change_type(rb_df, NA_SFLOAT);
  df = NA_PTR_TYPE(rb_df, real*);
  if (!NA_IsNArray(rb_e))
    rb_raise(rb_eArgError, "e (3th argument) must be NArray");
  if (NA_RANK(rb_e) != 1)
    rb_raise(rb_eArgError, "rank of e (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_e) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of e must be %d", n-1);
  if (NA_TYPE(rb_e) != NA_SCOMPLEX)
    rb_e = na_change_type(rb_e, NA_SCOMPLEX);
  e = NA_PTR_TYPE(rb_e, complex*);
  if (!NA_IsNArray(rb_ef))
    rb_raise(rb_eArgError, "ef (5th argument) must be NArray");
  if (NA_RANK(rb_ef) != 1)
    rb_raise(rb_eArgError, "rank of ef (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_ef) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of ef must be %d", n-1);
  if (NA_TYPE(rb_ef) != NA_SCOMPLEX)
    rb_ef = na_change_type(rb_ef, NA_SCOMPLEX);
  ef = NA_PTR_TYPE(rb_ef, complex*);
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
  work = ALLOC_N(complex, (n));
  rwork = ALLOC_N(real, (n));

  cptrfs_(&uplo, &n, &nrhs, d, e, df, ef, b, &ldb, x, &ldx, ferr, berr, work, rwork, &info);

  free(work);
  free(rwork);
  rb_info = INT2NUM(info);
  return rb_ary_new3(4, rb_ferr, rb_berr, rb_info, rb_x);
}

void
init_lapack_cptrfs(VALUE mLapack){
  rb_define_module_function(mLapack, "cptrfs", rb_cptrfs, -1);
}
