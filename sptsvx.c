#include "rb_lapack.h"

extern VOID sptsvx_(char *fact, integer *n, integer *nrhs, real *d, real *e, real *df, real *ef, real *b, integer *ldb, real *x, integer *ldx, real *rcond, real *ferr, real *berr, real *work, integer *info);

static VALUE
rb_sptsvx(int argc, VALUE *argv, VALUE self){
  VALUE rb_fact;
  char fact; 
  VALUE rb_d;
  real *d; 
  VALUE rb_e;
  real *e; 
  VALUE rb_df;
  real *df; 
  VALUE rb_ef;
  real *ef; 
  VALUE rb_b;
  real *b; 
  VALUE rb_x;
  real *x; 
  VALUE rb_rcond;
  real rcond; 
  VALUE rb_ferr;
  real *ferr; 
  VALUE rb_berr;
  real *berr; 
  VALUE rb_info;
  integer info; 
  VALUE rb_df_out__;
  real *df_out__;
  VALUE rb_ef_out__;
  real *ef_out__;
  real *work;

  integer n;
  integer ldb;
  integer nrhs;
  integer ldx;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  x, rcond, ferr, berr, info, df, ef = NumRu::Lapack.sptsvx( fact, d, e, df, ef, b)\n    or\n  NumRu::Lapack.sptsvx  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rb_fact = argv[0];
  rb_d = argv[1];
  rb_e = argv[2];
  rb_df = argv[3];
  rb_ef = argv[4];
  rb_b = argv[5];

  if (!NA_IsNArray(rb_df))
    rb_raise(rb_eArgError, "df (4th argument) must be NArray");
  if (NA_RANK(rb_df) != 1)
    rb_raise(rb_eArgError, "rank of df (4th argument) must be %d", 1);
  n = NA_SHAPE0(rb_df);
  if (NA_TYPE(rb_df) != NA_SFLOAT)
    rb_df = na_change_type(rb_df, NA_SFLOAT);
  df = NA_PTR_TYPE(rb_df, real*);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (2th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (2th argument) must be %d", 1);
  if (NA_SHAPE0(rb_d) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of d must be the same as shape 0 of df");
  if (NA_TYPE(rb_d) != NA_SFLOAT)
    rb_d = na_change_type(rb_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rb_d, real*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (6th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (6th argument) must be %d", 2);
  nrhs = NA_SHAPE1(rb_b);
  ldb = NA_SHAPE0(rb_b);
  if (NA_TYPE(rb_b) != NA_SFLOAT)
    rb_b = na_change_type(rb_b, NA_SFLOAT);
  b = NA_PTR_TYPE(rb_b, real*);
  fact = StringValueCStr(rb_fact)[0];
  ldx = MAX(1,n);
  if (!NA_IsNArray(rb_ef))
    rb_raise(rb_eArgError, "ef (5th argument) must be NArray");
  if (NA_RANK(rb_ef) != 1)
    rb_raise(rb_eArgError, "rank of ef (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_ef) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of ef must be %d", n-1);
  if (NA_TYPE(rb_ef) != NA_SFLOAT)
    rb_ef = na_change_type(rb_ef, NA_SFLOAT);
  ef = NA_PTR_TYPE(rb_ef, real*);
  if (!NA_IsNArray(rb_e))
    rb_raise(rb_eArgError, "e (3th argument) must be NArray");
  if (NA_RANK(rb_e) != 1)
    rb_raise(rb_eArgError, "rank of e (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_e) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of e must be %d", n-1);
  if (NA_TYPE(rb_e) != NA_SFLOAT)
    rb_e = na_change_type(rb_e, NA_SFLOAT);
  e = NA_PTR_TYPE(rb_e, real*);
  {
    int shape[2];
    shape[0] = ldx;
    shape[1] = nrhs;
    rb_x = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  x = NA_PTR_TYPE(rb_x, real*);
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
    shape[0] = n;
    rb_df_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  df_out__ = NA_PTR_TYPE(rb_df_out__, real*);
  MEMCPY(df_out__, df, real, NA_TOTAL(rb_df));
  rb_df = rb_df_out__;
  df = df_out__;
  {
    int shape[1];
    shape[0] = n-1;
    rb_ef_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  ef_out__ = NA_PTR_TYPE(rb_ef_out__, real*);
  MEMCPY(ef_out__, ef, real, NA_TOTAL(rb_ef));
  rb_ef = rb_ef_out__;
  ef = ef_out__;
  work = ALLOC_N(real, (2*n));

  sptsvx_(&fact, &n, &nrhs, d, e, df, ef, b, &ldb, x, &ldx, &rcond, ferr, berr, work, &info);

  free(work);
  rb_rcond = rb_float_new((double)rcond);
  rb_info = INT2NUM(info);
  return rb_ary_new3(7, rb_x, rb_rcond, rb_ferr, rb_berr, rb_info, rb_df, rb_ef);
}

void
init_lapack_sptsvx(VALUE mLapack){
  rb_define_module_function(mLapack, "sptsvx", rb_sptsvx, -1);
}
