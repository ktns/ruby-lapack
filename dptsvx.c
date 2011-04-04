#include "rb_lapack.h"

extern VOID dptsvx_(char *fact, integer *n, integer *nrhs, doublereal *d, doublereal *e, doublereal *df, doublereal *ef, doublereal *b, integer *ldb, doublereal *x, integer *ldx, doublereal *rcond, doublereal *ferr, doublereal *berr, doublereal *work, integer *info);

static VALUE
rb_dptsvx(int argc, VALUE *argv, VALUE self){
  VALUE rb_fact;
  char fact; 
  VALUE rb_d;
  doublereal *d; 
  VALUE rb_e;
  doublereal *e; 
  VALUE rb_df;
  doublereal *df; 
  VALUE rb_ef;
  doublereal *ef; 
  VALUE rb_b;
  doublereal *b; 
  VALUE rb_x;
  doublereal *x; 
  VALUE rb_rcond;
  doublereal rcond; 
  VALUE rb_ferr;
  doublereal *ferr; 
  VALUE rb_berr;
  doublereal *berr; 
  VALUE rb_info;
  integer info; 
  VALUE rb_df_out__;
  doublereal *df_out__;
  VALUE rb_ef_out__;
  doublereal *ef_out__;
  doublereal *work;

  integer n;
  integer ldb;
  integer nrhs;
  integer ldx;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  x, rcond, ferr, berr, info, df, ef = NumRu::Lapack.dptsvx( fact, d, e, df, ef, b)\n    or\n  NumRu::Lapack.dptsvx  # print help\n\n\nFORTRAN MANUAL\n\n");
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
  if (NA_TYPE(rb_df) != NA_DFLOAT)
    rb_df = na_change_type(rb_df, NA_DFLOAT);
  df = NA_PTR_TYPE(rb_df, doublereal*);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (2th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (2th argument) must be %d", 1);
  if (NA_SHAPE0(rb_d) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of d must be the same as shape 0 of df");
  if (NA_TYPE(rb_d) != NA_DFLOAT)
    rb_d = na_change_type(rb_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rb_d, doublereal*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (6th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (6th argument) must be %d", 2);
  nrhs = NA_SHAPE1(rb_b);
  ldb = NA_SHAPE0(rb_b);
  if (NA_TYPE(rb_b) != NA_DFLOAT)
    rb_b = na_change_type(rb_b, NA_DFLOAT);
  b = NA_PTR_TYPE(rb_b, doublereal*);
  fact = StringValueCStr(rb_fact)[0];
  ldx = MAX(1,n);
  if (!NA_IsNArray(rb_ef))
    rb_raise(rb_eArgError, "ef (5th argument) must be NArray");
  if (NA_RANK(rb_ef) != 1)
    rb_raise(rb_eArgError, "rank of ef (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_ef) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of ef must be %d", n-1);
  if (NA_TYPE(rb_ef) != NA_DFLOAT)
    rb_ef = na_change_type(rb_ef, NA_DFLOAT);
  ef = NA_PTR_TYPE(rb_ef, doublereal*);
  if (!NA_IsNArray(rb_e))
    rb_raise(rb_eArgError, "e (3th argument) must be NArray");
  if (NA_RANK(rb_e) != 1)
    rb_raise(rb_eArgError, "rank of e (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_e) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of e must be %d", n-1);
  if (NA_TYPE(rb_e) != NA_DFLOAT)
    rb_e = na_change_type(rb_e, NA_DFLOAT);
  e = NA_PTR_TYPE(rb_e, doublereal*);
  {
    int shape[2];
    shape[0] = ldx;
    shape[1] = nrhs;
    rb_x = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  x = NA_PTR_TYPE(rb_x, doublereal*);
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
    shape[0] = n;
    rb_df_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  df_out__ = NA_PTR_TYPE(rb_df_out__, doublereal*);
  MEMCPY(df_out__, df, doublereal, NA_TOTAL(rb_df));
  rb_df = rb_df_out__;
  df = df_out__;
  {
    int shape[1];
    shape[0] = n-1;
    rb_ef_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  ef_out__ = NA_PTR_TYPE(rb_ef_out__, doublereal*);
  MEMCPY(ef_out__, ef, doublereal, NA_TOTAL(rb_ef));
  rb_ef = rb_ef_out__;
  ef = ef_out__;
  work = ALLOC_N(doublereal, (2*n));

  dptsvx_(&fact, &n, &nrhs, d, e, df, ef, b, &ldb, x, &ldx, &rcond, ferr, berr, work, &info);

  free(work);
  rb_rcond = rb_float_new((double)rcond);
  rb_info = INT2NUM(info);
  return rb_ary_new3(7, rb_x, rb_rcond, rb_ferr, rb_berr, rb_info, rb_df, rb_ef);
}

void
init_lapack_dptsvx(VALUE mLapack){
  rb_define_module_function(mLapack, "dptsvx", rb_dptsvx, -1);
}
