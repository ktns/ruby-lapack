#include "rb_lapack.h"

extern VOID slaebz_(integer *ijob, integer *nitmax, integer *n, integer *mmax, integer *minp, integer *nbmin, real *abstol, real *reltol, real *pivmin, real *d, real *e, real *e2, integer *nval, real *ab, real *c, integer *mout, integer *nab, real *work, integer *iwork, integer *info);

static VALUE
rb_slaebz(int argc, VALUE *argv, VALUE self){
  VALUE rb_ijob;
  integer ijob; 
  VALUE rb_nitmax;
  integer nitmax; 
  VALUE rb_minp;
  integer minp; 
  VALUE rb_nbmin;
  integer nbmin; 
  VALUE rb_abstol;
  real abstol; 
  VALUE rb_reltol;
  real reltol; 
  VALUE rb_pivmin;
  real pivmin; 
  VALUE rb_d;
  real *d; 
  VALUE rb_e;
  real *e; 
  VALUE rb_e2;
  real *e2; 
  VALUE rb_nval;
  integer *nval; 
  VALUE rb_ab;
  real *ab; 
  VALUE rb_c;
  real *c; 
  VALUE rb_nab;
  integer *nab; 
  VALUE rb_mout;
  integer mout; 
  VALUE rb_info;
  integer info; 
  VALUE rb_nval_out__;
  integer *nval_out__;
  VALUE rb_ab_out__;
  real *ab_out__;
  VALUE rb_c_out__;
  real *c_out__;
  VALUE rb_nab_out__;
  integer *nab_out__;
  real *work;
  integer *iwork;

  integer n;
  integer mmax;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  mout, info, nval, ab, c, nab = NumRu::Lapack.slaebz( ijob, nitmax, minp, nbmin, abstol, reltol, pivmin, d, e, e2, nval, ab, c, nab)\n    or\n  NumRu::Lapack.slaebz  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 14)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 14)", argc);
  rb_ijob = argv[0];
  rb_nitmax = argv[1];
  rb_minp = argv[2];
  rb_nbmin = argv[3];
  rb_abstol = argv[4];
  rb_reltol = argv[5];
  rb_pivmin = argv[6];
  rb_d = argv[7];
  rb_e = argv[8];
  rb_e2 = argv[9];
  rb_nval = argv[10];
  rb_ab = argv[11];
  rb_c = argv[12];
  rb_nab = argv[13];

  abstol = (real)NUM2DBL(rb_abstol);
  ijob = NUM2INT(rb_ijob);
  if (!NA_IsNArray(rb_ab))
    rb_raise(rb_eArgError, "ab (12th argument) must be NArray");
  if (NA_RANK(rb_ab) != 2)
    rb_raise(rb_eArgError, "rank of ab (12th argument) must be %d", 2);
  if (NA_SHAPE1(rb_ab) != (2))
    rb_raise(rb_eRuntimeError, "shape 1 of ab must be %d", 2);
  mmax = NA_SHAPE0(rb_ab);
  if (NA_TYPE(rb_ab) != NA_SFLOAT)
    rb_ab = na_change_type(rb_ab, NA_SFLOAT);
  ab = NA_PTR_TYPE(rb_ab, real*);
  if (!NA_IsNArray(rb_e2))
    rb_raise(rb_eArgError, "e2 (10th argument) must be NArray");
  if (NA_RANK(rb_e2) != 1)
    rb_raise(rb_eArgError, "rank of e2 (10th argument) must be %d", 1);
  n = NA_SHAPE0(rb_e2);
  if (NA_TYPE(rb_e2) != NA_SFLOAT)
    rb_e2 = na_change_type(rb_e2, NA_SFLOAT);
  e2 = NA_PTR_TYPE(rb_e2, real*);
  nitmax = NUM2INT(rb_nitmax);
  pivmin = (real)NUM2DBL(rb_pivmin);
  if (!NA_IsNArray(rb_nab))
    rb_raise(rb_eArgError, "nab (14th argument) must be NArray");
  if (NA_RANK(rb_nab) != 2)
    rb_raise(rb_eArgError, "rank of nab (14th argument) must be %d", 2);
  if (NA_SHAPE1(rb_nab) != (2))
    rb_raise(rb_eRuntimeError, "shape 1 of nab must be %d", 2);
  if (NA_SHAPE0(rb_nab) != mmax)
    rb_raise(rb_eRuntimeError, "shape 0 of nab must be the same as shape 0 of ab");
  if (NA_TYPE(rb_nab) != NA_LINT)
    rb_nab = na_change_type(rb_nab, NA_LINT);
  nab = NA_PTR_TYPE(rb_nab, integer*);
  nbmin = NUM2INT(rb_nbmin);
  reltol = (real)NUM2DBL(rb_reltol);
  if (!NA_IsNArray(rb_e))
    rb_raise(rb_eArgError, "e (9th argument) must be NArray");
  if (NA_RANK(rb_e) != 1)
    rb_raise(rb_eArgError, "rank of e (9th argument) must be %d", 1);
  if (NA_SHAPE0(rb_e) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of e must be the same as shape 0 of e2");
  if (NA_TYPE(rb_e) != NA_SFLOAT)
    rb_e = na_change_type(rb_e, NA_SFLOAT);
  e = NA_PTR_TYPE(rb_e, real*);
  minp = NUM2INT(rb_minp);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (8th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (8th argument) must be %d", 1);
  if (NA_SHAPE0(rb_d) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of d must be the same as shape 0 of e2");
  if (NA_TYPE(rb_d) != NA_SFLOAT)
    rb_d = na_change_type(rb_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rb_d, real*);
  if (!NA_IsNArray(rb_nval))
    rb_raise(rb_eArgError, "nval (11th argument) must be NArray");
  if (NA_RANK(rb_nval) != 1)
    rb_raise(rb_eArgError, "rank of nval (11th argument) must be %d", 1);
  if (NA_SHAPE0(rb_nval) != ((ijob==1||ijob==2) ? 0 : ijob==3 ? minp : 0))
    rb_raise(rb_eRuntimeError, "shape 0 of nval must be %d", (ijob==1||ijob==2) ? 0 : ijob==3 ? minp : 0);
  if (NA_TYPE(rb_nval) != NA_LINT)
    rb_nval = na_change_type(rb_nval, NA_LINT);
  nval = NA_PTR_TYPE(rb_nval, integer*);
  if (!NA_IsNArray(rb_c))
    rb_raise(rb_eArgError, "c (13th argument) must be NArray");
  if (NA_RANK(rb_c) != 1)
    rb_raise(rb_eArgError, "rank of c (13th argument) must be %d", 1);
  if (NA_SHAPE0(rb_c) != (ijob==1 ? 0 : (ijob==2||ijob==3) ? mmax : 0))
    rb_raise(rb_eRuntimeError, "shape 0 of c must be %d", ijob==1 ? 0 : (ijob==2||ijob==3) ? mmax : 0);
  if (NA_TYPE(rb_c) != NA_SFLOAT)
    rb_c = na_change_type(rb_c, NA_SFLOAT);
  c = NA_PTR_TYPE(rb_c, real*);
  {
    int shape[1];
    shape[0] = (ijob==1||ijob==2) ? 0 : ijob==3 ? minp : 0;
    rb_nval_out__ = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  nval_out__ = NA_PTR_TYPE(rb_nval_out__, integer*);
  MEMCPY(nval_out__, nval, integer, NA_TOTAL(rb_nval));
  rb_nval = rb_nval_out__;
  nval = nval_out__;
  {
    int shape[2];
    shape[0] = mmax;
    shape[1] = 2;
    rb_ab_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  ab_out__ = NA_PTR_TYPE(rb_ab_out__, real*);
  MEMCPY(ab_out__, ab, real, NA_TOTAL(rb_ab));
  rb_ab = rb_ab_out__;
  ab = ab_out__;
  {
    int shape[1];
    shape[0] = ijob==1 ? 0 : (ijob==2||ijob==3) ? mmax : 0;
    rb_c_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  c_out__ = NA_PTR_TYPE(rb_c_out__, real*);
  MEMCPY(c_out__, c, real, NA_TOTAL(rb_c));
  rb_c = rb_c_out__;
  c = c_out__;
  {
    int shape[2];
    shape[0] = mmax;
    shape[1] = 2;
    rb_nab_out__ = na_make_object(NA_LINT, 2, shape, cNArray);
  }
  nab_out__ = NA_PTR_TYPE(rb_nab_out__, integer*);
  MEMCPY(nab_out__, nab, integer, NA_TOTAL(rb_nab));
  rb_nab = rb_nab_out__;
  nab = nab_out__;
  work = ALLOC_N(real, (mmax));
  iwork = ALLOC_N(integer, (mmax));

  slaebz_(&ijob, &nitmax, &n, &mmax, &minp, &nbmin, &abstol, &reltol, &pivmin, d, e, e2, nval, ab, c, &mout, nab, work, iwork, &info);

  free(work);
  free(iwork);
  rb_mout = INT2NUM(mout);
  rb_info = INT2NUM(info);
  return rb_ary_new3(6, rb_mout, rb_info, rb_nval, rb_ab, rb_c, rb_nab);
}

void
init_lapack_slaebz(VALUE mLapack){
  rb_define_module_function(mLapack, "slaebz", rb_slaebz, -1);
}
