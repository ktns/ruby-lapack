#include "rb_lapack.h"

extern VOID cpbsvx_(char *fact, char *uplo, integer *n, integer *kd, integer *nrhs, complex *ab, integer *ldab, complex *afb, integer *ldafb, char *equed, real *s, complex *b, integer *ldb, complex *x, integer *ldx, real *rcond, real *ferr, real *berr, complex *work, real *rwork, integer *info);

static VALUE
rb_cpbsvx(int argc, VALUE *argv, VALUE self){
  VALUE rb_fact;
  char fact; 
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_kd;
  integer kd; 
  VALUE rb_ab;
  complex *ab; 
  VALUE rb_afb;
  complex *afb; 
  VALUE rb_equed;
  char equed; 
  VALUE rb_s;
  real *s; 
  VALUE rb_b;
  complex *b; 
  VALUE rb_x;
  complex *x; 
  VALUE rb_rcond;
  real rcond; 
  VALUE rb_ferr;
  real *ferr; 
  VALUE rb_berr;
  real *berr; 
  VALUE rb_info;
  integer info; 
  VALUE rb_ab_out__;
  complex *ab_out__;
  VALUE rb_afb_out__;
  complex *afb_out__;
  VALUE rb_s_out__;
  real *s_out__;
  VALUE rb_b_out__;
  complex *b_out__;
  complex *work;
  real *rwork;

  integer ldab;
  integer n;
  integer ldafb;
  integer ldb;
  integer nrhs;
  integer ldx;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  x, rcond, ferr, berr, info, ab, afb, equed, s, b = NumRu::Lapack.cpbsvx( fact, uplo, kd, ab, afb, equed, s, b)\n    or\n  NumRu::Lapack.cpbsvx  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 8)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 8)", argc);
  rb_fact = argv[0];
  rb_uplo = argv[1];
  rb_kd = argv[2];
  rb_ab = argv[3];
  rb_afb = argv[4];
  rb_equed = argv[5];
  rb_s = argv[6];
  rb_b = argv[7];

  uplo = StringValueCStr(rb_uplo)[0];
  if (!NA_IsNArray(rb_ab))
    rb_raise(rb_eArgError, "ab (4th argument) must be NArray");
  if (NA_RANK(rb_ab) != 2)
    rb_raise(rb_eArgError, "rank of ab (4th argument) must be %d", 2);
  n = NA_SHAPE1(rb_ab);
  ldab = NA_SHAPE0(rb_ab);
  if (NA_TYPE(rb_ab) != NA_SCOMPLEX)
    rb_ab = na_change_type(rb_ab, NA_SCOMPLEX);
  ab = NA_PTR_TYPE(rb_ab, complex*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (8th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (8th argument) must be %d", 2);
  nrhs = NA_SHAPE1(rb_b);
  ldb = NA_SHAPE0(rb_b);
  if (NA_TYPE(rb_b) != NA_SCOMPLEX)
    rb_b = na_change_type(rb_b, NA_SCOMPLEX);
  b = NA_PTR_TYPE(rb_b, complex*);
  kd = NUM2INT(rb_kd);
  equed = StringValueCStr(rb_equed)[0];
  if (!NA_IsNArray(rb_afb))
    rb_raise(rb_eArgError, "afb (5th argument) must be NArray");
  if (NA_RANK(rb_afb) != 2)
    rb_raise(rb_eArgError, "rank of afb (5th argument) must be %d", 2);
  if (NA_SHAPE1(rb_afb) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of afb must be the same as shape 1 of ab");
  ldafb = NA_SHAPE0(rb_afb);
  if (NA_TYPE(rb_afb) != NA_SCOMPLEX)
    rb_afb = na_change_type(rb_afb, NA_SCOMPLEX);
  afb = NA_PTR_TYPE(rb_afb, complex*);
  if (!NA_IsNArray(rb_s))
    rb_raise(rb_eArgError, "s (7th argument) must be NArray");
  if (NA_RANK(rb_s) != 1)
    rb_raise(rb_eArgError, "rank of s (7th argument) must be %d", 1);
  if (NA_SHAPE0(rb_s) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of s must be the same as shape 1 of ab");
  if (NA_TYPE(rb_s) != NA_SFLOAT)
    rb_s = na_change_type(rb_s, NA_SFLOAT);
  s = NA_PTR_TYPE(rb_s, real*);
  fact = StringValueCStr(rb_fact)[0];
  ldx = MAX(1,n);
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
    shape[0] = ldab;
    shape[1] = n;
    rb_ab_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  ab_out__ = NA_PTR_TYPE(rb_ab_out__, complex*);
  MEMCPY(ab_out__, ab, complex, NA_TOTAL(rb_ab));
  rb_ab = rb_ab_out__;
  ab = ab_out__;
  {
    int shape[2];
    shape[0] = ldafb;
    shape[1] = n;
    rb_afb_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  afb_out__ = NA_PTR_TYPE(rb_afb_out__, complex*);
  MEMCPY(afb_out__, afb, complex, NA_TOTAL(rb_afb));
  rb_afb = rb_afb_out__;
  afb = afb_out__;
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
    shape[0] = ldb;
    shape[1] = nrhs;
    rb_b_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rb_b_out__, complex*);
  MEMCPY(b_out__, b, complex, NA_TOTAL(rb_b));
  rb_b = rb_b_out__;
  b = b_out__;
  work = ALLOC_N(complex, (2*n));
  rwork = ALLOC_N(real, (n));

  cpbsvx_(&fact, &uplo, &n, &kd, &nrhs, ab, &ldab, afb, &ldafb, &equed, s, b, &ldb, x, &ldx, &rcond, ferr, berr, work, rwork, &info);

  free(work);
  free(rwork);
  rb_rcond = rb_float_new((double)rcond);
  rb_info = INT2NUM(info);
  rb_equed = rb_str_new(&equed,1);
  return rb_ary_new3(10, rb_x, rb_rcond, rb_ferr, rb_berr, rb_info, rb_ab, rb_afb, rb_equed, rb_s, rb_b);
}

void
init_lapack_cpbsvx(VALUE mLapack){
  rb_define_module_function(mLapack, "cpbsvx", rb_cpbsvx, -1);
}
