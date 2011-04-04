#include "rb_lapack.h"

extern VOID cgbsvx_(char *fact, char *trans, integer *n, integer *kl, integer *ku, integer *nrhs, complex *ab, integer *ldab, complex *afb, integer *ldafb, integer *ipiv, char *equed, real *r, real *c, complex *b, integer *ldb, complex *x, integer *ldx, real *rcond, real *ferr, real *berr, complex *work, real *rwork, integer *info);

static VALUE
rb_cgbsvx(int argc, VALUE *argv, VALUE self){
  VALUE rb_fact;
  char fact; 
  VALUE rb_trans;
  char trans; 
  VALUE rb_kl;
  integer kl; 
  VALUE rb_ku;
  integer ku; 
  VALUE rb_ab;
  complex *ab; 
  VALUE rb_afb;
  complex *afb; 
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
  VALUE rb_x;
  complex *x; 
  VALUE rb_rcond;
  real rcond; 
  VALUE rb_ferr;
  real *ferr; 
  VALUE rb_berr;
  real *berr; 
  VALUE rb_rwork;
  real *rwork; 
  VALUE rb_info;
  integer info; 
  VALUE rb_ab_out__;
  complex *ab_out__;
  VALUE rb_afb_out__;
  complex *afb_out__;
  VALUE rb_ipiv_out__;
  integer *ipiv_out__;
  VALUE rb_r_out__;
  real *r_out__;
  VALUE rb_c_out__;
  real *c_out__;
  VALUE rb_b_out__;
  complex *b_out__;
  complex *work;

  integer ldab;
  integer n;
  integer ldafb;
  integer ldb;
  integer nrhs;
  integer ldx;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  x, rcond, ferr, berr, rwork, info, ab, afb, ipiv, equed, r, c, b = NumRu::Lapack.cgbsvx( fact, trans, kl, ku, ab, afb, ipiv, equed, r, c, b)\n    or\n  NumRu::Lapack.cgbsvx  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 11)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 11)", argc);
  rb_fact = argv[0];
  rb_trans = argv[1];
  rb_kl = argv[2];
  rb_ku = argv[3];
  rb_ab = argv[4];
  rb_afb = argv[5];
  rb_ipiv = argv[6];
  rb_equed = argv[7];
  rb_r = argv[8];
  rb_c = argv[9];
  rb_b = argv[10];

  trans = StringValueCStr(rb_trans)[0];
  if (!NA_IsNArray(rb_ipiv))
    rb_raise(rb_eArgError, "ipiv (7th argument) must be NArray");
  if (NA_RANK(rb_ipiv) != 1)
    rb_raise(rb_eArgError, "rank of ipiv (7th argument) must be %d", 1);
  n = NA_SHAPE0(rb_ipiv);
  if (NA_TYPE(rb_ipiv) != NA_LINT)
    rb_ipiv = na_change_type(rb_ipiv, NA_LINT);
  ipiv = NA_PTR_TYPE(rb_ipiv, integer*);
  if (!NA_IsNArray(rb_ab))
    rb_raise(rb_eArgError, "ab (5th argument) must be NArray");
  if (NA_RANK(rb_ab) != 2)
    rb_raise(rb_eArgError, "rank of ab (5th argument) must be %d", 2);
  if (NA_SHAPE1(rb_ab) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of ab must be the same as shape 0 of ipiv");
  ldab = NA_SHAPE0(rb_ab);
  if (NA_TYPE(rb_ab) != NA_SCOMPLEX)
    rb_ab = na_change_type(rb_ab, NA_SCOMPLEX);
  ab = NA_PTR_TYPE(rb_ab, complex*);
  kl = NUM2INT(rb_kl);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (11th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (11th argument) must be %d", 2);
  nrhs = NA_SHAPE1(rb_b);
  ldb = NA_SHAPE0(rb_b);
  if (NA_TYPE(rb_b) != NA_SCOMPLEX)
    rb_b = na_change_type(rb_b, NA_SCOMPLEX);
  b = NA_PTR_TYPE(rb_b, complex*);
  if (!NA_IsNArray(rb_c))
    rb_raise(rb_eArgError, "c (10th argument) must be NArray");
  if (NA_RANK(rb_c) != 1)
    rb_raise(rb_eArgError, "rank of c (10th argument) must be %d", 1);
  if (NA_SHAPE0(rb_c) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of c must be the same as shape 0 of ipiv");
  if (NA_TYPE(rb_c) != NA_SFLOAT)
    rb_c = na_change_type(rb_c, NA_SFLOAT);
  c = NA_PTR_TYPE(rb_c, real*);
  ku = NUM2INT(rb_ku);
  if (!NA_IsNArray(rb_r))
    rb_raise(rb_eArgError, "r (9th argument) must be NArray");
  if (NA_RANK(rb_r) != 1)
    rb_raise(rb_eArgError, "rank of r (9th argument) must be %d", 1);
  if (NA_SHAPE0(rb_r) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of r must be the same as shape 0 of ipiv");
  if (NA_TYPE(rb_r) != NA_SFLOAT)
    rb_r = na_change_type(rb_r, NA_SFLOAT);
  r = NA_PTR_TYPE(rb_r, real*);
  if (!NA_IsNArray(rb_afb))
    rb_raise(rb_eArgError, "afb (6th argument) must be NArray");
  if (NA_RANK(rb_afb) != 2)
    rb_raise(rb_eArgError, "rank of afb (6th argument) must be %d", 2);
  if (NA_SHAPE1(rb_afb) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of afb must be the same as shape 0 of ipiv");
  ldafb = NA_SHAPE0(rb_afb);
  if (NA_TYPE(rb_afb) != NA_SCOMPLEX)
    rb_afb = na_change_type(rb_afb, NA_SCOMPLEX);
  afb = NA_PTR_TYPE(rb_afb, complex*);
  equed = StringValueCStr(rb_equed)[0];
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
    int shape[1];
    shape[0] = n;
    rb_rwork = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  rwork = NA_PTR_TYPE(rb_rwork, real*);
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
  work = ALLOC_N(complex, (2*n));

  cgbsvx_(&fact, &trans, &n, &kl, &ku, &nrhs, ab, &ldab, afb, &ldafb, ipiv, &equed, r, c, b, &ldb, x, &ldx, &rcond, ferr, berr, work, rwork, &info);

  free(work);
  rb_rcond = rb_float_new((double)rcond);
  rb_info = INT2NUM(info);
  rb_equed = rb_str_new(&equed,1);
  return rb_ary_new3(13, rb_x, rb_rcond, rb_ferr, rb_berr, rb_rwork, rb_info, rb_ab, rb_afb, rb_ipiv, rb_equed, rb_r, rb_c, rb_b);
}

void
init_lapack_cgbsvx(VALUE mLapack){
  rb_define_module_function(mLapack, "cgbsvx", rb_cgbsvx, -1);
}
