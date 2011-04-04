#include "rb_lapack.h"

extern VOID clalsd_(char *uplo, integer *smlsiz, integer *n, integer *nrhs, real *d, real *e, complex *b, integer *ldb, real *rcond, integer *rank, complex *work, real *rwork, integer *iwork, integer *info);

static VALUE
rb_clalsd(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_smlsiz;
  integer smlsiz; 
  VALUE rb_d;
  real *d; 
  VALUE rb_e;
  real *e; 
  VALUE rb_b;
  complex *b; 
  VALUE rb_rcond;
  real rcond; 
  VALUE rb_rank;
  integer rank; 
  VALUE rb_info;
  integer info; 
  VALUE rb_d_out__;
  real *d_out__;
  VALUE rb_e_out__;
  real *e_out__;
  VALUE rb_b_out__;
  complex *b_out__;
  integer nlvl;
  complex *work;
  real *rwork;
  integer *iwork;

  integer n;
  integer ldb;
  integer nrhs;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  rank, info, d, e, b = NumRu::Lapack.clalsd( uplo, smlsiz, d, e, b, rcond)\n    or\n  NumRu::Lapack.clalsd  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rb_uplo = argv[0];
  rb_smlsiz = argv[1];
  rb_d = argv[2];
  rb_e = argv[3];
  rb_b = argv[4];
  rb_rcond = argv[5];

  rcond = (real)NUM2DBL(rb_rcond);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (5th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (5th argument) must be %d", 2);
  nrhs = NA_SHAPE1(rb_b);
  ldb = NA_SHAPE0(rb_b);
  if (NA_TYPE(rb_b) != NA_SCOMPLEX)
    rb_b = na_change_type(rb_b, NA_SCOMPLEX);
  b = NA_PTR_TYPE(rb_b, complex*);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (3th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (3th argument) must be %d", 1);
  n = NA_SHAPE0(rb_d);
  if (NA_TYPE(rb_d) != NA_SFLOAT)
    rb_d = na_change_type(rb_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rb_d, real*);
  smlsiz = NUM2INT(rb_smlsiz);
  uplo = StringValueCStr(rb_uplo)[0];
  nlvl = ( (int)( log(((double)n)/(smlsiz+1))/log(2.0) ) ) + 1;
  if (!NA_IsNArray(rb_e))
    rb_raise(rb_eArgError, "e (4th argument) must be NArray");
  if (NA_RANK(rb_e) != 1)
    rb_raise(rb_eArgError, "rank of e (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_e) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of e must be %d", n-1);
  if (NA_TYPE(rb_e) != NA_SFLOAT)
    rb_e = na_change_type(rb_e, NA_SFLOAT);
  e = NA_PTR_TYPE(rb_e, real*);
  {
    int shape[1];
    shape[0] = n;
    rb_d_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  d_out__ = NA_PTR_TYPE(rb_d_out__, real*);
  MEMCPY(d_out__, d, real, NA_TOTAL(rb_d));
  rb_d = rb_d_out__;
  d = d_out__;
  {
    int shape[1];
    shape[0] = n-1;
    rb_e_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  e_out__ = NA_PTR_TYPE(rb_e_out__, real*);
  MEMCPY(e_out__, e, real, NA_TOTAL(rb_e));
  rb_e = rb_e_out__;
  e = e_out__;
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
  work = ALLOC_N(complex, (n * nrhs));
  rwork = ALLOC_N(real, (9*n+2*n*smlsiz+8*n*nlvl+3*smlsiz*nrhs+(smlsiz+1)*(smlsiz+1)));
  iwork = ALLOC_N(integer, (3*n*nlvl + 11*n));

  clalsd_(&uplo, &smlsiz, &n, &nrhs, d, e, b, &ldb, &rcond, &rank, work, rwork, iwork, &info);

  free(work);
  free(rwork);
  free(iwork);
  rb_rank = INT2NUM(rank);
  rb_info = INT2NUM(info);
  return rb_ary_new3(5, rb_rank, rb_info, rb_d, rb_e, rb_b);
}

void
init_lapack_clalsd(VALUE mLapack){
  rb_define_module_function(mLapack, "clalsd", rb_clalsd, -1);
}
