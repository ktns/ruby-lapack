#include "rb_lapack.h"

extern VOID dlalsd_(char *uplo, integer *smlsiz, integer *n, integer *nrhs, doublereal *d, doublereal *e, doublereal *b, integer *ldb, doublereal *rcond, integer *rank, doublereal *work, integer *iwork, integer *info);

static VALUE
rb_dlalsd(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_smlsiz;
  integer smlsiz; 
  VALUE rb_d;
  doublereal *d; 
  VALUE rb_e;
  doublereal *e; 
  VALUE rb_b;
  doublereal *b; 
  VALUE rb_rcond;
  doublereal rcond; 
  VALUE rb_rank;
  integer rank; 
  VALUE rb_info;
  integer info; 
  VALUE rb_d_out__;
  doublereal *d_out__;
  VALUE rb_e_out__;
  doublereal *e_out__;
  VALUE rb_b_out__;
  doublereal *b_out__;
  integer nlvl;
  doublereal *work;
  integer *iwork;

  integer n;
  integer ldb;
  integer nrhs;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  rank, info, d, e, b = NumRu::Lapack.dlalsd( uplo, smlsiz, d, e, b, rcond)\n    or\n  NumRu::Lapack.dlalsd  # print help\n\n\nFORTRAN MANUAL\n\n");
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

  rcond = NUM2DBL(rb_rcond);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (5th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (5th argument) must be %d", 2);
  nrhs = NA_SHAPE1(rb_b);
  ldb = NA_SHAPE0(rb_b);
  if (NA_TYPE(rb_b) != NA_DFLOAT)
    rb_b = na_change_type(rb_b, NA_DFLOAT);
  b = NA_PTR_TYPE(rb_b, doublereal*);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (3th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (3th argument) must be %d", 1);
  n = NA_SHAPE0(rb_d);
  if (NA_TYPE(rb_d) != NA_DFLOAT)
    rb_d = na_change_type(rb_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rb_d, doublereal*);
  smlsiz = NUM2INT(rb_smlsiz);
  uplo = StringValueCStr(rb_uplo)[0];
  nlvl = MAX(0, (int)(1.0/log(2.0)*log((double)n/(smlsiz+1))) + 1);
  if (!NA_IsNArray(rb_e))
    rb_raise(rb_eArgError, "e (4th argument) must be NArray");
  if (NA_RANK(rb_e) != 1)
    rb_raise(rb_eArgError, "rank of e (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_e) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of e must be %d", n-1);
  if (NA_TYPE(rb_e) != NA_DFLOAT)
    rb_e = na_change_type(rb_e, NA_DFLOAT);
  e = NA_PTR_TYPE(rb_e, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_d_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  d_out__ = NA_PTR_TYPE(rb_d_out__, doublereal*);
  MEMCPY(d_out__, d, doublereal, NA_TOTAL(rb_d));
  rb_d = rb_d_out__;
  d = d_out__;
  {
    int shape[1];
    shape[0] = n-1;
    rb_e_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  e_out__ = NA_PTR_TYPE(rb_e_out__, doublereal*);
  MEMCPY(e_out__, e, doublereal, NA_TOTAL(rb_e));
  rb_e = rb_e_out__;
  e = e_out__;
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = nrhs;
    rb_b_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rb_b_out__, doublereal*);
  MEMCPY(b_out__, b, doublereal, NA_TOTAL(rb_b));
  rb_b = rb_b_out__;
  b = b_out__;
  work = ALLOC_N(doublereal, (9*n + 2*n*smlsiz + 8*n*nlvl + n*nrhs + pow(smlsiz+1,2)));
  iwork = ALLOC_N(integer, (3*n*nlvl + 11*n));

  dlalsd_(&uplo, &smlsiz, &n, &nrhs, d, e, b, &ldb, &rcond, &rank, work, iwork, &info);

  free(work);
  free(iwork);
  rb_rank = INT2NUM(rank);
  rb_info = INT2NUM(info);
  return rb_ary_new3(5, rb_rank, rb_info, rb_d, rb_e, rb_b);
}

void
init_lapack_dlalsd(VALUE mLapack){
  rb_define_module_function(mLapack, "dlalsd", rb_dlalsd, -1);
}
