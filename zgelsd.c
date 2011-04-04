#include "rb_lapack.h"

extern VOID zgelsd_(integer *m, integer *n, integer *nrhs, doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, doublereal *s, doublereal *rcond, integer *rank, doublecomplex *work, integer *lwork, doublereal *rwork, integer *iwork, integer *info);

static VALUE
rb_zgelsd(int argc, VALUE *argv, VALUE self){
  VALUE rb_m;
  integer m; 
  VALUE rb_a;
  doublecomplex *a; 
  VALUE rb_b;
  doublecomplex *b; 
  VALUE rb_rcond;
  doublereal rcond; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_s;
  doublereal *s; 
  VALUE rb_rank;
  integer rank; 
  VALUE rb_work;
  doublecomplex *work; 
  VALUE rb_info;
  integer info; 
  VALUE rb_b_out__;
  doublecomplex *b_out__;
  doublereal *rwork;
  integer *iwork;

  integer lda;
  integer n;
  integer ldb;
  integer nrhs;
  integer c__9;
  integer c__0;
  integer lrwork;
  integer liwork;
  integer smlsiz;
  integer nlvl;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  s, rank, work, info, b = NumRu::Lapack.zgelsd( m, a, b, rcond, lwork)\n    or\n  NumRu::Lapack.zgelsd  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_m = argv[0];
  rb_a = argv[1];
  rb_b = argv[2];
  rb_rcond = argv[3];
  rb_lwork = argv[4];

  c__9 = 9;
  rcond = NUM2DBL(rb_rcond);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_DCOMPLEX)
    rb_a = na_change_type(rb_a, NA_DCOMPLEX);
  a = NA_PTR_TYPE(rb_a, doublecomplex*);
  c__0 = 0;
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (3th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (3th argument) must be %d", 2);
  nrhs = NA_SHAPE1(rb_b);
  ldb = NA_SHAPE0(rb_b);
  if (NA_TYPE(rb_b) != NA_DCOMPLEX)
    rb_b = na_change_type(rb_b, NA_DCOMPLEX);
  b = NA_PTR_TYPE(rb_b, doublecomplex*);
  lwork = NUM2INT(rb_lwork);
  m = NUM2INT(rb_m);
  smlsiz = ilaenv_(&c__9,"ZGELSD"," ",&c__0,&c__0,&c__0,&c__0);
  nlvl = MAX(0,(int)(log(1.0*MIN(m,n)/(smlsiz+1))/log(2.0)));
  liwork = MAX(1,3*(MIN(m,n))*nlvl+11*(MIN(m,n)));
  lrwork = m>=n ? 10*n+2*n*smlsiz+8*n*nlvl+3*smlsiz*nrhs+(smlsiz+1)*(smlsiz+1) : 10*m+2*m*smlsiz+8*m*nlvl+2*smlsiz*nrhs+(smlsiz+1)*(smlsiz+1);
  {
    int shape[1];
    shape[0] = MIN(m,n);
    rb_s = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  s = NA_PTR_TYPE(rb_s, doublereal*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rb_work = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rb_work, doublecomplex*);
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = nrhs;
    rb_b_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rb_b_out__, doublecomplex*);
  MEMCPY(b_out__, b, doublecomplex, NA_TOTAL(rb_b));
  rb_b = rb_b_out__;
  b = b_out__;
  rwork = ALLOC_N(doublereal, (MAX(1,lrwork)));
  iwork = ALLOC_N(integer, (MAX(1,liwork)));

  zgelsd_(&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, &rank, work, &lwork, rwork, iwork, &info);

  free(rwork);
  free(iwork);
  rb_rank = INT2NUM(rank);
  rb_info = INT2NUM(info);
  return rb_ary_new3(5, rb_s, rb_rank, rb_work, rb_info, rb_b);
}

void
init_lapack_zgelsd(VALUE mLapack){
  rb_define_module_function(mLapack, "zgelsd", rb_zgelsd, -1);
}
