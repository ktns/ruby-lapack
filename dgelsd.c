#include "rb_lapack.h"

extern VOID dgelsd_(integer *m, integer *n, integer *nrhs, doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *s, doublereal *rcond, integer *rank, doublereal *work, integer *lwork, integer *iwork, integer *info);

static VALUE
rb_dgelsd(int argc, VALUE *argv, VALUE self){
  VALUE rb_m;
  integer m; 
  VALUE rb_a;
  doublereal *a; 
  VALUE rb_b;
  doublereal *b; 
  VALUE rb_rcond;
  doublereal rcond; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_s;
  doublereal *s; 
  VALUE rb_rank;
  integer rank; 
  VALUE rb_work;
  doublereal *work; 
  VALUE rb_info;
  integer info; 
  VALUE rb_b_out__;
  doublereal *b_out__;
  integer *iwork;

  integer lda;
  integer n;
  integer ldb;
  integer nrhs;
  integer c__9;
  integer c__0;
  integer liwork;
  integer smlsiz;
  integer nlvl;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  s, rank, work, info, b = NumRu::Lapack.dgelsd( m, a, b, rcond, lwork)\n    or\n  NumRu::Lapack.dgelsd  # print help\n\n\nFORTRAN MANUAL\n\n");
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
  if (NA_TYPE(rb_a) != NA_DFLOAT)
    rb_a = na_change_type(rb_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rb_a, doublereal*);
  c__0 = 0;
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (3th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (3th argument) must be %d", 2);
  nrhs = NA_SHAPE1(rb_b);
  ldb = NA_SHAPE0(rb_b);
  if (NA_TYPE(rb_b) != NA_DFLOAT)
    rb_b = na_change_type(rb_b, NA_DFLOAT);
  b = NA_PTR_TYPE(rb_b, doublereal*);
  lwork = NUM2INT(rb_lwork);
  m = NUM2INT(rb_m);
  smlsiz = ilaenv_(&c__9,"DGELSD"," ",&c__0,&c__0,&c__0,&c__0);
  nlvl = MAX(0,((int)(log(((double)(MIN(m,n)))/(smlsiz+1))/log(2.0))+1));
  liwork = 3*(MIN(m,n))*nlvl+11*(MIN(m,n));
  {
    int shape[1];
    shape[0] = MIN(m,n);
    rb_s = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  s = NA_PTR_TYPE(rb_s, doublereal*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rb_work = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rb_work, doublereal*);
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
  iwork = ALLOC_N(integer, (MAX(1,liwork)));

  dgelsd_(&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, &rank, work, &lwork, iwork, &info);

  free(iwork);
  rb_rank = INT2NUM(rank);
  rb_info = INT2NUM(info);
  return rb_ary_new3(5, rb_s, rb_rank, rb_work, rb_info, rb_b);
}

void
init_lapack_dgelsd(VALUE mLapack){
  rb_define_module_function(mLapack, "dgelsd", rb_dgelsd, -1);
}
