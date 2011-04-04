#include "rb_lapack.h"

extern VOID dtgsen_(integer *ijob, logical *wantq, logical *wantz, logical *select, integer *n, doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *alphar, doublereal *alphai, doublereal *beta, doublereal *q, integer *ldq, doublereal *z, integer *ldz, integer *m, doublereal *pl, doublereal *pr, doublereal *dif, doublereal *work, integer *lwork, integer *iwork, integer *liwork, integer *info);

static VALUE
rb_dtgsen(int argc, VALUE *argv, VALUE self){
  VALUE rb_ijob;
  integer ijob; 
  VALUE rb_wantq;
  logical wantq; 
  VALUE rb_wantz;
  logical wantz; 
  VALUE rb_select;
  logical *select; 
  VALUE rb_a;
  doublereal *a; 
  VALUE rb_b;
  doublereal *b; 
  VALUE rb_q;
  doublereal *q; 
  VALUE rb_z;
  doublereal *z; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_liwork;
  integer liwork; 
  VALUE rb_alphar;
  doublereal *alphar; 
  VALUE rb_alphai;
  doublereal *alphai; 
  VALUE rb_beta;
  doublereal *beta; 
  VALUE rb_m;
  integer m; 
  VALUE rb_pl;
  doublereal pl; 
  VALUE rb_pr;
  doublereal pr; 
  VALUE rb_dif;
  doublereal *dif; 
  VALUE rb_work;
  doublereal *work; 
  VALUE rb_iwork;
  integer *iwork; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  doublereal *a_out__;
  VALUE rb_b_out__;
  doublereal *b_out__;
  VALUE rb_q_out__;
  doublereal *q_out__;
  VALUE rb_z_out__;
  doublereal *z_out__;

  integer n;
  integer lda;
  integer ldb;
  integer ldq;
  integer ldz;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  alphar, alphai, beta, m, pl, pr, dif, work, iwork, info, a, b, q, z = NumRu::Lapack.dtgsen( ijob, wantq, wantz, select, a, b, q, z, lwork, liwork)\n    or\n  NumRu::Lapack.dtgsen  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 10)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 10)", argc);
  rb_ijob = argv[0];
  rb_wantq = argv[1];
  rb_wantz = argv[2];
  rb_select = argv[3];
  rb_a = argv[4];
  rb_b = argv[5];
  rb_q = argv[6];
  rb_z = argv[7];
  rb_lwork = argv[8];
  rb_liwork = argv[9];

  ijob = NUM2INT(rb_ijob);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (5th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (5th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_DFLOAT)
    rb_a = na_change_type(rb_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rb_a, doublereal*);
  wantz = (rb_wantz == Qtrue);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (6th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (6th argument) must be %d", 2);
  if (NA_SHAPE1(rb_b) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of b must be the same as shape 1 of a");
  ldb = NA_SHAPE0(rb_b);
  if (NA_TYPE(rb_b) != NA_DFLOAT)
    rb_b = na_change_type(rb_b, NA_DFLOAT);
  b = NA_PTR_TYPE(rb_b, doublereal*);
  liwork = NUM2INT(rb_liwork);
  wantq = (rb_wantq == Qtrue);
  if (!NA_IsNArray(rb_z))
    rb_raise(rb_eArgError, "z (8th argument) must be NArray");
  if (NA_RANK(rb_z) != 2)
    rb_raise(rb_eArgError, "rank of z (8th argument) must be %d", 2);
  if (NA_SHAPE1(rb_z) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of z must be the same as shape 1 of a");
  ldz = NA_SHAPE0(rb_z);
  if (NA_TYPE(rb_z) != NA_DFLOAT)
    rb_z = na_change_type(rb_z, NA_DFLOAT);
  z = NA_PTR_TYPE(rb_z, doublereal*);
  lwork = NUM2INT(rb_lwork);
  if (!NA_IsNArray(rb_q))
    rb_raise(rb_eArgError, "q (7th argument) must be NArray");
  if (NA_RANK(rb_q) != 2)
    rb_raise(rb_eArgError, "rank of q (7th argument) must be %d", 2);
  if (NA_SHAPE1(rb_q) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of q must be the same as shape 1 of a");
  ldq = NA_SHAPE0(rb_q);
  if (NA_TYPE(rb_q) != NA_DFLOAT)
    rb_q = na_change_type(rb_q, NA_DFLOAT);
  q = NA_PTR_TYPE(rb_q, doublereal*);
  if (!NA_IsNArray(rb_select))
    rb_raise(rb_eArgError, "select (4th argument) must be NArray");
  if (NA_RANK(rb_select) != 1)
    rb_raise(rb_eArgError, "rank of select (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_select) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of select must be the same as shape 1 of a");
  if (NA_TYPE(rb_select) != NA_LINT)
    rb_select = na_change_type(rb_select, NA_LINT);
  select = NA_PTR_TYPE(rb_select, logical*);
  {
    int shape[1];
    shape[0] = n;
    rb_alphar = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  alphar = NA_PTR_TYPE(rb_alphar, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_alphai = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  alphai = NA_PTR_TYPE(rb_alphai, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_beta = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  beta = NA_PTR_TYPE(rb_beta, doublereal*);
  {
    int shape[1];
    shape[0] = 2;
    rb_dif = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  dif = NA_PTR_TYPE(rb_dif, doublereal*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rb_work = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rb_work, doublereal*);
  {
    int shape[1];
    shape[0] = ijob==0 ? 0 : MAX(1,liwork);
    rb_iwork = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  iwork = NA_PTR_TYPE(rb_iwork, integer*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rb_a_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rb_a_out__, doublereal*);
  MEMCPY(a_out__, a, doublereal, NA_TOTAL(rb_a));
  rb_a = rb_a_out__;
  a = a_out__;
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = n;
    rb_b_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rb_b_out__, doublereal*);
  MEMCPY(b_out__, b, doublereal, NA_TOTAL(rb_b));
  rb_b = rb_b_out__;
  b = b_out__;
  {
    int shape[2];
    shape[0] = ldq;
    shape[1] = n;
    rb_q_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  q_out__ = NA_PTR_TYPE(rb_q_out__, doublereal*);
  MEMCPY(q_out__, q, doublereal, NA_TOTAL(rb_q));
  rb_q = rb_q_out__;
  q = q_out__;
  {
    int shape[2];
    shape[0] = ldz;
    shape[1] = n;
    rb_z_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  z_out__ = NA_PTR_TYPE(rb_z_out__, doublereal*);
  MEMCPY(z_out__, z, doublereal, NA_TOTAL(rb_z));
  rb_z = rb_z_out__;
  z = z_out__;

  dtgsen_(&ijob, &wantq, &wantz, select, &n, a, &lda, b, &ldb, alphar, alphai, beta, q, &ldq, z, &ldz, &m, &pl, &pr, dif, work, &lwork, iwork, &liwork, &info);

  rb_m = INT2NUM(m);
  rb_pl = rb_float_new((double)pl);
  rb_pr = rb_float_new((double)pr);
  rb_info = INT2NUM(info);
  return rb_ary_new3(14, rb_alphar, rb_alphai, rb_beta, rb_m, rb_pl, rb_pr, rb_dif, rb_work, rb_iwork, rb_info, rb_a, rb_b, rb_q, rb_z);
}

void
init_lapack_dtgsen(VALUE mLapack){
  rb_define_module_function(mLapack, "dtgsen", rb_dtgsen, -1);
}
