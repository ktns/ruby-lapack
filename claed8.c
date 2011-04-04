#include "rb_lapack.h"

extern VOID claed8_(integer *k, integer *n, integer *qsiz, complex *q, integer *ldq, real *d, real *rho, integer *cutpnt, real *z, real *dlamda, complex *q2, integer *ldq2, real *w, integer *indxp, integer *indx, integer *indxq, integer *perm, integer *givptr, integer *givcol, real *givnum, integer *info);

static VALUE
rb_claed8(int argc, VALUE *argv, VALUE self){
  VALUE rb_qsiz;
  integer qsiz; 
  VALUE rb_q;
  complex *q; 
  VALUE rb_d;
  real *d; 
  VALUE rb_rho;
  real rho; 
  VALUE rb_cutpnt;
  integer cutpnt; 
  VALUE rb_z;
  real *z; 
  VALUE rb_indxq;
  integer *indxq; 
  VALUE rb_k;
  integer k; 
  VALUE rb_dlamda;
  real *dlamda; 
  VALUE rb_q2;
  complex *q2; 
  VALUE rb_w;
  real *w; 
  VALUE rb_perm;
  integer *perm; 
  VALUE rb_givptr;
  integer givptr; 
  VALUE rb_givcol;
  integer *givcol; 
  VALUE rb_givnum;
  real *givnum; 
  VALUE rb_info;
  integer info; 
  VALUE rb_q_out__;
  complex *q_out__;
  VALUE rb_d_out__;
  real *d_out__;
  integer *indxp;
  integer *indx;

  integer ldq;
  integer n;
  integer ldq2;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  k, dlamda, q2, w, perm, givptr, givcol, givnum, info, q, d, rho = NumRu::Lapack.claed8( qsiz, q, d, rho, cutpnt, z, indxq)\n    or\n  NumRu::Lapack.claed8  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 7)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 7)", argc);
  rb_qsiz = argv[0];
  rb_q = argv[1];
  rb_d = argv[2];
  rb_rho = argv[3];
  rb_cutpnt = argv[4];
  rb_z = argv[5];
  rb_indxq = argv[6];

  qsiz = NUM2INT(rb_qsiz);
  if (!NA_IsNArray(rb_indxq))
    rb_raise(rb_eArgError, "indxq (7th argument) must be NArray");
  if (NA_RANK(rb_indxq) != 1)
    rb_raise(rb_eArgError, "rank of indxq (7th argument) must be %d", 1);
  n = NA_SHAPE0(rb_indxq);
  if (NA_TYPE(rb_indxq) != NA_LINT)
    rb_indxq = na_change_type(rb_indxq, NA_LINT);
  indxq = NA_PTR_TYPE(rb_indxq, integer*);
  cutpnt = NUM2INT(rb_cutpnt);
  rho = (real)NUM2DBL(rb_rho);
  if (!NA_IsNArray(rb_q))
    rb_raise(rb_eArgError, "q (2th argument) must be NArray");
  if (NA_RANK(rb_q) != 2)
    rb_raise(rb_eArgError, "rank of q (2th argument) must be %d", 2);
  if (NA_SHAPE1(rb_q) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of q must be the same as shape 0 of indxq");
  ldq = NA_SHAPE0(rb_q);
  if (NA_TYPE(rb_q) != NA_SCOMPLEX)
    rb_q = na_change_type(rb_q, NA_SCOMPLEX);
  q = NA_PTR_TYPE(rb_q, complex*);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (3th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_d) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of d must be the same as shape 0 of indxq");
  if (NA_TYPE(rb_d) != NA_SFLOAT)
    rb_d = na_change_type(rb_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rb_d, real*);
  if (!NA_IsNArray(rb_z))
    rb_raise(rb_eArgError, "z (6th argument) must be NArray");
  if (NA_RANK(rb_z) != 1)
    rb_raise(rb_eArgError, "rank of z (6th argument) must be %d", 1);
  if (NA_SHAPE0(rb_z) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of z must be the same as shape 0 of indxq");
  if (NA_TYPE(rb_z) != NA_SFLOAT)
    rb_z = na_change_type(rb_z, NA_SFLOAT);
  z = NA_PTR_TYPE(rb_z, real*);
  ldq2 = MAX( 1, n );
  {
    int shape[1];
    shape[0] = n;
    rb_dlamda = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  dlamda = NA_PTR_TYPE(rb_dlamda, real*);
  {
    int shape[2];
    shape[0] = ldq2;
    shape[1] = n;
    rb_q2 = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  q2 = NA_PTR_TYPE(rb_q2, complex*);
  {
    int shape[1];
    shape[0] = n;
    rb_w = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  w = NA_PTR_TYPE(rb_w, real*);
  {
    int shape[1];
    shape[0] = n;
    rb_perm = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  perm = NA_PTR_TYPE(rb_perm, integer*);
  {
    int shape[2];
    shape[0] = 2;
    shape[1] = n;
    rb_givcol = na_make_object(NA_LINT, 2, shape, cNArray);
  }
  givcol = NA_PTR_TYPE(rb_givcol, integer*);
  {
    int shape[2];
    shape[0] = 2;
    shape[1] = n;
    rb_givnum = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  givnum = NA_PTR_TYPE(rb_givnum, real*);
  {
    int shape[2];
    shape[0] = ldq;
    shape[1] = n;
    rb_q_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  q_out__ = NA_PTR_TYPE(rb_q_out__, complex*);
  MEMCPY(q_out__, q, complex, NA_TOTAL(rb_q));
  rb_q = rb_q_out__;
  q = q_out__;
  {
    int shape[1];
    shape[0] = n;
    rb_d_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  d_out__ = NA_PTR_TYPE(rb_d_out__, real*);
  MEMCPY(d_out__, d, real, NA_TOTAL(rb_d));
  rb_d = rb_d_out__;
  d = d_out__;
  indxp = ALLOC_N(integer, (n));
  indx = ALLOC_N(integer, (n));

  claed8_(&k, &n, &qsiz, q, &ldq, d, &rho, &cutpnt, z, dlamda, q2, &ldq2, w, indxp, indx, indxq, perm, &givptr, givcol, givnum, &info);

  free(indxp);
  free(indx);
  rb_k = INT2NUM(k);
  rb_givptr = INT2NUM(givptr);
  rb_info = INT2NUM(info);
  rb_rho = rb_float_new((double)rho);
  return rb_ary_new3(12, rb_k, rb_dlamda, rb_q2, rb_w, rb_perm, rb_givptr, rb_givcol, rb_givnum, rb_info, rb_q, rb_d, rb_rho);
}

void
init_lapack_claed8(VALUE mLapack){
  rb_define_module_function(mLapack, "claed8", rb_claed8, -1);
}
