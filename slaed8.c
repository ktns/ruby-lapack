#include "rb_lapack.h"

extern VOID slaed8_(integer *icompq, integer *k, integer *n, integer *qsiz, real *d, real *q, integer *ldq, integer *indxq, real *rho, integer *cutpnt, real *z, real *dlamda, real *q2, integer *ldq2, real *w, integer *perm, integer *givptr, integer *givcol, real *givnum, integer *indxp, integer *indx, integer *info);

static VALUE
rb_slaed8(int argc, VALUE *argv, VALUE self){
  VALUE rb_icompq;
  integer icompq; 
  VALUE rb_qsiz;
  integer qsiz; 
  VALUE rb_d;
  real *d; 
  VALUE rb_q;
  real *q; 
  VALUE rb_ldq;
  integer ldq; 
  VALUE rb_indxq;
  integer *indxq; 
  VALUE rb_rho;
  real rho; 
  VALUE rb_cutpnt;
  integer cutpnt; 
  VALUE rb_z;
  real *z; 
  VALUE rb_ldq2;
  integer ldq2; 
  VALUE rb_k;
  integer k; 
  VALUE rb_dlamda;
  real *dlamda; 
  VALUE rb_q2;
  real *q2; 
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
  VALUE rb_d_out__;
  real *d_out__;
  VALUE rb_q_out__;
  real *q_out__;
  integer *indxp;
  integer *indx;

  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  k, dlamda, q2, w, perm, givptr, givcol, givnum, info, d, q, rho = NumRu::Lapack.slaed8( icompq, qsiz, d, q, ldq, indxq, rho, cutpnt, z, ldq2)\n    or\n  NumRu::Lapack.slaed8  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 10)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 10)", argc);
  rb_icompq = argv[0];
  rb_qsiz = argv[1];
  rb_d = argv[2];
  rb_q = argv[3];
  rb_ldq = argv[4];
  rb_indxq = argv[5];
  rb_rho = argv[6];
  rb_cutpnt = argv[7];
  rb_z = argv[8];
  rb_ldq2 = argv[9];

  qsiz = NUM2INT(rb_qsiz);
  cutpnt = NUM2INT(rb_cutpnt);
  if (!NA_IsNArray(rb_indxq))
    rb_raise(rb_eArgError, "indxq (6th argument) must be NArray");
  if (NA_RANK(rb_indxq) != 1)
    rb_raise(rb_eArgError, "rank of indxq (6th argument) must be %d", 1);
  n = NA_SHAPE0(rb_indxq);
  if (NA_TYPE(rb_indxq) != NA_LINT)
    rb_indxq = na_change_type(rb_indxq, NA_LINT);
  indxq = NA_PTR_TYPE(rb_indxq, integer*);
  rho = (real)NUM2DBL(rb_rho);
  ldq = NUM2INT(rb_ldq);
  icompq = NUM2INT(rb_icompq);
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
    rb_raise(rb_eArgError, "z (9th argument) must be NArray");
  if (NA_RANK(rb_z) != 1)
    rb_raise(rb_eArgError, "rank of z (9th argument) must be %d", 1);
  if (NA_SHAPE0(rb_z) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of z must be the same as shape 0 of indxq");
  if (NA_TYPE(rb_z) != NA_SFLOAT)
    rb_z = na_change_type(rb_z, NA_SFLOAT);
  z = NA_PTR_TYPE(rb_z, real*);
  if (!NA_IsNArray(rb_q))
    rb_raise(rb_eArgError, "q (4th argument) must be NArray");
  if (NA_RANK(rb_q) != 2)
    rb_raise(rb_eArgError, "rank of q (4th argument) must be %d", 2);
  if (NA_SHAPE1(rb_q) != (icompq==0 ? 0 : n))
    rb_raise(rb_eRuntimeError, "shape 1 of q must be %d", icompq==0 ? 0 : n);
  if (NA_SHAPE0(rb_q) != (icompq==0 ? 0 : ldq))
    rb_raise(rb_eRuntimeError, "shape 0 of q must be %d", icompq==0 ? 0 : ldq);
  if (NA_TYPE(rb_q) != NA_SFLOAT)
    rb_q = na_change_type(rb_q, NA_SFLOAT);
  q = NA_PTR_TYPE(rb_q, real*);
  ldq2 = MAX(1,n);
  {
    int shape[1];
    shape[0] = n;
    rb_dlamda = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  dlamda = NA_PTR_TYPE(rb_dlamda, real*);
  {
    int shape[2];
    shape[0] = icompq==0 ? 0 : ldq2;
    shape[1] = icompq==0 ? 0 : n;
    rb_q2 = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  q2 = NA_PTR_TYPE(rb_q2, real*);
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
    int shape[1];
    shape[0] = n;
    rb_d_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  d_out__ = NA_PTR_TYPE(rb_d_out__, real*);
  MEMCPY(d_out__, d, real, NA_TOTAL(rb_d));
  rb_d = rb_d_out__;
  d = d_out__;
  {
    int shape[2];
    shape[0] = icompq==0 ? 0 : ldq;
    shape[1] = icompq==0 ? 0 : n;
    rb_q_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  q_out__ = NA_PTR_TYPE(rb_q_out__, real*);
  MEMCPY(q_out__, q, real, NA_TOTAL(rb_q));
  rb_q = rb_q_out__;
  q = q_out__;
  indxp = ALLOC_N(integer, (n));
  indx = ALLOC_N(integer, (n));

  slaed8_(&icompq, &k, &n, &qsiz, d, q, &ldq, indxq, &rho, &cutpnt, z, dlamda, q2, &ldq2, w, perm, &givptr, givcol, givnum, indxp, indx, &info);

  free(indxp);
  free(indx);
  rb_k = INT2NUM(k);
  rb_givptr = INT2NUM(givptr);
  rb_info = INT2NUM(info);
  rb_rho = rb_float_new((double)rho);
  return rb_ary_new3(12, rb_k, rb_dlamda, rb_q2, rb_w, rb_perm, rb_givptr, rb_givcol, rb_givnum, rb_info, rb_d, rb_q, rb_rho);
}

void
init_lapack_slaed8(VALUE mLapack){
  rb_define_module_function(mLapack, "slaed8", rb_slaed8, -1);
}
