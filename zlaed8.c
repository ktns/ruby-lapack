#include "rb_lapack.h"

extern VOID zlaed8_(integer *k, integer *n, integer *qsiz, doublecomplex *q, integer *ldq, doublereal *d, doublereal *rho, integer *cutpnt, doublereal *z, doublereal *dlamda, doublecomplex *q2, integer *ldq2, doublereal *w, integer *indxp, integer *indx, integer *indxq, integer *perm, integer *givptr, integer *givcol, doublereal *givnum, integer *info);

static VALUE
rb_zlaed8(int argc, VALUE *argv, VALUE self){
  VALUE rb_qsiz;
  integer qsiz; 
  VALUE rb_q;
  doublecomplex *q; 
  VALUE rb_d;
  doublereal *d; 
  VALUE rb_rho;
  doublereal rho; 
  VALUE rb_cutpnt;
  integer cutpnt; 
  VALUE rb_z;
  doublereal *z; 
  VALUE rb_indxq;
  integer *indxq; 
  VALUE rb_k;
  integer k; 
  VALUE rb_dlamda;
  doublereal *dlamda; 
  VALUE rb_q2;
  doublecomplex *q2; 
  VALUE rb_w;
  doublereal *w; 
  VALUE rb_perm;
  integer *perm; 
  VALUE rb_givptr;
  integer givptr; 
  VALUE rb_givcol;
  integer *givcol; 
  VALUE rb_givnum;
  doublereal *givnum; 
  VALUE rb_info;
  integer info; 
  VALUE rb_q_out__;
  doublecomplex *q_out__;
  VALUE rb_d_out__;
  doublereal *d_out__;
  integer *indxp;
  integer *indx;

  integer ldq;
  integer n;
  integer ldq2;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  k, dlamda, q2, w, perm, givptr, givcol, givnum, info, q, d, rho = NumRu::Lapack.zlaed8( qsiz, q, d, rho, cutpnt, z, indxq)\n    or\n  NumRu::Lapack.zlaed8  # print help\n\n\nFORTRAN MANUAL\n\n");
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
  rho = NUM2DBL(rb_rho);
  if (!NA_IsNArray(rb_q))
    rb_raise(rb_eArgError, "q (2th argument) must be NArray");
  if (NA_RANK(rb_q) != 2)
    rb_raise(rb_eArgError, "rank of q (2th argument) must be %d", 2);
  if (NA_SHAPE1(rb_q) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of q must be the same as shape 0 of indxq");
  ldq = NA_SHAPE0(rb_q);
  if (NA_TYPE(rb_q) != NA_DCOMPLEX)
    rb_q = na_change_type(rb_q, NA_DCOMPLEX);
  q = NA_PTR_TYPE(rb_q, doublecomplex*);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (3th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_d) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of d must be the same as shape 0 of indxq");
  if (NA_TYPE(rb_d) != NA_DFLOAT)
    rb_d = na_change_type(rb_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rb_d, doublereal*);
  if (!NA_IsNArray(rb_z))
    rb_raise(rb_eArgError, "z (6th argument) must be NArray");
  if (NA_RANK(rb_z) != 1)
    rb_raise(rb_eArgError, "rank of z (6th argument) must be %d", 1);
  if (NA_SHAPE0(rb_z) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of z must be the same as shape 0 of indxq");
  if (NA_TYPE(rb_z) != NA_DFLOAT)
    rb_z = na_change_type(rb_z, NA_DFLOAT);
  z = NA_PTR_TYPE(rb_z, doublereal*);
  ldq2 = MAX( 1, n );
  {
    int shape[1];
    shape[0] = n;
    rb_dlamda = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  dlamda = NA_PTR_TYPE(rb_dlamda, doublereal*);
  {
    int shape[2];
    shape[0] = ldq2;
    shape[1] = n;
    rb_q2 = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  q2 = NA_PTR_TYPE(rb_q2, doublecomplex*);
  {
    int shape[1];
    shape[0] = n;
    rb_w = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  w = NA_PTR_TYPE(rb_w, doublereal*);
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
    rb_givnum = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  givnum = NA_PTR_TYPE(rb_givnum, doublereal*);
  {
    int shape[2];
    shape[0] = ldq;
    shape[1] = n;
    rb_q_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  q_out__ = NA_PTR_TYPE(rb_q_out__, doublecomplex*);
  MEMCPY(q_out__, q, doublecomplex, NA_TOTAL(rb_q));
  rb_q = rb_q_out__;
  q = q_out__;
  {
    int shape[1];
    shape[0] = n;
    rb_d_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  d_out__ = NA_PTR_TYPE(rb_d_out__, doublereal*);
  MEMCPY(d_out__, d, doublereal, NA_TOTAL(rb_d));
  rb_d = rb_d_out__;
  d = d_out__;
  indxp = ALLOC_N(integer, (n));
  indx = ALLOC_N(integer, (n));

  zlaed8_(&k, &n, &qsiz, q, &ldq, d, &rho, &cutpnt, z, dlamda, q2, &ldq2, w, indxp, indx, indxq, perm, &givptr, givcol, givnum, &info);

  free(indxp);
  free(indx);
  rb_k = INT2NUM(k);
  rb_givptr = INT2NUM(givptr);
  rb_info = INT2NUM(info);
  rb_rho = rb_float_new((double)rho);
  return rb_ary_new3(12, rb_k, rb_dlamda, rb_q2, rb_w, rb_perm, rb_givptr, rb_givcol, rb_givnum, rb_info, rb_q, rb_d, rb_rho);
}

void
init_lapack_zlaed8(VALUE mLapack){
  rb_define_module_function(mLapack, "zlaed8", rb_zlaed8, -1);
}
