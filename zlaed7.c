#include "rb_lapack.h"

extern VOID zlaed7_(integer *n, integer *cutpnt, integer *qsiz, integer *tlvls, integer *curlvl, integer *curpbm, doublereal *d, doublecomplex *q, integer *ldq, doublereal *rho, integer *indxq, doublereal *qstore, integer *qptr, integer *prmptr, integer *perm, integer *givptr, integer *givcol, doublereal *givnum, doublecomplex *work, doublereal *rwork, integer *iwork, integer *info);

static VALUE
rb_zlaed7(int argc, VALUE *argv, VALUE self){
  VALUE rb_cutpnt;
  integer cutpnt; 
  VALUE rb_qsiz;
  integer qsiz; 
  VALUE rb_tlvls;
  integer tlvls; 
  VALUE rb_curlvl;
  integer curlvl; 
  VALUE rb_curpbm;
  integer curpbm; 
  VALUE rb_d;
  doublereal *d; 
  VALUE rb_q;
  doublecomplex *q; 
  VALUE rb_rho;
  doublereal rho; 
  VALUE rb_qstore;
  doublereal *qstore; 
  VALUE rb_qptr;
  integer *qptr; 
  VALUE rb_prmptr;
  integer *prmptr; 
  VALUE rb_perm;
  integer *perm; 
  VALUE rb_givptr;
  integer *givptr; 
  VALUE rb_givcol;
  integer *givcol; 
  VALUE rb_givnum;
  doublereal *givnum; 
  VALUE rb_indxq;
  integer *indxq; 
  VALUE rb_info;
  integer info; 
  VALUE rb_d_out__;
  doublereal *d_out__;
  VALUE rb_q_out__;
  doublecomplex *q_out__;
  VALUE rb_qstore_out__;
  doublereal *qstore_out__;
  VALUE rb_qptr_out__;
  integer *qptr_out__;
  doublecomplex *work;
  doublereal *rwork;
  integer *iwork;

  integer n;
  integer ldq;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  indxq, info, d, q, qstore, qptr = NumRu::Lapack.zlaed7( cutpnt, qsiz, tlvls, curlvl, curpbm, d, q, rho, qstore, qptr, prmptr, perm, givptr, givcol, givnum)\n    or\n  NumRu::Lapack.zlaed7  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 15)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 15)", argc);
  rb_cutpnt = argv[0];
  rb_qsiz = argv[1];
  rb_tlvls = argv[2];
  rb_curlvl = argv[3];
  rb_curpbm = argv[4];
  rb_d = argv[5];
  rb_q = argv[6];
  rb_rho = argv[7];
  rb_qstore = argv[8];
  rb_qptr = argv[9];
  rb_prmptr = argv[10];
  rb_perm = argv[11];
  rb_givptr = argv[12];
  rb_givcol = argv[13];
  rb_givnum = argv[14];

  qsiz = NUM2INT(rb_qsiz);
  cutpnt = NUM2INT(rb_cutpnt);
  tlvls = NUM2INT(rb_tlvls);
  rho = NUM2DBL(rb_rho);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (6th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (6th argument) must be %d", 1);
  n = NA_SHAPE0(rb_d);
  if (NA_TYPE(rb_d) != NA_DFLOAT)
    rb_d = na_change_type(rb_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rb_d, doublereal*);
  curlvl = NUM2INT(rb_curlvl);
  curpbm = NUM2INT(rb_curpbm);
  if (!NA_IsNArray(rb_q))
    rb_raise(rb_eArgError, "q (7th argument) must be NArray");
  if (NA_RANK(rb_q) != 2)
    rb_raise(rb_eArgError, "rank of q (7th argument) must be %d", 2);
  if (NA_SHAPE1(rb_q) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of q must be the same as shape 0 of d");
  ldq = NA_SHAPE0(rb_q);
  if (NA_TYPE(rb_q) != NA_DCOMPLEX)
    rb_q = na_change_type(rb_q, NA_DCOMPLEX);
  q = NA_PTR_TYPE(rb_q, doublecomplex*);
  if (!NA_IsNArray(rb_perm))
    rb_raise(rb_eArgError, "perm (12th argument) must be NArray");
  if (NA_RANK(rb_perm) != 1)
    rb_raise(rb_eArgError, "rank of perm (12th argument) must be %d", 1);
  if (NA_SHAPE0(rb_perm) != (n*LG(n)))
    rb_raise(rb_eRuntimeError, "shape 0 of perm must be %d", n*LG(n));
  if (NA_TYPE(rb_perm) != NA_LINT)
    rb_perm = na_change_type(rb_perm, NA_LINT);
  perm = NA_PTR_TYPE(rb_perm, integer*);
  if (!NA_IsNArray(rb_prmptr))
    rb_raise(rb_eArgError, "prmptr (11th argument) must be NArray");
  if (NA_RANK(rb_prmptr) != 1)
    rb_raise(rb_eArgError, "rank of prmptr (11th argument) must be %d", 1);
  if (NA_SHAPE0(rb_prmptr) != (n*LG(n)))
    rb_raise(rb_eRuntimeError, "shape 0 of prmptr must be %d", n*LG(n));
  if (NA_TYPE(rb_prmptr) != NA_LINT)
    rb_prmptr = na_change_type(rb_prmptr, NA_LINT);
  prmptr = NA_PTR_TYPE(rb_prmptr, integer*);
  if (!NA_IsNArray(rb_qstore))
    rb_raise(rb_eArgError, "qstore (9th argument) must be NArray");
  if (NA_RANK(rb_qstore) != 1)
    rb_raise(rb_eArgError, "rank of qstore (9th argument) must be %d", 1);
  if (NA_SHAPE0(rb_qstore) != (pow(n,2)+1))
    rb_raise(rb_eRuntimeError, "shape 0 of qstore must be %d", pow(n,2)+1);
  if (NA_TYPE(rb_qstore) != NA_DFLOAT)
    rb_qstore = na_change_type(rb_qstore, NA_DFLOAT);
  qstore = NA_PTR_TYPE(rb_qstore, doublereal*);
  if (!NA_IsNArray(rb_givptr))
    rb_raise(rb_eArgError, "givptr (13th argument) must be NArray");
  if (NA_RANK(rb_givptr) != 1)
    rb_raise(rb_eArgError, "rank of givptr (13th argument) must be %d", 1);
  if (NA_SHAPE0(rb_givptr) != (n*LG(n)))
    rb_raise(rb_eRuntimeError, "shape 0 of givptr must be %d", n*LG(n));
  if (NA_TYPE(rb_givptr) != NA_LINT)
    rb_givptr = na_change_type(rb_givptr, NA_LINT);
  givptr = NA_PTR_TYPE(rb_givptr, integer*);
  if (!NA_IsNArray(rb_givcol))
    rb_raise(rb_eArgError, "givcol (14th argument) must be NArray");
  if (NA_RANK(rb_givcol) != 2)
    rb_raise(rb_eArgError, "rank of givcol (14th argument) must be %d", 2);
  if (NA_SHAPE1(rb_givcol) != (n*LG(n)))
    rb_raise(rb_eRuntimeError, "shape 1 of givcol must be %d", n*LG(n));
  if (NA_SHAPE0(rb_givcol) != (2))
    rb_raise(rb_eRuntimeError, "shape 0 of givcol must be %d", 2);
  if (NA_TYPE(rb_givcol) != NA_LINT)
    rb_givcol = na_change_type(rb_givcol, NA_LINT);
  givcol = NA_PTR_TYPE(rb_givcol, integer*);
  if (!NA_IsNArray(rb_givnum))
    rb_raise(rb_eArgError, "givnum (15th argument) must be NArray");
  if (NA_RANK(rb_givnum) != 2)
    rb_raise(rb_eArgError, "rank of givnum (15th argument) must be %d", 2);
  if (NA_SHAPE1(rb_givnum) != (n*LG(n)))
    rb_raise(rb_eRuntimeError, "shape 1 of givnum must be %d", n*LG(n));
  if (NA_SHAPE0(rb_givnum) != (2))
    rb_raise(rb_eRuntimeError, "shape 0 of givnum must be %d", 2);
  if (NA_TYPE(rb_givnum) != NA_DFLOAT)
    rb_givnum = na_change_type(rb_givnum, NA_DFLOAT);
  givnum = NA_PTR_TYPE(rb_givnum, doublereal*);
  if (!NA_IsNArray(rb_qptr))
    rb_raise(rb_eArgError, "qptr (10th argument) must be NArray");
  if (NA_RANK(rb_qptr) != 1)
    rb_raise(rb_eArgError, "rank of qptr (10th argument) must be %d", 1);
  if (NA_SHAPE0(rb_qptr) != (n+2))
    rb_raise(rb_eRuntimeError, "shape 0 of qptr must be %d", n+2);
  if (NA_TYPE(rb_qptr) != NA_LINT)
    rb_qptr = na_change_type(rb_qptr, NA_LINT);
  qptr = NA_PTR_TYPE(rb_qptr, integer*);
  {
    int shape[1];
    shape[0] = n;
    rb_indxq = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  indxq = NA_PTR_TYPE(rb_indxq, integer*);
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
    shape[0] = pow(n,2)+1;
    rb_qstore_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  qstore_out__ = NA_PTR_TYPE(rb_qstore_out__, doublereal*);
  MEMCPY(qstore_out__, qstore, doublereal, NA_TOTAL(rb_qstore));
  rb_qstore = rb_qstore_out__;
  qstore = qstore_out__;
  {
    int shape[1];
    shape[0] = n+2;
    rb_qptr_out__ = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  qptr_out__ = NA_PTR_TYPE(rb_qptr_out__, integer*);
  MEMCPY(qptr_out__, qptr, integer, NA_TOTAL(rb_qptr));
  rb_qptr = rb_qptr_out__;
  qptr = qptr_out__;
  work = ALLOC_N(doublecomplex, (qsiz*n));
  rwork = ALLOC_N(doublereal, (3*n+2*qsiz*n));
  iwork = ALLOC_N(integer, (4*n));

  zlaed7_(&n, &cutpnt, &qsiz, &tlvls, &curlvl, &curpbm, d, q, &ldq, &rho, indxq, qstore, qptr, prmptr, perm, givptr, givcol, givnum, work, rwork, iwork, &info);

  free(work);
  free(rwork);
  free(iwork);
  rb_info = INT2NUM(info);
  return rb_ary_new3(6, rb_indxq, rb_info, rb_d, rb_q, rb_qstore, rb_qptr);
}

void
init_lapack_zlaed7(VALUE mLapack){
  rb_define_module_function(mLapack, "zlaed7", rb_zlaed7, -1);
}
