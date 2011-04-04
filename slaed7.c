#include "rb_lapack.h"

extern VOID slaed7_(integer *icompq, integer *n, integer *qsiz, integer *tlvls, integer *curlvl, integer *curpbm, real *d, real *q, integer *ldq, integer *indxq, real *rho, integer *cutpnt, real *qstore, integer *qptr, integer *prmptr, integer *perm, integer *givptr, integer *givcol, real *givnum, real *work, integer *iwork, integer *info);

static VALUE
rb_slaed7(int argc, VALUE *argv, VALUE self){
  VALUE rb_icompq;
  integer icompq; 
  VALUE rb_qsiz;
  integer qsiz; 
  VALUE rb_tlvls;
  integer tlvls; 
  VALUE rb_curlvl;
  integer curlvl; 
  VALUE rb_curpbm;
  integer curpbm; 
  VALUE rb_d;
  real *d; 
  VALUE rb_q;
  real *q; 
  VALUE rb_rho;
  real rho; 
  VALUE rb_cutpnt;
  integer cutpnt; 
  VALUE rb_qstore;
  real *qstore; 
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
  real *givnum; 
  VALUE rb_indxq;
  integer *indxq; 
  VALUE rb_info;
  integer info; 
  VALUE rb_d_out__;
  real *d_out__;
  VALUE rb_q_out__;
  real *q_out__;
  VALUE rb_qstore_out__;
  real *qstore_out__;
  VALUE rb_qptr_out__;
  integer *qptr_out__;
  real *work;
  integer *iwork;

  integer n;
  integer ldq;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  indxq, info, d, q, qstore, qptr = NumRu::Lapack.slaed7( icompq, qsiz, tlvls, curlvl, curpbm, d, q, rho, cutpnt, qstore, qptr, prmptr, perm, givptr, givcol, givnum)\n    or\n  NumRu::Lapack.slaed7  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 16)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 16)", argc);
  rb_icompq = argv[0];
  rb_qsiz = argv[1];
  rb_tlvls = argv[2];
  rb_curlvl = argv[3];
  rb_curpbm = argv[4];
  rb_d = argv[5];
  rb_q = argv[6];
  rb_rho = argv[7];
  rb_cutpnt = argv[8];
  rb_qstore = argv[9];
  rb_qptr = argv[10];
  rb_prmptr = argv[11];
  rb_perm = argv[12];
  rb_givptr = argv[13];
  rb_givcol = argv[14];
  rb_givnum = argv[15];

  qsiz = NUM2INT(rb_qsiz);
  cutpnt = NUM2INT(rb_cutpnt);
  tlvls = NUM2INT(rb_tlvls);
  rho = (real)NUM2DBL(rb_rho);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (6th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (6th argument) must be %d", 1);
  n = NA_SHAPE0(rb_d);
  if (NA_TYPE(rb_d) != NA_SFLOAT)
    rb_d = na_change_type(rb_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rb_d, real*);
  curlvl = NUM2INT(rb_curlvl);
  icompq = NUM2INT(rb_icompq);
  curpbm = NUM2INT(rb_curpbm);
  if (!NA_IsNArray(rb_q))
    rb_raise(rb_eArgError, "q (7th argument) must be NArray");
  if (NA_RANK(rb_q) != 2)
    rb_raise(rb_eArgError, "rank of q (7th argument) must be %d", 2);
  if (NA_SHAPE1(rb_q) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of q must be the same as shape 0 of d");
  ldq = NA_SHAPE0(rb_q);
  if (NA_TYPE(rb_q) != NA_SFLOAT)
    rb_q = na_change_type(rb_q, NA_SFLOAT);
  q = NA_PTR_TYPE(rb_q, real*);
  if (!NA_IsNArray(rb_perm))
    rb_raise(rb_eArgError, "perm (13th argument) must be NArray");
  if (NA_RANK(rb_perm) != 1)
    rb_raise(rb_eArgError, "rank of perm (13th argument) must be %d", 1);
  if (NA_SHAPE0(rb_perm) != (n*LG(n)))
    rb_raise(rb_eRuntimeError, "shape 0 of perm must be %d", n*LG(n));
  if (NA_TYPE(rb_perm) != NA_LINT)
    rb_perm = na_change_type(rb_perm, NA_LINT);
  perm = NA_PTR_TYPE(rb_perm, integer*);
  if (!NA_IsNArray(rb_prmptr))
    rb_raise(rb_eArgError, "prmptr (12th argument) must be NArray");
  if (NA_RANK(rb_prmptr) != 1)
    rb_raise(rb_eArgError, "rank of prmptr (12th argument) must be %d", 1);
  if (NA_SHAPE0(rb_prmptr) != (n*LG(n)))
    rb_raise(rb_eRuntimeError, "shape 0 of prmptr must be %d", n*LG(n));
  if (NA_TYPE(rb_prmptr) != NA_LINT)
    rb_prmptr = na_change_type(rb_prmptr, NA_LINT);
  prmptr = NA_PTR_TYPE(rb_prmptr, integer*);
  if (!NA_IsNArray(rb_qstore))
    rb_raise(rb_eArgError, "qstore (10th argument) must be NArray");
  if (NA_RANK(rb_qstore) != 1)
    rb_raise(rb_eArgError, "rank of qstore (10th argument) must be %d", 1);
  if (NA_SHAPE0(rb_qstore) != (pow(n,2)+1))
    rb_raise(rb_eRuntimeError, "shape 0 of qstore must be %d", pow(n,2)+1);
  if (NA_TYPE(rb_qstore) != NA_SFLOAT)
    rb_qstore = na_change_type(rb_qstore, NA_SFLOAT);
  qstore = NA_PTR_TYPE(rb_qstore, real*);
  if (!NA_IsNArray(rb_givptr))
    rb_raise(rb_eArgError, "givptr (14th argument) must be NArray");
  if (NA_RANK(rb_givptr) != 1)
    rb_raise(rb_eArgError, "rank of givptr (14th argument) must be %d", 1);
  if (NA_SHAPE0(rb_givptr) != (n*LG(n)))
    rb_raise(rb_eRuntimeError, "shape 0 of givptr must be %d", n*LG(n));
  if (NA_TYPE(rb_givptr) != NA_LINT)
    rb_givptr = na_change_type(rb_givptr, NA_LINT);
  givptr = NA_PTR_TYPE(rb_givptr, integer*);
  if (!NA_IsNArray(rb_givcol))
    rb_raise(rb_eArgError, "givcol (15th argument) must be NArray");
  if (NA_RANK(rb_givcol) != 2)
    rb_raise(rb_eArgError, "rank of givcol (15th argument) must be %d", 2);
  if (NA_SHAPE1(rb_givcol) != (n*LG(n)))
    rb_raise(rb_eRuntimeError, "shape 1 of givcol must be %d", n*LG(n));
  if (NA_SHAPE0(rb_givcol) != (2))
    rb_raise(rb_eRuntimeError, "shape 0 of givcol must be %d", 2);
  if (NA_TYPE(rb_givcol) != NA_LINT)
    rb_givcol = na_change_type(rb_givcol, NA_LINT);
  givcol = NA_PTR_TYPE(rb_givcol, integer*);
  if (!NA_IsNArray(rb_givnum))
    rb_raise(rb_eArgError, "givnum (16th argument) must be NArray");
  if (NA_RANK(rb_givnum) != 2)
    rb_raise(rb_eArgError, "rank of givnum (16th argument) must be %d", 2);
  if (NA_SHAPE1(rb_givnum) != (n*LG(n)))
    rb_raise(rb_eRuntimeError, "shape 1 of givnum must be %d", n*LG(n));
  if (NA_SHAPE0(rb_givnum) != (2))
    rb_raise(rb_eRuntimeError, "shape 0 of givnum must be %d", 2);
  if (NA_TYPE(rb_givnum) != NA_SFLOAT)
    rb_givnum = na_change_type(rb_givnum, NA_SFLOAT);
  givnum = NA_PTR_TYPE(rb_givnum, real*);
  if (!NA_IsNArray(rb_qptr))
    rb_raise(rb_eArgError, "qptr (11th argument) must be NArray");
  if (NA_RANK(rb_qptr) != 1)
    rb_raise(rb_eArgError, "rank of qptr (11th argument) must be %d", 1);
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
    rb_d_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  d_out__ = NA_PTR_TYPE(rb_d_out__, real*);
  MEMCPY(d_out__, d, real, NA_TOTAL(rb_d));
  rb_d = rb_d_out__;
  d = d_out__;
  {
    int shape[2];
    shape[0] = ldq;
    shape[1] = n;
    rb_q_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  q_out__ = NA_PTR_TYPE(rb_q_out__, real*);
  MEMCPY(q_out__, q, real, NA_TOTAL(rb_q));
  rb_q = rb_q_out__;
  q = q_out__;
  {
    int shape[1];
    shape[0] = pow(n,2)+1;
    rb_qstore_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  qstore_out__ = NA_PTR_TYPE(rb_qstore_out__, real*);
  MEMCPY(qstore_out__, qstore, real, NA_TOTAL(rb_qstore));
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
  work = ALLOC_N(real, (3*n+qsiz*n));
  iwork = ALLOC_N(integer, (4*n));

  slaed7_(&icompq, &n, &qsiz, &tlvls, &curlvl, &curpbm, d, q, &ldq, indxq, &rho, &cutpnt, qstore, qptr, prmptr, perm, givptr, givcol, givnum, work, iwork, &info);

  free(work);
  free(iwork);
  rb_info = INT2NUM(info);
  return rb_ary_new3(6, rb_indxq, rb_info, rb_d, rb_q, rb_qstore, rb_qptr);
}

void
init_lapack_slaed7(VALUE mLapack){
  rb_define_module_function(mLapack, "slaed7", rb_slaed7, -1);
}
