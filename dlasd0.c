#include "rb_lapack.h"

extern VOID dlasd0_(integer *n, integer *sqre, doublereal *d, doublereal *e, doublereal *u, integer *ldu, doublereal *vt, integer *ldvt, integer *smlsiz, integer *iwork, doublereal *work, integer *info);

static VALUE
rb_dlasd0(int argc, VALUE *argv, VALUE self){
  VALUE rb_sqre;
  integer sqre; 
  VALUE rb_d;
  doublereal *d; 
  VALUE rb_e;
  doublereal *e; 
  VALUE rb_smlsiz;
  integer smlsiz; 
  VALUE rb_u;
  doublereal *u; 
  VALUE rb_vt;
  doublereal *vt; 
  VALUE rb_info;
  integer info; 
  VALUE rb_d_out__;
  doublereal *d_out__;
  integer *iwork;
  doublereal *work;

  integer n;
  integer ldu;
  integer ldvt;
  integer m;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  u, vt, info, d = NumRu::Lapack.dlasd0( sqre, d, e, smlsiz)\n    or\n  NumRu::Lapack.dlasd0  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLASD0( N, SQRE, D, E, U, LDU, VT, LDVT, SMLSIZ, IWORK, WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  Using a divide and conquer approach, DLASD0 computes the singular\n*  value decomposition (SVD) of a real upper bidiagonal N-by-M\n*  matrix B with diagonal D and offdiagonal E, where M = N + SQRE.\n*  The algorithm computes orthogonal matrices U and VT such that\n*  B = U * S * VT. The singular values S are overwritten on D.\n*\n*  A related subroutine, DLASDA, computes only the singular values,\n*  and optionally, the singular vectors in compact form.\n*\n\n*  Arguments\n*  =========\n*\n*  N      (input) INTEGER\n*         On entry, the row dimension of the upper bidiagonal matrix.\n*         This is also the dimension of the main diagonal array D.\n*\n*  SQRE   (input) INTEGER\n*         Specifies the column dimension of the bidiagonal matrix.\n*         = 0: The bidiagonal matrix has column dimension M = N;\n*         = 1: The bidiagonal matrix has column dimension M = N+1;\n*\n*  D      (input/output) DOUBLE PRECISION array, dimension (N)\n*         On entry D contains the main diagonal of the bidiagonal\n*         matrix.\n*         On exit D, if INFO = 0, contains its singular values.\n*\n*  E      (input) DOUBLE PRECISION array, dimension (M-1)\n*         Contains the subdiagonal entries of the bidiagonal matrix.\n*         On exit, E has been destroyed.\n*\n*  U      (output) DOUBLE PRECISION array, dimension at least (LDQ, N)\n*         On exit, U contains the left singular vectors.\n*\n*  LDU    (input) INTEGER\n*         On entry, leading dimension of U.\n*\n*  VT     (output) DOUBLE PRECISION array, dimension at least (LDVT, M)\n*         On exit, VT' contains the right singular vectors.\n*\n*  LDVT   (input) INTEGER\n*         On entry, leading dimension of VT.\n*\n*  SMLSIZ (input) INTEGER\n*         On entry, maximum size of the subproblems at the\n*         bottom of the computation tree.\n*\n*  IWORK  (workspace) INTEGER work array.\n*         Dimension must be at least (8 * N)\n*\n*  WORK   (workspace) DOUBLE PRECISION work array.\n*         Dimension must be at least (3 * M**2 + 2 * M)\n*\n*  INFO   (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*          > 0:  if INFO = 1, a singular value did not converge\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Ming Gu and Huan Ren, Computer Science Division, University of\n*     California at Berkeley, USA\n*\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      INTEGER            I, I1, IC, IDXQ, IDXQC, IM1, INODE, ITEMP, IWK,\n     $                   J, LF, LL, LVL, M, NCC, ND, NDB1, NDIML, NDIMR,\n     $                   NL, NLF, NLP1, NLVL, NR, NRF, NRP1, SQREI\n      DOUBLE PRECISION   ALPHA, BETA\n*     ..\n*     .. External Subroutines ..\n      EXTERNAL           DLASD1, DLASDQ, DLASDT, XERBLA\n*     ..\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_sqre = argv[0];
  rb_d = argv[1];
  rb_e = argv[2];
  rb_smlsiz = argv[3];

  smlsiz = NUM2INT(rb_smlsiz);
  sqre = NUM2INT(rb_sqre);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (2th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (2th argument) must be %d", 1);
  n = NA_SHAPE0(rb_d);
  if (NA_TYPE(rb_d) != NA_DFLOAT)
    rb_d = na_change_type(rb_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rb_d, doublereal*);
  ldu = n;
  ldvt = n;
  m = sqre == 0 ? n : sqre == 1 ? n+1 : 0;
  if (!NA_IsNArray(rb_e))
    rb_raise(rb_eArgError, "e (3th argument) must be NArray");
  if (NA_RANK(rb_e) != 1)
    rb_raise(rb_eArgError, "rank of e (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_e) != (m-1))
    rb_raise(rb_eRuntimeError, "shape 0 of e must be %d", m-1);
  if (NA_TYPE(rb_e) != NA_DFLOAT)
    rb_e = na_change_type(rb_e, NA_DFLOAT);
  e = NA_PTR_TYPE(rb_e, doublereal*);
  {
    int shape[2];
    shape[0] = ldu;
    shape[1] = n;
    rb_u = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  u = NA_PTR_TYPE(rb_u, doublereal*);
  {
    int shape[2];
    shape[0] = ldvt;
    shape[1] = m;
    rb_vt = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  vt = NA_PTR_TYPE(rb_vt, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_d_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  d_out__ = NA_PTR_TYPE(rb_d_out__, doublereal*);
  MEMCPY(d_out__, d, doublereal, NA_TOTAL(rb_d));
  rb_d = rb_d_out__;
  d = d_out__;
  iwork = ALLOC_N(integer, ((8 * n)));
  work = ALLOC_N(doublereal, ((3 * pow(m,2) + 2 * m)));

  dlasd0_(&n, &sqre, d, e, u, &ldu, vt, &ldvt, &smlsiz, iwork, work, &info);

  free(iwork);
  free(work);
  rb_info = INT2NUM(info);
  return rb_ary_new3(4, rb_u, rb_vt, rb_info, rb_d);
}

void
init_lapack_dlasd0(VALUE mLapack){
  rb_define_module_function(mLapack, "dlasd0", rb_dlasd0, -1);
}
