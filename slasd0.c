#include "rb_lapack.h"

static VALUE
rb_slasd0(int argc, VALUE *argv, VALUE self){
  VALUE rb_sqre;
  integer sqre; 
  VALUE rb_d;
  real *d; 
  VALUE rb_e;
  real *e; 
  VALUE rb_smlsiz;
  integer smlsiz; 
  VALUE rb_u;
  real *u; 
  VALUE rb_vt;
  real *vt; 
  VALUE rb_info;
  integer info; 
  VALUE rb_d_out__;
  real *d_out__;
  integer *iwork;
  real *work;

  integer n;
  integer m;
  integer ldu;
  integer ldvt;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  u, vt, info, d = NumRu::Lapack.slasd0( sqre, d, e, smlsiz)\n    or\n  NumRu::Lapack.slasd0  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLASD0( N, SQRE, D, E, U, LDU, VT, LDVT, SMLSIZ, IWORK, WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  Using a divide and conquer approach, SLASD0 computes the singular\n*  value decomposition (SVD) of a real upper bidiagonal N-by-M\n*  matrix B with diagonal D and offdiagonal E, where M = N + SQRE.\n*  The algorithm computes orthogonal matrices U and VT such that\n*  B = U * S * VT. The singular values S are overwritten on D.\n*\n*  A related subroutine, SLASDA, computes only the singular values,\n*  and optionally, the singular vectors in compact form.\n*\n\n*  Arguments\n*  =========\n*\n*  N      (input) INTEGER\n*         On entry, the row dimension of the upper bidiagonal matrix.\n*         This is also the dimension of the main diagonal array D.\n*\n*  SQRE   (input) INTEGER\n*         Specifies the column dimension of the bidiagonal matrix.\n*         = 0: The bidiagonal matrix has column dimension M = N;\n*         = 1: The bidiagonal matrix has column dimension M = N+1;\n*\n*  D      (input/output) REAL array, dimension (N)\n*         On entry D contains the main diagonal of the bidiagonal\n*         matrix.\n*         On exit D, if INFO = 0, contains its singular values.\n*\n*  E      (input) REAL array, dimension (M-1)\n*         Contains the subdiagonal entries of the bidiagonal matrix.\n*         On exit, E has been destroyed.\n*\n*  U      (output) REAL array, dimension at least (LDQ, N)\n*         On exit, U contains the left singular vectors.\n*\n*  LDU    (input) INTEGER\n*         On entry, leading dimension of U.\n*\n*  VT     (output) REAL array, dimension at least (LDVT, M)\n*         On exit, VT' contains the right singular vectors.\n*\n*  LDVT   (input) INTEGER\n*         On entry, leading dimension of VT.\n*\n*  SMLSIZ (input) INTEGER\n*         On entry, maximum size of the subproblems at the\n*         bottom of the computation tree.\n*\n*  IWORK  (workspace) INTEGER array, dimension (8*N)\n*\n*  WORK   (workspace) REAL array, dimension (3*M**2+2*M)\n*\n*  INFO   (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*          > 0:  if INFO = 1, an singular value did not converge\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Ming Gu and Huan Ren, Computer Science Division, University of\n*     California at Berkeley, USA\n*\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      INTEGER            I, I1, IC, IDXQ, IDXQC, IM1, INODE, ITEMP, IWK,\n     $                   J, LF, LL, LVL, M, NCC, ND, NDB1, NDIML, NDIMR,\n     $                   NL, NLF, NLP1, NLVL, NR, NRF, NRP1, SQREI\n      REAL               ALPHA, BETA\n*     ..\n*     .. External Subroutines ..\n      EXTERNAL           SLASD1, SLASDQ, SLASDT, XERBLA\n*     ..\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_sqre = argv[0];
  rb_d = argv[1];
  rb_e = argv[2];
  rb_smlsiz = argv[3];

  sqre = NUM2INT(rb_sqre);
  smlsiz = NUM2INT(rb_smlsiz);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (2th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (2th argument) must be %d", 1);
  n = NA_SHAPE0(rb_d);
  if (NA_TYPE(rb_d) != NA_SFLOAT)
    rb_d = na_change_type(rb_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rb_d, real*);
  if (!NA_IsNArray(rb_e))
    rb_raise(rb_eArgError, "e (3th argument) must be NArray");
  if (NA_RANK(rb_e) != 1)
    rb_raise(rb_eArgError, "rank of e (3th argument) must be %d", 1);
  m = sqre == 0 ? n : sqre == 1 ? n+1 : 0;
  ldu = n;
  if (NA_SHAPE0(rb_e) != (m-1))
    rb_raise(rb_eRuntimeError, "shape 0 of e must be %d", m-1);
  if (NA_TYPE(rb_e) != NA_SFLOAT)
    rb_e = na_change_type(rb_e, NA_SFLOAT);
  e = NA_PTR_TYPE(rb_e, real*);
  {
    int shape[2];
    shape[0] = ldu;
    shape[1] = n;
    rb_u = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  u = NA_PTR_TYPE(rb_u, real*);
  ldvt = m;
  {
    int shape[2];
    shape[0] = ldvt;
    shape[1] = m;
    rb_vt = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  vt = NA_PTR_TYPE(rb_vt, real*);
  {
    int shape[1];
    shape[0] = n;
    rb_d_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  d_out__ = NA_PTR_TYPE(rb_d_out__, real*);
  MEMCPY(d_out__, d, real, NA_TOTAL(rb_d));
  rb_d = rb_d_out__;
  d = d_out__;
  iwork = ALLOC_N(integer, (8*n));
  work = ALLOC_N(real, (3*pow(m,2)+2*m));

  slasd0_(&n, &sqre, d, e, u, &ldu, vt, &ldvt, &smlsiz, iwork, work, &info);

  free(iwork);
  free(work);
  rb_info = INT2NUM(info);
  return rb_ary_new3(4, rb_u, rb_vt, rb_info, rb_d);
}

void
init_lapack_slasd0(VALUE mLapack){
  rb_define_module_function(mLapack, "slasd0", rb_slasd0, -1);
}
