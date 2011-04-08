#include "rb_lapack.h"

extern VOID slasd0_(integer *n, integer *sqre, real *d, real *e, real *u, integer *ldu, real *vt, integer *ldvt, integer *smlsiz, integer *iwork, real *work, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_slasd0(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_sqre;
  integer sqre; 
  VALUE rblapack_d;
  real *d; 
  VALUE rblapack_e;
  real *e; 
  VALUE rblapack_smlsiz;
  integer smlsiz; 
  VALUE rblapack_u;
  real *u; 
  VALUE rblapack_vt;
  real *vt; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_d_out__;
  real *d_out__;
  integer *iwork;
  real *work;

  integer n;
  integer ldu;
  integer ldvt;
  integer m;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  u, vt, info, d = NumRu::Lapack.slasd0( sqre, d, e, smlsiz, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLASD0( N, SQRE, D, E, U, LDU, VT, LDVT, SMLSIZ, IWORK, WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  Using a divide and conquer approach, SLASD0 computes the singular\n*  value decomposition (SVD) of a real upper bidiagonal N-by-M\n*  matrix B with diagonal D and offdiagonal E, where M = N + SQRE.\n*  The algorithm computes orthogonal matrices U and VT such that\n*  B = U * S * VT. The singular values S are overwritten on D.\n*\n*  A related subroutine, SLASDA, computes only the singular values,\n*  and optionally, the singular vectors in compact form.\n*\n\n*  Arguments\n*  =========\n*\n*  N      (input) INTEGER\n*         On entry, the row dimension of the upper bidiagonal matrix.\n*         This is also the dimension of the main diagonal array D.\n*\n*  SQRE   (input) INTEGER\n*         Specifies the column dimension of the bidiagonal matrix.\n*         = 0: The bidiagonal matrix has column dimension M = N;\n*         = 1: The bidiagonal matrix has column dimension M = N+1;\n*\n*  D      (input/output) REAL array, dimension (N)\n*         On entry D contains the main diagonal of the bidiagonal\n*         matrix.\n*         On exit D, if INFO = 0, contains its singular values.\n*\n*  E      (input) REAL array, dimension (M-1)\n*         Contains the subdiagonal entries of the bidiagonal matrix.\n*         On exit, E has been destroyed.\n*\n*  U      (output) REAL array, dimension at least (LDQ, N)\n*         On exit, U contains the left singular vectors.\n*\n*  LDU    (input) INTEGER\n*         On entry, leading dimension of U.\n*\n*  VT     (output) REAL array, dimension at least (LDVT, M)\n*         On exit, VT' contains the right singular vectors.\n*\n*  LDVT   (input) INTEGER\n*         On entry, leading dimension of VT.\n*\n*  SMLSIZ (input) INTEGER\n*         On entry, maximum size of the subproblems at the\n*         bottom of the computation tree.\n*\n*  IWORK  (workspace) INTEGER array, dimension (8*N)\n*\n*  WORK   (workspace) REAL array, dimension (3*M**2+2*M)\n*\n*  INFO   (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*          > 0:  if INFO = 1, a singular value did not converge\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Ming Gu and Huan Ren, Computer Science Division, University of\n*     California at Berkeley, USA\n*\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      INTEGER            I, I1, IC, IDXQ, IDXQC, IM1, INODE, ITEMP, IWK,\n     $                   J, LF, LL, LVL, M, NCC, ND, NDB1, NDIML, NDIMR,\n     $                   NL, NLF, NLP1, NLVL, NR, NRF, NRP1, SQREI\n      REAL               ALPHA, BETA\n*     ..\n*     .. External Subroutines ..\n      EXTERNAL           SLASD1, SLASDQ, SLASDT, XERBLA\n*     ..\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  u, vt, info, d = NumRu::Lapack.slasd0( sqre, d, e, smlsiz, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rblapack_sqre = argv[0];
  rblapack_d = argv[1];
  rblapack_e = argv[2];
  rblapack_smlsiz = argv[3];
  if (rb_options != Qnil) {
  }

  smlsiz = NUM2INT(rblapack_smlsiz);
  sqre = NUM2INT(rblapack_sqre);
  if (!NA_IsNArray(rblapack_d))
    rb_raise(rb_eArgError, "d (2th argument) must be NArray");
  if (NA_RANK(rblapack_d) != 1)
    rb_raise(rb_eArgError, "rank of d (2th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_d);
  if (NA_TYPE(rblapack_d) != NA_SFLOAT)
    rblapack_d = na_change_type(rblapack_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rblapack_d, real*);
  ldu = n;
  m = sqre == 0 ? n : sqre == 1 ? n+1 : 0;
  ldvt = m;
  if (!NA_IsNArray(rblapack_e))
    rb_raise(rb_eArgError, "e (3th argument) must be NArray");
  if (NA_RANK(rblapack_e) != 1)
    rb_raise(rb_eArgError, "rank of e (3th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_e) != (m-1))
    rb_raise(rb_eRuntimeError, "shape 0 of e must be %d", m-1);
  if (NA_TYPE(rblapack_e) != NA_SFLOAT)
    rblapack_e = na_change_type(rblapack_e, NA_SFLOAT);
  e = NA_PTR_TYPE(rblapack_e, real*);
  {
    int shape[2];
    shape[0] = ldu;
    shape[1] = n;
    rblapack_u = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  u = NA_PTR_TYPE(rblapack_u, real*);
  {
    int shape[2];
    shape[0] = ldvt;
    shape[1] = m;
    rblapack_vt = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  vt = NA_PTR_TYPE(rblapack_vt, real*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_d_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  d_out__ = NA_PTR_TYPE(rblapack_d_out__, real*);
  MEMCPY(d_out__, d, real, NA_TOTAL(rblapack_d));
  rblapack_d = rblapack_d_out__;
  d = d_out__;
  iwork = ALLOC_N(integer, (8*n));
  work = ALLOC_N(real, (3*pow(m,2)+2*m));

  slasd0_(&n, &sqre, d, e, u, &ldu, vt, &ldvt, &smlsiz, iwork, work, &info);

  free(iwork);
  free(work);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(4, rblapack_u, rblapack_vt, rblapack_info, rblapack_d);
}

void
init_lapack_slasd0(VALUE mLapack){
  rb_define_module_function(mLapack, "slasd0", rblapack_slasd0, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
