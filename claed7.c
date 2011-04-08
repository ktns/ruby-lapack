#include "rb_lapack.h"

extern VOID claed7_(integer *n, integer *cutpnt, integer *qsiz, integer *tlvls, integer *curlvl, integer *curpbm, real *d, complex *q, integer *ldq, real *rho, integer *indxq, real *qstore, integer *qptr, integer *prmptr, integer *perm, integer *givptr, integer *givcol, real *givnum, complex *work, real *rwork, integer *iwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_claed7(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_cutpnt;
  integer cutpnt; 
  VALUE rblapack_qsiz;
  integer qsiz; 
  VALUE rblapack_tlvls;
  integer tlvls; 
  VALUE rblapack_curlvl;
  integer curlvl; 
  VALUE rblapack_curpbm;
  integer curpbm; 
  VALUE rblapack_d;
  real *d; 
  VALUE rblapack_q;
  complex *q; 
  VALUE rblapack_rho;
  real rho; 
  VALUE rblapack_qstore;
  real *qstore; 
  VALUE rblapack_qptr;
  integer *qptr; 
  VALUE rblapack_prmptr;
  integer *prmptr; 
  VALUE rblapack_perm;
  integer *perm; 
  VALUE rblapack_givptr;
  integer *givptr; 
  VALUE rblapack_givcol;
  integer *givcol; 
  VALUE rblapack_givnum;
  real *givnum; 
  VALUE rblapack_indxq;
  integer *indxq; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_d_out__;
  real *d_out__;
  VALUE rblapack_q_out__;
  complex *q_out__;
  VALUE rblapack_qstore_out__;
  real *qstore_out__;
  VALUE rblapack_qptr_out__;
  integer *qptr_out__;
  complex *work;
  real *rwork;
  integer *iwork;

  integer n;
  integer ldq;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  indxq, info, d, q, qstore, qptr = NumRu::Lapack.claed7( cutpnt, qsiz, tlvls, curlvl, curpbm, d, q, rho, qstore, qptr, prmptr, perm, givptr, givcol, givnum, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE CLAED7( N, CUTPNT, QSIZ, TLVLS, CURLVL, CURPBM, D, Q, LDQ, RHO, INDXQ, QSTORE, QPTR, PRMPTR, PERM, GIVPTR, GIVCOL, GIVNUM, WORK, RWORK, IWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  CLAED7 computes the updated eigensystem of a diagonal\n*  matrix after modification by a rank-one symmetric matrix. This\n*  routine is used only for the eigenproblem which requires all\n*  eigenvalues and optionally eigenvectors of a dense or banded\n*  Hermitian matrix that has been reduced to tridiagonal form.\n*\n*    T = Q(in) ( D(in) + RHO * Z*Z' ) Q'(in) = Q(out) * D(out) * Q'(out)\n*\n*    where Z = Q'u, u is a vector of length N with ones in the\n*    CUTPNT and CUTPNT + 1 th elements and zeros elsewhere.\n*\n*     The eigenvectors of the original matrix are stored in Q, and the\n*     eigenvalues are in D.  The algorithm consists of three stages:\n*\n*        The first stage consists of deflating the size of the problem\n*        when there are multiple eigenvalues or if there is a zero in\n*        the Z vector.  For each such occurence the dimension of the\n*        secular equation problem is reduced by one.  This stage is\n*        performed by the routine SLAED2.\n*\n*        The second stage consists of calculating the updated\n*        eigenvalues. This is done by finding the roots of the secular\n*        equation via the routine SLAED4 (as called by SLAED3).\n*        This routine also calculates the eigenvectors of the current\n*        problem.\n*\n*        The final stage consists of computing the updated eigenvectors\n*        directly using the updated eigenvalues.  The eigenvectors for\n*        the current problem are multiplied with the eigenvectors from\n*        the overall problem.\n*\n\n*  Arguments\n*  =========\n*\n*  N      (input) INTEGER\n*         The dimension of the symmetric tridiagonal matrix.  N >= 0.\n*\n*  CUTPNT (input) INTEGER\n*         Contains the location of the last eigenvalue in the leading\n*         sub-matrix.  min(1,N) <= CUTPNT <= N.\n*\n*  QSIZ   (input) INTEGER\n*         The dimension of the unitary matrix used to reduce\n*         the full matrix to tridiagonal form.  QSIZ >= N.\n*\n*  TLVLS  (input) INTEGER\n*         The total number of merging levels in the overall divide and\n*         conquer tree.\n*\n*  CURLVL (input) INTEGER\n*         The current level in the overall merge routine,\n*         0 <= curlvl <= tlvls.\n*\n*  CURPBM (input) INTEGER\n*         The current problem in the current level in the overall\n*         merge routine (counting from upper left to lower right).\n*\n*  D      (input/output) REAL array, dimension (N)\n*         On entry, the eigenvalues of the rank-1-perturbed matrix.\n*         On exit, the eigenvalues of the repaired matrix.\n*\n*  Q      (input/output) COMPLEX array, dimension (LDQ,N)\n*         On entry, the eigenvectors of the rank-1-perturbed matrix.\n*         On exit, the eigenvectors of the repaired tridiagonal matrix.\n*\n*  LDQ    (input) INTEGER\n*         The leading dimension of the array Q.  LDQ >= max(1,N).\n*\n*  RHO    (input) REAL\n*         Contains the subdiagonal element used to create the rank-1\n*         modification.\n*\n*  INDXQ  (output) INTEGER array, dimension (N)\n*         This contains the permutation which will reintegrate the\n*         subproblem just solved back into sorted order,\n*         ie. D( INDXQ( I = 1, N ) ) will be in ascending order.\n*\n*  IWORK  (workspace) INTEGER array, dimension (4*N)\n*\n*  RWORK  (workspace) REAL array,\n*                                 dimension (3*N+2*QSIZ*N)\n*\n*  WORK   (workspace) COMPLEX array, dimension (QSIZ*N)\n*\n*  QSTORE (input/output) REAL array, dimension (N**2+1)\n*         Stores eigenvectors of submatrices encountered during\n*         divide and conquer, packed together. QPTR points to\n*         beginning of the submatrices.\n*\n*  QPTR   (input/output) INTEGER array, dimension (N+2)\n*         List of indices pointing to beginning of submatrices stored\n*         in QSTORE. The submatrices are numbered starting at the\n*         bottom left of the divide and conquer tree, from left to\n*         right and bottom to top.\n*\n*  PRMPTR (input) INTEGER array, dimension (N lg N)\n*         Contains a list of pointers which indicate where in PERM a\n*         level's permutation is stored.  PRMPTR(i+1) - PRMPTR(i)\n*         indicates the size of the permutation and also the size of\n*         the full, non-deflated problem.\n*\n*  PERM   (input) INTEGER array, dimension (N lg N)\n*         Contains the permutations (from deflation and sorting) to be\n*         applied to each eigenblock.\n*\n*  GIVPTR (input) INTEGER array, dimension (N lg N)\n*         Contains a list of pointers which indicate where in GIVCOL a\n*         level's Givens rotations are stored.  GIVPTR(i+1) - GIVPTR(i)\n*         indicates the number of Givens rotations.\n*\n*  GIVCOL (input) INTEGER array, dimension (2, N lg N)\n*         Each pair of numbers indicates a pair of columns to take place\n*         in a Givens rotation.\n*\n*  GIVNUM (input) REAL array, dimension (2, N lg N)\n*         Each number indicates the S value to be used in the\n*         corresponding Givens rotation.\n*\n*  INFO   (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*          > 0:  if INFO = 1, an eigenvalue did not converge\n*\n\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      INTEGER            COLTYP, CURR, I, IDLMDA, INDX,\n     $                   INDXC, INDXP, IQ, IW, IZ, K, N1, N2, PTR\n*     ..\n*     .. External Subroutines ..\n      EXTERNAL           CLACRM, CLAED8, SLAED9, SLAEDA, SLAMRG, XERBLA\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          MAX, MIN\n*     ..\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  indxq, info, d, q, qstore, qptr = NumRu::Lapack.claed7( cutpnt, qsiz, tlvls, curlvl, curpbm, d, q, rho, qstore, qptr, prmptr, perm, givptr, givcol, givnum, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 15)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 15)", argc);
  rblapack_cutpnt = argv[0];
  rblapack_qsiz = argv[1];
  rblapack_tlvls = argv[2];
  rblapack_curlvl = argv[3];
  rblapack_curpbm = argv[4];
  rblapack_d = argv[5];
  rblapack_q = argv[6];
  rblapack_rho = argv[7];
  rblapack_qstore = argv[8];
  rblapack_qptr = argv[9];
  rblapack_prmptr = argv[10];
  rblapack_perm = argv[11];
  rblapack_givptr = argv[12];
  rblapack_givcol = argv[13];
  rblapack_givnum = argv[14];
  if (rb_options != Qnil) {
  }

  qsiz = NUM2INT(rblapack_qsiz);
  cutpnt = NUM2INT(rblapack_cutpnt);
  tlvls = NUM2INT(rblapack_tlvls);
  rho = (real)NUM2DBL(rblapack_rho);
  if (!NA_IsNArray(rblapack_d))
    rb_raise(rb_eArgError, "d (6th argument) must be NArray");
  if (NA_RANK(rblapack_d) != 1)
    rb_raise(rb_eArgError, "rank of d (6th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_d);
  if (NA_TYPE(rblapack_d) != NA_SFLOAT)
    rblapack_d = na_change_type(rblapack_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rblapack_d, real*);
  curlvl = NUM2INT(rblapack_curlvl);
  curpbm = NUM2INT(rblapack_curpbm);
  if (!NA_IsNArray(rblapack_q))
    rb_raise(rb_eArgError, "q (7th argument) must be NArray");
  if (NA_RANK(rblapack_q) != 2)
    rb_raise(rb_eArgError, "rank of q (7th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_q) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of q must be the same as shape 0 of d");
  ldq = NA_SHAPE0(rblapack_q);
  if (NA_TYPE(rblapack_q) != NA_SCOMPLEX)
    rblapack_q = na_change_type(rblapack_q, NA_SCOMPLEX);
  q = NA_PTR_TYPE(rblapack_q, complex*);
  if (!NA_IsNArray(rblapack_perm))
    rb_raise(rb_eArgError, "perm (12th argument) must be NArray");
  if (NA_RANK(rblapack_perm) != 1)
    rb_raise(rb_eArgError, "rank of perm (12th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_perm) != (n*LG(n)))
    rb_raise(rb_eRuntimeError, "shape 0 of perm must be %d", n*LG(n));
  if (NA_TYPE(rblapack_perm) != NA_LINT)
    rblapack_perm = na_change_type(rblapack_perm, NA_LINT);
  perm = NA_PTR_TYPE(rblapack_perm, integer*);
  if (!NA_IsNArray(rblapack_prmptr))
    rb_raise(rb_eArgError, "prmptr (11th argument) must be NArray");
  if (NA_RANK(rblapack_prmptr) != 1)
    rb_raise(rb_eArgError, "rank of prmptr (11th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_prmptr) != (n*LG(n)))
    rb_raise(rb_eRuntimeError, "shape 0 of prmptr must be %d", n*LG(n));
  if (NA_TYPE(rblapack_prmptr) != NA_LINT)
    rblapack_prmptr = na_change_type(rblapack_prmptr, NA_LINT);
  prmptr = NA_PTR_TYPE(rblapack_prmptr, integer*);
  if (!NA_IsNArray(rblapack_qstore))
    rb_raise(rb_eArgError, "qstore (9th argument) must be NArray");
  if (NA_RANK(rblapack_qstore) != 1)
    rb_raise(rb_eArgError, "rank of qstore (9th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_qstore) != (pow(n,2)+1))
    rb_raise(rb_eRuntimeError, "shape 0 of qstore must be %d", pow(n,2)+1);
  if (NA_TYPE(rblapack_qstore) != NA_SFLOAT)
    rblapack_qstore = na_change_type(rblapack_qstore, NA_SFLOAT);
  qstore = NA_PTR_TYPE(rblapack_qstore, real*);
  if (!NA_IsNArray(rblapack_givptr))
    rb_raise(rb_eArgError, "givptr (13th argument) must be NArray");
  if (NA_RANK(rblapack_givptr) != 1)
    rb_raise(rb_eArgError, "rank of givptr (13th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_givptr) != (n*LG(n)))
    rb_raise(rb_eRuntimeError, "shape 0 of givptr must be %d", n*LG(n));
  if (NA_TYPE(rblapack_givptr) != NA_LINT)
    rblapack_givptr = na_change_type(rblapack_givptr, NA_LINT);
  givptr = NA_PTR_TYPE(rblapack_givptr, integer*);
  if (!NA_IsNArray(rblapack_givcol))
    rb_raise(rb_eArgError, "givcol (14th argument) must be NArray");
  if (NA_RANK(rblapack_givcol) != 2)
    rb_raise(rb_eArgError, "rank of givcol (14th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_givcol) != (n*LG(n)))
    rb_raise(rb_eRuntimeError, "shape 1 of givcol must be %d", n*LG(n));
  if (NA_SHAPE0(rblapack_givcol) != (2))
    rb_raise(rb_eRuntimeError, "shape 0 of givcol must be %d", 2);
  if (NA_TYPE(rblapack_givcol) != NA_LINT)
    rblapack_givcol = na_change_type(rblapack_givcol, NA_LINT);
  givcol = NA_PTR_TYPE(rblapack_givcol, integer*);
  if (!NA_IsNArray(rblapack_givnum))
    rb_raise(rb_eArgError, "givnum (15th argument) must be NArray");
  if (NA_RANK(rblapack_givnum) != 2)
    rb_raise(rb_eArgError, "rank of givnum (15th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_givnum) != (n*LG(n)))
    rb_raise(rb_eRuntimeError, "shape 1 of givnum must be %d", n*LG(n));
  if (NA_SHAPE0(rblapack_givnum) != (2))
    rb_raise(rb_eRuntimeError, "shape 0 of givnum must be %d", 2);
  if (NA_TYPE(rblapack_givnum) != NA_SFLOAT)
    rblapack_givnum = na_change_type(rblapack_givnum, NA_SFLOAT);
  givnum = NA_PTR_TYPE(rblapack_givnum, real*);
  if (!NA_IsNArray(rblapack_qptr))
    rb_raise(rb_eArgError, "qptr (10th argument) must be NArray");
  if (NA_RANK(rblapack_qptr) != 1)
    rb_raise(rb_eArgError, "rank of qptr (10th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_qptr) != (n+2))
    rb_raise(rb_eRuntimeError, "shape 0 of qptr must be %d", n+2);
  if (NA_TYPE(rblapack_qptr) != NA_LINT)
    rblapack_qptr = na_change_type(rblapack_qptr, NA_LINT);
  qptr = NA_PTR_TYPE(rblapack_qptr, integer*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_indxq = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  indxq = NA_PTR_TYPE(rblapack_indxq, integer*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_d_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  d_out__ = NA_PTR_TYPE(rblapack_d_out__, real*);
  MEMCPY(d_out__, d, real, NA_TOTAL(rblapack_d));
  rblapack_d = rblapack_d_out__;
  d = d_out__;
  {
    int shape[2];
    shape[0] = ldq;
    shape[1] = n;
    rblapack_q_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  q_out__ = NA_PTR_TYPE(rblapack_q_out__, complex*);
  MEMCPY(q_out__, q, complex, NA_TOTAL(rblapack_q));
  rblapack_q = rblapack_q_out__;
  q = q_out__;
  {
    int shape[1];
    shape[0] = pow(n,2)+1;
    rblapack_qstore_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  qstore_out__ = NA_PTR_TYPE(rblapack_qstore_out__, real*);
  MEMCPY(qstore_out__, qstore, real, NA_TOTAL(rblapack_qstore));
  rblapack_qstore = rblapack_qstore_out__;
  qstore = qstore_out__;
  {
    int shape[1];
    shape[0] = n+2;
    rblapack_qptr_out__ = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  qptr_out__ = NA_PTR_TYPE(rblapack_qptr_out__, integer*);
  MEMCPY(qptr_out__, qptr, integer, NA_TOTAL(rblapack_qptr));
  rblapack_qptr = rblapack_qptr_out__;
  qptr = qptr_out__;
  work = ALLOC_N(complex, (qsiz*n));
  rwork = ALLOC_N(real, (3*n+2*qsiz*n));
  iwork = ALLOC_N(integer, (4*n));

  claed7_(&n, &cutpnt, &qsiz, &tlvls, &curlvl, &curpbm, d, q, &ldq, &rho, indxq, qstore, qptr, prmptr, perm, givptr, givcol, givnum, work, rwork, iwork, &info);

  free(work);
  free(rwork);
  free(iwork);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(6, rblapack_indxq, rblapack_info, rblapack_d, rblapack_q, rblapack_qstore, rblapack_qptr);
}

void
init_lapack_claed7(VALUE mLapack){
  rb_define_module_function(mLapack, "claed7", rblapack_claed7, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
