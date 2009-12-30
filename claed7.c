#include "rb_lapack.h"

static VALUE
rb_claed7(int argc, VALUE *argv, VALUE self){
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
  real *d; 
  VALUE rb_q;
  complex *q; 
  VALUE rb_rho;
  real rho; 
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
  complex *q_out__;
  VALUE rb_qstore_out__;
  real *qstore_out__;
  VALUE rb_qptr_out__;
  integer *qptr_out__;
  complex *work;
  real *rwork;
  integer *iwork;

  integer n;
  integer ldq;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  indxq, info, d, q, qstore, qptr = NumRu::Lapack.claed7( cutpnt, qsiz, tlvls, curlvl, curpbm, d, q, rho, qstore, qptr, prmptr, perm, givptr, givcol, givnum)\n    or\n  NumRu::Lapack.claed7  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE CLAED7( N, CUTPNT, QSIZ, TLVLS, CURLVL, CURPBM, D, Q, LDQ, RHO, INDXQ, QSTORE, QPTR, PRMPTR, PERM, GIVPTR, GIVCOL, GIVNUM, WORK, RWORK, IWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  CLAED7 computes the updated eigensystem of a diagonal\n*  matrix after modification by a rank-one symmetric matrix. This\n*  routine is used only for the eigenproblem which requires all\n*  eigenvalues and optionally eigenvectors of a dense or banded\n*  Hermitian matrix that has been reduced to tridiagonal form.\n*\n*    T = Q(in) ( D(in) + RHO * Z*Z' ) Q'(in) = Q(out) * D(out) * Q'(out)\n*\n*    where Z = Q'u, u is a vector of length N with ones in the\n*    CUTPNT and CUTPNT + 1 th elements and zeros elsewhere.\n*\n*     The eigenvectors of the original matrix are stored in Q, and the\n*     eigenvalues are in D.  The algorithm consists of three stages:\n*\n*        The first stage consists of deflating the size of the problem\n*        when there are multiple eigenvalues or if there is a zero in\n*        the Z vector.  For each such occurence the dimension of the\n*        secular equation problem is reduced by one.  This stage is\n*        performed by the routine SLAED2.\n*\n*        The second stage consists of calculating the updated\n*        eigenvalues. This is done by finding the roots of the secular\n*        equation via the routine SLAED4 (as called by SLAED3).\n*        This routine also calculates the eigenvectors of the current\n*        problem.\n*\n*        The final stage consists of computing the updated eigenvectors\n*        directly using the updated eigenvalues.  The eigenvectors for\n*        the current problem are multiplied with the eigenvectors from\n*        the overall problem.\n*\n\n*  Arguments\n*  =========\n*\n*  N      (input) INTEGER\n*         The dimension of the symmetric tridiagonal matrix.  N >= 0.\n*\n*  CUTPNT (input) INTEGER\n*         Contains the location of the last eigenvalue in the leading\n*         sub-matrix.  min(1,N) <= CUTPNT <= N.\n*\n*  QSIZ   (input) INTEGER\n*         The dimension of the unitary matrix used to reduce\n*         the full matrix to tridiagonal form.  QSIZ >= N.\n*\n*  TLVLS  (input) INTEGER\n*         The total number of merging levels in the overall divide and\n*         conquer tree.\n*\n*  CURLVL (input) INTEGER\n*         The current level in the overall merge routine,\n*         0 <= curlvl <= tlvls.\n*\n*  CURPBM (input) INTEGER\n*         The current problem in the current level in the overall\n*         merge routine (counting from upper left to lower right).\n*\n*  D      (input/output) REAL array, dimension (N)\n*         On entry, the eigenvalues of the rank-1-perturbed matrix.\n*         On exit, the eigenvalues of the repaired matrix.\n*\n*  Q      (input/output) COMPLEX array, dimension (LDQ,N)\n*         On entry, the eigenvectors of the rank-1-perturbed matrix.\n*         On exit, the eigenvectors of the repaired tridiagonal matrix.\n*\n*  LDQ    (input) INTEGER\n*         The leading dimension of the array Q.  LDQ >= max(1,N).\n*\n*  RHO    (input) REAL\n*         Contains the subdiagonal element used to create the rank-1\n*         modification.\n*\n*  INDXQ  (output) INTEGER array, dimension (N)\n*         This contains the permutation which will reintegrate the\n*         subproblem just solved back into sorted order,\n*         ie. D( INDXQ( I = 1, N ) ) will be in ascending order.\n*\n*  IWORK  (workspace) INTEGER array, dimension (4*N)\n*\n*  RWORK  (workspace) REAL array,\n*                                 dimension (3*N+2*QSIZ*N)\n*\n*  WORK   (workspace) COMPLEX array, dimension (QSIZ*N)\n*\n*  QSTORE (input/output) REAL array, dimension (N**2+1)\n*         Stores eigenvectors of submatrices encountered during\n*         divide and conquer, packed together. QPTR points to\n*         beginning of the submatrices.\n*\n*  QPTR   (input/output) INTEGER array, dimension (N+2)\n*         List of indices pointing to beginning of submatrices stored\n*         in QSTORE. The submatrices are numbered starting at the\n*         bottom left of the divide and conquer tree, from left to\n*         right and bottom to top.\n*\n*  PRMPTR (input) INTEGER array, dimension (N lg N)\n*         Contains a list of pointers which indicate where in PERM a\n*         level's permutation is stored.  PRMPTR(i+1) - PRMPTR(i)\n*         indicates the size of the permutation and also the size of\n*         the full, non-deflated problem.\n*\n*  PERM   (input) INTEGER array, dimension (N lg N)\n*         Contains the permutations (from deflation and sorting) to be\n*         applied to each eigenblock.\n*\n*  GIVPTR (input) INTEGER array, dimension (N lg N)\n*         Contains a list of pointers which indicate where in GIVCOL a\n*         level's Givens rotations are stored.  GIVPTR(i+1) - GIVPTR(i)\n*         indicates the number of Givens rotations.\n*\n*  GIVCOL (input) INTEGER array, dimension (2, N lg N)\n*         Each pair of numbers indicates a pair of columns to take place\n*         in a Givens rotation.\n*\n*  GIVNUM (input) REAL array, dimension (2, N lg N)\n*         Each number indicates the S value to be used in the\n*         corresponding Givens rotation.\n*\n*  INFO   (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*          > 0:  if INFO = 1, an eigenvalue did not converge\n*\n\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      INTEGER            COLTYP, CURR, I, IDLMDA, INDX,\n     $                   INDXC, INDXP, IQ, IW, IZ, K, N1, N2, PTR\n*     ..\n*     .. External Subroutines ..\n      EXTERNAL           CLACRM, CLAED8, SLAED9, SLAEDA, SLAMRG, XERBLA\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          MAX, MIN\n*     ..\n\n");
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

  cutpnt = NUM2INT(rb_cutpnt);
  qsiz = NUM2INT(rb_qsiz);
  tlvls = NUM2INT(rb_tlvls);
  curlvl = NUM2INT(rb_curlvl);
  curpbm = NUM2INT(rb_curpbm);
  rho = (real)NUM2DBL(rb_rho);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (6th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (6th argument) must be %d", 1);
  n = NA_SHAPE0(rb_d);
  if (NA_TYPE(rb_d) != NA_SFLOAT)
    rb_d = na_change_type(rb_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rb_d, real*);
  if (!NA_IsNArray(rb_q))
    rb_raise(rb_eArgError, "q (7th argument) must be NArray");
  if (NA_RANK(rb_q) != 2)
    rb_raise(rb_eArgError, "rank of q (7th argument) must be %d", 2);
  ldq = NA_SHAPE0(rb_q);
  if (NA_SHAPE1(rb_q) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of q must be the same as shape 0 of d");
  if (NA_TYPE(rb_q) != NA_SCOMPLEX)
    rb_q = na_change_type(rb_q, NA_SCOMPLEX);
  q = NA_PTR_TYPE(rb_q, complex*);
  if (!NA_IsNArray(rb_qstore))
    rb_raise(rb_eArgError, "qstore (9th argument) must be NArray");
  if (NA_RANK(rb_qstore) != 1)
    rb_raise(rb_eArgError, "rank of qstore (9th argument) must be %d", 1);
  if (NA_SHAPE0(rb_qstore) != (pow(n,2)+1))
    rb_raise(rb_eRuntimeError, "shape 0 of qstore must be %d", pow(n,2)+1);
  if (NA_TYPE(rb_qstore) != NA_SFLOAT)
    rb_qstore = na_change_type(rb_qstore, NA_SFLOAT);
  qstore = NA_PTR_TYPE(rb_qstore, real*);
  if (!NA_IsNArray(rb_qptr))
    rb_raise(rb_eArgError, "qptr (10th argument) must be NArray");
  if (NA_RANK(rb_qptr) != 1)
    rb_raise(rb_eArgError, "rank of qptr (10th argument) must be %d", 1);
  if (NA_SHAPE0(rb_qptr) != (n+2))
    rb_raise(rb_eRuntimeError, "shape 0 of qptr must be %d", n+2);
  if (NA_TYPE(rb_qptr) != NA_LINT)
    rb_qptr = na_change_type(rb_qptr, NA_LINT);
  qptr = NA_PTR_TYPE(rb_qptr, integer*);
  if (!NA_IsNArray(rb_prmptr))
    rb_raise(rb_eArgError, "prmptr (11th argument) must be NArray");
  if (NA_RANK(rb_prmptr) != 1)
    rb_raise(rb_eArgError, "rank of prmptr (11th argument) must be %d", 1);
  if (NA_SHAPE0(rb_prmptr) != (n*LG(n)))
    rb_raise(rb_eRuntimeError, "shape 0 of prmptr must be %d", n*LG(n));
  if (NA_TYPE(rb_prmptr) != NA_LINT)
    rb_prmptr = na_change_type(rb_prmptr, NA_LINT);
  prmptr = NA_PTR_TYPE(rb_prmptr, integer*);
  if (!NA_IsNArray(rb_perm))
    rb_raise(rb_eArgError, "perm (12th argument) must be NArray");
  if (NA_RANK(rb_perm) != 1)
    rb_raise(rb_eArgError, "rank of perm (12th argument) must be %d", 1);
  if (NA_SHAPE0(rb_perm) != (n*LG(n)))
    rb_raise(rb_eRuntimeError, "shape 0 of perm must be %d", n*LG(n));
  if (NA_TYPE(rb_perm) != NA_LINT)
    rb_perm = na_change_type(rb_perm, NA_LINT);
  perm = NA_PTR_TYPE(rb_perm, integer*);
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
  if (NA_SHAPE0(rb_givcol) != (2))
    rb_raise(rb_eRuntimeError, "shape 0 of givcol must be %d", 2);
  if (NA_SHAPE1(rb_givcol) != (n*LG(n)))
    rb_raise(rb_eRuntimeError, "shape 1 of givcol must be %d", n*LG(n));
  if (NA_TYPE(rb_givcol) != NA_LINT)
    rb_givcol = na_change_type(rb_givcol, NA_LINT);
  givcol = NA_PTR_TYPE(rb_givcol, integer*);
  if (!NA_IsNArray(rb_givnum))
    rb_raise(rb_eArgError, "givnum (15th argument) must be NArray");
  if (NA_RANK(rb_givnum) != 2)
    rb_raise(rb_eArgError, "rank of givnum (15th argument) must be %d", 2);
  if (NA_SHAPE0(rb_givnum) != (2))
    rb_raise(rb_eRuntimeError, "shape 0 of givnum must be %d", 2);
  if (NA_SHAPE1(rb_givnum) != (n*LG(n)))
    rb_raise(rb_eRuntimeError, "shape 1 of givnum must be %d", n*LG(n));
  if (NA_TYPE(rb_givnum) != NA_SFLOAT)
    rb_givnum = na_change_type(rb_givnum, NA_SFLOAT);
  givnum = NA_PTR_TYPE(rb_givnum, real*);
  {
    int shape[1];
    shape[0] = DIM_LEN(n);
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
    rb_q_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  q_out__ = NA_PTR_TYPE(rb_q_out__, complex*);
  MEMCPY(q_out__, q, complex, NA_TOTAL(rb_q));
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
    shape[0] = DIM_LEN(n+2);
    rb_qptr_out__ = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  qptr_out__ = NA_PTR_TYPE(rb_qptr_out__, integer*);
  MEMCPY(qptr_out__, qptr, integer, NA_TOTAL(rb_qptr));
  rb_qptr = rb_qptr_out__;
  qptr = qptr_out__;
  work = ALLOC_N(complex, (qsiz*n));
  rwork = ALLOC_N(real, (3*n+2*qsiz*n));
  iwork = ALLOC_N(integer, (4*n));

  claed7_(&n, &cutpnt, &qsiz, &tlvls, &curlvl, &curpbm, d, q, &ldq, &rho, indxq, qstore, qptr, prmptr, perm, givptr, givcol, givnum, work, rwork, iwork, &info);

  free(work);
  free(rwork);
  free(iwork);
  rb_info = INT2NUM(info);
  return rb_ary_new3(6, rb_indxq, rb_info, rb_d, rb_q, rb_qstore, rb_qptr);
}

void
init_lapack_claed7(VALUE mLapack){
  rb_define_module_function(mLapack, "claed7", rb_claed7, -1);
}
