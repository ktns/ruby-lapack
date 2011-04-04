#include "rb_lapack.h"

extern VOID slaed1_(integer *n, real *d, real *q, integer *ldq, integer *indxq, real *rho, integer *cutpnt, real *work, integer *iwork, integer *info);

static VALUE
rb_slaed1(int argc, VALUE *argv, VALUE self){
  VALUE rb_d;
  real *d; 
  VALUE rb_q;
  real *q; 
  VALUE rb_indxq;
  integer *indxq; 
  VALUE rb_rho;
  real rho; 
  VALUE rb_cutpnt;
  integer cutpnt; 
  VALUE rb_info;
  integer info; 
  VALUE rb_d_out__;
  real *d_out__;
  VALUE rb_q_out__;
  real *q_out__;
  VALUE rb_indxq_out__;
  integer *indxq_out__;
  real *work;
  integer *iwork;

  integer n;
  integer ldq;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, d, q, indxq = NumRu::Lapack.slaed1( d, q, indxq, rho, cutpnt)\n    or\n  NumRu::Lapack.slaed1  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLAED1( N, D, Q, LDQ, INDXQ, RHO, CUTPNT, WORK, IWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  SLAED1 computes the updated eigensystem of a diagonal\n*  matrix after modification by a rank-one symmetric matrix.  This\n*  routine is used only for the eigenproblem which requires all\n*  eigenvalues and eigenvectors of a tridiagonal matrix.  SLAED7 handles\n*  the case in which eigenvalues only or eigenvalues and eigenvectors\n*  of a full symmetric matrix (which was reduced to tridiagonal form)\n*  are desired.\n*\n*    T = Q(in) ( D(in) + RHO * Z*Z' ) Q'(in) = Q(out) * D(out) * Q'(out)\n*\n*     where Z = Q'u, u is a vector of length N with ones in the\n*     CUTPNT and CUTPNT + 1 th elements and zeros elsewhere.\n*\n*     The eigenvectors of the original matrix are stored in Q, and the\n*     eigenvalues are in D.  The algorithm consists of three stages:\n*\n*        The first stage consists of deflating the size of the problem\n*        when there are multiple eigenvalues or if there is a zero in\n*        the Z vector.  For each such occurence the dimension of the\n*        secular equation problem is reduced by one.  This stage is\n*        performed by the routine SLAED2.\n*\n*        The second stage consists of calculating the updated\n*        eigenvalues. This is done by finding the roots of the secular\n*        equation via the routine SLAED4 (as called by SLAED3).\n*        This routine also calculates the eigenvectors of the current\n*        problem.\n*\n*        The final stage consists of computing the updated eigenvectors\n*        directly using the updated eigenvalues.  The eigenvectors for\n*        the current problem are multiplied with the eigenvectors from\n*        the overall problem.\n*\n\n*  Arguments\n*  =========\n*\n*  N      (input) INTEGER\n*         The dimension of the symmetric tridiagonal matrix.  N >= 0.\n*\n*  D      (input/output) REAL array, dimension (N)\n*         On entry, the eigenvalues of the rank-1-perturbed matrix.\n*         On exit, the eigenvalues of the repaired matrix.\n*\n*  Q      (input/output) REAL array, dimension (LDQ,N)\n*         On entry, the eigenvectors of the rank-1-perturbed matrix.\n*         On exit, the eigenvectors of the repaired tridiagonal matrix.\n*\n*  LDQ    (input) INTEGER\n*         The leading dimension of the array Q.  LDQ >= max(1,N).\n*\n*  INDXQ  (input/output) INTEGER array, dimension (N)\n*         On entry, the permutation which separately sorts the two\n*         subproblems in D into ascending order.\n*         On exit, the permutation which will reintegrate the\n*         subproblems back into sorted order,\n*         i.e. D( INDXQ( I = 1, N ) ) will be in ascending order.\n*\n*  RHO    (input) REAL\n*         The subdiagonal entry used to create the rank-1 modification.\n*\n*  CUTPNT (input) INTEGER\n*         The location of the last eigenvalue in the leading sub-matrix.\n*         min(1,N) <= CUTPNT <= N/2.\n*\n*  WORK   (workspace) REAL array, dimension (4*N + N**2)\n*\n*  IWORK  (workspace) INTEGER array, dimension (4*N)\n*\n*  INFO   (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*          > 0:  if INFO = 1, an eigenvalue did not converge\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Jeff Rutter, Computer Science Division, University of California\n*     at Berkeley, USA\n*  Modified by Francoise Tisseur, University of Tennessee.\n*\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      INTEGER            COLTYP, CPP1, I, IDLMDA, INDX, INDXC, INDXP,\n     $                   IQ2, IS, IW, IZ, K, N1, N2\n*     ..\n*     .. External Subroutines ..\n      EXTERNAL           SCOPY, SLAED2, SLAED3, SLAMRG, XERBLA\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          MAX, MIN\n*     ..\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_d = argv[0];
  rb_q = argv[1];
  rb_indxq = argv[2];
  rb_rho = argv[3];
  rb_cutpnt = argv[4];

  cutpnt = NUM2INT(rb_cutpnt);
  if (!NA_IsNArray(rb_indxq))
    rb_raise(rb_eArgError, "indxq (3th argument) must be NArray");
  if (NA_RANK(rb_indxq) != 1)
    rb_raise(rb_eArgError, "rank of indxq (3th argument) must be %d", 1);
  n = NA_SHAPE0(rb_indxq);
  if (NA_TYPE(rb_indxq) != NA_LINT)
    rb_indxq = na_change_type(rb_indxq, NA_LINT);
  indxq = NA_PTR_TYPE(rb_indxq, integer*);
  rho = (real)NUM2DBL(rb_rho);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (1th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (1th argument) must be %d", 1);
  if (NA_SHAPE0(rb_d) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of d must be the same as shape 0 of indxq");
  if (NA_TYPE(rb_d) != NA_SFLOAT)
    rb_d = na_change_type(rb_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rb_d, real*);
  if (!NA_IsNArray(rb_q))
    rb_raise(rb_eArgError, "q (2th argument) must be NArray");
  if (NA_RANK(rb_q) != 2)
    rb_raise(rb_eArgError, "rank of q (2th argument) must be %d", 2);
  if (NA_SHAPE1(rb_q) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of q must be the same as shape 0 of indxq");
  ldq = NA_SHAPE0(rb_q);
  if (NA_TYPE(rb_q) != NA_SFLOAT)
    rb_q = na_change_type(rb_q, NA_SFLOAT);
  q = NA_PTR_TYPE(rb_q, real*);
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
    shape[0] = n;
    rb_indxq_out__ = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  indxq_out__ = NA_PTR_TYPE(rb_indxq_out__, integer*);
  MEMCPY(indxq_out__, indxq, integer, NA_TOTAL(rb_indxq));
  rb_indxq = rb_indxq_out__;
  indxq = indxq_out__;
  work = ALLOC_N(real, (4*n + pow(n,2)));
  iwork = ALLOC_N(integer, (4*n));

  slaed1_(&n, d, q, &ldq, indxq, &rho, &cutpnt, work, iwork, &info);

  free(work);
  free(iwork);
  rb_info = INT2NUM(info);
  return rb_ary_new3(4, rb_info, rb_d, rb_q, rb_indxq);
}

void
init_lapack_slaed1(VALUE mLapack){
  rb_define_module_function(mLapack, "slaed1", rb_slaed1, -1);
}
