#include "rb_lapack.h"

extern VOID dgeevx_(char *balanc, char *jobvl, char *jobvr, char *sense, integer *n, doublereal *a, integer *lda, doublereal *wr, doublereal *wi, doublereal *vl, integer *ldvl, doublereal *vr, integer *ldvr, integer *ilo, integer *ihi, doublereal *scale, doublereal *abnrm, doublereal *rconde, doublereal *rcondv, doublereal *work, integer *lwork, integer *iwork, integer *info);

static VALUE
rb_dgeevx(int argc, VALUE *argv, VALUE self){
  VALUE rb_balanc;
  char balanc; 
  VALUE rb_jobvl;
  char jobvl; 
  VALUE rb_jobvr;
  char jobvr; 
  VALUE rb_sense;
  char sense; 
  VALUE rb_a;
  doublereal *a; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_wr;
  doublereal *wr; 
  VALUE rb_wi;
  doublereal *wi; 
  VALUE rb_vl;
  doublereal *vl; 
  VALUE rb_vr;
  doublereal *vr; 
  VALUE rb_ilo;
  integer ilo; 
  VALUE rb_ihi;
  integer ihi; 
  VALUE rb_scale;
  doublereal *scale; 
  VALUE rb_abnrm;
  doublereal abnrm; 
  VALUE rb_rconde;
  doublereal *rconde; 
  VALUE rb_rcondv;
  doublereal *rcondv; 
  VALUE rb_work;
  doublereal *work; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  doublereal *a_out__;
  integer *iwork;

  integer lda;
  integer n;
  integer ldvl;
  integer ldvr;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  wr, wi, vl, vr, ilo, ihi, scale, abnrm, rconde, rcondv, work, info, a = NumRu::Lapack.dgeevx( balanc, jobvl, jobvr, sense, a, lwork)\n    or\n  NumRu::Lapack.dgeevx  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DGEEVX( BALANC, JOBVL, JOBVR, SENSE, N, A, LDA, WR, WI, VL, LDVL, VR, LDVR, ILO, IHI, SCALE, ABNRM, RCONDE, RCONDV, WORK, LWORK, IWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  DGEEVX computes for an N-by-N real nonsymmetric matrix A, the\n*  eigenvalues and, optionally, the left and/or right eigenvectors.\n*\n*  Optionally also, it computes a balancing transformation to improve\n*  the conditioning of the eigenvalues and eigenvectors (ILO, IHI,\n*  SCALE, and ABNRM), reciprocal condition numbers for the eigenvalues\n*  (RCONDE), and reciprocal condition numbers for the right\n*  eigenvectors (RCONDV).\n*\n*  The right eigenvector v(j) of A satisfies\n*                   A * v(j) = lambda(j) * v(j)\n*  where lambda(j) is its eigenvalue.\n*  The left eigenvector u(j) of A satisfies\n*                u(j)**H * A = lambda(j) * u(j)**H\n*  where u(j)**H denotes the conjugate transpose of u(j).\n*\n*  The computed eigenvectors are normalized to have Euclidean norm\n*  equal to 1 and largest component real.\n*\n*  Balancing a matrix means permuting the rows and columns to make it\n*  more nearly upper triangular, and applying a diagonal similarity\n*  transformation D * A * D**(-1), where D is a diagonal matrix, to\n*  make its rows and columns closer in norm and the condition numbers\n*  of its eigenvalues and eigenvectors smaller.  The computed\n*  reciprocal condition numbers correspond to the balanced matrix.\n*  Permuting rows and columns will not change the condition numbers\n*  (in exact arithmetic) but diagonal scaling will.  For further\n*  explanation of balancing, see section 4.10.2 of the LAPACK\n*  Users' Guide.\n*\n\n*  Arguments\n*  =========\n*\n*  BALANC  (input) CHARACTER*1\n*          Indicates how the input matrix should be diagonally scaled\n*          and/or permuted to improve the conditioning of its\n*          eigenvalues.\n*          = 'N': Do not diagonally scale or permute;\n*          = 'P': Perform permutations to make the matrix more nearly\n*                 upper triangular. Do not diagonally scale;\n*          = 'S': Diagonally scale the matrix, i.e. replace A by\n*                 D*A*D**(-1), where D is a diagonal matrix chosen\n*                 to make the rows and columns of A more equal in\n*                 norm. Do not permute;\n*          = 'B': Both diagonally scale and permute A.\n*\n*          Computed reciprocal condition numbers will be for the matrix\n*          after balancing and/or permuting. Permuting does not change\n*          condition numbers (in exact arithmetic), but balancing does.\n*\n*  JOBVL   (input) CHARACTER*1\n*          = 'N': left eigenvectors of A are not computed;\n*          = 'V': left eigenvectors of A are computed.\n*          If SENSE = 'E' or 'B', JOBVL must = 'V'.\n*\n*  JOBVR   (input) CHARACTER*1\n*          = 'N': right eigenvectors of A are not computed;\n*          = 'V': right eigenvectors of A are computed.\n*          If SENSE = 'E' or 'B', JOBVR must = 'V'.\n*\n*  SENSE   (input) CHARACTER*1\n*          Determines which reciprocal condition numbers are computed.\n*          = 'N': None are computed;\n*          = 'E': Computed for eigenvalues only;\n*          = 'V': Computed for right eigenvectors only;\n*          = 'B': Computed for eigenvalues and right eigenvectors.\n*\n*          If SENSE = 'E' or 'B', both left and right eigenvectors\n*          must also be computed (JOBVL = 'V' and JOBVR = 'V').\n*\n*  N       (input) INTEGER\n*          The order of the matrix A. N >= 0.\n*\n*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)\n*          On entry, the N-by-N matrix A.\n*          On exit, A has been overwritten.  If JOBVL = 'V' or\n*          JOBVR = 'V', A contains the real Schur form of the balanced\n*          version of the input matrix A.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,N).\n*\n*  WR      (output) DOUBLE PRECISION array, dimension (N)\n*  WI      (output) DOUBLE PRECISION array, dimension (N)\n*          WR and WI contain the real and imaginary parts,\n*          respectively, of the computed eigenvalues.  Complex\n*          conjugate pairs of eigenvalues will appear consecutively\n*          with the eigenvalue having the positive imaginary part\n*          first.\n*\n*  VL      (output) DOUBLE PRECISION array, dimension (LDVL,N)\n*          If JOBVL = 'V', the left eigenvectors u(j) are stored one\n*          after another in the columns of VL, in the same order\n*          as their eigenvalues.\n*          If JOBVL = 'N', VL is not referenced.\n*          If the j-th eigenvalue is real, then u(j) = VL(:,j),\n*          the j-th column of VL.\n*          If the j-th and (j+1)-st eigenvalues form a complex\n*          conjugate pair, then u(j) = VL(:,j) + i*VL(:,j+1) and\n*          u(j+1) = VL(:,j) - i*VL(:,j+1).\n*\n*  LDVL    (input) INTEGER\n*          The leading dimension of the array VL.  LDVL >= 1; if\n*          JOBVL = 'V', LDVL >= N.\n*\n*  VR      (output) DOUBLE PRECISION array, dimension (LDVR,N)\n*          If JOBVR = 'V', the right eigenvectors v(j) are stored one\n*          after another in the columns of VR, in the same order\n*          as their eigenvalues.\n*          If JOBVR = 'N', VR is not referenced.\n*          If the j-th eigenvalue is real, then v(j) = VR(:,j),\n*          the j-th column of VR.\n*          If the j-th and (j+1)-st eigenvalues form a complex\n*          conjugate pair, then v(j) = VR(:,j) + i*VR(:,j+1) and\n*          v(j+1) = VR(:,j) - i*VR(:,j+1).\n*\n*  LDVR    (input) INTEGER\n*          The leading dimension of the array VR.  LDVR >= 1, and if\n*          JOBVR = 'V', LDVR >= N.\n*\n*  ILO     (output) INTEGER\n*  IHI     (output) INTEGER\n*          ILO and IHI are integer values determined when A was\n*          balanced.  The balanced A(i,j) = 0 if I > J and\n*          J = 1,...,ILO-1 or I = IHI+1,...,N.\n*\n*  SCALE   (output) DOUBLE PRECISION array, dimension (N)\n*          Details of the permutations and scaling factors applied\n*          when balancing A.  If P(j) is the index of the row and column\n*          interchanged with row and column j, and D(j) is the scaling\n*          factor applied to row and column j, then\n*          SCALE(J) = P(J),    for J = 1,...,ILO-1\n*                   = D(J),    for J = ILO,...,IHI\n*                   = P(J)     for J = IHI+1,...,N.\n*          The order in which the interchanges are made is N to IHI+1,\n*          then 1 to ILO-1.\n*\n*  ABNRM   (output) DOUBLE PRECISION\n*          The one-norm of the balanced matrix (the maximum\n*          of the sum of absolute values of elements of any column).\n*\n*  RCONDE  (output) DOUBLE PRECISION array, dimension (N)\n*          RCONDE(j) is the reciprocal condition number of the j-th\n*          eigenvalue.\n*\n*  RCONDV  (output) DOUBLE PRECISION array, dimension (N)\n*          RCONDV(j) is the reciprocal condition number of the j-th\n*          right eigenvector.\n*\n*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK.   If SENSE = 'N' or 'E',\n*          LWORK >= max(1,2*N), and if JOBVL = 'V' or JOBVR = 'V',\n*          LWORK >= 3*N.  If SENSE = 'V' or 'B', LWORK >= N*(N+6).\n*          For good performance, LWORK must generally be larger.\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  IWORK   (workspace) INTEGER array, dimension (2*N-2)\n*          If SENSE = 'N' or 'E', not referenced.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*          > 0:  if INFO = i, the QR algorithm failed to compute all the\n*                eigenvalues, and no eigenvectors or condition numbers\n*                have been computed; elements 1:ILO-1 and i+1:N of WR\n*                and WI contain eigenvalues which have converged.\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rb_balanc = argv[0];
  rb_jobvl = argv[1];
  rb_jobvr = argv[2];
  rb_sense = argv[3];
  rb_a = argv[4];
  rb_lwork = argv[5];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (5th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (5th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_DFLOAT)
    rb_a = na_change_type(rb_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rb_a, doublereal*);
  jobvr = StringValueCStr(rb_jobvr)[0];
  balanc = StringValueCStr(rb_balanc)[0];
  sense = StringValueCStr(rb_sense)[0];
  jobvl = StringValueCStr(rb_jobvl)[0];
  lwork = NUM2INT(rb_lwork);
  ldvl = lsame_(&jobvl,"V") ? n : 1;
  ldvr = lsame_(&jobvr,"V") ? n : 1;
  {
    int shape[1];
    shape[0] = n;
    rb_wr = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  wr = NA_PTR_TYPE(rb_wr, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_wi = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  wi = NA_PTR_TYPE(rb_wi, doublereal*);
  {
    int shape[2];
    shape[0] = ldvl;
    shape[1] = n;
    rb_vl = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  vl = NA_PTR_TYPE(rb_vl, doublereal*);
  {
    int shape[2];
    shape[0] = ldvr;
    shape[1] = n;
    rb_vr = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  vr = NA_PTR_TYPE(rb_vr, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_scale = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  scale = NA_PTR_TYPE(rb_scale, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_rconde = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  rconde = NA_PTR_TYPE(rb_rconde, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_rcondv = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  rcondv = NA_PTR_TYPE(rb_rcondv, doublereal*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rb_work = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rb_work, doublereal*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rb_a_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rb_a_out__, doublereal*);
  MEMCPY(a_out__, a, doublereal, NA_TOTAL(rb_a));
  rb_a = rb_a_out__;
  a = a_out__;
  iwork = ALLOC_N(integer, ((lsame_(&sense,"N")||lsame_(&sense,"E")) ? 0 : 2*n-2));

  dgeevx_(&balanc, &jobvl, &jobvr, &sense, &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr, &ilo, &ihi, scale, &abnrm, rconde, rcondv, work, &lwork, iwork, &info);

  free(iwork);
  rb_ilo = INT2NUM(ilo);
  rb_ihi = INT2NUM(ihi);
  rb_abnrm = rb_float_new((double)abnrm);
  rb_info = INT2NUM(info);
  return rb_ary_new3(13, rb_wr, rb_wi, rb_vl, rb_vr, rb_ilo, rb_ihi, rb_scale, rb_abnrm, rb_rconde, rb_rcondv, rb_work, rb_info, rb_a);
}

void
init_lapack_dgeevx(VALUE mLapack){
  rb_define_module_function(mLapack, "dgeevx", rb_dgeevx, -1);
}
