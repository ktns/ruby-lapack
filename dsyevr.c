#include "rb_lapack.h"

static VALUE
rb_dsyevr(int argc, VALUE *argv, VALUE self){
  VALUE rb_jobz;
  char jobz; 
  VALUE rb_range;
  char range; 
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_a;
  doublereal *a; 
  VALUE rb_vl;
  doublereal vl; 
  VALUE rb_vu;
  doublereal vu; 
  VALUE rb_il;
  integer il; 
  VALUE rb_iu;
  integer iu; 
  VALUE rb_abstol;
  doublereal abstol; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_liwork;
  integer liwork; 
  VALUE rb_m;
  integer m; 
  VALUE rb_w;
  doublereal *w; 
  VALUE rb_z;
  doublereal *z; 
  VALUE rb_isuppz;
  integer *isuppz; 
  VALUE rb_work;
  doublereal *work; 
  VALUE rb_iwork;
  integer *iwork; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  doublereal *a_out__;

  integer lda;
  integer n;
  integer ldz;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  m, w, z, isuppz, work, iwork, info, a = NumRu::Lapack.dsyevr( jobz, range, uplo, a, vl, vu, il, iu, abstol, lwork, liwork)\n    or\n  NumRu::Lapack.dsyevr  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DSYEVR( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK, IWORK, LIWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  DSYEVR computes selected eigenvalues and, optionally, eigenvectors\n*  of a real symmetric matrix A.  Eigenvalues and eigenvectors can be\n*  selected by specifying either a range of values or a range of\n*  indices for the desired eigenvalues.\n*\n*  DSYEVR first reduces the matrix A to tridiagonal form T with a call\n*  to DSYTRD.  Then, whenever possible, DSYEVR calls DSTEMR to compute\n*  the eigenspectrum using Relatively Robust Representations.  DSTEMR\n*  computes eigenvalues by the dqds algorithm, while orthogonal\n*  eigenvectors are computed from various \"good\" L D L^T representations\n*  (also known as Relatively Robust Representations). Gram-Schmidt\n*  orthogonalization is avoided as far as possible. More specifically,\n*  the various steps of the algorithm are as follows.\n*\n*  For each unreduced block (submatrix) of T,\n*     (a) Compute T - sigma I  = L D L^T, so that L and D\n*         define all the wanted eigenvalues to high relative accuracy.\n*         This means that small relative changes in the entries of D and L\n*         cause only small relative changes in the eigenvalues and\n*         eigenvectors. The standard (unfactored) representation of the\n*         tridiagonal matrix T does not have this property in general.\n*     (b) Compute the eigenvalues to suitable accuracy.\n*         If the eigenvectors are desired, the algorithm attains full\n*         accuracy of the computed eigenvalues only right before\n*         the corresponding vectors have to be computed, see steps c) and d).\n*     (c) For each cluster of close eigenvalues, select a new\n*         shift close to the cluster, find a new factorization, and refine\n*         the shifted eigenvalues to suitable accuracy.\n*     (d) For each eigenvalue with a large enough relative separation compute\n*         the corresponding eigenvector by forming a rank revealing twisted\n*         factorization. Go back to (c) for any clusters that remain.\n*\n*  The desired accuracy of the output can be specified by the input\n*  parameter ABSTOL.\n*\n*  For more details, see DSTEMR's documentation and:\n*  - Inderjit S. Dhillon and Beresford N. Parlett: \"Multiple representations\n*    to compute orthogonal eigenvectors of symmetric tridiagonal matrices,\"\n*    Linear Algebra and its Applications, 387(1), pp. 1-28, August 2004.\n*  - Inderjit Dhillon and Beresford Parlett: \"Orthogonal Eigenvectors and\n*    Relative Gaps,\" SIAM Journal on Matrix Analysis and Applications, Vol. 25,\n*    2004.  Also LAPACK Working Note 154.\n*  - Inderjit Dhillon: \"A new O(n^2) algorithm for the symmetric\n*    tridiagonal eigenvalue/eigenvector problem\",\n*    Computer Science Division Technical Report No. UCB/CSD-97-971,\n*    UC Berkeley, May 1997.\n*\n*\n*  Note 1 : DSYEVR calls DSTEMR when the full spectrum is requested\n*  on machines which conform to the ieee-754 floating point standard.\n*  DSYEVR calls DSTEBZ and SSTEIN on non-ieee machines and\n*  when partial spectrum requests are made.\n*\n*  Normal execution of DSTEMR may create NaNs and infinities and\n*  hence may abort due to a floating point exception in environments\n*  which do not handle NaNs and infinities in the ieee standard default\n*  manner.\n*\n\n*  Arguments\n*  =========\n*\n*  JOBZ    (input) CHARACTER*1\n*          = 'N':  Compute eigenvalues only;\n*          = 'V':  Compute eigenvalues and eigenvectors.\n*\n*  RANGE   (input) CHARACTER*1\n*          = 'A': all eigenvalues will be found.\n*          = 'V': all eigenvalues in the half-open interval (VL,VU]\n*                 will be found.\n*          = 'I': the IL-th through IU-th eigenvalues will be found.\n********** For RANGE = 'V' or 'I' and IU - IL < N - 1, DSTEBZ and\n********** DSTEIN are called\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U':  Upper triangle of A is stored;\n*          = 'L':  Lower triangle of A is stored.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)\n*          On entry, the symmetric matrix A.  If UPLO = 'U', the\n*          leading N-by-N upper triangular part of A contains the\n*          upper triangular part of the matrix A.  If UPLO = 'L',\n*          the leading N-by-N lower triangular part of A contains\n*          the lower triangular part of the matrix A.\n*          On exit, the lower triangle (if UPLO='L') or the upper\n*          triangle (if UPLO='U') of A, including the diagonal, is\n*          destroyed.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,N).\n*\n*  VL      (input) DOUBLE PRECISION\n*  VU      (input) DOUBLE PRECISION\n*          If RANGE='V', the lower and upper bounds of the interval to\n*          be searched for eigenvalues. VL < VU.\n*          Not referenced if RANGE = 'A' or 'I'.\n*\n*  IL      (input) INTEGER\n*  IU      (input) INTEGER\n*          If RANGE='I', the indices (in ascending order) of the\n*          smallest and largest eigenvalues to be returned.\n*          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.\n*          Not referenced if RANGE = 'A' or 'V'.\n*\n*  ABSTOL  (input) DOUBLE PRECISION\n*          The absolute error tolerance for the eigenvalues.\n*          An approximate eigenvalue is accepted as converged\n*          when it is determined to lie in an interval [a,b]\n*          of width less than or equal to\n*\n*                  ABSTOL + EPS *   max( |a|,|b| ) ,\n*\n*          where EPS is the machine precision.  If ABSTOL is less than\n*          or equal to zero, then  EPS*|T|  will be used in its place,\n*          where |T| is the 1-norm of the tridiagonal matrix obtained\n*          by reducing A to tridiagonal form.\n*\n*          See \"Computing Small Singular Values of Bidiagonal Matrices\n*          with Guaranteed High Relative Accuracy,\" by Demmel and\n*          Kahan, LAPACK Working Note #3.\n*\n*          If high relative accuracy is important, set ABSTOL to\n*          DLAMCH( 'Safe minimum' ).  Doing so will guarantee that\n*          eigenvalues are computed to high relative accuracy when\n*          possible in future releases.  The current code does not\n*          make any guarantees about high relative accuracy, but\n*          future releases will. See J. Barlow and J. Demmel,\n*          \"Computing Accurate Eigensystems of Scaled Diagonally\n*          Dominant Matrices\", LAPACK Working Note #7, for a discussion\n*          of which matrices define their eigenvalues to high relative\n*          accuracy.\n*\n*  M       (output) INTEGER\n*          The total number of eigenvalues found.  0 <= M <= N.\n*          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1.\n*\n*  W       (output) DOUBLE PRECISION array, dimension (N)\n*          The first M elements contain the selected eigenvalues in\n*          ascending order.\n*\n*  Z       (output) DOUBLE PRECISION array, dimension (LDZ, max(1,M))\n*          If JOBZ = 'V', then if INFO = 0, the first M columns of Z\n*          contain the orthonormal eigenvectors of the matrix A\n*          corresponding to the selected eigenvalues, with the i-th\n*          column of Z holding the eigenvector associated with W(i).\n*          If JOBZ = 'N', then Z is not referenced.\n*          Note: the user must ensure that at least max(1,M) columns are\n*          supplied in the array Z; if RANGE = 'V', the exact value of M\n*          is not known in advance and an upper bound must be used.\n*          Supplying N columns is always safe.\n*\n*  LDZ     (input) INTEGER\n*          The leading dimension of the array Z.  LDZ >= 1, and if\n*          JOBZ = 'V', LDZ >= max(1,N).\n*\n*  ISUPPZ  (output) INTEGER array, dimension ( 2*max(1,M) )\n*          The support of the eigenvectors in Z, i.e., the indices\n*          indicating the nonzero elements in Z. The i-th eigenvector\n*          is nonzero only in elements ISUPPZ( 2*i-1 ) through\n*          ISUPPZ( 2*i ).\n********** Implemented only for RANGE = 'A' or 'I' and IU - IL = N - 1\n*\n*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK.  LWORK >= max(1,26*N).\n*          For optimal efficiency, LWORK >= (NB+6)*N,\n*          where NB is the max of the blocksize for DSYTRD and DORMTR\n*          returned by ILAENV.\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  IWORK   (workspace/output) INTEGER array, dimension (MAX(1,LIWORK))\n*          On exit, if INFO = 0, IWORK(1) returns the optimal LWORK.\n*\n*  LIWORK  (input) INTEGER\n*          The dimension of the array IWORK.  LIWORK >= max(1,10*N).\n*\n*          If LIWORK = -1, then a workspace query is assumed; the\n*          routine only calculates the optimal size of the IWORK array,\n*          returns this value as the first entry of the IWORK array, and\n*          no error message related to LIWORK is issued by XERBLA.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*          > 0:  Internal error\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Inderjit Dhillon, IBM Almaden, USA\n*     Osni Marques, LBNL/NERSC, USA\n*     Ken Stanley, Computer Science Division, University of\n*       California at Berkeley, USA\n*     Jason Riedy, Computer Science Division, University of\n*       California at Berkeley, USA\n*\n* =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 11)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 11)", argc);
  rb_jobz = argv[0];
  rb_range = argv[1];
  rb_uplo = argv[2];
  rb_a = argv[3];
  rb_vl = argv[4];
  rb_vu = argv[5];
  rb_il = argv[6];
  rb_iu = argv[7];
  rb_abstol = argv[8];
  rb_lwork = argv[9];
  rb_liwork = argv[10];

  jobz = StringValueCStr(rb_jobz)[0];
  range = StringValueCStr(rb_range)[0];
  uplo = StringValueCStr(rb_uplo)[0];
  vl = NUM2DBL(rb_vl);
  vu = NUM2DBL(rb_vu);
  il = NUM2INT(rb_il);
  iu = NUM2INT(rb_iu);
  abstol = NUM2DBL(rb_abstol);
  lwork = NUM2INT(rb_lwork);
  liwork = NUM2INT(rb_liwork);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (4th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (4th argument) must be %d", 2);
  lda = NA_SHAPE0(rb_a);
  n = NA_SHAPE1(rb_a);
  if (NA_TYPE(rb_a) != NA_DFLOAT)
    rb_a = na_change_type(rb_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rb_a, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_w = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  w = NA_PTR_TYPE(rb_w, doublereal*);
  ldz = lsame_(&jobz,"V") ? MAX(1,n) : 1;
  m = lsame_(&range,"I") ? iu-il+1 : n;
  {
    int shape[2];
    shape[0] = ldz;
    shape[1] = MAX(1,m);
    rb_z = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  z = NA_PTR_TYPE(rb_z, doublereal*);
  {
    int shape[1];
    shape[0] = DIM_LEN(2*MAX(1,m));
    rb_isuppz = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  isuppz = NA_PTR_TYPE(rb_isuppz, integer*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rb_work = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rb_work, doublereal*);
  {
    int shape[1];
    shape[0] = DIM_LEN(MAX(1,liwork));
    rb_iwork = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  iwork = NA_PTR_TYPE(rb_iwork, integer*);
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

  dsyevr_(&jobz, &range, &uplo, &n, a, &lda, &vl, &vu, &il, &iu, &abstol, &m, w, z, &ldz, isuppz, work, &lwork, iwork, &liwork, &info);

  rb_m = INT2NUM(m);
  rb_info = INT2NUM(info);
  return rb_ary_new3(8, rb_m, rb_w, rb_z, rb_isuppz, rb_work, rb_iwork, rb_info, rb_a);
}

void
init_lapack_dsyevr(VALUE mLapack){
  rb_define_module_function(mLapack, "dsyevr", rb_dsyevr, -1);
}
