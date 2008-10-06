#include "rb_lapack.h"

static VALUE
rb_dsbgvx(int argc, VALUE *argv, VALUE self){
  VALUE rb_jobz;
  char jobz; 
  VALUE rb_range;
  char range; 
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_ka;
  integer ka; 
  VALUE rb_kb;
  integer kb; 
  VALUE rb_ab;
  doublereal *ab; 
  VALUE rb_bb;
  doublereal *bb; 
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
  VALUE rb_q;
  doublereal *q; 
  VALUE rb_m;
  integer m; 
  VALUE rb_w;
  doublereal *w; 
  VALUE rb_z;
  doublereal *z; 
  VALUE rb_work;
  doublereal *work; 
  VALUE rb_iwork;
  integer *iwork; 
  VALUE rb_ifail;
  integer *ifail; 
  VALUE rb_info;
  integer info; 
  VALUE rb_ab_out__;
  doublereal *ab_out__;
  VALUE rb_bb_out__;
  doublereal *bb_out__;

  integer ldab;
  integer n;
  integer ldbb;
  integer ldq;
  integer ldz;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  q, m, w, z, work, iwork, ifail, info, ab, bb = NumRu::Lapack.dsbgvx( jobz, range, uplo, ka, kb, ab, bb, vl, vu, il, iu, abstol)\n    or\n  NumRu::Lapack.dsbgvx  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DSBGVX( JOBZ, RANGE, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, Q, LDQ, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK, IWORK, IFAIL, INFO )\n\n*  Purpose\n*  =======\n*\n*  DSBGVX computes selected eigenvalues, and optionally, eigenvectors\n*  of a real generalized symmetric-definite banded eigenproblem, of\n*  the form A*x=(lambda)*B*x.  Here A and B are assumed to be symmetric\n*  and banded, and B is also positive definite.  Eigenvalues and\n*  eigenvectors can be selected by specifying either all eigenvalues,\n*  a range of values or a range of indices for the desired eigenvalues.\n*\n\n*  Arguments\n*  =========\n*\n*  JOBZ    (input) CHARACTER*1\n*          = 'N':  Compute eigenvalues only;\n*          = 'V':  Compute eigenvalues and eigenvectors.\n*\n*  RANGE   (input) CHARACTER*1\n*          = 'A': all eigenvalues will be found.\n*          = 'V': all eigenvalues in the half-open interval (VL,VU]\n*                 will be found.\n*          = 'I': the IL-th through IU-th eigenvalues will be found.\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U':  Upper triangles of A and B are stored;\n*          = 'L':  Lower triangles of A and B are stored.\n*\n*  N       (input) INTEGER\n*          The order of the matrices A and B.  N >= 0.\n*\n*  KA      (input) INTEGER\n*          The number of superdiagonals of the matrix A if UPLO = 'U',\n*          or the number of subdiagonals if UPLO = 'L'.  KA >= 0.\n*\n*  KB      (input) INTEGER\n*          The number of superdiagonals of the matrix B if UPLO = 'U',\n*          or the number of subdiagonals if UPLO = 'L'.  KB >= 0.\n*\n*  AB      (input/output) DOUBLE PRECISION array, dimension (LDAB, N)\n*          On entry, the upper or lower triangle of the symmetric band\n*          matrix A, stored in the first ka+1 rows of the array.  The\n*          j-th column of A is stored in the j-th column of the array AB\n*          as follows:\n*          if UPLO = 'U', AB(ka+1+i-j,j) = A(i,j) for max(1,j-ka)<=i<=j;\n*          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+ka).\n*\n*          On exit, the contents of AB are destroyed.\n*\n*  LDAB    (input) INTEGER\n*          The leading dimension of the array AB.  LDAB >= KA+1.\n*\n*  BB      (input/output) DOUBLE PRECISION array, dimension (LDBB, N)\n*          On entry, the upper or lower triangle of the symmetric band\n*          matrix B, stored in the first kb+1 rows of the array.  The\n*          j-th column of B is stored in the j-th column of the array BB\n*          as follows:\n*          if UPLO = 'U', BB(ka+1+i-j,j) = B(i,j) for max(1,j-kb)<=i<=j;\n*          if UPLO = 'L', BB(1+i-j,j)    = B(i,j) for j<=i<=min(n,j+kb).\n*\n*          On exit, the factor S from the split Cholesky factorization\n*          B = S**T*S, as returned by DPBSTF.\n*\n*  LDBB    (input) INTEGER\n*          The leading dimension of the array BB.  LDBB >= KB+1.\n*\n*  Q       (output) DOUBLE PRECISION array, dimension (LDQ, N)\n*          If JOBZ = 'V', the n-by-n matrix used in the reduction of\n*          A*x = (lambda)*B*x to standard form, i.e. C*x = (lambda)*x,\n*          and consequently C to tridiagonal form.\n*          If JOBZ = 'N', the array Q is not referenced.\n*\n*  LDQ     (input) INTEGER\n*          The leading dimension of the array Q.  If JOBZ = 'N',\n*          LDQ >= 1. If JOBZ = 'V', LDQ >= max(1,N).\n*\n*  VL      (input) DOUBLE PRECISION\n*  VU      (input) DOUBLE PRECISION\n*          If RANGE='V', the lower and upper bounds of the interval to\n*          be searched for eigenvalues. VL < VU.\n*          Not referenced if RANGE = 'A' or 'I'.\n*\n*  IL      (input) INTEGER\n*  IU      (input) INTEGER\n*          If RANGE='I', the indices (in ascending order) of the\n*          smallest and largest eigenvalues to be returned.\n*          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.\n*          Not referenced if RANGE = 'A' or 'V'.\n*\n*  ABSTOL  (input) DOUBLE PRECISION\n*          The absolute error tolerance for the eigenvalues.\n*          An approximate eigenvalue is accepted as converged\n*          when it is determined to lie in an interval [a,b]\n*          of width less than or equal to\n*\n*                  ABSTOL + EPS *   max( |a|,|b| ) ,\n*\n*          where EPS is the machine precision.  If ABSTOL is less than\n*          or equal to zero, then  EPS*|T|  will be used in its place,\n*          where |T| is the 1-norm of the tridiagonal matrix obtained\n*          by reducing A to tridiagonal form.\n*\n*          Eigenvalues will be computed most accurately when ABSTOL is\n*          set to twice the underflow threshold 2*DLAMCH('S'), not zero.\n*          If this routine returns with INFO>0, indicating that some\n*          eigenvectors did not converge, try setting ABSTOL to\n*          2*DLAMCH('S').\n*\n*  M       (output) INTEGER\n*          The total number of eigenvalues found.  0 <= M <= N.\n*          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1.\n*\n*  W       (output) DOUBLE PRECISION array, dimension (N)\n*          If INFO = 0, the eigenvalues in ascending order.\n*\n*  Z       (output) DOUBLE PRECISION array, dimension (LDZ, N)\n*          If JOBZ = 'V', then if INFO = 0, Z contains the matrix Z of\n*          eigenvectors, with the i-th column of Z holding the\n*          eigenvector associated with W(i).  The eigenvectors are\n*          normalized so Z**T*B*Z = I.\n*          If JOBZ = 'N', then Z is not referenced.\n*\n*  LDZ     (input) INTEGER\n*          The leading dimension of the array Z.  LDZ >= 1, and if\n*          JOBZ = 'V', LDZ >= max(1,N).\n*\n*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (7*N)\n*\n*  IWORK   (workspace/output) INTEGER array, dimension (5*N)\n*\n*  IFAIL   (output) INTEGER array, dimension (M)\n*          If JOBZ = 'V', then if INFO = 0, the first M elements of\n*          IFAIL are zero.  If INFO > 0, then IFAIL contains the\n*          indices of the eigenvalues that failed to converge.\n*          If JOBZ = 'N', then IFAIL is not referenced.\n*\n*  INFO    (output) INTEGER\n*          = 0 : successful exit\n*          < 0 : if INFO = -i, the i-th argument had an illegal value\n*          <= N: if INFO = i, then i eigenvectors failed to converge.\n*                  Their indices are stored in IFAIL.\n*          > N : DPBSTF returned an error code; i.e.,\n*                if INFO = N + i, for 1 <= i <= N, then the leading\n*                minor of order i of B is not positive definite.\n*                The factorization of B could not be completed and\n*                no eigenvalues or eigenvectors were computed.\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Mark Fahey, Department of Mathematics, Univ. of Kentucky, USA\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 12)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 12)", argc);
  rb_jobz = argv[0];
  rb_range = argv[1];
  rb_uplo = argv[2];
  rb_ka = argv[3];
  rb_kb = argv[4];
  rb_ab = argv[5];
  rb_bb = argv[6];
  rb_vl = argv[7];
  rb_vu = argv[8];
  rb_il = argv[9];
  rb_iu = argv[10];
  rb_abstol = argv[11];

  jobz = StringValueCStr(rb_jobz)[0];
  range = StringValueCStr(rb_range)[0];
  uplo = StringValueCStr(rb_uplo)[0];
  ka = NUM2INT(rb_ka);
  kb = NUM2INT(rb_kb);
  vl = NUM2DBL(rb_vl);
  vu = NUM2DBL(rb_vu);
  il = NUM2INT(rb_il);
  iu = NUM2INT(rb_iu);
  abstol = NUM2DBL(rb_abstol);
  if (!NA_IsNArray(rb_ab))
    rb_raise(rb_eArgError, "ab (6th argument) must be NArray");
  if (NA_RANK(rb_ab) != 2)
    rb_raise(rb_eArgError, "rank of ab (6th argument) must be %d", 2);
  ldab = NA_SHAPE0(rb_ab);
  n = NA_SHAPE1(rb_ab);
  if (NA_TYPE(rb_ab) != NA_DFLOAT)
    rb_ab = na_change_type(rb_ab, NA_DFLOAT);
  ab = NA_PTR_TYPE(rb_ab, doublereal*);
  if (!NA_IsNArray(rb_bb))
    rb_raise(rb_eArgError, "bb (7th argument) must be NArray");
  if (NA_RANK(rb_bb) != 2)
    rb_raise(rb_eArgError, "rank of bb (7th argument) must be %d", 2);
  ldbb = NA_SHAPE0(rb_bb);
  if (NA_SHAPE1(rb_bb) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of bb must be the same as shape 1 of ab");
  if (NA_TYPE(rb_bb) != NA_DFLOAT)
    rb_bb = na_change_type(rb_bb, NA_DFLOAT);
  bb = NA_PTR_TYPE(rb_bb, doublereal*);
  ldq = 1 ? jobz = 'n' : max(1,n) ? jobz = 'v' : 0;
  {
    int shape[2];
    shape[0] = ldq;
    shape[1] = n;
    rb_q = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  q = NA_PTR_TYPE(rb_q, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_w = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  w = NA_PTR_TYPE(rb_w, doublereal*);
  ldz = lsame_(&jobz,"V") ? MAX(1,n) : 1;
  {
    int shape[2];
    shape[0] = ldz;
    shape[1] = n;
    rb_z = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  z = NA_PTR_TYPE(rb_z, doublereal*);
  {
    int shape[1];
    shape[0] = 7*n;
    rb_work = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rb_work, doublereal*);
  {
    int shape[1];
    shape[0] = 5*n;
    rb_iwork = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  iwork = NA_PTR_TYPE(rb_iwork, integer*);
  m = lsame_(&range,"A") ? n : lsame_(&range,"I") ? iu-il+1 : 0;
  {
    int shape[1];
    shape[0] = m;
    rb_ifail = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  ifail = NA_PTR_TYPE(rb_ifail, integer*);
  {
    int shape[2];
    shape[0] = ldab;
    shape[1] = n;
    rb_ab_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  ab_out__ = NA_PTR_TYPE(rb_ab_out__, doublereal*);
  MEMCPY(ab_out__, ab, doublereal, NA_TOTAL(rb_ab));
  rb_ab = rb_ab_out__;
  ab = ab_out__;
  {
    int shape[2];
    shape[0] = ldbb;
    shape[1] = n;
    rb_bb_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  bb_out__ = NA_PTR_TYPE(rb_bb_out__, doublereal*);
  MEMCPY(bb_out__, bb, doublereal, NA_TOTAL(rb_bb));
  rb_bb = rb_bb_out__;
  bb = bb_out__;

  dsbgvx_(&jobz, &range, &uplo, &n, &ka, &kb, ab, &ldab, bb, &ldbb, q, &ldq, &vl, &vu, &il, &iu, &abstol, &m, w, z, &ldz, work, iwork, ifail, &info);

  rb_m = INT2NUM(m);
  rb_info = INT2NUM(info);
  return rb_ary_new3(10, rb_q, rb_m, rb_w, rb_z, rb_work, rb_iwork, rb_ifail, rb_info, rb_ab, rb_bb);
}

void
init_lapack_dsbgvx(VALUE mLapack){
  rb_define_module_function(mLapack, "dsbgvx", rb_dsbgvx, -1);
}
