#include "rb_lapack.h"

static VALUE
rb_cgesdd(int argc, VALUE *argv, VALUE self){
  VALUE rb_jobz;
  char jobz; 
  VALUE rb_m;
  integer m; 
  VALUE rb_a;
  complex *a; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_s;
  real *s; 
  VALUE rb_u;
  complex *u; 
  VALUE rb_vt;
  complex *vt; 
  VALUE rb_work;
  complex *work; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  complex *a_out__;
  real *rwork;
  integer *iwork;

  integer lda;
  integer n;
  integer ldu;
  integer ucol;
  integer ldvt;
  integer lrwork;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  s, u, vt, work, info, a = NumRu::Lapack.cgesdd( jobz, m, a, lwork)\n    or\n  NumRu::Lapack.cgesdd  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE CGESDD( JOBZ, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, RWORK, IWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  CGESDD computes the singular value decomposition (SVD) of a complex\n*  M-by-N matrix A, optionally computing the left and/or right singular\n*  vectors, by using divide-and-conquer method. The SVD is written\n*\n*       A = U * SIGMA * conjugate-transpose(V)\n*\n*  where SIGMA is an M-by-N matrix which is zero except for its\n*  min(m,n) diagonal elements, U is an M-by-M unitary matrix, and\n*  V is an N-by-N unitary matrix.  The diagonal elements of SIGMA\n*  are the singular values of A; they are real and non-negative, and\n*  are returned in descending order.  The first min(m,n) columns of\n*  U and V are the left and right singular vectors of A.\n*\n*  Note that the routine returns VT = V**H, not V.\n*\n*  The divide and conquer algorithm makes very mild assumptions about\n*  floating point arithmetic. It will work on machines with a guard\n*  digit in add/subtract, or on those binary machines without guard\n*  digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or\n*  Cray-2. It could conceivably fail on hexadecimal or decimal machines\n*  without guard digits, but we know of none.\n*\n\n*  Arguments\n*  =========\n*\n*  JOBZ    (input) CHARACTER*1\n*          Specifies options for computing all or part of the matrix U:\n*          = 'A':  all M columns of U and all N rows of V**H are\n*                  returned in the arrays U and VT;\n*          = 'S':  the first min(M,N) columns of U and the first\n*                  min(M,N) rows of V**H are returned in the arrays U\n*                  and VT;\n*          = 'O':  If M >= N, the first N columns of U are overwritten\n*                  in the array A and all rows of V**H are returned in\n*                  the array VT;\n*                  otherwise, all columns of U are returned in the\n*                  array U and the first M rows of V**H are overwritten\n*                  in the array A;\n*          = 'N':  no columns of U or rows of V**H are computed.\n*\n*  M       (input) INTEGER\n*          The number of rows of the input matrix A.  M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the input matrix A.  N >= 0.\n*\n*  A       (input/output) COMPLEX array, dimension (LDA,N)\n*          On entry, the M-by-N matrix A.\n*          On exit,\n*          if JOBZ = 'O',  A is overwritten with the first N columns\n*                          of U (the left singular vectors, stored\n*                          columnwise) if M >= N;\n*                          A is overwritten with the first M rows\n*                          of V**H (the right singular vectors, stored\n*                          rowwise) otherwise.\n*          if JOBZ .ne. 'O', the contents of A are destroyed.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,M).\n*\n*  S       (output) REAL array, dimension (min(M,N))\n*          The singular values of A, sorted so that S(i) >= S(i+1).\n*\n*  U       (output) COMPLEX array, dimension (LDU,UCOL)\n*          UCOL = M if JOBZ = 'A' or JOBZ = 'O' and M < N;\n*          UCOL = min(M,N) if JOBZ = 'S'.\n*          If JOBZ = 'A' or JOBZ = 'O' and M < N, U contains the M-by-M\n*          unitary matrix U;\n*          if JOBZ = 'S', U contains the first min(M,N) columns of U\n*          (the left singular vectors, stored columnwise);\n*          if JOBZ = 'O' and M >= N, or JOBZ = 'N', U is not referenced.\n*\n*  LDU     (input) INTEGER\n*          The leading dimension of the array U.  LDU >= 1; if\n*          JOBZ = 'S' or 'A' or JOBZ = 'O' and M < N, LDU >= M.\n*\n*  VT      (output) COMPLEX array, dimension (LDVT,N)\n*          If JOBZ = 'A' or JOBZ = 'O' and M >= N, VT contains the\n*          N-by-N unitary matrix V**H;\n*          if JOBZ = 'S', VT contains the first min(M,N) rows of\n*          V**H (the right singular vectors, stored rowwise);\n*          if JOBZ = 'O' and M < N, or JOBZ = 'N', VT is not referenced.\n*\n*  LDVT    (input) INTEGER\n*          The leading dimension of the array VT.  LDVT >= 1; if\n*          JOBZ = 'A' or JOBZ = 'O' and M >= N, LDVT >= N;\n*          if JOBZ = 'S', LDVT >= min(M,N).\n*\n*  WORK    (workspace/output) COMPLEX array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK. LWORK >= 1.\n*          if JOBZ = 'N', LWORK >= 2*min(M,N)+max(M,N).\n*          if JOBZ = 'O',\n*                LWORK >= 2*min(M,N)*min(M,N)+2*min(M,N)+max(M,N).\n*          if JOBZ = 'S' or 'A',\n*                LWORK >= min(M,N)*min(M,N)+2*min(M,N)+max(M,N).\n*          For good performance, LWORK should generally be larger.\n*\n*          If LWORK = -1, a workspace query is assumed.  The optimal\n*          size for the WORK array is calculated and stored in WORK(1),\n*          and no other work except argument checking is performed.\n*\n*  RWORK   (workspace) REAL array, dimension (MAX(1,LRWORK))\n*          If JOBZ = 'N', LRWORK >= 5*min(M,N).\n*          Otherwise, LRWORK >= 5*min(M,N)*min(M,N) + 7*min(M,N)\n*\n*  IWORK   (workspace) INTEGER array, dimension (8*min(M,N))\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*          > 0:  The updating process of SBDSDC did not converge.\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Ming Gu and Huan Ren, Computer Science Division, University of\n*     California at Berkeley, USA\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_jobz = argv[0];
  rb_m = argv[1];
  rb_a = argv[2];
  rb_lwork = argv[3];

  jobz = StringValueCStr(rb_jobz)[0];
  m = NUM2INT(rb_m);
  lwork = NUM2INT(rb_lwork);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  lda = NA_SHAPE0(rb_a);
  n = NA_SHAPE1(rb_a);
  if (NA_TYPE(rb_a) != NA_SCOMPLEX)
    rb_a = na_change_type(rb_a, NA_SCOMPLEX);
  a = NA_PTR_TYPE(rb_a, complex*);
  {
    int shape[1];
    shape[0] = MIN(m,n);
    rb_s = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  s = NA_PTR_TYPE(rb_s, real*);
  ucol = ((lsame_(&jobz,"A")) || (((lsame_(&jobz,"O")) && (m < n)))) ? m : lsame_(&jobz,"S") ? MIN(m,n) : 0;
  ldu = ((lsame_(&jobz,"S")) || ((('a') || (((lsame_(&jobz,"O")) && (m < n)))))) ? m : 1;
  {
    int shape[2];
    shape[0] = ldu;
    shape[1] = ucol;
    rb_u = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  u = NA_PTR_TYPE(rb_u, complex*);
  ldvt = ((lsame_(&jobz,"A")) || (((lsame_(&jobz,"O")) && (m == n)))) ? n : lsame_(&jobz,"S") ? MIN(m,n) : 1;
  {
    int shape[2];
    shape[0] = ldvt;
    shape[1] = n;
    rb_vt = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  vt = NA_PTR_TYPE(rb_vt, complex*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rb_work = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rb_work, complex*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rb_a_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rb_a_out__, complex*);
  MEMCPY(a_out__, a, complex, NA_TOTAL(rb_a));
  rb_a = rb_a_out__;
  a = a_out__;
  lrwork = lsame_(&jobz,"N") ? 5*MIN(m,n) : 5*MIN(m,n)*MIN(m,n) + 7*MIN(m,n);
  rwork = ALLOC_N(real, (MAX(1,lrwork)));
  iwork = ALLOC_N(integer, (8*MIN(m,n)));

  cgesdd_(&jobz, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, rwork, iwork, &info);

  free(rwork);
  free(iwork);
  rb_info = INT2NUM(info);
  return rb_ary_new3(6, rb_s, rb_u, rb_vt, rb_work, rb_info, rb_a);
}

void
init_lapack_cgesdd(VALUE mLapack){
  rb_define_module_function(mLapack, "cgesdd", rb_cgesdd, -1);
}
