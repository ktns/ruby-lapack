#include "rb_lapack.h"

extern VOID dgesvd_(char *jobu, char *jobvt, integer *m, integer *n, doublereal *a, integer *lda, doublereal *s, doublereal *u, integer *ldu, doublereal *vt, integer *ldvt, doublereal *work, integer *lwork, integer *info);

static VALUE
rb_dgesvd(int argc, VALUE *argv, VALUE self){
  VALUE rb_jobu;
  char jobu; 
  VALUE rb_jobvt;
  char jobvt; 
  VALUE rb_m;
  integer m; 
  VALUE rb_a;
  doublereal *a; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_s;
  doublereal *s; 
  VALUE rb_u;
  doublereal *u; 
  VALUE rb_vt;
  doublereal *vt; 
  VALUE rb_work;
  doublereal *work; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  doublereal *a_out__;

  integer lda;
  integer n;
  integer ldu;
  integer ldvt;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  s, u, vt, work, info, a = NumRu::Lapack.dgesvd( jobu, jobvt, m, a, lwork)\n    or\n  NumRu::Lapack.dgesvd  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  DGESVD computes the singular value decomposition (SVD) of a real\n*  M-by-N matrix A, optionally computing the left and/or right singular\n*  vectors. The SVD is written\n*\n*       A = U * SIGMA * transpose(V)\n*\n*  where SIGMA is an M-by-N matrix which is zero except for its\n*  min(m,n) diagonal elements, U is an M-by-M orthogonal matrix, and\n*  V is an N-by-N orthogonal matrix.  The diagonal elements of SIGMA\n*  are the singular values of A; they are real and non-negative, and\n*  are returned in descending order.  The first min(m,n) columns of\n*  U and V are the left and right singular vectors of A.\n*\n*  Note that the routine returns V**T, not V.\n*\n\n*  Arguments\n*  =========\n*\n*  JOBU    (input) CHARACTER*1\n*          Specifies options for computing all or part of the matrix U:\n*          = 'A':  all M columns of U are returned in array U:\n*          = 'S':  the first min(m,n) columns of U (the left singular\n*                  vectors) are returned in the array U;\n*          = 'O':  the first min(m,n) columns of U (the left singular\n*                  vectors) are overwritten on the array A;\n*          = 'N':  no columns of U (no left singular vectors) are\n*                  computed.\n*\n*  JOBVT   (input) CHARACTER*1\n*          Specifies options for computing all or part of the matrix\n*          V**T:\n*          = 'A':  all N rows of V**T are returned in the array VT;\n*          = 'S':  the first min(m,n) rows of V**T (the right singular\n*                  vectors) are returned in the array VT;\n*          = 'O':  the first min(m,n) rows of V**T (the right singular\n*                  vectors) are overwritten on the array A;\n*          = 'N':  no rows of V**T (no right singular vectors) are\n*                  computed.\n*\n*          JOBVT and JOBU cannot both be 'O'.\n*\n*  M       (input) INTEGER\n*          The number of rows of the input matrix A.  M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the input matrix A.  N >= 0.\n*\n*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)\n*          On entry, the M-by-N matrix A.\n*          On exit,\n*          if JOBU = 'O',  A is overwritten with the first min(m,n)\n*                          columns of U (the left singular vectors,\n*                          stored columnwise);\n*          if JOBVT = 'O', A is overwritten with the first min(m,n)\n*                          rows of V**T (the right singular vectors,\n*                          stored rowwise);\n*          if JOBU .ne. 'O' and JOBVT .ne. 'O', the contents of A\n*                          are destroyed.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,M).\n*\n*  S       (output) DOUBLE PRECISION array, dimension (min(M,N))\n*          The singular values of A, sorted so that S(i) >= S(i+1).\n*\n*  U       (output) DOUBLE PRECISION array, dimension (LDU,UCOL)\n*          (LDU,M) if JOBU = 'A' or (LDU,min(M,N)) if JOBU = 'S'.\n*          If JOBU = 'A', U contains the M-by-M orthogonal matrix U;\n*          if JOBU = 'S', U contains the first min(m,n) columns of U\n*          (the left singular vectors, stored columnwise);\n*          if JOBU = 'N' or 'O', U is not referenced.\n*\n*  LDU     (input) INTEGER\n*          The leading dimension of the array U.  LDU >= 1; if\n*          JOBU = 'S' or 'A', LDU >= M.\n*\n*  VT      (output) DOUBLE PRECISION array, dimension (LDVT,N)\n*          If JOBVT = 'A', VT contains the N-by-N orthogonal matrix\n*          V**T;\n*          if JOBVT = 'S', VT contains the first min(m,n) rows of\n*          V**T (the right singular vectors, stored rowwise);\n*          if JOBVT = 'N' or 'O', VT is not referenced.\n*\n*  LDVT    (input) INTEGER\n*          The leading dimension of the array VT.  LDVT >= 1; if\n*          JOBVT = 'A', LDVT >= N; if JOBVT = 'S', LDVT >= min(M,N).\n*\n*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK;\n*          if INFO > 0, WORK(2:MIN(M,N)) contains the unconverged\n*          superdiagonal elements of an upper bidiagonal matrix B\n*          whose diagonal is in S (not necessarily sorted). B\n*          satisfies A = U * B * VT, so it has the same singular values\n*          as A, and singular vectors related by U and VT.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK.\n*          LWORK >= MAX(1,3*MIN(M,N)+MAX(M,N),5*MIN(M,N)).\n*          For good performance, LWORK should generally be larger.\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*          > 0:  if DBDSQR did not converge, INFO specifies how many\n*                superdiagonals of an intermediate bidiagonal form B\n*                did not converge to zero. See the description of WORK\n*                above for details.\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_jobu = argv[0];
  rb_jobvt = argv[1];
  rb_m = argv[2];
  rb_a = argv[3];
  rb_lwork = argv[4];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (4th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (4th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_DFLOAT)
    rb_a = na_change_type(rb_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rb_a, doublereal*);
  m = NUM2INT(rb_m);
  jobvt = StringValueCStr(rb_jobvt)[0];
  lwork = NUM2INT(rb_lwork);
  jobu = StringValueCStr(rb_jobu)[0];
  ldvt = lsame_(&jobvt,"A") ? n : lsame_(&jobvt,"S") ? MIN(m,n) : 1;
  ldu = ((lsame_(&jobu,"S")) || (lsame_(&jobu,"A"))) ? m : 1;
  {
    int shape[1];
    shape[0] = MIN(m,n);
    rb_s = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  s = NA_PTR_TYPE(rb_s, doublereal*);
  {
    int shape[2];
    shape[0] = ldu;
    shape[1] = lsame_(&jobu,"A") ? m : lsame_(&jobu,"S") ? MIN(m,n) : 0;
    rb_u = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  u = NA_PTR_TYPE(rb_u, doublereal*);
  {
    int shape[2];
    shape[0] = ldvt;
    shape[1] = n;
    rb_vt = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  vt = NA_PTR_TYPE(rb_vt, doublereal*);
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

  dgesvd_(&jobu, &jobvt, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(6, rb_s, rb_u, rb_vt, rb_work, rb_info, rb_a);
}

void
init_lapack_dgesvd(VALUE mLapack){
  rb_define_module_function(mLapack, "dgesvd", rb_dgesvd, -1);
}
