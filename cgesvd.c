#include "rb_lapack.h"

extern VOID cgesvd_(char *jobu, char *jobvt, integer *m, integer *n, complex *a, integer *lda, real *s, complex *u, integer *ldu, complex *vt, integer *ldvt, complex *work, integer *lwork, real *rwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_cgesvd(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_jobu;
  char jobu; 
  VALUE rblapack_jobvt;
  char jobvt; 
  VALUE rblapack_m;
  integer m; 
  VALUE rblapack_a;
  complex *a; 
  VALUE rblapack_lwork;
  integer lwork; 
  VALUE rblapack_s;
  real *s; 
  VALUE rblapack_u;
  complex *u; 
  VALUE rblapack_vt;
  complex *vt; 
  VALUE rblapack_work;
  complex *work; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_a_out__;
  complex *a_out__;
  real *rwork;

  integer lda;
  integer n;
  integer ldu;
  integer ldvt;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  s, u, vt, work, info, a = NumRu::Lapack.cgesvd( jobu, jobvt, m, a, lwork, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE CGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, RWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  CGESVD computes the singular value decomposition (SVD) of a complex\n*  M-by-N matrix A, optionally computing the left and/or right singular\n*  vectors. The SVD is written\n*\n*       A = U * SIGMA * conjugate-transpose(V)\n*\n*  where SIGMA is an M-by-N matrix which is zero except for its\n*  min(m,n) diagonal elements, U is an M-by-M unitary matrix, and\n*  V is an N-by-N unitary matrix.  The diagonal elements of SIGMA\n*  are the singular values of A; they are real and non-negative, and\n*  are returned in descending order.  The first min(m,n) columns of\n*  U and V are the left and right singular vectors of A.\n*\n*  Note that the routine returns V**H, not V.\n*\n\n*  Arguments\n*  =========\n*\n*  JOBU    (input) CHARACTER*1\n*          Specifies options for computing all or part of the matrix U:\n*          = 'A':  all M columns of U are returned in array U:\n*          = 'S':  the first min(m,n) columns of U (the left singular\n*                  vectors) are returned in the array U;\n*          = 'O':  the first min(m,n) columns of U (the left singular\n*                  vectors) are overwritten on the array A;\n*          = 'N':  no columns of U (no left singular vectors) are\n*                  computed.\n*\n*  JOBVT   (input) CHARACTER*1\n*          Specifies options for computing all or part of the matrix\n*          V**H:\n*          = 'A':  all N rows of V**H are returned in the array VT;\n*          = 'S':  the first min(m,n) rows of V**H (the right singular\n*                  vectors) are returned in the array VT;\n*          = 'O':  the first min(m,n) rows of V**H (the right singular\n*                  vectors) are overwritten on the array A;\n*          = 'N':  no rows of V**H (no right singular vectors) are\n*                  computed.\n*\n*          JOBVT and JOBU cannot both be 'O'.\n*\n*  M       (input) INTEGER\n*          The number of rows of the input matrix A.  M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the input matrix A.  N >= 0.\n*\n*  A       (input/output) COMPLEX array, dimension (LDA,N)\n*          On entry, the M-by-N matrix A.\n*          On exit,\n*          if JOBU = 'O',  A is overwritten with the first min(m,n)\n*                          columns of U (the left singular vectors,\n*                          stored columnwise);\n*          if JOBVT = 'O', A is overwritten with the first min(m,n)\n*                          rows of V**H (the right singular vectors,\n*                          stored rowwise);\n*          if JOBU .ne. 'O' and JOBVT .ne. 'O', the contents of A\n*                          are destroyed.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,M).\n*\n*  S       (output) REAL array, dimension (min(M,N))\n*          The singular values of A, sorted so that S(i) >= S(i+1).\n*\n*  U       (output) COMPLEX array, dimension (LDU,UCOL)\n*          (LDU,M) if JOBU = 'A' or (LDU,min(M,N)) if JOBU = 'S'.\n*          If JOBU = 'A', U contains the M-by-M unitary matrix U;\n*          if JOBU = 'S', U contains the first min(m,n) columns of U\n*          (the left singular vectors, stored columnwise);\n*          if JOBU = 'N' or 'O', U is not referenced.\n*\n*  LDU     (input) INTEGER\n*          The leading dimension of the array U.  LDU >= 1; if\n*          JOBU = 'S' or 'A', LDU >= M.\n*\n*  VT      (output) COMPLEX array, dimension (LDVT,N)\n*          If JOBVT = 'A', VT contains the N-by-N unitary matrix\n*          V**H;\n*          if JOBVT = 'S', VT contains the first min(m,n) rows of\n*          V**H (the right singular vectors, stored rowwise);\n*          if JOBVT = 'N' or 'O', VT is not referenced.\n*\n*  LDVT    (input) INTEGER\n*          The leading dimension of the array VT.  LDVT >= 1; if\n*          JOBVT = 'A', LDVT >= N; if JOBVT = 'S', LDVT >= min(M,N).\n*\n*  WORK    (workspace/output) COMPLEX array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK.\n*          LWORK >=  MAX(1,2*MIN(M,N)+MAX(M,N)).\n*          For good performance, LWORK should generally be larger.\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  RWORK   (workspace) REAL array, dimension (5*min(M,N))\n*          On exit, if INFO > 0, RWORK(1:MIN(M,N)-1) contains the\n*          unconverged superdiagonal elements of an upper bidiagonal\n*          matrix B whose diagonal is in S (not necessarily sorted).\n*          B satisfies A = U * B * VT, so it has the same singular\n*          values as A, and singular vectors related by U and VT.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*          > 0:  if CBDSQR did not converge, INFO specifies how many\n*                superdiagonals of an intermediate bidiagonal form B\n*                did not converge to zero. See the description of RWORK\n*                above for details.\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  s, u, vt, work, info, a = NumRu::Lapack.cgesvd( jobu, jobvt, m, a, lwork, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rblapack_jobu = argv[0];
  rblapack_jobvt = argv[1];
  rblapack_m = argv[2];
  rblapack_a = argv[3];
  rblapack_lwork = argv[4];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (4th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (4th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_a);
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_SCOMPLEX)
    rblapack_a = na_change_type(rblapack_a, NA_SCOMPLEX);
  a = NA_PTR_TYPE(rblapack_a, complex*);
  m = NUM2INT(rblapack_m);
  jobvt = StringValueCStr(rblapack_jobvt)[0];
  lwork = NUM2INT(rblapack_lwork);
  jobu = StringValueCStr(rblapack_jobu)[0];
  ldvt = lsame_(&jobvt,"A") ? n : lsame_(&jobvt,"S") ? MIN(m,n) : 1;
  ldu = ((lsame_(&jobu,"S")) || (lsame_(&jobu,"A"))) ? m : 1;
  {
    int shape[1];
    shape[0] = MIN(m,n);
    rblapack_s = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  s = NA_PTR_TYPE(rblapack_s, real*);
  {
    int shape[2];
    shape[0] = ldu;
    shape[1] = lsame_(&jobu,"A") ? m : lsame_(&jobu,"S") ? MIN(m,n) : 0;
    rblapack_u = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  u = NA_PTR_TYPE(rblapack_u, complex*);
  {
    int shape[2];
    shape[0] = ldvt;
    shape[1] = n;
    rblapack_vt = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  vt = NA_PTR_TYPE(rblapack_vt, complex*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rblapack_work = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rblapack_work, complex*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rblapack_a_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rblapack_a_out__, complex*);
  MEMCPY(a_out__, a, complex, NA_TOTAL(rblapack_a));
  rblapack_a = rblapack_a_out__;
  a = a_out__;
  rwork = ALLOC_N(real, (5*MIN(m,n)));

  cgesvd_(&jobu, &jobvt, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, rwork, &info);

  free(rwork);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(6, rblapack_s, rblapack_u, rblapack_vt, rblapack_work, rblapack_info, rblapack_a);
}

void
init_lapack_cgesvd(VALUE mLapack){
  rb_define_module_function(mLapack, "cgesvd", rblapack_cgesvd, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
