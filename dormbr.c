#include "rb_lapack.h"

extern VOID dormbr_(char *vect, char *side, char *trans, integer *m, integer *n, integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *c, integer *ldc, doublereal *work, integer *lwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dormbr(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_vect;
  char vect; 
  VALUE rblapack_side;
  char side; 
  VALUE rblapack_trans;
  char trans; 
  VALUE rblapack_m;
  integer m; 
  VALUE rblapack_k;
  integer k; 
  VALUE rblapack_a;
  doublereal *a; 
  VALUE rblapack_tau;
  doublereal *tau; 
  VALUE rblapack_c;
  doublereal *c; 
  VALUE rblapack_lwork;
  integer lwork; 
  VALUE rblapack_work;
  doublereal *work; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_c_out__;
  doublereal *c_out__;

  integer lda;
  integer ldc;
  integer n;
  integer nq;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  work, info, c = NumRu::Lapack.dormbr( vect, side, trans, m, k, a, tau, c, lwork, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DORMBR( VECT, SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, LWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  If VECT = 'Q', DORMBR overwrites the general real M-by-N matrix C\n*  with\n*                  SIDE = 'L'     SIDE = 'R'\n*  TRANS = 'N':      Q * C          C * Q\n*  TRANS = 'T':      Q**T * C       C * Q**T\n*\n*  If VECT = 'P', DORMBR overwrites the general real M-by-N matrix C\n*  with\n*                  SIDE = 'L'     SIDE = 'R'\n*  TRANS = 'N':      P * C          C * P\n*  TRANS = 'T':      P**T * C       C * P**T\n*\n*  Here Q and P**T are the orthogonal matrices determined by DGEBRD when\n*  reducing a real matrix A to bidiagonal form: A = Q * B * P**T. Q and\n*  P**T are defined as products of elementary reflectors H(i) and G(i)\n*  respectively.\n*\n*  Let nq = m if SIDE = 'L' and nq = n if SIDE = 'R'. Thus nq is the\n*  order of the orthogonal matrix Q or P**T that is applied.\n*\n*  If VECT = 'Q', A is assumed to have been an NQ-by-K matrix:\n*  if nq >= k, Q = H(1) H(2) . . . H(k);\n*  if nq < k, Q = H(1) H(2) . . . H(nq-1).\n*\n*  If VECT = 'P', A is assumed to have been a K-by-NQ matrix:\n*  if k < nq, P = G(1) G(2) . . . G(k);\n*  if k >= nq, P = G(1) G(2) . . . G(nq-1).\n*\n\n*  Arguments\n*  =========\n*\n*  VECT    (input) CHARACTER*1\n*          = 'Q': apply Q or Q**T;\n*          = 'P': apply P or P**T.\n*\n*  SIDE    (input) CHARACTER*1\n*          = 'L': apply Q, Q**T, P or P**T from the Left;\n*          = 'R': apply Q, Q**T, P or P**T from the Right.\n*\n*  TRANS   (input) CHARACTER*1\n*          = 'N':  No transpose, apply Q  or P;\n*          = 'T':  Transpose, apply Q**T or P**T.\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix C. M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix C. N >= 0.\n*\n*  K       (input) INTEGER\n*          If VECT = 'Q', the number of columns in the original\n*          matrix reduced by DGEBRD.\n*          If VECT = 'P', the number of rows in the original\n*          matrix reduced by DGEBRD.\n*          K >= 0.\n*\n*  A       (input) DOUBLE PRECISION array, dimension\n*                                (LDA,min(nq,K)) if VECT = 'Q'\n*                                (LDA,nq)        if VECT = 'P'\n*          The vectors which define the elementary reflectors H(i) and\n*          G(i), whose products determine the matrices Q and P, as\n*          returned by DGEBRD.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.\n*          If VECT = 'Q', LDA >= max(1,nq);\n*          if VECT = 'P', LDA >= max(1,min(nq,K)).\n*\n*  TAU     (input) DOUBLE PRECISION array, dimension (min(nq,K))\n*          TAU(i) must contain the scalar factor of the elementary\n*          reflector H(i) or G(i) which determines Q or P, as returned\n*          by DGEBRD in the array argument TAUQ or TAUP.\n*\n*  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)\n*          On entry, the M-by-N matrix C.\n*          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q\n*          or P*C or P**T*C or C*P or C*P**T.\n*\n*  LDC     (input) INTEGER\n*          The leading dimension of the array C. LDC >= max(1,M).\n*\n*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK.\n*          If SIDE = 'L', LWORK >= max(1,N);\n*          if SIDE = 'R', LWORK >= max(1,M).\n*          For optimum performance LWORK >= N*NB if SIDE = 'L', and\n*          LWORK >= M*NB if SIDE = 'R', where NB is the optimal\n*          blocksize.\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*\n\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      LOGICAL            APPLYQ, LEFT, LQUERY, NOTRAN\n      CHARACTER          TRANST\n      INTEGER            I1, I2, IINFO, LWKOPT, MI, NB, NI, NQ, NW\n*     ..\n*     .. External Functions ..\n      LOGICAL            LSAME\n      INTEGER            ILAENV\n      EXTERNAL           LSAME, ILAENV\n*     ..\n*     .. External Subroutines ..\n      EXTERNAL           DORMLQ, DORMQR, XERBLA\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          MAX, MIN\n*     ..\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  work, info, c = NumRu::Lapack.dormbr( vect, side, trans, m, k, a, tau, c, lwork, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 9)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 9)", argc);
  rblapack_vect = argv[0];
  rblapack_side = argv[1];
  rblapack_trans = argv[2];
  rblapack_m = argv[3];
  rblapack_k = argv[4];
  rblapack_a = argv[5];
  rblapack_tau = argv[6];
  rblapack_c = argv[7];
  rblapack_lwork = argv[8];
  if (rb_options != Qnil) {
  }

  k = NUM2INT(rblapack_k);
  lwork = NUM2INT(rblapack_lwork);
  side = StringValueCStr(rblapack_side)[0];
  m = NUM2INT(rblapack_m);
  vect = StringValueCStr(rblapack_vect)[0];
  if (!NA_IsNArray(rblapack_c))
    rb_raise(rb_eArgError, "c (8th argument) must be NArray");
  if (NA_RANK(rblapack_c) != 2)
    rb_raise(rb_eArgError, "rank of c (8th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_c);
  ldc = NA_SHAPE0(rblapack_c);
  if (NA_TYPE(rblapack_c) != NA_DFLOAT)
    rblapack_c = na_change_type(rblapack_c, NA_DFLOAT);
  c = NA_PTR_TYPE(rblapack_c, doublereal*);
  trans = StringValueCStr(rblapack_trans)[0];
  nq = lsame_(&side,"L") ? m : lsame_(&side,"R") ? n : 0;
  if (!NA_IsNArray(rblapack_tau))
    rb_raise(rb_eArgError, "tau (7th argument) must be NArray");
  if (NA_RANK(rblapack_tau) != 1)
    rb_raise(rb_eArgError, "rank of tau (7th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_tau) != (MIN(nq,k)))
    rb_raise(rb_eRuntimeError, "shape 0 of tau must be %d", MIN(nq,k));
  if (NA_TYPE(rblapack_tau) != NA_DFLOAT)
    rblapack_tau = na_change_type(rblapack_tau, NA_DFLOAT);
  tau = NA_PTR_TYPE(rblapack_tau, doublereal*);
  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (6th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (6th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_a) != (MIN(nq,k)))
    rb_raise(rb_eRuntimeError, "shape 1 of a must be %d", MIN(nq,k));
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_DFLOAT)
    rblapack_a = na_change_type(rblapack_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rblapack_a, doublereal*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rblapack_work = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rblapack_work, doublereal*);
  {
    int shape[2];
    shape[0] = ldc;
    shape[1] = n;
    rblapack_c_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  c_out__ = NA_PTR_TYPE(rblapack_c_out__, doublereal*);
  MEMCPY(c_out__, c, doublereal, NA_TOTAL(rblapack_c));
  rblapack_c = rblapack_c_out__;
  c = c_out__;

  dormbr_(&vect, &side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork, &info);

  rblapack_info = INT2NUM(info);
  return rb_ary_new3(3, rblapack_work, rblapack_info, rblapack_c);
}

void
init_lapack_dormbr(VALUE mLapack){
  rb_define_module_function(mLapack, "dormbr", rblapack_dormbr, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
