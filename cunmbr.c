#include "rb_lapack.h"

extern VOID cunmbr_(char *vect, char *side, char *trans, integer *m, integer *n, integer *k, complex *a, integer *lda, complex *tau, complex *c, integer *ldc, complex *work, integer *lwork, integer *info);

static VALUE
rb_cunmbr(int argc, VALUE *argv, VALUE self){
  VALUE rb_vect;
  char vect; 
  VALUE rb_side;
  char side; 
  VALUE rb_trans;
  char trans; 
  VALUE rb_m;
  integer m; 
  VALUE rb_k;
  integer k; 
  VALUE rb_a;
  complex *a; 
  VALUE rb_tau;
  complex *tau; 
  VALUE rb_c;
  complex *c; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_work;
  complex *work; 
  VALUE rb_info;
  integer info; 
  VALUE rb_c_out__;
  complex *c_out__;
  integer nq;

  integer lda;
  integer ldc;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  work, info, c = NumRu::Lapack.cunmbr( vect, side, trans, m, k, a, tau, c, lwork)\n    or\n  NumRu::Lapack.cunmbr  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE CUNMBR( VECT, SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, LWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  If VECT = 'Q', CUNMBR overwrites the general complex M-by-N matrix C\n*  with\n*                  SIDE = 'L'     SIDE = 'R'\n*  TRANS = 'N':      Q * C          C * Q\n*  TRANS = 'C':      Q**H * C       C * Q**H\n*\n*  If VECT = 'P', CUNMBR overwrites the general complex M-by-N matrix C\n*  with\n*                  SIDE = 'L'     SIDE = 'R'\n*  TRANS = 'N':      P * C          C * P\n*  TRANS = 'C':      P**H * C       C * P**H\n*\n*  Here Q and P**H are the unitary matrices determined by CGEBRD when\n*  reducing a complex matrix A to bidiagonal form: A = Q * B * P**H. Q\n*  and P**H are defined as products of elementary reflectors H(i) and\n*  G(i) respectively.\n*\n*  Let nq = m if SIDE = 'L' and nq = n if SIDE = 'R'. Thus nq is the\n*  order of the unitary matrix Q or P**H that is applied.\n*\n*  If VECT = 'Q', A is assumed to have been an NQ-by-K matrix:\n*  if nq >= k, Q = H(1) H(2) . . . H(k);\n*  if nq < k, Q = H(1) H(2) . . . H(nq-1).\n*\n*  If VECT = 'P', A is assumed to have been a K-by-NQ matrix:\n*  if k < nq, P = G(1) G(2) . . . G(k);\n*  if k >= nq, P = G(1) G(2) . . . G(nq-1).\n*\n\n*  Arguments\n*  =========\n*\n*  VECT    (input) CHARACTER*1\n*          = 'Q': apply Q or Q**H;\n*          = 'P': apply P or P**H.\n*\n*  SIDE    (input) CHARACTER*1\n*          = 'L': apply Q, Q**H, P or P**H from the Left;\n*          = 'R': apply Q, Q**H, P or P**H from the Right.\n*\n*  TRANS   (input) CHARACTER*1\n*          = 'N':  No transpose, apply Q or P;\n*          = 'C':  Conjugate transpose, apply Q**H or P**H.\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix C. M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix C. N >= 0.\n*\n*  K       (input) INTEGER\n*          If VECT = 'Q', the number of columns in the original\n*          matrix reduced by CGEBRD.\n*          If VECT = 'P', the number of rows in the original\n*          matrix reduced by CGEBRD.\n*          K >= 0.\n*\n*  A       (input) COMPLEX array, dimension\n*                                (LDA,min(nq,K)) if VECT = 'Q'\n*                                (LDA,nq)        if VECT = 'P'\n*          The vectors which define the elementary reflectors H(i) and\n*          G(i), whose products determine the matrices Q and P, as\n*          returned by CGEBRD.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.\n*          If VECT = 'Q', LDA >= max(1,nq);\n*          if VECT = 'P', LDA >= max(1,min(nq,K)).\n*\n*  TAU     (input) COMPLEX array, dimension (min(nq,K))\n*          TAU(i) must contain the scalar factor of the elementary\n*          reflector H(i) or G(i) which determines Q or P, as returned\n*          by CGEBRD in the array argument TAUQ or TAUP.\n*\n*  C       (input/output) COMPLEX array, dimension (LDC,N)\n*          On entry, the M-by-N matrix C.\n*          On exit, C is overwritten by Q*C or Q**H*C or C*Q**H or C*Q\n*          or P*C or P**H*C or C*P or C*P**H.\n*\n*  LDC     (input) INTEGER\n*          The leading dimension of the array C. LDC >= max(1,M).\n*\n*  WORK    (workspace/output) COMPLEX array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK.\n*          If SIDE = 'L', LWORK >= max(1,N);\n*          if SIDE = 'R', LWORK >= max(1,M);\n*          if N = 0 or M = 0, LWORK >= 1.\n*          For optimum performance LWORK >= max(1,N*NB) if SIDE = 'L',\n*          and LWORK >= max(1,M*NB) if SIDE = 'R', where NB is the\n*          optimal blocksize. (NB = 0 if M = 0 or N = 0.)\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*\n\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      LOGICAL            APPLYQ, LEFT, LQUERY, NOTRAN\n      CHARACTER          TRANST\n      INTEGER            I1, I2, IINFO, LWKOPT, MI, NB, NI, NQ, NW\n*     ..\n*     .. External Functions ..\n      LOGICAL            LSAME\n      INTEGER            ILAENV\n      EXTERNAL           ILAENV, LSAME\n*     ..\n*     .. External Subroutines ..\n      EXTERNAL           CUNMLQ, CUNMQR, XERBLA\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          MAX, MIN\n*     ..\n\n");
    return Qnil;
  }
  if (argc != 9)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 9)", argc);
  rb_vect = argv[0];
  rb_side = argv[1];
  rb_trans = argv[2];
  rb_m = argv[3];
  rb_k = argv[4];
  rb_a = argv[5];
  rb_tau = argv[6];
  rb_c = argv[7];
  rb_lwork = argv[8];

  k = NUM2INT(rb_k);
  lwork = NUM2INT(rb_lwork);
  side = StringValueCStr(rb_side)[0];
  m = NUM2INT(rb_m);
  vect = StringValueCStr(rb_vect)[0];
  if (!NA_IsNArray(rb_c))
    rb_raise(rb_eArgError, "c (8th argument) must be NArray");
  if (NA_RANK(rb_c) != 2)
    rb_raise(rb_eArgError, "rank of c (8th argument) must be %d", 2);
  n = NA_SHAPE1(rb_c);
  ldc = NA_SHAPE0(rb_c);
  if (NA_TYPE(rb_c) != NA_SCOMPLEX)
    rb_c = na_change_type(rb_c, NA_SCOMPLEX);
  c = NA_PTR_TYPE(rb_c, complex*);
  trans = StringValueCStr(rb_trans)[0];
  nq = lsame_(&side,"L") ? m : lsame_(&side,"R") ? n : 0;
  if (!NA_IsNArray(rb_tau))
    rb_raise(rb_eArgError, "tau (7th argument) must be NArray");
  if (NA_RANK(rb_tau) != 1)
    rb_raise(rb_eArgError, "rank of tau (7th argument) must be %d", 1);
  if (NA_SHAPE0(rb_tau) != (MIN(nq,k)))
    rb_raise(rb_eRuntimeError, "shape 0 of tau must be %d", MIN(nq,k));
  if (NA_TYPE(rb_tau) != NA_SCOMPLEX)
    rb_tau = na_change_type(rb_tau, NA_SCOMPLEX);
  tau = NA_PTR_TYPE(rb_tau, complex*);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (6th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (6th argument) must be %d", 2);
  if (NA_SHAPE1(rb_a) != (MIN(nq,k)))
    rb_raise(rb_eRuntimeError, "shape 1 of a must be %d", MIN(nq,k));
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SCOMPLEX)
    rb_a = na_change_type(rb_a, NA_SCOMPLEX);
  a = NA_PTR_TYPE(rb_a, complex*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rb_work = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rb_work, complex*);
  {
    int shape[2];
    shape[0] = ldc;
    shape[1] = n;
    rb_c_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  c_out__ = NA_PTR_TYPE(rb_c_out__, complex*);
  MEMCPY(c_out__, c, complex, NA_TOTAL(rb_c));
  rb_c = rb_c_out__;
  c = c_out__;

  cunmbr_(&vect, &side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(3, rb_work, rb_info, rb_c);
}

void
init_lapack_cunmbr(VALUE mLapack){
  rb_define_module_function(mLapack, "cunmbr", rb_cunmbr, -1);
}
