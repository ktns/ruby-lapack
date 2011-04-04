#include "rb_lapack.h"

extern VOID zunmtr_(char *side, char *uplo, char *trans, integer *m, integer *n, doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *c, integer *ldc, doublecomplex *work, integer *lwork, integer *info);

static VALUE
rb_zunmtr(int argc, VALUE *argv, VALUE self){
  VALUE rb_side;
  char side; 
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_trans;
  char trans; 
  VALUE rb_a;
  doublecomplex *a; 
  VALUE rb_tau;
  doublecomplex *tau; 
  VALUE rb_c;
  doublecomplex *c; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_work;
  doublecomplex *work; 
  VALUE rb_info;
  integer info; 
  VALUE rb_c_out__;
  doublecomplex *c_out__;

  integer lda;
  integer m;
  integer ldc;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  work, info, c = NumRu::Lapack.zunmtr( side, uplo, trans, a, tau, c, lwork)\n    or\n  NumRu::Lapack.zunmtr  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZUNMTR( SIDE, UPLO, TRANS, M, N, A, LDA, TAU, C, LDC, WORK, LWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  ZUNMTR overwrites the general complex M-by-N matrix C with\n*\n*                  SIDE = 'L'     SIDE = 'R'\n*  TRANS = 'N':      Q * C          C * Q\n*  TRANS = 'C':      Q**H * C       C * Q**H\n*\n*  where Q is a complex unitary matrix of order nq, with nq = m if\n*  SIDE = 'L' and nq = n if SIDE = 'R'. Q is defined as the product of\n*  nq-1 elementary reflectors, as returned by ZHETRD:\n*\n*  if UPLO = 'U', Q = H(nq-1) . . . H(2) H(1);\n*\n*  if UPLO = 'L', Q = H(1) H(2) . . . H(nq-1).\n*\n\n*  Arguments\n*  =========\n*\n*  SIDE    (input) CHARACTER*1\n*          = 'L': apply Q or Q**H from the Left;\n*          = 'R': apply Q or Q**H from the Right.\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U': Upper triangle of A contains elementary reflectors\n*                 from ZHETRD;\n*          = 'L': Lower triangle of A contains elementary reflectors\n*                 from ZHETRD.\n*\n*  TRANS   (input) CHARACTER*1\n*          = 'N':  No transpose, apply Q;\n*          = 'C':  Conjugate transpose, apply Q**H.\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix C. M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix C. N >= 0.\n*\n*  A       (input) COMPLEX*16 array, dimension\n*                               (LDA,M) if SIDE = 'L'\n*                               (LDA,N) if SIDE = 'R'\n*          The vectors which define the elementary reflectors, as\n*          returned by ZHETRD.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.\n*          LDA >= max(1,M) if SIDE = 'L'; LDA >= max(1,N) if SIDE = 'R'.\n*\n*  TAU     (input) COMPLEX*16 array, dimension\n*                               (M-1) if SIDE = 'L'\n*                               (N-1) if SIDE = 'R'\n*          TAU(i) must contain the scalar factor of the elementary\n*          reflector H(i), as returned by ZHETRD.\n*\n*  C       (input/output) COMPLEX*16 array, dimension (LDC,N)\n*          On entry, the M-by-N matrix C.\n*          On exit, C is overwritten by Q*C or Q**H*C or C*Q**H or C*Q.\n*\n*  LDC     (input) INTEGER\n*          The leading dimension of the array C. LDC >= max(1,M).\n*\n*  WORK    (workspace/output) COMPLEX*16 array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK.\n*          If SIDE = 'L', LWORK >= max(1,N);\n*          if SIDE = 'R', LWORK >= max(1,M).\n*          For optimum performance LWORK >= N*NB if SIDE = 'L', and\n*          LWORK >=M*NB if SIDE = 'R', where NB is the optimal\n*          blocksize.\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*\n\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      LOGICAL            LEFT, LQUERY, UPPER\n      INTEGER            I1, I2, IINFO, LWKOPT, MI, NB, NI, NQ, NW\n*     ..\n*     .. External Functions ..\n      LOGICAL            LSAME\n      INTEGER            ILAENV\n      EXTERNAL           LSAME, ILAENV\n*     ..\n*     .. External Subroutines ..\n      EXTERNAL           XERBLA, ZUNMQL, ZUNMQR\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          MAX\n*     ..\n\n");
    return Qnil;
  }
  if (argc != 7)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 7)", argc);
  rb_side = argv[0];
  rb_uplo = argv[1];
  rb_trans = argv[2];
  rb_a = argv[3];
  rb_tau = argv[4];
  rb_c = argv[5];
  rb_lwork = argv[6];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (4th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (4th argument) must be %d", 2);
  m = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_DCOMPLEX)
    rb_a = na_change_type(rb_a, NA_DCOMPLEX);
  a = NA_PTR_TYPE(rb_a, doublecomplex*);
  side = StringValueCStr(rb_side)[0];
  if (!NA_IsNArray(rb_c))
    rb_raise(rb_eArgError, "c (6th argument) must be NArray");
  if (NA_RANK(rb_c) != 2)
    rb_raise(rb_eArgError, "rank of c (6th argument) must be %d", 2);
  n = NA_SHAPE1(rb_c);
  ldc = NA_SHAPE0(rb_c);
  if (NA_TYPE(rb_c) != NA_DCOMPLEX)
    rb_c = na_change_type(rb_c, NA_DCOMPLEX);
  c = NA_PTR_TYPE(rb_c, doublecomplex*);
  trans = StringValueCStr(rb_trans)[0];
  uplo = StringValueCStr(rb_uplo)[0];
  lwork = NUM2INT(rb_lwork);
  if (!NA_IsNArray(rb_tau))
    rb_raise(rb_eArgError, "tau (5th argument) must be NArray");
  if (NA_RANK(rb_tau) != 1)
    rb_raise(rb_eArgError, "rank of tau (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_tau) != (m-1))
    rb_raise(rb_eRuntimeError, "shape 0 of tau must be %d", m-1);
  if (NA_TYPE(rb_tau) != NA_DCOMPLEX)
    rb_tau = na_change_type(rb_tau, NA_DCOMPLEX);
  tau = NA_PTR_TYPE(rb_tau, doublecomplex*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rb_work = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rb_work, doublecomplex*);
  {
    int shape[2];
    shape[0] = ldc;
    shape[1] = n;
    rb_c_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  c_out__ = NA_PTR_TYPE(rb_c_out__, doublecomplex*);
  MEMCPY(c_out__, c, doublecomplex, NA_TOTAL(rb_c));
  rb_c = rb_c_out__;
  c = c_out__;

  zunmtr_(&side, &uplo, &trans, &m, &n, a, &lda, tau, c, &ldc, work, &lwork, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(3, rb_work, rb_info, rb_c);
}

void
init_lapack_zunmtr(VALUE mLapack){
  rb_define_module_function(mLapack, "zunmtr", rb_zunmtr, -1);
}
