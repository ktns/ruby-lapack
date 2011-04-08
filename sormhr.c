#include "rb_lapack.h"

extern VOID sormhr_(char *side, char *trans, integer *m, integer *n, integer *ilo, integer *ihi, real *a, integer *lda, real *tau, real *c, integer *ldc, real *work, integer *lwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_sormhr(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_side;
  char side; 
  VALUE rblapack_trans;
  char trans; 
  VALUE rblapack_ilo;
  integer ilo; 
  VALUE rblapack_ihi;
  integer ihi; 
  VALUE rblapack_a;
  real *a; 
  VALUE rblapack_tau;
  real *tau; 
  VALUE rblapack_c;
  real *c; 
  VALUE rblapack_lwork;
  integer lwork; 
  VALUE rblapack_work;
  real *work; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_c_out__;
  real *c_out__;

  integer lda;
  integer m;
  integer ldc;
  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  work, info, c = NumRu::Lapack.sormhr( side, trans, ilo, ihi, a, tau, c, lwork, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE SORMHR( SIDE, TRANS, M, N, ILO, IHI, A, LDA, TAU, C, LDC, WORK, LWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  SORMHR overwrites the general real M-by-N matrix C with\n*\n*                  SIDE = 'L'     SIDE = 'R'\n*  TRANS = 'N':      Q * C          C * Q\n*  TRANS = 'T':      Q**T * C       C * Q**T\n*\n*  where Q is a real orthogonal matrix of order nq, with nq = m if\n*  SIDE = 'L' and nq = n if SIDE = 'R'. Q is defined as the product of\n*  IHI-ILO elementary reflectors, as returned by SGEHRD:\n*\n*  Q = H(ilo) H(ilo+1) . . . H(ihi-1).\n*\n\n*  Arguments\n*  =========\n*\n*  SIDE    (input) CHARACTER*1\n*          = 'L': apply Q or Q**T from the Left;\n*          = 'R': apply Q or Q**T from the Right.\n*\n*  TRANS   (input) CHARACTER*1\n*          = 'N':  No transpose, apply Q;\n*          = 'T':  Transpose, apply Q**T.\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix C. M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix C. N >= 0.\n*\n*  ILO     (input) INTEGER\n*  IHI     (input) INTEGER\n*          ILO and IHI must have the same values as in the previous call\n*          of SGEHRD. Q is equal to the unit matrix except in the\n*          submatrix Q(ilo+1:ihi,ilo+1:ihi).\n*          If SIDE = 'L', then 1 <= ILO <= IHI <= M, if M > 0, and\n*          ILO = 1 and IHI = 0, if M = 0;\n*          if SIDE = 'R', then 1 <= ILO <= IHI <= N, if N > 0, and\n*          ILO = 1 and IHI = 0, if N = 0.\n*\n*  A       (input) REAL array, dimension\n*                               (LDA,M) if SIDE = 'L'\n*                               (LDA,N) if SIDE = 'R'\n*          The vectors which define the elementary reflectors, as\n*          returned by SGEHRD.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.\n*          LDA >= max(1,M) if SIDE = 'L'; LDA >= max(1,N) if SIDE = 'R'.\n*\n*  TAU     (input) REAL array, dimension\n*                               (M-1) if SIDE = 'L'\n*                               (N-1) if SIDE = 'R'\n*          TAU(i) must contain the scalar factor of the elementary\n*          reflector H(i), as returned by SGEHRD.\n*\n*  C       (input/output) REAL array, dimension (LDC,N)\n*          On entry, the M-by-N matrix C.\n*          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.\n*\n*  LDC     (input) INTEGER\n*          The leading dimension of the array C. LDC >= max(1,M).\n*\n*  WORK    (workspace/output) REAL array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK.\n*          If SIDE = 'L', LWORK >= max(1,N);\n*          if SIDE = 'R', LWORK >= max(1,M).\n*          For optimum performance LWORK >= N*NB if SIDE = 'L', and\n*          LWORK >= M*NB if SIDE = 'R', where NB is the optimal\n*          blocksize.\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*\n\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      LOGICAL            LEFT, LQUERY\n      INTEGER            I1, I2, IINFO, LWKOPT, MI, NB, NH, NI, NQ, NW\n*     ..\n*     .. External Functions ..\n      LOGICAL            LSAME\n      INTEGER            ILAENV\n      EXTERNAL           ILAENV, LSAME\n*     ..\n*     .. External Subroutines ..\n      EXTERNAL           SORMQR, XERBLA\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          MAX, MIN\n*     ..\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  work, info, c = NumRu::Lapack.sormhr( side, trans, ilo, ihi, a, tau, c, lwork, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 8)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 8)", argc);
  rblapack_side = argv[0];
  rblapack_trans = argv[1];
  rblapack_ilo = argv[2];
  rblapack_ihi = argv[3];
  rblapack_a = argv[4];
  rblapack_tau = argv[5];
  rblapack_c = argv[6];
  rblapack_lwork = argv[7];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (5th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (5th argument) must be %d", 2);
  m = NA_SHAPE1(rblapack_a);
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_SFLOAT)
    rblapack_a = na_change_type(rblapack_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rblapack_a, real*);
  ilo = NUM2INT(rblapack_ilo);
  side = StringValueCStr(rblapack_side)[0];
  if (!NA_IsNArray(rblapack_c))
    rb_raise(rb_eArgError, "c (7th argument) must be NArray");
  if (NA_RANK(rblapack_c) != 2)
    rb_raise(rb_eArgError, "rank of c (7th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_c);
  ldc = NA_SHAPE0(rblapack_c);
  if (NA_TYPE(rblapack_c) != NA_SFLOAT)
    rblapack_c = na_change_type(rblapack_c, NA_SFLOAT);
  c = NA_PTR_TYPE(rblapack_c, real*);
  ihi = NUM2INT(rblapack_ihi);
  trans = StringValueCStr(rblapack_trans)[0];
  lwork = NUM2INT(rblapack_lwork);
  if (!NA_IsNArray(rblapack_tau))
    rb_raise(rb_eArgError, "tau (6th argument) must be NArray");
  if (NA_RANK(rblapack_tau) != 1)
    rb_raise(rb_eArgError, "rank of tau (6th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_tau) != (m-1))
    rb_raise(rb_eRuntimeError, "shape 0 of tau must be %d", m-1);
  if (NA_TYPE(rblapack_tau) != NA_SFLOAT)
    rblapack_tau = na_change_type(rblapack_tau, NA_SFLOAT);
  tau = NA_PTR_TYPE(rblapack_tau, real*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rblapack_work = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rblapack_work, real*);
  {
    int shape[2];
    shape[0] = ldc;
    shape[1] = n;
    rblapack_c_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  c_out__ = NA_PTR_TYPE(rblapack_c_out__, real*);
  MEMCPY(c_out__, c, real, NA_TOTAL(rblapack_c));
  rblapack_c = rblapack_c_out__;
  c = c_out__;

  sormhr_(&side, &trans, &m, &n, &ilo, &ihi, a, &lda, tau, c, &ldc, work, &lwork, &info);

  rblapack_info = INT2NUM(info);
  return rb_ary_new3(3, rblapack_work, rblapack_info, rblapack_c);
}

void
init_lapack_sormhr(VALUE mLapack){
  rb_define_module_function(mLapack, "sormhr", rblapack_sormhr, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
