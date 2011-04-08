#include "rb_lapack.h"

extern VOID ztrexc_(char *compq, integer *n, doublecomplex *t, integer *ldt, doublecomplex *q, integer *ldq, integer *ifst, integer *ilst, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_ztrexc(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_compq;
  char compq; 
  VALUE rblapack_t;
  doublecomplex *t; 
  VALUE rblapack_q;
  doublecomplex *q; 
  VALUE rblapack_ifst;
  integer ifst; 
  VALUE rblapack_ilst;
  integer ilst; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_t_out__;
  doublecomplex *t_out__;
  VALUE rblapack_q_out__;
  doublecomplex *q_out__;

  integer ldt;
  integer n;
  integer ldq;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, t, q = NumRu::Lapack.ztrexc( compq, t, q, ifst, ilst, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZTREXC( COMPQ, N, T, LDT, Q, LDQ, IFST, ILST, INFO )\n\n*  Purpose\n*  =======\n*\n*  ZTREXC reorders the Schur factorization of a complex matrix\n*  A = Q*T*Q**H, so that the diagonal element of T with row index IFST\n*  is moved to row ILST.\n*\n*  The Schur form T is reordered by a unitary similarity transformation\n*  Z**H*T*Z, and optionally the matrix Q of Schur vectors is updated by\n*  postmultplying it with Z.\n*\n\n*  Arguments\n*  =========\n*\n*  COMPQ   (input) CHARACTER*1\n*          = 'V':  update the matrix Q of Schur vectors;\n*          = 'N':  do not update Q.\n*\n*  N       (input) INTEGER\n*          The order of the matrix T. N >= 0.\n*\n*  T       (input/output) COMPLEX*16 array, dimension (LDT,N)\n*          On entry, the upper triangular matrix T.\n*          On exit, the reordered upper triangular matrix.\n*\n*  LDT     (input) INTEGER\n*          The leading dimension of the array T. LDT >= max(1,N).\n*\n*  Q       (input/output) COMPLEX*16 array, dimension (LDQ,N)\n*          On entry, if COMPQ = 'V', the matrix Q of Schur vectors.\n*          On exit, if COMPQ = 'V', Q has been postmultiplied by the\n*          unitary transformation matrix Z which reorders T.\n*          If COMPQ = 'N', Q is not referenced.\n*\n*  LDQ     (input) INTEGER\n*          The leading dimension of the array Q.  LDQ >= max(1,N).\n*\n*  IFST    (input) INTEGER\n*  ILST    (input) INTEGER\n*          Specify the reordering of the diagonal elements of T:\n*          The element with row index IFST is moved to row ILST by a\n*          sequence of transpositions between adjacent elements.\n*          1 <= IFST <= N; 1 <= ILST <= N.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*\n\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      LOGICAL            WANTQ\n      INTEGER            K, M1, M2, M3\n      DOUBLE PRECISION   CS\n      COMPLEX*16         SN, T11, T22, TEMP\n*     ..\n*     .. External Functions ..\n      LOGICAL            LSAME\n      EXTERNAL           LSAME\n*     ..\n*     .. External Subroutines ..\n      EXTERNAL           XERBLA, ZLARTG, ZROT\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          DCONJG, MAX\n*     ..\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, t, q = NumRu::Lapack.ztrexc( compq, t, q, ifst, ilst, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rblapack_compq = argv[0];
  rblapack_t = argv[1];
  rblapack_q = argv[2];
  rblapack_ifst = argv[3];
  rblapack_ilst = argv[4];
  if (rb_options != Qnil) {
  }

  compq = StringValueCStr(rblapack_compq)[0];
  ilst = NUM2INT(rblapack_ilst);
  if (!NA_IsNArray(rblapack_q))
    rb_raise(rb_eArgError, "q (3th argument) must be NArray");
  if (NA_RANK(rblapack_q) != 2)
    rb_raise(rb_eArgError, "rank of q (3th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_q);
  ldq = NA_SHAPE0(rblapack_q);
  if (NA_TYPE(rblapack_q) != NA_DCOMPLEX)
    rblapack_q = na_change_type(rblapack_q, NA_DCOMPLEX);
  q = NA_PTR_TYPE(rblapack_q, doublecomplex*);
  if (!NA_IsNArray(rblapack_t))
    rb_raise(rb_eArgError, "t (2th argument) must be NArray");
  if (NA_RANK(rblapack_t) != 2)
    rb_raise(rb_eArgError, "rank of t (2th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_t) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of t must be the same as shape 1 of q");
  ldt = NA_SHAPE0(rblapack_t);
  if (NA_TYPE(rblapack_t) != NA_DCOMPLEX)
    rblapack_t = na_change_type(rblapack_t, NA_DCOMPLEX);
  t = NA_PTR_TYPE(rblapack_t, doublecomplex*);
  ifst = NUM2INT(rblapack_ifst);
  {
    int shape[2];
    shape[0] = ldt;
    shape[1] = n;
    rblapack_t_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  t_out__ = NA_PTR_TYPE(rblapack_t_out__, doublecomplex*);
  MEMCPY(t_out__, t, doublecomplex, NA_TOTAL(rblapack_t));
  rblapack_t = rblapack_t_out__;
  t = t_out__;
  {
    int shape[2];
    shape[0] = ldq;
    shape[1] = n;
    rblapack_q_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  q_out__ = NA_PTR_TYPE(rblapack_q_out__, doublecomplex*);
  MEMCPY(q_out__, q, doublecomplex, NA_TOTAL(rblapack_q));
  rblapack_q = rblapack_q_out__;
  q = q_out__;

  ztrexc_(&compq, &n, t, &ldt, q, &ldq, &ifst, &ilst, &info);

  rblapack_info = INT2NUM(info);
  return rb_ary_new3(3, rblapack_info, rblapack_t, rblapack_q);
}

void
init_lapack_ztrexc(VALUE mLapack){
  rb_define_module_function(mLapack, "ztrexc", rblapack_ztrexc, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
