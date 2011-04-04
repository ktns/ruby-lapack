#include "rb_lapack.h"

extern VOID ctrexc_(char *compq, integer *n, complex *t, integer *ldt, complex *q, integer *ldq, integer *ifst, integer *ilst, integer *info);

static VALUE
rb_ctrexc(int argc, VALUE *argv, VALUE self){
  VALUE rb_compq;
  char compq; 
  VALUE rb_t;
  complex *t; 
  VALUE rb_q;
  complex *q; 
  VALUE rb_ifst;
  integer ifst; 
  VALUE rb_ilst;
  integer ilst; 
  VALUE rb_info;
  integer info; 
  VALUE rb_t_out__;
  complex *t_out__;
  VALUE rb_q_out__;
  complex *q_out__;

  integer ldt;
  integer n;
  integer ldq;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, t, q = NumRu::Lapack.ctrexc( compq, t, q, ifst, ilst)\n    or\n  NumRu::Lapack.ctrexc  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE CTREXC( COMPQ, N, T, LDT, Q, LDQ, IFST, ILST, INFO )\n\n*  Purpose\n*  =======\n*\n*  CTREXC reorders the Schur factorization of a complex matrix\n*  A = Q*T*Q**H, so that the diagonal element of T with row index IFST\n*  is moved to row ILST.\n*\n*  The Schur form T is reordered by a unitary similarity transformation\n*  Z**H*T*Z, and optionally the matrix Q of Schur vectors is updated by\n*  postmultplying it with Z.\n*\n\n*  Arguments\n*  =========\n*\n*  COMPQ   (input) CHARACTER*1\n*          = 'V':  update the matrix Q of Schur vectors;\n*          = 'N':  do not update Q.\n*\n*  N       (input) INTEGER\n*          The order of the matrix T. N >= 0.\n*\n*  T       (input/output) COMPLEX array, dimension (LDT,N)\n*          On entry, the upper triangular matrix T.\n*          On exit, the reordered upper triangular matrix.\n*\n*  LDT     (input) INTEGER\n*          The leading dimension of the array T. LDT >= max(1,N).\n*\n*  Q       (input/output) COMPLEX array, dimension (LDQ,N)\n*          On entry, if COMPQ = 'V', the matrix Q of Schur vectors.\n*          On exit, if COMPQ = 'V', Q has been postmultiplied by the\n*          unitary transformation matrix Z which reorders T.\n*          If COMPQ = 'N', Q is not referenced.\n*\n*  LDQ     (input) INTEGER\n*          The leading dimension of the array Q.  LDQ >= max(1,N).\n*\n*  IFST    (input) INTEGER\n*  ILST    (input) INTEGER\n*          Specify the reordering of the diagonal elements of T:\n*          The element with row index IFST is moved to row ILST by a\n*          sequence of transpositions between adjacent elements.\n*          1 <= IFST <= N; 1 <= ILST <= N.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*\n\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      LOGICAL            WANTQ\n      INTEGER            K, M1, M2, M3\n      REAL               CS\n      COMPLEX            SN, T11, T22, TEMP\n*     ..\n*     .. External Functions ..\n      LOGICAL            LSAME\n      EXTERNAL           LSAME\n*     ..\n*     .. External Subroutines ..\n      EXTERNAL           CLARTG, CROT, XERBLA\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          CONJG, MAX\n*     ..\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_compq = argv[0];
  rb_t = argv[1];
  rb_q = argv[2];
  rb_ifst = argv[3];
  rb_ilst = argv[4];

  compq = StringValueCStr(rb_compq)[0];
  ilst = NUM2INT(rb_ilst);
  if (!NA_IsNArray(rb_q))
    rb_raise(rb_eArgError, "q (3th argument) must be NArray");
  if (NA_RANK(rb_q) != 2)
    rb_raise(rb_eArgError, "rank of q (3th argument) must be %d", 2);
  n = NA_SHAPE1(rb_q);
  ldq = NA_SHAPE0(rb_q);
  if (NA_TYPE(rb_q) != NA_SCOMPLEX)
    rb_q = na_change_type(rb_q, NA_SCOMPLEX);
  q = NA_PTR_TYPE(rb_q, complex*);
  if (!NA_IsNArray(rb_t))
    rb_raise(rb_eArgError, "t (2th argument) must be NArray");
  if (NA_RANK(rb_t) != 2)
    rb_raise(rb_eArgError, "rank of t (2th argument) must be %d", 2);
  if (NA_SHAPE1(rb_t) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of t must be the same as shape 1 of q");
  ldt = NA_SHAPE0(rb_t);
  if (NA_TYPE(rb_t) != NA_SCOMPLEX)
    rb_t = na_change_type(rb_t, NA_SCOMPLEX);
  t = NA_PTR_TYPE(rb_t, complex*);
  ifst = NUM2INT(rb_ifst);
  {
    int shape[2];
    shape[0] = ldt;
    shape[1] = n;
    rb_t_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  t_out__ = NA_PTR_TYPE(rb_t_out__, complex*);
  MEMCPY(t_out__, t, complex, NA_TOTAL(rb_t));
  rb_t = rb_t_out__;
  t = t_out__;
  {
    int shape[2];
    shape[0] = ldq;
    shape[1] = n;
    rb_q_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  q_out__ = NA_PTR_TYPE(rb_q_out__, complex*);
  MEMCPY(q_out__, q, complex, NA_TOTAL(rb_q));
  rb_q = rb_q_out__;
  q = q_out__;

  ctrexc_(&compq, &n, t, &ldt, q, &ldq, &ifst, &ilst, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(3, rb_info, rb_t, rb_q);
}

void
init_lapack_ctrexc(VALUE mLapack){
  rb_define_module_function(mLapack, "ctrexc", rb_ctrexc, -1);
}
