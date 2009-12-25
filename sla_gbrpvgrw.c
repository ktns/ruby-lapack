#include "rb_lapack.h"

extern VOID sla_gbrpvgrw_(real *__out__, integer *n, integer *kl, integer *ku, integer *ncols, real *ab, integer *ldab, real *afb, integer *ldafb);
static VALUE
rb_sla_gbrpvgrw(int argc, VALUE *argv, VALUE self){
  VALUE rb_kl;
  integer kl; 
  VALUE rb_ku;
  integer ku; 
  VALUE rb_ncols;
  integer ncols; 
  VALUE rb_ab;
  real *ab; 
  VALUE rb_afb;
  real *afb; 
  VALUE rb___out__;
  real __out__; 

  integer ldab;
  integer n;
  integer ldafb;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.sla_gbrpvgrw( kl, ku, ncols, ab, afb)\n    or\n  NumRu::Lapack.sla_gbrpvgrw  # print help\n\n\nFORTRAN MANUAL\n      REAL FUNCTION SLA_GBRPVGRW( N, KL, KU, NCOLS, AB, LDAB, AFB, LDAFB )\n\n*  Purpose\n*  =======\n*\n*  SLA_GBRPVGRW computes the reciprocal pivot growth factor\n*  norm(A)/norm(U). The \"max absolute element\" norm is used. If this is\n*  much less than 1, the stability of the LU factorization of the\n*  (equilibrated) matrix A could be poor. This also means that the\n*  solution X, estimated condition numbers, and error bounds could be\n*  unreliable.\n*\n\n*  Arguments\n*  =========\n*\n*     N       (input) INTEGER\n*     The number of linear equations, i.e., the order of the\n*     matrix A.  N >= 0.\n*\n*     KL      (input) INTEGER\n*     The number of subdiagonals within the band of A.  KL >= 0.\n*\n*     KU      (input) INTEGER\n*     The number of superdiagonals within the band of A.  KU >= 0.\n*\n*     NCOLS   (input) INTEGER\n*     The number of columns of the matrix A.  NCOLS >= 0.\n*\n*     AB      (input) REAL array, dimension (LDAB,N)\n*     On entry, the matrix A in band storage, in rows 1 to KL+KU+1.\n*     The j-th column of A is stored in the j-th column of the\n*     array AB as follows:\n*     AB(KU+1+i-j,j) = A(i,j) for max(1,j-KU)<=i<=min(N,j+kl)\n*\n*     LDAB    (input) INTEGER\n*     The leading dimension of the array AB.  LDAB >= KL+KU+1.\n*\n*     AFB     (input) REAL array, dimension (LDAFB,N)\n*     Details of the LU factorization of the band matrix A, as\n*     computed by SGBTRF.  U is stored as an upper triangular\n*     band matrix with KL+KU superdiagonals in rows 1 to KL+KU+1,\n*     and the multipliers used during the factorization are stored\n*     in rows KL+KU+2 to 2*KL+KU+1.\n*\n*     LDAFB   (input) INTEGER\n*     The leading dimension of the array AFB.  LDAFB >= 2*KL+KU+1.\n*\n\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      INTEGER            I, J, KD\n      REAL               AMAX, UMAX, RPVGRW\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          ABS, MAX, MIN\n*     ..\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_kl = argv[0];
  rb_ku = argv[1];
  rb_ncols = argv[2];
  rb_ab = argv[3];
  rb_afb = argv[4];

  kl = NUM2INT(rb_kl);
  ku = NUM2INT(rb_ku);
  ncols = NUM2INT(rb_ncols);
  if (!NA_IsNArray(rb_ab))
    rb_raise(rb_eArgError, "ab (4th argument) must be NArray");
  if (NA_RANK(rb_ab) != 2)
    rb_raise(rb_eArgError, "rank of ab (4th argument) must be %d", 2);
  ldab = NA_SHAPE0(rb_ab);
  n = NA_SHAPE1(rb_ab);
  if (NA_TYPE(rb_ab) != NA_SFLOAT)
    rb_ab = na_change_type(rb_ab, NA_SFLOAT);
  ab = NA_PTR_TYPE(rb_ab, real*);
  if (!NA_IsNArray(rb_afb))
    rb_raise(rb_eArgError, "afb (5th argument) must be NArray");
  if (NA_RANK(rb_afb) != 2)
    rb_raise(rb_eArgError, "rank of afb (5th argument) must be %d", 2);
  ldafb = NA_SHAPE0(rb_afb);
  if (NA_SHAPE1(rb_afb) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of afb must be the same as shape 1 of ab");
  if (NA_TYPE(rb_afb) != NA_SFLOAT)
    rb_afb = na_change_type(rb_afb, NA_SFLOAT);
  afb = NA_PTR_TYPE(rb_afb, real*);

  sla_gbrpvgrw_(&__out__, &n, &kl, &ku, &ncols, ab, &ldab, afb, &ldafb);

  rb___out__ = rb_float_new((double)__out__);
  return rb___out__;
}

void
init_lapack_sla_gbrpvgrw(VALUE mLapack){
  rb_define_module_function(mLapack, "sla_gbrpvgrw", rb_sla_gbrpvgrw, -1);
}
