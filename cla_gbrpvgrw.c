#include "rb_lapack.h"

extern real cla_gbrpvgrw_(integer *n, integer *kl, integer *ku, integer *ncols, complex *ab, integer *ldab, complex *afb, integer *ldafb);

static VALUE sHelp, sUsage;

static VALUE
rblapack_cla_gbrpvgrw(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_kl;
  integer kl; 
  VALUE rblapack_ku;
  integer ku; 
  VALUE rblapack_ncols;
  integer ncols; 
  VALUE rblapack_ab;
  complex *ab; 
  VALUE rblapack_afb;
  complex *afb; 
  VALUE rblapack___out__;
  real __out__; 

  integer ldab;
  integer n;
  integer ldafb;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.cla_gbrpvgrw( kl, ku, ncols, ab, afb, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      REAL FUNCTION CLA_GBRPVGRW( N, KL, KU, NCOLS, AB, LDAB, AFB, LDAFB )\n\n*  Purpose\n*  =======\n*\n*  CLA_GBRPVGRW computes the reciprocal pivot growth factor\n*  norm(A)/norm(U). The \"max absolute element\" norm is used. If this is\n*  much less than 1, the stability of the LU factorization of the\n*  (equilibrated) matrix A could be poor. This also means that the\n*  solution X, estimated condition numbers, and error bounds could be\n*  unreliable.\n*\n\n*  Arguments\n*  =========\n*\n*     N       (input) INTEGER\n*     The number of linear equations, i.e., the order of the\n*     matrix A.  N >= 0.\n*\n*     KL      (input) INTEGER\n*     The number of subdiagonals within the band of A.  KL >= 0.\n*\n*     KU      (input) INTEGER\n*     The number of superdiagonals within the band of A.  KU >= 0.\n*\n*     NCOLS   (input) INTEGER\n*     The number of columns of the matrix A.  NCOLS >= 0.\n*\n*     AB      (input) COMPLEX array, dimension (LDAB,N)\n*     On entry, the matrix A in band storage, in rows 1 to KL+KU+1.\n*     The j-th column of A is stored in the j-th column of the\n*     array AB as follows:\n*     AB(KU+1+i-j,j) = A(i,j) for max(1,j-KU)<=i<=min(N,j+kl)\n*\n*     LDAB    (input) INTEGER\n*     The leading dimension of the array AB.  LDAB >= KL+KU+1.\n*\n*     AFB     (input) COMPLEX array, dimension (LDAFB,N)\n*     Details of the LU factorization of the band matrix A, as\n*     computed by CGBTRF.  U is stored as an upper triangular\n*     band matrix with KL+KU superdiagonals in rows 1 to KL+KU+1,\n*     and the multipliers used during the factorization are stored\n*     in rows KL+KU+2 to 2*KL+KU+1.\n*\n*     LDAFB   (input) INTEGER\n*     The leading dimension of the array AFB.  LDAFB >= 2*KL+KU+1.\n*\n\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      INTEGER            I, J, KD\n      REAL               AMAX, UMAX, RPVGRW\n      COMPLEX            ZDUM\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          ABS, MAX, MIN, REAL, AIMAG\n*     ..\n*     .. Statement Functions ..\n      REAL               CABS1\n*     ..\n*     .. Statement Function Definitions ..\n      CABS1( ZDUM ) = ABS( REAL( ZDUM ) ) + ABS( AIMAG( ZDUM ) )\n*     ..\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.cla_gbrpvgrw( kl, ku, ncols, ab, afb, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rblapack_kl = argv[0];
  rblapack_ku = argv[1];
  rblapack_ncols = argv[2];
  rblapack_ab = argv[3];
  rblapack_afb = argv[4];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_ab))
    rb_raise(rb_eArgError, "ab (4th argument) must be NArray");
  if (NA_RANK(rblapack_ab) != 2)
    rb_raise(rb_eArgError, "rank of ab (4th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_ab);
  ldab = NA_SHAPE0(rblapack_ab);
  if (NA_TYPE(rblapack_ab) != NA_SCOMPLEX)
    rblapack_ab = na_change_type(rblapack_ab, NA_SCOMPLEX);
  ab = NA_PTR_TYPE(rblapack_ab, complex*);
  kl = NUM2INT(rblapack_kl);
  ncols = NUM2INT(rblapack_ncols);
  if (!NA_IsNArray(rblapack_afb))
    rb_raise(rb_eArgError, "afb (5th argument) must be NArray");
  if (NA_RANK(rblapack_afb) != 2)
    rb_raise(rb_eArgError, "rank of afb (5th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_afb) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of afb must be the same as shape 1 of ab");
  ldafb = NA_SHAPE0(rblapack_afb);
  if (NA_TYPE(rblapack_afb) != NA_SCOMPLEX)
    rblapack_afb = na_change_type(rblapack_afb, NA_SCOMPLEX);
  afb = NA_PTR_TYPE(rblapack_afb, complex*);
  ku = NUM2INT(rblapack_ku);

  __out__ = cla_gbrpvgrw_(&n, &kl, &ku, &ncols, ab, &ldab, afb, &ldafb);

  rblapack___out__ = rb_float_new((double)__out__);
  return rblapack___out__;
}

void
init_lapack_cla_gbrpvgrw(VALUE mLapack){
  rb_define_module_function(mLapack, "cla_gbrpvgrw", rblapack_cla_gbrpvgrw, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
