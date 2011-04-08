#include "rb_lapack.h"

extern VOID dsbtrd_(char *vect, char *uplo, integer *n, integer *kd, doublereal *ab, integer *ldab, doublereal *d, doublereal *e, doublereal *q, integer *ldq, doublereal *work, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dsbtrd(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_vect;
  char vect; 
  VALUE rblapack_uplo;
  char uplo; 
  VALUE rblapack_kd;
  integer kd; 
  VALUE rblapack_ab;
  doublereal *ab; 
  VALUE rblapack_q;
  doublereal *q; 
  VALUE rblapack_d;
  doublereal *d; 
  VALUE rblapack_e;
  doublereal *e; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_ab_out__;
  doublereal *ab_out__;
  VALUE rblapack_q_out__;
  doublereal *q_out__;
  doublereal *work;

  integer ldab;
  integer n;
  integer ldq;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  d, e, info, ab, q = NumRu::Lapack.dsbtrd( vect, uplo, kd, ab, q, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DSBTRD( VECT, UPLO, N, KD, AB, LDAB, D, E, Q, LDQ, WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  DSBTRD reduces a real symmetric band matrix A to symmetric\n*  tridiagonal form T by an orthogonal similarity transformation:\n*  Q**T * A * Q = T.\n*\n\n*  Arguments\n*  =========\n*\n*  VECT    (input) CHARACTER*1\n*          = 'N':  do not form Q;\n*          = 'V':  form Q;\n*          = 'U':  update a matrix X, by forming X*Q.\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U':  Upper triangle of A is stored;\n*          = 'L':  Lower triangle of A is stored.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  KD      (input) INTEGER\n*          The number of superdiagonals of the matrix A if UPLO = 'U',\n*          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.\n*\n*  AB      (input/output) DOUBLE PRECISION array, dimension (LDAB,N)\n*          On entry, the upper or lower triangle of the symmetric band\n*          matrix A, stored in the first KD+1 rows of the array.  The\n*          j-th column of A is stored in the j-th column of the array AB\n*          as follows:\n*          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;\n*          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).\n*          On exit, the diagonal elements of AB are overwritten by the\n*          diagonal elements of the tridiagonal matrix T; if KD > 0, the\n*          elements on the first superdiagonal (if UPLO = 'U') or the\n*          first subdiagonal (if UPLO = 'L') are overwritten by the\n*          off-diagonal elements of T; the rest of AB is overwritten by\n*          values generated during the reduction.\n*\n*  LDAB    (input) INTEGER\n*          The leading dimension of the array AB.  LDAB >= KD+1.\n*\n*  D       (output) DOUBLE PRECISION array, dimension (N)\n*          The diagonal elements of the tridiagonal matrix T.\n*\n*  E       (output) DOUBLE PRECISION array, dimension (N-1)\n*          The off-diagonal elements of the tridiagonal matrix T:\n*          E(i) = T(i,i+1) if UPLO = 'U'; E(i) = T(i+1,i) if UPLO = 'L'.\n*\n*  Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,N)\n*          On entry, if VECT = 'U', then Q must contain an N-by-N\n*          matrix X; if VECT = 'N' or 'V', then Q need not be set.\n*\n*          On exit:\n*          if VECT = 'V', Q contains the N-by-N orthogonal matrix Q;\n*          if VECT = 'U', Q contains the product X*Q;\n*          if VECT = 'N', the array Q is not referenced.\n*\n*  LDQ     (input) INTEGER\n*          The leading dimension of the array Q.\n*          LDQ >= 1, and LDQ >= N if VECT = 'V' or 'U'.\n*\n*  WORK    (workspace) DOUBLE PRECISION array, dimension (N)\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*\n\n*  Further Details\n*  ===============\n*\n*  Modified by Linda Kaufman, Bell Labs.\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  d, e, info, ab, q = NumRu::Lapack.dsbtrd( vect, uplo, kd, ab, q, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rblapack_vect = argv[0];
  rblapack_uplo = argv[1];
  rblapack_kd = argv[2];
  rblapack_ab = argv[3];
  rblapack_q = argv[4];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_ab))
    rb_raise(rb_eArgError, "ab (4th argument) must be NArray");
  if (NA_RANK(rblapack_ab) != 2)
    rb_raise(rb_eArgError, "rank of ab (4th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_ab);
  ldab = NA_SHAPE0(rblapack_ab);
  if (NA_TYPE(rblapack_ab) != NA_DFLOAT)
    rblapack_ab = na_change_type(rblapack_ab, NA_DFLOAT);
  ab = NA_PTR_TYPE(rblapack_ab, doublereal*);
  vect = StringValueCStr(rblapack_vect)[0];
  kd = NUM2INT(rblapack_kd);
  if (!NA_IsNArray(rblapack_q))
    rb_raise(rb_eArgError, "q (5th argument) must be NArray");
  if (NA_RANK(rblapack_q) != 2)
    rb_raise(rb_eArgError, "rank of q (5th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_q) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of q must be the same as shape 1 of ab");
  ldq = NA_SHAPE0(rblapack_q);
  if (NA_TYPE(rblapack_q) != NA_DFLOAT)
    rblapack_q = na_change_type(rblapack_q, NA_DFLOAT);
  q = NA_PTR_TYPE(rblapack_q, doublereal*);
  uplo = StringValueCStr(rblapack_uplo)[0];
  {
    int shape[1];
    shape[0] = n;
    rblapack_d = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  d = NA_PTR_TYPE(rblapack_d, doublereal*);
  {
    int shape[1];
    shape[0] = n-1;
    rblapack_e = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  e = NA_PTR_TYPE(rblapack_e, doublereal*);
  {
    int shape[2];
    shape[0] = ldab;
    shape[1] = n;
    rblapack_ab_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  ab_out__ = NA_PTR_TYPE(rblapack_ab_out__, doublereal*);
  MEMCPY(ab_out__, ab, doublereal, NA_TOTAL(rblapack_ab));
  rblapack_ab = rblapack_ab_out__;
  ab = ab_out__;
  {
    int shape[2];
    shape[0] = ldq;
    shape[1] = n;
    rblapack_q_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  q_out__ = NA_PTR_TYPE(rblapack_q_out__, doublereal*);
  MEMCPY(q_out__, q, doublereal, NA_TOTAL(rblapack_q));
  rblapack_q = rblapack_q_out__;
  q = q_out__;
  work = ALLOC_N(doublereal, (n));

  dsbtrd_(&vect, &uplo, &n, &kd, ab, &ldab, d, e, q, &ldq, work, &info);

  free(work);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(5, rblapack_d, rblapack_e, rblapack_info, rblapack_ab, rblapack_q);
}

void
init_lapack_dsbtrd(VALUE mLapack){
  rb_define_module_function(mLapack, "dsbtrd", rblapack_dsbtrd, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
