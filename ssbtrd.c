#include "rb_lapack.h"

static VALUE
rb_ssbtrd(int argc, VALUE *argv, VALUE self){
  VALUE rb_vect;
  char vect; 
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_kd;
  integer kd; 
  VALUE rb_ab;
  real *ab; 
  VALUE rb_q;
  real *q; 
  VALUE rb_d;
  real *d; 
  VALUE rb_e;
  real *e; 
  VALUE rb_info;
  integer info; 
  VALUE rb_ab_out__;
  real *ab_out__;
  VALUE rb_q_out__;
  real *q_out__;
  real *work;

  integer ldab;
  integer n;
  integer ldq;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  d, e, info, ab, q = NumRu::Lapack.ssbtrd( vect, uplo, kd, ab, q)\n    or\n  NumRu::Lapack.ssbtrd  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SSBTRD( VECT, UPLO, N, KD, AB, LDAB, D, E, Q, LDQ, WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  SSBTRD reduces a real symmetric band matrix A to symmetric\n*  tridiagonal form T by an orthogonal similarity transformation:\n*  Q**T * A * Q = T.\n*\n\n*  Arguments\n*  =========\n*\n*  VECT    (input) CHARACTER*1\n*          = 'N':  do not form Q;\n*          = 'V':  form Q;\n*          = 'U':  update a matrix X, by forming X*Q.\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U':  Upper triangle of A is stored;\n*          = 'L':  Lower triangle of A is stored.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  KD      (input) INTEGER\n*          The number of superdiagonals of the matrix A if UPLO = 'U',\n*          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.\n*\n*  AB      (input/output) REAL array, dimension (LDAB,N)\n*          On entry, the upper or lower triangle of the symmetric band\n*          matrix A, stored in the first KD+1 rows of the array.  The\n*          j-th column of A is stored in the j-th column of the array AB\n*          as follows:\n*          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;\n*          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).\n*          On exit, the diagonal elements of AB are overwritten by the\n*          diagonal elements of the tridiagonal matrix T; if KD > 0, the\n*          elements on the first superdiagonal (if UPLO = 'U') or the\n*          first subdiagonal (if UPLO = 'L') are overwritten by the\n*          off-diagonal elements of T; the rest of AB is overwritten by\n*          values generated during the reduction.\n*\n*  LDAB    (input) INTEGER\n*          The leading dimension of the array AB.  LDAB >= KD+1.\n*\n*  D       (output) REAL array, dimension (N)\n*          The diagonal elements of the tridiagonal matrix T.\n*\n*  E       (output) REAL array, dimension (N-1)\n*          The off-diagonal elements of the tridiagonal matrix T:\n*          E(i) = T(i,i+1) if UPLO = 'U'; E(i) = T(i+1,i) if UPLO = 'L'.\n*\n*  Q       (input/output) REAL array, dimension (LDQ,N)\n*          On entry, if VECT = 'U', then Q must contain an N-by-N\n*          matrix X; if VECT = 'N' or 'V', then Q need not be set.\n*\n*          On exit:\n*          if VECT = 'V', Q contains the N-by-N orthogonal matrix Q;\n*          if VECT = 'U', Q contains the product X*Q;\n*          if VECT = 'N', the array Q is not referenced.\n*\n*  LDQ     (input) INTEGER\n*          The leading dimension of the array Q.\n*          LDQ >= 1, and LDQ >= N if VECT = 'V' or 'U'.\n*\n*  WORK    (workspace) REAL array, dimension (N)\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*\n\n*  Further Details\n*  ===============\n*\n*  Modified by Linda Kaufman, Bell Labs.\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_vect = argv[0];
  rb_uplo = argv[1];
  rb_kd = argv[2];
  rb_ab = argv[3];
  rb_q = argv[4];

  vect = StringValueCStr(rb_vect)[0];
  uplo = StringValueCStr(rb_uplo)[0];
  kd = NUM2INT(rb_kd);
  if (!NA_IsNArray(rb_ab))
    rb_raise(rb_eArgError, "ab (4th argument) must be NArray");
  if (NA_RANK(rb_ab) != 2)
    rb_raise(rb_eArgError, "rank of ab (4th argument) must be %d", 2);
  ldab = NA_SHAPE0(rb_ab);
  n = NA_SHAPE1(rb_ab);
  if (NA_TYPE(rb_ab) != NA_SFLOAT)
    rb_ab = na_change_type(rb_ab, NA_SFLOAT);
  ab = NA_PTR_TYPE(rb_ab, real*);
  if (!NA_IsNArray(rb_q))
    rb_raise(rb_eArgError, "q (5th argument) must be NArray");
  if (NA_RANK(rb_q) != 2)
    rb_raise(rb_eArgError, "rank of q (5th argument) must be %d", 2);
  ldq = NA_SHAPE0(rb_q);
  if (NA_SHAPE1(rb_q) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of q must be the same as shape 1 of ab");
  if (NA_TYPE(rb_q) != NA_SFLOAT)
    rb_q = na_change_type(rb_q, NA_SFLOAT);
  q = NA_PTR_TYPE(rb_q, real*);
  {
    int shape[1];
    shape[0] = n;
    rb_d = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  d = NA_PTR_TYPE(rb_d, real*);
  {
    int shape[1];
    shape[0] = n-1;
    rb_e = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  e = NA_PTR_TYPE(rb_e, real*);
  {
    int shape[2];
    shape[0] = ldab;
    shape[1] = n;
    rb_ab_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  ab_out__ = NA_PTR_TYPE(rb_ab_out__, real*);
  MEMCPY(ab_out__, ab, real, NA_TOTAL(rb_ab));
  rb_ab = rb_ab_out__;
  ab = ab_out__;
  {
    int shape[2];
    shape[0] = ldq;
    shape[1] = n;
    rb_q_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  q_out__ = NA_PTR_TYPE(rb_q_out__, real*);
  MEMCPY(q_out__, q, real, NA_TOTAL(rb_q));
  rb_q = rb_q_out__;
  q = q_out__;
  work = ALLOC_N(real, (n));

  ssbtrd_(&vect, &uplo, &n, &kd, ab, &ldab, d, e, q, &ldq, work, &info);

  free(work);
  rb_info = INT2NUM(info);
  return rb_ary_new3(5, rb_d, rb_e, rb_info, rb_ab, rb_q);
}

void
init_lapack_ssbtrd(VALUE mLapack){
  rb_define_module_function(mLapack, "ssbtrd", rb_ssbtrd, -1);
}
