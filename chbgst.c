#include "rb_lapack.h"

extern VOID chbgst_(char *vect, char *uplo, integer *n, integer *ka, integer *kb, complex *ab, integer *ldab, complex *bb, integer *ldbb, complex *x, integer *ldx, complex *work, real *rwork, integer *info);

static VALUE
rb_chbgst(int argc, VALUE *argv, VALUE self){
  VALUE rb_vect;
  char vect; 
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_ka;
  integer ka; 
  VALUE rb_kb;
  integer kb; 
  VALUE rb_ab;
  complex *ab; 
  VALUE rb_bb;
  complex *bb; 
  VALUE rb_x;
  complex *x; 
  VALUE rb_info;
  integer info; 
  VALUE rb_ab_out__;
  complex *ab_out__;
  complex *work;
  real *rwork;

  integer ldab;
  integer n;
  integer ldbb;
  integer ldx;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  x, info, ab = NumRu::Lapack.chbgst( vect, uplo, ka, kb, ab, bb)\n    or\n  NumRu::Lapack.chbgst  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE CHBGST( VECT, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, X, LDX, WORK, RWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  CHBGST reduces a complex Hermitian-definite banded generalized\n*  eigenproblem  A*x = lambda*B*x  to standard form  C*y = lambda*y,\n*  such that C has the same bandwidth as A.\n*\n*  B must have been previously factorized as S**H*S by CPBSTF, using a\n*  split Cholesky factorization. A is overwritten by C = X**H*A*X, where\n*  X = S**(-1)*Q and Q is a unitary matrix chosen to preserve the\n*  bandwidth of A.\n*\n\n*  Arguments\n*  =========\n*\n*  VECT    (input) CHARACTER*1\n*          = 'N':  do not form the transformation matrix X;\n*          = 'V':  form X.\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U':  Upper triangle of A is stored;\n*          = 'L':  Lower triangle of A is stored.\n*\n*  N       (input) INTEGER\n*          The order of the matrices A and B.  N >= 0.\n*\n*  KA      (input) INTEGER\n*          The number of superdiagonals of the matrix A if UPLO = 'U',\n*          or the number of subdiagonals if UPLO = 'L'.  KA >= 0.\n*\n*  KB      (input) INTEGER\n*          The number of superdiagonals of the matrix B if UPLO = 'U',\n*          or the number of subdiagonals if UPLO = 'L'.  KA >= KB >= 0.\n*\n*  AB      (input/output) COMPLEX array, dimension (LDAB,N)\n*          On entry, the upper or lower triangle of the Hermitian band\n*          matrix A, stored in the first ka+1 rows of the array.  The\n*          j-th column of A is stored in the j-th column of the array AB\n*          as follows:\n*          if UPLO = 'U', AB(ka+1+i-j,j) = A(i,j) for max(1,j-ka)<=i<=j;\n*          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+ka).\n*\n*          On exit, the transformed matrix X**H*A*X, stored in the same\n*          format as A.\n*\n*  LDAB    (input) INTEGER\n*          The leading dimension of the array AB.  LDAB >= KA+1.\n*\n*  BB      (input) COMPLEX array, dimension (LDBB,N)\n*          The banded factor S from the split Cholesky factorization of\n*          B, as returned by CPBSTF, stored in the first kb+1 rows of\n*          the array.\n*\n*  LDBB    (input) INTEGER\n*          The leading dimension of the array BB.  LDBB >= KB+1.\n*\n*  X       (output) COMPLEX array, dimension (LDX,N)\n*          If VECT = 'V', the n-by-n matrix X.\n*          If VECT = 'N', the array X is not referenced.\n*\n*  LDX     (input) INTEGER\n*          The leading dimension of the array X.\n*          LDX >= max(1,N) if VECT = 'V'; LDX >= 1 otherwise.\n*\n*  WORK    (workspace) COMPLEX array, dimension (N)\n*\n*  RWORK   (workspace) REAL array, dimension (N)\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rb_vect = argv[0];
  rb_uplo = argv[1];
  rb_ka = argv[2];
  rb_kb = argv[3];
  rb_ab = argv[4];
  rb_bb = argv[5];

  uplo = StringValueCStr(rb_uplo)[0];
  if (!NA_IsNArray(rb_bb))
    rb_raise(rb_eArgError, "bb (6th argument) must be NArray");
  if (NA_RANK(rb_bb) != 2)
    rb_raise(rb_eArgError, "rank of bb (6th argument) must be %d", 2);
  n = NA_SHAPE1(rb_bb);
  ldbb = NA_SHAPE0(rb_bb);
  if (NA_TYPE(rb_bb) != NA_SCOMPLEX)
    rb_bb = na_change_type(rb_bb, NA_SCOMPLEX);
  bb = NA_PTR_TYPE(rb_bb, complex*);
  ka = NUM2INT(rb_ka);
  vect = StringValueCStr(rb_vect)[0];
  kb = NUM2INT(rb_kb);
  if (!NA_IsNArray(rb_ab))
    rb_raise(rb_eArgError, "ab (5th argument) must be NArray");
  if (NA_RANK(rb_ab) != 2)
    rb_raise(rb_eArgError, "rank of ab (5th argument) must be %d", 2);
  if (NA_SHAPE1(rb_ab) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of ab must be the same as shape 1 of bb");
  ldab = NA_SHAPE0(rb_ab);
  if (NA_TYPE(rb_ab) != NA_SCOMPLEX)
    rb_ab = na_change_type(rb_ab, NA_SCOMPLEX);
  ab = NA_PTR_TYPE(rb_ab, complex*);
  ldx = lsame_(&vect,"V") ? MAX(1,n) : 1;
  {
    int shape[2];
    shape[0] = ldx;
    shape[1] = n;
    rb_x = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  x = NA_PTR_TYPE(rb_x, complex*);
  {
    int shape[2];
    shape[0] = ldab;
    shape[1] = n;
    rb_ab_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  ab_out__ = NA_PTR_TYPE(rb_ab_out__, complex*);
  MEMCPY(ab_out__, ab, complex, NA_TOTAL(rb_ab));
  rb_ab = rb_ab_out__;
  ab = ab_out__;
  work = ALLOC_N(complex, (n));
  rwork = ALLOC_N(real, (n));

  chbgst_(&vect, &uplo, &n, &ka, &kb, ab, &ldab, bb, &ldbb, x, &ldx, work, rwork, &info);

  free(work);
  free(rwork);
  rb_info = INT2NUM(info);
  return rb_ary_new3(3, rb_x, rb_info, rb_ab);
}

void
init_lapack_chbgst(VALUE mLapack){
  rb_define_module_function(mLapack, "chbgst", rb_chbgst, -1);
}
