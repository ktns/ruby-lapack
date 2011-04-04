#include "rb_lapack.h"

extern VOID dpbstf_(char *uplo, integer *n, integer *kd, doublereal *ab, integer *ldab, integer *info);

static VALUE
rb_dpbstf(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_kd;
  integer kd; 
  VALUE rb_ab;
  doublereal *ab; 
  VALUE rb_info;
  integer info; 
  VALUE rb_ab_out__;
  doublereal *ab_out__;

  integer ldab;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, ab = NumRu::Lapack.dpbstf( uplo, kd, ab)\n    or\n  NumRu::Lapack.dpbstf  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DPBSTF( UPLO, N, KD, AB, LDAB, INFO )\n\n*  Purpose\n*  =======\n*\n*  DPBSTF computes a split Cholesky factorization of a real\n*  symmetric positive definite band matrix A.\n*\n*  This routine is designed to be used in conjunction with DSBGST.\n*\n*  The factorization has the form  A = S**T*S  where S is a band matrix\n*  of the same bandwidth as A and the following structure:\n*\n*    S = ( U    )\n*        ( M  L )\n*\n*  where U is upper triangular of order m = (n+kd)/2, and L is lower\n*  triangular of order n-m.\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U':  Upper triangle of A is stored;\n*          = 'L':  Lower triangle of A is stored.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  KD      (input) INTEGER\n*          The number of superdiagonals of the matrix A if UPLO = 'U',\n*          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.\n*\n*  AB      (input/output) DOUBLE PRECISION array, dimension (LDAB,N)\n*          On entry, the upper or lower triangle of the symmetric band\n*          matrix A, stored in the first kd+1 rows of the array.  The\n*          j-th column of A is stored in the j-th column of the array AB\n*          as follows:\n*          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;\n*          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).\n*\n*          On exit, if INFO = 0, the factor S from the split Cholesky\n*          factorization A = S**T*S. See Further Details.\n*\n*  LDAB    (input) INTEGER\n*          The leading dimension of the array AB.  LDAB >= KD+1.\n*\n*  INFO    (output) INTEGER\n*          = 0: successful exit\n*          < 0: if INFO = -i, the i-th argument had an illegal value\n*          > 0: if INFO = i, the factorization could not be completed,\n*               because the updated element a(i,i) was negative; the\n*               matrix A is not positive definite.\n*\n\n*  Further Details\n*  ===============\n*\n*  The band storage scheme is illustrated by the following example, when\n*  N = 7, KD = 2:\n*\n*  S = ( s11  s12  s13                     )\n*      (      s22  s23  s24                )\n*      (           s33  s34                )\n*      (                s44                )\n*      (           s53  s54  s55           )\n*      (                s64  s65  s66      )\n*      (                     s75  s76  s77 )\n*\n*  If UPLO = 'U', the array AB holds:\n*\n*  on entry:                          on exit:\n*\n*   *    *   a13  a24  a35  a46  a57   *    *   s13  s24  s53  s64  s75\n*   *   a12  a23  a34  a45  a56  a67   *   s12  s23  s34  s54  s65  s76\n*  a11  a22  a33  a44  a55  a66  a77  s11  s22  s33  s44  s55  s66  s77\n*\n*  If UPLO = 'L', the array AB holds:\n*\n*  on entry:                          on exit:\n*\n*  a11  a22  a33  a44  a55  a66  a77  s11  s22  s33  s44  s55  s66  s77\n*  a21  a32  a43  a54  a65  a76   *   s12  s23  s34  s54  s65  s76   *\n*  a31  a42  a53  a64  a64   *    *   s13  s24  s53  s64  s75   *    *\n*\n*  Array elements marked * are not used by the routine.\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_uplo = argv[0];
  rb_kd = argv[1];
  rb_ab = argv[2];

  if (!NA_IsNArray(rb_ab))
    rb_raise(rb_eArgError, "ab (3th argument) must be NArray");
  if (NA_RANK(rb_ab) != 2)
    rb_raise(rb_eArgError, "rank of ab (3th argument) must be %d", 2);
  n = NA_SHAPE1(rb_ab);
  ldab = NA_SHAPE0(rb_ab);
  if (NA_TYPE(rb_ab) != NA_DFLOAT)
    rb_ab = na_change_type(rb_ab, NA_DFLOAT);
  ab = NA_PTR_TYPE(rb_ab, doublereal*);
  kd = NUM2INT(rb_kd);
  uplo = StringValueCStr(rb_uplo)[0];
  {
    int shape[2];
    shape[0] = ldab;
    shape[1] = n;
    rb_ab_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  ab_out__ = NA_PTR_TYPE(rb_ab_out__, doublereal*);
  MEMCPY(ab_out__, ab, doublereal, NA_TOTAL(rb_ab));
  rb_ab = rb_ab_out__;
  ab = ab_out__;

  dpbstf_(&uplo, &n, &kd, ab, &ldab, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_info, rb_ab);
}

void
init_lapack_dpbstf(VALUE mLapack){
  rb_define_module_function(mLapack, "dpbstf", rb_dpbstf, -1);
}
