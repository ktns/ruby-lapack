#include "rb_lapack.h"

extern VOID strttf_(char *transr, char *uplo, integer *n, real *a, integer *lda, real *arf, integer *info);

static VALUE
rb_strttf(int argc, VALUE *argv, VALUE self){
  VALUE rb_transr;
  char transr; 
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_a;
  real *a; 
  VALUE rb_arf;
  real *arf; 
  VALUE rb_info;
  integer info; 

  integer lda;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  arf, info = NumRu::Lapack.strttf( transr, uplo, a)\n    or\n  NumRu::Lapack.strttf  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE STRTTF( TRANSR, UPLO, N, A, LDA, ARF, INFO )\n\n*  Purpose\n*  =======\n*\n*  STRTTF copies a triangular matrix A from standard full format (TR)\n*  to rectangular full packed format (TF) .\n*\n\n*  Arguments\n*  =========\n*\n*  TRANSR   (input) CHARACTER*1\n*          = 'N':  ARF in Normal form is wanted;\n*          = 'T':  ARF in Transpose form is wanted.\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U':  Upper triangle of A is stored;\n*          = 'L':  Lower triangle of A is stored.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A. N >= 0.\n*\n*  A       (input) REAL array, dimension (LDA,N).\n*          On entry, the triangular matrix A.  If UPLO = 'U', the\n*          leading N-by-N upper triangular part of the array A contains\n*          the upper triangular matrix, and the strictly lower\n*          triangular part of A is not referenced.  If UPLO = 'L', the\n*          leading N-by-N lower triangular part of the array A contains\n*          the lower triangular matrix, and the strictly upper\n*          triangular part of A is not referenced.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the matrix A. LDA >= max(1,N).\n*\n*  ARF     (output) REAL array, dimension (NT).\n*          NT=N*(N+1)/2. On exit, the triangular matrix A in RFP format.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*\n\n*  Further Details\n*  ===============\n*\n*  We first consider Rectangular Full Packed (RFP) Format when N is\n*  even. We give an example where N = 6.\n*\n*      AP is Upper             AP is Lower\n*\n*   00 01 02 03 04 05       00\n*      11 12 13 14 15       10 11\n*         22 23 24 25       20 21 22\n*            33 34 35       30 31 32 33\n*               44 45       40 41 42 43 44\n*                  55       50 51 52 53 54 55\n*\n*\n*  Let TRANSR = 'N'. RFP holds AP as follows:\n*  For UPLO = 'U' the upper trapezoid A(0:5,0:2) consists of the last\n*  three columns of AP upper. The lower triangle A(4:6,0:2) consists of\n*  the transpose of the first three columns of AP upper.\n*  For UPLO = 'L' the lower trapezoid A(1:6,0:2) consists of the first\n*  three columns of AP lower. The upper triangle A(0:2,0:2) consists of\n*  the transpose of the last three columns of AP lower.\n*  This covers the case N even and TRANSR = 'N'.\n*\n*         RFP A                   RFP A\n*\n*        03 04 05                33 43 53\n*        13 14 15                00 44 54\n*        23 24 25                10 11 55\n*        33 34 35                20 21 22\n*        00 44 45                30 31 32\n*        01 11 55                40 41 42\n*        02 12 22                50 51 52\n*\n*  Now let TRANSR = 'T'. RFP A in both UPLO cases is just the\n*  transpose of RFP A above. One therefore gets:\n*\n*\n*           RFP A                   RFP A\n*\n*     03 13 23 33 00 01 02    33 00 10 20 30 40 50\n*     04 14 24 34 44 11 12    43 44 11 21 31 41 51\n*     05 15 25 35 45 55 22    53 54 55 22 32 42 52\n*\n*\n*  We then consider Rectangular Full Packed (RFP) Format when N is\n*  odd. We give an example where N = 5.\n*\n*     AP is Upper                 AP is Lower\n*\n*   00 01 02 03 04              00\n*      11 12 13 14              10 11\n*         22 23 24              20 21 22\n*            33 34              30 31 32 33\n*               44              40 41 42 43 44\n*\n*\n*  Let TRANSR = 'N'. RFP holds AP as follows:\n*  For UPLO = 'U' the upper trapezoid A(0:4,0:2) consists of the last\n*  three columns of AP upper. The lower triangle A(3:4,0:1) consists of\n*  the transpose of the first two columns of AP upper.\n*  For UPLO = 'L' the lower trapezoid A(0:4,0:2) consists of the first\n*  three columns of AP lower. The upper triangle A(0:1,1:2) consists of\n*  the transpose of the last two columns of AP lower.\n*  This covers the case N odd and TRANSR = 'N'.\n*\n*         RFP A                   RFP A\n*\n*        02 03 04                00 33 43\n*        12 13 14                10 11 44\n*        22 23 24                20 21 22\n*        00 33 34                30 31 32\n*        01 11 44                40 41 42\n*\n*  Now let TRANSR = 'T'. RFP A in both UPLO cases is just the\n*  transpose of RFP A above. One therefore gets:\n*\n*           RFP A                   RFP A\n*\n*     02 12 22 00 01             00 10 20 30 40 50\n*     03 13 23 33 11             33 11 21 31 41 51\n*     04 14 24 34 44             43 44 22 32 42 52\n*\n*  Reference\n*  =========\n*\n*  =====================================================================\n*\n*     ..\n*     .. Local Scalars ..\n      LOGICAL            LOWER, NISODD, NORMALTRANSR\n      INTEGER            I, IJ, J, K, L, N1, N2, NT, NX2, NP1X2\n*     ..\n*     .. External Functions ..\n      LOGICAL            LSAME\n      EXTERNAL           LSAME\n*     ..\n*     .. External Subroutines ..\n      EXTERNAL           XERBLA\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          MAX, MOD\n*     ..\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_transr = argv[0];
  rb_uplo = argv[1];
  rb_a = argv[2];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);
  uplo = StringValueCStr(rb_uplo)[0];
  transr = StringValueCStr(rb_transr)[0];
  {
    int shape[1];
    shape[0] = n*(n+1)/2;
    rb_arf = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  arf = NA_PTR_TYPE(rb_arf, real*);

  strttf_(&transr, &uplo, &n, a, &lda, arf, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_arf, rb_info);
}

void
init_lapack_strttf(VALUE mLapack){
  rb_define_module_function(mLapack, "strttf", rb_strttf, -1);
}
