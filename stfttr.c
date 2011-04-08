#include "rb_lapack.h"

extern VOID stfttr_(char *transr, char *uplo, integer *n, real *arf, real *a, integer *lda, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_stfttr(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_transr;
  char transr; 
  VALUE rblapack_uplo;
  char uplo; 
  VALUE rblapack_arf;
  real *arf; 
  VALUE rblapack_a;
  real *a; 
  VALUE rblapack_info;
  integer info; 

  integer ldarf;
  integer lda;
  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  a, info = NumRu::Lapack.stfttr( transr, uplo, arf, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE STFTTR( TRANSR, UPLO, N, ARF, A, LDA, INFO )\n\n*  Purpose\n*  =======\n*\n*  STFTTR copies a triangular matrix A from rectangular full packed\n*  format (TF) to standard full format (TR).\n*\n\n*  Arguments\n*  =========\n*\n*  TRANSR   (input) CHARACTER*1\n*          = 'N':  ARF is in Normal format;\n*          = 'T':  ARF is in Transpose format.\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U':  A is upper triangular;\n*          = 'L':  A is lower triangular.\n*\n*  N       (input) INTEGER\n*          The order of the matrices ARF and A. N >= 0.\n*\n*  ARF     (input) REAL array, dimension (N*(N+1)/2).\n*          On entry, the upper (if UPLO = 'U') or lower (if UPLO = 'L')\n*          matrix A in RFP format. See the \"Notes\" below for more\n*          details.\n*\n*  A       (output) REAL array, dimension (LDA,N)\n*          On exit, the triangular matrix A.  If UPLO = 'U', the\n*          leading N-by-N upper triangular part of the array A contains\n*          the upper triangular matrix, and the strictly lower\n*          triangular part of A is not referenced.  If UPLO = 'L', the\n*          leading N-by-N lower triangular part of the array A contains\n*          the lower triangular matrix, and the strictly upper\n*          triangular part of A is not referenced.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,N).\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*\n\n*  Further Details\n*  ===============\n*\n*  We first consider Rectangular Full Packed (RFP) Format when N is\n*  even. We give an example where N = 6.\n*\n*      AP is Upper             AP is Lower\n*\n*   00 01 02 03 04 05       00\n*      11 12 13 14 15       10 11\n*         22 23 24 25       20 21 22\n*            33 34 35       30 31 32 33\n*               44 45       40 41 42 43 44\n*                  55       50 51 52 53 54 55\n*\n*\n*  Let TRANSR = 'N'. RFP holds AP as follows:\n*  For UPLO = 'U' the upper trapezoid A(0:5,0:2) consists of the last\n*  three columns of AP upper. The lower triangle A(4:6,0:2) consists of\n*  the transpose of the first three columns of AP upper.\n*  For UPLO = 'L' the lower trapezoid A(1:6,0:2) consists of the first\n*  three columns of AP lower. The upper triangle A(0:2,0:2) consists of\n*  the transpose of the last three columns of AP lower.\n*  This covers the case N even and TRANSR = 'N'.\n*\n*         RFP A                   RFP A\n*\n*        03 04 05                33 43 53\n*        13 14 15                00 44 54\n*        23 24 25                10 11 55\n*        33 34 35                20 21 22\n*        00 44 45                30 31 32\n*        01 11 55                40 41 42\n*        02 12 22                50 51 52\n*\n*  Now let TRANSR = 'T'. RFP A in both UPLO cases is just the\n*  transpose of RFP A above. One therefore gets:\n*\n*\n*           RFP A                   RFP A\n*\n*     03 13 23 33 00 01 02    33 00 10 20 30 40 50\n*     04 14 24 34 44 11 12    43 44 11 21 31 41 51\n*     05 15 25 35 45 55 22    53 54 55 22 32 42 52\n*\n*\n*  We then consider Rectangular Full Packed (RFP) Format when N is\n*  odd. We give an example where N = 5.\n*\n*     AP is Upper                 AP is Lower\n*\n*   00 01 02 03 04              00\n*      11 12 13 14              10 11\n*         22 23 24              20 21 22\n*            33 34              30 31 32 33\n*               44              40 41 42 43 44\n*\n*\n*  Let TRANSR = 'N'. RFP holds AP as follows:\n*  For UPLO = 'U' the upper trapezoid A(0:4,0:2) consists of the last\n*  three columns of AP upper. The lower triangle A(3:4,0:1) consists of\n*  the transpose of the first two columns of AP upper.\n*  For UPLO = 'L' the lower trapezoid A(0:4,0:2) consists of the first\n*  three columns of AP lower. The upper triangle A(0:1,1:2) consists of\n*  the transpose of the last two columns of AP lower.\n*  This covers the case N odd and TRANSR = 'N'.\n*\n*         RFP A                   RFP A\n*\n*        02 03 04                00 33 43\n*        12 13 14                10 11 44\n*        22 23 24                20 21 22\n*        00 33 34                30 31 32\n*        01 11 44                40 41 42\n*\n*  Now let TRANSR = 'T'. RFP A in both UPLO cases is just the\n*  transpose of RFP A above. One therefore gets:\n*\n*           RFP A                   RFP A\n*\n*     02 12 22 00 01             00 10 20 30 40 50\n*     03 13 23 33 11             33 11 21 31 41 51\n*     04 14 24 34 44             43 44 22 32 42 52\n*\n*  Reference\n*  =========\n*\n*  =====================================================================\n*\n*     ..\n*     .. Local Scalars ..\n      LOGICAL            LOWER, NISODD, NORMALTRANSR\n      INTEGER            N1, N2, K, NT, NX2, NP1X2\n      INTEGER            I, J, L, IJ\n*     ..\n*     .. External Functions ..\n      LOGICAL            LSAME\n      EXTERNAL           LSAME\n*     ..\n*     .. External Subroutines ..\n      EXTERNAL           XERBLA\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          MAX, MOD\n*     ..\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  a, info = NumRu::Lapack.stfttr( transr, uplo, arf, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rblapack_transr = argv[0];
  rblapack_uplo = argv[1];
  rblapack_arf = argv[2];
  if (rb_options != Qnil) {
  }

  transr = StringValueCStr(rblapack_transr)[0];
  uplo = StringValueCStr(rblapack_uplo)[0];
  if (!NA_IsNArray(rblapack_arf))
    rb_raise(rb_eArgError, "arf (3th argument) must be NArray");
  if (NA_RANK(rblapack_arf) != 1)
    rb_raise(rb_eArgError, "rank of arf (3th argument) must be %d", 1);
  ldarf = NA_SHAPE0(rblapack_arf);
  if (NA_TYPE(rblapack_arf) != NA_SFLOAT)
    rblapack_arf = na_change_type(rblapack_arf, NA_SFLOAT);
  arf = NA_PTR_TYPE(rblapack_arf, real*);
  n = ((int)sqrtf(ldarf*8+1.0f)-1)/2;
  lda = MAX(1,n);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rblapack_a = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  a = NA_PTR_TYPE(rblapack_a, real*);

  stfttr_(&transr, &uplo, &n, arf, a, &lda, &info);

  rblapack_info = INT2NUM(info);
  return rb_ary_new3(2, rblapack_a, rblapack_info);
}

void
init_lapack_stfttr(VALUE mLapack){
  rb_define_module_function(mLapack, "stfttr", rblapack_stfttr, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
