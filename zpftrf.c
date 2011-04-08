#include "rb_lapack.h"

extern VOID zpftrf_(char *transr, char *uplo, integer *n, complex *a, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_zpftrf(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_transr;
  char transr; 
  VALUE rblapack_uplo;
  char uplo; 
  VALUE rblapack_n;
  integer n; 
  VALUE rblapack_a;
  complex *a; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_a_out__;
  complex *a_out__;


  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, a = NumRu::Lapack.zpftrf( transr, uplo, n, a, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZPFTRF( TRANSR, UPLO, N, A, INFO )\n\n*  Purpose\n*  =======\n*\n*  ZPFTRF computes the Cholesky factorization of a complex Hermitian\n*  positive definite matrix A.\n*\n*  The factorization has the form\n*     A = U**H * U,  if UPLO = 'U', or\n*     A = L  * L**H,  if UPLO = 'L',\n*  where U is an upper triangular matrix and L is lower triangular.\n*\n*  This is the block version of the algorithm, calling Level 3 BLAS.\n*\n\n*  Arguments\n*  =========\n*\n*  TRANSR    (input) CHARACTER*1\n*          = 'N':  The Normal TRANSR of RFP A is stored;\n*          = 'C':  The Conjugate-transpose TRANSR of RFP A is stored.\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U':  Upper triangle of RFP A is stored;\n*          = 'L':  Lower triangle of RFP A is stored.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  A       (input/output) COMPLEX array, dimension ( N*(N+1)/2 );\n*          On entry, the Hermitian matrix A in RFP format. RFP format is\n*          described by TRANSR, UPLO, and N as follows: If TRANSR = 'N'\n*          then RFP A is (0:N,0:k-1) when N is even; k=N/2. RFP A is\n*          (0:N-1,0:k) when N is odd; k=N/2. IF TRANSR = 'C' then RFP is\n*          the Conjugate-transpose of RFP A as defined when\n*          TRANSR = 'N'. The contents of RFP A are defined by UPLO as\n*          follows: If UPLO = 'U' the RFP A contains the nt elements of\n*          upper packed A. If UPLO = 'L' the RFP A contains the elements\n*          of lower packed A. The LDA of RFP A is (N+1)/2 when TRANSR =\n*          'C'. When TRANSR is 'N' the LDA is N+1 when N is even and N\n*          is odd. See the Note below for more details.\n*\n*          On exit, if INFO = 0, the factor U or L from the Cholesky\n*          factorization RFP A = U**H*U or RFP A = L*L**H.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*          > 0:  if INFO = i, the leading minor of order i is not\n*                positive definite, and the factorization could not be\n*                completed.\n*\n*  Further Notes on RFP Format:\n*  ============================\n*\n*  We first consider Standard Packed Format when N is even.\n*  We give an example where N = 6.\n*\n*     AP is Upper             AP is Lower\n*\n*   00 01 02 03 04 05       00\n*      11 12 13 14 15       10 11\n*         22 23 24 25       20 21 22\n*            33 34 35       30 31 32 33\n*               44 45       40 41 42 43 44\n*                  55       50 51 52 53 54 55\n*\n*\n*  Let TRANSR = 'N'. RFP holds AP as follows:\n*  For UPLO = 'U' the upper trapezoid A(0:5,0:2) consists of the last\n*  three columns of AP upper. The lower triangle A(4:6,0:2) consists of\n*  conjugate-transpose of the first three columns of AP upper.\n*  For UPLO = 'L' the lower trapezoid A(1:6,0:2) consists of the first\n*  three columns of AP lower. The upper triangle A(0:2,0:2) consists of\n*  conjugate-transpose of the last three columns of AP lower.\n*  To denote conjugate we place -- above the element. This covers the\n*  case N even and TRANSR = 'N'.\n*\n*         RFP A                   RFP A\n*\n*                                -- -- --\n*        03 04 05                33 43 53\n*                                   -- --\n*        13 14 15                00 44 54\n*                                      --\n*        23 24 25                10 11 55\n*\n*        33 34 35                20 21 22\n*        --\n*        00 44 45                30 31 32\n*        -- --\n*        01 11 55                40 41 42\n*        -- -- --\n*        02 12 22                50 51 52\n*\n*  Now let TRANSR = 'C'. RFP A in both UPLO cases is just the conjugate-\n*  transpose of RFP A above. One therefore gets:\n*\n*\n*           RFP A                   RFP A\n*\n*     -- -- -- --                -- -- -- -- -- --\n*     03 13 23 33 00 01 02    33 00 10 20 30 40 50\n*     -- -- -- -- --                -- -- -- -- --\n*     04 14 24 34 44 11 12    43 44 11 21 31 41 51\n*     -- -- -- -- -- --                -- -- -- --\n*     05 15 25 35 45 55 22    53 54 55 22 32 42 52\n*\n*\n*  We next  consider Standard Packed Format when N is odd.\n*  We give an example where N = 5.\n*\n*     AP is Upper                 AP is Lower\n*\n*   00 01 02 03 04              00\n*      11 12 13 14              10 11\n*         22 23 24              20 21 22\n*            33 34              30 31 32 33\n*               44              40 41 42 43 44\n*\n*\n*  Let TRANSR = 'N'. RFP holds AP as follows:\n*  For UPLO = 'U' the upper trapezoid A(0:4,0:2) consists of the last\n*  three columns of AP upper. The lower triangle A(3:4,0:1) consists of\n*  conjugate-transpose of the first two   columns of AP upper.\n*  For UPLO = 'L' the lower trapezoid A(0:4,0:2) consists of the first\n*  three columns of AP lower. The upper triangle A(0:1,1:2) consists of\n*  conjugate-transpose of the last two   columns of AP lower.\n*  To denote conjugate we place -- above the element. This covers the\n*  case N odd  and TRANSR = 'N'.\n*\n*         RFP A                   RFP A\n*\n*                                   -- --\n*        02 03 04                00 33 43\n*                                      --\n*        12 13 14                10 11 44\n*\n*        22 23 24                20 21 22\n*        --\n*        00 33 34                30 31 32\n*        -- --\n*        01 11 44                40 41 42\n*\n*  Now let TRANSR = 'C'. RFP A in both UPLO cases is just the conjugate-\n*  transpose of RFP A above. One therefore gets:\n*\n*\n*           RFP A                   RFP A\n*\n*     -- -- --                   -- -- -- -- -- --\n*     02 12 22 00 01             00 10 20 30 40 50\n*     -- -- -- --                   -- -- -- -- --\n*     03 13 23 33 11             33 11 21 31 41 51\n*     -- -- -- -- --                   -- -- -- --\n*     04 14 24 34 44             43 44 22 32 42 52\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, a = NumRu::Lapack.zpftrf( transr, uplo, n, a, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rblapack_transr = argv[0];
  rblapack_uplo = argv[1];
  rblapack_n = argv[2];
  rblapack_a = argv[3];
  if (rb_options != Qnil) {
  }

  transr = StringValueCStr(rblapack_transr)[0];
  n = NUM2INT(rblapack_n);
  uplo = StringValueCStr(rblapack_uplo)[0];
  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (4th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 1)
    rb_raise(rb_eArgError, "rank of a (4th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_a) != (n*(n+1)/2))
    rb_raise(rb_eRuntimeError, "shape 0 of a must be %d", n*(n+1)/2);
  if (NA_TYPE(rblapack_a) != NA_SCOMPLEX)
    rblapack_a = na_change_type(rblapack_a, NA_SCOMPLEX);
  a = NA_PTR_TYPE(rblapack_a, complex*);
  {
    int shape[1];
    shape[0] = n*(n+1)/2;
    rblapack_a_out__ = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rblapack_a_out__, complex*);
  MEMCPY(a_out__, a, complex, NA_TOTAL(rblapack_a));
  rblapack_a = rblapack_a_out__;
  a = a_out__;

  zpftrf_(&transr, &uplo, &n, a, &info);

  rblapack_info = INT2NUM(info);
  return rb_ary_new3(2, rblapack_info, rblapack_a);
}

void
init_lapack_zpftrf(VALUE mLapack){
  rb_define_module_function(mLapack, "zpftrf", rblapack_zpftrf, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
