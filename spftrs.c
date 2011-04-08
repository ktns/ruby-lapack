#include "rb_lapack.h"

extern VOID spftrs_(char *transr, char *uplo, integer *n, integer *nrhs, real *a, real *b, integer *ldb, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_spftrs(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_transr;
  char transr; 
  VALUE rblapack_uplo;
  char uplo; 
  VALUE rblapack_n;
  integer n; 
  VALUE rblapack_a;
  real *a; 
  VALUE rblapack_b;
  real *b; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_b_out__;
  real *b_out__;

  integer ldb;
  integer nrhs;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, b = NumRu::Lapack.spftrs( transr, uplo, n, a, b, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE SPFTRS( TRANSR, UPLO, N, NRHS, A, B, LDB, INFO )\n\n*  Purpose\n*  =======\n*\n*  SPFTRS solves a system of linear equations A*X = B with a symmetric\n*  positive definite matrix A using the Cholesky factorization\n*  A = U**T*U or A = L*L**T computed by SPFTRF.\n*\n\n*  Arguments\n*  =========\n*\n*  TRANSR    (input) CHARACTER*1\n*          = 'N':  The Normal TRANSR of RFP A is stored;\n*          = 'T':  The Transpose TRANSR of RFP A is stored.\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U':  Upper triangle of RFP A is stored;\n*          = 'L':  Lower triangle of RFP A is stored.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  NRHS    (input) INTEGER\n*          The number of right hand sides, i.e., the number of columns\n*          of the matrix B.  NRHS >= 0.\n*\n*  A       (input) REAL array, dimension ( N*(N+1)/2 )\n*          The triangular factor U or L from the Cholesky factorization\n*          of RFP A = U**H*U or RFP A = L*L**T, as computed by SPFTRF.\n*          See note below for more details about RFP A.\n*\n*  B       (input/output) REAL array, dimension (LDB,NRHS)\n*          On entry, the right hand side matrix B.\n*          On exit, the solution matrix X.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B.  LDB >= max(1,N).\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*\n\n*  Further Details\n*  ===============\n*\n*  We first consider Rectangular Full Packed (RFP) Format when N is\n*  even. We give an example where N = 6.\n*\n*      AP is Upper             AP is Lower\n*\n*   00 01 02 03 04 05       00\n*      11 12 13 14 15       10 11\n*         22 23 24 25       20 21 22\n*            33 34 35       30 31 32 33\n*               44 45       40 41 42 43 44\n*                  55       50 51 52 53 54 55\n*\n*\n*  Let TRANSR = 'N'. RFP holds AP as follows:\n*  For UPLO = 'U' the upper trapezoid A(0:5,0:2) consists of the last\n*  three columns of AP upper. The lower triangle A(4:6,0:2) consists of\n*  the transpose of the first three columns of AP upper.\n*  For UPLO = 'L' the lower trapezoid A(1:6,0:2) consists of the first\n*  three columns of AP lower. The upper triangle A(0:2,0:2) consists of\n*  the transpose of the last three columns of AP lower.\n*  This covers the case N even and TRANSR = 'N'.\n*\n*         RFP A                   RFP A\n*\n*        03 04 05                33 43 53\n*        13 14 15                00 44 54\n*        23 24 25                10 11 55\n*        33 34 35                20 21 22\n*        00 44 45                30 31 32\n*        01 11 55                40 41 42\n*        02 12 22                50 51 52\n*\n*  Now let TRANSR = 'T'. RFP A in both UPLO cases is just the\n*  transpose of RFP A above. One therefore gets:\n*\n*\n*           RFP A                   RFP A\n*\n*     03 13 23 33 00 01 02    33 00 10 20 30 40 50\n*     04 14 24 34 44 11 12    43 44 11 21 31 41 51\n*     05 15 25 35 45 55 22    53 54 55 22 32 42 52\n*\n*\n*  We then consider Rectangular Full Packed (RFP) Format when N is\n*  odd. We give an example where N = 5.\n*\n*     AP is Upper                 AP is Lower\n*\n*   00 01 02 03 04              00\n*      11 12 13 14              10 11\n*         22 23 24              20 21 22\n*            33 34              30 31 32 33\n*               44              40 41 42 43 44\n*\n*\n*  Let TRANSR = 'N'. RFP holds AP as follows:\n*  For UPLO = 'U' the upper trapezoid A(0:4,0:2) consists of the last\n*  three columns of AP upper. The lower triangle A(3:4,0:1) consists of\n*  the transpose of the first two columns of AP upper.\n*  For UPLO = 'L' the lower trapezoid A(0:4,0:2) consists of the first\n*  three columns of AP lower. The upper triangle A(0:1,1:2) consists of\n*  the transpose of the last two columns of AP lower.\n*  This covers the case N odd and TRANSR = 'N'.\n*\n*         RFP A                   RFP A\n*\n*        02 03 04                00 33 43\n*        12 13 14                10 11 44\n*        22 23 24                20 21 22\n*        00 33 34                30 31 32\n*        01 11 44                40 41 42\n*\n*  Now let TRANSR = 'T'. RFP A in both UPLO cases is just the\n*  transpose of RFP A above. One therefore gets:\n*\n*           RFP A                   RFP A\n*\n*     02 12 22 00 01             00 10 20 30 40 50\n*     03 13 23 33 11             33 11 21 31 41 51\n*     04 14 24 34 44             43 44 22 32 42 52\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, b = NumRu::Lapack.spftrs( transr, uplo, n, a, b, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rblapack_transr = argv[0];
  rblapack_uplo = argv[1];
  rblapack_n = argv[2];
  rblapack_a = argv[3];
  rblapack_b = argv[4];
  if (rb_options != Qnil) {
  }

  transr = StringValueCStr(rblapack_transr)[0];
  if (!NA_IsNArray(rblapack_b))
    rb_raise(rb_eArgError, "b (5th argument) must be NArray");
  if (NA_RANK(rblapack_b) != 2)
    rb_raise(rb_eArgError, "rank of b (5th argument) must be %d", 2);
  nrhs = NA_SHAPE1(rblapack_b);
  ldb = NA_SHAPE0(rblapack_b);
  if (NA_TYPE(rblapack_b) != NA_SFLOAT)
    rblapack_b = na_change_type(rblapack_b, NA_SFLOAT);
  b = NA_PTR_TYPE(rblapack_b, real*);
  uplo = StringValueCStr(rblapack_uplo)[0];
  n = NUM2INT(rblapack_n);
  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (4th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 1)
    rb_raise(rb_eArgError, "rank of a (4th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_a) != (n*(n+1)/2))
    rb_raise(rb_eRuntimeError, "shape 0 of a must be %d", n*(n+1)/2);
  if (NA_TYPE(rblapack_a) != NA_SFLOAT)
    rblapack_a = na_change_type(rblapack_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rblapack_a, real*);
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = nrhs;
    rblapack_b_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rblapack_b_out__, real*);
  MEMCPY(b_out__, b, real, NA_TOTAL(rblapack_b));
  rblapack_b = rblapack_b_out__;
  b = b_out__;

  spftrs_(&transr, &uplo, &n, &nrhs, a, b, &ldb, &info);

  rblapack_info = INT2NUM(info);
  return rb_ary_new3(2, rblapack_info, rblapack_b);
}

void
init_lapack_spftrs(VALUE mLapack){
  rb_define_module_function(mLapack, "spftrs", rblapack_spftrs, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
