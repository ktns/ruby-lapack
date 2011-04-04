#include "rb_lapack.h"

extern VOID ctftri_(char *transr, char *uplo, char *diag, integer *n, complex *a, integer *info);

static VALUE
rb_ctftri(int argc, VALUE *argv, VALUE self){
  VALUE rb_transr;
  char transr; 
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_diag;
  char diag; 
  VALUE rb_n;
  integer n; 
  VALUE rb_a;
  complex *a; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  complex *a_out__;


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, a = NumRu::Lapack.ctftri( transr, uplo, diag, n, a)\n    or\n  NumRu::Lapack.ctftri  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE CTFTRI( TRANSR, UPLO, DIAG, N, A, INFO )\n\n*  Purpose\n*  =======\n*\n*  CTFTRI computes the inverse of a triangular matrix A stored in RFP\n*  format.\n*\n*  This is a Level 3 BLAS version of the algorithm.\n*\n\n*  Arguments\n*  =========\n*\n*  TRANSR    (input) CHARACTER*1\n*          = 'N':  The Normal TRANSR of RFP A is stored;\n*          = 'C':  The Conjugate-transpose TRANSR of RFP A is stored.\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U':  A is upper triangular;\n*          = 'L':  A is lower triangular.\n*\n*  DIAG    (input) CHARACTER*1\n*          = 'N':  A is non-unit triangular;\n*          = 'U':  A is unit triangular.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  A       (input/output) COMPLEX array, dimension ( N*(N+1)/2 );\n*          On entry, the triangular matrix A in RFP format. RFP format\n*          is described by TRANSR, UPLO, and N as follows: If TRANSR =\n*          'N' then RFP A is (0:N,0:k-1) when N is even; k=N/2. RFP A is\n*          (0:N-1,0:k) when N is odd; k=N/2. IF TRANSR = 'C' then RFP is\n*          the Conjugate-transpose of RFP A as defined when\n*          TRANSR = 'N'. The contents of RFP A are defined by UPLO as\n*          follows: If UPLO = 'U' the RFP A contains the nt elements of\n*          upper packed A; If UPLO = 'L' the RFP A contains the nt\n*          elements of lower packed A. The LDA of RFP A is (N+1)/2 when\n*          TRANSR = 'C'. When TRANSR is 'N' the LDA is N+1 when N is\n*          even and N is odd. See the Note below for more details.\n*\n*          On exit, the (triangular) inverse of the original matrix, in\n*          the same storage format.\n*\n*  INFO    (output) INTEGER\n*          = 0: successful exit\n*          < 0: if INFO = -i, the i-th argument had an illegal value\n*          > 0: if INFO = i, A(i,i) is exactly zero.  The triangular\n*               matrix is singular and its inverse can not be computed.\n*\n\n*  Further Details\n*  ===============\n*\n*  We first consider Standard Packed Format when N is even.\n*  We give an example where N = 6.\n*\n*      AP is Upper             AP is Lower\n*\n*   00 01 02 03 04 05       00\n*      11 12 13 14 15       10 11\n*         22 23 24 25       20 21 22\n*            33 34 35       30 31 32 33\n*               44 45       40 41 42 43 44\n*                  55       50 51 52 53 54 55\n*\n*\n*  Let TRANSR = 'N'. RFP holds AP as follows:\n*  For UPLO = 'U' the upper trapezoid A(0:5,0:2) consists of the last\n*  three columns of AP upper. The lower triangle A(4:6,0:2) consists of\n*  conjugate-transpose of the first three columns of AP upper.\n*  For UPLO = 'L' the lower trapezoid A(1:6,0:2) consists of the first\n*  three columns of AP lower. The upper triangle A(0:2,0:2) consists of\n*  conjugate-transpose of the last three columns of AP lower.\n*  To denote conjugate we place -- above the element. This covers the\n*  case N even and TRANSR = 'N'.\n*\n*         RFP A                   RFP A\n*\n*                                -- -- --\n*        03 04 05                33 43 53\n*                                   -- --\n*        13 14 15                00 44 54\n*                                      --\n*        23 24 25                10 11 55\n*\n*        33 34 35                20 21 22\n*        --\n*        00 44 45                30 31 32\n*        -- --\n*        01 11 55                40 41 42\n*        -- -- --\n*        02 12 22                50 51 52\n*\n*  Now let TRANSR = 'C'. RFP A in both UPLO cases is just the conjugate-\n*  transpose of RFP A above. One therefore gets:\n*\n*\n*           RFP A                   RFP A\n*\n*     -- -- -- --                -- -- -- -- -- --\n*     03 13 23 33 00 01 02    33 00 10 20 30 40 50\n*     -- -- -- -- --                -- -- -- -- --\n*     04 14 24 34 44 11 12    43 44 11 21 31 41 51\n*     -- -- -- -- -- --                -- -- -- --\n*     05 15 25 35 45 55 22    53 54 55 22 32 42 52\n*\n*\n*  We next  consider Standard Packed Format when N is odd.\n*  We give an example where N = 5.\n*\n*     AP is Upper                 AP is Lower\n*\n*   00 01 02 03 04              00\n*      11 12 13 14              10 11\n*         22 23 24              20 21 22\n*            33 34              30 31 32 33\n*               44              40 41 42 43 44\n*\n*\n*  Let TRANSR = 'N'. RFP holds AP as follows:\n*  For UPLO = 'U' the upper trapezoid A(0:4,0:2) consists of the last\n*  three columns of AP upper. The lower triangle A(3:4,0:1) consists of\n*  conjugate-transpose of the first two   columns of AP upper.\n*  For UPLO = 'L' the lower trapezoid A(0:4,0:2) consists of the first\n*  three columns of AP lower. The upper triangle A(0:1,1:2) consists of\n*  conjugate-transpose of the last two   columns of AP lower.\n*  To denote conjugate we place -- above the element. This covers the\n*  case N odd  and TRANSR = 'N'.\n*\n*         RFP A                   RFP A\n*\n*                                   -- --\n*        02 03 04                00 33 43\n*                                      --\n*        12 13 14                10 11 44\n*\n*        22 23 24                20 21 22\n*        --\n*        00 33 34                30 31 32\n*        -- --\n*        01 11 44                40 41 42\n*\n*  Now let TRANSR = 'C'. RFP A in both UPLO cases is just the conjugate-\n*  transpose of RFP A above. One therefore gets:\n*\n*\n*           RFP A                   RFP A\n*\n*     -- -- --                   -- -- -- -- -- --\n*     02 12 22 00 01             00 10 20 30 40 50\n*     -- -- -- --                   -- -- -- -- --\n*     03 13 23 33 11             33 11 21 31 41 51\n*     -- -- -- -- --                   -- -- -- --\n*     04 14 24 34 44             43 44 22 32 42 52\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_transr = argv[0];
  rb_uplo = argv[1];
  rb_diag = argv[2];
  rb_n = argv[3];
  rb_a = argv[4];

  transr = StringValueCStr(rb_transr)[0];
  diag = StringValueCStr(rb_diag)[0];
  uplo = StringValueCStr(rb_uplo)[0];
  n = NUM2INT(rb_n);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (5th argument) must be NArray");
  if (NA_RANK(rb_a) != 1)
    rb_raise(rb_eArgError, "rank of a (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_a) != (n*(n+1)/2))
    rb_raise(rb_eRuntimeError, "shape 0 of a must be %d", n*(n+1)/2);
  if (NA_TYPE(rb_a) != NA_SCOMPLEX)
    rb_a = na_change_type(rb_a, NA_SCOMPLEX);
  a = NA_PTR_TYPE(rb_a, complex*);
  {
    int shape[1];
    shape[0] = n*(n+1)/2;
    rb_a_out__ = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rb_a_out__, complex*);
  MEMCPY(a_out__, a, complex, NA_TOTAL(rb_a));
  rb_a = rb_a_out__;
  a = a_out__;

  ctftri_(&transr, &uplo, &diag, &n, a, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_info, rb_a);
}

void
init_lapack_ctftri(VALUE mLapack){
  rb_define_module_function(mLapack, "ctftri", rb_ctftri, -1);
}
