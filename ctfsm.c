#include "rb_lapack.h"

extern VOID ctfsm_(char *transr, char *side, char *uplo, char *trans, char *diag, integer *m, integer *n, complex *alpha, complex *a, complex *b, integer *ldb);

static VALUE
rb_ctfsm(int argc, VALUE *argv, VALUE self){
  VALUE rb_transr;
  char transr; 
  VALUE rb_side;
  char side; 
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_trans;
  char trans; 
  VALUE rb_diag;
  char diag; 
  VALUE rb_m;
  integer m; 
  VALUE rb_alpha;
  complex alpha; 
  VALUE rb_a;
  complex *a; 
  VALUE rb_b;
  complex *b; 
  VALUE rb_b_out__;
  complex *b_out__;

  integer ldb;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  b = NumRu::Lapack.ctfsm( transr, side, uplo, trans, diag, m, alpha, a, b)\n    or\n  NumRu::Lapack.ctfsm  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE CTFSM( TRANSR, SIDE, UPLO, TRANS, DIAG, M, N, ALPHA, A, B, LDB )\n\n*  Purpose\n*  =======\n*\n*  Level 3 BLAS like routine for A in RFP Format.\n*\n*  CTFSM solves the matrix equation\n*\n*     op( A )*X = alpha*B  or  X*op( A ) = alpha*B\n*\n*  where alpha is a scalar, X and B are m by n matrices, A is a unit, or\n*  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of\n*\n*     op( A ) = A   or   op( A ) = conjg( A' ).\n*\n*  A is in Rectangular Full Packed (RFP) Format.\n*\n*  The matrix X is overwritten on B.\n*\n\n*  Arguments\n*  ==========\n*\n*  TRANSR  (input) CHARACTER*1\n*          = 'N':  The Normal Form of RFP A is stored;\n*          = 'C':  The Conjugate-transpose Form of RFP A is stored.\n*\n*  SIDE    (input) CHARACTER*1\n*           On entry, SIDE specifies whether op( A ) appears on the left\n*           or right of X as follows:\n*\n*              SIDE = 'L' or 'l'   op( A )*X = alpha*B.\n*\n*              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.\n*\n*           Unchanged on exit.\n*\n*  UPLO    (input) CHARACTER*1\n*           On entry, UPLO specifies whether the RFP matrix A came from\n*           an upper or lower triangular matrix as follows:\n*           UPLO = 'U' or 'u' RFP A came from an upper triangular matrix\n*           UPLO = 'L' or 'l' RFP A came from a  lower triangular matrix\n*\n*           Unchanged on exit.\n*\n*  TRANS   (input) CHARACTER*1\n*           On entry, TRANS  specifies the form of op( A ) to be used\n*           in the matrix multiplication as follows:\n*\n*              TRANS  = 'N' or 'n'   op( A ) = A.\n*\n*              TRANS  = 'C' or 'c'   op( A ) = conjg( A' ).\n*\n*           Unchanged on exit.\n*\n*  DIAG    (input) CHARACTER*1\n*           On entry, DIAG specifies whether or not RFP A is unit\n*           triangular as follows:\n*\n*              DIAG = 'U' or 'u'   A is assumed to be unit triangular.\n*\n*              DIAG = 'N' or 'n'   A is not assumed to be unit\n*                                  triangular.\n*\n*           Unchanged on exit.\n*\n*  M       (input) INTEGER\n*           On entry, M specifies the number of rows of B. M must be at\n*           least zero.\n*           Unchanged on exit.\n*\n*  N       (input) INTEGER\n*           On entry, N specifies the number of columns of B.  N must be\n*           at least zero.\n*           Unchanged on exit.\n*\n*  ALPHA   (input) COMPLEX\n*           On entry,  ALPHA specifies the scalar  alpha. When  alpha is\n*           zero then  A is not referenced and  B need not be set before\n*           entry.\n*           Unchanged on exit.\n*\n*  A       (input) COMPLEX array, dimension (N*(N+1)/2)\n*           NT = N*(N+1)/2. On entry, the matrix A in RFP Format.\n*           RFP Format is described by TRANSR, UPLO and N as follows:\n*           If TRANSR='N' then RFP A is (0:N,0:K-1) when N is even;\n*           K=N/2. RFP A is (0:N-1,0:K) when N is odd; K=N/2. If\n*           TRANSR = 'C' then RFP is the Conjugate-transpose of RFP A as\n*           defined when TRANSR = 'N'. The contents of RFP A are defined\n*           by UPLO as follows: If UPLO = 'U' the RFP A contains the NT\n*           elements of upper packed A either in normal or\n*           conjugate-transpose Format. If UPLO = 'L' the RFP A contains\n*           the NT elements of lower packed A either in normal or\n*           conjugate-transpose Format. The LDA of RFP A is (N+1)/2 when\n*           TRANSR = 'C'. When TRANSR is 'N' the LDA is N+1 when N is\n*           even and is N when is odd.\n*           See the Note below for more details. Unchanged on exit.\n*\n*  B       (input/output) COMPLEX array,  dimension (LDB,N)\n*           Before entry,  the leading  m by n part of the array  B must\n*           contain  the  right-hand  side  matrix  B,  and  on exit  is\n*           overwritten by the solution matrix  X.\n*\n*  LDB     (input) INTEGER\n*           On entry, LDB specifies the first dimension of B as declared\n*           in  the  calling  (sub)  program.   LDB  must  be  at  least\n*           max( 1, m ).\n*           Unchanged on exit.\n*\n\n*  Further Details\n*  ===============\n*\n*  We first consider Standard Packed Format when N is even.\n*  We give an example where N = 6.\n*\n*      AP is Upper             AP is Lower\n*\n*   00 01 02 03 04 05       00\n*      11 12 13 14 15       10 11\n*         22 23 24 25       20 21 22\n*            33 34 35       30 31 32 33\n*               44 45       40 41 42 43 44\n*                  55       50 51 52 53 54 55\n*\n*\n*  Let TRANSR = 'N'. RFP holds AP as follows:\n*  For UPLO = 'U' the upper trapezoid A(0:5,0:2) consists of the last\n*  three columns of AP upper. The lower triangle A(4:6,0:2) consists of\n*  conjugate-transpose of the first three columns of AP upper.\n*  For UPLO = 'L' the lower trapezoid A(1:6,0:2) consists of the first\n*  three columns of AP lower. The upper triangle A(0:2,0:2) consists of\n*  conjugate-transpose of the last three columns of AP lower.\n*  To denote conjugate we place -- above the element. This covers the\n*  case N even and TRANSR = 'N'.\n*\n*         RFP A                   RFP A\n*\n*                                -- -- --\n*        03 04 05                33 43 53\n*                                   -- --\n*        13 14 15                00 44 54\n*                                      --\n*        23 24 25                10 11 55\n*\n*        33 34 35                20 21 22\n*        --\n*        00 44 45                30 31 32\n*        -- --\n*        01 11 55                40 41 42\n*        -- -- --\n*        02 12 22                50 51 52\n*\n*  Now let TRANSR = 'C'. RFP A in both UPLO cases is just the conjugate-\n*  transpose of RFP A above. One therefore gets:\n*\n*\n*           RFP A                   RFP A\n*\n*     -- -- -- --                -- -- -- -- -- --\n*     03 13 23 33 00 01 02    33 00 10 20 30 40 50\n*     -- -- -- -- --                -- -- -- -- --\n*     04 14 24 34 44 11 12    43 44 11 21 31 41 51\n*     -- -- -- -- -- --                -- -- -- --\n*     05 15 25 35 45 55 22    53 54 55 22 32 42 52\n*\n*\n*  We next  consider Standard Packed Format when N is odd.\n*  We give an example where N = 5.\n*\n*     AP is Upper                 AP is Lower\n*\n*   00 01 02 03 04              00\n*      11 12 13 14              10 11\n*         22 23 24              20 21 22\n*            33 34              30 31 32 33\n*               44              40 41 42 43 44\n*\n*\n*  Let TRANSR = 'N'. RFP holds AP as follows:\n*  For UPLO = 'U' the upper trapezoid A(0:4,0:2) consists of the last\n*  three columns of AP upper. The lower triangle A(3:4,0:1) consists of\n*  conjugate-transpose of the first two   columns of AP upper.\n*  For UPLO = 'L' the lower trapezoid A(0:4,0:2) consists of the first\n*  three columns of AP lower. The upper triangle A(0:1,1:2) consists of\n*  conjugate-transpose of the last two   columns of AP lower.\n*  To denote conjugate we place -- above the element. This covers the\n*  case N odd  and TRANSR = 'N'.\n*\n*         RFP A                   RFP A\n*\n*                                   -- --\n*        02 03 04                00 33 43\n*                                      --\n*        12 13 14                10 11 44\n*\n*        22 23 24                20 21 22\n*        --\n*        00 33 34                30 31 32\n*        -- --\n*        01 11 44                40 41 42\n*\n*  Now let TRANSR = 'C'. RFP A in both UPLO cases is just the conjugate-\n*  transpose of RFP A above. One therefore gets:\n*\n*\n*           RFP A                   RFP A\n*\n*     -- -- --                   -- -- -- -- -- --\n*     02 12 22 00 01             00 10 20 30 40 50\n*     -- -- -- --                   -- -- -- -- --\n*     03 13 23 33 11             33 11 21 31 41 51\n*     -- -- -- -- --                   -- -- -- --\n*     04 14 24 34 44             43 44 22 32 42 52\n*\n*     ..\n\n");
    return Qnil;
  }
  if (argc != 9)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 9)", argc);
  rb_transr = argv[0];
  rb_side = argv[1];
  rb_uplo = argv[2];
  rb_trans = argv[3];
  rb_diag = argv[4];
  rb_m = argv[5];
  rb_alpha = argv[6];
  rb_a = argv[7];
  rb_b = argv[8];

  transr = StringValueCStr(rb_transr)[0];
  side = StringValueCStr(rb_side)[0];
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (9th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (9th argument) must be %d", 2);
  n = NA_SHAPE1(rb_b);
  ldb = NA_SHAPE0(rb_b);
  if (NA_TYPE(rb_b) != NA_SCOMPLEX)
    rb_b = na_change_type(rb_b, NA_SCOMPLEX);
  b = NA_PTR_TYPE(rb_b, complex*);
  m = NUM2INT(rb_m);
  uplo = StringValueCStr(rb_uplo)[0];
  alpha.r = (real)NUM2DBL(rb_funcall(rb_alpha, rb_intern("real"), 0));
  alpha.i = (real)NUM2DBL(rb_funcall(rb_alpha, rb_intern("imag"), 0));
  trans = StringValueCStr(rb_trans)[0];
  diag = StringValueCStr(rb_diag)[0];
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (8th argument) must be NArray");
  if (NA_RANK(rb_a) != 1)
    rb_raise(rb_eArgError, "rank of a (8th argument) must be %d", 1);
  if (NA_SHAPE0(rb_a) != (n*(n+1)/2))
    rb_raise(rb_eRuntimeError, "shape 0 of a must be %d", n*(n+1)/2);
  if (NA_TYPE(rb_a) != NA_SCOMPLEX)
    rb_a = na_change_type(rb_a, NA_SCOMPLEX);
  a = NA_PTR_TYPE(rb_a, complex*);
  {
    int shape[2];
    shape[0] = ldb;
    shape[1] = n;
    rb_b_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  b_out__ = NA_PTR_TYPE(rb_b_out__, complex*);
  MEMCPY(b_out__, b, complex, NA_TOTAL(rb_b));
  rb_b = rb_b_out__;
  b = b_out__;

  ctfsm_(&transr, &side, &uplo, &trans, &diag, &m, &n, &alpha, a, b, &ldb);

  return rb_b;
}

void
init_lapack_ctfsm(VALUE mLapack){
  rb_define_module_function(mLapack, "ctfsm", rb_ctfsm, -1);
}
