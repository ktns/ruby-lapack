#include "rb_lapack.h"

static VALUE
rb_dsfrk(int argc, VALUE *argv, VALUE self){
  VALUE rb_transr;
  char transr; 
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_trans;
  char trans; 
  VALUE rb_n;
  integer n; 
  VALUE rb_k;
  integer k; 
  VALUE rb_alpha;
  doublereal alpha; 
  VALUE rb_a;
  doublereal a; 
  VALUE rb_lda;
  integer lda; 
  VALUE rb_beta;
  doublereal beta; 
  VALUE rb_c;
  doublereal c; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n   = NumRu::Lapack.dsfrk( transr, uplo, trans, n, k, alpha, a, lda, beta, c)\n    or\n  NumRu::Lapack.dsfrk  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DSFRK( TRANSR, UPLO, TRANS, N, K, ALPHA, A, LDA, BETA, C )\n\n*  Purpose\n*  =======\n*\n*  Level 3 BLAS like routine for C in RFP Format.\n*\n*  DSFRK performs one of the symmetric rank--k operations\n*\n*     C := alpha*A*A' + beta*C,\n*\n*  or\n*\n*     C := alpha*A'*A + beta*C,\n*\n*  where alpha and beta are real scalars, C is an n--by--n symmetric\n*  matrix and A is an n--by--k matrix in the first case and a k--by--n\n*  matrix in the second case.\n*\n\n*  Arguments\n*  ==========\n*\n*  TRANSR    (input) CHARACTER\n*          = 'N':  The Normal Form of RFP A is stored;\n*          = 'T':  The Transpose Form of RFP A is stored.\n*\n*  UPLO   - (input) CHARACTER\n*           On  entry, UPLO specifies whether the upper or lower\n*           triangular part of the array C is to be referenced as\n*           follows:\n*\n*              UPLO = 'U' or 'u'   Only the upper triangular part of C\n*                                  is to be referenced.\n*\n*              UPLO = 'L' or 'l'   Only the lower triangular part of C\n*                                  is to be referenced.\n*\n*           Unchanged on exit.\n*\n*  TRANS  - (input) CHARACTER\n*           On entry, TRANS specifies the operation to be performed as\n*           follows:\n*\n*              TRANS = 'N' or 'n'   C := alpha*A*A' + beta*C.\n*\n*              TRANS = 'T' or 't'   C := alpha*A'*A + beta*C.\n*\n*           Unchanged on exit.\n*\n*  N      - (input) INTEGER.\n*           On entry, N specifies the order of the matrix C. N must be\n*           at least zero.\n*           Unchanged on exit.\n*\n*  K      - (input) INTEGER.\n*           On entry with TRANS = 'N' or 'n', K specifies the number\n*           of  columns of the matrix A, and on entry with TRANS = 'T'\n*           or 't', K specifies the number of rows of the matrix A. K\n*           must be at least zero.\n*           Unchanged on exit.\n*\n*  ALPHA  - (input) DOUBLE PRECISION.\n*           On entry, ALPHA specifies the scalar alpha.\n*           Unchanged on exit.\n*\n*  A      - (input) DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where KA\n*           is K  when TRANS = 'N' or 'n', and is N otherwise. Before\n*           entry with TRANS = 'N' or 'n', the leading N--by--K part of\n*           the array A must contain the matrix A, otherwise the leading\n*           K--by--N part of the array A must contain the matrix A.\n*           Unchanged on exit.\n*\n*  LDA    - (input) INTEGER.\n*           On entry, LDA specifies the first dimension of A as declared\n*           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'\n*           then  LDA must be at least  max( 1, n ), otherwise  LDA must\n*           be at least  max( 1, k ).\n*           Unchanged on exit.\n*\n*  BETA   - (input) DOUBLE PRECISION.\n*           On entry, BETA specifies the scalar beta.\n*           Unchanged on exit.\n*\n*\n*  C      - (input/output) DOUBLE PRECISION array, dimension ( NT );\n*           NT = N*(N+1)/2. On entry, the symmetric matrix C in RFP\n*           Format. RFP Format is described by TRANSR, UPLO and N.\n*\n*  Arguments\n*  ==========\n*\n*     ..\n\n");
    return Qnil;
  }
  if (argc != 10)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 10)", argc);
  rb_transr = argv[0];
  rb_uplo = argv[1];
  rb_trans = argv[2];
  rb_n = argv[3];
  rb_k = argv[4];
  rb_alpha = argv[5];
  rb_a = argv[6];
  rb_lda = argv[7];
  rb_beta = argv[8];
  rb_c = argv[9];

  transr = StringValueCStr(rb_transr)[0];
  uplo = StringValueCStr(rb_uplo)[0];
  trans = StringValueCStr(rb_trans)[0];
  n = NUM2INT(rb_n);
  k = NUM2INT(rb_k);
  alpha = NUM2DBL(rb_alpha);
  a = NUM2DBL(rb_a);
  lda = NUM2INT(rb_lda);
  beta = NUM2DBL(rb_beta);
  c = NUM2DBL(rb_c);

  dsfrk_(&transr, &uplo, &trans, &n, &k, &alpha, &a, &lda, &beta, &c);

  return Qnil;
}

void
init_lapack_dsfrk(VALUE mLapack){
  rb_define_module_function(mLapack, "dsfrk", rb_dsfrk, -1);
}
