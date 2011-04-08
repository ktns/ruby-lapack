#include "rb_lapack.h"

extern VOID ssfrk_(char *transr, char *uplo, char *trans, integer *n, integer *k, real *alpha, real *a, integer *lda, real *beta, real *c);

static VALUE sHelp, sUsage;

static VALUE
rblapack_ssfrk(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_transr;
  char transr; 
  VALUE rblapack_uplo;
  char uplo; 
  VALUE rblapack_trans;
  char trans; 
  VALUE rblapack_n;
  integer n; 
  VALUE rblapack_k;
  integer k; 
  VALUE rblapack_alpha;
  real alpha; 
  VALUE rblapack_a;
  real *a; 
  VALUE rblapack_beta;
  real beta; 
  VALUE rblapack_c;
  real *c; 
  VALUE rblapack_c_out__;
  real *c_out__;

  integer lda;
  integer nt;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  c = NumRu::Lapack.ssfrk( transr, uplo, trans, n, k, alpha, a, beta, c, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE SSFRK( TRANSR, UPLO, TRANS, N, K, ALPHA, A, LDA, BETA, C )\n\n*  Purpose\n*  =======\n*\n*  Level 3 BLAS like routine for C in RFP Format.\n*\n*  SSFRK performs one of the symmetric rank--k operations\n*\n*     C := alpha*A*A' + beta*C,\n*\n*  or\n*\n*     C := alpha*A'*A + beta*C,\n*\n*  where alpha and beta are real scalars, C is an n--by--n symmetric\n*  matrix and A is an n--by--k matrix in the first case and a k--by--n\n*  matrix in the second case.\n*\n\n*  Arguments\n*  ==========\n*\n*  TRANSR   (input) CHARACTER*1\n*          = 'N':  The Normal Form of RFP A is stored;\n*          = 'T':  The Transpose Form of RFP A is stored.\n*\n*  UPLO     (input) CHARACTER*1\n*           On  entry, UPLO specifies whether the upper or lower\n*           triangular part of the array C is to be referenced as\n*           follows:\n*\n*              UPLO = 'U' or 'u'   Only the upper triangular part of C\n*                                  is to be referenced.\n*\n*              UPLO = 'L' or 'l'   Only the lower triangular part of C\n*                                  is to be referenced.\n*\n*           Unchanged on exit.\n*\n*  TRANS    (input) CHARACTER*1\n*           On entry, TRANS specifies the operation to be performed as\n*           follows:\n*\n*              TRANS = 'N' or 'n'   C := alpha*A*A' + beta*C.\n*\n*              TRANS = 'T' or 't'   C := alpha*A'*A + beta*C.\n*\n*           Unchanged on exit.\n*\n*  N       (input) INTEGER\n*           On entry, N specifies the order of the matrix C. N must be\n*           at least zero.\n*           Unchanged on exit.\n*\n*  K       (input) INTEGER\n*           On entry with TRANS = 'N' or 'n', K specifies the number\n*           of  columns of the matrix A, and on entry with TRANS = 'T'\n*           or 't', K specifies the number of rows of the matrix A. K\n*           must be at least zero.\n*           Unchanged on exit.\n*\n*  ALPHA   (input) REAL\n*           On entry, ALPHA specifies the scalar alpha.\n*           Unchanged on exit.\n*\n*  A       (input) REAL array of DIMENSION (LDA,ka)\n*           where KA\n*           is K  when TRANS = 'N' or 'n', and is N otherwise. Before\n*           entry with TRANS = 'N' or 'n', the leading N--by--K part of\n*           the array A must contain the matrix A, otherwise the leading\n*           K--by--N part of the array A must contain the matrix A.\n*           Unchanged on exit.\n*\n*  LDA     (input) INTEGER\n*           On entry, LDA specifies the first dimension of A as declared\n*           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'\n*           then  LDA must be at least  max( 1, n ), otherwise  LDA must\n*           be at least  max( 1, k ).\n*           Unchanged on exit.\n*\n*  BETA    (input) REAL\n*           On entry, BETA specifies the scalar beta.\n*           Unchanged on exit.\n*\n*\n*  C       (input/output) REAL array, dimension (NT)\n*           NT = N*(N+1)/2. On entry, the symmetric matrix C in RFP\n*           Format. RFP Format is described by TRANSR, UPLO and N.\n*\n*  Arguments\n*  ==========\n*\n*     ..\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  c = NumRu::Lapack.ssfrk( transr, uplo, trans, n, k, alpha, a, beta, c, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 9)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 9)", argc);
  rblapack_transr = argv[0];
  rblapack_uplo = argv[1];
  rblapack_trans = argv[2];
  rblapack_n = argv[3];
  rblapack_k = argv[4];
  rblapack_alpha = argv[5];
  rblapack_a = argv[6];
  rblapack_beta = argv[7];
  rblapack_c = argv[8];
  if (rb_options != Qnil) {
  }

  k = NUM2INT(rblapack_k);
  transr = StringValueCStr(rblapack_transr)[0];
  if (!NA_IsNArray(rblapack_c))
    rb_raise(rb_eArgError, "c (9th argument) must be NArray");
  if (NA_RANK(rblapack_c) != 1)
    rb_raise(rb_eArgError, "rank of c (9th argument) must be %d", 1);
  nt = NA_SHAPE0(rblapack_c);
  if (NA_TYPE(rblapack_c) != NA_SFLOAT)
    rblapack_c = na_change_type(rblapack_c, NA_SFLOAT);
  c = NA_PTR_TYPE(rblapack_c, real*);
  n = NUM2INT(rblapack_n);
  uplo = StringValueCStr(rblapack_uplo)[0];
  alpha = (real)NUM2DBL(rblapack_alpha);
  trans = StringValueCStr(rblapack_trans)[0];
  beta = (real)NUM2DBL(rblapack_beta);
  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (7th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (7th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_a) != ((lsame_(&trans,"N") || lsame_(&trans,"n")) ? k : n))
    rb_raise(rb_eRuntimeError, "shape 1 of a must be %d", (lsame_(&trans,"N") || lsame_(&trans,"n")) ? k : n);
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_SFLOAT)
    rblapack_a = na_change_type(rblapack_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rblapack_a, real*);
  {
    int shape[1];
    shape[0] = nt;
    rblapack_c_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  c_out__ = NA_PTR_TYPE(rblapack_c_out__, real*);
  MEMCPY(c_out__, c, real, NA_TOTAL(rblapack_c));
  rblapack_c = rblapack_c_out__;
  c = c_out__;

  ssfrk_(&transr, &uplo, &trans, &n, &k, &alpha, a, &lda, &beta, c);

  return rblapack_c;
}

void
init_lapack_ssfrk(VALUE mLapack){
  rb_define_module_function(mLapack, "ssfrk", rblapack_ssfrk, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
