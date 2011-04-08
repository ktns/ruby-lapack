#include "rb_lapack.h"

extern VOID slatrs_(char *uplo, char *trans, char *diag, char *normin, integer *n, real *a, integer *lda, real *x, real *scale, real *cnorm, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_slatrs(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_uplo;
  char uplo; 
  VALUE rblapack_trans;
  char trans; 
  VALUE rblapack_diag;
  char diag; 
  VALUE rblapack_normin;
  char normin; 
  VALUE rblapack_a;
  real *a; 
  VALUE rblapack_x;
  real *x; 
  VALUE rblapack_cnorm;
  real *cnorm; 
  VALUE rblapack_scale;
  real scale; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_x_out__;
  real *x_out__;
  VALUE rblapack_cnorm_out__;
  real *cnorm_out__;

  integer lda;
  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  scale, info, x, cnorm = NumRu::Lapack.slatrs( uplo, trans, diag, normin, a, x, cnorm, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLATRS( UPLO, TRANS, DIAG, NORMIN, N, A, LDA, X, SCALE, CNORM, INFO )\n\n*  Purpose\n*  =======\n*\n*  SLATRS solves one of the triangular systems\n*\n*     A *x = s*b  or  A'*x = s*b\n*\n*  with scaling to prevent overflow.  Here A is an upper or lower\n*  triangular matrix, A' denotes the transpose of A, x and b are\n*  n-element vectors, and s is a scaling factor, usually less than\n*  or equal to 1, chosen so that the components of x will be less than\n*  the overflow threshold.  If the unscaled problem will not cause\n*  overflow, the Level 2 BLAS routine STRSV is called.  If the matrix A\n*  is singular (A(j,j) = 0 for some j), then s is set to 0 and a\n*  non-trivial solution to A*x = 0 is returned.\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO    (input) CHARACTER*1\n*          Specifies whether the matrix A is upper or lower triangular.\n*          = 'U':  Upper triangular\n*          = 'L':  Lower triangular\n*\n*  TRANS   (input) CHARACTER*1\n*          Specifies the operation applied to A.\n*          = 'N':  Solve A * x = s*b  (No transpose)\n*          = 'T':  Solve A'* x = s*b  (Transpose)\n*          = 'C':  Solve A'* x = s*b  (Conjugate transpose = Transpose)\n*\n*  DIAG    (input) CHARACTER*1\n*          Specifies whether or not the matrix A is unit triangular.\n*          = 'N':  Non-unit triangular\n*          = 'U':  Unit triangular\n*\n*  NORMIN  (input) CHARACTER*1\n*          Specifies whether CNORM has been set or not.\n*          = 'Y':  CNORM contains the column norms on entry\n*          = 'N':  CNORM is not set on entry.  On exit, the norms will\n*                  be computed and stored in CNORM.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  A       (input) REAL array, dimension (LDA,N)\n*          The triangular matrix A.  If UPLO = 'U', the leading n by n\n*          upper triangular part of the array A contains the upper\n*          triangular matrix, and the strictly lower triangular part of\n*          A is not referenced.  If UPLO = 'L', the leading n by n lower\n*          triangular part of the array A contains the lower triangular\n*          matrix, and the strictly upper triangular part of A is not\n*          referenced.  If DIAG = 'U', the diagonal elements of A are\n*          also not referenced and are assumed to be 1.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max (1,N).\n*\n*  X       (input/output) REAL array, dimension (N)\n*          On entry, the right hand side b of the triangular system.\n*          On exit, X is overwritten by the solution vector x.\n*\n*  SCALE   (output) REAL\n*          The scaling factor s for the triangular system\n*             A * x = s*b  or  A'* x = s*b.\n*          If SCALE = 0, the matrix A is singular or badly scaled, and\n*          the vector x is an exact or approximate solution to A*x = 0.\n*\n*  CNORM   (input or output) REAL array, dimension (N)\n*\n*          If NORMIN = 'Y', CNORM is an input argument and CNORM(j)\n*          contains the norm of the off-diagonal part of the j-th column\n*          of A.  If TRANS = 'N', CNORM(j) must be greater than or equal\n*          to the infinity-norm, and if TRANS = 'T' or 'C', CNORM(j)\n*          must be greater than or equal to the 1-norm.\n*\n*          If NORMIN = 'N', CNORM is an output argument and CNORM(j)\n*          returns the 1-norm of the offdiagonal part of the j-th column\n*          of A.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -k, the k-th argument had an illegal value\n*\n\n*  Further Details\n*  ======= =======\n*\n*  A rough bound on x is computed; if that is less than overflow, STRSV\n*  is called, otherwise, specific code is used which checks for possible\n*  overflow or divide-by-zero at every operation.\n*\n*  A columnwise scheme is used for solving A*x = b.  The basic algorithm\n*  if A is lower triangular is\n*\n*       x[1:n] := b[1:n]\n*       for j = 1, ..., n\n*            x(j) := x(j) / A(j,j)\n*            x[j+1:n] := x[j+1:n] - x(j) * A[j+1:n,j]\n*       end\n*\n*  Define bounds on the components of x after j iterations of the loop:\n*     M(j) = bound on x[1:j]\n*     G(j) = bound on x[j+1:n]\n*  Initially, let M(0) = 0 and G(0) = max{x(i), i=1,...,n}.\n*\n*  Then for iteration j+1 we have\n*     M(j+1) <= G(j) / | A(j+1,j+1) |\n*     G(j+1) <= G(j) + M(j+1) * | A[j+2:n,j+1] |\n*            <= G(j) ( 1 + CNORM(j+1) / | A(j+1,j+1) | )\n*\n*  where CNORM(j+1) is greater than or equal to the infinity-norm of\n*  column j+1 of A, not counting the diagonal.  Hence\n*\n*     G(j) <= G(0) product ( 1 + CNORM(i) / | A(i,i) | )\n*                  1<=i<=j\n*  and\n*\n*     |x(j)| <= ( G(0) / |A(j,j)| ) product ( 1 + CNORM(i) / |A(i,i)| )\n*                                   1<=i< j\n*\n*  Since |x(j)| <= M(j), we use the Level 2 BLAS routine STRSV if the\n*  reciprocal of the largest M(j), j=1,..,n, is larger than\n*  max(underflow, 1/overflow).\n*\n*  The bound on x(j) is also used to determine when a step in the\n*  columnwise method can be performed without fear of overflow.  If\n*  the computed bound is greater than a large constant, x is scaled to\n*  prevent overflow, but if the bound overflows, x is set to 0, x(j) to\n*  1, and scale to 0, and a non-trivial solution to A*x = 0 is found.\n*\n*  Similarly, a row-wise scheme is used to solve A'*x = b.  The basic\n*  algorithm for A upper triangular is\n*\n*       for j = 1, ..., n\n*            x(j) := ( b(j) - A[1:j-1,j]' * x[1:j-1] ) / A(j,j)\n*       end\n*\n*  We simultaneously compute two bounds\n*       G(j) = bound on ( b(i) - A[1:i-1,i]' * x[1:i-1] ), 1<=i<=j\n*       M(j) = bound on x(i), 1<=i<=j\n*\n*  The initial values are G(0) = 0, M(0) = max{b(i), i=1,..,n}, and we\n*  add the constraint G(j) >= G(j-1) and M(j) >= M(j-1) for j >= 1.\n*  Then the bound on x(j) is\n*\n*       M(j) <= M(j-1) * ( 1 + CNORM(j) ) / | A(j,j) |\n*\n*            <= M(0) * product ( ( 1 + CNORM(i) ) / |A(i,i)| )\n*                      1<=i<=j\n*\n*  and we can safely call STRSV if 1/M(n) and 1/G(n) are both greater\n*  than max(underflow, 1/overflow).\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  scale, info, x, cnorm = NumRu::Lapack.slatrs( uplo, trans, diag, normin, a, x, cnorm, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 7)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 7)", argc);
  rblapack_uplo = argv[0];
  rblapack_trans = argv[1];
  rblapack_diag = argv[2];
  rblapack_normin = argv[3];
  rblapack_a = argv[4];
  rblapack_x = argv[5];
  rblapack_cnorm = argv[6];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_cnorm))
    rb_raise(rb_eArgError, "cnorm (7th argument) must be NArray");
  if (NA_RANK(rblapack_cnorm) != 1)
    rb_raise(rb_eArgError, "rank of cnorm (7th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_cnorm);
  if (NA_TYPE(rblapack_cnorm) != NA_SFLOAT)
    rblapack_cnorm = na_change_type(rblapack_cnorm, NA_SFLOAT);
  cnorm = NA_PTR_TYPE(rblapack_cnorm, real*);
  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (5th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (5th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_a) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of a must be the same as shape 0 of cnorm");
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_SFLOAT)
    rblapack_a = na_change_type(rblapack_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rblapack_a, real*);
  if (!NA_IsNArray(rblapack_x))
    rb_raise(rb_eArgError, "x (6th argument) must be NArray");
  if (NA_RANK(rblapack_x) != 1)
    rb_raise(rb_eArgError, "rank of x (6th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_x) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of x must be the same as shape 0 of cnorm");
  if (NA_TYPE(rblapack_x) != NA_SFLOAT)
    rblapack_x = na_change_type(rblapack_x, NA_SFLOAT);
  x = NA_PTR_TYPE(rblapack_x, real*);
  diag = StringValueCStr(rblapack_diag)[0];
  normin = StringValueCStr(rblapack_normin)[0];
  trans = StringValueCStr(rblapack_trans)[0];
  uplo = StringValueCStr(rblapack_uplo)[0];
  {
    int shape[1];
    shape[0] = n;
    rblapack_x_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  x_out__ = NA_PTR_TYPE(rblapack_x_out__, real*);
  MEMCPY(x_out__, x, real, NA_TOTAL(rblapack_x));
  rblapack_x = rblapack_x_out__;
  x = x_out__;
  {
    int shape[1];
    shape[0] = n;
    rblapack_cnorm_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  cnorm_out__ = NA_PTR_TYPE(rblapack_cnorm_out__, real*);
  MEMCPY(cnorm_out__, cnorm, real, NA_TOTAL(rblapack_cnorm));
  rblapack_cnorm = rblapack_cnorm_out__;
  cnorm = cnorm_out__;

  slatrs_(&uplo, &trans, &diag, &normin, &n, a, &lda, x, &scale, cnorm, &info);

  rblapack_scale = rb_float_new((double)scale);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(4, rblapack_scale, rblapack_info, rblapack_x, rblapack_cnorm);
}

void
init_lapack_slatrs(VALUE mLapack){
  rb_define_module_function(mLapack, "slatrs", rblapack_slatrs, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
