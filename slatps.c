#include "rb_lapack.h"

extern VOID slatps_(char *uplo, char *trans, char *diag, char *normin, integer *n, real *ap, real *x, real *scale, real *cnorm, integer *info);

static VALUE
rb_slatps(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_trans;
  char trans; 
  VALUE rb_diag;
  char diag; 
  VALUE rb_normin;
  char normin; 
  VALUE rb_ap;
  real *ap; 
  VALUE rb_x;
  real *x; 
  VALUE rb_cnorm;
  real *cnorm; 
  VALUE rb_scale;
  real scale; 
  VALUE rb_info;
  integer info; 
  VALUE rb_x_out__;
  real *x_out__;
  VALUE rb_cnorm_out__;
  real *cnorm_out__;

  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  scale, info, x, cnorm = NumRu::Lapack.slatps( uplo, trans, diag, normin, ap, x, cnorm)\n    or\n  NumRu::Lapack.slatps  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLATPS( UPLO, TRANS, DIAG, NORMIN, N, AP, X, SCALE, CNORM, INFO )\n\n*  Purpose\n*  =======\n*\n*  SLATPS solves one of the triangular systems\n*\n*     A *x = s*b  or  A'*x = s*b\n*\n*  with scaling to prevent overflow, where A is an upper or lower\n*  triangular matrix stored in packed form.  Here A' denotes the\n*  transpose of A, x and b are n-element vectors, and s is a scaling\n*  factor, usually less than or equal to 1, chosen so that the\n*  components of x will be less than the overflow threshold.  If the\n*  unscaled problem will not cause overflow, the Level 2 BLAS routine\n*  STPSV is called. If the matrix A is singular (A(j,j) = 0 for some j),\n*  then s is set to 0 and a non-trivial solution to A*x = 0 is returned.\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO    (input) CHARACTER*1\n*          Specifies whether the matrix A is upper or lower triangular.\n*          = 'U':  Upper triangular\n*          = 'L':  Lower triangular\n*\n*  TRANS   (input) CHARACTER*1\n*          Specifies the operation applied to A.\n*          = 'N':  Solve A * x = s*b  (No transpose)\n*          = 'T':  Solve A'* x = s*b  (Transpose)\n*          = 'C':  Solve A'* x = s*b  (Conjugate transpose = Transpose)\n*\n*  DIAG    (input) CHARACTER*1\n*          Specifies whether or not the matrix A is unit triangular.\n*          = 'N':  Non-unit triangular\n*          = 'U':  Unit triangular\n*\n*  NORMIN  (input) CHARACTER*1\n*          Specifies whether CNORM has been set or not.\n*          = 'Y':  CNORM contains the column norms on entry\n*          = 'N':  CNORM is not set on entry.  On exit, the norms will\n*                  be computed and stored in CNORM.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  AP      (input) REAL array, dimension (N*(N+1)/2)\n*          The upper or lower triangular matrix A, packed columnwise in\n*          a linear array.  The j-th column of A is stored in the array\n*          AP as follows:\n*          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;\n*          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.\n*\n*  X       (input/output) REAL array, dimension (N)\n*          On entry, the right hand side b of the triangular system.\n*          On exit, X is overwritten by the solution vector x.\n*\n*  SCALE   (output) REAL\n*          The scaling factor s for the triangular system\n*             A * x = s*b  or  A'* x = s*b.\n*          If SCALE = 0, the matrix A is singular or badly scaled, and\n*          the vector x is an exact or approximate solution to A*x = 0.\n*\n*  CNORM   (input or output) REAL array, dimension (N)\n*\n*          If NORMIN = 'Y', CNORM is an input argument and CNORM(j)\n*          contains the norm of the off-diagonal part of the j-th column\n*          of A.  If TRANS = 'N', CNORM(j) must be greater than or equal\n*          to the infinity-norm, and if TRANS = 'T' or 'C', CNORM(j)\n*          must be greater than or equal to the 1-norm.\n*\n*          If NORMIN = 'N', CNORM is an output argument and CNORM(j)\n*          returns the 1-norm of the offdiagonal part of the j-th column\n*          of A.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -k, the k-th argument had an illegal value\n*\n\n*  Further Details\n*  ======= =======\n*\n*  A rough bound on x is computed; if that is less than overflow, STPSV\n*  is called, otherwise, specific code is used which checks for possible\n*  overflow or divide-by-zero at every operation.\n*\n*  A columnwise scheme is used for solving A*x = b.  The basic algorithm\n*  if A is lower triangular is\n*\n*       x[1:n] := b[1:n]\n*       for j = 1, ..., n\n*            x(j) := x(j) / A(j,j)\n*            x[j+1:n] := x[j+1:n] - x(j) * A[j+1:n,j]\n*       end\n*\n*  Define bounds on the components of x after j iterations of the loop:\n*     M(j) = bound on x[1:j]\n*     G(j) = bound on x[j+1:n]\n*  Initially, let M(0) = 0 and G(0) = max{x(i), i=1,...,n}.\n*\n*  Then for iteration j+1 we have\n*     M(j+1) <= G(j) / | A(j+1,j+1) |\n*     G(j+1) <= G(j) + M(j+1) * | A[j+2:n,j+1] |\n*            <= G(j) ( 1 + CNORM(j+1) / | A(j+1,j+1) | )\n*\n*  where CNORM(j+1) is greater than or equal to the infinity-norm of\n*  column j+1 of A, not counting the diagonal.  Hence\n*\n*     G(j) <= G(0) product ( 1 + CNORM(i) / | A(i,i) | )\n*                  1<=i<=j\n*  and\n*\n*     |x(j)| <= ( G(0) / |A(j,j)| ) product ( 1 + CNORM(i) / |A(i,i)| )\n*                                   1<=i< j\n*\n*  Since |x(j)| <= M(j), we use the Level 2 BLAS routine STPSV if the\n*  reciprocal of the largest M(j), j=1,..,n, is larger than\n*  max(underflow, 1/overflow).\n*\n*  The bound on x(j) is also used to determine when a step in the\n*  columnwise method can be performed without fear of overflow.  If\n*  the computed bound is greater than a large constant, x is scaled to\n*  prevent overflow, but if the bound overflows, x is set to 0, x(j) to\n*  1, and scale to 0, and a non-trivial solution to A*x = 0 is found.\n*\n*  Similarly, a row-wise scheme is used to solve A'*x = b.  The basic\n*  algorithm for A upper triangular is\n*\n*       for j = 1, ..., n\n*            x(j) := ( b(j) - A[1:j-1,j]' * x[1:j-1] ) / A(j,j)\n*       end\n*\n*  We simultaneously compute two bounds\n*       G(j) = bound on ( b(i) - A[1:i-1,i]' * x[1:i-1] ), 1<=i<=j\n*       M(j) = bound on x(i), 1<=i<=j\n*\n*  The initial values are G(0) = 0, M(0) = max{b(i), i=1,..,n}, and we\n*  add the constraint G(j) >= G(j-1) and M(j) >= M(j-1) for j >= 1.\n*  Then the bound on x(j) is\n*\n*       M(j) <= M(j-1) * ( 1 + CNORM(j) ) / | A(j,j) |\n*\n*            <= M(0) * product ( ( 1 + CNORM(i) ) / |A(i,i)| )\n*                      1<=i<=j\n*\n*  and we can safely call STPSV if 1/M(n) and 1/G(n) are both greater\n*  than max(underflow, 1/overflow).\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 7)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 7)", argc);
  rb_uplo = argv[0];
  rb_trans = argv[1];
  rb_diag = argv[2];
  rb_normin = argv[3];
  rb_ap = argv[4];
  rb_x = argv[5];
  rb_cnorm = argv[6];

  if (!NA_IsNArray(rb_cnorm))
    rb_raise(rb_eArgError, "cnorm (7th argument) must be NArray");
  if (NA_RANK(rb_cnorm) != 1)
    rb_raise(rb_eArgError, "rank of cnorm (7th argument) must be %d", 1);
  n = NA_SHAPE0(rb_cnorm);
  if (NA_TYPE(rb_cnorm) != NA_SFLOAT)
    rb_cnorm = na_change_type(rb_cnorm, NA_SFLOAT);
  cnorm = NA_PTR_TYPE(rb_cnorm, real*);
  if (!NA_IsNArray(rb_x))
    rb_raise(rb_eArgError, "x (6th argument) must be NArray");
  if (NA_RANK(rb_x) != 1)
    rb_raise(rb_eArgError, "rank of x (6th argument) must be %d", 1);
  if (NA_SHAPE0(rb_x) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of x must be the same as shape 0 of cnorm");
  if (NA_TYPE(rb_x) != NA_SFLOAT)
    rb_x = na_change_type(rb_x, NA_SFLOAT);
  x = NA_PTR_TYPE(rb_x, real*);
  diag = StringValueCStr(rb_diag)[0];
  trans = StringValueCStr(rb_trans)[0];
  normin = StringValueCStr(rb_normin)[0];
  uplo = StringValueCStr(rb_uplo)[0];
  if (!NA_IsNArray(rb_ap))
    rb_raise(rb_eArgError, "ap (5th argument) must be NArray");
  if (NA_RANK(rb_ap) != 1)
    rb_raise(rb_eArgError, "rank of ap (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_ap) != (n*(n+1)/2))
    rb_raise(rb_eRuntimeError, "shape 0 of ap must be %d", n*(n+1)/2);
  if (NA_TYPE(rb_ap) != NA_SFLOAT)
    rb_ap = na_change_type(rb_ap, NA_SFLOAT);
  ap = NA_PTR_TYPE(rb_ap, real*);
  {
    int shape[1];
    shape[0] = n;
    rb_x_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  x_out__ = NA_PTR_TYPE(rb_x_out__, real*);
  MEMCPY(x_out__, x, real, NA_TOTAL(rb_x));
  rb_x = rb_x_out__;
  x = x_out__;
  {
    int shape[1];
    shape[0] = n;
    rb_cnorm_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  cnorm_out__ = NA_PTR_TYPE(rb_cnorm_out__, real*);
  MEMCPY(cnorm_out__, cnorm, real, NA_TOTAL(rb_cnorm));
  rb_cnorm = rb_cnorm_out__;
  cnorm = cnorm_out__;

  slatps_(&uplo, &trans, &diag, &normin, &n, ap, x, &scale, cnorm, &info);

  rb_scale = rb_float_new((double)scale);
  rb_info = INT2NUM(info);
  return rb_ary_new3(4, rb_scale, rb_info, rb_x, rb_cnorm);
}

void
init_lapack_slatps(VALUE mLapack){
  rb_define_module_function(mLapack, "slatps", rb_slatps, -1);
}
