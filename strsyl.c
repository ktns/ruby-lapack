#include "rb_lapack.h"

static VALUE
rb_strsyl(int argc, VALUE *argv, VALUE self){
  VALUE rb_trana;
  char trana; 
  VALUE rb_tranb;
  char tranb; 
  VALUE rb_isgn;
  integer isgn; 
  VALUE rb_a;
  real *a; 
  VALUE rb_b;
  real *b; 
  VALUE rb_c;
  real *c; 
  VALUE rb_scale;
  real scale; 
  VALUE rb_info;
  integer info; 
  VALUE rb_c_out__;
  real *c_out__;

  integer lda;
  integer m;
  integer ldb;
  integer n;
  integer ldc;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  scale, info, c = NumRu::Lapack.strsyl( trana, tranb, isgn, a, b, c)\n    or\n  NumRu::Lapack.strsyl  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE STRSYL( TRANA, TRANB, ISGN, M, N, A, LDA, B, LDB, C, LDC, SCALE, INFO )\n\n*  Purpose\n*  =======\n*\n*  STRSYL solves the real Sylvester matrix equation:\n*\n*     op(A)*X + X*op(B) = scale*C or\n*     op(A)*X - X*op(B) = scale*C,\n*\n*  where op(A) = A or A**T, and  A and B are both upper quasi-\n*  triangular. A is M-by-M and B is N-by-N; the right hand side C and\n*  the solution X are M-by-N; and scale is an output scale factor, set\n*  <= 1 to avoid overflow in X.\n*\n*  A and B must be in Schur canonical form (as returned by SHSEQR), that\n*  is, block upper triangular with 1-by-1 and 2-by-2 diagonal blocks;\n*  each 2-by-2 diagonal block has its diagonal elements equal and its\n*  off-diagonal elements of opposite sign.\n*\n\n*  Arguments\n*  =========\n*\n*  TRANA   (input) CHARACTER*1\n*          Specifies the option op(A):\n*          = 'N': op(A) = A    (No transpose)\n*          = 'T': op(A) = A**T (Transpose)\n*          = 'C': op(A) = A**H (Conjugate transpose = Transpose)\n*\n*  TRANB   (input) CHARACTER*1\n*          Specifies the option op(B):\n*          = 'N': op(B) = B    (No transpose)\n*          = 'T': op(B) = B**T (Transpose)\n*          = 'C': op(B) = B**H (Conjugate transpose = Transpose)\n*\n*  ISGN    (input) INTEGER\n*          Specifies the sign in the equation:\n*          = +1: solve op(A)*X + X*op(B) = scale*C\n*          = -1: solve op(A)*X - X*op(B) = scale*C\n*\n*  M       (input) INTEGER\n*          The order of the matrix A, and the number of rows in the\n*          matrices X and C. M >= 0.\n*\n*  N       (input) INTEGER\n*          The order of the matrix B, and the number of columns in the\n*          matrices X and C. N >= 0.\n*\n*  A       (input) REAL array, dimension (LDA,M)\n*          The upper quasi-triangular matrix A, in Schur canonical form.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A. LDA >= max(1,M).\n*\n*  B       (input) REAL array, dimension (LDB,N)\n*          The upper quasi-triangular matrix B, in Schur canonical form.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B. LDB >= max(1,N).\n*\n*  C       (input/output) REAL array, dimension (LDC,N)\n*          On entry, the M-by-N right hand side matrix C.\n*          On exit, C is overwritten by the solution matrix X.\n*\n*  LDC     (input) INTEGER\n*          The leading dimension of the array C. LDC >= max(1,M)\n*\n*  SCALE   (output) REAL\n*          The scale factor, scale, set <= 1 to avoid overflow in X.\n*\n*  INFO    (output) INTEGER\n*          = 0: successful exit\n*          < 0: if INFO = -i, the i-th argument had an illegal value\n*          = 1: A and B have common or very close eigenvalues; perturbed\n*               values were used to solve the equation (but the matrices\n*               A and B are unchanged).\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rb_trana = argv[0];
  rb_tranb = argv[1];
  rb_isgn = argv[2];
  rb_a = argv[3];
  rb_b = argv[4];
  rb_c = argv[5];

  trana = StringValueCStr(rb_trana)[0];
  tranb = StringValueCStr(rb_tranb)[0];
  isgn = NUM2INT(rb_isgn);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (4th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (4th argument) must be %d", 2);
  lda = NA_SHAPE0(rb_a);
  m = NA_SHAPE1(rb_a);
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (5th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (5th argument) must be %d", 2);
  ldb = NA_SHAPE0(rb_b);
  n = NA_SHAPE1(rb_b);
  if (NA_TYPE(rb_b) != NA_SFLOAT)
    rb_b = na_change_type(rb_b, NA_SFLOAT);
  b = NA_PTR_TYPE(rb_b, real*);
  if (!NA_IsNArray(rb_c))
    rb_raise(rb_eArgError, "c (6th argument) must be NArray");
  if (NA_RANK(rb_c) != 2)
    rb_raise(rb_eArgError, "rank of c (6th argument) must be %d", 2);
  ldc = NA_SHAPE0(rb_c);
  if (NA_SHAPE1(rb_c) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of c must be the same as shape 1 of b");
  if (NA_TYPE(rb_c) != NA_SFLOAT)
    rb_c = na_change_type(rb_c, NA_SFLOAT);
  c = NA_PTR_TYPE(rb_c, real*);
  {
    int shape[2];
    shape[0] = ldc;
    shape[1] = n;
    rb_c_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  c_out__ = NA_PTR_TYPE(rb_c_out__, real*);
  MEMCPY(c_out__, c, real, NA_TOTAL(rb_c));
  rb_c = rb_c_out__;
  c = c_out__;

  strsyl_(&trana, &tranb, &isgn, &m, &n, a, &lda, b, &ldb, c, &ldc, &scale, &info);

  rb_scale = rb_float_new((double)scale);
  rb_info = INT2NUM(info);
  return rb_ary_new3(3, rb_scale, rb_info, rb_c);
}

void
init_lapack_strsyl(VALUE mLapack){
  rb_define_module_function(mLapack, "strsyl", rb_strsyl, -1);
}
