#include "rb_lapack.h"

extern VOID strsyl_(char *trana, char *tranb, integer *isgn, integer *m, integer *n, real *a, integer *lda, real *b, integer *ldb, real *c, integer *ldc, real *scale, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_strsyl(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_trana;
  char trana; 
  VALUE rblapack_tranb;
  char tranb; 
  VALUE rblapack_isgn;
  integer isgn; 
  VALUE rblapack_a;
  real *a; 
  VALUE rblapack_b;
  real *b; 
  VALUE rblapack_c;
  real *c; 
  VALUE rblapack_scale;
  real scale; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_c_out__;
  real *c_out__;

  integer lda;
  integer m;
  integer ldb;
  integer n;
  integer ldc;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  scale, info, c = NumRu::Lapack.strsyl( trana, tranb, isgn, a, b, c, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE STRSYL( TRANA, TRANB, ISGN, M, N, A, LDA, B, LDB, C, LDC, SCALE, INFO )\n\n*  Purpose\n*  =======\n*\n*  STRSYL solves the real Sylvester matrix equation:\n*\n*     op(A)*X + X*op(B) = scale*C or\n*     op(A)*X - X*op(B) = scale*C,\n*\n*  where op(A) = A or A**T, and  A and B are both upper quasi-\n*  triangular. A is M-by-M and B is N-by-N; the right hand side C and\n*  the solution X are M-by-N; and scale is an output scale factor, set\n*  <= 1 to avoid overflow in X.\n*\n*  A and B must be in Schur canonical form (as returned by SHSEQR), that\n*  is, block upper triangular with 1-by-1 and 2-by-2 diagonal blocks;\n*  each 2-by-2 diagonal block has its diagonal elements equal and its\n*  off-diagonal elements of opposite sign.\n*\n\n*  Arguments\n*  =========\n*\n*  TRANA   (input) CHARACTER*1\n*          Specifies the option op(A):\n*          = 'N': op(A) = A    (No transpose)\n*          = 'T': op(A) = A**T (Transpose)\n*          = 'C': op(A) = A**H (Conjugate transpose = Transpose)\n*\n*  TRANB   (input) CHARACTER*1\n*          Specifies the option op(B):\n*          = 'N': op(B) = B    (No transpose)\n*          = 'T': op(B) = B**T (Transpose)\n*          = 'C': op(B) = B**H (Conjugate transpose = Transpose)\n*\n*  ISGN    (input) INTEGER\n*          Specifies the sign in the equation:\n*          = +1: solve op(A)*X + X*op(B) = scale*C\n*          = -1: solve op(A)*X - X*op(B) = scale*C\n*\n*  M       (input) INTEGER\n*          The order of the matrix A, and the number of rows in the\n*          matrices X and C. M >= 0.\n*\n*  N       (input) INTEGER\n*          The order of the matrix B, and the number of columns in the\n*          matrices X and C. N >= 0.\n*\n*  A       (input) REAL array, dimension (LDA,M)\n*          The upper quasi-triangular matrix A, in Schur canonical form.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A. LDA >= max(1,M).\n*\n*  B       (input) REAL array, dimension (LDB,N)\n*          The upper quasi-triangular matrix B, in Schur canonical form.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the array B. LDB >= max(1,N).\n*\n*  C       (input/output) REAL array, dimension (LDC,N)\n*          On entry, the M-by-N right hand side matrix C.\n*          On exit, C is overwritten by the solution matrix X.\n*\n*  LDC     (input) INTEGER\n*          The leading dimension of the array C. LDC >= max(1,M)\n*\n*  SCALE   (output) REAL\n*          The scale factor, scale, set <= 1 to avoid overflow in X.\n*\n*  INFO    (output) INTEGER\n*          = 0: successful exit\n*          < 0: if INFO = -i, the i-th argument had an illegal value\n*          = 1: A and B have common or very close eigenvalues; perturbed\n*               values were used to solve the equation (but the matrices\n*               A and B are unchanged).\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  scale, info, c = NumRu::Lapack.strsyl( trana, tranb, isgn, a, b, c, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rblapack_trana = argv[0];
  rblapack_tranb = argv[1];
  rblapack_isgn = argv[2];
  rblapack_a = argv[3];
  rblapack_b = argv[4];
  rblapack_c = argv[5];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (4th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (4th argument) must be %d", 2);
  m = NA_SHAPE1(rblapack_a);
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_SFLOAT)
    rblapack_a = na_change_type(rblapack_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rblapack_a, real*);
  if (!NA_IsNArray(rblapack_b))
    rb_raise(rb_eArgError, "b (5th argument) must be NArray");
  if (NA_RANK(rblapack_b) != 2)
    rb_raise(rb_eArgError, "rank of b (5th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_b);
  ldb = NA_SHAPE0(rblapack_b);
  if (NA_TYPE(rblapack_b) != NA_SFLOAT)
    rblapack_b = na_change_type(rblapack_b, NA_SFLOAT);
  b = NA_PTR_TYPE(rblapack_b, real*);
  trana = StringValueCStr(rblapack_trana)[0];
  if (!NA_IsNArray(rblapack_c))
    rb_raise(rb_eArgError, "c (6th argument) must be NArray");
  if (NA_RANK(rblapack_c) != 2)
    rb_raise(rb_eArgError, "rank of c (6th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_c) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of c must be the same as shape 1 of b");
  ldc = NA_SHAPE0(rblapack_c);
  if (NA_TYPE(rblapack_c) != NA_SFLOAT)
    rblapack_c = na_change_type(rblapack_c, NA_SFLOAT);
  c = NA_PTR_TYPE(rblapack_c, real*);
  tranb = StringValueCStr(rblapack_tranb)[0];
  isgn = NUM2INT(rblapack_isgn);
  {
    int shape[2];
    shape[0] = ldc;
    shape[1] = n;
    rblapack_c_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  c_out__ = NA_PTR_TYPE(rblapack_c_out__, real*);
  MEMCPY(c_out__, c, real, NA_TOTAL(rblapack_c));
  rblapack_c = rblapack_c_out__;
  c = c_out__;

  strsyl_(&trana, &tranb, &isgn, &m, &n, a, &lda, b, &ldb, c, &ldc, &scale, &info);

  rblapack_scale = rb_float_new((double)scale);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(3, rblapack_scale, rblapack_info, rblapack_c);
}

void
init_lapack_strsyl(VALUE mLapack){
  rb_define_module_function(mLapack, "strsyl", rblapack_strsyl, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
