#include "rb_lapack.h"

extern VOID zgesc2_(integer *n, doublecomplex *a, integer *lda, doublecomplex *rhs, integer *ipiv, integer *jpiv, doublereal *scale);

static VALUE sHelp, sUsage;

static VALUE
rblapack_zgesc2(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_a;
  doublecomplex *a; 
  VALUE rblapack_rhs;
  doublecomplex *rhs; 
  VALUE rblapack_ipiv;
  integer *ipiv; 
  VALUE rblapack_jpiv;
  integer *jpiv; 
  VALUE rblapack_scale;
  doublereal scale; 
  VALUE rblapack_rhs_out__;
  doublecomplex *rhs_out__;

  integer lda;
  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  scale, rhs = NumRu::Lapack.zgesc2( a, rhs, ipiv, jpiv, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZGESC2( N, A, LDA, RHS, IPIV, JPIV, SCALE )\n\n*  Purpose\n*  =======\n*\n*  ZGESC2 solves a system of linear equations\n*\n*            A * X = scale* RHS\n*\n*  with a general N-by-N matrix A using the LU factorization with\n*  complete pivoting computed by ZGETC2.\n*\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix A.\n*\n*  A       (input) COMPLEX*16 array, dimension (LDA, N)\n*          On entry, the  LU part of the factorization of the n-by-n\n*          matrix A computed by ZGETC2:  A = P * L * U * Q\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1, N).\n*\n*  RHS     (input/output) COMPLEX*16 array, dimension N.\n*          On entry, the right hand side vector b.\n*          On exit, the solution vector X.\n*\n*  IPIV    (input) INTEGER array, dimension (N).\n*          The pivot indices; for 1 <= i <= N, row i of the\n*          matrix has been interchanged with row IPIV(i).\n*\n*  JPIV    (input) INTEGER array, dimension (N).\n*          The pivot indices; for 1 <= j <= N, column j of the\n*          matrix has been interchanged with column JPIV(j).\n*\n*  SCALE    (output) DOUBLE PRECISION\n*           On exit, SCALE contains the scale factor. SCALE is chosen\n*           0 <= SCALE <= 1 to prevent owerflow in the solution.\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Bo Kagstrom and Peter Poromaa, Department of Computing Science,\n*     Umea University, S-901 87 Umea, Sweden.\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  scale, rhs = NumRu::Lapack.zgesc2( a, rhs, ipiv, jpiv, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rblapack_a = argv[0];
  rblapack_rhs = argv[1];
  rblapack_ipiv = argv[2];
  rblapack_jpiv = argv[3];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_ipiv))
    rb_raise(rb_eArgError, "ipiv (3th argument) must be NArray");
  if (NA_RANK(rblapack_ipiv) != 1)
    rb_raise(rb_eArgError, "rank of ipiv (3th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_ipiv);
  if (NA_TYPE(rblapack_ipiv) != NA_LINT)
    rblapack_ipiv = na_change_type(rblapack_ipiv, NA_LINT);
  ipiv = NA_PTR_TYPE(rblapack_ipiv, integer*);
  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (1th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (1th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_a) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of a must be the same as shape 0 of ipiv");
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_DCOMPLEX)
    rblapack_a = na_change_type(rblapack_a, NA_DCOMPLEX);
  a = NA_PTR_TYPE(rblapack_a, doublecomplex*);
  if (!NA_IsNArray(rblapack_rhs))
    rb_raise(rb_eArgError, "rhs (2th argument) must be NArray");
  if (NA_RANK(rblapack_rhs) != 1)
    rb_raise(rb_eArgError, "rank of rhs (2th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_rhs) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of rhs must be the same as shape 0 of ipiv");
  if (NA_TYPE(rblapack_rhs) != NA_DCOMPLEX)
    rblapack_rhs = na_change_type(rblapack_rhs, NA_DCOMPLEX);
  rhs = NA_PTR_TYPE(rblapack_rhs, doublecomplex*);
  if (!NA_IsNArray(rblapack_jpiv))
    rb_raise(rb_eArgError, "jpiv (4th argument) must be NArray");
  if (NA_RANK(rblapack_jpiv) != 1)
    rb_raise(rb_eArgError, "rank of jpiv (4th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_jpiv) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of jpiv must be the same as shape 0 of ipiv");
  if (NA_TYPE(rblapack_jpiv) != NA_LINT)
    rblapack_jpiv = na_change_type(rblapack_jpiv, NA_LINT);
  jpiv = NA_PTR_TYPE(rblapack_jpiv, integer*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_rhs_out__ = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  rhs_out__ = NA_PTR_TYPE(rblapack_rhs_out__, doublecomplex*);
  MEMCPY(rhs_out__, rhs, doublecomplex, NA_TOTAL(rblapack_rhs));
  rblapack_rhs = rblapack_rhs_out__;
  rhs = rhs_out__;

  zgesc2_(&n, a, &lda, rhs, ipiv, jpiv, &scale);

  rblapack_scale = rb_float_new((double)scale);
  return rb_ary_new3(2, rblapack_scale, rblapack_rhs);
}

void
init_lapack_zgesc2(VALUE mLapack){
  rb_define_module_function(mLapack, "zgesc2", rblapack_zgesc2, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
