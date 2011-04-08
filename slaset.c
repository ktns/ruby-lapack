#include "rb_lapack.h"

extern VOID slaset_(char *uplo, integer *m, integer *n, real *alpha, real *beta, real *a, integer *lda);

static VALUE sHelp, sUsage;

static VALUE
rblapack_slaset(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_uplo;
  char uplo; 
  VALUE rblapack_m;
  integer m; 
  VALUE rblapack_alpha;
  real alpha; 
  VALUE rblapack_beta;
  real beta; 
  VALUE rblapack_a;
  real *a; 
  VALUE rblapack_a_out__;
  real *a_out__;

  integer lda;
  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  a = NumRu::Lapack.slaset( uplo, m, alpha, beta, a, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLASET( UPLO, M, N, ALPHA, BETA, A, LDA )\n\n*  Purpose\n*  =======\n*\n*  SLASET initializes an m-by-n matrix A to BETA on the diagonal and\n*  ALPHA on the offdiagonals.\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO    (input) CHARACTER*1\n*          Specifies the part of the matrix A to be set.\n*          = 'U':      Upper triangular part is set; the strictly lower\n*                      triangular part of A is not changed.\n*          = 'L':      Lower triangular part is set; the strictly upper\n*                      triangular part of A is not changed.\n*          Otherwise:  All of the matrix A is set.\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix A.  M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix A.  N >= 0.\n*\n*  ALPHA   (input) REAL\n*          The constant to which the offdiagonal elements are to be set.\n*\n*  BETA    (input) REAL\n*          The constant to which the diagonal elements are to be set.\n*\n*  A       (input/output) REAL array, dimension (LDA,N)\n*          On exit, the leading m-by-n submatrix of A is set as follows:\n*\n*          if UPLO = 'U', A(i,j) = ALPHA, 1<=i<=j-1, 1<=j<=n,\n*          if UPLO = 'L', A(i,j) = ALPHA, j+1<=i<=m, 1<=j<=n,\n*          otherwise,     A(i,j) = ALPHA, 1<=i<=m, 1<=j<=n, i.ne.j,\n*\n*          and, for all UPLO, A(i,i) = BETA, 1<=i<=min(m,n).\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,M).\n*\n\n* =====================================================================\n*\n*     .. Local Scalars ..\n      INTEGER            I, J\n*     ..\n*     .. External Functions ..\n      LOGICAL            LSAME\n      EXTERNAL           LSAME\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          MIN\n*     ..\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  a = NumRu::Lapack.slaset( uplo, m, alpha, beta, a, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rblapack_uplo = argv[0];
  rblapack_m = argv[1];
  rblapack_alpha = argv[2];
  rblapack_beta = argv[3];
  rblapack_a = argv[4];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (5th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (5th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_a);
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_SFLOAT)
    rblapack_a = na_change_type(rblapack_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rblapack_a, real*);
  m = NUM2INT(rblapack_m);
  beta = (real)NUM2DBL(rblapack_beta);
  alpha = (real)NUM2DBL(rblapack_alpha);
  uplo = StringValueCStr(rblapack_uplo)[0];
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rblapack_a_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rblapack_a_out__, real*);
  MEMCPY(a_out__, a, real, NA_TOTAL(rblapack_a));
  rblapack_a = rblapack_a_out__;
  a = a_out__;

  slaset_(&uplo, &m, &n, &alpha, &beta, a, &lda);

  return rblapack_a;
}

void
init_lapack_slaset(VALUE mLapack){
  rb_define_module_function(mLapack, "slaset", rblapack_slaset, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
