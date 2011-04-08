#include "rb_lapack.h"

extern VOID zsytri2_(char *uplo, integer *n, doublecomplex *a, integer *lda, integer *ipiv, doublecomplex *work, integer *lwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_zsytri2(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_uplo;
  char uplo; 
  VALUE rblapack_a;
  doublecomplex *a; 
  VALUE rblapack_ipiv;
  integer *ipiv; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_a_out__;
  doublecomplex *a_out__;
  integer c__1;
  integer nb;
  integer c__m1;
  doublecomplex *work;

  integer lda;
  integer n;
  integer lwork;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, a = NumRu::Lapack.zsytri2( uplo, a, ipiv, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZSYTRI2( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  ZSYTRI2 computes the inverse of a complex symmetric indefinite matrix\n*  A using the factorization A = U*D*U**T or A = L*D*L**T computed by\n*  ZSYTRF. ZSYTRI2 sets the LEADING DIMENSION of the workspace\n*  before calling ZSYTRI2X that actually computes the inverse.\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO    (input) CHARACTER*1\n*          Specifies whether the details of the factorization are stored\n*          as an upper or lower triangular matrix.\n*          = 'U':  Upper triangular, form is A = U*D*U**T;\n*          = 'L':  Lower triangular, form is A = L*D*L**T.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  A       (input/output) DOUBLE COMPLEX array, dimension (LDA,N)\n*          On entry, the NB diagonal matrix D and the multipliers\n*          used to obtain the factor U or L as computed by ZSYTRF.\n*\n*          On exit, if INFO = 0, the (symmetric) inverse of the original\n*          matrix.  If UPLO = 'U', the upper triangular part of the\n*          inverse is formed and the part of A below the diagonal is not\n*          referenced; if UPLO = 'L' the lower triangular part of the\n*          inverse is formed and the part of A above the diagonal is\n*          not referenced.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,N).\n*\n*  IPIV    (input) INTEGER array, dimension (N)\n*          Details of the interchanges and the NB structure of D\n*          as determined by ZSYTRF.\n*\n*  WORK    (workspace) DOUBLE COMPLEX array, dimension (N+NB+1)*(NB+3)\n*\n*  LWORK   (input) INTEGER\n*          The dimension of the array WORK.\n*          WORK is size >= (N+NB+1)*(NB+3)\n*          If LDWORK = -1, then a workspace query is assumed; the routine\n*           calculates:\n*              - the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array,\n*              - and no error message related to LDWORK is issued by XERBLA.\n*\n*  INFO    (output) INTEGER\n*          = 0: successful exit\n*          < 0: if INFO = -i, the i-th argument had an illegal value\n*          > 0: if INFO = i, D(i,i) = 0; the matrix is singular and its\n*               inverse could not be computed.\n*\n\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      LOGICAL            UPPER, LQUERY\n      INTEGER            MINSIZE, NBMAX\n*     ..\n*     .. External Functions ..\n      LOGICAL            LSAME\n      INTEGER            ILAENV\n      EXTERNAL           LSAME, ILAENV\n*     ..\n*     .. External Subroutines ..\n      EXTERNAL           ZSYTRI2X\n*     ..\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, a = NumRu::Lapack.zsytri2( uplo, a, ipiv, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rblapack_uplo = argv[0];
  rblapack_a = argv[1];
  rblapack_ipiv = argv[2];
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
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_a) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of a must be the same as shape 0 of ipiv");
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_DCOMPLEX)
    rblapack_a = na_change_type(rblapack_a, NA_DCOMPLEX);
  a = NA_PTR_TYPE(rblapack_a, doublecomplex*);
  c__1 = 1;
  c__m1 = -1;
  uplo = StringValueCStr(rblapack_uplo)[0];
  nb = ilaenv_(&c__1, "ZSYTRF", &uplo, &n, &c__m1, &c__m1, &c__m1);
  lwork = (n+nb+1)*(nb+3);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rblapack_a_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rblapack_a_out__, doublecomplex*);
  MEMCPY(a_out__, a, doublecomplex, NA_TOTAL(rblapack_a));
  rblapack_a = rblapack_a_out__;
  a = a_out__;
  work = ALLOC_N(doublecomplex, (lwork));

  zsytri2_(&uplo, &n, a, &lda, ipiv, work, &lwork, &info);

  free(work);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(2, rblapack_info, rblapack_a);
}

void
init_lapack_zsytri2(VALUE mLapack){
  rb_define_module_function(mLapack, "zsytri2", rblapack_zsytri2, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
