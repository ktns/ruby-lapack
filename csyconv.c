#include "rb_lapack.h"

extern VOID csyconv_(char *uplo, char *way, integer *n, complex *a, integer *lda, integer *ipiv, complex *work, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_csyconv(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_uplo;
  char uplo; 
  VALUE rblapack_way;
  char way; 
  VALUE rblapack_a;
  complex *a; 
  VALUE rblapack_ipiv;
  integer *ipiv; 
  VALUE rblapack_info;
  integer info; 
  complex *work;

  integer lda;
  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  info = NumRu::Lapack.csyconv( uplo, way, a, ipiv, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE CSYCONV( UPLO, WAY, N, A, LDA, IPIV, WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  CSYCONV convert A given by TRF into L and D and vice-versa.\n*  Get Non-diag elements of D (returned in workspace) and \n*  apply or reverse permutation done in TRF.\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO    (input) CHARACTER*1\n*          Specifies whether the details of the factorization are stored\n*          as an upper or lower triangular matrix.\n*          = 'U':  Upper triangular, form is A = U*D*U**T;\n*          = 'L':  Lower triangular, form is A = L*D*L**T.\n* \n*  WAY     (input) CHARACTER*1\n*          = 'C': Convert \n*          = 'R': Revert\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  A       (input) COMPLEX array, dimension (LDA,N)\n*          The block diagonal matrix D and the multipliers used to\n*          obtain the factor U or L as computed by CSYTRF.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,N).\n*\n*  IPIV    (input) INTEGER array, dimension (N)\n*          Details of the interchanges and the block structure of D\n*          as determined by CSYTRF.\n*\n* WORK     (workspace) COMPLEX array, dimension (N)\n*\n* LWORK    (input) INTEGER\n*          The length of WORK.  LWORK >=1. \n*          LWORK = N\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  info = NumRu::Lapack.csyconv( uplo, way, a, ipiv, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rblapack_uplo = argv[0];
  rblapack_way = argv[1];
  rblapack_a = argv[2];
  rblapack_ipiv = argv[3];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_ipiv))
    rb_raise(rb_eArgError, "ipiv (4th argument) must be NArray");
  if (NA_RANK(rblapack_ipiv) != 1)
    rb_raise(rb_eArgError, "rank of ipiv (4th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_ipiv);
  if (NA_TYPE(rblapack_ipiv) != NA_LINT)
    rblapack_ipiv = na_change_type(rblapack_ipiv, NA_LINT);
  ipiv = NA_PTR_TYPE(rblapack_ipiv, integer*);
  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_a) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of a must be the same as shape 0 of ipiv");
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_SCOMPLEX)
    rblapack_a = na_change_type(rblapack_a, NA_SCOMPLEX);
  a = NA_PTR_TYPE(rblapack_a, complex*);
  way = StringValueCStr(rblapack_way)[0];
  uplo = StringValueCStr(rblapack_uplo)[0];
  work = ALLOC_N(complex, (MAX(1,n)));

  csyconv_(&uplo, &way, &n, a, &lda, ipiv, work, &info);

  free(work);
  rblapack_info = INT2NUM(info);
  return rblapack_info;
}

void
init_lapack_csyconv(VALUE mLapack){
  rb_define_module_function(mLapack, "csyconv", rblapack_csyconv, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
