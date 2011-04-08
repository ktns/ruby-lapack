#include "rb_lapack.h"

extern VOID slasyf_(char *uplo, integer *n, integer *nb, integer *kb, real *a, integer *lda, integer *ipiv, real *w, integer *ldw, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_slasyf(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_uplo;
  char uplo; 
  VALUE rblapack_nb;
  integer nb; 
  VALUE rblapack_a;
  real *a; 
  VALUE rblapack_kb;
  integer kb; 
  VALUE rblapack_ipiv;
  integer *ipiv; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_a_out__;
  real *a_out__;
  real *w;

  integer lda;
  integer n;
  integer ldw;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  kb, ipiv, info, a = NumRu::Lapack.slasyf( uplo, nb, a, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLASYF( UPLO, N, NB, KB, A, LDA, IPIV, W, LDW, INFO )\n\n*  Purpose\n*  =======\n*\n*  SLASYF computes a partial factorization of a real symmetric matrix A\n*  using the Bunch-Kaufman diagonal pivoting method. The partial\n*  factorization has the form:\n*\n*  A  =  ( I  U12 ) ( A11  0  ) (  I    0   )  if UPLO = 'U', or:\n*        ( 0  U22 ) (  0   D  ) ( U12' U22' )\n*\n*  A  =  ( L11  0 ) (  D   0  ) ( L11' L21' )  if UPLO = 'L'\n*        ( L21  I ) (  0  A22 ) (  0    I   )\n*\n*  where the order of D is at most NB. The actual order is returned in\n*  the argument KB, and is either NB or NB-1, or N if N <= NB.\n*\n*  SLASYF is an auxiliary routine called by SSYTRF. It uses blocked code\n*  (calling Level 3 BLAS) to update the submatrix A11 (if UPLO = 'U') or\n*  A22 (if UPLO = 'L').\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO    (input) CHARACTER*1\n*          Specifies whether the upper or lower triangular part of the\n*          symmetric matrix A is stored:\n*          = 'U':  Upper triangular\n*          = 'L':  Lower triangular\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  NB      (input) INTEGER\n*          The maximum number of columns of the matrix A that should be\n*          factored.  NB should be at least 2 to allow for 2-by-2 pivot\n*          blocks.\n*\n*  KB      (output) INTEGER\n*          The number of columns of A that were actually factored.\n*          KB is either NB-1 or NB, or N if N <= NB.\n*\n*  A       (input/output) REAL array, dimension (LDA,N)\n*          On entry, the symmetric matrix A.  If UPLO = 'U', the leading\n*          n-by-n upper triangular part of A contains the upper\n*          triangular part of the matrix A, and the strictly lower\n*          triangular part of A is not referenced.  If UPLO = 'L', the\n*          leading n-by-n lower triangular part of A contains the lower\n*          triangular part of the matrix A, and the strictly upper\n*          triangular part of A is not referenced.\n*          On exit, A contains details of the partial factorization.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,N).\n*\n*  IPIV    (output) INTEGER array, dimension (N)\n*          Details of the interchanges and the block structure of D.\n*          If UPLO = 'U', only the last KB elements of IPIV are set;\n*          if UPLO = 'L', only the first KB elements are set.\n*\n*          If IPIV(k) > 0, then rows and columns k and IPIV(k) were\n*          interchanged and D(k,k) is a 1-by-1 diagonal block.\n*          If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and\n*          columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k)\n*          is a 2-by-2 diagonal block.  If UPLO = 'L' and IPIV(k) =\n*          IPIV(k+1) < 0, then rows and columns k+1 and -IPIV(k) were\n*          interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block.\n*\n*  W       (workspace) REAL array, dimension (LDW,NB)\n*\n*  LDW     (input) INTEGER\n*          The leading dimension of the array W.  LDW >= max(1,N).\n*\n*  INFO    (output) INTEGER\n*          = 0: successful exit\n*          > 0: if INFO = k, D(k,k) is exactly zero.  The factorization\n*               has been completed, but the block diagonal matrix D is\n*               exactly singular.\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  kb, ipiv, info, a = NumRu::Lapack.slasyf( uplo, nb, a, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rblapack_uplo = argv[0];
  rblapack_nb = argv[1];
  rblapack_a = argv[2];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_a);
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_SFLOAT)
    rblapack_a = na_change_type(rblapack_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rblapack_a, real*);
  uplo = StringValueCStr(rblapack_uplo)[0];
  nb = NUM2INT(rblapack_nb);
  ldw = MAX(1,n);
  {
    int shape[1];
    shape[0] = n;
    rblapack_ipiv = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  ipiv = NA_PTR_TYPE(rblapack_ipiv, integer*);
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
  w = ALLOC_N(real, (ldw)*(MAX(1,nb)));

  slasyf_(&uplo, &n, &nb, &kb, a, &lda, ipiv, w, &ldw, &info);

  free(w);
  rblapack_kb = INT2NUM(kb);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(4, rblapack_kb, rblapack_ipiv, rblapack_info, rblapack_a);
}

void
init_lapack_slasyf(VALUE mLapack){
  rb_define_module_function(mLapack, "slasyf", rblapack_slasyf, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
