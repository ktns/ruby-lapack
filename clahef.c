#include "rb_lapack.h"

static VALUE
rb_clahef(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_nb;
  integer nb; 
  VALUE rb_a;
  complex *a; 
  VALUE rb_kb;
  integer kb; 
  VALUE rb_ipiv;
  integer *ipiv; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  complex *a_out__;
  complex *w;

  integer lda;
  integer n;
  integer ldw;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  kb, ipiv, info, a = NumRu::Lapack.clahef( uplo, nb, a)\n    or\n  NumRu::Lapack.clahef  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE CLAHEF( UPLO, N, NB, KB, A, LDA, IPIV, W, LDW, INFO )\n\n*  Purpose\n*  =======\n*\n*  CLAHEF computes a partial factorization of a complex Hermitian\n*  matrix A using the Bunch-Kaufman diagonal pivoting method. The\n*  partial factorization has the form:\n*\n*  A  =  ( I  U12 ) ( A11  0  ) (  I    0   )  if UPLO = 'U', or:\n*        ( 0  U22 ) (  0   D  ) ( U12' U22' )\n*\n*  A  =  ( L11  0 ) (  D   0  ) ( L11' L21' )  if UPLO = 'L'\n*        ( L21  I ) (  0  A22 ) (  0    I   )\n*\n*  where the order of D is at most NB. The actual order is returned in\n*  the argument KB, and is either NB or NB-1, or N if N <= NB.\n*  Note that U' denotes the conjugate transpose of U.\n*\n*  CLAHEF is an auxiliary routine called by CHETRF. It uses blocked code\n*  (calling Level 3 BLAS) to update the submatrix A11 (if UPLO = 'U') or\n*  A22 (if UPLO = 'L').\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO    (input) CHARACTER*1\n*          Specifies whether the upper or lower triangular part of the\n*          Hermitian matrix A is stored:\n*          = 'U':  Upper triangular\n*          = 'L':  Lower triangular\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  NB      (input) INTEGER\n*          The maximum number of columns of the matrix A that should be\n*          factored.  NB should be at least 2 to allow for 2-by-2 pivot\n*          blocks.\n*\n*  KB      (output) INTEGER\n*          The number of columns of A that were actually factored.\n*          KB is either NB-1 or NB, or N if N <= NB.\n*\n*  A       (input/output) COMPLEX array, dimension (LDA,N)\n*          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading\n*          n-by-n upper triangular part of A contains the upper\n*          triangular part of the matrix A, and the strictly lower\n*          triangular part of A is not referenced.  If UPLO = 'L', the\n*          leading n-by-n lower triangular part of A contains the lower\n*          triangular part of the matrix A, and the strictly upper\n*          triangular part of A is not referenced.\n*          On exit, A contains details of the partial factorization.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,N).\n*\n*  IPIV    (output) INTEGER array, dimension (N)\n*          Details of the interchanges and the block structure of D.\n*          If UPLO = 'U', only the last KB elements of IPIV are set;\n*          if UPLO = 'L', only the first KB elements are set.\n*\n*          If IPIV(k) > 0, then rows and columns k and IPIV(k) were\n*          interchanged and D(k,k) is a 1-by-1 diagonal block.\n*          If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and\n*          columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k)\n*          is a 2-by-2 diagonal block.  If UPLO = 'L' and IPIV(k) =\n*          IPIV(k+1) < 0, then rows and columns k+1 and -IPIV(k) were\n*          interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block.\n*\n*  W       (workspace) COMPLEX array, dimension (LDW,NB)\n*\n*  LDW     (input) INTEGER\n*          The leading dimension of the array W.  LDW >= max(1,N).\n*\n*  INFO    (output) INTEGER\n*          = 0: successful exit\n*          > 0: if INFO = k, D(k,k) is exactly zero.  The factorization\n*               has been completed, but the block diagonal matrix D is\n*               exactly singular.\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_uplo = argv[0];
  rb_nb = argv[1];
  rb_a = argv[2];

  uplo = StringValueCStr(rb_uplo)[0];
  nb = NUM2INT(rb_nb);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  lda = NA_SHAPE0(rb_a);
  n = NA_SHAPE1(rb_a);
  if (NA_TYPE(rb_a) != NA_SCOMPLEX)
    rb_a = na_change_type(rb_a, NA_SCOMPLEX);
  a = NA_PTR_TYPE(rb_a, complex*);
  {
    int shape[1];
    shape[0] = DIM_LEN(n);
    rb_ipiv = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  ipiv = NA_PTR_TYPE(rb_ipiv, integer*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rb_a_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rb_a_out__, complex*);
  MEMCPY(a_out__, a, complex, NA_TOTAL(rb_a));
  rb_a = rb_a_out__;
  a = a_out__;
  ldw = MAX(1,n);
  w = ALLOC_N(complex, (ldw)*(MAX(n,nb)));

  clahef_(&uplo, &n, &nb, &kb, a, &lda, ipiv, w, &ldw, &info);

  free(w);
  rb_kb = INT2NUM(kb);
  rb_info = INT2NUM(info);
  return rb_ary_new3(4, rb_kb, rb_ipiv, rb_info, rb_a);
}

void
init_lapack_clahef(VALUE mLapack){
  rb_define_module_function(mLapack, "clahef", rb_clahef, -1);
}
