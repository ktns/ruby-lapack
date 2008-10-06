#include "rb_lapack.h"

static VALUE
rb_clatrd(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_nb;
  integer nb; 
  VALUE rb_a;
  complex *a; 
  VALUE rb_e;
  real *e; 
  VALUE rb_tau;
  complex *tau; 
  VALUE rb_w;
  complex *w; 
  VALUE rb_a_out__;
  complex *a_out__;

  integer lda;
  integer n;
  integer ldw;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  e, tau, w, a = NumRu::Lapack.clatrd( uplo, nb, a)\n    or\n  NumRu::Lapack.clatrd  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE CLATRD( UPLO, N, NB, A, LDA, E, TAU, W, LDW )\n\n*  Purpose\n*  =======\n*\n*  CLATRD reduces NB rows and columns of a complex Hermitian matrix A to\n*  Hermitian tridiagonal form by a unitary similarity\n*  transformation Q' * A * Q, and returns the matrices V and W which are\n*  needed to apply the transformation to the unreduced part of A.\n*\n*  If UPLO = 'U', CLATRD reduces the last NB rows and columns of a\n*  matrix, of which the upper triangle is supplied;\n*  if UPLO = 'L', CLATRD reduces the first NB rows and columns of a\n*  matrix, of which the lower triangle is supplied.\n*\n*  This is an auxiliary routine called by CHETRD.\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO    (input) CHARACTER*1\n*          Specifies whether the upper or lower triangular part of the\n*          Hermitian matrix A is stored:\n*          = 'U': Upper triangular\n*          = 'L': Lower triangular\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.\n*\n*  NB      (input) INTEGER\n*          The number of rows and columns to be reduced.\n*\n*  A       (input/output) COMPLEX array, dimension (LDA,N)\n*          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading\n*          n-by-n upper triangular part of A contains the upper\n*          triangular part of the matrix A, and the strictly lower\n*          triangular part of A is not referenced.  If UPLO = 'L', the\n*          leading n-by-n lower triangular part of A contains the lower\n*          triangular part of the matrix A, and the strictly upper\n*          triangular part of A is not referenced.\n*          On exit:\n*          if UPLO = 'U', the last NB columns have been reduced to\n*            tridiagonal form, with the diagonal elements overwriting\n*            the diagonal elements of A; the elements above the diagonal\n*            with the array TAU, represent the unitary matrix Q as a\n*            product of elementary reflectors;\n*          if UPLO = 'L', the first NB columns have been reduced to\n*            tridiagonal form, with the diagonal elements overwriting\n*            the diagonal elements of A; the elements below the diagonal\n*            with the array TAU, represent the  unitary matrix Q as a\n*            product of elementary reflectors.\n*          See Further Details.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,N).\n*\n*  E       (output) REAL array, dimension (N-1)\n*          If UPLO = 'U', E(n-nb:n-1) contains the superdiagonal\n*          elements of the last NB columns of the reduced matrix;\n*          if UPLO = 'L', E(1:nb) contains the subdiagonal elements of\n*          the first NB columns of the reduced matrix.\n*\n*  TAU     (output) COMPLEX array, dimension (N-1)\n*          The scalar factors of the elementary reflectors, stored in\n*          TAU(n-nb:n-1) if UPLO = 'U', and in TAU(1:nb) if UPLO = 'L'.\n*          See Further Details.\n*\n*  W       (output) COMPLEX array, dimension (LDW,NB)\n*          The n-by-nb matrix W required to update the unreduced part\n*          of A.\n*\n*  LDW     (input) INTEGER\n*          The leading dimension of the array W. LDW >= max(1,N).\n*\n\n*  Further Details\n*  ===============\n*\n*  If UPLO = 'U', the matrix Q is represented as a product of elementary\n*  reflectors\n*\n*     Q = H(n) H(n-1) . . . H(n-nb+1).\n*\n*  Each H(i) has the form\n*\n*     H(i) = I - tau * v * v'\n*\n*  where tau is a complex scalar, and v is a complex vector with\n*  v(i:n) = 0 and v(i-1) = 1; v(1:i-1) is stored on exit in A(1:i-1,i),\n*  and tau in TAU(i-1).\n*\n*  If UPLO = 'L', the matrix Q is represented as a product of elementary\n*  reflectors\n*\n*     Q = H(1) H(2) . . . H(nb).\n*\n*  Each H(i) has the form\n*\n*     H(i) = I - tau * v * v'\n*\n*  where tau is a complex scalar, and v is a complex vector with\n*  v(1:i) = 0 and v(i+1) = 1; v(i+1:n) is stored on exit in A(i+1:n,i),\n*  and tau in TAU(i).\n*\n*  The elements of the vectors v together form the n-by-nb matrix V\n*  which is needed, with W, to apply the transformation to the unreduced\n*  part of the matrix, using a Hermitian rank-2k update of the form:\n*  A := A - V*W' - W*V'.\n*\n*  The contents of A on exit are illustrated by the following examples\n*  with n = 5 and nb = 2:\n*\n*  if UPLO = 'U':                       if UPLO = 'L':\n*\n*    (  a   a   a   v4  v5 )              (  d                  )\n*    (      a   a   v4  v5 )              (  1   d              )\n*    (          a   1   v5 )              (  v1  1   a          )\n*    (              d   1  )              (  v1  v2  a   a      )\n*    (                  d  )              (  v1  v2  a   a   a  )\n*\n*  where d denotes a diagonal element of the reduced matrix, a denotes\n*  an element of the original matrix that is unchanged, and vi denotes\n*  an element of the vector defining H(i).\n*\n*  =====================================================================\n*\n\n");
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
    shape[0] = n-1;
    rb_e = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  e = NA_PTR_TYPE(rb_e, real*);
  {
    int shape[1];
    shape[0] = n-1;
    rb_tau = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  tau = NA_PTR_TYPE(rb_tau, complex*);
  ldw = MAX(1,n);
  {
    int shape[2];
    shape[0] = ldw;
    shape[1] = MAX(n,nb);
    rb_w = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  w = NA_PTR_TYPE(rb_w, complex*);
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

  clatrd_(&uplo, &n, &nb, a, &lda, e, tau, w, &ldw);

  return rb_ary_new3(4, rb_e, rb_tau, rb_w, rb_a);
}

void
init_lapack_clatrd(VALUE mLapack){
  rb_define_module_function(mLapack, "clatrd", rb_clatrd, -1);
}
