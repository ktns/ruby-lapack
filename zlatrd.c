#include "rb_lapack.h"

extern VOID zlatrd_(char *uplo, integer *n, integer *nb, doublecomplex *a, integer *lda, doublereal *e, doublecomplex *tau, doublecomplex *w, integer *ldw);

static VALUE sHelp, sUsage;

static VALUE
rblapack_zlatrd(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_uplo;
  char uplo; 
  VALUE rblapack_nb;
  integer nb; 
  VALUE rblapack_a;
  doublecomplex *a; 
  VALUE rblapack_e;
  doublereal *e; 
  VALUE rblapack_tau;
  doublecomplex *tau; 
  VALUE rblapack_w;
  doublecomplex *w; 
  VALUE rblapack_a_out__;
  doublecomplex *a_out__;

  integer lda;
  integer n;
  integer ldw;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  e, tau, w, a = NumRu::Lapack.zlatrd( uplo, nb, a, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZLATRD( UPLO, N, NB, A, LDA, E, TAU, W, LDW )\n\n*  Purpose\n*  =======\n*\n*  ZLATRD reduces NB rows and columns of a complex Hermitian matrix A to\n*  Hermitian tridiagonal form by a unitary similarity\n*  transformation Q' * A * Q, and returns the matrices V and W which are\n*  needed to apply the transformation to the unreduced part of A.\n*\n*  If UPLO = 'U', ZLATRD reduces the last NB rows and columns of a\n*  matrix, of which the upper triangle is supplied;\n*  if UPLO = 'L', ZLATRD reduces the first NB rows and columns of a\n*  matrix, of which the lower triangle is supplied.\n*\n*  This is an auxiliary routine called by ZHETRD.\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO    (input) CHARACTER*1\n*          Specifies whether the upper or lower triangular part of the\n*          Hermitian matrix A is stored:\n*          = 'U': Upper triangular\n*          = 'L': Lower triangular\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.\n*\n*  NB      (input) INTEGER\n*          The number of rows and columns to be reduced.\n*\n*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)\n*          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading\n*          n-by-n upper triangular part of A contains the upper\n*          triangular part of the matrix A, and the strictly lower\n*          triangular part of A is not referenced.  If UPLO = 'L', the\n*          leading n-by-n lower triangular part of A contains the lower\n*          triangular part of the matrix A, and the strictly upper\n*          triangular part of A is not referenced.\n*          On exit:\n*          if UPLO = 'U', the last NB columns have been reduced to\n*            tridiagonal form, with the diagonal elements overwriting\n*            the diagonal elements of A; the elements above the diagonal\n*            with the array TAU, represent the unitary matrix Q as a\n*            product of elementary reflectors;\n*          if UPLO = 'L', the first NB columns have been reduced to\n*            tridiagonal form, with the diagonal elements overwriting\n*            the diagonal elements of A; the elements below the diagonal\n*            with the array TAU, represent the  unitary matrix Q as a\n*            product of elementary reflectors.\n*          See Further Details.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,N).\n*\n*  E       (output) DOUBLE PRECISION array, dimension (N-1)\n*          If UPLO = 'U', E(n-nb:n-1) contains the superdiagonal\n*          elements of the last NB columns of the reduced matrix;\n*          if UPLO = 'L', E(1:nb) contains the subdiagonal elements of\n*          the first NB columns of the reduced matrix.\n*\n*  TAU     (output) COMPLEX*16 array, dimension (N-1)\n*          The scalar factors of the elementary reflectors, stored in\n*          TAU(n-nb:n-1) if UPLO = 'U', and in TAU(1:nb) if UPLO = 'L'.\n*          See Further Details.\n*\n*  W       (output) COMPLEX*16 array, dimension (LDW,NB)\n*          The n-by-nb matrix W required to update the unreduced part\n*          of A.\n*\n*  LDW     (input) INTEGER\n*          The leading dimension of the array W. LDW >= max(1,N).\n*\n\n*  Further Details\n*  ===============\n*\n*  If UPLO = 'U', the matrix Q is represented as a product of elementary\n*  reflectors\n*\n*     Q = H(n) H(n-1) . . . H(n-nb+1).\n*\n*  Each H(i) has the form\n*\n*     H(i) = I - tau * v * v'\n*\n*  where tau is a complex scalar, and v is a complex vector with\n*  v(i:n) = 0 and v(i-1) = 1; v(1:i-1) is stored on exit in A(1:i-1,i),\n*  and tau in TAU(i-1).\n*\n*  If UPLO = 'L', the matrix Q is represented as a product of elementary\n*  reflectors\n*\n*     Q = H(1) H(2) . . . H(nb).\n*\n*  Each H(i) has the form\n*\n*     H(i) = I - tau * v * v'\n*\n*  where tau is a complex scalar, and v is a complex vector with\n*  v(1:i) = 0 and v(i+1) = 1; v(i+1:n) is stored on exit in A(i+1:n,i),\n*  and tau in TAU(i).\n*\n*  The elements of the vectors v together form the n-by-nb matrix V\n*  which is needed, with W, to apply the transformation to the unreduced\n*  part of the matrix, using a Hermitian rank-2k update of the form:\n*  A := A - V*W' - W*V'.\n*\n*  The contents of A on exit are illustrated by the following examples\n*  with n = 5 and nb = 2:\n*\n*  if UPLO = 'U':                       if UPLO = 'L':\n*\n*    (  a   a   a   v4  v5 )              (  d                  )\n*    (      a   a   v4  v5 )              (  1   d              )\n*    (          a   1   v5 )              (  v1  1   a          )\n*    (              d   1  )              (  v1  v2  a   a      )\n*    (                  d  )              (  v1  v2  a   a   a  )\n*\n*  where d denotes a diagonal element of the reduced matrix, a denotes\n*  an element of the original matrix that is unchanged, and vi denotes\n*  an element of the vector defining H(i).\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  e, tau, w, a = NumRu::Lapack.zlatrd( uplo, nb, a, [:usage => usage, :help => help])\n");
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
  if (NA_TYPE(rblapack_a) != NA_DCOMPLEX)
    rblapack_a = na_change_type(rblapack_a, NA_DCOMPLEX);
  a = NA_PTR_TYPE(rblapack_a, doublecomplex*);
  uplo = StringValueCStr(rblapack_uplo)[0];
  nb = NUM2INT(rblapack_nb);
  ldw = MAX(1,n);
  {
    int shape[1];
    shape[0] = n-1;
    rblapack_e = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  e = NA_PTR_TYPE(rblapack_e, doublereal*);
  {
    int shape[1];
    shape[0] = n-1;
    rblapack_tau = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  tau = NA_PTR_TYPE(rblapack_tau, doublecomplex*);
  {
    int shape[2];
    shape[0] = ldw;
    shape[1] = MAX(n,nb);
    rblapack_w = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  w = NA_PTR_TYPE(rblapack_w, doublecomplex*);
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

  zlatrd_(&uplo, &n, &nb, a, &lda, e, tau, w, &ldw);

  return rb_ary_new3(4, rblapack_e, rblapack_tau, rblapack_w, rblapack_a);
}

void
init_lapack_zlatrd(VALUE mLapack){
  rb_define_module_function(mLapack, "zlatrd", rblapack_zlatrd, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
