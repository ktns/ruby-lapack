#include "rb_lapack.h"

static VALUE
rb_clahrd(int argc, VALUE *argv, VALUE self){
  VALUE rb_n;
  integer n; 
  VALUE rb_k;
  integer k; 
  VALUE rb_nb;
  integer nb; 
  VALUE rb_a;
  complex *a; 
  VALUE rb_tau;
  complex *tau; 
  VALUE rb_t;
  complex *t; 
  VALUE rb_y;
  complex *y; 
  VALUE rb_a_out__;
  complex *a_out__;

  integer lda;
  integer ldt;
  integer ldy;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  tau, t, y, a = NumRu::Lapack.clahrd( n, k, nb, a)\n    or\n  NumRu::Lapack.clahrd  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE CLAHRD( N, K, NB, A, LDA, TAU, T, LDT, Y, LDY )\n\n*  Purpose\n*  =======\n*\n*  CLAHRD reduces the first NB columns of a complex general n-by-(n-k+1)\n*  matrix A so that elements below the k-th subdiagonal are zero. The\n*  reduction is performed by a unitary similarity transformation\n*  Q' * A * Q. The routine returns the matrices V and T which determine\n*  Q as a block reflector I - V*T*V', and also the matrix Y = A * V * T.\n*\n*  This is an OBSOLETE auxiliary routine. \n*  This routine will be 'deprecated' in a  future release.\n*  Please use the new routine CLAHR2 instead.\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.\n*\n*  K       (input) INTEGER\n*          The offset for the reduction. Elements below the k-th\n*          subdiagonal in the first NB columns are reduced to zero.\n*\n*  NB      (input) INTEGER\n*          The number of columns to be reduced.\n*\n*  A       (input/output) COMPLEX array, dimension (LDA,N-K+1)\n*          On entry, the n-by-(n-k+1) general matrix A.\n*          On exit, the elements on and above the k-th subdiagonal in\n*          the first NB columns are overwritten with the corresponding\n*          elements of the reduced matrix; the elements below the k-th\n*          subdiagonal, with the array TAU, represent the matrix Q as a\n*          product of elementary reflectors. The other columns of A are\n*          unchanged. See Further Details.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,N).\n*\n*  TAU     (output) COMPLEX array, dimension (NB)\n*          The scalar factors of the elementary reflectors. See Further\n*          Details.\n*\n*  T       (output) COMPLEX array, dimension (LDT,NB)\n*          The upper triangular matrix T.\n*\n*  LDT     (input) INTEGER\n*          The leading dimension of the array T.  LDT >= NB.\n*\n*  Y       (output) COMPLEX array, dimension (LDY,NB)\n*          The n-by-nb matrix Y.\n*\n*  LDY     (input) INTEGER\n*          The leading dimension of the array Y. LDY >= max(1,N).\n*\n\n*  Further Details\n*  ===============\n*\n*  The matrix Q is represented as a product of nb elementary reflectors\n*\n*     Q = H(1) H(2) . . . H(nb).\n*\n*  Each H(i) has the form\n*\n*     H(i) = I - tau * v * v'\n*\n*  where tau is a complex scalar, and v is a complex vector with\n*  v(1:i+k-1) = 0, v(i+k) = 1; v(i+k+1:n) is stored on exit in\n*  A(i+k+1:n,i), and tau in TAU(i).\n*\n*  The elements of the vectors v together form the (n-k+1)-by-nb matrix\n*  V which is needed, with T and Y, to apply the transformation to the\n*  unreduced part of the matrix, using an update of the form:\n*  A := (I - V*T*V') * (A - Y*V').\n*\n*  The contents of A on exit are illustrated by the following example\n*  with n = 7, k = 3 and nb = 2:\n*\n*     ( a   h   a   a   a )\n*     ( a   h   a   a   a )\n*     ( a   h   a   a   a )\n*     ( h   h   a   a   a )\n*     ( v1  h   a   a   a )\n*     ( v1  v2  a   a   a )\n*     ( v1  v2  a   a   a )\n*\n*  where a denotes an element of the original matrix A, h denotes a\n*  modified element of the upper Hessenberg matrix H, and vi denotes an\n*  element of the vector defining H(i).\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_n = argv[0];
  rb_k = argv[1];
  rb_nb = argv[2];
  rb_a = argv[3];

  n = NUM2INT(rb_n);
  k = NUM2INT(rb_k);
  nb = NUM2INT(rb_nb);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (4th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (4th argument) must be %d", 2);
  lda = NA_SHAPE0(rb_a);
  if (NA_SHAPE1(rb_a) != (n-k+1))
    rb_raise(rb_eRuntimeError, "shape 1 of a must be %d", n-k+1);
  if (NA_TYPE(rb_a) != NA_SCOMPLEX)
    rb_a = na_change_type(rb_a, NA_SCOMPLEX);
  a = NA_PTR_TYPE(rb_a, complex*);
  {
    int shape[1];
    shape[0] = MAX(1,nb);
    rb_tau = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  tau = NA_PTR_TYPE(rb_tau, complex*);
  ldt = nb;
  {
    int shape[2];
    shape[0] = ldt;
    shape[1] = MAX(1,nb);
    rb_t = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  t = NA_PTR_TYPE(rb_t, complex*);
  ldy = MAX(1,n);
  {
    int shape[2];
    shape[0] = ldy;
    shape[1] = MAX(1,nb);
    rb_y = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  y = NA_PTR_TYPE(rb_y, complex*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n-k+1;
    rb_a_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rb_a_out__, complex*);
  MEMCPY(a_out__, a, complex, NA_TOTAL(rb_a));
  rb_a = rb_a_out__;
  a = a_out__;

  clahrd_(&n, &k, &nb, a, &lda, tau, t, &ldt, y, &ldy);

  return rb_ary_new3(4, rb_tau, rb_t, rb_y, rb_a);
}

void
init_lapack_clahrd(VALUE mLapack){
  rb_define_module_function(mLapack, "clahrd", rb_clahrd, -1);
}
