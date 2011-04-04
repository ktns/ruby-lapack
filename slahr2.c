#include "rb_lapack.h"

extern VOID slahr2_(integer *n, integer *k, integer *nb, real *a, integer *lda, real *tau, real *t, integer *ldt, real *y, integer *ldy);

static VALUE
rb_slahr2(int argc, VALUE *argv, VALUE self){
  VALUE rb_n;
  integer n; 
  VALUE rb_k;
  integer k; 
  VALUE rb_nb;
  integer nb; 
  VALUE rb_a;
  real *a; 
  VALUE rb_tau;
  real *tau; 
  VALUE rb_t;
  real *t; 
  VALUE rb_y;
  real *y; 
  VALUE rb_a_out__;
  real *a_out__;

  integer lda;
  integer ldt;
  integer ldy;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  tau, t, y, a = NumRu::Lapack.slahr2( n, k, nb, a)\n    or\n  NumRu::Lapack.slahr2  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLAHR2( N, K, NB, A, LDA, TAU, T, LDT, Y, LDY )\n\n*  Purpose\n*  =======\n*\n*  SLAHR2 reduces the first NB columns of A real general n-BY-(n-k+1)\n*  matrix A so that elements below the k-th subdiagonal are zero. The\n*  reduction is performed by an orthogonal similarity transformation\n*  Q' * A * Q. The routine returns the matrices V and T which determine\n*  Q as a block reflector I - V*T*V', and also the matrix Y = A * V * T.\n*\n*  This is an auxiliary routine called by SGEHRD.\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.\n*\n*  K       (input) INTEGER\n*          The offset for the reduction. Elements below the k-th\n*          subdiagonal in the first NB columns are reduced to zero.\n*          K < N.\n*\n*  NB      (input) INTEGER\n*          The number of columns to be reduced.\n*\n*  A       (input/output) REAL array, dimension (LDA,N-K+1)\n*          On entry, the n-by-(n-k+1) general matrix A.\n*          On exit, the elements on and above the k-th subdiagonal in\n*          the first NB columns are overwritten with the corresponding\n*          elements of the reduced matrix; the elements below the k-th\n*          subdiagonal, with the array TAU, represent the matrix Q as a\n*          product of elementary reflectors. The other columns of A are\n*          unchanged. See Further Details.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,N).\n*\n*  TAU     (output) REAL array, dimension (NB)\n*          The scalar factors of the elementary reflectors. See Further\n*          Details.\n*\n*  T       (output) REAL array, dimension (LDT,NB)\n*          The upper triangular matrix T.\n*\n*  LDT     (input) INTEGER\n*          The leading dimension of the array T.  LDT >= NB.\n*\n*  Y       (output) REAL array, dimension (LDY,NB)\n*          The n-by-nb matrix Y.\n*\n*  LDY     (input) INTEGER\n*          The leading dimension of the array Y. LDY >= N.\n*\n\n*  Further Details\n*  ===============\n*\n*  The matrix Q is represented as a product of nb elementary reflectors\n*\n*     Q = H(1) H(2) . . . H(nb).\n*\n*  Each H(i) has the form\n*\n*     H(i) = I - tau * v * v'\n*\n*  where tau is a real scalar, and v is a real vector with\n*  v(1:i+k-1) = 0, v(i+k) = 1; v(i+k+1:n) is stored on exit in\n*  A(i+k+1:n,i), and tau in TAU(i).\n*\n*  The elements of the vectors v together form the (n-k+1)-by-nb matrix\n*  V which is needed, with T and Y, to apply the transformation to the\n*  unreduced part of the matrix, using an update of the form:\n*  A := (I - V*T*V') * (A - Y*V').\n*\n*  The contents of A on exit are illustrated by the following example\n*  with n = 7, k = 3 and nb = 2:\n*\n*     ( a   a   a   a   a )\n*     ( a   a   a   a   a )\n*     ( a   a   a   a   a )\n*     ( h   h   a   a   a )\n*     ( v1  h   a   a   a )\n*     ( v1  v2  a   a   a )\n*     ( v1  v2  a   a   a )\n*\n*  where a denotes an element of the original matrix A, h denotes a\n*  modified element of the upper Hessenberg matrix H, and vi denotes an\n*  element of the vector defining H(i).\n*\n*  This subroutine is a slight modification of LAPACK-3.0's DLAHRD\n*  incorporating improvements proposed by Quintana-Orti and Van de\n*  Gejin. Note that the entries of A(1:K,2:NB) differ from those\n*  returned by the original LAPACK-3.0's DLAHRD routine. (This\n*  subroutine is not backward compatible with LAPACK-3.0's DLAHRD.)\n*\n*  References\n*  ==========\n*\n*  Gregorio Quintana-Orti and Robert van de Geijn, \"Improving the\n*  performance of reduction to Hessenberg form,\" ACM Transactions on\n*  Mathematical Software, 32(2):180-194, June 2006.\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_n = argv[0];
  rb_k = argv[1];
  rb_nb = argv[2];
  rb_a = argv[3];

  k = NUM2INT(rb_k);
  nb = NUM2INT(rb_nb);
  n = NUM2INT(rb_n);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (4th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (4th argument) must be %d", 2);
  if (NA_SHAPE1(rb_a) != (n-k+1))
    rb_raise(rb_eRuntimeError, "shape 1 of a must be %d", n-k+1);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);
  ldt = nb;
  ldy = n;
  {
    int shape[1];
    shape[0] = MAX(1,nb);
    rb_tau = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  tau = NA_PTR_TYPE(rb_tau, real*);
  {
    int shape[2];
    shape[0] = ldt;
    shape[1] = MAX(1,nb);
    rb_t = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  t = NA_PTR_TYPE(rb_t, real*);
  {
    int shape[2];
    shape[0] = ldy;
    shape[1] = MAX(1,nb);
    rb_y = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  y = NA_PTR_TYPE(rb_y, real*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n-k+1;
    rb_a_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rb_a_out__, real*);
  MEMCPY(a_out__, a, real, NA_TOTAL(rb_a));
  rb_a = rb_a_out__;
  a = a_out__;

  slahr2_(&n, &k, &nb, a, &lda, tau, t, &ldt, y, &ldy);

  return rb_ary_new3(4, rb_tau, rb_t, rb_y, rb_a);
}

void
init_lapack_slahr2(VALUE mLapack){
  rb_define_module_function(mLapack, "slahr2", rb_slahr2, -1);
}
