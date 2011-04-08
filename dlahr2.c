#include "rb_lapack.h"

extern VOID dlahr2_(integer *n, integer *k, integer *nb, doublereal *a, integer *lda, doublereal *tau, doublereal *t, integer *ldt, doublereal *y, integer *ldy);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dlahr2(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_n;
  integer n; 
  VALUE rblapack_k;
  integer k; 
  VALUE rblapack_nb;
  integer nb; 
  VALUE rblapack_a;
  doublereal *a; 
  VALUE rblapack_tau;
  doublereal *tau; 
  VALUE rblapack_t;
  doublereal *t; 
  VALUE rblapack_y;
  doublereal *y; 
  VALUE rblapack_a_out__;
  doublereal *a_out__;

  integer lda;
  integer ldt;
  integer ldy;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  tau, t, y, a = NumRu::Lapack.dlahr2( n, k, nb, a, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLAHR2( N, K, NB, A, LDA, TAU, T, LDT, Y, LDY )\n\n*  Purpose\n*  =======\n*\n*  DLAHR2 reduces the first NB columns of A real general n-BY-(n-k+1)\n*  matrix A so that elements below the k-th subdiagonal are zero. The\n*  reduction is performed by an orthogonal similarity transformation\n*  Q' * A * Q. The routine returns the matrices V and T which determine\n*  Q as a block reflector I - V*T*V', and also the matrix Y = A * V * T.\n*\n*  This is an auxiliary routine called by DGEHRD.\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.\n*\n*  K       (input) INTEGER\n*          The offset for the reduction. Elements below the k-th\n*          subdiagonal in the first NB columns are reduced to zero.\n*          K < N.\n*\n*  NB      (input) INTEGER\n*          The number of columns to be reduced.\n*\n*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N-K+1)\n*          On entry, the n-by-(n-k+1) general matrix A.\n*          On exit, the elements on and above the k-th subdiagonal in\n*          the first NB columns are overwritten with the corresponding\n*          elements of the reduced matrix; the elements below the k-th\n*          subdiagonal, with the array TAU, represent the matrix Q as a\n*          product of elementary reflectors. The other columns of A are\n*          unchanged. See Further Details.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,N).\n*\n*  TAU     (output) DOUBLE PRECISION array, dimension (NB)\n*          The scalar factors of the elementary reflectors. See Further\n*          Details.\n*\n*  T       (output) DOUBLE PRECISION array, dimension (LDT,NB)\n*          The upper triangular matrix T.\n*\n*  LDT     (input) INTEGER\n*          The leading dimension of the array T.  LDT >= NB.\n*\n*  Y       (output) DOUBLE PRECISION array, dimension (LDY,NB)\n*          The n-by-nb matrix Y.\n*\n*  LDY     (input) INTEGER\n*          The leading dimension of the array Y. LDY >= N.\n*\n\n*  Further Details\n*  ===============\n*\n*  The matrix Q is represented as a product of nb elementary reflectors\n*\n*     Q = H(1) H(2) . . . H(nb).\n*\n*  Each H(i) has the form\n*\n*     H(i) = I - tau * v * v'\n*\n*  where tau is a real scalar, and v is a real vector with\n*  v(1:i+k-1) = 0, v(i+k) = 1; v(i+k+1:n) is stored on exit in\n*  A(i+k+1:n,i), and tau in TAU(i).\n*\n*  The elements of the vectors v together form the (n-k+1)-by-nb matrix\n*  V which is needed, with T and Y, to apply the transformation to the\n*  unreduced part of the matrix, using an update of the form:\n*  A := (I - V*T*V') * (A - Y*V').\n*\n*  The contents of A on exit are illustrated by the following example\n*  with n = 7, k = 3 and nb = 2:\n*\n*     ( a   a   a   a   a )\n*     ( a   a   a   a   a )\n*     ( a   a   a   a   a )\n*     ( h   h   a   a   a )\n*     ( v1  h   a   a   a )\n*     ( v1  v2  a   a   a )\n*     ( v1  v2  a   a   a )\n*\n*  where a denotes an element of the original matrix A, h denotes a\n*  modified element of the upper Hessenberg matrix H, and vi denotes an\n*  element of the vector defining H(i).\n*\n*  This subroutine is a slight modification of LAPACK-3.0's DLAHRD\n*  incorporating improvements proposed by Quintana-Orti and Van de\n*  Gejin. Note that the entries of A(1:K,2:NB) differ from those\n*  returned by the original LAPACK-3.0's DLAHRD routine. (This\n*  subroutine is not backward compatible with LAPACK-3.0's DLAHRD.)\n*\n*  References\n*  ==========\n*\n*  Gregorio Quintana-Orti and Robert van de Geijn, \"Improving the\n*  performance of reduction to Hessenberg form,\" ACM Transactions on\n*  Mathematical Software, 32(2):180-194, June 2006.\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  tau, t, y, a = NumRu::Lapack.dlahr2( n, k, nb, a, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rblapack_n = argv[0];
  rblapack_k = argv[1];
  rblapack_nb = argv[2];
  rblapack_a = argv[3];
  if (rb_options != Qnil) {
  }

  k = NUM2INT(rblapack_k);
  nb = NUM2INT(rblapack_nb);
  n = NUM2INT(rblapack_n);
  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (4th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (4th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_a) != (n-k+1))
    rb_raise(rb_eRuntimeError, "shape 1 of a must be %d", n-k+1);
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_DFLOAT)
    rblapack_a = na_change_type(rblapack_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rblapack_a, doublereal*);
  ldt = nb;
  ldy = n;
  {
    int shape[1];
    shape[0] = MAX(1,nb);
    rblapack_tau = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  tau = NA_PTR_TYPE(rblapack_tau, doublereal*);
  {
    int shape[2];
    shape[0] = ldt;
    shape[1] = MAX(1,nb);
    rblapack_t = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  t = NA_PTR_TYPE(rblapack_t, doublereal*);
  {
    int shape[2];
    shape[0] = ldy;
    shape[1] = MAX(1,nb);
    rblapack_y = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  y = NA_PTR_TYPE(rblapack_y, doublereal*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n-k+1;
    rblapack_a_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rblapack_a_out__, doublereal*);
  MEMCPY(a_out__, a, doublereal, NA_TOTAL(rblapack_a));
  rblapack_a = rblapack_a_out__;
  a = a_out__;

  dlahr2_(&n, &k, &nb, a, &lda, tau, t, &ldt, y, &ldy);

  return rb_ary_new3(4, rblapack_tau, rblapack_t, rblapack_y, rblapack_a);
}

void
init_lapack_dlahr2(VALUE mLapack){
  rb_define_module_function(mLapack, "dlahr2", rblapack_dlahr2, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
