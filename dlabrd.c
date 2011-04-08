#include "rb_lapack.h"

extern VOID dlabrd_(integer *m, integer *n, integer *nb, doublereal *a, integer *lda, doublereal *d, doublereal *e, doublereal *tauq, doublereal *taup, doublereal *x, integer *ldx, doublereal *y, integer *ldy);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dlabrd(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_m;
  integer m; 
  VALUE rblapack_nb;
  integer nb; 
  VALUE rblapack_a;
  doublereal *a; 
  VALUE rblapack_d;
  doublereal *d; 
  VALUE rblapack_e;
  doublereal *e; 
  VALUE rblapack_tauq;
  doublereal *tauq; 
  VALUE rblapack_taup;
  doublereal *taup; 
  VALUE rblapack_x;
  doublereal *x; 
  VALUE rblapack_y;
  doublereal *y; 
  VALUE rblapack_a_out__;
  doublereal *a_out__;

  integer lda;
  integer n;
  integer ldx;
  integer ldy;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  d, e, tauq, taup, x, y, a = NumRu::Lapack.dlabrd( m, nb, a, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLABRD( M, N, NB, A, LDA, D, E, TAUQ, TAUP, X, LDX, Y, LDY )\n\n*  Purpose\n*  =======\n*\n*  DLABRD reduces the first NB rows and columns of a real general\n*  m by n matrix A to upper or lower bidiagonal form by an orthogonal\n*  transformation Q' * A * P, and returns the matrices X and Y which\n*  are needed to apply the transformation to the unreduced part of A.\n*\n*  If m >= n, A is reduced to upper bidiagonal form; if m < n, to lower\n*  bidiagonal form.\n*\n*  This is an auxiliary routine called by DGEBRD\n*\n\n*  Arguments\n*  =========\n*\n*  M       (input) INTEGER\n*          The number of rows in the matrix A.\n*\n*  N       (input) INTEGER\n*          The number of columns in the matrix A.\n*\n*  NB      (input) INTEGER\n*          The number of leading rows and columns of A to be reduced.\n*\n*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)\n*          On entry, the m by n general matrix to be reduced.\n*          On exit, the first NB rows and columns of the matrix are\n*          overwritten; the rest of the array is unchanged.\n*          If m >= n, elements on and below the diagonal in the first NB\n*            columns, with the array TAUQ, represent the orthogonal\n*            matrix Q as a product of elementary reflectors; and\n*            elements above the diagonal in the first NB rows, with the\n*            array TAUP, represent the orthogonal matrix P as a product\n*            of elementary reflectors.\n*          If m < n, elements below the diagonal in the first NB\n*            columns, with the array TAUQ, represent the orthogonal\n*            matrix Q as a product of elementary reflectors, and\n*            elements on and above the diagonal in the first NB rows,\n*            with the array TAUP, represent the orthogonal matrix P as\n*            a product of elementary reflectors.\n*          See Further Details.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,M).\n*\n*  D       (output) DOUBLE PRECISION array, dimension (NB)\n*          The diagonal elements of the first NB rows and columns of\n*          the reduced matrix.  D(i) = A(i,i).\n*\n*  E       (output) DOUBLE PRECISION array, dimension (NB)\n*          The off-diagonal elements of the first NB rows and columns of\n*          the reduced matrix.\n*\n*  TAUQ    (output) DOUBLE PRECISION array dimension (NB)\n*          The scalar factors of the elementary reflectors which\n*          represent the orthogonal matrix Q. See Further Details.\n*\n*  TAUP    (output) DOUBLE PRECISION array, dimension (NB)\n*          The scalar factors of the elementary reflectors which\n*          represent the orthogonal matrix P. See Further Details.\n*\n*  X       (output) DOUBLE PRECISION array, dimension (LDX,NB)\n*          The m-by-nb matrix X required to update the unreduced part\n*          of A.\n*\n*  LDX     (input) INTEGER\n*          The leading dimension of the array X. LDX >= M.\n*\n*  Y       (output) DOUBLE PRECISION array, dimension (LDY,NB)\n*          The n-by-nb matrix Y required to update the unreduced part\n*          of A.\n*\n*  LDY     (input) INTEGER\n*          The leading dimension of the array Y. LDY >= N.\n*\n\n*  Further Details\n*  ===============\n*\n*  The matrices Q and P are represented as products of elementary\n*  reflectors:\n*\n*     Q = H(1) H(2) . . . H(nb)  and  P = G(1) G(2) . . . G(nb)\n*\n*  Each H(i) and G(i) has the form:\n*\n*     H(i) = I - tauq * v * v'  and G(i) = I - taup * u * u'\n*\n*  where tauq and taup are real scalars, and v and u are real vectors.\n*\n*  If m >= n, v(1:i-1) = 0, v(i) = 1, and v(i:m) is stored on exit in\n*  A(i:m,i); u(1:i) = 0, u(i+1) = 1, and u(i+1:n) is stored on exit in\n*  A(i,i+1:n); tauq is stored in TAUQ(i) and taup in TAUP(i).\n*\n*  If m < n, v(1:i) = 0, v(i+1) = 1, and v(i+1:m) is stored on exit in\n*  A(i+2:m,i); u(1:i-1) = 0, u(i) = 1, and u(i:n) is stored on exit in\n*  A(i,i+1:n); tauq is stored in TAUQ(i) and taup in TAUP(i).\n*\n*  The elements of the vectors v and u together form the m-by-nb matrix\n*  V and the nb-by-n matrix U' which are needed, with X and Y, to apply\n*  the transformation to the unreduced part of the matrix, using a block\n*  update of the form:  A := A - V*Y' - X*U'.\n*\n*  The contents of A on exit are illustrated by the following examples\n*  with nb = 2:\n*\n*  m = 6 and n = 5 (m > n):          m = 5 and n = 6 (m < n):\n*\n*    (  1   1   u1  u1  u1 )           (  1   u1  u1  u1  u1  u1 )\n*    (  v1  1   1   u2  u2 )           (  1   1   u2  u2  u2  u2 )\n*    (  v1  v2  a   a   a  )           (  v1  1   a   a   a   a  )\n*    (  v1  v2  a   a   a  )           (  v1  v2  a   a   a   a  )\n*    (  v1  v2  a   a   a  )           (  v1  v2  a   a   a   a  )\n*    (  v1  v2  a   a   a  )\n*\n*  where a denotes an element of the original matrix which is unchanged,\n*  vi denotes an element of the vector defining H(i), and ui an element\n*  of the vector defining G(i).\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  d, e, tauq, taup, x, y, a = NumRu::Lapack.dlabrd( m, nb, a, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rblapack_m = argv[0];
  rblapack_nb = argv[1];
  rblapack_a = argv[2];
  if (rb_options != Qnil) {
  }

  nb = NUM2INT(rblapack_nb);
  m = NUM2INT(rblapack_m);
  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_a);
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_DFLOAT)
    rblapack_a = na_change_type(rblapack_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rblapack_a, doublereal*);
  ldy = n;
  ldx = m;
  {
    int shape[1];
    shape[0] = MAX(1,nb);
    rblapack_d = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  d = NA_PTR_TYPE(rblapack_d, doublereal*);
  {
    int shape[1];
    shape[0] = MAX(1,nb);
    rblapack_e = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  e = NA_PTR_TYPE(rblapack_e, doublereal*);
  {
    int shape[1];
    shape[0] = MAX(1,nb);
    rblapack_tauq = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  tauq = NA_PTR_TYPE(rblapack_tauq, doublereal*);
  {
    int shape[1];
    shape[0] = MAX(1,nb);
    rblapack_taup = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  taup = NA_PTR_TYPE(rblapack_taup, doublereal*);
  {
    int shape[2];
    shape[0] = ldx;
    shape[1] = MAX(1,nb);
    rblapack_x = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  x = NA_PTR_TYPE(rblapack_x, doublereal*);
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
    shape[1] = n;
    rblapack_a_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rblapack_a_out__, doublereal*);
  MEMCPY(a_out__, a, doublereal, NA_TOTAL(rblapack_a));
  rblapack_a = rblapack_a_out__;
  a = a_out__;

  dlabrd_(&m, &n, &nb, a, &lda, d, e, tauq, taup, x, &ldx, y, &ldy);

  return rb_ary_new3(7, rblapack_d, rblapack_e, rblapack_tauq, rblapack_taup, rblapack_x, rblapack_y, rblapack_a);
}

void
init_lapack_dlabrd(VALUE mLapack){
  rb_define_module_function(mLapack, "dlabrd", rblapack_dlabrd, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
