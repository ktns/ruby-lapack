#include "rb_lapack.h"

extern VOID dgebrd_(integer *m, integer *n, doublereal *a, integer *lda, doublereal *d, doublereal *e, doublereal *tauq, doublereal *taup, doublereal *work, integer *lwork, integer *info);

static VALUE
rb_dgebrd(int argc, VALUE *argv, VALUE self){
  VALUE rb_m;
  integer m; 
  VALUE rb_a;
  doublereal *a; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_d;
  doublereal *d; 
  VALUE rb_e;
  doublereal *e; 
  VALUE rb_tauq;
  doublereal *tauq; 
  VALUE rb_taup;
  doublereal *taup; 
  VALUE rb_work;
  doublereal *work; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  doublereal *a_out__;

  integer lda;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  d, e, tauq, taup, work, info, a = NumRu::Lapack.dgebrd( m, a, lwork)\n    or\n  NumRu::Lapack.dgebrd  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DGEBRD( M, N, A, LDA, D, E, TAUQ, TAUP, WORK, LWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  DGEBRD reduces a general real M-by-N matrix A to upper or lower\n*  bidiagonal form B by an orthogonal transformation: Q**T * A * P = B.\n*\n*  If m >= n, B is upper bidiagonal; if m < n, B is lower bidiagonal.\n*\n\n*  Arguments\n*  =========\n*\n*  M       (input) INTEGER\n*          The number of rows in the matrix A.  M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns in the matrix A.  N >= 0.\n*\n*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)\n*          On entry, the M-by-N general matrix to be reduced.\n*          On exit,\n*          if m >= n, the diagonal and the first superdiagonal are\n*            overwritten with the upper bidiagonal matrix B; the\n*            elements below the diagonal, with the array TAUQ, represent\n*            the orthogonal matrix Q as a product of elementary\n*            reflectors, and the elements above the first superdiagonal,\n*            with the array TAUP, represent the orthogonal matrix P as\n*            a product of elementary reflectors;\n*          if m < n, the diagonal and the first subdiagonal are\n*            overwritten with the lower bidiagonal matrix B; the\n*            elements below the first subdiagonal, with the array TAUQ,\n*            represent the orthogonal matrix Q as a product of\n*            elementary reflectors, and the elements above the diagonal,\n*            with the array TAUP, represent the orthogonal matrix P as\n*            a product of elementary reflectors.\n*          See Further Details.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,M).\n*\n*  D       (output) DOUBLE PRECISION array, dimension (min(M,N))\n*          The diagonal elements of the bidiagonal matrix B:\n*          D(i) = A(i,i).\n*\n*  E       (output) DOUBLE PRECISION array, dimension (min(M,N)-1)\n*          The off-diagonal elements of the bidiagonal matrix B:\n*          if m >= n, E(i) = A(i,i+1) for i = 1,2,...,n-1;\n*          if m < n, E(i) = A(i+1,i) for i = 1,2,...,m-1.\n*\n*  TAUQ    (output) DOUBLE PRECISION array dimension (min(M,N))\n*          The scalar factors of the elementary reflectors which\n*          represent the orthogonal matrix Q. See Further Details.\n*\n*  TAUP    (output) DOUBLE PRECISION array, dimension (min(M,N))\n*          The scalar factors of the elementary reflectors which\n*          represent the orthogonal matrix P. See Further Details.\n*\n*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The length of the array WORK.  LWORK >= max(1,M,N).\n*          For optimum performance LWORK >= (M+N)*NB, where NB\n*          is the optimal blocksize.\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*\n\n*  Further Details\n*  ===============\n*\n*  The matrices Q and P are represented as products of elementary\n*  reflectors:\n*\n*  If m >= n,\n*\n*     Q = H(1) H(2) . . . H(n)  and  P = G(1) G(2) . . . G(n-1)\n*\n*  Each H(i) and G(i) has the form:\n*\n*     H(i) = I - tauq * v * v'  and G(i) = I - taup * u * u'\n*\n*  where tauq and taup are real scalars, and v and u are real vectors;\n*  v(1:i-1) = 0, v(i) = 1, and v(i+1:m) is stored on exit in A(i+1:m,i);\n*  u(1:i) = 0, u(i+1) = 1, and u(i+2:n) is stored on exit in A(i,i+2:n);\n*  tauq is stored in TAUQ(i) and taup in TAUP(i).\n*\n*  If m < n,\n*\n*     Q = H(1) H(2) . . . H(m-1)  and  P = G(1) G(2) . . . G(m)\n*\n*  Each H(i) and G(i) has the form:\n*\n*     H(i) = I - tauq * v * v'  and G(i) = I - taup * u * u'\n*\n*  where tauq and taup are real scalars, and v and u are real vectors;\n*  v(1:i) = 0, v(i+1) = 1, and v(i+2:m) is stored on exit in A(i+2:m,i);\n*  u(1:i-1) = 0, u(i) = 1, and u(i+1:n) is stored on exit in A(i,i+1:n);\n*  tauq is stored in TAUQ(i) and taup in TAUP(i).\n*\n*  The contents of A on exit are illustrated by the following examples:\n*\n*  m = 6 and n = 5 (m > n):          m = 5 and n = 6 (m < n):\n*\n*    (  d   e   u1  u1  u1 )           (  d   u1  u1  u1  u1  u1 )\n*    (  v1  d   e   u2  u2 )           (  e   d   u2  u2  u2  u2 )\n*    (  v1  v2  d   e   u3 )           (  v1  e   d   u3  u3  u3 )\n*    (  v1  v2  v3  d   e  )           (  v1  v2  e   d   u4  u4 )\n*    (  v1  v2  v3  v4  d  )           (  v1  v2  v3  e   d   u5 )\n*    (  v1  v2  v3  v4  v5 )\n*\n*  where d and e denote diagonal and off-diagonal elements of B, vi\n*  denotes an element of the vector defining H(i), and ui an element of\n*  the vector defining G(i).\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_m = argv[0];
  rb_a = argv[1];
  rb_lwork = argv[2];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_DFLOAT)
    rb_a = na_change_type(rb_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rb_a, doublereal*);
  m = NUM2INT(rb_m);
  lwork = NUM2INT(rb_lwork);
  {
    int shape[1];
    shape[0] = MIN(m,n);
    rb_d = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  d = NA_PTR_TYPE(rb_d, doublereal*);
  {
    int shape[1];
    shape[0] = MIN(m,n)-1;
    rb_e = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  e = NA_PTR_TYPE(rb_e, doublereal*);
  {
    int shape[1];
    shape[0] = MIN(m,n);
    rb_tauq = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  tauq = NA_PTR_TYPE(rb_tauq, doublereal*);
  {
    int shape[1];
    shape[0] = MIN(m,n);
    rb_taup = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  taup = NA_PTR_TYPE(rb_taup, doublereal*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rb_work = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rb_work, doublereal*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rb_a_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rb_a_out__, doublereal*);
  MEMCPY(a_out__, a, doublereal, NA_TOTAL(rb_a));
  rb_a = rb_a_out__;
  a = a_out__;

  dgebrd_(&m, &n, a, &lda, d, e, tauq, taup, work, &lwork, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(7, rb_d, rb_e, rb_tauq, rb_taup, rb_work, rb_info, rb_a);
}

void
init_lapack_dgebrd(VALUE mLapack){
  rb_define_module_function(mLapack, "dgebrd", rb_dgebrd, -1);
}
