#include "rb_lapack.h"

extern VOID zgebrd_(integer *m, integer *n, doublecomplex *a, integer *lda, doublereal *d, doublereal *e, doublecomplex *tauq, doublecomplex *taup, doublecomplex *work, integer *lwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_zgebrd(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_m;
  integer m; 
  VALUE rblapack_a;
  doublecomplex *a; 
  VALUE rblapack_lwork;
  integer lwork; 
  VALUE rblapack_d;
  doublereal *d; 
  VALUE rblapack_e;
  doublereal *e; 
  VALUE rblapack_tauq;
  doublecomplex *tauq; 
  VALUE rblapack_taup;
  doublecomplex *taup; 
  VALUE rblapack_work;
  doublecomplex *work; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_a_out__;
  doublecomplex *a_out__;

  integer lda;
  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  d, e, tauq, taup, work, info, a = NumRu::Lapack.zgebrd( m, a, lwork, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZGEBRD( M, N, A, LDA, D, E, TAUQ, TAUP, WORK, LWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  ZGEBRD reduces a general complex M-by-N matrix A to upper or lower\n*  bidiagonal form B by a unitary transformation: Q**H * A * P = B.\n*\n*  If m >= n, B is upper bidiagonal; if m < n, B is lower bidiagonal.\n*\n\n*  Arguments\n*  =========\n*\n*  M       (input) INTEGER\n*          The number of rows in the matrix A.  M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns in the matrix A.  N >= 0.\n*\n*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)\n*          On entry, the M-by-N general matrix to be reduced.\n*          On exit,\n*          if m >= n, the diagonal and the first superdiagonal are\n*            overwritten with the upper bidiagonal matrix B; the\n*            elements below the diagonal, with the array TAUQ, represent\n*            the unitary matrix Q as a product of elementary\n*            reflectors, and the elements above the first superdiagonal,\n*            with the array TAUP, represent the unitary matrix P as\n*            a product of elementary reflectors;\n*          if m < n, the diagonal and the first subdiagonal are\n*            overwritten with the lower bidiagonal matrix B; the\n*            elements below the first subdiagonal, with the array TAUQ,\n*            represent the unitary matrix Q as a product of\n*            elementary reflectors, and the elements above the diagonal,\n*            with the array TAUP, represent the unitary matrix P as\n*            a product of elementary reflectors.\n*          See Further Details.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,M).\n*\n*  D       (output) DOUBLE PRECISION array, dimension (min(M,N))\n*          The diagonal elements of the bidiagonal matrix B:\n*          D(i) = A(i,i).\n*\n*  E       (output) DOUBLE PRECISION array, dimension (min(M,N)-1)\n*          The off-diagonal elements of the bidiagonal matrix B:\n*          if m >= n, E(i) = A(i,i+1) for i = 1,2,...,n-1;\n*          if m < n, E(i) = A(i+1,i) for i = 1,2,...,m-1.\n*\n*  TAUQ    (output) COMPLEX*16 array dimension (min(M,N))\n*          The scalar factors of the elementary reflectors which\n*          represent the unitary matrix Q. See Further Details.\n*\n*  TAUP    (output) COMPLEX*16 array, dimension (min(M,N))\n*          The scalar factors of the elementary reflectors which\n*          represent the unitary matrix P. See Further Details.\n*\n*  WORK    (workspace/output) COMPLEX*16 array, dimension (MAX(1,LWORK))\n*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.\n*\n*  LWORK   (input) INTEGER\n*          The length of the array WORK.  LWORK >= max(1,M,N).\n*          For optimum performance LWORK >= (M+N)*NB, where NB\n*          is the optimal blocksize.\n*\n*          If LWORK = -1, then a workspace query is assumed; the routine\n*          only calculates the optimal size of the WORK array, returns\n*          this value as the first entry of the WORK array, and no error\n*          message related to LWORK is issued by XERBLA.\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*\n\n*  Further Details\n*  ===============\n*\n*  The matrices Q and P are represented as products of elementary\n*  reflectors:\n*\n*  If m >= n,\n*\n*     Q = H(1) H(2) . . . H(n)  and  P = G(1) G(2) . . . G(n-1)\n*\n*  Each H(i) and G(i) has the form:\n*\n*     H(i) = I - tauq * v * v'  and G(i) = I - taup * u * u'\n*\n*  where tauq and taup are complex scalars, and v and u are complex\n*  vectors; v(1:i-1) = 0, v(i) = 1, and v(i+1:m) is stored on exit in\n*  A(i+1:m,i); u(1:i) = 0, u(i+1) = 1, and u(i+2:n) is stored on exit in\n*  A(i,i+2:n); tauq is stored in TAUQ(i) and taup in TAUP(i).\n*\n*  If m < n,\n*\n*     Q = H(1) H(2) . . . H(m-1)  and  P = G(1) G(2) . . . G(m)\n*\n*  Each H(i) and G(i) has the form:\n*\n*     H(i) = I - tauq * v * v'  and G(i) = I - taup * u * u'\n*\n*  where tauq and taup are complex scalars, and v and u are complex\n*  vectors; v(1:i) = 0, v(i+1) = 1, and v(i+2:m) is stored on exit in\n*  A(i+2:m,i); u(1:i-1) = 0, u(i) = 1, and u(i+1:n) is stored on exit in\n*  A(i,i+1:n); tauq is stored in TAUQ(i) and taup in TAUP(i).\n*\n*  The contents of A on exit are illustrated by the following examples:\n*\n*  m = 6 and n = 5 (m > n):          m = 5 and n = 6 (m < n):\n*\n*    (  d   e   u1  u1  u1 )           (  d   u1  u1  u1  u1  u1 )\n*    (  v1  d   e   u2  u2 )           (  e   d   u2  u2  u2  u2 )\n*    (  v1  v2  d   e   u3 )           (  v1  e   d   u3  u3  u3 )\n*    (  v1  v2  v3  d   e  )           (  v1  v2  e   d   u4  u4 )\n*    (  v1  v2  v3  v4  d  )           (  v1  v2  v3  e   d   u5 )\n*    (  v1  v2  v3  v4  v5 )\n*\n*  where d and e denote diagonal and off-diagonal elements of B, vi\n*  denotes an element of the vector defining H(i), and ui an element of\n*  the vector defining G(i).\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  d, e, tauq, taup, work, info, a = NumRu::Lapack.zgebrd( m, a, lwork, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rblapack_m = argv[0];
  rblapack_a = argv[1];
  rblapack_lwork = argv[2];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_a);
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_DCOMPLEX)
    rblapack_a = na_change_type(rblapack_a, NA_DCOMPLEX);
  a = NA_PTR_TYPE(rblapack_a, doublecomplex*);
  m = NUM2INT(rblapack_m);
  lwork = NUM2INT(rblapack_lwork);
  {
    int shape[1];
    shape[0] = MIN(m,n);
    rblapack_d = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  d = NA_PTR_TYPE(rblapack_d, doublereal*);
  {
    int shape[1];
    shape[0] = MIN(m,n)-1;
    rblapack_e = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  e = NA_PTR_TYPE(rblapack_e, doublereal*);
  {
    int shape[1];
    shape[0] = MIN(m,n);
    rblapack_tauq = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  tauq = NA_PTR_TYPE(rblapack_tauq, doublecomplex*);
  {
    int shape[1];
    shape[0] = MIN(m,n);
    rblapack_taup = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  taup = NA_PTR_TYPE(rblapack_taup, doublecomplex*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rblapack_work = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rblapack_work, doublecomplex*);
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

  zgebrd_(&m, &n, a, &lda, d, e, tauq, taup, work, &lwork, &info);

  rblapack_info = INT2NUM(info);
  return rb_ary_new3(7, rblapack_d, rblapack_e, rblapack_tauq, rblapack_taup, rblapack_work, rblapack_info, rblapack_a);
}

void
init_lapack_zgebrd(VALUE mLapack){
  rb_define_module_function(mLapack, "zgebrd", rblapack_zgebrd, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
