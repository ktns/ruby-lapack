#include "rb_lapack.h"

extern VOID dsytd2_(char *uplo, integer *n, doublereal *a, integer *lda, doublereal *d, doublereal *e, doublereal *tau, integer *info);

static VALUE
rb_dsytd2(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_a;
  doublereal *a; 
  VALUE rb_d;
  doublereal *d; 
  VALUE rb_e;
  doublereal *e; 
  VALUE rb_tau;
  doublereal *tau; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  doublereal *a_out__;

  integer lda;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  d, e, tau, info, a = NumRu::Lapack.dsytd2( uplo, a)\n    or\n  NumRu::Lapack.dsytd2  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DSYTD2( UPLO, N, A, LDA, D, E, TAU, INFO )\n\n*  Purpose\n*  =======\n*\n*  DSYTD2 reduces a real symmetric matrix A to symmetric tridiagonal\n*  form T by an orthogonal similarity transformation: Q' * A * Q = T.\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO    (input) CHARACTER*1\n*          Specifies whether the upper or lower triangular part of the\n*          symmetric matrix A is stored:\n*          = 'U':  Upper triangular\n*          = 'L':  Lower triangular\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)\n*          On entry, the symmetric matrix A.  If UPLO = 'U', the leading\n*          n-by-n upper triangular part of A contains the upper\n*          triangular part of the matrix A, and the strictly lower\n*          triangular part of A is not referenced.  If UPLO = 'L', the\n*          leading n-by-n lower triangular part of A contains the lower\n*          triangular part of the matrix A, and the strictly upper\n*          triangular part of A is not referenced.\n*          On exit, if UPLO = 'U', the diagonal and first superdiagonal\n*          of A are overwritten by the corresponding elements of the\n*          tridiagonal matrix T, and the elements above the first\n*          superdiagonal, with the array TAU, represent the orthogonal\n*          matrix Q as a product of elementary reflectors; if UPLO\n*          = 'L', the diagonal and first subdiagonal of A are over-\n*          written by the corresponding elements of the tridiagonal\n*          matrix T, and the elements below the first subdiagonal, with\n*          the array TAU, represent the orthogonal matrix Q as a product\n*          of elementary reflectors. See Further Details.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,N).\n*\n*  D       (output) DOUBLE PRECISION array, dimension (N)\n*          The diagonal elements of the tridiagonal matrix T:\n*          D(i) = A(i,i).\n*\n*  E       (output) DOUBLE PRECISION array, dimension (N-1)\n*          The off-diagonal elements of the tridiagonal matrix T:\n*          E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'.\n*\n*  TAU     (output) DOUBLE PRECISION array, dimension (N-1)\n*          The scalar factors of the elementary reflectors (see Further\n*          Details).\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*\n\n*  Further Details\n*  ===============\n*\n*  If UPLO = 'U', the matrix Q is represented as a product of elementary\n*  reflectors\n*\n*     Q = H(n-1) . . . H(2) H(1).\n*\n*  Each H(i) has the form\n*\n*     H(i) = I - tau * v * v'\n*\n*  where tau is a real scalar, and v is a real vector with\n*  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in\n*  A(1:i-1,i+1), and tau in TAU(i).\n*\n*  If UPLO = 'L', the matrix Q is represented as a product of elementary\n*  reflectors\n*\n*     Q = H(1) H(2) . . . H(n-1).\n*\n*  Each H(i) has the form\n*\n*     H(i) = I - tau * v * v'\n*\n*  where tau is a real scalar, and v is a real vector with\n*  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i),\n*  and tau in TAU(i).\n*\n*  The contents of A on exit are illustrated by the following examples\n*  with n = 5:\n*\n*  if UPLO = 'U':                       if UPLO = 'L':\n*\n*    (  d   e   v2  v3  v4 )              (  d                  )\n*    (      d   e   v3  v4 )              (  e   d              )\n*    (          d   e   v4 )              (  v1  e   d          )\n*    (              d   e  )              (  v1  v2  e   d      )\n*    (                  d  )              (  v1  v2  v3  e   d  )\n*\n*  where d and e denote diagonal and off-diagonal elements of T, and vi\n*  denotes an element of the vector defining H(i).\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rb_uplo = argv[0];
  rb_a = argv[1];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_DFLOAT)
    rb_a = na_change_type(rb_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rb_a, doublereal*);
  uplo = StringValueCStr(rb_uplo)[0];
  {
    int shape[1];
    shape[0] = n;
    rb_d = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  d = NA_PTR_TYPE(rb_d, doublereal*);
  {
    int shape[1];
    shape[0] = n-1;
    rb_e = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  e = NA_PTR_TYPE(rb_e, doublereal*);
  {
    int shape[1];
    shape[0] = n-1;
    rb_tau = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  tau = NA_PTR_TYPE(rb_tau, doublereal*);
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

  dsytd2_(&uplo, &n, a, &lda, d, e, tau, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(5, rb_d, rb_e, rb_tau, rb_info, rb_a);
}

void
init_lapack_dsytd2(VALUE mLapack){
  rb_define_module_function(mLapack, "dsytd2", rb_dsytd2, -1);
}
