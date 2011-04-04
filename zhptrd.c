#include "rb_lapack.h"

extern VOID zhptrd_(char *uplo, integer *n, doublecomplex *ap, doublereal *d, doublereal *e, doublecomplex *tau, integer *info);

static VALUE
rb_zhptrd(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_ap;
  doublecomplex *ap; 
  VALUE rb_d;
  doublereal *d; 
  VALUE rb_e;
  doublereal *e; 
  VALUE rb_tau;
  doublecomplex *tau; 
  VALUE rb_info;
  integer info; 
  VALUE rb_ap_out__;
  doublecomplex *ap_out__;

  integer ldap;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  d, e, tau, info, ap = NumRu::Lapack.zhptrd( uplo, ap)\n    or\n  NumRu::Lapack.zhptrd  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZHPTRD( UPLO, N, AP, D, E, TAU, INFO )\n\n*  Purpose\n*  =======\n*\n*  ZHPTRD reduces a complex Hermitian matrix A stored in packed form to\n*  real symmetric tridiagonal form T by a unitary similarity\n*  transformation: Q**H * A * Q = T.\n*\n\n*  Arguments\n*  =========\n*\n*  UPLO    (input) CHARACTER*1\n*          = 'U':  Upper triangle of A is stored;\n*          = 'L':  Lower triangle of A is stored.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  AP      (input/output) COMPLEX*16 array, dimension (N*(N+1)/2)\n*          On entry, the upper or lower triangle of the Hermitian matrix\n*          A, packed columnwise in a linear array.  The j-th column of A\n*          is stored in the array AP as follows:\n*          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;\n*          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.\n*          On exit, if UPLO = 'U', the diagonal and first superdiagonal\n*          of A are overwritten by the corresponding elements of the\n*          tridiagonal matrix T, and the elements above the first\n*          superdiagonal, with the array TAU, represent the unitary\n*          matrix Q as a product of elementary reflectors; if UPLO\n*          = 'L', the diagonal and first subdiagonal of A are over-\n*          written by the corresponding elements of the tridiagonal\n*          matrix T, and the elements below the first subdiagonal, with\n*          the array TAU, represent the unitary matrix Q as a product\n*          of elementary reflectors. See Further Details.\n*\n*  D       (output) DOUBLE PRECISION array, dimension (N)\n*          The diagonal elements of the tridiagonal matrix T:\n*          D(i) = A(i,i).\n*\n*  E       (output) DOUBLE PRECISION array, dimension (N-1)\n*          The off-diagonal elements of the tridiagonal matrix T:\n*          E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'.\n*\n*  TAU     (output) COMPLEX*16 array, dimension (N-1)\n*          The scalar factors of the elementary reflectors (see Further\n*          Details).\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value\n*\n\n*  Further Details\n*  ===============\n*\n*  If UPLO = 'U', the matrix Q is represented as a product of elementary\n*  reflectors\n*\n*     Q = H(n-1) . . . H(2) H(1).\n*\n*  Each H(i) has the form\n*\n*     H(i) = I - tau * v * v'\n*\n*  where tau is a complex scalar, and v is a complex vector with\n*  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in AP,\n*  overwriting A(1:i-1,i+1), and tau is stored in TAU(i).\n*\n*  If UPLO = 'L', the matrix Q is represented as a product of elementary\n*  reflectors\n*\n*     Q = H(1) H(2) . . . H(n-1).\n*\n*  Each H(i) has the form\n*\n*     H(i) = I - tau * v * v'\n*\n*  where tau is a complex scalar, and v is a complex vector with\n*  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in AP,\n*  overwriting A(i+2:n,i), and tau is stored in TAU(i).\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rb_uplo = argv[0];
  rb_ap = argv[1];

  uplo = StringValueCStr(rb_uplo)[0];
  if (!NA_IsNArray(rb_ap))
    rb_raise(rb_eArgError, "ap (2th argument) must be NArray");
  if (NA_RANK(rb_ap) != 1)
    rb_raise(rb_eArgError, "rank of ap (2th argument) must be %d", 1);
  ldap = NA_SHAPE0(rb_ap);
  if (NA_TYPE(rb_ap) != NA_DCOMPLEX)
    rb_ap = na_change_type(rb_ap, NA_DCOMPLEX);
  ap = NA_PTR_TYPE(rb_ap, doublecomplex*);
  n = ((int)sqrtf(ldap*8+1.0f)-1)/2;
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
    rb_tau = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  tau = NA_PTR_TYPE(rb_tau, doublecomplex*);
  {
    int shape[1];
    shape[0] = ldap;
    rb_ap_out__ = na_make_object(NA_DCOMPLEX, 1, shape, cNArray);
  }
  ap_out__ = NA_PTR_TYPE(rb_ap_out__, doublecomplex*);
  MEMCPY(ap_out__, ap, doublecomplex, NA_TOTAL(rb_ap));
  rb_ap = rb_ap_out__;
  ap = ap_out__;

  zhptrd_(&uplo, &n, ap, d, e, tau, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(5, rb_d, rb_e, rb_tau, rb_info, rb_ap);
}

void
init_lapack_zhptrd(VALUE mLapack){
  rb_define_module_function(mLapack, "zhptrd", rb_zhptrd, -1);
}
