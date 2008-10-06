#include "rb_lapack.h"

static VALUE
rb_zgbtrf(int argc, VALUE *argv, VALUE self){
  VALUE rb_m;
  integer m; 
  VALUE rb_kl;
  integer kl; 
  VALUE rb_ku;
  integer ku; 
  VALUE rb_ab;
  doublecomplex *ab; 
  VALUE rb_ipiv;
  integer *ipiv; 
  VALUE rb_info;
  integer info; 
  VALUE rb_ab_out__;
  doublecomplex *ab_out__;

  integer ldab;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  ipiv, info, ab = NumRu::Lapack.zgbtrf( m, kl, ku, ab)\n    or\n  NumRu::Lapack.zgbtrf  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZGBTRF( M, N, KL, KU, AB, LDAB, IPIV, INFO )\n\n*  Purpose\n*  =======\n*\n*  ZGBTRF computes an LU factorization of a complex m-by-n band matrix A\n*  using partial pivoting with row interchanges.\n*\n*  This is the blocked version of the algorithm, calling Level 3 BLAS.\n*\n\n*  Arguments\n*  =========\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix A.  M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix A.  N >= 0.\n*\n*  KL      (input) INTEGER\n*          The number of subdiagonals within the band of A.  KL >= 0.\n*\n*  KU      (input) INTEGER\n*          The number of superdiagonals within the band of A.  KU >= 0.\n*\n*  AB      (input/output) COMPLEX*16 array, dimension (LDAB,N)\n*          On entry, the matrix A in band storage, in rows KL+1 to\n*          2*KL+KU+1; rows 1 to KL of the array need not be set.\n*          The j-th column of A is stored in the j-th column of the\n*          array AB as follows:\n*          AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)\n*\n*          On exit, details of the factorization: U is stored as an\n*          upper triangular band matrix with KL+KU superdiagonals in\n*          rows 1 to KL+KU+1, and the multipliers used during the\n*          factorization are stored in rows KL+KU+2 to 2*KL+KU+1.\n*          See below for further details.\n*\n*  LDAB    (input) INTEGER\n*          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.\n*\n*  IPIV    (output) INTEGER array, dimension (min(M,N))\n*          The pivot indices; for 1 <= i <= min(M,N), row i of the\n*          matrix was interchanged with row IPIV(i).\n*\n*  INFO    (output) INTEGER\n*          = 0: successful exit\n*          < 0: if INFO = -i, the i-th argument had an illegal value\n*          > 0: if INFO = +i, U(i,i) is exactly zero. The factorization\n*               has been completed, but the factor U is exactly\n*               singular, and division by zero will occur if it is used\n*               to solve a system of equations.\n*\n\n*  Further Details\n*  ===============\n*\n*  The band storage scheme is illustrated by the following example, when\n*  M = N = 6, KL = 2, KU = 1:\n*\n*  On entry:                       On exit:\n*\n*      *    *    *    +    +    +       *    *    *   u14  u25  u36\n*      *    *    +    +    +    +       *    *   u13  u24  u35  u46\n*      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56\n*     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66\n*     a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   *\n*     a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    *\n*\n*  Array elements marked * are not used by the routine; elements marked\n*  + need not be set on entry, but are required by the routine to store\n*  elements of U because of fill-in resulting from the row interchanges.\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_m = argv[0];
  rb_kl = argv[1];
  rb_ku = argv[2];
  rb_ab = argv[3];

  m = NUM2INT(rb_m);
  kl = NUM2INT(rb_kl);
  ku = NUM2INT(rb_ku);
  if (!NA_IsNArray(rb_ab))
    rb_raise(rb_eArgError, "ab (4th argument) must be NArray");
  if (NA_RANK(rb_ab) != 2)
    rb_raise(rb_eArgError, "rank of ab (4th argument) must be %d", 2);
  ldab = NA_SHAPE0(rb_ab);
  n = NA_SHAPE1(rb_ab);
  if (NA_TYPE(rb_ab) != NA_DCOMPLEX)
    rb_ab = na_change_type(rb_ab, NA_DCOMPLEX);
  ab = NA_PTR_TYPE(rb_ab, doublecomplex*);
  {
    int shape[1];
    shape[0] = MIN(m,n);
    rb_ipiv = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  ipiv = NA_PTR_TYPE(rb_ipiv, integer*);
  {
    int shape[2];
    shape[0] = ldab;
    shape[1] = n;
    rb_ab_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  ab_out__ = NA_PTR_TYPE(rb_ab_out__, doublecomplex*);
  MEMCPY(ab_out__, ab, doublecomplex, NA_TOTAL(rb_ab));
  rb_ab = rb_ab_out__;
  ab = ab_out__;

  zgbtrf_(&m, &n, &kl, &ku, ab, &ldab, ipiv, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(3, rb_ipiv, rb_info, rb_ab);
}

void
init_lapack_zgbtrf(VALUE mLapack){
  rb_define_module_function(mLapack, "zgbtrf", rb_zgbtrf, -1);
}
