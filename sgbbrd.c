#include "rb_lapack.h"

static VALUE
rb_sgbbrd(int argc, VALUE *argv, VALUE self){
  VALUE rb_vect;
  char vect; 
  VALUE rb_kl;
  integer kl; 
  VALUE rb_ku;
  integer ku; 
  VALUE rb_ab;
  real *ab; 
  VALUE rb_c;
  real *c; 
  VALUE rb_d;
  real *d; 
  VALUE rb_e;
  real *e; 
  VALUE rb_q;
  real *q; 
  VALUE rb_pt;
  real *pt; 
  VALUE rb_info;
  integer info; 
  VALUE rb_ab_out__;
  real *ab_out__;
  VALUE rb_c_out__;
  real *c_out__;
  real *work;

  integer ldab;
  integer n;
  integer ldc;
  integer ncc;
  integer m;
  integer ldq;
  integer ldpt;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  d, e, q, pt, info, ab, c = NumRu::Lapack.sgbbrd( vect, kl, ku, ab, c)\n    or\n  NumRu::Lapack.sgbbrd  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SGBBRD( VECT, M, N, NCC, KL, KU, AB, LDAB, D, E, Q, LDQ, PT, LDPT, C, LDC, WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  SGBBRD reduces a real general m-by-n band matrix A to upper\n*  bidiagonal form B by an orthogonal transformation: Q' * A * P = B.\n*\n*  The routine computes B, and optionally forms Q or P', or computes\n*  Q'*C for a given matrix C.\n*\n\n*  Arguments\n*  =========\n*\n*  VECT    (input) CHARACTER*1\n*          Specifies whether or not the matrices Q and P' are to be\n*          formed.\n*          = 'N': do not form Q or P';\n*          = 'Q': form Q only;\n*          = 'P': form P' only;\n*          = 'B': form both.\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix A.  M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix A.  N >= 0.\n*\n*  NCC     (input) INTEGER\n*          The number of columns of the matrix C.  NCC >= 0.\n*\n*  KL      (input) INTEGER\n*          The number of subdiagonals of the matrix A. KL >= 0.\n*\n*  KU      (input) INTEGER\n*          The number of superdiagonals of the matrix A. KU >= 0.\n*\n*  AB      (input/output) REAL array, dimension (LDAB,N)\n*          On entry, the m-by-n band matrix A, stored in rows 1 to\n*          KL+KU+1. The j-th column of A is stored in the j-th column of\n*          the array AB as follows:\n*          AB(ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl).\n*          On exit, A is overwritten by values generated during the\n*          reduction.\n*\n*  LDAB    (input) INTEGER\n*          The leading dimension of the array A. LDAB >= KL+KU+1.\n*\n*  D       (output) REAL array, dimension (min(M,N))\n*          The diagonal elements of the bidiagonal matrix B.\n*\n*  E       (output) REAL array, dimension (min(M,N)-1)\n*          The superdiagonal elements of the bidiagonal matrix B.\n*\n*  Q       (output) REAL array, dimension (LDQ,M)\n*          If VECT = 'Q' or 'B', the m-by-m orthogonal matrix Q.\n*          If VECT = 'N' or 'P', the array Q is not referenced.\n*\n*  LDQ     (input) INTEGER\n*          The leading dimension of the array Q.\n*          LDQ >= max(1,M) if VECT = 'Q' or 'B'; LDQ >= 1 otherwise.\n*\n*  PT      (output) REAL array, dimension (LDPT,N)\n*          If VECT = 'P' or 'B', the n-by-n orthogonal matrix P'.\n*          If VECT = 'N' or 'Q', the array PT is not referenced.\n*\n*  LDPT    (input) INTEGER\n*          The leading dimension of the array PT.\n*          LDPT >= max(1,N) if VECT = 'P' or 'B'; LDPT >= 1 otherwise.\n*\n*  C       (input/output) REAL array, dimension (LDC,NCC)\n*          On entry, an m-by-ncc matrix C.\n*          On exit, C is overwritten by Q'*C.\n*          C is not referenced if NCC = 0.\n*\n*  LDC     (input) INTEGER\n*          The leading dimension of the array C.\n*          LDC >= max(1,M) if NCC > 0; LDC >= 1 if NCC = 0.\n*\n*  WORK    (workspace) REAL array, dimension (2*max(M,N))\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_vect = argv[0];
  rb_kl = argv[1];
  rb_ku = argv[2];
  rb_ab = argv[3];
  rb_c = argv[4];

  vect = StringValueCStr(rb_vect)[0];
  kl = NUM2INT(rb_kl);
  ku = NUM2INT(rb_ku);
  if (!NA_IsNArray(rb_ab))
    rb_raise(rb_eArgError, "ab (4th argument) must be NArray");
  if (NA_RANK(rb_ab) != 2)
    rb_raise(rb_eArgError, "rank of ab (4th argument) must be %d", 2);
  ldab = NA_SHAPE0(rb_ab);
  n = NA_SHAPE1(rb_ab);
  if (NA_TYPE(rb_ab) != NA_SFLOAT)
    rb_ab = na_change_type(rb_ab, NA_SFLOAT);
  ab = NA_PTR_TYPE(rb_ab, real*);
  if (!NA_IsNArray(rb_c))
    rb_raise(rb_eArgError, "c (5th argument) must be NArray");
  if (NA_RANK(rb_c) != 2)
    rb_raise(rb_eArgError, "rank of c (5th argument) must be %d", 2);
  ldc = NA_SHAPE0(rb_c);
  ncc = NA_SHAPE1(rb_c);
  if (NA_TYPE(rb_c) != NA_SFLOAT)
    rb_c = na_change_type(rb_c, NA_SFLOAT);
  c = NA_PTR_TYPE(rb_c, real*);
  m = ldab;
  {
    int shape[1];
    shape[0] = MIN(m,n);
    rb_d = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  d = NA_PTR_TYPE(rb_d, real*);
  {
    int shape[1];
    shape[0] = MIN(m,n)-1;
    rb_e = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  e = NA_PTR_TYPE(rb_e, real*);
  ldq = ((lsame_(&vect,"Q")) || (lsame_(&vect,"B"))) ? MAX(1,m) : 1;
  {
    int shape[2];
    shape[0] = ldq;
    shape[1] = m;
    rb_q = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  q = NA_PTR_TYPE(rb_q, real*);
  ldpt = ((lsame_(&vect,"P")) || (lsame_(&vect,"B"))) ? MAX(1,n) : 1;
  {
    int shape[2];
    shape[0] = ldpt;
    shape[1] = n;
    rb_pt = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  pt = NA_PTR_TYPE(rb_pt, real*);
  {
    int shape[2];
    shape[0] = ldab;
    shape[1] = n;
    rb_ab_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  ab_out__ = NA_PTR_TYPE(rb_ab_out__, real*);
  MEMCPY(ab_out__, ab, real, NA_TOTAL(rb_ab));
  rb_ab = rb_ab_out__;
  ab = ab_out__;
  {
    int shape[2];
    shape[0] = ldc;
    shape[1] = ncc;
    rb_c_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  c_out__ = NA_PTR_TYPE(rb_c_out__, real*);
  MEMCPY(c_out__, c, real, NA_TOTAL(rb_c));
  rb_c = rb_c_out__;
  c = c_out__;
  work = ALLOC_N(real, (2*MAX(m,n)));

  sgbbrd_(&vect, &m, &n, &ncc, &kl, &ku, ab, &ldab, d, e, q, &ldq, pt, &ldpt, c, &ldc, work, &info);

  free(work);
  rb_info = INT2NUM(info);
  return rb_ary_new3(7, rb_d, rb_e, rb_q, rb_pt, rb_info, rb_ab, rb_c);
}

void
init_lapack_sgbbrd(VALUE mLapack){
  rb_define_module_function(mLapack, "sgbbrd", rb_sgbbrd, -1);
}
