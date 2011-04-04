#include "rb_lapack.h"

extern VOID dgbbrd_(char *vect, integer *m, integer *n, integer *ncc, integer *kl, integer *ku, doublereal *ab, integer *ldab, doublereal *d, doublereal *e, doublereal *q, integer *ldq, doublereal *pt, integer *ldpt, doublereal *c, integer *ldc, doublereal *work, integer *info);

static VALUE
rb_dgbbrd(int argc, VALUE *argv, VALUE self){
  VALUE rb_vect;
  char vect; 
  VALUE rb_kl;
  integer kl; 
  VALUE rb_ku;
  integer ku; 
  VALUE rb_ab;
  doublereal *ab; 
  VALUE rb_c;
  doublereal *c; 
  VALUE rb_d;
  doublereal *d; 
  VALUE rb_e;
  doublereal *e; 
  VALUE rb_q;
  doublereal *q; 
  VALUE rb_pt;
  doublereal *pt; 
  VALUE rb_info;
  integer info; 
  VALUE rb_ab_out__;
  doublereal *ab_out__;
  VALUE rb_c_out__;
  doublereal *c_out__;
  doublereal *work;

  integer ldab;
  integer n;
  integer ldc;
  integer ncc;
  integer ldq;
  integer m;
  integer ldpt;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  d, e, q, pt, info, ab, c = NumRu::Lapack.dgbbrd( vect, kl, ku, ab, c)\n    or\n  NumRu::Lapack.dgbbrd  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DGBBRD( VECT, M, N, NCC, KL, KU, AB, LDAB, D, E, Q, LDQ, PT, LDPT, C, LDC, WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  DGBBRD reduces a real general m-by-n band matrix A to upper\n*  bidiagonal form B by an orthogonal transformation: Q' * A * P = B.\n*\n*  The routine computes B, and optionally forms Q or P', or computes\n*  Q'*C for a given matrix C.\n*\n\n*  Arguments\n*  =========\n*\n*  VECT    (input) CHARACTER*1\n*          Specifies whether or not the matrices Q and P' are to be\n*          formed.\n*          = 'N': do not form Q or P';\n*          = 'Q': form Q only;\n*          = 'P': form P' only;\n*          = 'B': form both.\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix A.  M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix A.  N >= 0.\n*\n*  NCC     (input) INTEGER\n*          The number of columns of the matrix C.  NCC >= 0.\n*\n*  KL      (input) INTEGER\n*          The number of subdiagonals of the matrix A. KL >= 0.\n*\n*  KU      (input) INTEGER\n*          The number of superdiagonals of the matrix A. KU >= 0.\n*\n*  AB      (input/output) DOUBLE PRECISION array, dimension (LDAB,N)\n*          On entry, the m-by-n band matrix A, stored in rows 1 to\n*          KL+KU+1. The j-th column of A is stored in the j-th column of\n*          the array AB as follows:\n*          AB(ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl).\n*          On exit, A is overwritten by values generated during the\n*          reduction.\n*\n*  LDAB    (input) INTEGER\n*          The leading dimension of the array A. LDAB >= KL+KU+1.\n*\n*  D       (output) DOUBLE PRECISION array, dimension (min(M,N))\n*          The diagonal elements of the bidiagonal matrix B.\n*\n*  E       (output) DOUBLE PRECISION array, dimension (min(M,N)-1)\n*          The superdiagonal elements of the bidiagonal matrix B.\n*\n*  Q       (output) DOUBLE PRECISION array, dimension (LDQ,M)\n*          If VECT = 'Q' or 'B', the m-by-m orthogonal matrix Q.\n*          If VECT = 'N' or 'P', the array Q is not referenced.\n*\n*  LDQ     (input) INTEGER\n*          The leading dimension of the array Q.\n*          LDQ >= max(1,M) if VECT = 'Q' or 'B'; LDQ >= 1 otherwise.\n*\n*  PT      (output) DOUBLE PRECISION array, dimension (LDPT,N)\n*          If VECT = 'P' or 'B', the n-by-n orthogonal matrix P'.\n*          If VECT = 'N' or 'Q', the array PT is not referenced.\n*\n*  LDPT    (input) INTEGER\n*          The leading dimension of the array PT.\n*          LDPT >= max(1,N) if VECT = 'P' or 'B'; LDPT >= 1 otherwise.\n*\n*  C       (input/output) DOUBLE PRECISION array, dimension (LDC,NCC)\n*          On entry, an m-by-ncc matrix C.\n*          On exit, C is overwritten by Q'*C.\n*          C is not referenced if NCC = 0.\n*\n*  LDC     (input) INTEGER\n*          The leading dimension of the array C.\n*          LDC >= max(1,M) if NCC > 0; LDC >= 1 if NCC = 0.\n*\n*  WORK    (workspace) DOUBLE PRECISION array, dimension (2*max(M,N))\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_vect = argv[0];
  rb_kl = argv[1];
  rb_ku = argv[2];
  rb_ab = argv[3];
  rb_c = argv[4];

  if (!NA_IsNArray(rb_ab))
    rb_raise(rb_eArgError, "ab (4th argument) must be NArray");
  if (NA_RANK(rb_ab) != 2)
    rb_raise(rb_eArgError, "rank of ab (4th argument) must be %d", 2);
  n = NA_SHAPE1(rb_ab);
  ldab = NA_SHAPE0(rb_ab);
  if (NA_TYPE(rb_ab) != NA_DFLOAT)
    rb_ab = na_change_type(rb_ab, NA_DFLOAT);
  ab = NA_PTR_TYPE(rb_ab, doublereal*);
  kl = NUM2INT(rb_kl);
  if (!NA_IsNArray(rb_c))
    rb_raise(rb_eArgError, "c (5th argument) must be NArray");
  if (NA_RANK(rb_c) != 2)
    rb_raise(rb_eArgError, "rank of c (5th argument) must be %d", 2);
  ncc = NA_SHAPE1(rb_c);
  ldc = NA_SHAPE0(rb_c);
  if (NA_TYPE(rb_c) != NA_DFLOAT)
    rb_c = na_change_type(rb_c, NA_DFLOAT);
  c = NA_PTR_TYPE(rb_c, doublereal*);
  vect = StringValueCStr(rb_vect)[0];
  ku = NUM2INT(rb_ku);
  m = ldab;
  ldpt = ((lsame_(&vect,"P")) || (lsame_(&vect,"B"))) ? MAX(1,n) : 1;
  ldq = ((lsame_(&vect,"Q")) || (lsame_(&vect,"B"))) ? MAX(1,m) : 1;
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
    int shape[2];
    shape[0] = ldq;
    shape[1] = m;
    rb_q = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  q = NA_PTR_TYPE(rb_q, doublereal*);
  {
    int shape[2];
    shape[0] = ldpt;
    shape[1] = n;
    rb_pt = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  pt = NA_PTR_TYPE(rb_pt, doublereal*);
  {
    int shape[2];
    shape[0] = ldab;
    shape[1] = n;
    rb_ab_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  ab_out__ = NA_PTR_TYPE(rb_ab_out__, doublereal*);
  MEMCPY(ab_out__, ab, doublereal, NA_TOTAL(rb_ab));
  rb_ab = rb_ab_out__;
  ab = ab_out__;
  {
    int shape[2];
    shape[0] = ldc;
    shape[1] = ncc;
    rb_c_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  c_out__ = NA_PTR_TYPE(rb_c_out__, doublereal*);
  MEMCPY(c_out__, c, doublereal, NA_TOTAL(rb_c));
  rb_c = rb_c_out__;
  c = c_out__;
  work = ALLOC_N(doublereal, (2*MAX(m,n)));

  dgbbrd_(&vect, &m, &n, &ncc, &kl, &ku, ab, &ldab, d, e, q, &ldq, pt, &ldpt, c, &ldc, work, &info);

  free(work);
  rb_info = INT2NUM(info);
  return rb_ary_new3(7, rb_d, rb_e, rb_q, rb_pt, rb_info, rb_ab, rb_c);
}

void
init_lapack_dgbbrd(VALUE mLapack){
  rb_define_module_function(mLapack, "dgbbrd", rb_dgbbrd, -1);
}
