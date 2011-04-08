#include "rb_lapack.h"

extern VOID cgbbrd_(char *vect, integer *m, integer *n, integer *ncc, integer *kl, integer *ku, complex *ab, integer *ldab, real *d, real *e, complex *q, integer *ldq, complex *pt, integer *ldpt, complex *c, integer *ldc, complex *work, real *rwork, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_cgbbrd(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_vect;
  char vect; 
  VALUE rblapack_kl;
  integer kl; 
  VALUE rblapack_ku;
  integer ku; 
  VALUE rblapack_ab;
  complex *ab; 
  VALUE rblapack_c;
  complex *c; 
  VALUE rblapack_d;
  real *d; 
  VALUE rblapack_e;
  real *e; 
  VALUE rblapack_q;
  complex *q; 
  VALUE rblapack_pt;
  complex *pt; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_ab_out__;
  complex *ab_out__;
  VALUE rblapack_c_out__;
  complex *c_out__;
  complex *work;
  real *rwork;

  integer ldab;
  integer n;
  integer ldc;
  integer ncc;
  integer ldq;
  integer m;
  integer ldpt;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  d, e, q, pt, info, ab, c = NumRu::Lapack.cgbbrd( vect, kl, ku, ab, c, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE CGBBRD( VECT, M, N, NCC, KL, KU, AB, LDAB, D, E, Q, LDQ, PT, LDPT, C, LDC, WORK, RWORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  CGBBRD reduces a complex general m-by-n band matrix A to real upper\n*  bidiagonal form B by a unitary transformation: Q' * A * P = B.\n*\n*  The routine computes B, and optionally forms Q or P', or computes\n*  Q'*C for a given matrix C.\n*\n\n*  Arguments\n*  =========\n*\n*  VECT    (input) CHARACTER*1\n*          Specifies whether or not the matrices Q and P' are to be\n*          formed.\n*          = 'N': do not form Q or P';\n*          = 'Q': form Q only;\n*          = 'P': form P' only;\n*          = 'B': form both.\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix A.  M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix A.  N >= 0.\n*\n*  NCC     (input) INTEGER\n*          The number of columns of the matrix C.  NCC >= 0.\n*\n*  KL      (input) INTEGER\n*          The number of subdiagonals of the matrix A. KL >= 0.\n*\n*  KU      (input) INTEGER\n*          The number of superdiagonals of the matrix A. KU >= 0.\n*\n*  AB      (input/output) COMPLEX array, dimension (LDAB,N)\n*          On entry, the m-by-n band matrix A, stored in rows 1 to\n*          KL+KU+1. The j-th column of A is stored in the j-th column of\n*          the array AB as follows:\n*          AB(ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl).\n*          On exit, A is overwritten by values generated during the\n*          reduction.\n*\n*  LDAB    (input) INTEGER\n*          The leading dimension of the array A. LDAB >= KL+KU+1.\n*\n*  D       (output) REAL array, dimension (min(M,N))\n*          The diagonal elements of the bidiagonal matrix B.\n*\n*  E       (output) REAL array, dimension (min(M,N)-1)\n*          The superdiagonal elements of the bidiagonal matrix B.\n*\n*  Q       (output) COMPLEX array, dimension (LDQ,M)\n*          If VECT = 'Q' or 'B', the m-by-m unitary matrix Q.\n*          If VECT = 'N' or 'P', the array Q is not referenced.\n*\n*  LDQ     (input) INTEGER\n*          The leading dimension of the array Q.\n*          LDQ >= max(1,M) if VECT = 'Q' or 'B'; LDQ >= 1 otherwise.\n*\n*  PT      (output) COMPLEX array, dimension (LDPT,N)\n*          If VECT = 'P' or 'B', the n-by-n unitary matrix P'.\n*          If VECT = 'N' or 'Q', the array PT is not referenced.\n*\n*  LDPT    (input) INTEGER\n*          The leading dimension of the array PT.\n*          LDPT >= max(1,N) if VECT = 'P' or 'B'; LDPT >= 1 otherwise.\n*\n*  C       (input/output) COMPLEX array, dimension (LDC,NCC)\n*          On entry, an m-by-ncc matrix C.\n*          On exit, C is overwritten by Q'*C.\n*          C is not referenced if NCC = 0.\n*\n*  LDC     (input) INTEGER\n*          The leading dimension of the array C.\n*          LDC >= max(1,M) if NCC > 0; LDC >= 1 if NCC = 0.\n*\n*  WORK    (workspace) COMPLEX array, dimension (max(M,N))\n*\n*  RWORK   (workspace) REAL array, dimension (max(M,N))\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  d, e, q, pt, info, ab, c = NumRu::Lapack.cgbbrd( vect, kl, ku, ab, c, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rblapack_vect = argv[0];
  rblapack_kl = argv[1];
  rblapack_ku = argv[2];
  rblapack_ab = argv[3];
  rblapack_c = argv[4];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_ab))
    rb_raise(rb_eArgError, "ab (4th argument) must be NArray");
  if (NA_RANK(rblapack_ab) != 2)
    rb_raise(rb_eArgError, "rank of ab (4th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_ab);
  ldab = NA_SHAPE0(rblapack_ab);
  if (NA_TYPE(rblapack_ab) != NA_SCOMPLEX)
    rblapack_ab = na_change_type(rblapack_ab, NA_SCOMPLEX);
  ab = NA_PTR_TYPE(rblapack_ab, complex*);
  kl = NUM2INT(rblapack_kl);
  if (!NA_IsNArray(rblapack_c))
    rb_raise(rb_eArgError, "c (5th argument) must be NArray");
  if (NA_RANK(rblapack_c) != 2)
    rb_raise(rb_eArgError, "rank of c (5th argument) must be %d", 2);
  ncc = NA_SHAPE1(rblapack_c);
  ldc = NA_SHAPE0(rblapack_c);
  if (NA_TYPE(rblapack_c) != NA_SCOMPLEX)
    rblapack_c = na_change_type(rblapack_c, NA_SCOMPLEX);
  c = NA_PTR_TYPE(rblapack_c, complex*);
  vect = StringValueCStr(rblapack_vect)[0];
  ku = NUM2INT(rblapack_ku);
  m = ldab;
  ldpt = ((lsame_(&vect,"P")) || (lsame_(&vect,"B"))) ? MAX(1,n) : 1;
  ldq = ((lsame_(&vect,"Q")) || (lsame_(&vect,"B"))) ? MAX(1,m) : 1;
  {
    int shape[1];
    shape[0] = MIN(m,n);
    rblapack_d = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  d = NA_PTR_TYPE(rblapack_d, real*);
  {
    int shape[1];
    shape[0] = MIN(m,n)-1;
    rblapack_e = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  e = NA_PTR_TYPE(rblapack_e, real*);
  {
    int shape[2];
    shape[0] = ldq;
    shape[1] = m;
    rblapack_q = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  q = NA_PTR_TYPE(rblapack_q, complex*);
  {
    int shape[2];
    shape[0] = ldpt;
    shape[1] = n;
    rblapack_pt = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  pt = NA_PTR_TYPE(rblapack_pt, complex*);
  {
    int shape[2];
    shape[0] = ldab;
    shape[1] = n;
    rblapack_ab_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  ab_out__ = NA_PTR_TYPE(rblapack_ab_out__, complex*);
  MEMCPY(ab_out__, ab, complex, NA_TOTAL(rblapack_ab));
  rblapack_ab = rblapack_ab_out__;
  ab = ab_out__;
  {
    int shape[2];
    shape[0] = ldc;
    shape[1] = ncc;
    rblapack_c_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  c_out__ = NA_PTR_TYPE(rblapack_c_out__, complex*);
  MEMCPY(c_out__, c, complex, NA_TOTAL(rblapack_c));
  rblapack_c = rblapack_c_out__;
  c = c_out__;
  work = ALLOC_N(complex, (MAX(m,n)));
  rwork = ALLOC_N(real, (MAX(m,n)));

  cgbbrd_(&vect, &m, &n, &ncc, &kl, &ku, ab, &ldab, d, e, q, &ldq, pt, &ldpt, c, &ldc, work, rwork, &info);

  free(work);
  free(rwork);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(7, rblapack_d, rblapack_e, rblapack_q, rblapack_pt, rblapack_info, rblapack_ab, rblapack_c);
}

void
init_lapack_cgbbrd(VALUE mLapack){
  rb_define_module_function(mLapack, "cgbbrd", rblapack_cgbbrd, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
