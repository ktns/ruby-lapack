#include "rb_lapack.h"

extern VOID dgehd2_(integer *n, integer *ilo, integer *ihi, doublereal *a, integer *lda, doublereal *tau, doublereal *work, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dgehd2(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_ilo;
  integer ilo; 
  VALUE rblapack_ihi;
  integer ihi; 
  VALUE rblapack_a;
  doublereal *a; 
  VALUE rblapack_tau;
  doublereal *tau; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_a_out__;
  doublereal *a_out__;
  doublereal *work;

  integer lda;
  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  tau, info, a = NumRu::Lapack.dgehd2( ilo, ihi, a, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DGEHD2( N, ILO, IHI, A, LDA, TAU, WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  DGEHD2 reduces a real general matrix A to upper Hessenberg form H by\n*  an orthogonal similarity transformation:  Q' * A * Q = H .\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  ILO     (input) INTEGER\n*  IHI     (input) INTEGER\n*          It is assumed that A is already upper triangular in rows\n*          and columns 1:ILO-1 and IHI+1:N. ILO and IHI are normally\n*          set by a previous call to DGEBAL; otherwise they should be\n*          set to 1 and N respectively. See Further Details.\n*          1 <= ILO <= IHI <= max(1,N).\n*\n*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)\n*          On entry, the n by n general matrix to be reduced.\n*          On exit, the upper triangle and the first subdiagonal of A\n*          are overwritten with the upper Hessenberg matrix H, and the\n*          elements below the first subdiagonal, with the array TAU,\n*          represent the orthogonal matrix Q as a product of elementary\n*          reflectors. See Further Details.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,N).\n*\n*  TAU     (output) DOUBLE PRECISION array, dimension (N-1)\n*          The scalar factors of the elementary reflectors (see Further\n*          Details).\n*\n*  WORK    (workspace) DOUBLE PRECISION array, dimension (N)\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*\n\n*  Further Details\n*  ===============\n*\n*  The matrix Q is represented as a product of (ihi-ilo) elementary\n*  reflectors\n*\n*     Q = H(ilo) H(ilo+1) . . . H(ihi-1).\n*\n*  Each H(i) has the form\n*\n*     H(i) = I - tau * v * v'\n*\n*  where tau is a real scalar, and v is a real vector with\n*  v(1:i) = 0, v(i+1) = 1 and v(ihi+1:n) = 0; v(i+2:ihi) is stored on\n*  exit in A(i+2:ihi,i), and tau in TAU(i).\n*\n*  The contents of A are illustrated by the following example, with\n*  n = 7, ilo = 2 and ihi = 6:\n*\n*  on entry,                        on exit,\n*\n*  ( a   a   a   a   a   a   a )    (  a   a   h   h   h   h   a )\n*  (     a   a   a   a   a   a )    (      a   h   h   h   h   a )\n*  (     a   a   a   a   a   a )    (      h   h   h   h   h   h )\n*  (     a   a   a   a   a   a )    (      v2  h   h   h   h   h )\n*  (     a   a   a   a   a   a )    (      v2  v3  h   h   h   h )\n*  (     a   a   a   a   a   a )    (      v2  v3  v4  h   h   h )\n*  (                         a )    (                          a )\n*\n*  where a denotes an element of the original matrix A, h denotes a\n*  modified element of the upper Hessenberg matrix H, and vi denotes an\n*  element of the vector defining H(i).\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  tau, info, a = NumRu::Lapack.dgehd2( ilo, ihi, a, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rblapack_ilo = argv[0];
  rblapack_ihi = argv[1];
  rblapack_a = argv[2];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_a);
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_DFLOAT)
    rblapack_a = na_change_type(rblapack_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rblapack_a, doublereal*);
  ilo = NUM2INT(rblapack_ilo);
  ihi = NUM2INT(rblapack_ihi);
  {
    int shape[1];
    shape[0] = n-1;
    rblapack_tau = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  tau = NA_PTR_TYPE(rblapack_tau, doublereal*);
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
  work = ALLOC_N(doublereal, (n));

  dgehd2_(&n, &ilo, &ihi, a, &lda, tau, work, &info);

  free(work);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(3, rblapack_tau, rblapack_info, rblapack_a);
}

void
init_lapack_dgehd2(VALUE mLapack){
  rb_define_module_function(mLapack, "dgehd2", rblapack_dgehd2, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
