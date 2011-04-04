#include "rb_lapack.h"

extern VOID cgehd2_(integer *n, integer *ilo, integer *ihi, complex *a, integer *lda, complex *tau, complex *work, integer *info);

static VALUE
rb_cgehd2(int argc, VALUE *argv, VALUE self){
  VALUE rb_ilo;
  integer ilo; 
  VALUE rb_ihi;
  integer ihi; 
  VALUE rb_a;
  complex *a; 
  VALUE rb_tau;
  complex *tau; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  complex *a_out__;
  complex *work;

  integer lda;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  tau, info, a = NumRu::Lapack.cgehd2( ilo, ihi, a)\n    or\n  NumRu::Lapack.cgehd2  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE CGEHD2( N, ILO, IHI, A, LDA, TAU, WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  CGEHD2 reduces a complex general matrix A to upper Hessenberg form H\n*  by a unitary similarity transformation:  Q' * A * Q = H .\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.\n*\n*  ILO     (input) INTEGER\n*  IHI     (input) INTEGER\n*          It is assumed that A is already upper triangular in rows\n*          and columns 1:ILO-1 and IHI+1:N. ILO and IHI are normally\n*          set by a previous call to CGEBAL; otherwise they should be\n*          set to 1 and N respectively. See Further Details.\n*          1 <= ILO <= IHI <= max(1,N).\n*\n*  A       (input/output) COMPLEX array, dimension (LDA,N)\n*          On entry, the n by n general matrix to be reduced.\n*          On exit, the upper triangle and the first subdiagonal of A\n*          are overwritten with the upper Hessenberg matrix H, and the\n*          elements below the first subdiagonal, with the array TAU,\n*          represent the unitary matrix Q as a product of elementary\n*          reflectors. See Further Details.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,N).\n*\n*  TAU     (output) COMPLEX array, dimension (N-1)\n*          The scalar factors of the elementary reflectors (see Further\n*          Details).\n*\n*  WORK    (workspace) COMPLEX array, dimension (N)\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*\n\n*  Further Details\n*  ===============\n*\n*  The matrix Q is represented as a product of (ihi-ilo) elementary\n*  reflectors\n*\n*     Q = H(ilo) H(ilo+1) . . . H(ihi-1).\n*\n*  Each H(i) has the form\n*\n*     H(i) = I - tau * v * v'\n*\n*  where tau is a complex scalar, and v is a complex vector with\n*  v(1:i) = 0, v(i+1) = 1 and v(ihi+1:n) = 0; v(i+2:ihi) is stored on\n*  exit in A(i+2:ihi,i), and tau in TAU(i).\n*\n*  The contents of A are illustrated by the following example, with\n*  n = 7, ilo = 2 and ihi = 6:\n*\n*  on entry,                        on exit,\n*\n*  ( a   a   a   a   a   a   a )    (  a   a   h   h   h   h   a )\n*  (     a   a   a   a   a   a )    (      a   h   h   h   h   a )\n*  (     a   a   a   a   a   a )    (      h   h   h   h   h   h )\n*  (     a   a   a   a   a   a )    (      v2  h   h   h   h   h )\n*  (     a   a   a   a   a   a )    (      v2  v3  h   h   h   h )\n*  (     a   a   a   a   a   a )    (      v2  v3  v4  h   h   h )\n*  (                         a )    (                          a )\n*\n*  where a denotes an element of the original matrix A, h denotes a\n*  modified element of the upper Hessenberg matrix H, and vi denotes an\n*  element of the vector defining H(i).\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_ilo = argv[0];
  rb_ihi = argv[1];
  rb_a = argv[2];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SCOMPLEX)
    rb_a = na_change_type(rb_a, NA_SCOMPLEX);
  a = NA_PTR_TYPE(rb_a, complex*);
  ilo = NUM2INT(rb_ilo);
  ihi = NUM2INT(rb_ihi);
  {
    int shape[1];
    shape[0] = n-1;
    rb_tau = na_make_object(NA_SCOMPLEX, 1, shape, cNArray);
  }
  tau = NA_PTR_TYPE(rb_tau, complex*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rb_a_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rb_a_out__, complex*);
  MEMCPY(a_out__, a, complex, NA_TOTAL(rb_a));
  rb_a = rb_a_out__;
  a = a_out__;
  work = ALLOC_N(complex, (n));

  cgehd2_(&n, &ilo, &ihi, a, &lda, tau, work, &info);

  free(work);
  rb_info = INT2NUM(info);
  return rb_ary_new3(3, rb_tau, rb_info, rb_a);
}

void
init_lapack_cgehd2(VALUE mLapack){
  rb_define_module_function(mLapack, "cgehd2", rb_cgehd2, -1);
}
