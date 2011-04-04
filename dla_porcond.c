#include "rb_lapack.h"

extern doublereal dla_porcond_(char *uplo, integer *n, doublereal *a, integer *lda, doublereal *af, integer *ldaf, integer *cmode, doublereal *c, integer *info, doublereal *work, integer *iwork);

static VALUE
rb_dla_porcond(int argc, VALUE *argv, VALUE self){
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_a;
  doublereal *a; 
  VALUE rb_af;
  doublereal *af; 
  VALUE rb_cmode;
  integer cmode; 
  VALUE rb_c;
  doublereal *c; 
  VALUE rb_work;
  doublereal *work; 
  VALUE rb_iwork;
  integer *iwork; 
  VALUE rb_info;
  integer info; 
  VALUE rb___out__;
  doublereal __out__; 

  integer lda;
  integer n;
  integer ldaf;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, __out__ = NumRu::Lapack.dla_porcond( uplo, a, af, cmode, c, work, iwork)\n    or\n  NumRu::Lapack.dla_porcond  # print help\n\n\nFORTRAN MANUAL\n      DOUBLE PRECISION FUNCTION DLA_PORCOND( UPLO, N, A, LDA, AF, LDAF, CMODE, C, INFO, WORK, IWORK )\n\n*  Purpose\n*  =======\n*\n*     DLA_PORCOND Estimates the Skeel condition number of  op(A) * op2(C)\n*     where op2 is determined by CMODE as follows\n*     CMODE =  1    op2(C) = C\n*     CMODE =  0    op2(C) = I\n*     CMODE = -1    op2(C) = inv(C)\n*     The Skeel condition number  cond(A) = norminf( |inv(A)||A| )\n*     is computed by computing scaling factors R such that\n*     diag(R)*A*op2(C) is row equilibrated and computing the standard\n*     infinity-norm condition number.\n*\n\n*  Arguments\n*  ==========\n*\n*     UPLO    (input) CHARACTER*1\n*       = 'U':  Upper triangle of A is stored;\n*       = 'L':  Lower triangle of A is stored.\n*\n*     N       (input) INTEGER\n*     The number of linear equations, i.e., the order of the\n*     matrix A.  N >= 0.\n*\n*     A       (input) DOUBLE PRECISION array, dimension (LDA,N)\n*     On entry, the N-by-N matrix A.\n*\n*     LDA     (input) INTEGER\n*     The leading dimension of the array A.  LDA >= max(1,N).\n*\n*     AF      (input) DOUBLE PRECISION array, dimension (LDAF,N)\n*     The triangular factor U or L from the Cholesky factorization\n*     A = U**T*U or A = L*L**T, as computed by DPOTRF.\n*\n*     LDAF    (input) INTEGER\n*     The leading dimension of the array AF.  LDAF >= max(1,N).\n*\n*     CMODE   (input) INTEGER\n*     Determines op2(C) in the formula op(A) * op2(C) as follows:\n*     CMODE =  1    op2(C) = C\n*     CMODE =  0    op2(C) = I\n*     CMODE = -1    op2(C) = inv(C)\n*\n*     C       (input) DOUBLE PRECISION array, dimension (N)\n*     The vector C in the formula op(A) * op2(C).\n*\n*     INFO    (output) INTEGER\n*       = 0:  Successful exit.\n*     i > 0:  The ith argument is invalid.\n*\n*     WORK    (input) DOUBLE PRECISION array, dimension (3*N).\n*     Workspace.\n*\n*     IWORK   (input) INTEGER array, dimension (N).\n*     Workspace.\n*\n\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      INTEGER            KASE, I, J\n      DOUBLE PRECISION   AINVNM, TMP\n      LOGICAL            UP\n*     ..\n*     .. Array Arguments ..\n      INTEGER            ISAVE( 3 )\n*     ..\n*     .. External Functions ..\n      LOGICAL            LSAME\n      INTEGER            IDAMAX\n      EXTERNAL           LSAME, IDAMAX\n*     ..\n*     .. External Subroutines ..\n      EXTERNAL           DLACN2, DPOTRS, XERBLA\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          ABS, MAX\n*     ..\n\n");
    return Qnil;
  }
  if (argc != 7)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 7)", argc);
  rb_uplo = argv[0];
  rb_a = argv[1];
  rb_af = argv[2];
  rb_cmode = argv[3];
  rb_c = argv[4];
  rb_work = argv[5];
  rb_iwork = argv[6];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_DFLOAT)
    rb_a = na_change_type(rb_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rb_a, doublereal*);
  if (!NA_IsNArray(rb_c))
    rb_raise(rb_eArgError, "c (5th argument) must be NArray");
  if (NA_RANK(rb_c) != 1)
    rb_raise(rb_eArgError, "rank of c (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_c) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of c must be the same as shape 1 of a");
  if (NA_TYPE(rb_c) != NA_DFLOAT)
    rb_c = na_change_type(rb_c, NA_DFLOAT);
  c = NA_PTR_TYPE(rb_c, doublereal*);
  if (!NA_IsNArray(rb_af))
    rb_raise(rb_eArgError, "af (3th argument) must be NArray");
  if (NA_RANK(rb_af) != 2)
    rb_raise(rb_eArgError, "rank of af (3th argument) must be %d", 2);
  if (NA_SHAPE1(rb_af) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of af must be the same as shape 1 of a");
  ldaf = NA_SHAPE0(rb_af);
  if (NA_TYPE(rb_af) != NA_DFLOAT)
    rb_af = na_change_type(rb_af, NA_DFLOAT);
  af = NA_PTR_TYPE(rb_af, doublereal*);
  uplo = StringValueCStr(rb_uplo)[0];
  cmode = NUM2INT(rb_cmode);
  if (!NA_IsNArray(rb_iwork))
    rb_raise(rb_eArgError, "iwork (7th argument) must be NArray");
  if (NA_RANK(rb_iwork) != 1)
    rb_raise(rb_eArgError, "rank of iwork (7th argument) must be %d", 1);
  if (NA_SHAPE0(rb_iwork) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of iwork must be the same as shape 1 of a");
  if (NA_TYPE(rb_iwork) != NA_LINT)
    rb_iwork = na_change_type(rb_iwork, NA_LINT);
  iwork = NA_PTR_TYPE(rb_iwork, integer*);
  if (!NA_IsNArray(rb_work))
    rb_raise(rb_eArgError, "work (6th argument) must be NArray");
  if (NA_RANK(rb_work) != 1)
    rb_raise(rb_eArgError, "rank of work (6th argument) must be %d", 1);
  if (NA_SHAPE0(rb_work) != (3*n))
    rb_raise(rb_eRuntimeError, "shape 0 of work must be %d", 3*n);
  if (NA_TYPE(rb_work) != NA_DFLOAT)
    rb_work = na_change_type(rb_work, NA_DFLOAT);
  work = NA_PTR_TYPE(rb_work, doublereal*);

  __out__ = dla_porcond_(&uplo, &n, a, &lda, af, &ldaf, &cmode, c, &info, work, iwork);

  rb_info = INT2NUM(info);
  rb___out__ = rb_float_new((double)__out__);
  return rb_ary_new3(2, rb_info, rb___out__);
}

void
init_lapack_dla_porcond(VALUE mLapack){
  rb_define_module_function(mLapack, "dla_porcond", rb_dla_porcond, -1);
}
