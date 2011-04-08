#include "rb_lapack.h"

extern doublereal dla_syrcond_(char *uplo, integer *n, doublereal *a, integer *lda, doublereal *af, integer *ldaf, integer *ipiv, integer *cmode, doublereal *c, integer *info, doublereal *work, integer *iwork);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dla_syrcond(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_uplo;
  char uplo; 
  VALUE rblapack_a;
  doublereal *a; 
  VALUE rblapack_af;
  doublereal *af; 
  VALUE rblapack_ipiv;
  integer *ipiv; 
  VALUE rblapack_cmode;
  integer cmode; 
  VALUE rblapack_c;
  doublereal *c; 
  VALUE rblapack_work;
  doublereal *work; 
  VALUE rblapack_iwork;
  integer *iwork; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack___out__;
  doublereal __out__; 

  integer lda;
  integer n;
  integer ldaf;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, __out__ = NumRu::Lapack.dla_syrcond( uplo, a, af, ipiv, cmode, c, work, iwork, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      DOUBLE PRECISION FUNCTION DLA_SYRCOND( UPLO, N, A, LDA, AF, LDAF,  IPIV, CMODE, C, INFO, WORK, IWORK )\n\n*  Purpose\n*  =======\n*\n*     DLA_SYRCOND estimates the Skeel condition number of  op(A) * op2(C)\n*     where op2 is determined by CMODE as follows\n*     CMODE =  1    op2(C) = C\n*     CMODE =  0    op2(C) = I\n*     CMODE = -1    op2(C) = inv(C)\n*     The Skeel condition number cond(A) = norminf( |inv(A)||A| )\n*     is computed by computing scaling factors R such that\n*     diag(R)*A*op2(C) is row equilibrated and computing the standard\n*     infinity-norm condition number.\n*\n\n*  Arguments\n*  ==========\n*\n*     UPLO    (input) CHARACTER*1\n*       = 'U':  Upper triangle of A is stored;\n*       = 'L':  Lower triangle of A is stored.\n*\n*     N       (input) INTEGER\n*     The number of linear equations, i.e., the order of the\n*     matrix A.  N >= 0.\n*\n*     A       (input) DOUBLE PRECISION array, dimension (LDA,N)\n*     On entry, the N-by-N matrix A.\n*\n*     LDA     (input) INTEGER\n*     The leading dimension of the array A.  LDA >= max(1,N).\n*\n*     AF      (input) DOUBLE PRECISION array, dimension (LDAF,N)\n*     The block diagonal matrix D and the multipliers used to\n*     obtain the factor U or L as computed by DSYTRF.\n*\n*     LDAF    (input) INTEGER\n*     The leading dimension of the array AF.  LDAF >= max(1,N).\n*\n*     IPIV    (input) INTEGER array, dimension (N)\n*     Details of the interchanges and the block structure of D\n*     as determined by DSYTRF.\n*\n*     CMODE   (input) INTEGER\n*     Determines op2(C) in the formula op(A) * op2(C) as follows:\n*     CMODE =  1    op2(C) = C\n*     CMODE =  0    op2(C) = I\n*     CMODE = -1    op2(C) = inv(C)\n*\n*     C       (input) DOUBLE PRECISION array, dimension (N)\n*     The vector C in the formula op(A) * op2(C).\n*\n*     INFO    (output) INTEGER\n*       = 0:  Successful exit.\n*     i > 0:  The ith argument is invalid.\n*\n*     WORK    (input) DOUBLE PRECISION array, dimension (3*N).\n*     Workspace.\n*\n*     IWORK   (input) INTEGER array, dimension (N).\n*     Workspace.\n*\n\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      CHARACTER          NORMIN\n      INTEGER            KASE, I, J\n      DOUBLE PRECISION   AINVNM, SMLNUM, TMP\n      LOGICAL            UP\n*     ..\n*     .. Local Arrays ..\n      INTEGER            ISAVE( 3 )\n*     ..\n*     .. External Functions ..\n      LOGICAL            LSAME\n      INTEGER            IDAMAX\n      DOUBLE PRECISION   DLAMCH\n      EXTERNAL           LSAME, IDAMAX, DLAMCH\n*     ..\n*     .. External Subroutines ..\n      EXTERNAL           DLACN2, DLATRS, DRSCL, XERBLA, DSYTRS\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          ABS, MAX\n*     ..\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  info, __out__ = NumRu::Lapack.dla_syrcond( uplo, a, af, ipiv, cmode, c, work, iwork, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 8)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 8)", argc);
  rblapack_uplo = argv[0];
  rblapack_a = argv[1];
  rblapack_af = argv[2];
  rblapack_ipiv = argv[3];
  rblapack_cmode = argv[4];
  rblapack_c = argv[5];
  rblapack_work = argv[6];
  rblapack_iwork = argv[7];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_ipiv))
    rb_raise(rb_eArgError, "ipiv (4th argument) must be NArray");
  if (NA_RANK(rblapack_ipiv) != 1)
    rb_raise(rb_eArgError, "rank of ipiv (4th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_ipiv);
  if (NA_TYPE(rblapack_ipiv) != NA_LINT)
    rblapack_ipiv = na_change_type(rblapack_ipiv, NA_LINT);
  ipiv = NA_PTR_TYPE(rblapack_ipiv, integer*);
  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (2th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 2)
    rb_raise(rb_eArgError, "rank of a (2th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_a) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of a must be the same as shape 0 of ipiv");
  lda = NA_SHAPE0(rblapack_a);
  if (NA_TYPE(rblapack_a) != NA_DFLOAT)
    rblapack_a = na_change_type(rblapack_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rblapack_a, doublereal*);
  if (!NA_IsNArray(rblapack_c))
    rb_raise(rb_eArgError, "c (6th argument) must be NArray");
  if (NA_RANK(rblapack_c) != 1)
    rb_raise(rb_eArgError, "rank of c (6th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_c) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of c must be the same as shape 0 of ipiv");
  if (NA_TYPE(rblapack_c) != NA_DFLOAT)
    rblapack_c = na_change_type(rblapack_c, NA_DFLOAT);
  c = NA_PTR_TYPE(rblapack_c, doublereal*);
  if (!NA_IsNArray(rblapack_af))
    rb_raise(rb_eArgError, "af (3th argument) must be NArray");
  if (NA_RANK(rblapack_af) != 2)
    rb_raise(rb_eArgError, "rank of af (3th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_af) != n)
    rb_raise(rb_eRuntimeError, "shape 1 of af must be the same as shape 0 of ipiv");
  ldaf = NA_SHAPE0(rblapack_af);
  if (NA_TYPE(rblapack_af) != NA_DFLOAT)
    rblapack_af = na_change_type(rblapack_af, NA_DFLOAT);
  af = NA_PTR_TYPE(rblapack_af, doublereal*);
  uplo = StringValueCStr(rblapack_uplo)[0];
  cmode = NUM2INT(rblapack_cmode);
  if (!NA_IsNArray(rblapack_iwork))
    rb_raise(rb_eArgError, "iwork (8th argument) must be NArray");
  if (NA_RANK(rblapack_iwork) != 1)
    rb_raise(rb_eArgError, "rank of iwork (8th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_iwork) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of iwork must be the same as shape 0 of ipiv");
  if (NA_TYPE(rblapack_iwork) != NA_LINT)
    rblapack_iwork = na_change_type(rblapack_iwork, NA_LINT);
  iwork = NA_PTR_TYPE(rblapack_iwork, integer*);
  if (!NA_IsNArray(rblapack_work))
    rb_raise(rb_eArgError, "work (7th argument) must be NArray");
  if (NA_RANK(rblapack_work) != 1)
    rb_raise(rb_eArgError, "rank of work (7th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_work) != (3*n))
    rb_raise(rb_eRuntimeError, "shape 0 of work must be %d", 3*n);
  if (NA_TYPE(rblapack_work) != NA_DFLOAT)
    rblapack_work = na_change_type(rblapack_work, NA_DFLOAT);
  work = NA_PTR_TYPE(rblapack_work, doublereal*);

  __out__ = dla_syrcond_(&uplo, &n, a, &lda, af, &ldaf, ipiv, &cmode, c, &info, work, iwork);

  rblapack_info = INT2NUM(info);
  rblapack___out__ = rb_float_new((double)__out__);
  return rb_ary_new3(2, rblapack_info, rblapack___out__);
}

void
init_lapack_dla_syrcond(VALUE mLapack){
  rb_define_module_function(mLapack, "dla_syrcond", rblapack_dla_syrcond, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
