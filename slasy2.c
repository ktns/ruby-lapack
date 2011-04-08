#include "rb_lapack.h"

extern VOID slasy2_(logical *ltranl, logical *ltranr, integer *isgn, integer *n1, integer *n2, real *tl, integer *ldtl, real *tr, integer *ldtr, real *b, integer *ldb, real *scale, real *x, integer *ldx, real *xnorm, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_slasy2(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_ltranl;
  logical ltranl; 
  VALUE rblapack_ltranr;
  logical ltranr; 
  VALUE rblapack_isgn;
  integer isgn; 
  VALUE rblapack_n1;
  integer n1; 
  VALUE rblapack_n2;
  integer n2; 
  VALUE rblapack_tl;
  real *tl; 
  VALUE rblapack_tr;
  real *tr; 
  VALUE rblapack_b;
  real *b; 
  VALUE rblapack_scale;
  real scale; 
  VALUE rblapack_x;
  real *x; 
  VALUE rblapack_xnorm;
  real xnorm; 
  VALUE rblapack_info;
  integer info; 

  integer ldtl;
  integer ldtr;
  integer ldb;
  integer ldx;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  scale, x, xnorm, info = NumRu::Lapack.slasy2( ltranl, ltranr, isgn, n1, n2, tl, tr, b, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLASY2( LTRANL, LTRANR, ISGN, N1, N2, TL, LDTL, TR, LDTR, B, LDB, SCALE, X, LDX, XNORM, INFO )\n\n*  Purpose\n*  =======\n*\n*  SLASY2 solves for the N1 by N2 matrix X, 1 <= N1,N2 <= 2, in\n*\n*         op(TL)*X + ISGN*X*op(TR) = SCALE*B,\n*\n*  where TL is N1 by N1, TR is N2 by N2, B is N1 by N2, and ISGN = 1 or\n*  -1.  op(T) = T or T', where T' denotes the transpose of T.\n*\n\n*  Arguments\n*  =========\n*\n*  LTRANL  (input) LOGICAL\n*          On entry, LTRANL specifies the op(TL):\n*             = .FALSE., op(TL) = TL,\n*             = .TRUE., op(TL) = TL'.\n*\n*  LTRANR  (input) LOGICAL\n*          On entry, LTRANR specifies the op(TR):\n*            = .FALSE., op(TR) = TR,\n*            = .TRUE., op(TR) = TR'.\n*\n*  ISGN    (input) INTEGER\n*          On entry, ISGN specifies the sign of the equation\n*          as described before. ISGN may only be 1 or -1.\n*\n*  N1      (input) INTEGER\n*          On entry, N1 specifies the order of matrix TL.\n*          N1 may only be 0, 1 or 2.\n*\n*  N2      (input) INTEGER\n*          On entry, N2 specifies the order of matrix TR.\n*          N2 may only be 0, 1 or 2.\n*\n*  TL      (input) REAL array, dimension (LDTL,2)\n*          On entry, TL contains an N1 by N1 matrix.\n*\n*  LDTL    (input) INTEGER\n*          The leading dimension of the matrix TL. LDTL >= max(1,N1).\n*\n*  TR      (input) REAL array, dimension (LDTR,2)\n*          On entry, TR contains an N2 by N2 matrix.\n*\n*  LDTR    (input) INTEGER\n*          The leading dimension of the matrix TR. LDTR >= max(1,N2).\n*\n*  B       (input) REAL array, dimension (LDB,2)\n*          On entry, the N1 by N2 matrix B contains the right-hand\n*          side of the equation.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the matrix B. LDB >= max(1,N1).\n*\n*  SCALE   (output) REAL\n*          On exit, SCALE contains the scale factor. SCALE is chosen\n*          less than or equal to 1 to prevent the solution overflowing.\n*\n*  X       (output) REAL array, dimension (LDX,2)\n*          On exit, X contains the N1 by N2 solution.\n*\n*  LDX     (input) INTEGER\n*          The leading dimension of the matrix X. LDX >= max(1,N1).\n*\n*  XNORM   (output) REAL\n*          On exit, XNORM is the infinity-norm of the solution.\n*\n*  INFO    (output) INTEGER\n*          On exit, INFO is set to\n*             0: successful exit.\n*             1: TL and TR have too close eigenvalues, so TL or\n*                TR is perturbed to get a nonsingular equation.\n*          NOTE: In the interests of speed, this routine does not\n*                check the inputs for errors.\n*\n\n* =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  scale, x, xnorm, info = NumRu::Lapack.slasy2( ltranl, ltranr, isgn, n1, n2, tl, tr, b, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 8)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 8)", argc);
  rblapack_ltranl = argv[0];
  rblapack_ltranr = argv[1];
  rblapack_isgn = argv[2];
  rblapack_n1 = argv[3];
  rblapack_n2 = argv[4];
  rblapack_tl = argv[5];
  rblapack_tr = argv[6];
  rblapack_b = argv[7];
  if (rb_options != Qnil) {
  }

  ltranl = (rblapack_ltranl == Qtrue);
  if (!NA_IsNArray(rblapack_b))
    rb_raise(rb_eArgError, "b (8th argument) must be NArray");
  if (NA_RANK(rblapack_b) != 2)
    rb_raise(rb_eArgError, "rank of b (8th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_b) != (2))
    rb_raise(rb_eRuntimeError, "shape 1 of b must be %d", 2);
  ldb = NA_SHAPE0(rblapack_b);
  if (NA_TYPE(rblapack_b) != NA_SFLOAT)
    rblapack_b = na_change_type(rblapack_b, NA_SFLOAT);
  b = NA_PTR_TYPE(rblapack_b, real*);
  if (!NA_IsNArray(rblapack_tl))
    rb_raise(rb_eArgError, "tl (6th argument) must be NArray");
  if (NA_RANK(rblapack_tl) != 2)
    rb_raise(rb_eArgError, "rank of tl (6th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_tl) != (2))
    rb_raise(rb_eRuntimeError, "shape 1 of tl must be %d", 2);
  ldtl = NA_SHAPE0(rblapack_tl);
  if (NA_TYPE(rblapack_tl) != NA_SFLOAT)
    rblapack_tl = na_change_type(rblapack_tl, NA_SFLOAT);
  tl = NA_PTR_TYPE(rblapack_tl, real*);
  n1 = NUM2INT(rblapack_n1);
  isgn = NUM2INT(rblapack_isgn);
  ltranr = (rblapack_ltranr == Qtrue);
  if (!NA_IsNArray(rblapack_tr))
    rb_raise(rb_eArgError, "tr (7th argument) must be NArray");
  if (NA_RANK(rblapack_tr) != 2)
    rb_raise(rb_eArgError, "rank of tr (7th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_tr) != (2))
    rb_raise(rb_eRuntimeError, "shape 1 of tr must be %d", 2);
  ldtr = NA_SHAPE0(rblapack_tr);
  if (NA_TYPE(rblapack_tr) != NA_SFLOAT)
    rblapack_tr = na_change_type(rblapack_tr, NA_SFLOAT);
  tr = NA_PTR_TYPE(rblapack_tr, real*);
  n2 = NUM2INT(rblapack_n2);
  ldx = MAX(1,n1);
  {
    int shape[2];
    shape[0] = ldx;
    shape[1] = 2;
    rblapack_x = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  x = NA_PTR_TYPE(rblapack_x, real*);

  slasy2_(&ltranl, &ltranr, &isgn, &n1, &n2, tl, &ldtl, tr, &ldtr, b, &ldb, &scale, x, &ldx, &xnorm, &info);

  rblapack_scale = rb_float_new((double)scale);
  rblapack_xnorm = rb_float_new((double)xnorm);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(4, rblapack_scale, rblapack_x, rblapack_xnorm, rblapack_info);
}

void
init_lapack_slasy2(VALUE mLapack){
  rb_define_module_function(mLapack, "slasy2", rblapack_slasy2, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
