#include "rb_lapack.h"

static VALUE
rb_dlasy2(int argc, VALUE *argv, VALUE self){
  VALUE rb_ltranl;
  logical ltranl; 
  VALUE rb_ltranr;
  logical ltranr; 
  VALUE rb_isgn;
  integer isgn; 
  VALUE rb_n1;
  integer n1; 
  VALUE rb_n2;
  integer n2; 
  VALUE rb_tl;
  doublereal *tl; 
  VALUE rb_tr;
  doublereal *tr; 
  VALUE rb_b;
  doublereal *b; 
  VALUE rb_scale;
  doublereal scale; 
  VALUE rb_x;
  doublereal *x; 
  VALUE rb_xnorm;
  doublereal xnorm; 
  VALUE rb_info;
  integer info; 

  integer ldtl;
  integer ldtr;
  integer ldb;
  integer ldx;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  scale, x, xnorm, info = NumRu::Lapack.dlasy2( ltranl, ltranr, isgn, n1, n2, tl, tr, b)\n    or\n  NumRu::Lapack.dlasy2  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLASY2( LTRANL, LTRANR, ISGN, N1, N2, TL, LDTL, TR, LDTR, B, LDB, SCALE, X, LDX, XNORM, INFO )\n\n*  Purpose\n*  =======\n*\n*  DLASY2 solves for the N1 by N2 matrix X, 1 <= N1,N2 <= 2, in\n*\n*         op(TL)*X + ISGN*X*op(TR) = SCALE*B,\n*\n*  where TL is N1 by N1, TR is N2 by N2, B is N1 by N2, and ISGN = 1 or\n*  -1.  op(T) = T or T', where T' denotes the transpose of T.\n*\n\n*  Arguments\n*  =========\n*\n*  LTRANL  (input) LOGICAL\n*          On entry, LTRANL specifies the op(TL):\n*             = .FALSE., op(TL) = TL,\n*             = .TRUE., op(TL) = TL'.\n*\n*  LTRANR  (input) LOGICAL\n*          On entry, LTRANR specifies the op(TR):\n*            = .FALSE., op(TR) = TR,\n*            = .TRUE., op(TR) = TR'.\n*\n*  ISGN    (input) INTEGER\n*          On entry, ISGN specifies the sign of the equation\n*          as described before. ISGN may only be 1 or -1.\n*\n*  N1      (input) INTEGER\n*          On entry, N1 specifies the order of matrix TL.\n*          N1 may only be 0, 1 or 2.\n*\n*  N2      (input) INTEGER\n*          On entry, N2 specifies the order of matrix TR.\n*          N2 may only be 0, 1 or 2.\n*\n*  TL      (input) DOUBLE PRECISION array, dimension (LDTL,2)\n*          On entry, TL contains an N1 by N1 matrix.\n*\n*  LDTL    (input) INTEGER\n*          The leading dimension of the matrix TL. LDTL >= max(1,N1).\n*\n*  TR      (input) DOUBLE PRECISION array, dimension (LDTR,2)\n*          On entry, TR contains an N2 by N2 matrix.\n*\n*  LDTR    (input) INTEGER\n*          The leading dimension of the matrix TR. LDTR >= max(1,N2).\n*\n*  B       (input) DOUBLE PRECISION array, dimension (LDB,2)\n*          On entry, the N1 by N2 matrix B contains the right-hand\n*          side of the equation.\n*\n*  LDB     (input) INTEGER\n*          The leading dimension of the matrix B. LDB >= max(1,N1).\n*\n*  SCALE   (output) DOUBLE PRECISION\n*          On exit, SCALE contains the scale factor. SCALE is chosen\n*          less than or equal to 1 to prevent the solution overflowing.\n*\n*  X       (output) DOUBLE PRECISION array, dimension (LDX,2)\n*          On exit, X contains the N1 by N2 solution.\n*\n*  LDX     (input) INTEGER\n*          The leading dimension of the matrix X. LDX >= max(1,N1).\n*\n*  XNORM   (output) DOUBLE PRECISION\n*          On exit, XNORM is the infinity-norm of the solution.\n*\n*  INFO    (output) INTEGER\n*          On exit, INFO is set to\n*             0: successful exit.\n*             1: TL and TR have too close eigenvalues, so TL or\n*                TR is perturbed to get a nonsingular equation.\n*          NOTE: In the interests of speed, this routine does not\n*                check the inputs for errors.\n*\n\n* =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 8)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 8)", argc);
  rb_ltranl = argv[0];
  rb_ltranr = argv[1];
  rb_isgn = argv[2];
  rb_n1 = argv[3];
  rb_n2 = argv[4];
  rb_tl = argv[5];
  rb_tr = argv[6];
  rb_b = argv[7];

  ltranl = (rb_ltranl == Qtrue);
  ltranr = (rb_ltranr == Qtrue);
  isgn = NUM2INT(rb_isgn);
  n1 = NUM2INT(rb_n1);
  n2 = NUM2INT(rb_n2);
  if (!NA_IsNArray(rb_tl))
    rb_raise(rb_eArgError, "tl (6th argument) must be NArray");
  if (NA_RANK(rb_tl) != 2)
    rb_raise(rb_eArgError, "rank of tl (6th argument) must be %d", 2);
  ldtl = NA_SHAPE0(rb_tl);
  if (NA_SHAPE1(rb_tl) != (2))
    rb_raise(rb_eRuntimeError, "shape 1 of tl must be %d", 2);
  if (NA_TYPE(rb_tl) != NA_DFLOAT)
    rb_tl = na_change_type(rb_tl, NA_DFLOAT);
  tl = NA_PTR_TYPE(rb_tl, doublereal*);
  if (!NA_IsNArray(rb_tr))
    rb_raise(rb_eArgError, "tr (7th argument) must be NArray");
  if (NA_RANK(rb_tr) != 2)
    rb_raise(rb_eArgError, "rank of tr (7th argument) must be %d", 2);
  ldtr = NA_SHAPE0(rb_tr);
  if (NA_SHAPE1(rb_tr) != (2))
    rb_raise(rb_eRuntimeError, "shape 1 of tr must be %d", 2);
  if (NA_TYPE(rb_tr) != NA_DFLOAT)
    rb_tr = na_change_type(rb_tr, NA_DFLOAT);
  tr = NA_PTR_TYPE(rb_tr, doublereal*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (8th argument) must be NArray");
  if (NA_RANK(rb_b) != 2)
    rb_raise(rb_eArgError, "rank of b (8th argument) must be %d", 2);
  ldb = NA_SHAPE0(rb_b);
  if (NA_SHAPE1(rb_b) != (2))
    rb_raise(rb_eRuntimeError, "shape 1 of b must be %d", 2);
  if (NA_TYPE(rb_b) != NA_DFLOAT)
    rb_b = na_change_type(rb_b, NA_DFLOAT);
  b = NA_PTR_TYPE(rb_b, doublereal*);
  ldx = MAX(1,n1);
  {
    int shape[2];
    shape[0] = ldx;
    shape[1] = 2;
    rb_x = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  x = NA_PTR_TYPE(rb_x, doublereal*);

  dlasy2_(&ltranl, &ltranr, &isgn, &n1, &n2, tl, &ldtl, tr, &ldtr, b, &ldb, &scale, x, &ldx, &xnorm, &info);

  rb_scale = rb_float_new((double)scale);
  rb_xnorm = rb_float_new((double)xnorm);
  rb_info = INT2NUM(info);
  return rb_ary_new3(4, rb_scale, rb_x, rb_xnorm, rb_info);
}

void
init_lapack_dlasy2(VALUE mLapack){
  rb_define_module_function(mLapack, "dlasy2", rb_dlasy2, -1);
}
