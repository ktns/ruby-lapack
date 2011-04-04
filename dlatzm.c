#include "rb_lapack.h"

extern VOID dlatzm_(char *side, integer *m, integer *n, doublereal *v, integer *incv, doublereal *tau, doublereal *c1, doublereal *c2, integer *ldc, doublereal *work);

static VALUE
rb_dlatzm(int argc, VALUE *argv, VALUE self){
  VALUE rb_side;
  char side; 
  VALUE rb_m;
  integer m; 
  VALUE rb_n;
  integer n; 
  VALUE rb_v;
  doublereal *v; 
  VALUE rb_incv;
  integer incv; 
  VALUE rb_tau;
  doublereal tau; 
  VALUE rb_c1;
  doublereal *c1; 
  VALUE rb_c2;
  doublereal *c2; 
  VALUE rb_c1_out__;
  doublereal *c1_out__;
  VALUE rb_c2_out__;
  doublereal *c2_out__;
  doublereal *work;

  integer ldc;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  c1, c2 = NumRu::Lapack.dlatzm( side, m, n, v, incv, tau, c1, c2)\n    or\n  NumRu::Lapack.dlatzm  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLATZM( SIDE, M, N, V, INCV, TAU, C1, C2, LDC, WORK )\n\n*  Purpose\n*  =======\n*\n*  This routine is deprecated and has been replaced by routine DORMRZ.\n*\n*  DLATZM applies a Householder matrix generated by DTZRQF to a matrix.\n*\n*  Let P = I - tau*u*u',   u = ( 1 ),\n*                              ( v )\n*  where v is an (m-1) vector if SIDE = 'L', or a (n-1) vector if\n*  SIDE = 'R'.\n*\n*  If SIDE equals 'L', let\n*         C = [ C1 ] 1\n*             [ C2 ] m-1\n*               n\n*  Then C is overwritten by P*C.\n*\n*  If SIDE equals 'R', let\n*         C = [ C1, C2 ] m\n*                1  n-1\n*  Then C is overwritten by C*P.\n*\n\n*  Arguments\n*  =========\n*\n*  SIDE    (input) CHARACTER*1\n*          = 'L': form P * C\n*          = 'R': form C * P\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix C.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix C.\n*\n*  V       (input) DOUBLE PRECISION array, dimension\n*                  (1 + (M-1)*abs(INCV)) if SIDE = 'L'\n*                  (1 + (N-1)*abs(INCV)) if SIDE = 'R'\n*          The vector v in the representation of P. V is not used\n*          if TAU = 0.\n*\n*  INCV    (input) INTEGER\n*          The increment between elements of v. INCV <> 0\n*\n*  TAU     (input) DOUBLE PRECISION\n*          The value tau in the representation of P.\n*\n*  C1      (input/output) DOUBLE PRECISION array, dimension\n*                         (LDC,N) if SIDE = 'L'\n*                         (M,1)   if SIDE = 'R'\n*          On entry, the n-vector C1 if SIDE = 'L', or the m-vector C1\n*          if SIDE = 'R'.\n*\n*          On exit, the first row of P*C if SIDE = 'L', or the first\n*          column of C*P if SIDE = 'R'.\n*\n*  C2      (input/output) DOUBLE PRECISION array, dimension\n*                         (LDC, N)   if SIDE = 'L'\n*                         (LDC, N-1) if SIDE = 'R'\n*          On entry, the (m - 1) x n matrix C2 if SIDE = 'L', or the\n*          m x (n - 1) matrix C2 if SIDE = 'R'.\n*\n*          On exit, rows 2:m of P*C if SIDE = 'L', or columns 2:m of C*P\n*          if SIDE = 'R'.\n*\n*  LDC     (input) INTEGER\n*          The leading dimension of the arrays C1 and C2. LDC >= (1,M).\n*\n*  WORK    (workspace) DOUBLE PRECISION array, dimension\n*                      (N) if SIDE = 'L'\n*                      (M) if SIDE = 'R'\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 8)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 8)", argc);
  rb_side = argv[0];
  rb_m = argv[1];
  rb_n = argv[2];
  rb_v = argv[3];
  rb_incv = argv[4];
  rb_tau = argv[5];
  rb_c1 = argv[6];
  rb_c2 = argv[7];

  tau = NUM2DBL(rb_tau);
  side = StringValueCStr(rb_side)[0];
  m = NUM2INT(rb_m);
  incv = NUM2INT(rb_incv);
  n = NUM2INT(rb_n);
  if (!NA_IsNArray(rb_c2))
    rb_raise(rb_eArgError, "c2 (8th argument) must be NArray");
  if (NA_RANK(rb_c2) != 2)
    rb_raise(rb_eArgError, "rank of c2 (8th argument) must be %d", 2);
  if (NA_SHAPE1(rb_c2) != (lsame_(&side,"L") ? n : lsame_(&side,"R") ? n-1 : 0))
    rb_raise(rb_eRuntimeError, "shape 1 of c2 must be %d", lsame_(&side,"L") ? n : lsame_(&side,"R") ? n-1 : 0);
  ldc = NA_SHAPE0(rb_c2);
  if (NA_TYPE(rb_c2) != NA_DFLOAT)
    rb_c2 = na_change_type(rb_c2, NA_DFLOAT);
  c2 = NA_PTR_TYPE(rb_c2, doublereal*);
  if (!NA_IsNArray(rb_c1))
    rb_raise(rb_eArgError, "c1 (7th argument) must be NArray");
  if (NA_RANK(rb_c1) != 2)
    rb_raise(rb_eArgError, "rank of c1 (7th argument) must be %d", 2);
  if (NA_SHAPE1(rb_c1) != (lsame_(&side,"L") ? n : lsame_(&side,"R") ? 1 : 0))
    rb_raise(rb_eRuntimeError, "shape 1 of c1 must be %d", lsame_(&side,"L") ? n : lsame_(&side,"R") ? 1 : 0);
  if (NA_SHAPE0(rb_c1) != (lsame_(&side,"L") ? ldc : lsame_(&side,"R") ? m : 0))
    rb_raise(rb_eRuntimeError, "shape 0 of c1 must be %d", lsame_(&side,"L") ? ldc : lsame_(&side,"R") ? m : 0);
  if (NA_TYPE(rb_c1) != NA_DFLOAT)
    rb_c1 = na_change_type(rb_c1, NA_DFLOAT);
  c1 = NA_PTR_TYPE(rb_c1, doublereal*);
  if (!NA_IsNArray(rb_v))
    rb_raise(rb_eArgError, "v (4th argument) must be NArray");
  if (NA_RANK(rb_v) != 1)
    rb_raise(rb_eArgError, "rank of v (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_v) != (1 + (m-1)*abs(incv)))
    rb_raise(rb_eRuntimeError, "shape 0 of v must be %d", 1 + (m-1)*abs(incv));
  if (NA_TYPE(rb_v) != NA_DFLOAT)
    rb_v = na_change_type(rb_v, NA_DFLOAT);
  v = NA_PTR_TYPE(rb_v, doublereal*);
  {
    int shape[2];
    shape[0] = lsame_(&side,"L") ? ldc : lsame_(&side,"R") ? m : 0;
    shape[1] = lsame_(&side,"L") ? n : lsame_(&side,"R") ? 1 : 0;
    rb_c1_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  c1_out__ = NA_PTR_TYPE(rb_c1_out__, doublereal*);
  MEMCPY(c1_out__, c1, doublereal, NA_TOTAL(rb_c1));
  rb_c1 = rb_c1_out__;
  c1 = c1_out__;
  {
    int shape[2];
    shape[0] = ldc;
    shape[1] = lsame_(&side,"L") ? n : lsame_(&side,"R") ? n-1 : 0;
    rb_c2_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  c2_out__ = NA_PTR_TYPE(rb_c2_out__, doublereal*);
  MEMCPY(c2_out__, c2, doublereal, NA_TOTAL(rb_c2));
  rb_c2 = rb_c2_out__;
  c2 = c2_out__;
  work = ALLOC_N(doublereal, (lsame_(&side,"L") ? n : lsame_(&side,"R") ? m : 0));

  dlatzm_(&side, &m, &n, v, &incv, &tau, c1, c2, &ldc, work);

  free(work);
  return rb_ary_new3(2, rb_c1, rb_c2);
}

void
init_lapack_dlatzm(VALUE mLapack){
  rb_define_module_function(mLapack, "dlatzm", rb_dlatzm, -1);
}
