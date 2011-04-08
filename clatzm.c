#include "rb_lapack.h"

extern VOID clatzm_(char *side, integer *m, integer *n, complex *v, integer *incv, complex *tau, complex *c1, complex *c2, integer *ldc, complex *work);

static VALUE sHelp, sUsage;

static VALUE
rblapack_clatzm(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_side;
  char side; 
  VALUE rblapack_m;
  integer m; 
  VALUE rblapack_n;
  integer n; 
  VALUE rblapack_v;
  complex *v; 
  VALUE rblapack_incv;
  integer incv; 
  VALUE rblapack_tau;
  complex tau; 
  VALUE rblapack_c1;
  complex *c1; 
  VALUE rblapack_c2;
  complex *c2; 
  VALUE rblapack_c1_out__;
  complex *c1_out__;
  VALUE rblapack_c2_out__;
  complex *c2_out__;
  complex *work;

  integer ldc;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  c1, c2 = NumRu::Lapack.clatzm( side, m, n, v, incv, tau, c1, c2, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE CLATZM( SIDE, M, N, V, INCV, TAU, C1, C2, LDC, WORK )\n\n*  Purpose\n*  =======\n*\n*  This routine is deprecated and has been replaced by routine CUNMRZ.\n*\n*  CLATZM applies a Householder matrix generated by CTZRQF to a matrix.\n*\n*  Let P = I - tau*u*u',   u = ( 1 ),\n*                              ( v )\n*  where v is an (m-1) vector if SIDE = 'L', or a (n-1) vector if\n*  SIDE = 'R'.\n*\n*  If SIDE equals 'L', let\n*         C = [ C1 ] 1\n*             [ C2 ] m-1\n*               n\n*  Then C is overwritten by P*C.\n*\n*  If SIDE equals 'R', let\n*         C = [ C1, C2 ] m\n*                1  n-1\n*  Then C is overwritten by C*P.\n*\n\n*  Arguments\n*  =========\n*\n*  SIDE    (input) CHARACTER*1\n*          = 'L': form P * C\n*          = 'R': form C * P\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix C.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix C.\n*\n*  V       (input) COMPLEX array, dimension\n*                  (1 + (M-1)*abs(INCV)) if SIDE = 'L'\n*                  (1 + (N-1)*abs(INCV)) if SIDE = 'R'\n*          The vector v in the representation of P. V is not used\n*          if TAU = 0.\n*\n*  INCV    (input) INTEGER\n*          The increment between elements of v. INCV <> 0\n*\n*  TAU     (input) COMPLEX\n*          The value tau in the representation of P.\n*\n*  C1      (input/output) COMPLEX array, dimension\n*                         (LDC,N) if SIDE = 'L'\n*                         (M,1)   if SIDE = 'R'\n*          On entry, the n-vector C1 if SIDE = 'L', or the m-vector C1\n*          if SIDE = 'R'.\n*\n*          On exit, the first row of P*C if SIDE = 'L', or the first\n*          column of C*P if SIDE = 'R'.\n*\n*  C2      (input/output) COMPLEX array, dimension\n*                         (LDC, N)   if SIDE = 'L'\n*                         (LDC, N-1) if SIDE = 'R'\n*          On entry, the (m - 1) x n matrix C2 if SIDE = 'L', or the\n*          m x (n - 1) matrix C2 if SIDE = 'R'.\n*\n*          On exit, rows 2:m of P*C if SIDE = 'L', or columns 2:m of C*P\n*          if SIDE = 'R'.\n*\n*  LDC     (input) INTEGER\n*          The leading dimension of the arrays C1 and C2.\n*          LDC >= max(1,M).\n*\n*  WORK    (workspace) COMPLEX array, dimension\n*                      (N) if SIDE = 'L'\n*                      (M) if SIDE = 'R'\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  c1, c2 = NumRu::Lapack.clatzm( side, m, n, v, incv, tau, c1, c2, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 8)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 8)", argc);
  rblapack_side = argv[0];
  rblapack_m = argv[1];
  rblapack_n = argv[2];
  rblapack_v = argv[3];
  rblapack_incv = argv[4];
  rblapack_tau = argv[5];
  rblapack_c1 = argv[6];
  rblapack_c2 = argv[7];
  if (rb_options != Qnil) {
  }

  tau.r = (real)NUM2DBL(rb_funcall(rblapack_tau, rb_intern("real"), 0));
  tau.i = (real)NUM2DBL(rb_funcall(rblapack_tau, rb_intern("imag"), 0));
  side = StringValueCStr(rblapack_side)[0];
  m = NUM2INT(rblapack_m);
  incv = NUM2INT(rblapack_incv);
  n = NUM2INT(rblapack_n);
  if (!NA_IsNArray(rblapack_c2))
    rb_raise(rb_eArgError, "c2 (8th argument) must be NArray");
  if (NA_RANK(rblapack_c2) != 2)
    rb_raise(rb_eArgError, "rank of c2 (8th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_c2) != (lsame_(&side,"L") ? n : lsame_(&side,"R") ? n-1 : 0))
    rb_raise(rb_eRuntimeError, "shape 1 of c2 must be %d", lsame_(&side,"L") ? n : lsame_(&side,"R") ? n-1 : 0);
  ldc = NA_SHAPE0(rblapack_c2);
  if (NA_TYPE(rblapack_c2) != NA_SCOMPLEX)
    rblapack_c2 = na_change_type(rblapack_c2, NA_SCOMPLEX);
  c2 = NA_PTR_TYPE(rblapack_c2, complex*);
  if (!NA_IsNArray(rblapack_c1))
    rb_raise(rb_eArgError, "c1 (7th argument) must be NArray");
  if (NA_RANK(rblapack_c1) != 2)
    rb_raise(rb_eArgError, "rank of c1 (7th argument) must be %d", 2);
  if (NA_SHAPE1(rblapack_c1) != (lsame_(&side,"L") ? n : lsame_(&side,"R") ? 1 : 0))
    rb_raise(rb_eRuntimeError, "shape 1 of c1 must be %d", lsame_(&side,"L") ? n : lsame_(&side,"R") ? 1 : 0);
  if (NA_SHAPE0(rblapack_c1) != (lsame_(&side,"L") ? ldc : lsame_(&side,"R") ? m : 0))
    rb_raise(rb_eRuntimeError, "shape 0 of c1 must be %d", lsame_(&side,"L") ? ldc : lsame_(&side,"R") ? m : 0);
  if (NA_TYPE(rblapack_c1) != NA_SCOMPLEX)
    rblapack_c1 = na_change_type(rblapack_c1, NA_SCOMPLEX);
  c1 = NA_PTR_TYPE(rblapack_c1, complex*);
  if (!NA_IsNArray(rblapack_v))
    rb_raise(rb_eArgError, "v (4th argument) must be NArray");
  if (NA_RANK(rblapack_v) != 1)
    rb_raise(rb_eArgError, "rank of v (4th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_v) != (1 + (m-1)*abs(incv)))
    rb_raise(rb_eRuntimeError, "shape 0 of v must be %d", 1 + (m-1)*abs(incv));
  if (NA_TYPE(rblapack_v) != NA_SCOMPLEX)
    rblapack_v = na_change_type(rblapack_v, NA_SCOMPLEX);
  v = NA_PTR_TYPE(rblapack_v, complex*);
  {
    int shape[2];
    shape[0] = lsame_(&side,"L") ? ldc : lsame_(&side,"R") ? m : 0;
    shape[1] = lsame_(&side,"L") ? n : lsame_(&side,"R") ? 1 : 0;
    rblapack_c1_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  c1_out__ = NA_PTR_TYPE(rblapack_c1_out__, complex*);
  MEMCPY(c1_out__, c1, complex, NA_TOTAL(rblapack_c1));
  rblapack_c1 = rblapack_c1_out__;
  c1 = c1_out__;
  {
    int shape[2];
    shape[0] = ldc;
    shape[1] = lsame_(&side,"L") ? n : lsame_(&side,"R") ? n-1 : 0;
    rblapack_c2_out__ = na_make_object(NA_SCOMPLEX, 2, shape, cNArray);
  }
  c2_out__ = NA_PTR_TYPE(rblapack_c2_out__, complex*);
  MEMCPY(c2_out__, c2, complex, NA_TOTAL(rblapack_c2));
  rblapack_c2 = rblapack_c2_out__;
  c2 = c2_out__;
  work = ALLOC_N(complex, (lsame_(&side,"L") ? n : lsame_(&side,"R") ? m : 0));

  clatzm_(&side, &m, &n, v, &incv, &tau, c1, c2, &ldc, work);

  free(work);
  return rb_ary_new3(2, rblapack_c1, rblapack_c2);
}

void
init_lapack_clatzm(VALUE mLapack){
  rb_define_module_function(mLapack, "clatzm", rblapack_clatzm, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
