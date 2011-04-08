#include "rb_lapack.h"

extern VOID slaqtr_(logical *ltran, logical *lreal, integer *n, real *t, integer *ldt, real *b, real *w, real *scale, real *x, real *work, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_slaqtr(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_ltran;
  logical ltran; 
  VALUE rblapack_lreal;
  logical lreal; 
  VALUE rblapack_t;
  real *t; 
  VALUE rblapack_b;
  real *b; 
  VALUE rblapack_w;
  real w; 
  VALUE rblapack_x;
  real *x; 
  VALUE rblapack_scale;
  real scale; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_x_out__;
  real *x_out__;
  real *work;

  integer ldt;
  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  scale, info, x = NumRu::Lapack.slaqtr( ltran, lreal, t, b, w, x, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLAQTR( LTRAN, LREAL, N, T, LDT, B, W, SCALE, X, WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  SLAQTR solves the real quasi-triangular system\n*\n*               op(T)*p = scale*c,               if LREAL = .TRUE.\n*\n*  or the complex quasi-triangular systems\n*\n*             op(T + iB)*(p+iq) = scale*(c+id),  if LREAL = .FALSE.\n*\n*  in real arithmetic, where T is upper quasi-triangular.\n*  If LREAL = .FALSE., then the first diagonal block of T must be\n*  1 by 1, B is the specially structured matrix\n*\n*                 B = [ b(1) b(2) ... b(n) ]\n*                     [       w            ]\n*                     [           w        ]\n*                     [              .     ]\n*                     [                 w  ]\n*\n*  op(A) = A or A', A' denotes the conjugate transpose of\n*  matrix A.\n*\n*  On input, X = [ c ].  On output, X = [ p ].\n*                [ d ]                  [ q ]\n*\n*  This subroutine is designed for the condition number estimation\n*  in routine STRSNA.\n*\n\n*  Arguments\n*  =========\n*\n*  LTRAN   (input) LOGICAL\n*          On entry, LTRAN specifies the option of conjugate transpose:\n*             = .FALSE.,    op(T+i*B) = T+i*B,\n*             = .TRUE.,     op(T+i*B) = (T+i*B)'.\n*\n*  LREAL   (input) LOGICAL\n*          On entry, LREAL specifies the input matrix structure:\n*             = .FALSE.,    the input is complex\n*             = .TRUE.,     the input is real\n*\n*  N       (input) INTEGER\n*          On entry, N specifies the order of T+i*B. N >= 0.\n*\n*  T       (input) REAL array, dimension (LDT,N)\n*          On entry, T contains a matrix in Schur canonical form.\n*          If LREAL = .FALSE., then the first diagonal block of T must\n*          be 1 by 1.\n*\n*  LDT     (input) INTEGER\n*          The leading dimension of the matrix T. LDT >= max(1,N).\n*\n*  B       (input) REAL array, dimension (N)\n*          On entry, B contains the elements to form the matrix\n*          B as described above.\n*          If LREAL = .TRUE., B is not referenced.\n*\n*  W       (input) REAL\n*          On entry, W is the diagonal element of the matrix B.\n*          If LREAL = .TRUE., W is not referenced.\n*\n*  SCALE   (output) REAL\n*          On exit, SCALE is the scale factor.\n*\n*  X       (input/output) REAL array, dimension (2*N)\n*          On entry, X contains the right hand side of the system.\n*          On exit, X is overwritten by the solution.\n*\n*  WORK    (workspace) REAL array, dimension (N)\n*\n*  INFO    (output) INTEGER\n*          On exit, INFO is set to\n*             0: successful exit.\n*               1: the some diagonal 1 by 1 block has been perturbed by\n*                  a small number SMIN to keep nonsingularity.\n*               2: the some diagonal 2 by 2 block has been perturbed by\n*                  a small number in SLALN2 to keep nonsingularity.\n*          NOTE: In the interests of speed, this routine does not\n*                check the inputs for errors.\n*\n\n* =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  scale, info, x = NumRu::Lapack.slaqtr( ltran, lreal, t, b, w, x, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rblapack_ltran = argv[0];
  rblapack_lreal = argv[1];
  rblapack_t = argv[2];
  rblapack_b = argv[3];
  rblapack_w = argv[4];
  rblapack_x = argv[5];
  if (rb_options != Qnil) {
  }

  w = (real)NUM2DBL(rblapack_w);
  lreal = (rblapack_lreal == Qtrue);
  if (!NA_IsNArray(rblapack_t))
    rb_raise(rb_eArgError, "t (3th argument) must be NArray");
  if (NA_RANK(rblapack_t) != 2)
    rb_raise(rb_eArgError, "rank of t (3th argument) must be %d", 2);
  n = NA_SHAPE1(rblapack_t);
  ldt = NA_SHAPE0(rblapack_t);
  if (NA_TYPE(rblapack_t) != NA_SFLOAT)
    rblapack_t = na_change_type(rblapack_t, NA_SFLOAT);
  t = NA_PTR_TYPE(rblapack_t, real*);
  ltran = (rblapack_ltran == Qtrue);
  if (!NA_IsNArray(rblapack_b))
    rb_raise(rb_eArgError, "b (4th argument) must be NArray");
  if (NA_RANK(rblapack_b) != 1)
    rb_raise(rb_eArgError, "rank of b (4th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_b) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of b must be the same as shape 1 of t");
  if (NA_TYPE(rblapack_b) != NA_SFLOAT)
    rblapack_b = na_change_type(rblapack_b, NA_SFLOAT);
  b = NA_PTR_TYPE(rblapack_b, real*);
  if (!NA_IsNArray(rblapack_x))
    rb_raise(rb_eArgError, "x (6th argument) must be NArray");
  if (NA_RANK(rblapack_x) != 1)
    rb_raise(rb_eArgError, "rank of x (6th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_x) != (2*n))
    rb_raise(rb_eRuntimeError, "shape 0 of x must be %d", 2*n);
  if (NA_TYPE(rblapack_x) != NA_SFLOAT)
    rblapack_x = na_change_type(rblapack_x, NA_SFLOAT);
  x = NA_PTR_TYPE(rblapack_x, real*);
  {
    int shape[1];
    shape[0] = 2*n;
    rblapack_x_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  x_out__ = NA_PTR_TYPE(rblapack_x_out__, real*);
  MEMCPY(x_out__, x, real, NA_TOTAL(rblapack_x));
  rblapack_x = rblapack_x_out__;
  x = x_out__;
  work = ALLOC_N(real, (n));

  slaqtr_(&ltran, &lreal, &n, t, &ldt, b, &w, &scale, x, work, &info);

  free(work);
  rblapack_scale = rb_float_new((double)scale);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(3, rblapack_scale, rblapack_info, rblapack_x);
}

void
init_lapack_slaqtr(VALUE mLapack){
  rb_define_module_function(mLapack, "slaqtr", rblapack_slaqtr, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
