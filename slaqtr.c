#include "rb_lapack.h"

static VALUE
rb_slaqtr(int argc, VALUE *argv, VALUE self){
  VALUE rb_ltran;
  logical ltran; 
  VALUE rb_lreal;
  logical lreal; 
  VALUE rb_t;
  real *t; 
  VALUE rb_b;
  real *b; 
  VALUE rb_w;
  real w; 
  VALUE rb_x;
  real *x; 
  VALUE rb_scale;
  real scale; 
  VALUE rb_info;
  integer info; 
  VALUE rb_x_out__;
  real *x_out__;
  real *work;

  integer ldt;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  scale, info, x = NumRu::Lapack.slaqtr( ltran, lreal, t, b, w, x)\n    or\n  NumRu::Lapack.slaqtr  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLAQTR( LTRAN, LREAL, N, T, LDT, B, W, SCALE, X, WORK, INFO )\n\n*  Purpose\n*  =======\n*\n*  SLAQTR solves the real quasi-triangular system\n*\n*               op(T)*p = scale*c,               if LREAL = .TRUE.\n*\n*  or the complex quasi-triangular systems\n*\n*             op(T + iB)*(p+iq) = scale*(c+id),  if LREAL = .FALSE.\n*\n*  in real arithmetic, where T is upper quasi-triangular.\n*  If LREAL = .FALSE., then the first diagonal block of T must be\n*  1 by 1, B is the specially structured matrix\n*\n*                 B = [ b(1) b(2) ... b(n) ]\n*                     [       w            ]\n*                     [           w        ]\n*                     [              .     ]\n*                     [                 w  ]\n*\n*  op(A) = A or A', A' denotes the conjugate transpose of\n*  matrix A.\n*\n*  On input, X = [ c ].  On output, X = [ p ].\n*                [ d ]                  [ q ]\n*\n*  This subroutine is designed for the condition number estimation\n*  in routine STRSNA.\n*\n\n*  Arguments\n*  =========\n*\n*  LTRAN   (input) LOGICAL\n*          On entry, LTRAN specifies the option of conjugate transpose:\n*             = .FALSE.,    op(T+i*B) = T+i*B,\n*             = .TRUE.,     op(T+i*B) = (T+i*B)'.\n*\n*  LREAL   (input) LOGICAL\n*          On entry, LREAL specifies the input matrix structure:\n*             = .FALSE.,    the input is complex\n*             = .TRUE.,     the input is real\n*\n*  N       (input) INTEGER\n*          On entry, N specifies the order of T+i*B. N >= 0.\n*\n*  T       (input) REAL array, dimension (LDT,N)\n*          On entry, T contains a matrix in Schur canonical form.\n*          If LREAL = .FALSE., then the first diagonal block of T must\n*          be 1 by 1.\n*\n*  LDT     (input) INTEGER\n*          The leading dimension of the matrix T. LDT >= max(1,N).\n*\n*  B       (input) REAL array, dimension (N)\n*          On entry, B contains the elements to form the matrix\n*          B as described above.\n*          If LREAL = .TRUE., B is not referenced.\n*\n*  W       (input) REAL\n*          On entry, W is the diagonal element of the matrix B.\n*          If LREAL = .TRUE., W is not referenced.\n*\n*  SCALE   (output) REAL\n*          On exit, SCALE is the scale factor.\n*\n*  X       (input/output) REAL array, dimension (2*N)\n*          On entry, X contains the right hand side of the system.\n*          On exit, X is overwritten by the solution.\n*\n*  WORK    (workspace) REAL array, dimension (N)\n*\n*  INFO    (output) INTEGER\n*          On exit, INFO is set to\n*             0: successful exit.\n*               1: the some diagonal 1 by 1 block has been perturbed by\n*                  a small number SMIN to keep nonsingularity.\n*               2: the some diagonal 2 by 2 block has been perturbed by\n*                  a small number in SLALN2 to keep nonsingularity.\n*          NOTE: In the interests of speed, this routine does not\n*                check the inputs for errors.\n*\n\n* =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rb_ltran = argv[0];
  rb_lreal = argv[1];
  rb_t = argv[2];
  rb_b = argv[3];
  rb_w = argv[4];
  rb_x = argv[5];

  ltran = (rb_ltran == Qtrue);
  lreal = (rb_lreal == Qtrue);
  w = (real)NUM2DBL(rb_w);
  if (!NA_IsNArray(rb_t))
    rb_raise(rb_eArgError, "t (3th argument) must be NArray");
  if (NA_RANK(rb_t) != 2)
    rb_raise(rb_eArgError, "rank of t (3th argument) must be %d", 2);
  ldt = NA_SHAPE0(rb_t);
  n = NA_SHAPE1(rb_t);
  if (NA_TYPE(rb_t) != NA_SFLOAT)
    rb_t = na_change_type(rb_t, NA_SFLOAT);
  t = NA_PTR_TYPE(rb_t, real*);
  if (!NA_IsNArray(rb_b))
    rb_raise(rb_eArgError, "b (4th argument) must be NArray");
  if (NA_RANK(rb_b) != 1)
    rb_raise(rb_eArgError, "rank of b (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_b) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of b must be the same as shape 1 of t");
  if (NA_TYPE(rb_b) != NA_SFLOAT)
    rb_b = na_change_type(rb_b, NA_SFLOAT);
  b = NA_PTR_TYPE(rb_b, real*);
  if (!NA_IsNArray(rb_x))
    rb_raise(rb_eArgError, "x (6th argument) must be NArray");
  if (NA_RANK(rb_x) != 1)
    rb_raise(rb_eArgError, "rank of x (6th argument) must be %d", 1);
  if (NA_SHAPE0(rb_x) != (2*n))
    rb_raise(rb_eRuntimeError, "shape 0 of x must be %d", 2*n);
  if (NA_TYPE(rb_x) != NA_SFLOAT)
    rb_x = na_change_type(rb_x, NA_SFLOAT);
  x = NA_PTR_TYPE(rb_x, real*);
  {
    int shape[1];
    shape[0] = 2*n;
    rb_x_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  x_out__ = NA_PTR_TYPE(rb_x_out__, real*);
  MEMCPY(x_out__, x, real, NA_TOTAL(rb_x));
  rb_x = rb_x_out__;
  x = x_out__;
  work = ALLOC_N(real, (n));

  slaqtr_(&ltran, &lreal, &n, t, &ldt, b, &w, &scale, x, work, &info);

  free(work);
  rb_scale = rb_float_new((double)scale);
  rb_info = INT2NUM(info);
  return rb_ary_new3(3, rb_scale, rb_info, rb_x);
}

void
init_lapack_slaqtr(VALUE mLapack){
  rb_define_module_function(mLapack, "slaqtr", rb_slaqtr, -1);
}
