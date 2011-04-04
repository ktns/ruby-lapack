#include "rb_lapack.h"

extern VOID slasr_(char *side, char *pivot, char *direct, integer *m, integer *n, real *c, real *s, real *a, integer *lda);

static VALUE
rb_slasr(int argc, VALUE *argv, VALUE self){
  VALUE rb_side;
  char side; 
  VALUE rb_pivot;
  char pivot; 
  VALUE rb_direct;
  char direct; 
  VALUE rb_m;
  integer m; 
  VALUE rb_c;
  real *c; 
  VALUE rb_s;
  real *s; 
  VALUE rb_a;
  real *a; 
  VALUE rb_a_out__;
  real *a_out__;

  integer lda;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  a = NumRu::Lapack.slasr( side, pivot, direct, m, c, s, a)\n    or\n  NumRu::Lapack.slasr  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLASR( SIDE, PIVOT, DIRECT, M, N, C, S, A, LDA )\n\n*  Purpose\n*  =======\n*\n*  SLASR applies a sequence of plane rotations to a real matrix A,\n*  from either the left or the right.\n*  \n*  When SIDE = 'L', the transformation takes the form\n*  \n*     A := P*A\n*  \n*  and when SIDE = 'R', the transformation takes the form\n*  \n*     A := A*P**T\n*  \n*  where P is an orthogonal matrix consisting of a sequence of z plane\n*  rotations, with z = M when SIDE = 'L' and z = N when SIDE = 'R',\n*  and P**T is the transpose of P.\n*  \n*  When DIRECT = 'F' (Forward sequence), then\n*  \n*     P = P(z-1) * ... * P(2) * P(1)\n*  \n*  and when DIRECT = 'B' (Backward sequence), then\n*  \n*     P = P(1) * P(2) * ... * P(z-1)\n*  \n*  where P(k) is a plane rotation matrix defined by the 2-by-2 rotation\n*  \n*     R(k) = (  c(k)  s(k) )\n*          = ( -s(k)  c(k) ).\n*  \n*  When PIVOT = 'V' (Variable pivot), the rotation is performed\n*  for the plane (k,k+1), i.e., P(k) has the form\n*  \n*     P(k) = (  1                                            )\n*            (       ...                                     )\n*            (              1                                )\n*            (                   c(k)  s(k)                  )\n*            (                  -s(k)  c(k)                  )\n*            (                                1              )\n*            (                                     ...       )\n*            (                                            1  )\n*  \n*  where R(k) appears as a rank-2 modification to the identity matrix in\n*  rows and columns k and k+1.\n*  \n*  When PIVOT = 'T' (Top pivot), the rotation is performed for the\n*  plane (1,k+1), so P(k) has the form\n*  \n*     P(k) = (  c(k)                    s(k)                 )\n*            (         1                                     )\n*            (              ...                              )\n*            (                     1                         )\n*            ( -s(k)                    c(k)                 )\n*            (                                 1             )\n*            (                                      ...      )\n*            (                                             1 )\n*  \n*  where R(k) appears in rows and columns 1 and k+1.\n*  \n*  Similarly, when PIVOT = 'B' (Bottom pivot), the rotation is\n*  performed for the plane (k,z), giving P(k) the form\n*  \n*     P(k) = ( 1                                             )\n*            (      ...                                      )\n*            (             1                                 )\n*            (                  c(k)                    s(k) )\n*            (                         1                     )\n*            (                              ...              )\n*            (                                     1         )\n*            (                 -s(k)                    c(k) )\n*  \n*  where R(k) appears in rows and columns k and z.  The rotations are\n*  performed without ever forming P(k) explicitly.\n*\n\n*  Arguments\n*  =========\n*\n*  SIDE    (input) CHARACTER*1\n*          Specifies whether the plane rotation matrix P is applied to\n*          A on the left or the right.\n*          = 'L':  Left, compute A := P*A\n*          = 'R':  Right, compute A:= A*P**T\n*\n*  PIVOT   (input) CHARACTER*1\n*          Specifies the plane for which P(k) is a plane rotation\n*          matrix.\n*          = 'V':  Variable pivot, the plane (k,k+1)\n*          = 'T':  Top pivot, the plane (1,k+1)\n*          = 'B':  Bottom pivot, the plane (k,z)\n*\n*  DIRECT  (input) CHARACTER*1\n*          Specifies whether P is a forward or backward sequence of\n*          plane rotations.\n*          = 'F':  Forward, P = P(z-1)*...*P(2)*P(1)\n*          = 'B':  Backward, P = P(1)*P(2)*...*P(z-1)\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix A.  If m <= 1, an immediate\n*          return is effected.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix A.  If n <= 1, an\n*          immediate return is effected.\n*\n*  C       (input) REAL array, dimension\n*                  (M-1) if SIDE = 'L'\n*                  (N-1) if SIDE = 'R'\n*          The cosines c(k) of the plane rotations.\n*\n*  S       (input) REAL array, dimension\n*                  (M-1) if SIDE = 'L'\n*                  (N-1) if SIDE = 'R'\n*          The sines s(k) of the plane rotations.  The 2-by-2 plane\n*          rotation part of the matrix P(k), R(k), has the form\n*          R(k) = (  c(k)  s(k) )\n*                 ( -s(k)  c(k) ).\n*\n*  A       (input/output) REAL array, dimension (LDA,N)\n*          The M-by-N matrix A.  On exit, A is overwritten by P*A if\n*          SIDE = 'R' or by A*P**T if SIDE = 'L'.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(1,M).\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 7)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 7)", argc);
  rb_side = argv[0];
  rb_pivot = argv[1];
  rb_direct = argv[2];
  rb_m = argv[3];
  rb_c = argv[4];
  rb_s = argv[5];
  rb_a = argv[6];

  direct = StringValueCStr(rb_direct)[0];
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (7th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (7th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_SFLOAT)
    rb_a = na_change_type(rb_a, NA_SFLOAT);
  a = NA_PTR_TYPE(rb_a, real*);
  side = StringValueCStr(rb_side)[0];
  pivot = StringValueCStr(rb_pivot)[0];
  m = NUM2INT(rb_m);
  if (!NA_IsNArray(rb_c))
    rb_raise(rb_eArgError, "c (5th argument) must be NArray");
  if (NA_RANK(rb_c) != 1)
    rb_raise(rb_eArgError, "rank of c (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_c) != (m-1))
    rb_raise(rb_eRuntimeError, "shape 0 of c must be %d", m-1);
  if (NA_TYPE(rb_c) != NA_SFLOAT)
    rb_c = na_change_type(rb_c, NA_SFLOAT);
  c = NA_PTR_TYPE(rb_c, real*);
  if (!NA_IsNArray(rb_s))
    rb_raise(rb_eArgError, "s (6th argument) must be NArray");
  if (NA_RANK(rb_s) != 1)
    rb_raise(rb_eArgError, "rank of s (6th argument) must be %d", 1);
  if (NA_SHAPE0(rb_s) != (m-1))
    rb_raise(rb_eRuntimeError, "shape 0 of s must be %d", m-1);
  if (NA_TYPE(rb_s) != NA_SFLOAT)
    rb_s = na_change_type(rb_s, NA_SFLOAT);
  s = NA_PTR_TYPE(rb_s, real*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rb_a_out__ = na_make_object(NA_SFLOAT, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rb_a_out__, real*);
  MEMCPY(a_out__, a, real, NA_TOTAL(rb_a));
  rb_a = rb_a_out__;
  a = a_out__;

  slasr_(&side, &pivot, &direct, &m, &n, c, s, a, &lda);

  return rb_a;
}

void
init_lapack_slasr(VALUE mLapack){
  rb_define_module_function(mLapack, "slasr", rb_slasr, -1);
}
