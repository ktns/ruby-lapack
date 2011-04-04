#include "rb_lapack.h"

extern VOID zlaqge_(integer *m, integer *n, doublecomplex *a, integer *lda, doublereal *r, doublereal *c, doublereal *rowcnd, doublereal *colcnd, doublereal *amax, char *equed);

static VALUE
rb_zlaqge(int argc, VALUE *argv, VALUE self){
  VALUE rb_a;
  doublecomplex *a; 
  VALUE rb_r;
  doublereal *r; 
  VALUE rb_c;
  doublereal *c; 
  VALUE rb_rowcnd;
  doublereal rowcnd; 
  VALUE rb_colcnd;
  doublereal colcnd; 
  VALUE rb_amax;
  doublereal amax; 
  VALUE rb_equed;
  char equed; 
  VALUE rb_a_out__;
  doublecomplex *a_out__;

  integer lda;
  integer n;
  integer m;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  equed, a = NumRu::Lapack.zlaqge( a, r, c, rowcnd, colcnd, amax)\n    or\n  NumRu::Lapack.zlaqge  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZLAQGE( M, N, A, LDA, R, C, ROWCND, COLCND, AMAX, EQUED )\n\n*  Purpose\n*  =======\n*\n*  ZLAQGE equilibrates a general M by N matrix A using the row and\n*  column scaling factors in the vectors R and C.\n*\n\n*  Arguments\n*  =========\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix A.  M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix A.  N >= 0.\n*\n*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)\n*          On entry, the M by N matrix A.\n*          On exit, the equilibrated matrix.  See EQUED for the form of\n*          the equilibrated matrix.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(M,1).\n*\n*  R       (input) DOUBLE PRECISION array, dimension (M)\n*          The row scale factors for A.\n*\n*  C       (input) DOUBLE PRECISION array, dimension (N)\n*          The column scale factors for A.\n*\n*  ROWCND  (input) DOUBLE PRECISION\n*          Ratio of the smallest R(i) to the largest R(i).\n*\n*  COLCND  (input) DOUBLE PRECISION\n*          Ratio of the smallest C(i) to the largest C(i).\n*\n*  AMAX    (input) DOUBLE PRECISION\n*          Absolute value of largest matrix entry.\n*\n*  EQUED   (output) CHARACTER*1\n*          Specifies the form of equilibration that was done.\n*          = 'N':  No equilibration\n*          = 'R':  Row equilibration, i.e., A has been premultiplied by\n*                  diag(R).\n*          = 'C':  Column equilibration, i.e., A has been postmultiplied\n*                  by diag(C).\n*          = 'B':  Both row and column equilibration, i.e., A has been\n*                  replaced by diag(R) * A * diag(C).\n*\n*  Internal Parameters\n*  ===================\n*\n*  THRESH is a threshold value used to decide if row or column scaling\n*  should be done based on the ratio of the row or column scaling\n*  factors.  If ROWCND < THRESH, row scaling is done, and if\n*  COLCND < THRESH, column scaling is done.\n*\n*  LARGE and SMALL are threshold values used to decide if row scaling\n*  should be done based on the absolute size of the largest matrix\n*  element.  If AMAX > LARGE or AMAX < SMALL, row scaling is done.\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rb_a = argv[0];
  rb_r = argv[1];
  rb_c = argv[2];
  rb_rowcnd = argv[3];
  rb_colcnd = argv[4];
  rb_amax = argv[5];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (1th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (1th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_DCOMPLEX)
    rb_a = na_change_type(rb_a, NA_DCOMPLEX);
  a = NA_PTR_TYPE(rb_a, doublecomplex*);
  if (!NA_IsNArray(rb_c))
    rb_raise(rb_eArgError, "c (3th argument) must be NArray");
  if (NA_RANK(rb_c) != 1)
    rb_raise(rb_eArgError, "rank of c (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_c) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of c must be the same as shape 1 of a");
  if (NA_TYPE(rb_c) != NA_DFLOAT)
    rb_c = na_change_type(rb_c, NA_DFLOAT);
  c = NA_PTR_TYPE(rb_c, doublereal*);
  amax = NUM2DBL(rb_amax);
  colcnd = NUM2DBL(rb_colcnd);
  if (!NA_IsNArray(rb_r))
    rb_raise(rb_eArgError, "r (2th argument) must be NArray");
  if (NA_RANK(rb_r) != 1)
    rb_raise(rb_eArgError, "rank of r (2th argument) must be %d", 1);
  m = NA_SHAPE0(rb_r);
  if (NA_TYPE(rb_r) != NA_DFLOAT)
    rb_r = na_change_type(rb_r, NA_DFLOAT);
  r = NA_PTR_TYPE(rb_r, doublereal*);
  rowcnd = NUM2DBL(rb_rowcnd);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rb_a_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rb_a_out__, doublecomplex*);
  MEMCPY(a_out__, a, doublecomplex, NA_TOTAL(rb_a));
  rb_a = rb_a_out__;
  a = a_out__;

  zlaqge_(&m, &n, a, &lda, r, c, &rowcnd, &colcnd, &amax, &equed);

  rb_equed = rb_str_new(&equed,1);
  return rb_ary_new3(2, rb_equed, rb_a);
}

void
init_lapack_zlaqge(VALUE mLapack){
  rb_define_module_function(mLapack, "zlaqge", rb_zlaqge, -1);
}
