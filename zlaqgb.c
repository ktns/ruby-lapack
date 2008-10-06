#include "rb_lapack.h"

static VALUE
rb_zlaqgb(int argc, VALUE *argv, VALUE self){
  VALUE rb_kl;
  integer kl; 
  VALUE rb_ku;
  integer ku; 
  VALUE rb_ab;
  doublecomplex *ab; 
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
  VALUE rb_ab_out__;
  doublecomplex *ab_out__;

  integer ldab;
  integer n;
  integer m;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  equed, ab = NumRu::Lapack.zlaqgb( kl, ku, ab, r, c, rowcnd, colcnd, amax)\n    or\n  NumRu::Lapack.zlaqgb  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZLAQGB( M, N, KL, KU, AB, LDAB, R, C, ROWCND, COLCND, AMAX, EQUED )\n\n*  Purpose\n*  =======\n*\n*  ZLAQGB equilibrates a general M by N band matrix A with KL\n*  subdiagonals and KU superdiagonals using the row and scaling factors\n*  in the vectors R and C.\n*\n\n*  Arguments\n*  =========\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix A.  M >= 0.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix A.  N >= 0.\n*\n*  KL      (input) INTEGER\n*          The number of subdiagonals within the band of A.  KL >= 0.\n*\n*  KU      (input) INTEGER\n*          The number of superdiagonals within the band of A.  KU >= 0.\n*\n*  AB      (input/output) COMPLEX*16 array, dimension (LDAB,N)\n*          On entry, the matrix A in band storage, in rows 1 to KL+KU+1.\n*          The j-th column of A is stored in the j-th column of the\n*          array AB as follows:\n*          AB(ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)\n*\n*          On exit, the equilibrated matrix, in the same storage format\n*          as A.  See EQUED for the form of the equilibrated matrix.\n*\n*  LDAB    (input) INTEGER\n*          The leading dimension of the array AB.  LDA >= KL+KU+1.\n*\n*  R       (input) DOUBLE PRECISION array, dimension (M)\n*          The row scale factors for A.\n*\n*  C       (input) DOUBLE PRECISION array, dimension (N)\n*          The column scale factors for A.\n*\n*  ROWCND  (input) DOUBLE PRECISION\n*          Ratio of the smallest R(i) to the largest R(i).\n*\n*  COLCND  (input) DOUBLE PRECISION\n*          Ratio of the smallest C(i) to the largest C(i).\n*\n*  AMAX    (input) DOUBLE PRECISION\n*          Absolute value of largest matrix entry.\n*\n*  EQUED   (output) CHARACTER*1\n*          Specifies the form of equilibration that was done.\n*          = 'N':  No equilibration\n*          = 'R':  Row equilibration, i.e., A has been premultiplied by\n*                  diag(R).\n*          = 'C':  Column equilibration, i.e., A has been postmultiplied\n*                  by diag(C).\n*          = 'B':  Both row and column equilibration, i.e., A has been\n*                  replaced by diag(R) * A * diag(C).\n*\n*  Internal Parameters\n*  ===================\n*\n*  THRESH is a threshold value used to decide if row or column scaling\n*  should be done based on the ratio of the row or column scaling\n*  factors.  If ROWCND < THRESH, row scaling is done, and if\n*  COLCND < THRESH, column scaling is done.\n*\n*  LARGE and SMALL are threshold values used to decide if row scaling\n*  should be done based on the absolute size of the largest matrix\n*  element.  If AMAX > LARGE or AMAX < SMALL, row scaling is done.\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 8)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 8)", argc);
  rb_kl = argv[0];
  rb_ku = argv[1];
  rb_ab = argv[2];
  rb_r = argv[3];
  rb_c = argv[4];
  rb_rowcnd = argv[5];
  rb_colcnd = argv[6];
  rb_amax = argv[7];

  kl = NUM2INT(rb_kl);
  ku = NUM2INT(rb_ku);
  rowcnd = NUM2DBL(rb_rowcnd);
  colcnd = NUM2DBL(rb_colcnd);
  amax = NUM2DBL(rb_amax);
  if (!NA_IsNArray(rb_ab))
    rb_raise(rb_eArgError, "ab (3th argument) must be NArray");
  if (NA_RANK(rb_ab) != 2)
    rb_raise(rb_eArgError, "rank of ab (3th argument) must be %d", 2);
  ldab = NA_SHAPE0(rb_ab);
  n = NA_SHAPE1(rb_ab);
  if (NA_TYPE(rb_ab) != NA_DCOMPLEX)
    rb_ab = na_change_type(rb_ab, NA_DCOMPLEX);
  ab = NA_PTR_TYPE(rb_ab, doublecomplex*);
  if (!NA_IsNArray(rb_r))
    rb_raise(rb_eArgError, "r (4th argument) must be NArray");
  if (NA_RANK(rb_r) != 1)
    rb_raise(rb_eArgError, "rank of r (4th argument) must be %d", 1);
  m = NA_SHAPE0(rb_r);
  if (NA_TYPE(rb_r) != NA_DFLOAT)
    rb_r = na_change_type(rb_r, NA_DFLOAT);
  r = NA_PTR_TYPE(rb_r, doublereal*);
  if (!NA_IsNArray(rb_c))
    rb_raise(rb_eArgError, "c (5th argument) must be NArray");
  if (NA_RANK(rb_c) != 1)
    rb_raise(rb_eArgError, "rank of c (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_c) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of c must be the same as shape 1 of ab");
  if (NA_TYPE(rb_c) != NA_DFLOAT)
    rb_c = na_change_type(rb_c, NA_DFLOAT);
  c = NA_PTR_TYPE(rb_c, doublereal*);
  {
    int shape[2];
    shape[0] = ldab;
    shape[1] = n;
    rb_ab_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  ab_out__ = NA_PTR_TYPE(rb_ab_out__, doublecomplex*);
  MEMCPY(ab_out__, ab, doublecomplex, NA_TOTAL(rb_ab));
  rb_ab = rb_ab_out__;
  ab = ab_out__;

  zlaqgb_(&m, &n, &kl, &ku, ab, &ldab, r, c, &rowcnd, &colcnd, &amax, &equed);

  rb_equed = rb_str_new(&equed,1);
  return rb_ary_new3(2, rb_equed, rb_ab);
}

void
init_lapack_zlaqgb(VALUE mLapack){
  rb_define_module_function(mLapack, "zlaqgb", rb_zlaqgb, -1);
}
