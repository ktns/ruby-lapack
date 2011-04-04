#include "rb_lapack.h"

extern doublereal zlange_(char *norm, integer *m, integer *n, doublecomplex *a, integer *lda, doublereal *work);

static VALUE
rb_zlange(int argc, VALUE *argv, VALUE self){
  VALUE rb_norm;
  char norm; 
  VALUE rb_m;
  integer m; 
  VALUE rb_a;
  doublecomplex *a; 
  VALUE rb___out__;
  doublereal __out__; 
  doublereal *work;

  integer lda;
  integer n;
  integer lwork;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.zlange( norm, m, a)\n    or\n  NumRu::Lapack.zlange  # print help\n\n\nFORTRAN MANUAL\n      DOUBLE PRECISION FUNCTION ZLANGE( NORM, M, N, A, LDA, WORK )\n\n*  Purpose\n*  =======\n*\n*  ZLANGE  returns the value of the one norm,  or the Frobenius norm, or\n*  the  infinity norm,  or the  element of  largest absolute value  of a\n*  complex matrix A.\n*\n*  Description\n*  ===========\n*\n*  ZLANGE returns the value\n*\n*     ZLANGE = ( max(abs(A(i,j))), NORM = 'M' or 'm'\n*              (\n*              ( norm1(A),         NORM = '1', 'O' or 'o'\n*              (\n*              ( normI(A),         NORM = 'I' or 'i'\n*              (\n*              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'\n*\n*  where  norm1  denotes the  one norm of a matrix (maximum column sum),\n*  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and\n*  normF  denotes the  Frobenius norm of a matrix (square root of sum of\n*  squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.\n*\n\n*  Arguments\n*  =========\n*\n*  NORM    (input) CHARACTER*1\n*          Specifies the value to be returned in ZLANGE as described\n*          above.\n*\n*  M       (input) INTEGER\n*          The number of rows of the matrix A.  M >= 0.  When M = 0,\n*          ZLANGE is set to zero.\n*\n*  N       (input) INTEGER\n*          The number of columns of the matrix A.  N >= 0.  When N = 0,\n*          ZLANGE is set to zero.\n*\n*  A       (input) COMPLEX*16 array, dimension (LDA,N)\n*          The m by n matrix A.\n*\n*  LDA     (input) INTEGER\n*          The leading dimension of the array A.  LDA >= max(M,1).\n*\n*  WORK    (workspace) DOUBLE PRECISION array, dimension (MAX(1,LWORK)),\n*          where LWORK >= M when NORM = 'I'; otherwise, WORK is not\n*          referenced.\n*\n\n* =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_norm = argv[0];
  rb_m = argv[1];
  rb_a = argv[2];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_DCOMPLEX)
    rb_a = na_change_type(rb_a, NA_DCOMPLEX);
  a = NA_PTR_TYPE(rb_a, doublecomplex*);
  m = NUM2INT(rb_m);
  norm = StringValueCStr(rb_norm)[0];
  lwork = lsame_(&norm,"I") ? m : 0;
  work = ALLOC_N(doublereal, (MAX(1,lwork)));

  __out__ = zlange_(&norm, &m, &n, a, &lda, work);

  free(work);
  rb___out__ = rb_float_new((double)__out__);
  return rb___out__;
}

void
init_lapack_zlange(VALUE mLapack){
  rb_define_module_function(mLapack, "zlange", rb_zlange, -1);
}
