#include "rb_lapack.h"

extern doublereal dlangb_(char *norm, integer *n, integer *kl, integer *ku, doublereal *ab, integer *ldab, doublereal *work);

static VALUE
rb_dlangb(int argc, VALUE *argv, VALUE self){
  VALUE rb_norm;
  char norm; 
  VALUE rb_kl;
  integer kl; 
  VALUE rb_ku;
  integer ku; 
  VALUE rb_ab;
  doublereal *ab; 
  VALUE rb___out__;
  doublereal __out__; 
  doublereal *work;

  integer ldab;
  integer n;
  integer lwork;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.dlangb( norm, kl, ku, ab)\n    or\n  NumRu::Lapack.dlangb  # print help\n\n\nFORTRAN MANUAL\n      DOUBLE PRECISION FUNCTION DLANGB( NORM, N, KL, KU, AB, LDAB, WORK )\n\n*  Purpose\n*  =======\n*\n*  DLANGB  returns the value of the one norm,  or the Frobenius norm, or\n*  the  infinity norm,  or the element of  largest absolute value  of an\n*  n by n band matrix  A,  with kl sub-diagonals and ku super-diagonals.\n*\n*  Description\n*  ===========\n*\n*  DLANGB returns the value\n*\n*     DLANGB = ( max(abs(A(i,j))), NORM = 'M' or 'm'\n*              (\n*              ( norm1(A),         NORM = '1', 'O' or 'o'\n*              (\n*              ( normI(A),         NORM = 'I' or 'i'\n*              (\n*              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'\n*\n*  where  norm1  denotes the  one norm of a matrix (maximum column sum),\n*  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and\n*  normF  denotes the  Frobenius norm of a matrix (square root of sum of\n*  squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.\n*\n\n*  Arguments\n*  =========\n*\n*  NORM    (input) CHARACTER*1\n*          Specifies the value to be returned in DLANGB as described\n*          above.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.  When N = 0, DLANGB is\n*          set to zero.\n*\n*  KL      (input) INTEGER\n*          The number of sub-diagonals of the matrix A.  KL >= 0.\n*\n*  KU      (input) INTEGER\n*          The number of super-diagonals of the matrix A.  KU >= 0.\n*\n*  AB      (input) DOUBLE PRECISION array, dimension (LDAB,N)\n*          The band matrix A, stored in rows 1 to KL+KU+1.  The j-th\n*          column of A is stored in the j-th column of the array AB as\n*          follows:\n*          AB(ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(n,j+kl).\n*\n*  LDAB    (input) INTEGER\n*          The leading dimension of the array AB.  LDAB >= KL+KU+1.\n*\n*  WORK    (workspace) DOUBLE PRECISION array, dimension (MAX(1,LWORK)),\n*          where LWORK >= N when NORM = 'I'; otherwise, WORK is not\n*          referenced.\n*\n\n* =====================================================================\n*\n*\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_norm = argv[0];
  rb_kl = argv[1];
  rb_ku = argv[2];
  rb_ab = argv[3];

  if (!NA_IsNArray(rb_ab))
    rb_raise(rb_eArgError, "ab (4th argument) must be NArray");
  if (NA_RANK(rb_ab) != 2)
    rb_raise(rb_eArgError, "rank of ab (4th argument) must be %d", 2);
  n = NA_SHAPE1(rb_ab);
  ldab = NA_SHAPE0(rb_ab);
  if (NA_TYPE(rb_ab) != NA_DFLOAT)
    rb_ab = na_change_type(rb_ab, NA_DFLOAT);
  ab = NA_PTR_TYPE(rb_ab, doublereal*);
  kl = NUM2INT(rb_kl);
  norm = StringValueCStr(rb_norm)[0];
  ku = NUM2INT(rb_ku);
  lwork = lsame_(&norm,"I") ? n : 0;
  work = ALLOC_N(doublereal, (MAX(1,lwork)));

  __out__ = dlangb_(&norm, &n, &kl, &ku, ab, &ldab, work);

  free(work);
  rb___out__ = rb_float_new((double)__out__);
  return rb___out__;
}

void
init_lapack_dlangb(VALUE mLapack){
  rb_define_module_function(mLapack, "dlangb", rb_dlangb, -1);
}
