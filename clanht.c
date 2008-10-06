#include "rb_lapack.h"

extern VOID clanht_(real *__out__, char *norm, integer *n, real *d, complex *e);
static VALUE
rb_clanht(int argc, VALUE *argv, VALUE self){
  VALUE rb_norm;
  char norm; 
  VALUE rb_d;
  real *d; 
  VALUE rb_e;
  complex *e; 
  VALUE rb___out__;
  real __out__; 

  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.clanht( norm, d, e)\n    or\n  NumRu::Lapack.clanht  # print help\n\n\nFORTRAN MANUAL\n      REAL             FUNCTION CLANHT( NORM, N, D, E )\n\n*  Purpose\n*  =======\n*\n*  CLANHT  returns the value of the one norm,  or the Frobenius norm, or\n*  the  infinity norm,  or the  element of  largest absolute value  of a\n*  complex Hermitian tridiagonal matrix A.\n*\n*  Description\n*  ===========\n*\n*  CLANHT returns the value\n*\n*     CLANHT = ( max(abs(A(i,j))), NORM = 'M' or 'm'\n*              (\n*              ( norm1(A),         NORM = '1', 'O' or 'o'\n*              (\n*              ( normI(A),         NORM = 'I' or 'i'\n*              (\n*              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'\n*\n*  where  norm1  denotes the  one norm of a matrix (maximum column sum),\n*  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and\n*  normF  denotes the  Frobenius norm of a matrix (square root of sum of\n*  squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.\n*\n\n*  Arguments\n*  =========\n*\n*  NORM    (input) CHARACTER*1\n*          Specifies the value to be returned in CLANHT as described\n*          above.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.  When N = 0, CLANHT is\n*          set to zero.\n*\n*  D       (input) REAL array, dimension (N)\n*          The diagonal elements of A.\n*\n*  E       (input) COMPLEX array, dimension (N-1)\n*          The (n-1) sub-diagonal or super-diagonal elements of A.\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rb_norm = argv[0];
  rb_d = argv[1];
  rb_e = argv[2];

  norm = StringValueCStr(rb_norm)[0];
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (2th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (2th argument) must be %d", 1);
  n = NA_SHAPE0(rb_d);
  if (NA_TYPE(rb_d) != NA_SFLOAT)
    rb_d = na_change_type(rb_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rb_d, real*);
  if (!NA_IsNArray(rb_e))
    rb_raise(rb_eArgError, "e (3th argument) must be NArray");
  if (NA_RANK(rb_e) != 1)
    rb_raise(rb_eArgError, "rank of e (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_e) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of e must be %d", n-1);
  if (NA_TYPE(rb_e) != NA_SCOMPLEX)
    rb_e = na_change_type(rb_e, NA_SCOMPLEX);
  e = NA_PTR_TYPE(rb_e, complex*);

  clanht_(&__out__, &norm, &n, d, e);

  rb___out__ = rb_float_new((double)__out__);
  return rb___out__;
}

void
init_lapack_clanht(VALUE mLapack){
  rb_define_module_function(mLapack, "clanht", rb_clanht, -1);
}
