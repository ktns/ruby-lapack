#include "rb_lapack.h"

extern VOID slangt_(real *__out__, char *norm, integer *n, real *dl, real *d, real *du);
static VALUE
rb_slangt(int argc, VALUE *argv, VALUE self){
  VALUE rb_norm;
  char norm; 
  VALUE rb_dl;
  real *dl; 
  VALUE rb_d;
  real *d; 
  VALUE rb_du;
  real *du; 
  VALUE rb___out__;
  real __out__; 

  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.slangt( norm, dl, d, du)\n    or\n  NumRu::Lapack.slangt  # print help\n\n\nFORTRAN MANUAL\n      REAL             FUNCTION SLANGT( NORM, N, DL, D, DU )\n\n*  Purpose\n*  =======\n*\n*  SLANGT  returns the value of the one norm,  or the Frobenius norm, or\n*  the  infinity norm,  or the  element of  largest absolute value  of a\n*  real tridiagonal matrix A.\n*\n*  Description\n*  ===========\n*\n*  SLANGT returns the value\n*\n*     SLANGT = ( max(abs(A(i,j))), NORM = 'M' or 'm'\n*              (\n*              ( norm1(A),         NORM = '1', 'O' or 'o'\n*              (\n*              ( normI(A),         NORM = 'I' or 'i'\n*              (\n*              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'\n*\n*  where  norm1  denotes the  one norm of a matrix (maximum column sum),\n*  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and\n*  normF  denotes the  Frobenius norm of a matrix (square root of sum of\n*  squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.\n*\n\n*  Arguments\n*  =========\n*\n*  NORM    (input) CHARACTER*1\n*          Specifies the value to be returned in SLANGT as described\n*          above.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.  When N = 0, SLANGT is\n*          set to zero.\n*\n*  DL      (input) REAL array, dimension (N-1)\n*          The (n-1) sub-diagonal elements of A.\n*\n*  D       (input) REAL array, dimension (N)\n*          The diagonal elements of A.\n*\n*  DU      (input) REAL array, dimension (N-1)\n*          The (n-1) super-diagonal elements of A.\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rb_norm = argv[0];
  rb_dl = argv[1];
  rb_d = argv[2];
  rb_du = argv[3];

  norm = StringValueCStr(rb_norm)[0];
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (2th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (2th argument) must be %d", 1);
  n = NA_SHAPE0(rb_d);
  if (NA_TYPE(rb_d) != NA_SFLOAT)
    rb_d = na_change_type(rb_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rb_d, real*);
  if (!NA_IsNArray(rb_dl))
    rb_raise(rb_eArgError, "dl (3th argument) must be NArray");
  if (NA_RANK(rb_dl) != 1)
    rb_raise(rb_eArgError, "rank of dl (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_dl) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of dl must be %d", n-1);
  if (NA_TYPE(rb_dl) != NA_SFLOAT)
    rb_dl = na_change_type(rb_dl, NA_SFLOAT);
  dl = NA_PTR_TYPE(rb_dl, real*);
  if (!NA_IsNArray(rb_du))
    rb_raise(rb_eArgError, "du (4th argument) must be NArray");
  if (NA_RANK(rb_du) != 1)
    rb_raise(rb_eArgError, "rank of du (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_du) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of du must be %d", n-1);
  if (NA_TYPE(rb_du) != NA_SFLOAT)
    rb_du = na_change_type(rb_du, NA_SFLOAT);
  du = NA_PTR_TYPE(rb_du, real*);

  slangt_(&__out__, &norm, &n, dl, d, du);

  rb___out__ = rb_float_new((double)__out__);
  return rb___out__;
}

void
init_lapack_slangt(VALUE mLapack){
  rb_define_module_function(mLapack, "slangt", rb_slangt, -1);
}
