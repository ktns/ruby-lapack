#include "rb_lapack.h"

extern VOID zlangt_(doublereal *__out__, char *norm, integer *n, doublecomplex *dl, doublecomplex *d, doublecomplex *du);
static VALUE
rb_zlangt(int argc, VALUE *argv, VALUE self){
  VALUE rb_norm;
  char norm; 
  VALUE rb_dl;
  doublecomplex *dl; 
  VALUE rb_d;
  doublecomplex *d; 
  VALUE rb_du;
  doublecomplex *du; 
  VALUE rb___out__;
  doublereal __out__; 

  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.zlangt( norm, dl, d, du)\n    or\n  NumRu::Lapack.zlangt  # print help\n\n\nFORTRAN MANUAL\n      DOUBLE PRECISION FUNCTION ZLANGT( NORM, N, DL, D, DU )\n\n*  Purpose\n*  =======\n*\n*  ZLANGT  returns the value of the one norm,  or the Frobenius norm, or\n*  the  infinity norm,  or the  element of  largest absolute value  of a\n*  complex tridiagonal matrix A.\n*\n*  Description\n*  ===========\n*\n*  ZLANGT returns the value\n*\n*     ZLANGT = ( max(abs(A(i,j))), NORM = 'M' or 'm'\n*              (\n*              ( norm1(A),         NORM = '1', 'O' or 'o'\n*              (\n*              ( normI(A),         NORM = 'I' or 'i'\n*              (\n*              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'\n*\n*  where  norm1  denotes the  one norm of a matrix (maximum column sum),\n*  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and\n*  normF  denotes the  Frobenius norm of a matrix (square root of sum of\n*  squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.\n*\n\n*  Arguments\n*  =========\n*\n*  NORM    (input) CHARACTER*1\n*          Specifies the value to be returned in ZLANGT as described\n*          above.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.  When N = 0, ZLANGT is\n*          set to zero.\n*\n*  DL      (input) COMPLEX*16 array, dimension (N-1)\n*          The (n-1) sub-diagonal elements of A.\n*\n*  D       (input) COMPLEX*16 array, dimension (N)\n*          The diagonal elements of A.\n*\n*  DU      (input) COMPLEX*16 array, dimension (N-1)\n*          The (n-1) super-diagonal elements of A.\n*\n\n*  =====================================================================\n*\n\n");
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
  if (NA_TYPE(rb_d) != NA_DCOMPLEX)
    rb_d = na_change_type(rb_d, NA_DCOMPLEX);
  d = NA_PTR_TYPE(rb_d, doublecomplex*);
  if (!NA_IsNArray(rb_dl))
    rb_raise(rb_eArgError, "dl (3th argument) must be NArray");
  if (NA_RANK(rb_dl) != 1)
    rb_raise(rb_eArgError, "rank of dl (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_dl) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of dl must be %d", n-1);
  if (NA_TYPE(rb_dl) != NA_DCOMPLEX)
    rb_dl = na_change_type(rb_dl, NA_DCOMPLEX);
  dl = NA_PTR_TYPE(rb_dl, doublecomplex*);
  if (!NA_IsNArray(rb_du))
    rb_raise(rb_eArgError, "du (4th argument) must be NArray");
  if (NA_RANK(rb_du) != 1)
    rb_raise(rb_eArgError, "rank of du (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_du) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of du must be %d", n-1);
  if (NA_TYPE(rb_du) != NA_DCOMPLEX)
    rb_du = na_change_type(rb_du, NA_DCOMPLEX);
  du = NA_PTR_TYPE(rb_du, doublecomplex*);

  zlangt_(&__out__, &norm, &n, dl, d, du);

  rb___out__ = rb_float_new((double)__out__);
  return rb___out__;
}

void
init_lapack_zlangt(VALUE mLapack){
  rb_define_module_function(mLapack, "zlangt", rb_zlangt, -1);
}
