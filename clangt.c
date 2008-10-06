#include "rb_lapack.h"

extern VOID clangt_(real *__out__, char *norm, integer *n, complex *dl, complex *d, complex *du);
static VALUE
rb_clangt(int argc, VALUE *argv, VALUE self){
  VALUE rb_norm;
  char norm; 
  VALUE rb_dl;
  complex *dl; 
  VALUE rb_d;
  complex *d; 
  VALUE rb_du;
  complex *du; 
  VALUE rb___out__;
  real __out__; 

  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.clangt( norm, dl, d, du)\n    or\n  NumRu::Lapack.clangt  # print help\n\n\nFORTRAN MANUAL\n      REAL             FUNCTION CLANGT( NORM, N, DL, D, DU )\n\n*  Purpose\n*  =======\n*\n*  CLANGT  returns the value of the one norm,  or the Frobenius norm, or\n*  the  infinity norm,  or the  element of  largest absolute value  of a\n*  complex tridiagonal matrix A.\n*\n*  Description\n*  ===========\n*\n*  CLANGT returns the value\n*\n*     CLANGT = ( max(abs(A(i,j))), NORM = 'M' or 'm'\n*              (\n*              ( norm1(A),         NORM = '1', 'O' or 'o'\n*              (\n*              ( normI(A),         NORM = 'I' or 'i'\n*              (\n*              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'\n*\n*  where  norm1  denotes the  one norm of a matrix (maximum column sum),\n*  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and\n*  normF  denotes the  Frobenius norm of a matrix (square root of sum of\n*  squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.\n*\n\n*  Arguments\n*  =========\n*\n*  NORM    (input) CHARACTER*1\n*          Specifies the value to be returned in CLANGT as described\n*          above.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.  When N = 0, CLANGT is\n*          set to zero.\n*\n*  DL      (input) COMPLEX array, dimension (N-1)\n*          The (n-1) sub-diagonal elements of A.\n*\n*  D       (input) COMPLEX array, dimension (N)\n*          The diagonal elements of A.\n*\n*  DU      (input) COMPLEX array, dimension (N-1)\n*          The (n-1) super-diagonal elements of A.\n*\n\n*  =====================================================================\n*\n\n");
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
  if (NA_TYPE(rb_d) != NA_SCOMPLEX)
    rb_d = na_change_type(rb_d, NA_SCOMPLEX);
  d = NA_PTR_TYPE(rb_d, complex*);
  if (!NA_IsNArray(rb_dl))
    rb_raise(rb_eArgError, "dl (3th argument) must be NArray");
  if (NA_RANK(rb_dl) != 1)
    rb_raise(rb_eArgError, "rank of dl (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_dl) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of dl must be %d", n-1);
  if (NA_TYPE(rb_dl) != NA_SCOMPLEX)
    rb_dl = na_change_type(rb_dl, NA_SCOMPLEX);
  dl = NA_PTR_TYPE(rb_dl, complex*);
  if (!NA_IsNArray(rb_du))
    rb_raise(rb_eArgError, "du (4th argument) must be NArray");
  if (NA_RANK(rb_du) != 1)
    rb_raise(rb_eArgError, "rank of du (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_du) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of du must be %d", n-1);
  if (NA_TYPE(rb_du) != NA_SCOMPLEX)
    rb_du = na_change_type(rb_du, NA_SCOMPLEX);
  du = NA_PTR_TYPE(rb_du, complex*);

  clangt_(&__out__, &norm, &n, dl, d, du);

  rb___out__ = rb_float_new((double)__out__);
  return rb___out__;
}

void
init_lapack_clangt(VALUE mLapack){
  rb_define_module_function(mLapack, "clangt", rb_clangt, -1);
}
