#include "rb_lapack.h"

extern doublereal zlangt_(char *norm, integer *n, doublecomplex *dl, doublecomplex *d, doublecomplex *du);

static VALUE sHelp, sUsage;

static VALUE
rblapack_zlangt(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_norm;
  char norm; 
  VALUE rblapack_dl;
  doublecomplex *dl; 
  VALUE rblapack_d;
  doublecomplex *d; 
  VALUE rblapack_du;
  doublecomplex *du; 
  VALUE rblapack___out__;
  doublereal __out__; 

  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.zlangt( norm, dl, d, du, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      DOUBLE PRECISION FUNCTION ZLANGT( NORM, N, DL, D, DU )\n\n*  Purpose\n*  =======\n*\n*  ZLANGT  returns the value of the one norm,  or the Frobenius norm, or\n*  the  infinity norm,  or the  element of  largest absolute value  of a\n*  complex tridiagonal matrix A.\n*\n*  Description\n*  ===========\n*\n*  ZLANGT returns the value\n*\n*     ZLANGT = ( max(abs(A(i,j))), NORM = 'M' or 'm'\n*              (\n*              ( norm1(A),         NORM = '1', 'O' or 'o'\n*              (\n*              ( normI(A),         NORM = 'I' or 'i'\n*              (\n*              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'\n*\n*  where  norm1  denotes the  one norm of a matrix (maximum column sum),\n*  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and\n*  normF  denotes the  Frobenius norm of a matrix (square root of sum of\n*  squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.\n*\n\n*  Arguments\n*  =========\n*\n*  NORM    (input) CHARACTER*1\n*          Specifies the value to be returned in ZLANGT as described\n*          above.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.  When N = 0, ZLANGT is\n*          set to zero.\n*\n*  DL      (input) COMPLEX*16 array, dimension (N-1)\n*          The (n-1) sub-diagonal elements of A.\n*\n*  D       (input) COMPLEX*16 array, dimension (N)\n*          The diagonal elements of A.\n*\n*  DU      (input) COMPLEX*16 array, dimension (N-1)\n*          The (n-1) super-diagonal elements of A.\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.zlangt( norm, dl, d, du, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 4)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 4)", argc);
  rblapack_norm = argv[0];
  rblapack_dl = argv[1];
  rblapack_d = argv[2];
  rblapack_du = argv[3];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_d))
    rb_raise(rb_eArgError, "d (3th argument) must be NArray");
  if (NA_RANK(rblapack_d) != 1)
    rb_raise(rb_eArgError, "rank of d (3th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_d);
  if (NA_TYPE(rblapack_d) != NA_DCOMPLEX)
    rblapack_d = na_change_type(rblapack_d, NA_DCOMPLEX);
  d = NA_PTR_TYPE(rblapack_d, doublecomplex*);
  norm = StringValueCStr(rblapack_norm)[0];
  if (!NA_IsNArray(rblapack_du))
    rb_raise(rb_eArgError, "du (4th argument) must be NArray");
  if (NA_RANK(rblapack_du) != 1)
    rb_raise(rb_eArgError, "rank of du (4th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_du) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of du must be %d", n-1);
  if (NA_TYPE(rblapack_du) != NA_DCOMPLEX)
    rblapack_du = na_change_type(rblapack_du, NA_DCOMPLEX);
  du = NA_PTR_TYPE(rblapack_du, doublecomplex*);
  if (!NA_IsNArray(rblapack_dl))
    rb_raise(rb_eArgError, "dl (2th argument) must be NArray");
  if (NA_RANK(rblapack_dl) != 1)
    rb_raise(rb_eArgError, "rank of dl (2th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_dl) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of dl must be %d", n-1);
  if (NA_TYPE(rblapack_dl) != NA_DCOMPLEX)
    rblapack_dl = na_change_type(rblapack_dl, NA_DCOMPLEX);
  dl = NA_PTR_TYPE(rblapack_dl, doublecomplex*);

  __out__ = zlangt_(&norm, &n, dl, d, du);

  rblapack___out__ = rb_float_new((double)__out__);
  return rblapack___out__;
}

void
init_lapack_zlangt(VALUE mLapack){
  rb_define_module_function(mLapack, "zlangt", rblapack_zlangt, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
