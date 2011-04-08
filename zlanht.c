#include "rb_lapack.h"

extern doublereal zlanht_(char *norm, integer *n, doublereal *d, doublecomplex *e);

static VALUE sHelp, sUsage;

static VALUE
rblapack_zlanht(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_norm;
  char norm; 
  VALUE rblapack_d;
  doublereal *d; 
  VALUE rblapack_e;
  doublecomplex *e; 
  VALUE rblapack___out__;
  doublereal __out__; 

  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.zlanht( norm, d, e, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      DOUBLE PRECISION FUNCTION ZLANHT( NORM, N, D, E )\n\n*  Purpose\n*  =======\n*\n*  ZLANHT  returns the value of the one norm,  or the Frobenius norm, or\n*  the  infinity norm,  or the  element of  largest absolute value  of a\n*  complex Hermitian tridiagonal matrix A.\n*\n*  Description\n*  ===========\n*\n*  ZLANHT returns the value\n*\n*     ZLANHT = ( max(abs(A(i,j))), NORM = 'M' or 'm'\n*              (\n*              ( norm1(A),         NORM = '1', 'O' or 'o'\n*              (\n*              ( normI(A),         NORM = 'I' or 'i'\n*              (\n*              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'\n*\n*  where  norm1  denotes the  one norm of a matrix (maximum column sum),\n*  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and\n*  normF  denotes the  Frobenius norm of a matrix (square root of sum of\n*  squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.\n*\n\n*  Arguments\n*  =========\n*\n*  NORM    (input) CHARACTER*1\n*          Specifies the value to be returned in ZLANHT as described\n*          above.\n*\n*  N       (input) INTEGER\n*          The order of the matrix A.  N >= 0.  When N = 0, ZLANHT is\n*          set to zero.\n*\n*  D       (input) DOUBLE PRECISION array, dimension (N)\n*          The diagonal elements of A.\n*\n*  E       (input) COMPLEX*16 array, dimension (N-1)\n*          The (n-1) sub-diagonal or super-diagonal elements of A.\n*\n\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  __out__ = NumRu::Lapack.zlanht( norm, d, e, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 3)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 3)", argc);
  rblapack_norm = argv[0];
  rblapack_d = argv[1];
  rblapack_e = argv[2];
  if (rb_options != Qnil) {
  }

  if (!NA_IsNArray(rblapack_d))
    rb_raise(rb_eArgError, "d (2th argument) must be NArray");
  if (NA_RANK(rblapack_d) != 1)
    rb_raise(rb_eArgError, "rank of d (2th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_d);
  if (NA_TYPE(rblapack_d) != NA_DFLOAT)
    rblapack_d = na_change_type(rblapack_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rblapack_d, doublereal*);
  norm = StringValueCStr(rblapack_norm)[0];
  if (!NA_IsNArray(rblapack_e))
    rb_raise(rb_eArgError, "e (3th argument) must be NArray");
  if (NA_RANK(rblapack_e) != 1)
    rb_raise(rb_eArgError, "rank of e (3th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_e) != (n-1))
    rb_raise(rb_eRuntimeError, "shape 0 of e must be %d", n-1);
  if (NA_TYPE(rblapack_e) != NA_DCOMPLEX)
    rblapack_e = na_change_type(rblapack_e, NA_DCOMPLEX);
  e = NA_PTR_TYPE(rblapack_e, doublecomplex*);

  __out__ = zlanht_(&norm, &n, d, e);

  rblapack___out__ = rb_float_new((double)__out__);
  return rblapack___out__;
}

void
init_lapack_zlanht(VALUE mLapack){
  rb_define_module_function(mLapack, "zlanht", rblapack_zlanht, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
