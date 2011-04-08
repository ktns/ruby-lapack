#include "rb_lapack.h"

extern VOID dlamrg_(integer *n1, integer *n2, doublereal *a, integer *dtrd1, integer *dtrd2, integer *index);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dlamrg(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_n1;
  integer n1; 
  VALUE rblapack_n2;
  integer n2; 
  VALUE rblapack_a;
  doublereal *a; 
  VALUE rblapack_dtrd1;
  integer dtrd1; 
  VALUE rblapack_dtrd2;
  integer dtrd2; 
  VALUE rblapack_index;
  integer *index; 


  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  index = NumRu::Lapack.dlamrg( n1, n2, a, dtrd1, dtrd2, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLAMRG( N1, N2, A, DTRD1, DTRD2, INDEX )\n\n*  Purpose\n*  =======\n*\n*  DLAMRG will create a permutation list which will merge the elements\n*  of A (which is composed of two independently sorted sets) into a\n*  single set which is sorted in ascending order.\n*\n\n*  Arguments\n*  =========\n*\n*  N1     (input) INTEGER\n*  N2     (input) INTEGER\n*         These arguements contain the respective lengths of the two\n*         sorted lists to be merged.\n*\n*  A      (input) DOUBLE PRECISION array, dimension (N1+N2)\n*         The first N1 elements of A contain a list of numbers which\n*         are sorted in either ascending or descending order.  Likewise\n*         for the final N2 elements.\n*\n*  DTRD1  (input) INTEGER\n*  DTRD2  (input) INTEGER\n*         These are the strides to be taken through the array A.\n*         Allowable strides are 1 and -1.  They indicate whether a\n*         subset of A is sorted in ascending (DTRDx = 1) or descending\n*         (DTRDx = -1) order.\n*\n*  INDEX  (output) INTEGER array, dimension (N1+N2)\n*         On exit this array will contain a permutation such that\n*         if B( I ) = A( INDEX( I ) ) for I=1,N1+N2, then B will be\n*         sorted in ascending order.\n*\n\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      INTEGER            I, IND1, IND2, N1SV, N2SV\n*     ..\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  index = NumRu::Lapack.dlamrg( n1, n2, a, dtrd1, dtrd2, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rblapack_n1 = argv[0];
  rblapack_n2 = argv[1];
  rblapack_a = argv[2];
  rblapack_dtrd1 = argv[3];
  rblapack_dtrd2 = argv[4];
  if (rb_options != Qnil) {
  }

  dtrd2 = NUM2INT(rblapack_dtrd2);
  n1 = NUM2INT(rblapack_n1);
  dtrd1 = NUM2INT(rblapack_dtrd1);
  n2 = NUM2INT(rblapack_n2);
  if (!NA_IsNArray(rblapack_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rblapack_a) != 1)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_a) != (n1+n2))
    rb_raise(rb_eRuntimeError, "shape 0 of a must be %d", n1+n2);
  if (NA_TYPE(rblapack_a) != NA_DFLOAT)
    rblapack_a = na_change_type(rblapack_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rblapack_a, doublereal*);
  {
    int shape[1];
    shape[0] = n1+n2;
    rblapack_index = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  index = NA_PTR_TYPE(rblapack_index, integer*);

  dlamrg_(&n1, &n2, a, &dtrd1, &dtrd2, index);

  return rblapack_index;
}

void
init_lapack_dlamrg(VALUE mLapack){
  rb_define_module_function(mLapack, "dlamrg", rblapack_dlamrg, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
