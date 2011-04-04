#include "rb_lapack.h"

extern VOID dlamrg_(integer *n1, integer *n2, doublereal *a, integer *dtrd1, integer *dtrd2, integer *index);

static VALUE
rb_dlamrg(int argc, VALUE *argv, VALUE self){
  VALUE rb_n1;
  integer n1; 
  VALUE rb_n2;
  integer n2; 
  VALUE rb_a;
  doublereal *a; 
  VALUE rb_dtrd1;
  integer dtrd1; 
  VALUE rb_dtrd2;
  integer dtrd2; 
  VALUE rb_index;
  integer *index; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  index = NumRu::Lapack.dlamrg( n1, n2, a, dtrd1, dtrd2)\n    or\n  NumRu::Lapack.dlamrg  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLAMRG( N1, N2, A, DTRD1, DTRD2, INDEX )\n\n*  Purpose\n*  =======\n*\n*  DLAMRG will create a permutation list which will merge the elements\n*  of A (which is composed of two independently sorted sets) into a\n*  single set which is sorted in ascending order.\n*\n\n*  Arguments\n*  =========\n*\n*  N1     (input) INTEGER\n*  N2     (input) INTEGER\n*         These arguements contain the respective lengths of the two\n*         sorted lists to be merged.\n*\n*  A      (input) DOUBLE PRECISION array, dimension (N1+N2)\n*         The first N1 elements of A contain a list of numbers which\n*         are sorted in either ascending or descending order.  Likewise\n*         for the final N2 elements.\n*\n*  DTRD1  (input) INTEGER\n*  DTRD2  (input) INTEGER\n*         These are the strides to be taken through the array A.\n*         Allowable strides are 1 and -1.  They indicate whether a\n*         subset of A is sorted in ascending (DTRDx = 1) or descending\n*         (DTRDx = -1) order.\n*\n*  INDEX  (output) INTEGER array, dimension (N1+N2)\n*         On exit this array will contain a permutation such that\n*         if B( I ) = A( INDEX( I ) ) for I=1,N1+N2, then B will be\n*         sorted in ascending order.\n*\n\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      INTEGER            I, IND1, IND2, N1SV, N2SV\n*     ..\n\n");
    return Qnil;
  }
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rb_n1 = argv[0];
  rb_n2 = argv[1];
  rb_a = argv[2];
  rb_dtrd1 = argv[3];
  rb_dtrd2 = argv[4];

  dtrd2 = NUM2INT(rb_dtrd2);
  n1 = NUM2INT(rb_n1);
  dtrd1 = NUM2INT(rb_dtrd1);
  n2 = NUM2INT(rb_n2);
  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (3th argument) must be NArray");
  if (NA_RANK(rb_a) != 1)
    rb_raise(rb_eArgError, "rank of a (3th argument) must be %d", 1);
  if (NA_SHAPE0(rb_a) != (n1+n2))
    rb_raise(rb_eRuntimeError, "shape 0 of a must be %d", n1+n2);
  if (NA_TYPE(rb_a) != NA_DFLOAT)
    rb_a = na_change_type(rb_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rb_a, doublereal*);
  {
    int shape[1];
    shape[0] = n1+n2;
    rb_index = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  index = NA_PTR_TYPE(rb_index, integer*);

  dlamrg_(&n1, &n2, a, &dtrd1, &dtrd2, index);

  return rb_index;
}

void
init_lapack_dlamrg(VALUE mLapack){
  rb_define_module_function(mLapack, "dlamrg", rb_dlamrg, -1);
}
