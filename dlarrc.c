#include "rb_lapack.h"

static VALUE
rb_dlarrc(int argc, VALUE *argv, VALUE self){
  VALUE rb_jobt;
  char jobt; 
  VALUE rb_vl;
  doublereal vl; 
  VALUE rb_vu;
  doublereal vu; 
  VALUE rb_d;
  doublereal *d; 
  VALUE rb_e;
  doublereal *e; 
  VALUE rb_pivmin;
  doublereal pivmin; 
  VALUE rb_eigcnt;
  integer eigcnt; 
  VALUE rb_lcnt;
  integer lcnt; 
  VALUE rb_rcnt;
  integer rcnt; 
  VALUE rb_info;
  integer info; 

  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  eigcnt, lcnt, rcnt, info = NumRu::Lapack.dlarrc( jobt, vl, vu, d, e, pivmin)\n    or\n  NumRu::Lapack.dlarrc  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLARRC( JOBT, N, VL, VU, D, E, PIVMIN, EIGCNT, LCNT, RCNT, INFO )\n\n*  Purpose\n*  =======\n*\n*  Find the number of eigenvalues of the symmetric tridiagonal matrix T\n*  that are in the interval (VL,VU] if JOBT = 'T', and of L D L^T\n*  if JOBT = 'L'.\n*\n\n*  Arguments\n*  =========\n*\n*  JOBT    (input) CHARACTER*1\n*          = 'T':  Compute Sturm count for matrix T.\n*          = 'L':  Compute Sturm count for matrix L D L^T.\n*\n*  N       (input) INTEGER\n*          The order of the matrix. N > 0.\n*\n*  VL      (input) DOUBLE PRECISION\n*  VU      (input) DOUBLE PRECISION\n*          The lower and upper bounds for the eigenvalues.\n*\n*  D       (input) DOUBLE PRECISION array, dimension (N)\n*          JOBT = 'T': The N diagonal elements of the tridiagonal matrix T.\n*          JOBT = 'L': The N diagonal elements of the diagonal matrix D.\n*\n*  E       (input) DOUBLE PRECISION array, dimension (N)\n*          JOBT = 'T': The N-1 offdiagonal elements of the matrix T.\n*          JOBT = 'L': The N-1 offdiagonal elements of the matrix L.\n*\n*  PIVMIN  (input) DOUBLE PRECISION\n*          The minimum pivot in the Sturm sequence for T.\n*\n*  EIGCNT  (output) INTEGER\n*          The number of eigenvalues of the symmetric tridiagonal matrix T\n*          that are in the interval (VL,VU]\n*\n*  LCNT    (output) INTEGER\n*  RCNT    (output) INTEGER\n*          The left and right negcounts of the interval.\n*\n*  INFO    (output) INTEGER\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Beresford Parlett, University of California, Berkeley, USA\n*     Jim Demmel, University of California, Berkeley, USA\n*     Inderjit Dhillon, University of Texas, Austin, USA\n*     Osni Marques, LBNL/NERSC, USA\n*     Christof Voemel, University of California, Berkeley, USA\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rb_jobt = argv[0];
  rb_vl = argv[1];
  rb_vu = argv[2];
  rb_d = argv[3];
  rb_e = argv[4];
  rb_pivmin = argv[5];

  jobt = StringValueCStr(rb_jobt)[0];
  vl = NUM2DBL(rb_vl);
  vu = NUM2DBL(rb_vu);
  pivmin = NUM2DBL(rb_pivmin);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (4th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (4th argument) must be %d", 1);
  n = NA_SHAPE0(rb_d);
  if (NA_TYPE(rb_d) != NA_DFLOAT)
    rb_d = na_change_type(rb_d, NA_DFLOAT);
  d = NA_PTR_TYPE(rb_d, doublereal*);
  if (!NA_IsNArray(rb_e))
    rb_raise(rb_eArgError, "e (5th argument) must be NArray");
  if (NA_RANK(rb_e) != 1)
    rb_raise(rb_eArgError, "rank of e (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_e) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of e must be the same as shape 0 of d");
  if (NA_TYPE(rb_e) != NA_DFLOAT)
    rb_e = na_change_type(rb_e, NA_DFLOAT);
  e = NA_PTR_TYPE(rb_e, doublereal*);

  dlarrc_(&jobt, &n, &vl, &vu, d, e, &pivmin, &eigcnt, &lcnt, &rcnt, &info);

  rb_eigcnt = INT2NUM(eigcnt);
  rb_lcnt = INT2NUM(lcnt);
  rb_rcnt = INT2NUM(rcnt);
  rb_info = INT2NUM(info);
  return rb_ary_new3(4, rb_eigcnt, rb_lcnt, rb_rcnt, rb_info);
}

void
init_lapack_dlarrc(VALUE mLapack){
  rb_define_module_function(mLapack, "dlarrc", rb_dlarrc, -1);
}
