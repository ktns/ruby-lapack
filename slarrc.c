#include "rb_lapack.h"

extern VOID slarrc_(char *jobt, integer *n, real *vl, real *vu, real *d, real *e, real *pivmin, integer *eigcnt, integer *lcnt, integer *rcnt, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_slarrc(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_jobt;
  char jobt; 
  VALUE rblapack_vl;
  real vl; 
  VALUE rblapack_vu;
  real vu; 
  VALUE rblapack_d;
  real *d; 
  VALUE rblapack_e;
  real *e; 
  VALUE rblapack_pivmin;
  real pivmin; 
  VALUE rblapack_eigcnt;
  integer eigcnt; 
  VALUE rblapack_lcnt;
  integer lcnt; 
  VALUE rblapack_rcnt;
  integer rcnt; 
  VALUE rblapack_info;
  integer info; 

  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  eigcnt, lcnt, rcnt, info = NumRu::Lapack.slarrc( jobt, vl, vu, d, e, pivmin, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLARRC( JOBT, N, VL, VU, D, E, PIVMIN, EIGCNT, LCNT, RCNT, INFO )\n\n*  Purpose\n*  =======\n*\n*  Find the number of eigenvalues of the symmetric tridiagonal matrix T\n*  that are in the interval (VL,VU] if JOBT = 'T', and of L D L^T\n*  if JOBT = 'L'.\n*\n\n*  Arguments\n*  =========\n*\n*  JOBT    (input) CHARACTER*1\n*          = 'T':  Compute Sturm count for matrix T.\n*          = 'L':  Compute Sturm count for matrix L D L^T.\n*\n*  N       (input) INTEGER\n*          The order of the matrix. N > 0.\n*\n*  VL      (input) DOUBLE PRECISION\n*  VU      (input) DOUBLE PRECISION\n*          The lower and upper bounds for the eigenvalues.\n*\n*  D       (input) DOUBLE PRECISION array, dimension (N)\n*          JOBT = 'T': The N diagonal elements of the tridiagonal matrix T.\n*          JOBT = 'L': The N diagonal elements of the diagonal matrix D.\n*\n*  E       (input) DOUBLE PRECISION array, dimension (N)\n*          JOBT = 'T': The N-1 offdiagonal elements of the matrix T.\n*          JOBT = 'L': The N-1 offdiagonal elements of the matrix L.\n*\n*  PIVMIN  (input) REAL\n*          The minimum pivot in the Sturm sequence for T.\n*\n*  EIGCNT  (output) INTEGER\n*          The number of eigenvalues of the symmetric tridiagonal matrix T\n*          that are in the interval (VL,VU]\n*\n*  LCNT    (output) INTEGER\n*  RCNT    (output) INTEGER\n*          The left and right negcounts of the interval.\n*\n*  INFO    (output) INTEGER\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Beresford Parlett, University of California, Berkeley, USA\n*     Jim Demmel, University of California, Berkeley, USA\n*     Inderjit Dhillon, University of Texas, Austin, USA\n*     Osni Marques, LBNL/NERSC, USA\n*     Christof Voemel, University of California, Berkeley, USA\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  eigcnt, lcnt, rcnt, info = NumRu::Lapack.slarrc( jobt, vl, vu, d, e, pivmin, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rblapack_jobt = argv[0];
  rblapack_vl = argv[1];
  rblapack_vu = argv[2];
  rblapack_d = argv[3];
  rblapack_e = argv[4];
  rblapack_pivmin = argv[5];
  if (rb_options != Qnil) {
  }

  vl = (real)NUM2DBL(rblapack_vl);
  jobt = StringValueCStr(rblapack_jobt)[0];
  if (!NA_IsNArray(rblapack_d))
    rb_raise(rb_eArgError, "d (4th argument) must be NArray");
  if (NA_RANK(rblapack_d) != 1)
    rb_raise(rb_eArgError, "rank of d (4th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_d);
  if (NA_TYPE(rblapack_d) != NA_SFLOAT)
    rblapack_d = na_change_type(rblapack_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rblapack_d, real*);
  if (!NA_IsNArray(rblapack_e))
    rb_raise(rb_eArgError, "e (5th argument) must be NArray");
  if (NA_RANK(rblapack_e) != 1)
    rb_raise(rb_eArgError, "rank of e (5th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_e) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of e must be the same as shape 0 of d");
  if (NA_TYPE(rblapack_e) != NA_SFLOAT)
    rblapack_e = na_change_type(rblapack_e, NA_SFLOAT);
  e = NA_PTR_TYPE(rblapack_e, real*);
  vu = (real)NUM2DBL(rblapack_vu);
  pivmin = (real)NUM2DBL(rblapack_pivmin);

  slarrc_(&jobt, &n, &vl, &vu, d, e, &pivmin, &eigcnt, &lcnt, &rcnt, &info);

  rblapack_eigcnt = INT2NUM(eigcnt);
  rblapack_lcnt = INT2NUM(lcnt);
  rblapack_rcnt = INT2NUM(rcnt);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(4, rblapack_eigcnt, rblapack_lcnt, rblapack_rcnt, rblapack_info);
}

void
init_lapack_slarrc(VALUE mLapack){
  rb_define_module_function(mLapack, "slarrc", rblapack_slarrc, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
