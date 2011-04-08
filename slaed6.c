#include "rb_lapack.h"

extern VOID slaed6_(integer *kniter, logical *orgati, real *rho, real *d, real *z, real *finit, real *tau, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_slaed6(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_kniter;
  integer kniter; 
  VALUE rblapack_orgati;
  logical orgati; 
  VALUE rblapack_rho;
  real rho; 
  VALUE rblapack_d;
  real *d; 
  VALUE rblapack_z;
  real *z; 
  VALUE rblapack_finit;
  real finit; 
  VALUE rblapack_tau;
  real tau; 
  VALUE rblapack_info;
  integer info; 


  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  tau, info = NumRu::Lapack.slaed6( kniter, orgati, rho, d, z, finit, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLAED6( KNITER, ORGATI, RHO, D, Z, FINIT, TAU, INFO )\n\n*  Purpose\n*  =======\n*\n*  SLAED6 computes the positive or negative root (closest to the origin)\n*  of\n*                   z(1)        z(2)        z(3)\n*  f(x) =   rho + --------- + ---------- + ---------\n*                  d(1)-x      d(2)-x      d(3)-x\n*\n*  It is assumed that\n*\n*        if ORGATI = .true. the root is between d(2) and d(3);\n*        otherwise it is between d(1) and d(2)\n*\n*  This routine will be called by SLAED4 when necessary. In most cases,\n*  the root sought is the smallest in magnitude, though it might not be\n*  in some extremely rare situations.\n*\n\n*  Arguments\n*  =========\n*\n*  KNITER       (input) INTEGER\n*               Refer to SLAED4 for its significance.\n*\n*  ORGATI       (input) LOGICAL\n*               If ORGATI is true, the needed root is between d(2) and\n*               d(3); otherwise it is between d(1) and d(2).  See\n*               SLAED4 for further details.\n*\n*  RHO          (input) REAL            \n*               Refer to the equation f(x) above.\n*\n*  D            (input) REAL array, dimension (3)\n*               D satisfies d(1) < d(2) < d(3).\n*\n*  Z            (input) REAL array, dimension (3)\n*               Each of the elements in z must be positive.\n*\n*  FINIT        (input) REAL            \n*               The value of f at 0. It is more accurate than the one\n*               evaluated inside this routine (if someone wants to do\n*               so).\n*\n*  TAU          (output) REAL            \n*               The root of the equation f(x).\n*\n*  INFO         (output) INTEGER\n*               = 0: successful exit\n*               > 0: if INFO = 1, failure to converge\n*\n\n*  Further Details\n*  ===============\n*\n*  30/06/99: Based on contributions by\n*     Ren-Cang Li, Computer Science Division, University of California\n*     at Berkeley, USA\n*\n*  10/02/03: This version has a few statements commented out for thread safety\n*     (machine parameters are computed on each entry). SJH.\n*\n*  05/10/06: Modified from a new version of Ren-Cang Li, use\n*     Gragg-Thornton-Warner cubic convergent scheme for better stability.\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  tau, info = NumRu::Lapack.slaed6( kniter, orgati, rho, d, z, finit, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rblapack_kniter = argv[0];
  rblapack_orgati = argv[1];
  rblapack_rho = argv[2];
  rblapack_d = argv[3];
  rblapack_z = argv[4];
  rblapack_finit = argv[5];
  if (rb_options != Qnil) {
  }

  orgati = (rblapack_orgati == Qtrue);
  finit = (real)NUM2DBL(rblapack_finit);
  rho = (real)NUM2DBL(rblapack_rho);
  if (!NA_IsNArray(rblapack_z))
    rb_raise(rb_eArgError, "z (5th argument) must be NArray");
  if (NA_RANK(rblapack_z) != 1)
    rb_raise(rb_eArgError, "rank of z (5th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_z) != (3))
    rb_raise(rb_eRuntimeError, "shape 0 of z must be %d", 3);
  if (NA_TYPE(rblapack_z) != NA_SFLOAT)
    rblapack_z = na_change_type(rblapack_z, NA_SFLOAT);
  z = NA_PTR_TYPE(rblapack_z, real*);
  if (!NA_IsNArray(rblapack_d))
    rb_raise(rb_eArgError, "d (4th argument) must be NArray");
  if (NA_RANK(rblapack_d) != 1)
    rb_raise(rb_eArgError, "rank of d (4th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_d) != (3))
    rb_raise(rb_eRuntimeError, "shape 0 of d must be %d", 3);
  if (NA_TYPE(rblapack_d) != NA_SFLOAT)
    rblapack_d = na_change_type(rblapack_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rblapack_d, real*);
  kniter = NUM2INT(rblapack_kniter);

  slaed6_(&kniter, &orgati, &rho, d, z, &finit, &tau, &info);

  rblapack_tau = rb_float_new((double)tau);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(2, rblapack_tau, rblapack_info);
}

void
init_lapack_slaed6(VALUE mLapack){
  rb_define_module_function(mLapack, "slaed6", rblapack_slaed6, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
