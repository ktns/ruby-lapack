#include "rb_lapack.h"

static VALUE
rb_slaed6(int argc, VALUE *argv, VALUE self){
  VALUE rb_kniter;
  integer kniter; 
  VALUE rb_orgati;
  logical orgati; 
  VALUE rb_rho;
  real rho; 
  VALUE rb_d;
  real *d; 
  VALUE rb_z;
  real *z; 
  VALUE rb_finit;
  real finit; 
  VALUE rb_tau;
  real tau; 
  VALUE rb_info;
  integer info; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  tau, info = NumRu::Lapack.slaed6( kniter, orgati, rho, d, z, finit)\n    or\n  NumRu::Lapack.slaed6  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLAED6( KNITER, ORGATI, RHO, D, Z, FINIT, TAU, INFO )\n\n*  Purpose\n*  =======\n*\n*  SLAED6 computes the positive or negative root (closest to the origin)\n*  of\n*                   z(1)        z(2)        z(3)\n*  f(x) =   rho + --------- + ---------- + ---------\n*                  d(1)-x      d(2)-x      d(3)-x\n*\n*  It is assumed that\n*\n*        if ORGATI = .true. the root is between d(2) and d(3);\n*        otherwise it is between d(1) and d(2)\n*\n*  This routine will be called by SLAED4 when necessary. In most cases,\n*  the root sought is the smallest in magnitude, though it might not be\n*  in some extremely rare situations.\n*\n\n*  Arguments\n*  =========\n*\n*  KNITER       (input) INTEGER\n*               Refer to SLAED4 for its significance.\n*\n*  ORGATI       (input) LOGICAL\n*               If ORGATI is true, the needed root is between d(2) and\n*               d(3); otherwise it is between d(1) and d(2).  See\n*               SLAED4 for further details.\n*\n*  RHO          (input) REAL            \n*               Refer to the equation f(x) above.\n*\n*  D            (input) REAL array, dimension (3)\n*               D satisfies d(1) < d(2) < d(3).\n*\n*  Z            (input) REAL array, dimension (3)\n*               Each of the elements in z must be positive.\n*\n*  FINIT        (input) REAL            \n*               The value of f at 0. It is more accurate than the one\n*               evaluated inside this routine (if someone wants to do\n*               so).\n*\n*  TAU          (output) REAL            \n*               The root of the equation f(x).\n*\n*  INFO         (output) INTEGER\n*               = 0: successful exit\n*               > 0: if INFO = 1, failure to converge\n*\n\n*  Further Details\n*  ===============\n*\n*  30/06/99: Based on contributions by\n*     Ren-Cang Li, Computer Science Division, University of California\n*     at Berkeley, USA\n*\n*  10/02/03: This version has a few statements commented out for thread safety\n*     (machine parameters are computed on each entry). SJH.\n*\n*  05/10/06: Modified from a new version of Ren-Cang Li, use\n*     Gragg-Thornton-Warner cubic convergent scheme for better stability.\n*\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rb_kniter = argv[0];
  rb_orgati = argv[1];
  rb_rho = argv[2];
  rb_d = argv[3];
  rb_z = argv[4];
  rb_finit = argv[5];

  kniter = NUM2INT(rb_kniter);
  orgati = (rb_orgati == Qtrue);
  rho = (real)NUM2DBL(rb_rho);
  finit = (real)NUM2DBL(rb_finit);
  if (!NA_IsNArray(rb_d))
    rb_raise(rb_eArgError, "d (4th argument) must be NArray");
  if (NA_RANK(rb_d) != 1)
    rb_raise(rb_eArgError, "rank of d (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_d) != (3))
    rb_raise(rb_eRuntimeError, "shape 0 of d must be %d", 3);
  if (NA_TYPE(rb_d) != NA_SFLOAT)
    rb_d = na_change_type(rb_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rb_d, real*);
  if (!NA_IsNArray(rb_z))
    rb_raise(rb_eArgError, "z (5th argument) must be NArray");
  if (NA_RANK(rb_z) != 1)
    rb_raise(rb_eArgError, "rank of z (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_z) != (3))
    rb_raise(rb_eRuntimeError, "shape 0 of z must be %d", 3);
  if (NA_TYPE(rb_z) != NA_SFLOAT)
    rb_z = na_change_type(rb_z, NA_SFLOAT);
  z = NA_PTR_TYPE(rb_z, real*);

  slaed6_(&kniter, &orgati, &rho, d, z, &finit, &tau, &info);

  rb_tau = rb_float_new((double)tau);
  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_tau, rb_info);
}

void
init_lapack_slaed6(VALUE mLapack){
  rb_define_module_function(mLapack, "slaed6", rb_slaed6, -1);
}
