#include "rb_lapack.h"

extern VOID slarra_(integer *n, real *d, real *e, real *e2, real *spltol, real *tnrm, integer *nsplit, integer *isplit, integer *info);

static VALUE sHelp, sUsage;

static VALUE
rblapack_slarra(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_d;
  real *d; 
  VALUE rblapack_e;
  real *e; 
  VALUE rblapack_e2;
  real *e2; 
  VALUE rblapack_spltol;
  real spltol; 
  VALUE rblapack_tnrm;
  real tnrm; 
  VALUE rblapack_nsplit;
  integer nsplit; 
  VALUE rblapack_isplit;
  integer *isplit; 
  VALUE rblapack_info;
  integer info; 
  VALUE rblapack_e_out__;
  real *e_out__;
  VALUE rblapack_e2_out__;
  real *e2_out__;

  integer n;

  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  nsplit, isplit, info, e, e2 = NumRu::Lapack.slarra( d, e, e2, spltol, tnrm, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE SLARRA( N, D, E, E2, SPLTOL, TNRM, NSPLIT, ISPLIT, INFO )\n\n*  Purpose\n*  =======\n*\n*  Compute the splitting points with threshold SPLTOL.\n*  SLARRA sets any \"small\" off-diagonal elements to zero.\n*\n\n*  Arguments\n*  =========\n*\n*  N       (input) INTEGER\n*          The order of the matrix. N > 0.\n*\n*  D       (input) REAL             array, dimension (N)\n*          On entry, the N diagonal elements of the tridiagonal\n*          matrix T.\n*\n*  E       (input/output) REAL             array, dimension (N)\n*          On entry, the first (N-1) entries contain the subdiagonal\n*          elements of the tridiagonal matrix T; E(N) need not be set.\n*          On exit, the entries E( ISPLIT( I ) ), 1 <= I <= NSPLIT,\n*          are set to zero, the other entries of E are untouched.\n*\n*  E2      (input/output) REAL             array, dimension (N)\n*          On entry, the first (N-1) entries contain the SQUARES of the\n*          subdiagonal elements of the tridiagonal matrix T;\n*          E2(N) need not be set.\n*          On exit, the entries E2( ISPLIT( I ) ),\n*          1 <= I <= NSPLIT, have been set to zero\n*\n*  SPLTOL (input) REAL            \n*          The threshold for splitting. Two criteria can be used:\n*          SPLTOL<0 : criterion based on absolute off-diagonal value\n*          SPLTOL>0 : criterion that preserves relative accuracy\n*\n*  TNRM (input) REAL            \n*          The norm of the matrix.\n*\n*  NSPLIT  (output) INTEGER\n*          The number of blocks T splits into. 1 <= NSPLIT <= N.\n*\n*  ISPLIT  (output) INTEGER array, dimension (N)\n*          The splitting points, at which T breaks up into blocks.\n*          The first block consists of rows/columns 1 to ISPLIT(1),\n*          the second of rows/columns ISPLIT(1)+1 through ISPLIT(2),\n*          etc., and the NSPLIT-th consists of rows/columns\n*          ISPLIT(NSPLIT-1)+1 through ISPLIT(NSPLIT)=N.\n*\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Beresford Parlett, University of California, Berkeley, USA\n*     Jim Demmel, University of California, Berkeley, USA\n*     Inderjit Dhillon, University of Texas, Austin, USA\n*     Osni Marques, LBNL/NERSC, USA\n*     Christof Voemel, University of California, Berkeley, USA\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  nsplit, isplit, info, e, e2 = NumRu::Lapack.slarra( d, e, e2, spltol, tnrm, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 5)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 5)", argc);
  rblapack_d = argv[0];
  rblapack_e = argv[1];
  rblapack_e2 = argv[2];
  rblapack_spltol = argv[3];
  rblapack_tnrm = argv[4];
  if (rb_options != Qnil) {
  }

  tnrm = (real)NUM2DBL(rblapack_tnrm);
  if (!NA_IsNArray(rblapack_e2))
    rb_raise(rb_eArgError, "e2 (3th argument) must be NArray");
  if (NA_RANK(rblapack_e2) != 1)
    rb_raise(rb_eArgError, "rank of e2 (3th argument) must be %d", 1);
  n = NA_SHAPE0(rblapack_e2);
  if (NA_TYPE(rblapack_e2) != NA_SFLOAT)
    rblapack_e2 = na_change_type(rblapack_e2, NA_SFLOAT);
  e2 = NA_PTR_TYPE(rblapack_e2, real*);
  spltol = (real)NUM2DBL(rblapack_spltol);
  if (!NA_IsNArray(rblapack_d))
    rb_raise(rb_eArgError, "d (1th argument) must be NArray");
  if (NA_RANK(rblapack_d) != 1)
    rb_raise(rb_eArgError, "rank of d (1th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_d) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of d must be the same as shape 0 of e2");
  if (NA_TYPE(rblapack_d) != NA_SFLOAT)
    rblapack_d = na_change_type(rblapack_d, NA_SFLOAT);
  d = NA_PTR_TYPE(rblapack_d, real*);
  if (!NA_IsNArray(rblapack_e))
    rb_raise(rb_eArgError, "e (2th argument) must be NArray");
  if (NA_RANK(rblapack_e) != 1)
    rb_raise(rb_eArgError, "rank of e (2th argument) must be %d", 1);
  if (NA_SHAPE0(rblapack_e) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of e must be the same as shape 0 of e2");
  if (NA_TYPE(rblapack_e) != NA_SFLOAT)
    rblapack_e = na_change_type(rblapack_e, NA_SFLOAT);
  e = NA_PTR_TYPE(rblapack_e, real*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_isplit = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  isplit = NA_PTR_TYPE(rblapack_isplit, integer*);
  {
    int shape[1];
    shape[0] = n;
    rblapack_e_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  e_out__ = NA_PTR_TYPE(rblapack_e_out__, real*);
  MEMCPY(e_out__, e, real, NA_TOTAL(rblapack_e));
  rblapack_e = rblapack_e_out__;
  e = e_out__;
  {
    int shape[1];
    shape[0] = n;
    rblapack_e2_out__ = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  e2_out__ = NA_PTR_TYPE(rblapack_e2_out__, real*);
  MEMCPY(e2_out__, e2, real, NA_TOTAL(rblapack_e2));
  rblapack_e2 = rblapack_e2_out__;
  e2 = e2_out__;

  slarra_(&n, d, e, e2, &spltol, &tnrm, &nsplit, isplit, &info);

  rblapack_nsplit = INT2NUM(nsplit);
  rblapack_info = INT2NUM(info);
  return rb_ary_new3(5, rblapack_nsplit, rblapack_isplit, rblapack_info, rblapack_e, rblapack_e2);
}

void
init_lapack_slarra(VALUE mLapack){
  rb_define_module_function(mLapack, "slarra", rblapack_slarra, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
