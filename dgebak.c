#include "rb_lapack.h"

static VALUE
rb_dgebak(int argc, VALUE *argv, VALUE self){
  VALUE rb_job;
  char job; 
  VALUE rb_side;
  char side; 
  VALUE rb_ilo;
  integer ilo; 
  VALUE rb_ihi;
  integer ihi; 
  VALUE rb_scale;
  doublereal *scale; 
  VALUE rb_v;
  doublereal *v; 
  VALUE rb_info;
  integer info; 
  VALUE rb_v_out__;
  doublereal *v_out__;

  integer n;
  integer ldv;
  integer m;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, v = NumRu::Lapack.dgebak( job, side, ilo, ihi, scale, v)\n    or\n  NumRu::Lapack.dgebak  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE DGEBAK( JOB, SIDE, N, ILO, IHI, SCALE, M, V, LDV, INFO )\n\n*  Purpose\n*  =======\n*\n*  DGEBAK forms the right or left eigenvectors of a real general matrix\n*  by backward transformation on the computed eigenvectors of the\n*  balanced matrix output by DGEBAL.\n*\n\n*  Arguments\n*  =========\n*\n*  JOB     (input) CHARACTER*1\n*          Specifies the type of backward transformation required:\n*          = 'N', do nothing, return immediately;\n*          = 'P', do backward transformation for permutation only;\n*          = 'S', do backward transformation for scaling only;\n*          = 'B', do backward transformations for both permutation and\n*                 scaling.\n*          JOB must be the same as the argument JOB supplied to DGEBAL.\n*\n*  SIDE    (input) CHARACTER*1\n*          = 'R':  V contains right eigenvectors;\n*          = 'L':  V contains left eigenvectors.\n*\n*  N       (input) INTEGER\n*          The number of rows of the matrix V.  N >= 0.\n*\n*  ILO     (input) INTEGER\n*  IHI     (input) INTEGER\n*          The integers ILO and IHI determined by DGEBAL.\n*          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.\n*\n*  SCALE   (input) DOUBLE PRECISION array, dimension (N)\n*          Details of the permutation and scaling factors, as returned\n*          by DGEBAL.\n*\n*  M       (input) INTEGER\n*          The number of columns of the matrix V.  M >= 0.\n*\n*  V       (input/output) DOUBLE PRECISION array, dimension (LDV,M)\n*          On entry, the matrix of right or left eigenvectors to be\n*          transformed, as returned by DHSEIN or DTREVC.\n*          On exit, V is overwritten by the transformed eigenvectors.\n*\n*  LDV     (input) INTEGER\n*          The leading dimension of the array V. LDV >= max(1,N).\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*\n\n*  =====================================================================\n*\n\n");
    return Qnil;
  }
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rb_job = argv[0];
  rb_side = argv[1];
  rb_ilo = argv[2];
  rb_ihi = argv[3];
  rb_scale = argv[4];
  rb_v = argv[5];

  job = StringValueCStr(rb_job)[0];
  side = StringValueCStr(rb_side)[0];
  ilo = NUM2INT(rb_ilo);
  ihi = NUM2INT(rb_ihi);
  if (!NA_IsNArray(rb_scale))
    rb_raise(rb_eArgError, "scale (5th argument) must be NArray");
  if (NA_RANK(rb_scale) != 1)
    rb_raise(rb_eArgError, "rank of scale (5th argument) must be %d", 1);
  n = NA_SHAPE0(rb_scale);
  if (NA_TYPE(rb_scale) != NA_DFLOAT)
    rb_scale = na_change_type(rb_scale, NA_DFLOAT);
  scale = NA_PTR_TYPE(rb_scale, doublereal*);
  if (!NA_IsNArray(rb_v))
    rb_raise(rb_eArgError, "v (6th argument) must be NArray");
  if (NA_RANK(rb_v) != 2)
    rb_raise(rb_eArgError, "rank of v (6th argument) must be %d", 2);
  ldv = NA_SHAPE0(rb_v);
  m = NA_SHAPE1(rb_v);
  if (NA_TYPE(rb_v) != NA_DFLOAT)
    rb_v = na_change_type(rb_v, NA_DFLOAT);
  v = NA_PTR_TYPE(rb_v, doublereal*);
  {
    int shape[2];
    shape[0] = ldv;
    shape[1] = m;
    rb_v_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  v_out__ = NA_PTR_TYPE(rb_v_out__, doublereal*);
  MEMCPY(v_out__, v, doublereal, NA_TOTAL(rb_v));
  rb_v = rb_v_out__;
  v = v_out__;

  dgebak_(&job, &side, &n, &ilo, &ihi, scale, &m, v, &ldv, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_info, rb_v);
}

void
init_lapack_dgebak(VALUE mLapack){
  rb_define_module_function(mLapack, "dgebak", rb_dgebak, -1);
}
