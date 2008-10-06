#include "rb_lapack.h"

static VALUE
rb_zggbak(int argc, VALUE *argv, VALUE self){
  VALUE rb_job;
  char job; 
  VALUE rb_side;
  char side; 
  VALUE rb_ilo;
  integer ilo; 
  VALUE rb_ihi;
  integer ihi; 
  VALUE rb_lscale;
  doublereal *lscale; 
  VALUE rb_rscale;
  doublereal *rscale; 
  VALUE rb_v;
  doublecomplex *v; 
  VALUE rb_info;
  integer info; 
  VALUE rb_v_out__;
  doublecomplex *v_out__;

  integer n;
  integer ldv;
  integer m;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  info, v = NumRu::Lapack.zggbak( job, side, ilo, ihi, lscale, rscale, v)\n    or\n  NumRu::Lapack.zggbak  # print help\n\n\nFORTRAN MANUAL\n      SUBROUTINE ZGGBAK( JOB, SIDE, N, ILO, IHI, LSCALE, RSCALE, M, V, LDV, INFO )\n\n*  Purpose\n*  =======\n*\n*  ZGGBAK forms the right or left eigenvectors of a complex generalized\n*  eigenvalue problem A*x = lambda*B*x, by backward transformation on\n*  the computed eigenvectors of the balanced pair of matrices output by\n*  ZGGBAL.\n*\n\n*  Arguments\n*  =========\n*\n*  JOB     (input) CHARACTER*1\n*          Specifies the type of backward transformation required:\n*          = 'N':  do nothing, return immediately;\n*          = 'P':  do backward transformation for permutation only;\n*          = 'S':  do backward transformation for scaling only;\n*          = 'B':  do backward transformations for both permutation and\n*                  scaling.\n*          JOB must be the same as the argument JOB supplied to ZGGBAL.\n*\n*  SIDE    (input) CHARACTER*1\n*          = 'R':  V contains right eigenvectors;\n*          = 'L':  V contains left eigenvectors.\n*\n*  N       (input) INTEGER\n*          The number of rows of the matrix V.  N >= 0.\n*\n*  ILO     (input) INTEGER\n*  IHI     (input) INTEGER\n*          The integers ILO and IHI determined by ZGGBAL.\n*          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.\n*\n*  LSCALE  (input) DOUBLE PRECISION array, dimension (N)\n*          Details of the permutations and/or scaling factors applied\n*          to the left side of A and B, as returned by ZGGBAL.\n*\n*  RSCALE  (input) DOUBLE PRECISION array, dimension (N)\n*          Details of the permutations and/or scaling factors applied\n*          to the right side of A and B, as returned by ZGGBAL.\n*\n*  M       (input) INTEGER\n*          The number of columns of the matrix V.  M >= 0.\n*\n*  V       (input/output) COMPLEX*16 array, dimension (LDV,M)\n*          On entry, the matrix of right or left eigenvectors to be\n*          transformed, as returned by ZTGEVC.\n*          On exit, V is overwritten by the transformed eigenvectors.\n*\n*  LDV     (input) INTEGER\n*          The leading dimension of the matrix V. LDV >= max(1,N).\n*\n*  INFO    (output) INTEGER\n*          = 0:  successful exit.\n*          < 0:  if INFO = -i, the i-th argument had an illegal value.\n*\n\n*  Further Details\n*  ===============\n*\n*  See R.C. Ward, Balancing the generalized eigenvalue problem,\n*                 SIAM J. Sci. Stat. Comp. 2 (1981), 141-152.\n*\n*  =====================================================================\n*\n*     .. Local Scalars ..\n      LOGICAL            LEFTV, RIGHTV\n      INTEGER            I, K\n*     ..\n*     .. External Functions ..\n      LOGICAL            LSAME\n      EXTERNAL           LSAME\n*     ..\n*     .. External Subroutines ..\n      EXTERNAL           XERBLA, ZDSCAL, ZSWAP\n*     ..\n*     .. Intrinsic Functions ..\n      INTRINSIC          MAX\n*     ..\n\n");
    return Qnil;
  }
  if (argc != 7)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 7)", argc);
  rb_job = argv[0];
  rb_side = argv[1];
  rb_ilo = argv[2];
  rb_ihi = argv[3];
  rb_lscale = argv[4];
  rb_rscale = argv[5];
  rb_v = argv[6];

  job = StringValueCStr(rb_job)[0];
  side = StringValueCStr(rb_side)[0];
  ilo = NUM2INT(rb_ilo);
  ihi = NUM2INT(rb_ihi);
  if (!NA_IsNArray(rb_lscale))
    rb_raise(rb_eArgError, "lscale (5th argument) must be NArray");
  if (NA_RANK(rb_lscale) != 1)
    rb_raise(rb_eArgError, "rank of lscale (5th argument) must be %d", 1);
  n = NA_SHAPE0(rb_lscale);
  if (NA_TYPE(rb_lscale) != NA_DFLOAT)
    rb_lscale = na_change_type(rb_lscale, NA_DFLOAT);
  lscale = NA_PTR_TYPE(rb_lscale, doublereal*);
  if (!NA_IsNArray(rb_rscale))
    rb_raise(rb_eArgError, "rscale (6th argument) must be NArray");
  if (NA_RANK(rb_rscale) != 1)
    rb_raise(rb_eArgError, "rank of rscale (6th argument) must be %d", 1);
  if (NA_SHAPE0(rb_rscale) != n)
    rb_raise(rb_eRuntimeError, "shape 0 of rscale must be the same as shape 0 of lscale");
  if (NA_TYPE(rb_rscale) != NA_DFLOAT)
    rb_rscale = na_change_type(rb_rscale, NA_DFLOAT);
  rscale = NA_PTR_TYPE(rb_rscale, doublereal*);
  if (!NA_IsNArray(rb_v))
    rb_raise(rb_eArgError, "v (7th argument) must be NArray");
  if (NA_RANK(rb_v) != 2)
    rb_raise(rb_eArgError, "rank of v (7th argument) must be %d", 2);
  ldv = NA_SHAPE0(rb_v);
  m = NA_SHAPE1(rb_v);
  if (NA_TYPE(rb_v) != NA_DCOMPLEX)
    rb_v = na_change_type(rb_v, NA_DCOMPLEX);
  v = NA_PTR_TYPE(rb_v, doublecomplex*);
  {
    int shape[2];
    shape[0] = ldv;
    shape[1] = m;
    rb_v_out__ = na_make_object(NA_DCOMPLEX, 2, shape, cNArray);
  }
  v_out__ = NA_PTR_TYPE(rb_v_out__, doublecomplex*);
  MEMCPY(v_out__, v, doublecomplex, NA_TOTAL(rb_v));
  rb_v = rb_v_out__;
  v = v_out__;

  zggbak_(&job, &side, &n, &ilo, &ihi, lscale, rscale, &m, v, &ldv, &info);

  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_info, rb_v);
}

void
init_lapack_zggbak(VALUE mLapack){
  rb_define_module_function(mLapack, "zggbak", rb_zggbak, -1);
}
