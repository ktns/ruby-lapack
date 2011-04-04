#include "rb_lapack.h"

extern VOID dspgvx_(integer *itype, char *jobz, char *range, char *uplo, integer *n, doublereal *ap, doublereal *bp, doublereal *vl, doublereal *vu, integer *il, integer *iu, doublereal *abstol, integer *m, doublereal *w, doublereal *z, integer *ldz, doublereal *work, integer *iwork, integer *ifail, integer *info);

static VALUE
rb_dspgvx(int argc, VALUE *argv, VALUE self){
  VALUE rb_itype;
  integer itype; 
  VALUE rb_jobz;
  char jobz; 
  VALUE rb_range;
  char range; 
  VALUE rb_uplo;
  char uplo; 
  VALUE rb_ap;
  doublereal *ap; 
  VALUE rb_bp;
  doublereal *bp; 
  VALUE rb_vl;
  doublereal vl; 
  VALUE rb_vu;
  doublereal vu; 
  VALUE rb_il;
  integer il; 
  VALUE rb_iu;
  integer iu; 
  VALUE rb_abstol;
  doublereal abstol; 
  VALUE rb_ldz;
  integer ldz; 
  VALUE rb_m;
  integer m; 
  VALUE rb_w;
  doublereal *w; 
  VALUE rb_z;
  doublereal *z; 
  VALUE rb_ifail;
  integer *ifail; 
  VALUE rb_info;
  integer info; 
  VALUE rb_ap_out__;
  doublereal *ap_out__;
  VALUE rb_bp_out__;
  doublereal *bp_out__;
  doublereal *work;
  integer *iwork;

  integer ldap;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  m, w, z, ifail, info, ap, bp = NumRu::Lapack.dspgvx( itype, jobz, range, uplo, ap, bp, vl, vu, il, iu, abstol, ldz)\n    or\n  NumRu::Lapack.dspgvx  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 12)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 12)", argc);
  rb_itype = argv[0];
  rb_jobz = argv[1];
  rb_range = argv[2];
  rb_uplo = argv[3];
  rb_ap = argv[4];
  rb_bp = argv[5];
  rb_vl = argv[6];
  rb_vu = argv[7];
  rb_il = argv[8];
  rb_iu = argv[9];
  rb_abstol = argv[10];
  rb_ldz = argv[11];

  abstol = NUM2DBL(rb_abstol);
  vl = NUM2DBL(rb_vl);
  iu = NUM2INT(rb_iu);
  jobz = StringValueCStr(rb_jobz)[0];
  vu = NUM2DBL(rb_vu);
  range = StringValueCStr(rb_range)[0];
  il = NUM2INT(rb_il);
  if (!NA_IsNArray(rb_ap))
    rb_raise(rb_eArgError, "ap (5th argument) must be NArray");
  if (NA_RANK(rb_ap) != 1)
    rb_raise(rb_eArgError, "rank of ap (5th argument) must be %d", 1);
  ldap = NA_SHAPE0(rb_ap);
  if (NA_TYPE(rb_ap) != NA_DFLOAT)
    rb_ap = na_change_type(rb_ap, NA_DFLOAT);
  ap = NA_PTR_TYPE(rb_ap, doublereal*);
  itype = NUM2INT(rb_itype);
  uplo = StringValueCStr(rb_uplo)[0];
  n = ((int)sqrtf(ldap*8+1.0f)-1)/2;
  m = lsame_(&range,"A") ? n : lsame_(&range,"I") ? iu-il+1 : 0;
  ldz = lsame_(&jobz,"V") ? MAX(1,n) : 1;
  if (!NA_IsNArray(rb_bp))
    rb_raise(rb_eArgError, "bp (6th argument) must be NArray");
  if (NA_RANK(rb_bp) != 1)
    rb_raise(rb_eArgError, "rank of bp (6th argument) must be %d", 1);
  if (NA_SHAPE0(rb_bp) != (n*(n+1)/2))
    rb_raise(rb_eRuntimeError, "shape 0 of bp must be %d", n*(n+1)/2);
  if (NA_TYPE(rb_bp) != NA_DFLOAT)
    rb_bp = na_change_type(rb_bp, NA_DFLOAT);
  bp = NA_PTR_TYPE(rb_bp, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_w = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  w = NA_PTR_TYPE(rb_w, doublereal*);
  {
    int shape[2];
    shape[0] = lsame_(&jobz,"N") ? 0 : ldz;
    shape[1] = lsame_(&jobz,"N") ? 0 : MAX(1,m);
    rb_z = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  z = NA_PTR_TYPE(rb_z, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_ifail = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  ifail = NA_PTR_TYPE(rb_ifail, integer*);
  {
    int shape[1];
    shape[0] = ldap;
    rb_ap_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  ap_out__ = NA_PTR_TYPE(rb_ap_out__, doublereal*);
  MEMCPY(ap_out__, ap, doublereal, NA_TOTAL(rb_ap));
  rb_ap = rb_ap_out__;
  ap = ap_out__;
  {
    int shape[1];
    shape[0] = n*(n+1)/2;
    rb_bp_out__ = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  bp_out__ = NA_PTR_TYPE(rb_bp_out__, doublereal*);
  MEMCPY(bp_out__, bp, doublereal, NA_TOTAL(rb_bp));
  rb_bp = rb_bp_out__;
  bp = bp_out__;
  work = ALLOC_N(doublereal, (8*n));
  iwork = ALLOC_N(integer, (5*n));

  dspgvx_(&itype, &jobz, &range, &uplo, &n, ap, bp, &vl, &vu, &il, &iu, &abstol, &m, w, z, &ldz, work, iwork, ifail, &info);

  free(work);
  free(iwork);
  rb_m = INT2NUM(m);
  rb_info = INT2NUM(info);
  return rb_ary_new3(7, rb_m, rb_w, rb_z, rb_ifail, rb_info, rb_ap, rb_bp);
}

void
init_lapack_dspgvx(VALUE mLapack){
  rb_define_module_function(mLapack, "dspgvx", rb_dspgvx, -1);
}
