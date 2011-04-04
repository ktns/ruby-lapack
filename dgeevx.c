#include "rb_lapack.h"

extern VOID dgeevx_(char *balanc, char *jobvl, char *jobvr, char *sense, integer *n, doublereal *a, integer *lda, doublereal *wr, doublereal *wi, doublereal *vl, integer *ldvl, doublereal *vr, integer *ldvr, integer *ilo, integer *ihi, doublereal *scale, doublereal *abnrm, doublereal *rconde, doublereal *rcondv, doublereal *work, integer *lwork, integer *iwork, integer *info);

static VALUE
rb_dgeevx(int argc, VALUE *argv, VALUE self){
  VALUE rb_balanc;
  char balanc; 
  VALUE rb_jobvl;
  char jobvl; 
  VALUE rb_jobvr;
  char jobvr; 
  VALUE rb_sense;
  char sense; 
  VALUE rb_a;
  doublereal *a; 
  VALUE rb_lwork;
  integer lwork; 
  VALUE rb_wr;
  doublereal *wr; 
  VALUE rb_wi;
  doublereal *wi; 
  VALUE rb_vl;
  doublereal *vl; 
  VALUE rb_vr;
  doublereal *vr; 
  VALUE rb_ilo;
  integer ilo; 
  VALUE rb_ihi;
  integer ihi; 
  VALUE rb_scale;
  doublereal *scale; 
  VALUE rb_abnrm;
  doublereal abnrm; 
  VALUE rb_rconde;
  doublereal *rconde; 
  VALUE rb_rcondv;
  doublereal *rcondv; 
  VALUE rb_work;
  doublereal *work; 
  VALUE rb_info;
  integer info; 
  VALUE rb_a_out__;
  doublereal *a_out__;
  integer *iwork;

  integer lda;
  integer n;
  integer ldvl;
  integer ldvr;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  wr, wi, vl, vr, ilo, ihi, scale, abnrm, rconde, rcondv, work, info, a = NumRu::Lapack.dgeevx( balanc, jobvl, jobvr, sense, a, lwork)\n    or\n  NumRu::Lapack.dgeevx  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 6)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 6)", argc);
  rb_balanc = argv[0];
  rb_jobvl = argv[1];
  rb_jobvr = argv[2];
  rb_sense = argv[3];
  rb_a = argv[4];
  rb_lwork = argv[5];

  if (!NA_IsNArray(rb_a))
    rb_raise(rb_eArgError, "a (5th argument) must be NArray");
  if (NA_RANK(rb_a) != 2)
    rb_raise(rb_eArgError, "rank of a (5th argument) must be %d", 2);
  n = NA_SHAPE1(rb_a);
  lda = NA_SHAPE0(rb_a);
  if (NA_TYPE(rb_a) != NA_DFLOAT)
    rb_a = na_change_type(rb_a, NA_DFLOAT);
  a = NA_PTR_TYPE(rb_a, doublereal*);
  jobvr = StringValueCStr(rb_jobvr)[0];
  balanc = StringValueCStr(rb_balanc)[0];
  sense = StringValueCStr(rb_sense)[0];
  jobvl = StringValueCStr(rb_jobvl)[0];
  lwork = NUM2INT(rb_lwork);
  ldvl = lsame_(&jobvl,"V") ? n : 1;
  ldvr = lsame_(&jobvr,"V") ? n : 1;
  {
    int shape[1];
    shape[0] = n;
    rb_wr = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  wr = NA_PTR_TYPE(rb_wr, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_wi = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  wi = NA_PTR_TYPE(rb_wi, doublereal*);
  {
    int shape[2];
    shape[0] = ldvl;
    shape[1] = n;
    rb_vl = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  vl = NA_PTR_TYPE(rb_vl, doublereal*);
  {
    int shape[2];
    shape[0] = ldvr;
    shape[1] = n;
    rb_vr = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  vr = NA_PTR_TYPE(rb_vr, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_scale = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  scale = NA_PTR_TYPE(rb_scale, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_rconde = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  rconde = NA_PTR_TYPE(rb_rconde, doublereal*);
  {
    int shape[1];
    shape[0] = n;
    rb_rcondv = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  rcondv = NA_PTR_TYPE(rb_rcondv, doublereal*);
  {
    int shape[1];
    shape[0] = MAX(1,lwork);
    rb_work = na_make_object(NA_DFLOAT, 1, shape, cNArray);
  }
  work = NA_PTR_TYPE(rb_work, doublereal*);
  {
    int shape[2];
    shape[0] = lda;
    shape[1] = n;
    rb_a_out__ = na_make_object(NA_DFLOAT, 2, shape, cNArray);
  }
  a_out__ = NA_PTR_TYPE(rb_a_out__, doublereal*);
  MEMCPY(a_out__, a, doublereal, NA_TOTAL(rb_a));
  rb_a = rb_a_out__;
  a = a_out__;
  iwork = ALLOC_N(integer, ((lsame_(&sense,"N")||lsame_(&sense,"E")) ? 0 : 2*n-2));

  dgeevx_(&balanc, &jobvl, &jobvr, &sense, &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr, &ilo, &ihi, scale, &abnrm, rconde, rcondv, work, &lwork, iwork, &info);

  free(iwork);
  rb_ilo = INT2NUM(ilo);
  rb_ihi = INT2NUM(ihi);
  rb_abnrm = rb_float_new((double)abnrm);
  rb_info = INT2NUM(info);
  return rb_ary_new3(13, rb_wr, rb_wi, rb_vl, rb_vr, rb_ilo, rb_ihi, rb_scale, rb_abnrm, rb_rconde, rb_rcondv, rb_work, rb_info, rb_a);
}

void
init_lapack_dgeevx(VALUE mLapack){
  rb_define_module_function(mLapack, "dgeevx", rb_dgeevx, -1);
}
