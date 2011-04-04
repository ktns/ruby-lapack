#include "rb_lapack.h"

extern VOID slaeda_(integer *n, integer *tlvls, integer *curlvl, integer *curpbm, integer *prmptr, integer *perm, integer *givptr, integer *givcol, real *givnum, real *q, integer *qptr, real *z, real *ztemp, integer *info);

static VALUE
rb_slaeda(int argc, VALUE *argv, VALUE self){
  VALUE rb_tlvls;
  integer tlvls; 
  VALUE rb_curlvl;
  integer curlvl; 
  VALUE rb_curpbm;
  integer curpbm; 
  VALUE rb_prmptr;
  integer *prmptr; 
  VALUE rb_perm;
  integer *perm; 
  VALUE rb_givptr;
  integer *givptr; 
  VALUE rb_givcol;
  integer *givcol; 
  VALUE rb_givnum;
  real *givnum; 
  VALUE rb_q;
  real *q; 
  VALUE rb_qptr;
  integer *qptr; 
  VALUE rb_z;
  real *z; 
  VALUE rb_info;
  integer info; 
  real *ztemp;

  integer ldqptr;
  integer n;

  if (argc == 0) {
    printf("%s\n", "USAGE:\n  z, info = NumRu::Lapack.slaeda( tlvls, curlvl, curpbm, prmptr, perm, givptr, givcol, givnum, q, qptr)\n    or\n  NumRu::Lapack.slaeda  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 10)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 10)", argc);
  rb_tlvls = argv[0];
  rb_curlvl = argv[1];
  rb_curpbm = argv[2];
  rb_prmptr = argv[3];
  rb_perm = argv[4];
  rb_givptr = argv[5];
  rb_givcol = argv[6];
  rb_givnum = argv[7];
  rb_q = argv[8];
  rb_qptr = argv[9];

  curpbm = NUM2INT(rb_curpbm);
  if (!NA_IsNArray(rb_qptr))
    rb_raise(rb_eArgError, "qptr (10th argument) must be NArray");
  if (NA_RANK(rb_qptr) != 1)
    rb_raise(rb_eArgError, "rank of qptr (10th argument) must be %d", 1);
  ldqptr = NA_SHAPE0(rb_qptr);
  if (NA_TYPE(rb_qptr) != NA_LINT)
    rb_qptr = na_change_type(rb_qptr, NA_LINT);
  qptr = NA_PTR_TYPE(rb_qptr, integer*);
  tlvls = NUM2INT(rb_tlvls);
  curlvl = NUM2INT(rb_curlvl);
  n = ldqptr-2;
  if (!NA_IsNArray(rb_perm))
    rb_raise(rb_eArgError, "perm (5th argument) must be NArray");
  if (NA_RANK(rb_perm) != 1)
    rb_raise(rb_eArgError, "rank of perm (5th argument) must be %d", 1);
  if (NA_SHAPE0(rb_perm) != (n*LG(n)))
    rb_raise(rb_eRuntimeError, "shape 0 of perm must be %d", n*LG(n));
  if (NA_TYPE(rb_perm) != NA_LINT)
    rb_perm = na_change_type(rb_perm, NA_LINT);
  perm = NA_PTR_TYPE(rb_perm, integer*);
  if (!NA_IsNArray(rb_prmptr))
    rb_raise(rb_eArgError, "prmptr (4th argument) must be NArray");
  if (NA_RANK(rb_prmptr) != 1)
    rb_raise(rb_eArgError, "rank of prmptr (4th argument) must be %d", 1);
  if (NA_SHAPE0(rb_prmptr) != (n*LG(n)))
    rb_raise(rb_eRuntimeError, "shape 0 of prmptr must be %d", n*LG(n));
  if (NA_TYPE(rb_prmptr) != NA_LINT)
    rb_prmptr = na_change_type(rb_prmptr, NA_LINT);
  prmptr = NA_PTR_TYPE(rb_prmptr, integer*);
  if (!NA_IsNArray(rb_givptr))
    rb_raise(rb_eArgError, "givptr (6th argument) must be NArray");
  if (NA_RANK(rb_givptr) != 1)
    rb_raise(rb_eArgError, "rank of givptr (6th argument) must be %d", 1);
  if (NA_SHAPE0(rb_givptr) != (n*LG(n)))
    rb_raise(rb_eRuntimeError, "shape 0 of givptr must be %d", n*LG(n));
  if (NA_TYPE(rb_givptr) != NA_LINT)
    rb_givptr = na_change_type(rb_givptr, NA_LINT);
  givptr = NA_PTR_TYPE(rb_givptr, integer*);
  if (!NA_IsNArray(rb_q))
    rb_raise(rb_eArgError, "q (9th argument) must be NArray");
  if (NA_RANK(rb_q) != 1)
    rb_raise(rb_eArgError, "rank of q (9th argument) must be %d", 1);
  if (NA_SHAPE0(rb_q) != (pow(n,2)))
    rb_raise(rb_eRuntimeError, "shape 0 of q must be %d", pow(n,2));
  if (NA_TYPE(rb_q) != NA_SFLOAT)
    rb_q = na_change_type(rb_q, NA_SFLOAT);
  q = NA_PTR_TYPE(rb_q, real*);
  if (!NA_IsNArray(rb_givcol))
    rb_raise(rb_eArgError, "givcol (7th argument) must be NArray");
  if (NA_RANK(rb_givcol) != 2)
    rb_raise(rb_eArgError, "rank of givcol (7th argument) must be %d", 2);
  if (NA_SHAPE1(rb_givcol) != (n*LG(n)))
    rb_raise(rb_eRuntimeError, "shape 1 of givcol must be %d", n*LG(n));
  if (NA_SHAPE0(rb_givcol) != (2))
    rb_raise(rb_eRuntimeError, "shape 0 of givcol must be %d", 2);
  if (NA_TYPE(rb_givcol) != NA_LINT)
    rb_givcol = na_change_type(rb_givcol, NA_LINT);
  givcol = NA_PTR_TYPE(rb_givcol, integer*);
  if (!NA_IsNArray(rb_givnum))
    rb_raise(rb_eArgError, "givnum (8th argument) must be NArray");
  if (NA_RANK(rb_givnum) != 2)
    rb_raise(rb_eArgError, "rank of givnum (8th argument) must be %d", 2);
  if (NA_SHAPE1(rb_givnum) != (n*LG(n)))
    rb_raise(rb_eRuntimeError, "shape 1 of givnum must be %d", n*LG(n));
  if (NA_SHAPE0(rb_givnum) != (2))
    rb_raise(rb_eRuntimeError, "shape 0 of givnum must be %d", 2);
  if (NA_TYPE(rb_givnum) != NA_SFLOAT)
    rb_givnum = na_change_type(rb_givnum, NA_SFLOAT);
  givnum = NA_PTR_TYPE(rb_givnum, real*);
  {
    int shape[1];
    shape[0] = n;
    rb_z = na_make_object(NA_SFLOAT, 1, shape, cNArray);
  }
  z = NA_PTR_TYPE(rb_z, real*);
  ztemp = ALLOC_N(real, (n));

  slaeda_(&n, &tlvls, &curlvl, &curpbm, prmptr, perm, givptr, givcol, givnum, q, qptr, z, ztemp, &info);

  free(ztemp);
  rb_info = INT2NUM(info);
  return rb_ary_new3(2, rb_z, rb_info);
}

void
init_lapack_slaeda(VALUE mLapack){
  rb_define_module_function(mLapack, "slaeda", rb_slaeda, -1);
}
