#include "rb_lapack.h"

extern VOID dlasdt_(integer *n, integer *lvl, integer *nd, integer *inode, integer *ndiml, integer *ndimr, integer *msub);

static VALUE sHelp, sUsage;

static VALUE
rblapack_dlasdt(int argc, VALUE *argv, VALUE self){
  VALUE rblapack_n;
  integer n; 
  VALUE rblapack_msub;
  integer msub; 
  VALUE rblapack_lvl;
  integer lvl; 
  VALUE rblapack_nd;
  integer nd; 
  VALUE rblapack_inode;
  integer *inode; 
  VALUE rblapack_ndiml;
  integer *ndiml; 
  VALUE rblapack_ndimr;
  integer *ndimr; 


  VALUE rb_options;
  if (argc > 0 && TYPE(argv[argc-1]) == T_HASH) {
    argc--;
    rb_options = argv[argc];
    if (rb_hash_aref(rb_options, sHelp) == Qtrue) {
      printf("%s\n", "USAGE:\n  lvl, nd, inode, ndiml, ndimr = NumRu::Lapack.dlasdt( n, msub, [:usage => usage, :help => help])\n\n\nFORTRAN MANUAL\n      SUBROUTINE DLASDT( N, LVL, ND, INODE, NDIML, NDIMR, MSUB )\n\n*  Purpose\n*  =======\n*\n*  DLASDT creates a tree of subproblems for bidiagonal divide and\n*  conquer.\n*\n\n*  Arguments\n*  =========\n*\n*   N      (input) INTEGER\n*          On entry, the number of diagonal elements of the\n*          bidiagonal matrix.\n*\n*   LVL    (output) INTEGER\n*          On exit, the number of levels on the computation tree.\n*\n*   ND     (output) INTEGER\n*          On exit, the number of nodes on the tree.\n*\n*   INODE  (output) INTEGER array, dimension ( N )\n*          On exit, centers of subproblems.\n*\n*   NDIML  (output) INTEGER array, dimension ( N )\n*          On exit, row dimensions of left children.\n*\n*   NDIMR  (output) INTEGER array, dimension ( N )\n*          On exit, row dimensions of right children.\n*\n*   MSUB   (input) INTEGER\n*          On entry, the maximum row dimension each subproblem at the\n*          bottom of the tree can be of.\n*\n\n*  Further Details\n*  ===============\n*\n*  Based on contributions by\n*     Ming Gu and Huan Ren, Computer Science Division, University of\n*     California at Berkeley, USA\n*\n*  =====================================================================\n*\n\n");
      rb_exit(0);
    }
    if (rb_hash_aref(rb_options, sUsage) == Qtrue) {
      printf("%s\n", "USAGE:\n  lvl, nd, inode, ndiml, ndimr = NumRu::Lapack.dlasdt( n, msub, [:usage => usage, :help => help])\n");
      rb_exit(0);
    } 
  } else
    rb_options = Qnil;
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rblapack_n = argv[0];
  rblapack_msub = argv[1];
  if (rb_options != Qnil) {
  }

  n = NUM2INT(rblapack_n);
  msub = NUM2INT(rblapack_msub);
  {
    int shape[1];
    shape[0] = MAX(1,n);
    rblapack_inode = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  inode = NA_PTR_TYPE(rblapack_inode, integer*);
  {
    int shape[1];
    shape[0] = MAX(1,n);
    rblapack_ndiml = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  ndiml = NA_PTR_TYPE(rblapack_ndiml, integer*);
  {
    int shape[1];
    shape[0] = MAX(1,n);
    rblapack_ndimr = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  ndimr = NA_PTR_TYPE(rblapack_ndimr, integer*);

  dlasdt_(&n, &lvl, &nd, inode, ndiml, ndimr, &msub);

  rblapack_lvl = INT2NUM(lvl);
  rblapack_nd = INT2NUM(nd);
  return rb_ary_new3(5, rblapack_lvl, rblapack_nd, rblapack_inode, rblapack_ndiml, rblapack_ndimr);
}

void
init_lapack_dlasdt(VALUE mLapack){
  rb_define_module_function(mLapack, "dlasdt", rblapack_dlasdt, -1);
  sHelp = ID2SYM(rb_intern("help"));
  sUsage = ID2SYM(rb_intern("usage"));
}
