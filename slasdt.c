#include "rb_lapack.h"

extern VOID slasdt_(integer *n, integer *lvl, integer *nd, integer *inode, integer *ndiml, integer *ndimr, integer *msub);

static VALUE
rb_slasdt(int argc, VALUE *argv, VALUE self){
  VALUE rb_n;
  integer n; 
  VALUE rb_msub;
  integer msub; 
  VALUE rb_lvl;
  integer lvl; 
  VALUE rb_nd;
  integer nd; 
  VALUE rb_inode;
  integer *inode; 
  VALUE rb_ndiml;
  integer *ndiml; 
  VALUE rb_ndimr;
  integer *ndimr; 


  if (argc == 0) {
    printf("%s\n", "USAGE:\n  lvl, nd, inode, ndiml, ndimr = NumRu::Lapack.slasdt( n, msub)\n    or\n  NumRu::Lapack.slasdt  # print help\n\n\nFORTRAN MANUAL\n\n");
    return Qnil;
  }
  if (argc != 2)
    rb_raise(rb_eArgError,"wrong number of arguments (%d for 2)", argc);
  rb_n = argv[0];
  rb_msub = argv[1];

  n = NUM2INT(rb_n);
  msub = NUM2INT(rb_msub);
  {
    int shape[1];
    shape[0] = MAX(1,n);
    rb_inode = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  inode = NA_PTR_TYPE(rb_inode, integer*);
  {
    int shape[1];
    shape[0] = MAX(1,n);
    rb_ndiml = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  ndiml = NA_PTR_TYPE(rb_ndiml, integer*);
  {
    int shape[1];
    shape[0] = MAX(1,n);
    rb_ndimr = na_make_object(NA_LINT, 1, shape, cNArray);
  }
  ndimr = NA_PTR_TYPE(rb_ndimr, integer*);

  slasdt_(&n, &lvl, &nd, inode, ndiml, ndimr, &msub);

  rb_lvl = INT2NUM(lvl);
  rb_nd = INT2NUM(nd);
  return rb_ary_new3(5, rb_lvl, rb_nd, rb_inode, rb_ndiml, rb_ndimr);
}

void
init_lapack_slasdt(VALUE mLapack){
  rb_define_module_function(mLapack, "slasdt", rb_slasdt, -1);
}
