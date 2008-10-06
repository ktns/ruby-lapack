#include <string.h>
#include <math.h>
#include "ruby.h"
#include "narray.h"
#include "f2c.h"
#include "clapack.h"

#define MAX(a,b) a > b ? a : b
#define MIN(a,b) a < b ? a : b
#define LG(n) (int)ceil(log((double)n)/log(2.0))

extern logical lsame_(char *ca, char *cb);
extern logical ilaenv_(integer *ispec, char *name__, char *opts, integer *n1, integer *n2, integer *n3, integer *n4, ftnlen name_len, ftnlen opts_len);
extern int cunmtr_(char *side, char *uplo, char *trans, integer *m, integer *n, complex *a, integer *lda, complex *tau, complex *c__, integer *ldc, complex *work, integer *lwork, integer *info);
extern int cunmrz_(char *side, char *trans, integer *m, integer *n, integer *k, integer *l, complex *a, integer *lda, complex *tau, complex *c__, integer *ldc, complex *work, integer *lwork, integer *info);
