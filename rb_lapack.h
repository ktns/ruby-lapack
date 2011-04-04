#include <string.h>
#include <math.h>
#include "ruby.h"
#include "narray.h"
#include "f2c_minimal.h"

#define MAX(a,b) (a > b ? a : b)
#define MIN(a,b) (a < b ? a : b)
#define LG(n) ((int)ceil(log((double)n)/log(2.0)))

extern logical lsame_(char *ca, char *cb);
