#include "data.inc"

#define __INT__            integer
#define __REAL__           real*8
#define __COMPLEX__        complex*16
#define __C_STRING_S_      character(len=STRING_SIZE_S_)
#define __C_STRING_L_      character(len=STRING_SIZE_L_)
#define __C_STRING_XL_      character(len=STRING_SIZE_XL_)
#define __C_OUT_STRING_S_  __C_STRING_S_, intent(out) 
#define __C_OUT_STRING_L_  __C_STRING_L_, intent(out) 
#define __C_OUT_STRING_XL_  __C_STRING_XL_, intent(out) 
