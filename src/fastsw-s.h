#ifndef __FASTSW__
#define __FASTSW__
#include "limits.h"
#include "assert.h"
#include "string.h"
#include "common-s.h"
#include "global-s.h"
#include "math.h"
#include "filters-s.h"

Alignment Fast_SW(char *Ref,char *Read,char *Quality,int OFF,char Sign);
Alignment Fast_SWX(char *Ref,char *Read,int OFF);
#endif
