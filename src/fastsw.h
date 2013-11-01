#ifndef __FASTSW__
#define __FASTSW__
#include "limits.h"
#include "assert.h"
#include "string.h"
#include "common.h"
#include "global.h"
#include "math.h"
#include "filters.h"

Alignment Fast_SW(char *Ref,char *Read,char *Quality,int OFF,char Sign);
Alignment Fast_SWX(char *Ref,char *Read,int OFF);
void Mismatch_Scan_With_Score(char* Ref,char* Read,char* Qual,int Read_Len,int Max_Mis,int Current_Mis,Alignment & A);
#endif
