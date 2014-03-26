#ifndef __PRINT__
#define __PRINT__
#include "global.h"
#include "common.h"
#include "swroutines.h"
#include "filters.h"

extern const Alignment Default_Alignment;
bool Quality_Passed(Cigar_Info & C,int StringLength);
//void Cigar_Check_And_Print(Hit_Info &  H,BATREAD & Read,int StringLength,FILE* Single_File,READ & R,bool Dont_Check_Quality,int Quality_Score,Alignment A=Default_Alignment,int Clip_H=0,int Clip_T=0,char* CIG=NULL);
void Cigar_Check_And_Print(Hit_Info &  H,BATREAD & Read,int StringLength,Final_Hit & Single_File,READ & R,bool Dont_Check_Quality,int Quality_Score,Alignment A,int Clip_H,int Clip_T,char* CIG);
void Print_Sam(Final_Hit & Single_File,READ & R,Hit_Info & H,int StringLength,int Quality_Score,Alignment A,int TClip_H,int TClip_T,char* TCIG);
//void Print_Sam(FILE* Single_File,READ & R,Hit_Info & H,int StringLength,int Quality_Score,Alignment A=Default_Alignment,int TClip_H=0,int TClip_T=0,char* TCIG=NULL);
void Print_Sam(Final_Hit & Single_File,READ & R,Hit_Info & H,int StringLength,int Quality_Score,Alignment A,int TClip_H,int TClip_T,char* TCIG);
void Print_Mishit(READ & R,FILE* Mishit_File);
void Print_Blank_Line(FILE* Single_File,READ & R);
void Read2Bin(char* Dest,char* Source,int StringLength);
void Read2RevCBin(char* Dest,char* Source,int StringLength);
void Reverse_Tag(char* Dest,const READ & R,int StringLength);
void Print_Unmapped(FILE* Single_File,READ & R,bool Mate_Mapped,unsigned Paired,unsigned HT);
#endif
