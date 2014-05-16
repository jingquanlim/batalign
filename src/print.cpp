#include "print.h"
#include <cstdio>
#include "stdlib.h"
#include "assert.h"
const int SOFTCLIPLENGTH=10;
extern int INDELGAP;
extern bool STACK_LOWQ;
int Calc_MapQ(Hit_Info & H,Alignment & A,int Clip_Count);
extern int Top_Penalty;
extern int JUMP;
extern int SCOREGAP;
extern int SW_THRESHOLD;
extern std::string RGID;

int Find_Cigar(char* Cigar,Hit_Info & H,char* Current_Tag,int StringLength,READ & R,int & Clip_H,int & Clip_T)
{
	s_align* Aln;
	Ann_Info Ann;
	char Org_String[StringLength+50];
	Cigar_Info Cig_Info;

	int Jump=0;if(H.Sign=='-') Jump= 0+JUMP;
	s_profile* p = ssw_init((int8_t*)Current_Tag, StringLength, mata, n, 1);
	Get_Bases(H.Loc+Jump,StringLength+INDELGAP,Org_String);
	Aln=mengyao_ssw_core(Org_String,StringLength, Current_Tag,StringLength+INDELGAP,0,0/*DP*/, p);
	if(Aln->score1 >= ACC_SCORE)
	{
		H.SW_Score=Aln->score1;H.SW_Sub_Opt_Score=Aln->score2;
		Clip_H=Aln->read_begin1;Clip_T=0;
		if(Aln->read_end1!=StringLength-1) Clip_T=StringLength-1-Aln->read_end1;
		ssw_cigar_processQ(Aln,Cig_Info,Org_String,Aln->ref_begin1,Current_Tag,Aln->read_begin1,StringLength,R.Quality,NULL,Clip_H,Clip_T);
		ssw_cigar_print(Aln,Cigar,Cig_Info,Org_String+Aln->ref_begin1,Current_Tag+Aln->read_begin1,Clip_H,Clip_T,StringLength);
		H.Score= -Cig_Info.Score;
		H.QScore=Cig_Info.QScore;
		H.BQScore=Cig_Info.BQScore;
		H.Mismatch=Cig_Info.Mis;
		H.Indel=Cig_Info.Indel_Count;
		//A.Sign=H.Sign;

		H.Loc=H.Loc+Jump+Aln->ref_begin1-1;
		align_destroy(Aln);
		init_destroy(p); 
		return (Cig_Info.Length+Clip_T);
	}
	align_destroy(Aln);
	init_destroy(p); 
	return 0;
}

bool Quality_Passed(Cigar_Info & C,int StringLength)
{
	if(C.Mis>5) return false;
	if(C.Length<StringLength)
	{
		if(StringLength-C.Length+C.Mis >5) return false;
		//if(StringLength-C.Length>10) return false;

	}
	return true;
}

void Reverse_Tag(char* Dest,const READ & R,int StringLength)
{
	for (int i=StringLength-1;i>=0;i--)
	{
		*Dest=Char_To_CharC[R.Tag_Copy[i]];Dest++;
	}
	*Dest=0;
}

void Read2Bin(char* Dest,char* Source,int StringLength)
{
	for (int i=0;i<StringLength;i++)
	{
		*Dest=Char_To_Code[Source[i]];Dest++;
	}
	*Dest=0;
}

void Read2RevCBin(char* Dest,char* Source,int StringLength)
{
	for (int i=StringLength-1;i>=0;i--)
	{
		*Dest=Char_To_CodeC[Source[i]];Dest++;
	}
	*Dest=0;
}

void Cigar_Check_And_Print(Hit_Info &  H,BATREAD & Read,int StringLength,Final_Hit & Printed_Hit,READ & R,bool Dont_Check_Quality,int Quality_Score,Alignment A,int Clip_H,int Clip_T,char* CIG)
{
	H.Loc=H.Org_Loc;
	Print_Sam(Printed_Hit,R,H,StringLength,Quality_Score,A,Clip_H,Clip_T,CIG);
}

void Print_Sam(Final_Hit & Printed_Hit,READ & R,Hit_Info & H,int StringLength,int Quality_Score,Alignment A,int TClip_H,int TClip_T,char* TCIG)
{
	int Flag=0;
	int Skip=0;//StringLength;
	char Printed_CIG[MAX_SIGLEN];
	char* Qual=R.Quality,Rev_Qual[MAXTAG];char *CIG;
	char* Tag=R.Tag_Copy,Rev_Tag[MAXTAG];
	assert(H.Loc!=UINT_MAX && A.Score >=0);
	//int Real_Len=0;
	int Clip_H=TClip_H,Clip_T=TClip_T;
	if(TCIG)
	{
		if(TCIG[0])
		{
			CIG=TCIG;
			if(A.Loc)
				H.SW_Score=A.SW_Score;
			H.SW_Sub_Opt_Score=0;
		}
		else/*Should be bad Cigar*/
		{
			CIG=Printed_CIG;
			TCIG=NULL;
		}
	}
	else
	{
		assert(false);
		CIG=Printed_CIG;
	}

	int Sub_Opt_Score=0;

	if(R.Real_Len>=StringLength)
	{
		R.Tag_Copy[R.Real_Len]=0;
		R.Quality[R.Real_Len]=0;
		char Real_String[R.Real_Len];
		if(H.Sign=='+')
		{
			if(!TCIG)
			{
				Read2Bin(Real_String,R.Tag_Copy,R.Real_Len);
				Skip=Find_Cigar(CIG,H,Real_String,R.Real_Len,R,Clip_H,Clip_T);
			}

		}
		else
		{
			Reverse_Quality(Rev_Qual,R,R.Real_Len);
			Reverse_Tag(Rev_Tag,R,R.Real_Len);
			for(int i=0;i<R.Real_Len;i++)
			{
				R.Tag_Copy[i]=Rev_Tag[i];
				R.Quality[i]=Rev_Qual[i];
			}
			if(!TCIG)
			{
				Read2Bin(Real_String,R.Tag_Copy,R.Real_Len);
				H.Loc-=(R.Real_Len-StringLength)+INDELGAP-1;
				Skip=Find_Cigar(CIG,H,Real_String,R.Real_Len,R,Clip_H,Clip_T);
			}
		}

		if(TCIG)
		{
			if(A.Loc)
			{
				H.Score= A.Score;
				H.QScore=A.QScore;
				H.BQScore=A.BQScore;
				H.Mismatch=A.Mismatch;
				H.Indel=A.Indel;
			}
		}
		else
			H.Score= -H.Score;
	}
	else
	{
		assert(false);
		sprintf(CIG,"%dM",StringLength);
		R.Tag_Copy[StringLength]=0;
		R.Quality[StringLength]=0;
	}

	if(Quality_Score)
	{
		Quality_Score=Calc_MapQ(H,A,Clip_H+Clip_T);
		if(Quality_Score==1)
		{
			//if(A.SW_Score < 290)
			if(A.SW_Score < SW_THRESHOLD)
				Quality_Score=0;
			if(R.NCount >std::max(NCOUNT,int(5*StringLength/100)))
				Quality_Score=0;
		}
	}

	if(!CIG[0])
	{
		sprintf(CIG,"%dX",R.Real_Len);
		Quality_Score=0;
		fprintf (stdout,"\nCigar Error:%s\n",R.Description);
	}

	if(H.Sign=='-') 
	{
		Flag=16; 
		Qual=Rev_Qual;
		Tag=Rev_Tag;
	}
	else
	{
		Flag=0;
	}
	if(Sub_Opt_Score!=INT_MAX)
	{
		Printed_Hit.Loc=H.Loc;
		Printed_Hit.Skip=Skip;
		Printed_Hit.Flag=Flag;
		Printed_Hit.Quality_Score=Quality_Score;
		Printed_Hit.CIG=std::string(CIG);
		Printed_Hit.Tag=std::string(Tag);
		Printed_Hit.Qual=std::string(Qual);
		Printed_Hit.Mismatch=H.Mismatch;
		Printed_Hit.Score=H.Score;
		Printed_Hit.Sub_Opt_Score=Sub_Opt_Score;///////
		Printed_Hit.QScore=H.QScore;
		Printed_Hit.SW_Score=H.SW_Score;
		Printed_Hit.SW_Sub_Opt_Score=H.SW_Sub_Opt_Score;
	}
	else
	{
		Printed_Hit.Loc=H.Loc;
		Printed_Hit.Flag=Flag;
		Printed_Hit.Skip=Skip;
		Printed_Hit.Quality_Score=Quality_Score;
		Printed_Hit.CIG=std::string(CIG);
		Printed_Hit.Tag=std::string(Tag);
		Printed_Hit.Qual=std::string(Qual);
		Printed_Hit.Mismatch=H.Mismatch;
		Printed_Hit.Score=H.Score;
		Printed_Hit.Sub_Opt_Score=INT_MAX;///////
		Printed_Hit.QScore=H.QScore;
		Printed_Hit.SW_Score=H.SW_Score;
		Printed_Hit.SW_Sub_Opt_Score=H.SW_Sub_Opt_Score;
	}
}

void Print_Mishit(READ & R,FILE* Mishit_File)
{
	if(PRINT_MISHIT) fprintf(Mishit_File,"@%s\n%s\n+\n%s\n",R.Description+1,R.Tag_Copy,R.Quality);	
}

void Print_Blank_Line(FILE* Single_File,READ & R)
{
	R.Tag_Copy[R.Real_Len]=R.Quality[R.Real_Len]=0;
	fprintf(Single_File,"%s\t4\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\n",R.Description+1,R.Tag_Copy,R.Quality);
}

int Calc_MapQ(Hit_Info & H,Alignment & A,int Clip_Count)
{
	int Quality_Score;
	const int QUAL_START=60,Q_GAP=10;
	//int BOPEN=gap_open,BEXT=gap_extension;// BOPEN=6,BEXT=3,MATCH_BONUS=0;//2;
	//int BOPEN=gap_open,BEXT=gap_extension;// BOPEN=6,BEXT=3,MATCH_BONUS=0;//2;
	/*if(H.Status==PAIRED_SW)
	{
		Quality_Score=5;
	}
	else*/ if(H.Status==SW_RECOVERED || H.Status==PAIRED_SW)
	{
		int Sub_Opt_Score=A.Sub_Opt_Score;
		int MapQ=H.BQScore;
		if(Clip_Count)//>30)
		{
			int Score_Add=std::min(Clip_Count,gap_open);;//std::min((Clip_Count)*MN,BOPEN+BEXT*(Clip_Count-1));
			MapQ-=Score_Add;
		}	
		if(A.Score<A.Sub_Opt_Score-Q_GAP/3 )
		{
			Quality_Score=std::max(1,QUAL_START/2+MapQ);//-std::min(Top_Penalty,QUAL_START/2));
			if(!STACK_LOWQ && Quality_Score==1)
				Quality_Score=1;
			else
				Quality_Score+=5;
		}
		else
		{
			if(STACK_LOWQ)
				Quality_Score=1;
			else
				Quality_Score=1;
		}
		//Quality_Score=0;
	}
	//Unique hits..
	else if(H.Status==UNIQUEHIT)
	{
		int Sub_Opt_Score=INT_MAX;
		int MapQ=H.BQScore;
		int Offset=0;
		if(Clip_Count)//>30)
		{
			int Score_Add=std::min((Clip_Count)*MN,BOPEN+BEXT*(Clip_Count-1));//std::min(Clip_Count,gap_open);//0;//
			MapQ-=Score_Add;
		}	
		Quality_Score=std::max(1,QUAL_START-Offset+MapQ-std::min(Top_Penalty,QUAL_START/3));
		if(Quality_Score>1) 
			Quality_Score+=5;
		//if(!STACK_LOWQ && Quality_Score==1)
			//Quality_Score=0;

	}
	else if(H.Status==SHARP_UNIQUEHIT)
	{
		int Sub_Opt_Score=INT_MAX;
		int MapQ=H.BQScore;
		int Offset=0;
		if(Clip_Count)//>30)
		{
			int Score_Add=std::min((Clip_Count)*MN,BOPEN+BEXT*(Clip_Count-1));//std::min(Clip_Count,gap_open);//
			MapQ-=Score_Add;
			//Offset=5;
		}	
		//if(H.Indel) MapQ+=BOPEN;
		Quality_Score=std::max(1,QUAL_START-Offset+MapQ-std::min(Top_Penalty,QUAL_START/3));
		if(Quality_Score>1) 
			Quality_Score+=5;
		//if(!STACK_LOWQ && Quality_Score==1)
		//	Quality_Score=0;
	}
	//Multiple hits..
	else 
	{
		int Sub_Opt_Score=H.Sub_Opt_Score;
		assert(H.Status==MULTI_HIT);
		assert(H.Sub_Opt_Score!=INT_MAX);
		int MapQ=H.BQScore;
		if(Clip_Count)//>30)
		{
			int Score_Add=std::min((Clip_Count)*MN,BOPEN+BEXT*(Clip_Count-1));//std::min(Clip_Count,gap_open);
			MapQ-=Score_Add;
		}	
		if(A.Score<A.Sub_Opt_Score-SCOREGAP)//Q_GAP/3 )
		{
			Quality_Score=std::max(1,QUAL_START/2+MapQ);//-std::min(Top_Penalty,QUAL_START/2));
			if(Quality_Score==1)
			{
				Quality_Score=std::max(1,QUAL_START/2+5+MapQ);//-std::min(Top_Penalty,QUAL_START/2));
				if(Quality_Score==1)
				{
					Quality_Score=std::max(1,QUAL_START/2+5+5+MapQ);//-std::min(Top_Penalty,QUAL_START/2));
				}

				if(Quality_Score>5) Quality_Score=5;
				assert(Quality_Score<=5);
				return Quality_Score;
			}
			else
			{
				Quality_Score+=5;
				return Quality_Score;
			}
			/*if(!STACK_LOWQ && Quality_Score==1)
				Quality_Score=0;
			else
				Quality_Score+=5;*/
		}
		else
		{
			if(STACK_LOWQ)
				Quality_Score=1;
			else
				Quality_Score=1;
		}
		//Quality_Score=0;
	}
	return Quality_Score;
}

/*int Calc_MapQX(Hit_Info & H,Alignment & A,int Clip_Count)
{
	int Quality_Score;
	const int QUAL_START=210,Q_GAP=10;
	if(H.Status==SW_RECOVERED)
	{
		int Sub_Opt_Score=A.Sub_Opt_Score;
		int MapQ= -H.Score;
		if(Clip_Count)
		{
			int Score_Add=std::min((Clip_Count)*MN,BOPEN+BEXT*(Clip_Count-1));
			MapQ-=Score_Add;
		}	
		//MapQ= MapQ/2; 
		if(A.Score<A.Sub_Opt_Score-Q_GAP/2 )
		{
			Quality_Score=std::max(1,QUAL_START/3+MapQ);
		}
		else if(A.Score<A.Sub_Opt_Score-Q_GAP )
		{
			Quality_Score=std::max(1,QUAL_START/3+MapQ);
		}
		else
		{
			Quality_Score=0;
		}
		Quality_Score=0;
	}
	//Unique hits..
	else if(H.Status==UNIQUEHIT)
	{
		int Sub_Opt_Score=INT_MAX;
		int MapQ= -H.Score+H.BQScore;
		if(Clip_Count)
		{
			int Score_Add=std::min((Clip_Count)*MN,BOPEN+BEXT*(Clip_Count-1));
			MapQ+=Score_Add;
		}	
		if(H.Indel==1) MapQ+=BOPEN;
		//MapQ= MapQ/2; 
		Quality_Score=std::max(1,QUAL_START+MapQ);

	}
	else if(H.Status==SHARP_UNIQUEHIT)
	{
		int Sub_Opt_Score=INT_MAX;
		int MapQ= -H.Score+H.BQScore;
		if(Clip_Count)
		{
			int Score_Add=std::min((Clip_Count)*MN,BOPEN+BEXT*(Clip_Count-1));
			MapQ+=Score_Add;
		}	
		if(H.Indel==1) MapQ+=BOPEN;
		//MapQ= MapQ/2; 
		Quality_Score=std::max(1,QUAL_START+MapQ);

	}
	//Multiple hits..
	else 
	{
		int Sub_Opt_Score=H.Sub_Opt_Score;
		assert(H.Status==MULTI_HIT);
		assert(H.Sub_Opt_Score!=INT_MAX);
		int MapQ= -H.Score;
		if(Clip_Count)
		{
			int Score_Add=std::min((Clip_Count)*MN,BOPEN+BEXT*(Clip_Count-1));
			MapQ-=Score_Add;
		}	
		//MapQ= MapQ/2; 
		if(H.Score<H.Sub_Opt_Score-Q_GAP/2 )
		{
			Quality_Score=std::max(1,QUAL_START/3+MapQ);
		}
		else if(A.Score<A.Sub_Opt_Score-Q_GAP )
		{
			Quality_Score=std::max(1,QUAL_START/2+MapQ);
		}
		else
		{
			Quality_Score=0;
		}
		Quality_Score=0;
	}
	return Quality_Score;
}*/

void Print_Unmapped(FILE* Single_File,READ & R,bool Mate_Mapped,unsigned Paired,unsigned HT,int Read_Len)
{
	char* Qual=R.Quality;
	char* Tag=R.Tag_Copy;
	unsigned Flag=4;
	Flag=((Flag|Paired)|HT);
	if(!Mate_Mapped)
	{
		Flag|=8;
	}
	R.Tag_Copy[Read_Len]=R.Quality[Read_Len]=0;
	fprintf(Single_File,"%s\t%u\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\tRG:Z:%s\n",R.Description+1,Flag,R.Tag_Copy,R.Quality,RGID.c_str());
}
