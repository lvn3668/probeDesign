#include<stdio.h>
#include<stdlib.h>
#include<math.h>

// Author: Lalitha Viswanathan
// Org: Tata Consultancy Services
long cdna,tdna,size,qdna;
char *DNA,*CDNA,*QDNA,*TDNA;
float crossreactivity;

int main(int argc,char* argv[])
{
	long temp,line,i,start,stop,j,length;
	float len;
	int counter,new_counter,sub_counter,flag,probe_counter,in_counter;
	char* location,*possible_probe,*test_probe;
	FILE *filehandlefnafile,*filepointer1stfile,*filepointer2ndfile,*filepointerresultsfile;
	long strt,stp,strt2=-1,stp2=-1;
	char stri[1000],*str,strn;
	char gfile[100],pfile[100],tfile[100],ch;
	if(argc!=7)
	{
		printf("\n probedesign <fna file of genome> <pttfile of genome> <length of probe required> <%%crossreactivity allowed> <%% to be left on either side> <file containing gene for which probe is designed>\n");
		exit(1);
	}
	strcpy(gfile,argv[1]);
	strcpy(pfile,argv[2]);
	length=atoi(argv[3]);
	crossreactivity=atof(argv[4]);
	len=atof(argv[5]);
	strcpy(tfile,argv[6]);
	filepointerresultsfile=fopen("probes.txt","w+");
	possible_probe=(char *)malloc(sizeof(char)*len);
	test_probe=(char *)malloc(sizeof(char)*len);
	if(!(filehandlefnafile=fopen(gfile,"r")))
	{
			printf("File %s does not exist\n",gfile);
			exit(1);
	}
	 fseek(filehandlefnafile,0L,SEEK_END);
         temp=ftell(filehandlefnafile);
         DNA=(char *)malloc(sizeof(char)*temp);
         CDNA=(char *)malloc(sizeof(char)*temp);
         QDNA=(char *)malloc(sizeof(char)*temp);
         TDNA=(char *)malloc(sizeof(char)*temp);
	 size=1;
	 rewind(filehandlefnafile);
	 fgets(stri,500,filehandlefnafile);
	 while(!feof(filehandlefnafile))
	 {
		 ch=fgetc(filehandlefnafile);
		 if((isalpha(ch)) && ((ch=='A')||(ch=='T')||(ch=='G')||(ch=='C')))
				 DNA[size++]=ch;
	 }//end of while
	fclose(filehandlefnafile);
	printf("Size of the genome =%ld(bp)\n",size);
	//size is size of fna file	
	if(!(filehandlefnafile=fopen(tfile,"r")))
	{
			printf("File %s does not exist\n",tfile);
			exit(1);
	}
	 tdna=1;
	 rewind(filehandlefnafile);
	 fgets(stri,500,filehandlefnafile);
	 while(!feof(filehandlefnafile))
	 {
		 ch=fgetc(filehandlefnafile);
		 if(isalpha(ch))
			 if((ch=='A')||(ch=='T')||(ch=='G')||(ch=='C'))
				 TDNA[tdna++]=ch;
	}//end of while
	fclose(filehandlefnafile);
	printf("Size of the test gene =%ld(bp)\n",tdna);
	//beginning of function for checking for probe
	//len tells how much to leave 
	for(counter=0;counter<tdna-((len/100)*tdna);counter++)
	{
		new_counter=counter+((len/100)*tdna);
		probe_counter=0;
		for(in_counter=new_counter;in_counter<new_counter+length;in_counter++)
			possible_probe[probe_counter++]=TDNA[in_counter];
		possible_probe[probe_counter]='\0';
		if(!(filehandlefnafile=fopen(pfile,"r")))
		{
			printf("File %s does not exist\n",pfile);
			exit(1);
		}
		temp=1;
		line=1;
	 	while(1)
	 	{
		 if(feof(filehandlefnafile))
		 {
		  	 break;	
		 }
		 str=fgets(stri,1000,filehandlefnafile);
		 line++;
		 if(line>6)
		  {
		  cdna=0;qdna=0;
		  sscanf(stri,"%ld..%ld %c %*s %*s %*s %*s %*s %*s \n",&strt,&stp,&strn);
			if(strn=='+')
			{
			flag=0;
			if (strt<stp)
			{
			for(temp=strt;temp<=stp-3;temp++)
				CDNA[cdna++]=DNA[temp];
			for(sub_counter=0;sub_counter<cdna;sub_counter++)
			{
			probe_counter=0;
			for(in_counter=sub_counter;in_counter<sub_counter+length;in_counter++)
				test_probe[probe_counter++]=CDNA[in_counter];
			test_probe[probe_counter]='\0';
			//if theres a greater degree of match, then reject it
			if(!check_for_crossreactivity(possible_probe,test_probe))
			{
				flag=0;
				break;
			}
			else
				flag=1;
			}
			if(flag)
				fprintf(filepointerresultsfile,"%s \n",possible_probe);
			}}
			else if(strn=='-')
			{
				for(temp=strt+3;temp<=stp;temp++)
					QDNA[qdna++]=DNA[temp];
				for(sub_counter=0;sub_counter<cdna;sub_counter++)
				{
				probe_counter=0;
				for(in_counter=sub_counter;in_counter<sub_counter+length;in_counter++)
					test_probe[probe_counter++]=QDNA[in_counter];
					test_probe[probe_counter]='\0';
				//if theres a greater degree of match, then reject it
				if(!check_for_crossreactivity(possible_probe,test_probe))
				{
					flag=0;
					break;
				}
				else
				flag=1;
				}
				if(flag)
					fprintf(filepointerresultsfile,"%s \n",possible_probe);
			}
	   		stp2=stp;
			probe_counter=0;
  		}//end of temp>6
	}//end of while
	}
	free(DNA);
	free(CDNA);
	free(QDNA);
}
// check for cross reactivity 
int check_for_crossreactivity(char* sequence1, char* sequence2)
{
		int counter,matches=0,mismatches=0;
		float score;
		for(counter=0;counter<strlen(sequence1);counter++)
		{
			if(sequence1[counter]==sequence2[counter])
				matches++;
			else
				mismatches++;
		}	
		score=matches*2 - mismatches*1;
		if((score/tdna)*100>crossreactivity)
				return 0;
		else
				return 1;
}

