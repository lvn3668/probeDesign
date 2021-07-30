#include<stdio.h>
#include<stdlib.h>
#include<math.h>

// Author: Lalitha Viswanathan
// Org: Tata Consultancy Services
long cdnalength,testgenelength,sizeofDNAStringignoringjunkchars,qdna;
char *DNA,*CDNA,*QDNA,*TESTGENE;
float crossreactivity;
# define NUMARGS 7

int main(int argc,char* argv[])
{
	long sizeoffnafile,line,i,start,stop,j,length;
	float len;
	int counter,new_counter,sub_counter,flag,probe_counter,in_counter;
	char* location,*possible_probe,*test_probe;
	FILE *filehandlefnafile,*filepointer1stfile,*filepointer2ndfile,*filepointerresultsfile;
	long startpos,stoppos,start2=-1,stop2=-1;
	char stri[1000],*str,strn;
	char genefile[100],probefile[100],tfile[100],DNAChar;
	// A .ptt file is an NCBI Protein Table file, 
	// which is a tab delimited file containing a list of all the proteins for their genomes 
	// (ftp://ftp.ncbi.nih.gov/genomes/). 
	// It corresponds with the CDS annotations from the GenBank file and can be created by 
	// parsing the GenBank files and writing the appropriate output.

	// The columns are :
	// Location    Strand    Length    PID    Gene    Synonym    Code    COG    Product
	if(argc!=NUMARGS)
	{
		printf("\n probedesign <fna file of genome> <pttfile of genome> <length of probe required> <%%crossreactivity allowed> <%% to be left on either side> <file containing gene for which probe is designed>\n");
		exit(1);
	}
	strcpy(genefile,argv[1]);
	strcpy(probefile,argv[2]);
	length=atoi(argv[3]);
	crossreactivity=atof(argv[4]);
	len=atof(argv[5]);
	strcpy(tfile,argv[6]);
	filepointerresultsfile=fopen("probes.txt","w+");
	possible_probe=(char *)malloc(sizeof(char)*len);
	test_probe=(char *)malloc(sizeof(char)*len);
	if(!(filehandlefnafile=fopen(genefile,"r")))
	{
			printf("File %s does not exist\n",genefile);
			exit(1);
	}
	 fseek(filehandlefnafile,0L,SEEK_END);
         sizeoffnafile=ftell(filehandlefnafile);
         DNA=(char *)malloc(sizeof(char)*sizeoffnafile);
         CDNA=(char *)malloc(sizeof(char)*sizeoffnafile);
         QDNA=(char *)malloc(sizeof(char)*sizeoffnafile);
         TESTGENE=(char *)malloc(sizeof(char)*sizeoffnafile);
	 sizeofDNAStringignoringjunkchars=1;
	 rewind(filehandlefnafile);
	 fgets(stri,500,filehandlefnafile);
	 while(!feof(filehandlefnafile))
	 {
		 DNAChar=fgetc(filehandlefnafile);
		 // Ignore junk characters
		 if((isalpha(DNAChar)) && ((DNAChar=='A')||(DNAChar=='T')||(DNAChar=='G')||(DNAChar=='C')))
				 DNA[sizeofDNAStringignoringjunkchars++]=DNAChar;
	 }//end of while
	fclose(filehandlefnafile);
	printf("Size of the genome =%ld(bp)\n",sizeofDNAStringignoringjunkchars);
	//size is size of fna file	
	if(!(filehandlefnafile=fopen(tfile,"r")))
	{
			printf("File %s does not exist\n",tfile);
			exit(1);
	}
	 testgenelength=1;
	 rewind(filehandlefnafile);
	 fgets(stri,500,filehandlefnafile);
	 while(!feof(filehandlefnafile))
	 {
		 DNAChar=fgetc(filehandlefnafile);
		 if(isalpha(DNAChar))
			 if((DNAChar=='A')||(DNAChar=='T')||(DNAChar=='G')||(DNAChar=='C'))
				 TESTGENE[testgenelength++]=DNAChar;
	}//end of while
	fclose(filehandlefnafile);
	printf("Size of the test gene =%ld(bp)\n",testgenelength);
	//beginning of function for checking for probe
	//len tells how much to leave 
	for(counter=0;counter<testgenelength-((len/100)*testgenelength);counter++)
	{
		new_counter=counter+((len/100)*testgenelength);
		probe_counter=0;
		for(in_counter=new_counter;in_counter<new_counter+length;in_counter++)
			possible_probe[probe_counter++]=TESTGENE[in_counter];
		possible_probe[probe_counter]='\0';
		if(!(filehandlefnafile=fopen(probefile,"r")))
		{
			printf("File %s does not exist\n",probefile);
			exit(1);
		}
		sizeoffnafile=1;
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
		  cdnalength=0;qdna=0;
		  sscanf(stri,"%ld..%ld %c %*s %*s %*s %*s %*s %*s \n",&startpos,&stoppos,&strn);
			if(strn=='+')
			{
			flag=0;
			if (startpos<stoppos)
			{
			for(sizeoffnafile=startpos;sizeoffnafile<=stoppos-3;sizeoffnafile++)
				CDNA[cdnalength++]=DNA[sizeoffnafile];
			for(sub_counter=0;sub_counter<cdnalength;sub_counter++)
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
				for(sizeoffnafile=startpos+3;sizeoffnafile<=stoppos;sizeoffnafile++)
					QDNA[qdna++]=DNA[sizeoffnafile];
				for(sub_counter=0;sub_counter<cdnalength;sub_counter++)
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
	   		stop2=stoppos;
			probe_counter=0;
  		}//end of line >6
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
		// cross reactivity is defined as 2m - mm to uniqueness of cleavage site
		// to ensure that probe does not splice gene at multiple locations 
		// Defined here as 2matches - mismatches
		// 2*matches - mismatches
		score=matches*2 - mismatches*1;
		if((score/testgenelength)*100>crossreactivity)
				return 0;
		else
				return 1;
}

