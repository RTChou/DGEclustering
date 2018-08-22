/* GenClust
 * Copyright (C) 2002,2004 
 *  Vito Di Gesù, Raffaele Giancarlo,
 *  Giosuè Lo Bosco, Alessandra Raimondi
 *  and Davide Scaturro
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 * 02111-1307, USA.
 */

#include "genclustlb.h"
#include "string.h"

//this structure is defined to store information about the description and the name of each item
typedef struct
{
 char *geneName;
 char *descrizione;
}geneDescr;

// comma2point convert the ',' (comma) to '.' (point) in S
void comma2point(char *S)
{
 int i=0;
 while (S[i]!='\0')
 {
  if (S[i]==',') S[i]='.';
  i++;
 }
}

// point2comma convert the  '.' (point) to ',' (comma) in S
void point2comma(char *S)
{
 int i=0;
 while (S[i]!='\0')
 {
  if (S[i]=='.')
  { 
   if (S[i+1]!='\0') S[i]=',';
    else 
	   {
	    S[i]='\0';
        break;	   
	   }
  }
  i++;
 }
}

int main(int argc, char* argv[])
{
 if ((argc <= 3))
 {
  printf("Error: not enough program arguments (<=3) \n");
  printf("Usage : %s <file_input> <ncluster> <ngenerations> <fileoutput> <random_init> <output_type>\n",argv[0]);
 }
 else 
 {
  //filename is the name of the input file (dataset) the standard of this file
  //is as in input_sample.txt in the directory input_file_formats
  char *filename = argv[1];
  //npart is the number of classes
  uint nparti = atoi(argv[2]);
  if ((nparti<=0)||(nparti>255))
  {
   printf("Error: wrong ncluster value (0 < ncluster < 255)\n");
   printf("Usage : %s <file_input> <ncluster> <ngenerations> <fileoutput> <random_init> <output_type>\n",argv[0]);
  }
  else
  {
   ubyte npart = nparti;
   //nstepmax is the maximum number of generations
   uint nstepmax = atoi(argv[3]);
   //fileout is the output file, if not specified it will be considered the file with name
   //out.txt
   char *fileout = "out.txt";
   //randominit=0 random initialization of items in the partition,
   //otherwise the initialitation is what is specified in the file out.tmp  
   ubyte randominit = 1;
   //outputype=1 means that the final solution is the one with minimum internal variance found so far,
   //otherwise is the one found in the last generation of algorithm  
   ubyte outputtype = 0;
   FILE *f;
   //if the number of argument is >4 then fileout = argv[4]
   if (argc > 4) fileout = argv[4]; 
   //if the number of argument is >5 then randominit = argv[5]
   if (argc > 5) randominit = atoi(argv[5]);
    //if the number of argument is >6 then the type of output is argv[6]
   if (argc > 6)  outputtype = atoi(argv[6]);

   //open the input file
   f=fopen(filename,"rt");
   //S is the string which correspons to a line in the file
   if (f==NULL)
    printf("Error: Can't open input file %s\n",filename);
   else
   {
    char *S=(char *) malloc(255*sizeof(char));
	uint nrig,nfeatures,flags;
	//read nrig which is the first string value in the first row of the input file
    fscanf(f,"%s",S);
    nrig = atoi(S);
    //read the second string value in the first line of the input file which is the
    //number of features (or the number of columns)
    fscanf(f,"%s",S);
    nfeatures = atoi(S);
    //read the descritpion of the data and store it in flags, this is the third string
    //value in the first row of the input file
    fscanf(f,"%s",S);
    flags = atoi(S); 
    //flags =  0 there is no description of the items, there are no name of the items
    //flags =  10 there are name of the items, there is not descrition of the items
    //flags =  11 there are name of the items, there is descrition of the items
    //flags = 1  there are not name of the items, there is descrition of the items          
    if (nrig==0)
    {
	 printf("Error: The input file has zero rows, check if in the first line of the input file the number of rows of the dataset has been specified\n");
     printf("Usage : %s <file_input> <ncluster> <ngenerations> <fileoutput> <random_init> <output_type>\n",argv[0]);
 	}
	else
    if (nfeatures==0)
    {     
	 printf("Error: The input file has zero features, check if in the first line of the input file the number of featrues has been specified\n");
       printf("Usage : %s <file_input> <ncluster> <ngenerations> <fileoutput> <random_init> <output_type>\n",argv[0]);
	}
	else
   if (nrig<npart)
    {
	 printf("Error: The number of cluster is greater than the items in the dataset\n");
      printf("Usage : %s <file_input> <ncluster> <ngenerations> <fileoutput> <random_init> <output_type>\n",argv[0]);
	}
	else
    {
	 uint j,i,k,kk;
	 int jj;
	 int m,c,mm,n;
	 char *buff;
	 double *values,max=0,min=0,a;
	 //neuronNames is a vector which stores the names of the items (if present)
     geneDescr *neuronNames = (geneDescr *) malloc (nrig*sizeof(geneDescr));
     //featuresName is a vector which stores the names of the features
     char **featuresName = (char **) malloc(nfeatures*sizeof(char *));
     //buffer is a vector which allows to read the entire line
     char *buffer = (char *) malloc(nfeatures*255*sizeof(char));
     //read all the line and store it in buffer
     fscanf(f,"%s",buffer);
     m = strlen(buffer)+1;
     while((c=getc(f)!='\n')&&(feof(f) == 0))
	 {
      fscanf(f,"%s",&buffer[m]);
      m = m+strlen(&buffer[m])+1;
	 }
     // we know we have nfeatures number of features, so scan the buffer from the end 
     // and put the descriptions in the feature vector featuresName
      mm = m-1 ;
     for (jj=nfeatures-1;jj>=0;jj--)
	 {
      mm--;
	  while (buffer[mm]!=0) mm--;
      featuresName[jj]=&buffer[mm+1];
	 }
 	 // values is the vector which stores all the numerical values of the items 
     values = (double *) malloc(nrig*nfeatures*sizeof(double));
     // read all the values 
     i=0;
     while(fscanf(f,"%s",S)!=EOF)
	 {
      //buff is a vector which allows to read the entire line  
      buff = (char *) malloc(nfeatures*255*sizeof(char));
      m=0;
      n=0;
      //copy S in the buffer  
      strcpy(buff,S);
      m=strlen(S)+1;
      //read all the line and stores it in buuffer
      while((c=getc(f)!='\n')&&(feof(f) == 0))
	  {
       fscanf(f,"%s",&buff[m]);
       m = m+strlen(&buff[m])+1;
       n++;
	  }
      // we know we have nfeatures number of features, so scan buff from the end 
      // and store the values in the vector value
      mm = m-1 ;
 	  for (jj=nfeatures-1;jj>=0;jj--)
	  {
       mm--;
	   while (buff[mm]!=0) mm--;
       S=(char *) malloc(255*sizeof(char));;
	   strcpy(S,&buff[mm+1]);
       comma2point(S);
	   a = atof(S);
	   values[i*nfeatures+jj]=a;
       //while storing values find maximum and minimum values
	   if (a>max) max=a;
	   if (a<min) min=a; 
	  }
      //flags =  0 there is no description of the items, there are no name of the items
      //flags =  10 there are name of the items, there is not descrition of the items
      //flags =  11 there are name of the items, there is descrition of the items
      //flags = 1  there are not name of the items, there is descrition of the items          
      switch (flags)
	  {
       case 0: neuronNames[i].descrizione=0;
		    neuronNames[i].geneName=0;
            break;
       case 10: 
		     neuronNames[i].geneName = (char *) malloc(255*sizeof(char));
		     strcpy(neuronNames[i].geneName,&buff[0]);     
		     neuronNames[i].descrizione=0;
            break;
       case 11:
		     neuronNames[i].geneName = (char *) malloc(255*sizeof(char));
		     strcpy(neuronNames[i].geneName,&buff[0]);     
		     m=strlen(&buff[0])+1;
             for(j=0;j<n-nfeatures;j++) 
			 {
			  neuronNames[i].descrizione = (char *) malloc(255*sizeof(char));
			  strcpy(neuronNames[i].descrizione,&buff[m]);
			  m = m+strlen(&buff[m])+1;
			 }

			 break;
	   case 1: neuronNames[i].geneName=0;
		    m=strlen(&buff[0])+1;
             for(j=0;j<n-nfeatures;j++) 
			 {
			 neuronNames[i].descrizione = (char *) malloc(255*sizeof(char));
			  strcpy(neuronNames[i].descrizione,&buff[m]);
			  m = m+strlen(&buff[m])+1;
			 }
	  }    
      i++; 
      free(buff);
	 }
     fclose(f);
     if (i!=nrig)
 	  printf("Error: The number of items in the dataset file does not corresponds to what specified in the first line of the input file, check the input file %s\n",filename);
     else
	 {
       //parameters of the genetic algorithm
      double pc = 0.9;
      double pm = 0.01;
	  double *var=NULL,*varbest=NULL;
      ubyte *pinit=NULL,*notempty;
	  uint *neuronFreq;
      ubyte *result,l;
      FILE *fo;
	  ubyte flag;
	  //if randominit = 0 read the file INITFILE
      if (randominit == 0)
	  {
    	FILE *ff=fopen(INITFILE,"rt");
		if ((ff != NULL))
		{
 		 pinit = (ubyte *) malloc(nrig*sizeof(ubyte));
	  	 fgets(S,255,f);
		 fgets(S,255,f);
         i=0;		 
         while(feof(ff) == 0)
		 {
   		  fscanf(f,"%s",S);
	      kk=atoi(S); 
       	  for (j=0;j<kk;j++)
		  {
           fscanf(f,"%s",S);
	       k=atoi(S)-1;
           pinit[k]=i;
		  }
		  i++;
		 }
		 if ((i-1)!=npart)
         {
		 printf("Warning : Error reading initialization file %s, not enough rows, performing random inizialization\n",INITFILE);
		 free(pinit);
		 pinit = NULL;
		 randominit = 1;
		 }
		}
        else
		{
		 printf("Warning : Couldn't find initialization file %s, performing random inizialization\n",INITFILE);
		 randominit = 1;
		}
	  }
	  //inizialization
      result = genclust(nrig, nfeatures, nstepmax, npart, values, pc, pm, pinit,outputtype);
      //output the results
      fo=fopen(fileout,"wt");
      //neuronFreq is a vector wich stores how many times an item is present. This is used to verify
      //that the output is a partition
      neuronFreq = (uint *) malloc (sizeof(uint)*nrig);
	  notempty = (ubyte *) malloc(sizeof(ubyte)*npart);
      for (i=0;i<nrig;i++)
       neuronFreq[i]=0;
	  for (i=0;i<npart;i++)
       notempty[i]=0;
      for (l=0;l<npart;l++)
	  {
       fprintf(fo,"CLUSTER %d\n",l+1);
       for (j=0;j<nrig;j++)
	   {
	    if (result[j]==l)
		{
		 notempty[l]=1;
		 fprintf(fo," %s ",neuronNames[j].geneName);
   	     neuronFreq[j] = neuronFreq[j]+1;
	     for (k=0;k<nfeatures;k++) 
		 {
	      //the genetic clustering algorithm normalize all the values, so it needs to
		  //restore all the values to their original
		  sprintf(S,"%f",values[j*nfeatures+k]);
	      point2comma(S);
          fprintf(fo,"%s ",S);
		 }   
         putc('\n',fo);
	    }
	   }
	  }
	  flag = TRUE;
	  for (i=0;i<nrig;i++)
	   flag = (flag) & (neuronFreq[i]==1);
      for (i=0;i<npart;i++)
	   flag = (flag) & (notempty[i]==1);
      if (flag)
       fprintf(fo,"%s","- END -\n");
      else 
       fprintf(fo,"%s","- END WITH ERRORS -\n");
   	  printf(", output written in %s \n",fileout);
	 }
	}
   }   
  }
 }
 return 0;
}

