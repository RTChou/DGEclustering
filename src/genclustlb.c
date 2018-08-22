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

uint nrig,nfeatures;
partitiontype part;
uint *pn,*pb;
ubyte npart;
int nparti;
double *values;
double pc,pm; 
uint *ind2bin,*lab2bin,*bin2ind;
ubyte *bin2lab;

//random float numbers generator
double random_double()
{
 return ((double)rand())/RAND_MAX;
}

//random integer numbers generator
uint random_intera(uint n)
{
 double num;
 uint converted;
 if (n==0) num = n;
  else
   num = (n-1)*random_double();
 converted = (uint) num;
  
 return(converted);
}

/* Memory allocators*/

//Memory allocation for the partition
partitiontype create_partition(ubyte npartmax)
{
 partitiontype partition = (partitiontype) malloc(npartmax*sizeof(ptype));
 uint i,j;
 for (i=0;i<npartmax;i++)
 {
  partition[i].var = (double *) malloc(nfeatures*sizeof(double));
  partition[i].bar = (double *) malloc(nfeatures*sizeof(double));
  for (j=0;j<nfeatures;j++)
  {
   partition[i].var[j] = 0;
   partition[i].bar[j] = 0;
  }
  partition[i].np=0;
 }
 return(partition);
}

//Memory allocation for the population
uint *create_population(uint n)
{
 uint *p = (uint *) malloc(n*sizeof(uint));
 uint i;
 for (i=0;i<n;i++) p[i]=i;
 return(p);
}

/* Coding and decoding functions*/

// codify is the function which codes the chromosome
uint codify(uint i,uint n, ubyte ind, ubyte npartmax)
{
 //a generic chromosome is a 32 bit value, the first 24 codify the  
 //generic item index, the last 8 codify the label
 uint indice = ind2bin[i];
 uint l = lab2bin[ind]; 
 indice=(indice<<8)|l;
 return indice;
}

// label finds the label where a coded chromosome belongs to
ubyte label(uint code)
{
 uint c = 255;
 //get the last 8 ubyte from code
 uint parti=(code&c);
 //decode the value
 ubyte li = (ubyte) parti;
 ubyte l = bin2lab[li];
 return(l);
}

// indexes find the item position.
uint indexes(uint codice)
{
 //get the first 24 ubyte from code
 uint mask = (uint) 16777215 << 8;
 uint cromind=(codice&mask);  
 cromind=(cromind>>8);  
 return(bin2ind[cromind]);
}

// Update the centroids and the variances after adding an element 
void add_item(partitiontype partition, uint i, ubyte l)
{
 uint j;
 double xn,nelem,mp;
 partition[l].np=partition[l].np+1;
 nelem = (double) partition[l].np;
 //update the centroids and the variances
 for(j=0;j<nfeatures;j++)
 {
  //xn is the element to add to cluster l
  xn = values[i*nfeatures+j];
  mp = partition[l].bar[j];
  partition[l].bar[j]=partition[l].bar[j]*((nelem-1)/nelem)+(xn/nelem);
  partition[l].var[j]=((nelem-1)*(partition[l].var[j]+pow(mp,2))-nelem*(pow(partition[l].bar[j],2))+pow(xn,2))/nelem;
 }
}

// Update the centroids and the variances after removing an element 
void remove_item(partitiontype partition,uint i,ubyte lab)
{
 uint j;
 double xn,nelem,mp;
 //update the centroids and the variances
 nelem = (double) partition[lab].np;
 for(j=0;j<nfeatures;j++)
 {
  //xn is the element to remove from cluster lab
  xn = values[i*nfeatures+j];
  mp=partition[lab].bar[j];	 
  if (nelem-1 > 0)  
  {
   partition[lab].bar[j] = (partition[lab].bar[j]-(xn/nelem))*(nelem/(nelem-1));
   partition[lab].var[j] = (((nelem*partition[lab].var[j])+(nelem*(pow(mp,2)))-pow(xn,2))/(nelem-1))-pow(partition[lab].bar[j],2);
  }
  else 
  {	  
   partition[lab].bar[j]=0;
   partition[lab].var[j]=0;
  }
 }
 partition[lab].np=partition[lab].np-1;
}

//Create the partition associated to a population of chromosomes, 
partitiontype partitioning(uint *p,uint n, ubyte npartmax)
{
 //allocate partition
 partitiontype partition = create_partition(npartmax);	
 uint i;
 //put the items of the population in the parition data structure
 for (i=0;i<n;i++)
  add_item(partition,i,label(p[i]));
 return(partition);	
}

//Initialize the population of chromosomes 
uint *initialize_population(uint n,ubyte npartmax, ubyte *pinit)
{
 uint *p;
 uint i,j,k,rt,szc,nv,r;
 ubyte l;
 //memory allocation
 //assigned stores the items not yet assigned
 uint *assigned = (uint *) malloc (n*sizeof(uint));
 p=create_population(n);
 if (pinit == NULL)
 {
  //random generation of assignements
  for (i=0;i<n;i++)
   assigned[i]=i;
  //rt is the maximum number of items to assign for each cluster
  rt = (uint) ((double) n )/npartmax;
  //nv is the number of items not yet assigned
  nv = n;
  for (l=0;l<npartmax;l++)
  {
   //szc is the number of items to assign to cluster l
   szc = 1+random_intera(rt-1);
   for (j=0;j<szc;j++)
   {
    //r is a random in the range 0 nv-1
    r = random_intera(nv);
	//assign the label l to the item assigned[r]
    p[assigned[r]]=codify(assigned[r],n,l,npartmax);
    //update the list of assigned items
	for (k=r;k<nv-1;k++) assigned[k]=assigned[k+1];
	nv = nv-1;
   }
  }
  for (i=0;i<nv;i++)
  {
   l=(ubyte) random_intera(npartmax);
   p[assigned[i]]=codify(assigned[i],n,l,npartmax);
  }
  free(assigned);
 }
 else 
 {
  for (i=0;i<n;i++)
  { 
   //pinit is the vector with the initial value for each item
   l = pinit[i];
   //codify the chromosome
   p[i]=codify(i,n,l,npartmax);
  }
 }
 return p;
}

//signle point crossover of a crhomosome p[i] and another random cromosome, the offspring are 
//stored in population ps, pc is the crossover rate
void single_point_crossover(uint *p , uint *ps , uint i , double pc)
{
 ubyte cut;
 uint numero,n,h1,h2,g2;
 double c;
 uint c1,c2;
 c=random_double();
 //act the crossover if c<pc
 if (c<pc)
 {
  numero=1;
  cut=(ubyte)random_intera(32);
  n=(numero<<(32-cut))-1;
  h1=p[i];
  g2=random_intera(nrig-1);
  h2=p[g2];
  c1=(uint)((h1 & (~n))|(h2 & n));
  c2=(uint)((h2 & (~n))|(h1 & n));
  ps[indexes(c1)]=c1;
  ps[indexes(c2)]=c2;
 }
}

//bit mutation of a single chromosome p[i], pm is the mutation rate
void bit_mutation(uint * p , uint i , double pm)
{
 ubyte k;
 double m;
 unsigned short int n,numero=1;
 uint cm;
 numero=1;
 cm=p[i];
 for (k=0;k<32;k++)
 {
  m=random_double();
  if (m<pm)
  {
   n=(uint)(numero<<k);
   if (n & cm) cm&=(~n);
	else cm|=n;
  }
 }
 p[indexes(cm)]=cm;
}

//calculate the total variance of a partition
double total_variance(partitiontype partition,ubyte npartmax)
{
 ubyte ind,j;
 double tot=0;
 double result;
 for(ind=0;ind<npartmax;ind++)
 {
  result=0;
  //calculate the variance for each features 
  for (j=0;j<nfeatures;j++)
    result=result+partition[ind].var[j];
  tot = tot + result;
 }
 return(sqrt(tot));
}

//fitness calculation
double fitness(partitiontype part,ubyte lab,uint i)
{
 // i is the index item (on the dataset) 
 // values is the vector whith all the dataset values
 // this fitness is the euclidean distance between the item i and the cluster l 
 // suggested by the chromosome p[i]
 uint j;
 double a,result=0,max;
 for (j=0;j<nfeatures;j++)
 {
  a=pow((values[i*nfeatures+j]-part[lab].bar[j]),2);
  max = values[i*nfeatures+j];
  if (max < part[lab].bar[j]) max = part[lab].bar[j]; 
  if (max == 0) max = 1;
  a = a/(max*max);
  result = result + a;
 }
 result = sqrt(result/nfeatures);
 return(result);
}

//seletion procedure, it gets the best chromosomes of the population p1 and p2
void selection(uint * p1 ,uint * p2, partitiontype part1, uint n)
{
 uint i;
 for (i=0;i<n;i++)
 {
  //new label of item i (in p2)
  ubyte labn=label(p2[i]); 
  //old label of item i (in p1)
  ubyte labv=label(p1[i]);
  double e1,e2;
  //fitness of the item i with the old label
  e1 = fitness(part1,labv,i);
  if ((labn!=labv) && (part1[labv].np > 1))
  {
   //remove the item i from the cluster labv 
   remove_item(part1,i,labv);
   //add the item i to the cluster labn 
   add_item(part1,i,labn);
   //fitness of the item i with the new label
   e2 = fitness(part1,labn,i);
   if (e2<e1)
    p1[i]=p2[i];
   else
   {	  
    //remove the item i from the cluster labn 
    remove_item(part1,i,labn);
    //add the item i to the cluster labv 
    add_item(part1,i,labv);	  
   }
  }
 }
}

//the crossover procedure, 
uint * crossover(uint * p , uint n , ubyte npartmax,  double pc)
{
 uint * ps;
 uint i;
 //ps is the population after the crossover of all elements
 ps=create_population(n);
 //copy population p in ps
 for(i=0;i<n;i++)
  ps[i]=p[i];
 //do the single point crossover for each chromosome in p
 for (i=0;i<n;i++)
  single_point_crossover(p,ps,i,pc);
 return ps;
}

//the mutation procedure, it mutates group of bits in the chromosome
void mutation(uint * p , uint n , double pm)
{
 uint i;
 for (i=0;i<n;bit_mutation(p,i,pm),i++);
}

//the genetic process for a single generation, p is the population at generation t, it returns the
//new population
void generate_population(uint * p ,partitiontype part1,uint n , ubyte npartmax)
{
 uint *p_succ=crossover(p,n,npartmax,pc);
 mutation(p_succ,n,pm);
 selection(p ,p_succ,part1,n);
 free(p_succ);
}

// Set the starting population, 
void Initialize(ubyte *pinit, double *var, double *varbest)
{
 uint i,ind=0;
 ubyte lab = 0;
 uint N = 16777216;
 uint n = 256;
 uint f = (uint) floor(N/nrig);
 ubyte g = (ubyte) floor(n/npart);
 int r = N%nrig;	
 int s = n%npart;	
 // Create the functions to encode and decode chromosomes
 // ind2bin[i] convert the index i in the range 0..nrig-1 into a binary in the range 0..(2^24)-1
 // lab2bin[l] convert the label l in the range 0..npart-1 into a binary in the range 0..(2^8)-1
 ind2bin = (uint *) malloc((nrig+1)*sizeof(uint));
 lab2bin = (uint *) malloc((npart+1)*sizeof(uint));
 // bin2ind[b] convert the binary b in the range 0..(2^24)-1 into an index in the range 0..nrig-1
 // bin2lab[b] convert the binary b in the range 0..(2^8)-1 into a label l in the range 0..npart-1
 bin2ind = (uint *) malloc(N*sizeof(uint));
 bin2lab = (ubyte *) malloc(n*sizeof(ubyte));
 ind2bin[0]=0;
 lab2bin[0]=0;
 for (i=1;i<nrig;i++)
   ind2bin[i]=ind2bin[i-1]+f+((r--)>0);
 for (i=1;i<npart;i++)
   lab2bin[i]=lab2bin[i-1]+g+((s--)>0);
 ind2bin[nrig]=N;
 lab2bin[npart]=n;
 for (i=0;i<N;i++)
 {
  ind = ind + (i>=ind2bin[ind+1]);
  bin2ind[i]=ind;
 }
 for (i=0;i<n;i++)
 {
  lab = lab + (i>=lab2bin[lab+1]);
  bin2lab[i]=lab;
 }
 //Initialize population pn 
 pn=initialize_population(nrig, npart, pinit);
 //Initialize population pb
 pb=create_population(nrig);
 //Create the structure which stores the partition
 part=partitioning(pn,nrig,npart);
 //Copy pn in pb
 for (i=0;i<nrig;i++) pb[i]=pn[i];
 *var=total_variance(part,npart);
 *varbest=*var;
}

// This corresponds to what is done for each generation of the genetic algorithm
void Generation(double *var, double *varbest) 
{
 uint i;
 generate_population(pn ,part, nrig ,npart); 
 *(var+1)=total_variance(part,npart);
 if ((*(var+1))<(*varbest)) 
 {
  for (i=0;i<nrig;i++) pb[i]=pn[i];
  *(varbest+1)=*(var+1);
 }
 else *(varbest+1) = *varbest;
}

ubyte *genclust(uint nr, uint nfeat, uint nstepmax, uint nparti, double *vl, double pcc, double pmm, ubyte *pinit, ubyte bestout)
{
 ubyte *outc; 
 nrig = nr;
 nfeatures = nfeat;
 pc = pcc;
 pm = pmm;
 outc = NULL;
 // srand((uint) time (NULL));
 //npart is the number of classes
 if ((nparti<=0)||(nparti>255))
  printf("Error: Wrong number of cluster (ncluster) value (0 < ncluster < 255)\n");
 else
 {
  npart = nparti;
  //nstepmax is the maximum number of generation
  if (nrig==0)
   printf("Error: Zero number of rows in the dataset\n");
  else
  if (nfeatures==0)
   printf("Error: Zero number of features in the dataset\n");
  else
  if (nrig<npart)
   printf("Error: The number of cluster is greater than the items in the dataset\n");
  else
   {
    const char *R = (pinit == NULL) ? "Yes" : "No";
    const char *T = (bestout) ? "Minimum Variance" :"Last Generation" ;
    double max = 0, min= 0, f;
    uint i,j, niter=0;
  	char *p,ch;
	uint uv;
	int v;
    double vd,*var,*varbest;
	FILE *fv,*fvbest;
    char S[255]; 
	printf("Number of itmes in the dataset : %d\n",nrig);
    printf("Number of features in the dataset : %d\n",nfeatures);
    printf("Number of desired clusters : %d\n",npart);
    printf("Number of generations : %d\n",nstepmax);
    printf("Random Inizialization : %s\n",R);
    printf("Output type : %s\n",T);
    //searching maxmin
    for (i = 0;i<nrig;i++)
     for(j = 0;j<nfeatures;j++)
	 {
	  f=vl[i*nfeatures+j];    
	  if (f>max) max=f;
	  if (f<min) min=f;
	 }
    //var and varbest stores the values of variance and bestvariance at each generation
    var = (double *) malloc ((nstepmax+1)*sizeof(double));
    varbest = (double *) malloc ((nstepmax+1)*sizeof(double));
    values = (double *) malloc(nrig*nfeatures*sizeof(double));
	//normalization of the values
    for (i = 0;i<nrig;i++)
     for(j = 0;j<nfeatures;j++)
        values[i*nfeatures+j]=(vl[i*nfeatures+j]-min)/(max-min);
    //inizialization
    Initialize(pinit,var,varbest);
    p = (char *) malloc(6*sizeof(char));
    p[0]='|';p[1]='/';p[2]='-';p[3]='|';p[4]='-';p[5]='\\';
    ch = '%';
    printf("Running the algorithm : Completed ");
    //repeat for the number of generations niter 
    if (nstepmax >0)
	{
	 do 
	 {
      Generation(var+niter,varbest+niter);
      if (niter > 0)
	  {
	   for (i=0;i<4;i++) putchar(0x08);
	   if (v>=10) putchar(0x08);
	  }
      vd = niter;
	  vd = (vd/nstepmax)*100;
	  v=(int) vd;
	  sprintf(S,"%d",v);
	  printf(S);
	  putchar(ch);
      putchar(' ');
	  putchar(p[niter%6]);
      fflush(stdout);
      niter++;
	 }
     while(nstepmax>=niter);
    }
	free(p);
	putchar(0x08);
	putchar(0x08);
	outc = (ubyte *) malloc(nrig*sizeof(ubyte));
	for(i=0;i<nrig;i++)
    {
	 uv = (bestout) ? pb[i] : pn[i]; 
     outc[i]=label(uv);
	}
	fv=fopen("variance.txt","wt");
    fvbest=fopen("bestvariance.txt","wt");
	for (i=0;i<nstepmax;i++)
	  {
       fprintf(fv,"%f\n",var[i]);
       fprintf(fvbest,"%f\n",varbest[i]);
	  }
	fclose(fv);
	fclose(fvbest);
  }
 }   
 return outc;
}
