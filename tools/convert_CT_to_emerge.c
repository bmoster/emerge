///////////////////////////////////////////////////////////////////////////////////////////////////
// Emerge code tools - File convert_CT_to_emerge.c                                               //
// Converts a standard output file from consistent-trees (e.g. tree_0_0_0.dat) to the emerge     //
// tree format. Call as:                                                                         //
// convert_CT_to_emerge <input_file> <output_file>                                               //
///////////////////////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

///////////////////////////////////////////////////////////////////////////////////////////////////
// Enable LONGIDS if IDs are larger than 4e9
// #define LONGIDS
///////////////////////////////////////////////////////////////////////////////////////////////////


#define BLKLEN fwrite(&blklen, sizeof(int), 1, ofp);


///////////////////////////////////////////////////////////////////////////////////////////////////
// Check that all variables are at the right column!
///////////////////////////////////////////////////////////////////////////////////////////////////
// Variable     Column_Number
///////////////////////////////////////////////////////////////////////////////////////////////////
#define SCALE   0
#define HALOID  1
#define DESCID  3
#define NP      4
#define UPID    6
#define MVIR    10
#define RVIR    11
#define RS      12
#define MMP     14
#define POSX    17
#define POSY    18
#define POSZ    19
#define VELX    20
#define VELY    21
#define VELZ    22
#define LAMBDA  26
#define BFID    27

#ifndef LONGIDS
typedef unsigned int IDType;
#else
typedef unsigned long long IDType;
#endif

//Declarations
struct haloin
{
  unsigned short np;
  unsigned short mmp;
  
  long long haloid;
  long long descid;
  long long upid;
  long long bfid;

  float scale;
  float mvir;
  float rvir;
  float rs;
  float lambda;

  float  pos[3];
  float  vel[3];
};

struct haloout
{
  IDType haloid;
  IDType descid;
  IDType upid;
  
  unsigned short np;
  unsigned short mmp;
  
  float scale;
  float mvir;
  float rvir;
  float c;
  float lambda;
  
  float  pos[3];
  float  vel[3];
};

//Function declarations
int compare_scale(const void *a, const void *b);
void assign_value(char *value, struct haloin *h, int column);
void loadBar(int x, int n, int r, int w);


//Main code
int main(int argc, char **argv)
{

  int i,j,blklen;
  int itree, ihalo;
  int Ntrees,Nhalostot;
  int Nmaxhalo;
  int *Nhalos;
  
  IDType *TreeID;
  
  long long MaxID;
  
  FILE *ifp, *ofp;

  const int linesize = 1000;
  char line[linesize], symbol;
  char *infname, *outfname, *tmp;

  struct haloin *hread;
  struct haloout hout;

  // Read input and output files
  if (argc < 3)
  {
    printf("Enter input file name and output file name!\n");
    exit(0);
  }
  
  infname  = calloc(200,sizeof *infname);
  outfname = calloc(200,sizeof *outfname);

  sprintf(infname, "%s",argv[1]);
  sprintf(outfname,"%s",argv[2]);

  Ntrees=0;
  Nhalostot=0;
  Nmaxhalo=0;
  ihalo=0;
  
  MaxID = (1UL << (sizeof(unsigned int)*8)) - 1;
  
  // Check sizes of types on this machine
  if(sizeof(long long) != 8)
  {
    printf("Type `long long' is not 64 bit on this platform.\n");
    exit(0);
  }
  if(sizeof(int) != 4)
  {
    printf("Type `int' is not 32 bit on this platform.\n");
    exit(0);
  }
  if(sizeof(unsigned long long) != 8)
  {
    printf("Type `unsigned long long' is not 64 bit on this platform.\n");
    exit(0);
  }
  if(sizeof(unsigned int) != 4)
  {
    printf("Type `unsigned int' is not 32 bit on this platform.\n");
    exit(0);
  }
  if(sizeof(float) != 4)
  {
    printf("Type `float' is not 32 bit on this platform.\n");
    exit(0);
  }
  if(sizeof(double) != 8)
  {
    printf("Type `double' is not 64 bit on this platform.\n");
    exit(0);
  }

  tmp  = calloc(200,sizeof *tmp);
  sprintf(tmp, "%s to %s",infname,outfname);
  
  printf("\n**********************************************************************\n");
  printf("*                                                                    *\n");
  printf("* TreeConverter: Converts Consistent-Trees file into Emerge format   *\n");
  printf("*                                                                    *\n");
  printf("**********************************************************************\n");
  printf("*                                                                    *\n");
  printf("* Converting %-55s *\n",tmp);
  printf("*                                                                    *\n");
  printf("**********************************************************************\n");
  printf("*                                                                    *\n");
  printf("* Checking the number of haloes in the tree file...                  *\n");
  
  free(tmp);
  
  //Count all haloes in file
  ifp=fopen(infname,"r");
  
  //Skip Header
  i=0;
  while ((fgets(line,linesize,ifp))&&(line[0]=='#')) i++;
  
  //Read Number of trees
  sscanf(line,"%d\n",&Ntrees);
  
  //Allocate Nhalos
  if(!(Nhalos=calloc(Ntrees,sizeof *Nhalos)))
  {
    fprintf(stderr,"Failed to allocate memory for Nhalos.\n");
    exit(0);
  }

  //Allocate TreeID
  if(!(TreeID=calloc(Ntrees,sizeof *TreeID)))
  {
    fprintf(stderr,"Failed to allocate memory for TreeID.\n");
    exit(0);
  }
  
  //Go through file and count number of haloes per tree
  itree = ihalo = 0;
  while (fgets(line,linesize,ifp))
  {
    if (line[0]=='#'){
      if (ihalo>0) itree++;
#ifndef LONGIDS
      sscanf(line,"%*s %u\n",&TreeID[itree]);
#else
      sscanf(line,"%*s %llu\n",&TreeID[itree]);
#endif
      if (ihalo>Nmaxhalo) Nmaxhalo=ihalo;
      ihalo=0;
      loadBar(itree, Ntrees, 50, 59);
    }
    if (line[0]!='#')
    {
      Nhalostot++;
      ihalo++;
      Nhalos[itree] = ihalo;
    }
  }
  itree++;
  fclose(ifp);
  
  if (itree!=Ntrees){
    printf("* Error! The number of trees stated in the file is unequal to the actual number of trees!\n");
    exit(0);
  }

  printf("* ... done checking the number of haloes                             *\n");
  printf("*                                                                    *\n");
  printf("* The file contains %10d trees and %10d haloes!          *\n",Ntrees,Nhalostot);
  printf("*                                                                    *\n");
  printf("**********************************************************************\n");
  
  if(!(hread=malloc(Nmaxhalo*sizeof(struct haloin))))
  {
    fprintf(stderr,"Failed to allocate memory.\n");
    exit(0);
  }

  ofp=fopen(outfname,"w");
  
  //Write Header containing the number of trees...
  blklen=sizeof(int);
  BLKLEN;
  fwrite(&Ntrees, sizeof(int), 1, ofp);
  BLKLEN;
  
  //...the number of haloes for each tree...
  blklen=Ntrees*sizeof(int);
  BLKLEN;
  fwrite(Nhalos, sizeof(int), Ntrees, ofp);
  BLKLEN;

  //...and the Ids of the trees
  blklen=Ntrees*sizeof(int);
  BLKLEN;
  fwrite(TreeID, sizeof(IDType), Ntrees, ofp);
  BLKLEN;
  
  printf("*                                                                    *\n");
  printf("* Loading haloes...                                                  *\n");

  itree=0;

  ifp=fopen(infname,"r");
  free(infname);
  i=0;
  while ((fgets(line,linesize,ifp))&&(line[0]=='#')) i++;
  sscanf(line,"%d\n",&Ntrees);
  fgets(line,linesize,ifp);
  
  //Write number of bytes for all haloes
  blklen=Nhalostot*sizeof(struct haloout);
  BLKLEN;

  // Loop over all trees in file
  while(itree<Ntrees)
  {
      i=0;
      while ((fgets(line,linesize,ifp))&&(line[0]!='#'))
      {
        //Parse line
        tmp = strtok(line, " ");
        assign_value(tmp, &hread[i], 0);
        j=1;
        while (tmp = strtok(NULL, " "))
        {
          assign_value(tmp, &hread[i], j);
          j++;
        }
        i++;
      }
    
    //Number of haloes in this tree
    Nhalos[itree]=i;
    
    //Sort tree in reverse breath first order
    qsort(hread, Nhalos[itree], sizeof(struct haloin), compare_scale);
    
    for (i=0;i<Nhalos[itree];i++)
    {
#ifndef LONGIDS
      if (hread[i].haloid > MaxID)
      {
        printf("* Halo at position %d has an ID of %lld.\n",i,hread[i].haloid);
        printf("* The maximum allowed value for integers is %lld.\n",MaxID);
        printf("* Use the option LONGIDS instead (also for EMERGE)!\n");
        exit(0);
      }
#endif
      //Assign values to output halo
      hout.np     = hread[i].np;
      hout.mmp    = hread[i].mmp;
      hout.scale  = hread[i].scale;
      hout.mvir   = hread[i].mvir;
      hout.rvir   = hread[i].rvir;
      hout.c      = hread[i].rvir/hread[i].rs;
      hout.lambda = hread[i].lambda;
      hout.pos[0] = hread[i].pos[0];
      hout.pos[1] = hread[i].pos[1];
      hout.pos[2] = hread[i].pos[2];
      hout.vel[0] = hread[i].vel[0];
      hout.vel[1] = hread[i].vel[1];
      hout.vel[2] = hread[i].vel[2];
      hout.haloid = (IDType) (hread[i].haloid + 1); //Convert to unsigned (need to avoid -1)
      hout.descid = (IDType) (hread[i].descid + 1);
      hout.upid   = (IDType) (hread[i].upid + 1);
    
      //Write halo to output file
      fwrite(&hout, sizeof(struct haloout), 1, ofp);
    }

    //Next tree
    itree++;

    loadBar(itree, Ntrees, 50, 59);
    
  }
  
  BLKLEN;
  
  printf("* ... done loading haloes                                            *\n");
  printf("*                                                                    *\n");
  printf("**********************************************************************\n");
  printf("*                                                                    *\n");
  printf("* Written to file %40s           *\n",outfname);
  printf("*                                                                    *\n");
  printf("**********************************************************************\n");

  free(TreeID);
  free(Nhalos);
  free(outfname);
  fclose(ifp);
  fclose(ofp);
  free(hread);
  
}


// Process has done i out of n rounds,
// and we want a bar of width w and resolution r.
void loadBar(int x, int n, int r, int w)
{
  if (n<r) r = n;
  
  // Only update r times.
  if (( x % (n/r) != 0 ) && (x!=0)) return;
  
  // Calculuate the ratio of complete-to-incomplete.
  float ratio = x/(float)n;
  int   c     = ratio * w;
  
  // Show the percentage complete.
  printf("* %3d%% [", (int)(ratio*100) );
  
  // Show the load bar.
  for (x=0; x<c; x++)
    printf("=");
  
  for (x=c; x<w; x++)
    printf(" ");
  
  // ANSI Control codes to go back to the
  // previous line and clear it.
  printf("] *\n\033[F\033[J");
}


//Asigns the string in a given column to the correct value
void assign_value(char *value, struct haloin *h, int column)
{
  
  if (column == SCALE)  h->scale  = strtof(value, NULL);
  if (column == HALOID) h->haloid = strtoll(value, NULL, 0);
  if (column == DESCID) h->descid = strtoll(value, NULL, 0);
  if (column == NP)     h->np     = (unsigned short) strtoul(value, NULL, 0);
  if (column == UPID)   h->upid   = strtoll(value, NULL, 0);
  if (column == MVIR)   h->mvir   = strtof(value, NULL);
  if (column == RVIR)   h->rvir   = strtof(value, NULL);
  if (column == RS)     h->rs     = strtof(value, NULL);
  if (column == MMP)    h->mmp    = (unsigned short) strtoul(value, NULL, 0);
  if (column == POSX)   h->pos[0] = strtof(value, NULL);
  if (column == POSY)   h->pos[1] = strtof(value, NULL);
  if (column == POSZ)   h->pos[2] = strtof(value, NULL);
  if (column == VELX)   h->vel[0] = strtof(value, NULL);
  if (column == VELY)   h->vel[1] = strtof(value, NULL);
  if (column == VELZ)   h->vel[2] = strtof(value, NULL);
  if (column == LAMBDA) h->lambda = strtof(value, NULL);
  if (column == BFID)   h->bfid   = strtoll(value, NULL, 0);
  
}


//Compare the scale factors of two haloes
//If they are the same compare the BFID of the haloes
//This is done to make the sorting stable
int compare_scale(const void *a, const void *b)
{
  const struct haloin *c = a;
  const struct haloin *d = b;
  int scale = (c->scale > d->scale) - (c->scale < d->scale);
  return scale != 0 ? scale : d->bfid - c->bfid;
}



