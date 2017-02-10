/* FPG  - forward population genetic simulation  Jody Hey */


/******** INCLUDE FILES  *************** */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <time.h>
#include <ctype.h>
#include <signal.h>
#include <assert.h>
#ifdef  __MWERKS__
#include <console.h>
#include "sioux.h"
#endif


/******** DEBUG RELATED OPTIONS *************** */
/* debugging related definitions */
/*MVC++  defines _DEBUG  under debuggin, and NDEBUG  under Release version*/

/* can define these to add extra debugging routines,*/
/* to use the myassert function define USE_MYASSERT */
#ifdef _DEBUG
#define USE_MYASSERT 
#define MYDEBUG
#endif


#ifndef USE_MYASSERT
#define myassert(a)
#endif

/******** SIMPLE DEFINITIONS MACROS *************** */

/* small for debug     
MAXALLCHROMESEGMENTS is the length of all the chromosomes
MAXGAMSIZE is maximum number of gametes in population

#define MAXALLCHROMESEGMENTS 40
#define MAXGAMSIZE  200   
*/

/*in MVC++
if the memory model is very large, you need to make the stack larger
go to project->settings 
on the link page
in window select 'output'
in reserve window  enter '10000000' or even '30000000'
*important*  must do a complete rebuild for it to matter
This worked for the following values

#define MAXALLCHROMESEGMENTS 50000
#define MAXGAMSIZE  100
 */

/* MVC++ has an integer data type with 64 bits called _int64 
but have had some problems with it, avoid*/

/* if compiled under Microsoft Visual C++, then _MSC_VER  is defined */
/* to find info under MVC++ help for other compiler predefinitions search help for
'predefined macros' */
/* if compiled under Metrowerks  CodeWarrior  then __MWERKS__ is defined */
/* to find info under Metrowerks Codewarrior  search the Codewarrior Manuals 
(not the basic help)  for 'Predefined Symbols */

/* This section was meant to use __int64 when compiled under MVC++
however this actually slows things down a bit, so don't bother 
#ifdef _MSC_VER 
#define SEGMENT_BITS  64
#define  SEGMENT_TYPE  unsigned __int64
#define MAXALLCHROMESEGMENTS 500
#else 
#define SEGMENT_BITS  32
#define  SEGMENT_TYPE  unsigned long 
#define MAXALLCHROMESEGMENTS 1000 
#endif */

#define SEGMENT_BITS  32
#define  SEGMENT_TYPE  unsigned long 
#define MAXALLCHROMESEGMENTS 1000 

/*#define MAXALLCHROMESEGMENTS 1000 */
#define MAXGAMSIZE  1000

#define MAXNUMCHROMOSOMES MAXALLCHROMESEGMENTS

/* maximum number of distinct fitness values*/
#define MAXSGP MAXGAMSIZE /2 
#define MAXZSIZE MAXGAMSIZE/2 
#define MAXALLCHROMELENGTH MAXALLCHROMESEGMENTS * SEGMENT_BITS
#define MINBASESPERREC 10 /* minmum number of average base pairs (or segment bits) between recombinaton events */

#if (MAXGAMSIZE >= 500)
#define MAXBINS  50
#elif (MAXGAMSIZE <=20)
#define MAXBINS 2;
#else 
#define MAXBINS  MAXGAMSIZE/10
#endif

#define OPTION_CHAR  '-'   /* precedes options on command line */
#define FP fprintf(rfile,
#define FPLINE fprintf(rfile,"-------------------------------------------------------------------------\n")
#define MAXPOPS  10
#define MAXMIGRANTS 100

/*
#define NOREC - sets recombination rate to zero and skips a bunch of steps
this can be useful if you want a special compiled version that does not do recombination 
*/

#define LDREGIONS  5
#define _FNSIZE 50
#define Dfault  FLT_MAX
#define Dfault_check  (Dfault / 10.0)
#define MAXOUTOFROOM 10
#define MAXOPTIONS   25
#define MYFLT_MAX   1e10
#define MYFLT_MIN   1e-10

/********   FUNCTIONAL MACROS *************** */
/* these SET macros deal with values from 1 to SEGMENT_BITS, rather than 0 to SEGMENT_BITS-1  */
#define single(i)  ((SEGMENT_TYPE) 1 << (i-1))

/* true if i is in the set x */
#define element(i,x)   ((single(i) & (x) )> 0)

#define isempty(x)  ((x) == (SEGMENT_TYPE) 0)

#define round(X) (long) ((X)+0.5)

#define forall(i,x) \
   for ( (i) = 1; (i) <= SEGMENT_BITS; ++(i)) if(element((i),(x)))
  /* forall : permits an action to be applied to all members of a segment
        for example 
                int k;
                forall(k,z) printf("%d ",k);   */

/* need a uniform random number generator that will not return  0 or 1 
use ran1(idum) from numerical recipes */

/* for random positive long integers between 0 and val-1 */

#define randint(val)   (unsigned long) (ran1(idum) * (val))


/********  NEW TYPES *************** */
typedef enum bool {_false, _true} boolean;

struct genome {
       SEGMENT_TYPE b[MAXALLCHROMESEGMENTS];
       float fit;  /* fitness*/
};

/* POLYALL, measure LD between all polys; POLYPOS, measure LD only among favored alleles, POLYNEG - only among deleterious */
enum {POLYALL,POLYPOS,POLYNEG, POLYNEUT};

/* LDR = r, LDRS = r squared  LDD = D LDABS = |D|  LDDP = D prime */
enum {LDR,LDRS,LDD,LDABS, LDDP};

/*
SADD  - deploys additive fitness model
SMULTI - deploys multiplicative fitness model
SEPI  - deploys epistatic fitness model 
*/
enum fitnessscheme {SADD,SMULTI,SEPI} fitmod;

/********  GLOBAL VARIABLES *************** */

/*
Notes on indexing mutations.
Each chromosome has chromesegments* SEGMENT_BITS positions for a mutation.
They are numbered beginning with 1. 
Each mutation has a segment number and a position within that segment.
Segment numbers are numbered beginning with zero.
Positions within segments are numberd from 1 to SEGMENT_BITS
Thus mutation number 1 is in segment 0, at position 1.
Mutation number SEGMENT_BITS is in segment 0 at position SEGMENT_BITS.
mutation number SEGMENT_BITS + 1 is in segment 1 at position 0.
Here is the basic routine for returning the segment and the position in the segment of a mutation at rspot

seg = rspot / SEGMENT_BITS;
inseg = rspot % SEGMENT_BITS;
if ( inseg ==0 ){
     inseg = SEGMENT_BITS;
     seg--;
}

For a mutation in seg at position inseg,  the mutation spot is seg*SEGMENT_BITS + inseg. 

The mutation positions that are available (i.e. not segregating) are in mutarray
*/



/* pointers to sections of memory  allocated by calloc */
struct genome *gpheap, *zpheap;
void *tempheapptr;
struct genome *zptr_list[MAXZSIZE];
int  Wzpheap[MAXPOPS+1];
struct genome migrants[MAXZSIZE];

/* gpheap and zpheap are kept pointing to large sections of memory 
zptr_list in an array of pointers, each to a 'z' 
think of zptre_list[i] as a pointer to the pointer at position i
in order to get zptr_list[i] to point to a particular z
zptr_list[i] = &particular_z
*/

unsigned long gamsize,zsize;
unsigned long chromesegments,allchromesegments;    
unsigned long seed_for_ran1;   
float pop_u_s, pop_u_n, total_pop_u, total_u,s_prop;         
float popscoeff, scoeff, pscoeff, nscoeff, pscoeff_dom_adj,nscoeff_dom_adj;       
float del_prop;        
float poprecrate, recrate;      
float bp;
float epi;
int num_measurements;    
int countints; /* count of the number of times things were measured */
long timebetween; /* the time in between measurements in units of gamsize generations*/
long timestartup; /* the time to run before measurements in units of gamsize generations*/
int   num_chromosomes;  /* the number of separate chromosomes */
int chromosome_ends[MAXNUMCHROMOSOMES];
int isegment_bits = SEGMENT_BITS; 
float mrate;
int numpops, nowpop;
unsigned long Wgamsize,Wzsize;

/* some things that will mostly be constant */
/* startpopgens*gamsize is the number of generations between assessments of mean fitness during burn in phase */
/* measure_interval* gamsize is the number of generations between measurements, starting after the burn in phase has ended */
float startpopgens, measure_interval; 

int checkint; 
/* fitness values associated with each group */
float fitnesses[MAXSGP];

/* group sizes of fitness classes */
int countfits[MAXSGP]; 
/* numbers of offspring gametes */
int gpsizes[MAXSGP];

unsigned long maxcrecs;
SEGMENT_TYPE emptyseg, fullseg;
int toomanygps;
long numgen, nowgen;
char rname[_FNSIZE]="", dname[_FNSIZE];
unsigned long  chromelength,allchromelength, sgp; /*sgp is possible # of distinct fitness */
FILE  *rfile, *dfile;
char *p;
int datasetnum = 0;
struct genome emptygenome, fullgenome;
unsigned long genomesize;
unsigned int sqrt_allchromelength;
long *idum, idumval;  /* used by ran1 */

/* # of positions left for mutations to occur (n is neutral,s selected)*/
/* mutavail-1 is also the position in mutarray from which to draw the location 
of the next mutation*/
unsigned long mutavail, minmutavail, meanmutavail_counts, mutavail_getting_low;
unsigned long availablespace;
boolean mutroomok;
int countoutofroom;
/* mean number of positions available */
float  meanmutavail;
boolean  fitness_absolute, resetfitness_absolute;

/* arrays indicated which positions in the chromosomes have not had mutations  */
/* mutarray holds the positions of segregating mutations, lowest possible # is 1 */
unsigned long mutarray[MAXALLCHROMELENGTH]={0};

/* average number of fitness classes */
float mean_numfitgps;
/* maximum number of fitness classes */
int max_numfitgps ;


/* scoeff_array has the selection coefficient associated with each
segregating mutation - be careful of 0,1 indexing problems
Use this, so the program does not really have to keep track of what is in it. 
Values can be assigned selection coefficients, but they need not be reset to 0, if 
the corresponding mutation is lost or fixed.
*/
float scoeff_array[MAXALLCHROMESEGMENTS][SEGMENT_BITS];

/* measure how many mutations segregating in the population */
float dominance;

int begin_count[MAXSGP];
int numfitgps;
float popmeanfit, overallmeanfit;
int burnend;

/* binfithet contains the fitnesses and heterozygosities
of a groups of zygotes */
struct binfithet{
	float  fit;
	float  het[POLYNEUT+1];
};

/* fithet contains the fitness and the heterozygosity of a zygote, a
pointer to that zygote, and the zygote number before sorting.
an array of these will be sorted by fitness. for use in measurement */
struct fithet{
	float  fit;
	float  het[POLYNEUT+1];
	int   onum; 
	struct genome *zptr1, *zptr2;
};

float  meanhet[POLYNEUT+1], mean_s[POLYNEUT+1];
float a1, Wa1;
char polystr[POLYNEUT+1][12];
struct  binfithet absbins[MAXBINS];
/* fh is filled in measure, and used there and by measureLD */
struct fithet fh[MAXZSIZE];
char  polyruntime[6], LDruntime[6];


boolean DoLD;
int LDmeasure, polys;


boolean LDmeasurebool[LDDP+1];
boolean polysbool[POLYNEUT+1];
struct LDaccumulate {
	float  LD[LDDP+1][POLYNEUT+1];
	unsigned long n[LDDP+1][POLYNEUT+1];
};
struct genome polygenome[POLYNEUT+1];

struct LDaccumulate LDall;
struct LDaccumulate LD_by_regions[LDREGIONS+1];
struct LDaccumulate LD_by_regions_by_fitness[LDREGIONS+1][MAXBINS];
struct LDaccumulate  LD_mean_fitnesses[MAXBINS];


/* number of bins for fitness ranks */
int bins, Wbins; 
/* the rank of the least fit individual to go in the least fit bin
if # bins * # of individuals per bin   = all individuals, then this is zero. */
int lowbinbottom, Wlowbinbottom;

/* number of LD regions */
int regions;

float tajd_sum[POLYNEUT+1], tajd_by_fitness_sum[MAXBINS][POLYNEUT+1];
unsigned long tajd_n[POLYNEUT+1],tajd_by_fitness_n[MAXBINS][POLYNEUT+1];


/* print out a data set  there should be MAXBINS rows*/
/* each row should be  num_chromosomes * chromesegments*SEGMENT_BITS in length */
boolean print_pseudo_data;
char *pchars[MAXPOPS*MAXBINS];
int pcharlines;

unsigned long fixcountgen, fixcountstart[POLYNEUT+1];
unsigned long  totalmuts, mutsincelastcheckfix, mutcount[POLYNEUT +1], polycount[POLYNEUT +1], lostcount[POLYNEUT+1], fixcount[POLYNEUT+1];
float meanpolycount[POLYNEUT+1];
int numcheckfix;

enum {SEXCLUSIVE,SSHARED,SFIXED};
float Fst[POLYNEUT+1], Nm[POLYNEUT+1], Hw[POLYNEUT+1], Hb[POLYNEUT+1];
int numFst[POLYNEUT+1], numNm[POLYNEUT+1], numHw[POLYNEUT+1], numHb[POLYNEUT+1];
float countstypes[SFIXED+1][POLYNEUT+1];
float meanvarfit, covarfit, meanvarhet[POLYNEUT+1],covarhet[POLYNEUT+1];

/*debug struct LDaccumulate Dposneg[2],Dpposneg[2];
unsigned long ndp = 0, pdp = 0; */

/* debug mutation positions */
/*double meanmutspot = 0.0, sqrmeanmutspot = 0.0;
unsigned long nummutspot = 0; 
*/
#ifdef  MYDEBUG
double meansmuts = 0.0, mutscalled = 0.0, mutsadded=0.0;
int meansmutcalls = 0;
#endif


/*********** PROTOTYPES ****************************/

#ifdef USE_MYASSERT
void myassert(unsigned short isit);
#endif
void myFPEsig(int sig);
void mySEGVsig(int sig);
void establish_signal_handler(void);
void strdelete(char *s,int pos,int len);
int pickpoisson(float param);
#ifndef NOREC
void shuffle(unsigned long *front, unsigned long length);
#endif
void shell(unsigned *lptr, int length);
void shellfithet(struct fithet *lptr, int length);
int segpos(int mspot, int *inseg);
void add_inseg(struct genome *r, unsigned long rspot);
void remove_inseg(struct genome *r, unsigned long rspot);
int cardinality(SEGMENT_TYPE x);
int countdiffs(struct genome *zptr1,struct genome *zptr2,int polys);
float ran1(long *idum);
float normdev(float mean, float variance);
void start(void);
void checkfix(struct genome *cptr);
float assignfit(struct genome *zptr, int *numhom, int *numhet);
void fitgroups(void);
double probability_mode(long mode, long gametes, double freq_good);
long binom_distr(double freq_good, double randnum, long gametes);
void get_famsizes(void);
void pairs(void);
void add_scoeff(unsigned long mspot);
void produce_gametes(void);
void migrate(void);
void generation(void);
void getmeanfit(void);
void measure(void);
void measureW(void);
void measureLD(void);
void calcLD(unsigned int x[4], float *LD);
float dotajd(int n, int ss, float pi);
void process_options(void);
void process_filename(int argc, char *argv[]);
void setup(int argc, char *argv[]);
void write_results(char *p);
boolean mutroom(void);
void run(void);
/******** LITTLE UTILITY FUNCTIONS *****************/

#ifdef USE_MYASSERT
void myassert(unsigned short isit){
   if (!isit){
      fclose(rfile);
      exit(0);
   }
}

#endif

void myFPEsig(int sig){
   if (rfile) {
      fprintf(rfile,"\n floating point crash \n");
      fclose(rfile);
   }
   printf("floating point crash \n");
    exit(1);
}
void mySEGVsig(int sig){
   if (rfile) {
      fprintf(rfile,"\n memory crash \n");
      fclose(rfile);
   }
   printf("memory crash in %d\n",nowgen);
    exit(2);
}

void establish_signal_handler(void){
   signal(SIGFPE,myFPEsig);
   signal(SIGSEGV,mySEGVsig);
}   

/* Delete the substring of length "len" at index "pos" from "s".
   Delete less if out-of-range. */
void strdelete(char *s,int pos,int len)
{
    int slen;

    if (--pos < 0)
        return;
    slen = strlen(s) - pos;
    if (slen <= 0)
        return;
    s += pos;
    if (slen <= len) {
        *s = 0;
        return;
    }
    while ((*s = s[len])) s++;
}

/* gets a poisson random variable.  If param > minpp the normal distribution is used as an approximation*/
#define minpp  50.0
int pickpoisson(float param){
long double randnum, raised;
int   i;
if (param >= minpp){ 
	i = round(normdev(param,param)); 
	return(i);
}
else {
 raised = exp(-param);
   randnum =  ran1(idum) ;
  for (i=0; randnum > raised; i++){
        randnum *= ran1(idum) ;
        }
 return (i);
}
} /* pickpoisson */
#undef minpp

#ifndef NOREC
void shuffle(unsigned long *front, unsigned long length){
/* take an array of int with some length and produce a random shuffle */
/* the first position in the array that can be shuffled is at front
the last position is at front + (length-1)*/
/* note that front may not actually be the front of the array, it is just the first
position that can be involved in the shuffle*/
/* got this routine off the net, supposedly it produces a properly random order */
unsigned long i,j,t;
for (i=length-1 ; i>=1 ;i--){
   /*j  = rand() % (i+1); */
	j = randint(i+1);
   t = *(front + i);
   *(front+i) = *(front+j);
   *(front+j) = t;
}
} /* shuffle */
#endif

void shell(unsigned *lptr, int length){
  float aln = 1.442695022, tiny = 1.0e-5, t;
 static   int nn,m,lognb2,i,j,k,l;
  lognb2 = floor(log(length)*aln + tiny);
  m = length;
  for (nn=1; nn<= lognb2; nn++){
      m = m /2;
      k = length - m;
      for (j = 0; j <= k-1; j++){
          i = j;
          reloop:
          l = i+m;
          if (*(lptr+l)< *(lptr+i)){
             t = *(lptr+i);
             *(lptr+i) = *(lptr+l);
             *(lptr+l) = t;
             i = i - m;
             if (i>= 0) goto reloop;
             }
          }
      }
} /* shell */


void shellfithet(struct fithet *lptr, int length){
  float aln = 1.442695022, tiny = 1.0e-5;
  struct fithet t;
 static   int nn,m,lognb2,i,j,k,l;
  lognb2 = floor(log(length)*aln + tiny);
  m = length;
  for (nn=1; nn<= lognb2; nn++){
      m = m /2;
      k = length - m;
      for (j = 0; j <= k-1; j++){
          i = j;
          reloop:
          l = i+m;
          if ((lptr+l)->fit < (lptr+i)->fit){
             t = *(lptr+i);
             *(lptr+i) = *(lptr+l);
             *(lptr+l) = t;
             i = i - m;
             if (i>= 0) goto reloop;
             }
          }
      }
} /* shellfithet */

int segpos(int mspot, int *inseg){
	int seg;
	seg = mspot / SEGMENT_BITS;
	*inseg = mspot % SEGMENT_BITS;
	if ( *inseg ==0 ){
		*inseg = SEGMENT_BITS;
		seg--;
	}
	return(seg);
} /* segpos */

void add_inseg(struct genome *r, unsigned long rspot){
   int seg, inseg;
   SEGMENT_TYPE tempseg;

   seg = segpos(rspot, &inseg);
   tempseg = single(inseg);
   r->b[seg] |= tempseg;
}

void remove_inseg(struct genome *r, unsigned long rspot){
   int seg, inseg;
   SEGMENT_TYPE tempseg;

   seg = segpos(rspot, &inseg);
   tempseg = single(inseg);
   r->b[seg] &= ~tempseg;
}

int cardinality(SEGMENT_TYPE x){
        int count = 0;
        while (!(isempty(x))) {
                x ^= (x & -x); ++count;
                }
        return(count);
}  /* cardinality */

int countdiffs(struct genome *zptr1,struct genome *zptr2,int polys){
int j, count;
SEGMENT_TYPE am, af;
	count = 0;
	for (j=0;j <= allchromesegments-1;j++) { 
	  am = zptr1->b[j] & polygenome[polys].b[j];
	  af = zptr2->b[j] & polygenome[polys].b[j];
	  count += cardinality(am ^ af);
	}
	return(count);
} /* countdiffs */

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1e-6 /* 1.2e-38 is too small ! */
#define RNMX 0.999999 /*(1.0-EPS) */

float ran1(long *idum)
/*Minimal" random number generator of Park and Miller with Bays-Durham shuffe and added
safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of the endpoint
values). Call with idum a negative integer to initialize; thereafter, do not alter idum between
successive deviates in a sequence. RNMX should approximate the largest oating value that is
less than 1. */
{
int j;
long k;
static long iy=0;
static long iv[NTAB];
float temp;
if (*idum <= 0 || !iy){
	if (-(*idum) < 1) *idum=1;  
	else *idum = -(*idum);
	for (j=NTAB+7;j>=0;j--){
		k=(*idum)/IQ;
		*idum=IA*(*idum-k*IQ)-IR*k;
		if (*idum < 0) *idum += IM;
		if (j < NTAB) iv[j] = *idum;
		}
	iy=iv[0];
}
k=(*idum)/IQ; 
*idum=IA*(*idum-k*IQ)-IR*k; 
if (*idum < 0) *idum += IM;
j=iy/NDIV; 
iy=iv[j];
iv[j] = *idum;
if ((temp=AM*iy) > RNMX) return RNMX; 
else return temp;
}
#undef IA
#undef IM 
#undef AM 
#undef IQ 
#undef IR
#undef NTAB 
#undef NDIV 
#undef EPS
#undef RNMX 

float normdev(float mean, float variance)
{
        static int iset=0;
        static float gset;
        float fac,rsq,v1,v2;
		float rescale;

        if  (iset == 0) {
                do {
                        v1=2.0*ran1(idum)-1.0;
                        v2=2.0*ran1(idum)-1.0;
                        rsq=v1*v1+v2*v2;
                } while (rsq >= 1.0 || rsq == 0.0);
                fac=sqrt(-2.0*log(rsq)/rsq);
                gset=v1*fac;
                iset=1;
				rescale = v2*fac*sqrt(variance) + mean;
                return rescale;
        } else {
                iset=0;
				rescale = gset*sqrt(variance) + mean;
                return rescale;
        }
}


/********  MAJOR GENERATION-RELATED FUNCTIONS *************** */

void start(void){
SEGMENT_TYPE tempulong;
unsigned long  i,j;

/* seed random number generation */
	idumval = - seed_for_ran1;
	idum = &idumval;

/* initialize population */
	allchromesegments = num_chromosomes * chromesegments;
	chromelength = chromesegments * SEGMENT_BITS;
	for (j=0;j< num_chromosomes;j++) chromosome_ends[j] = ((j+1)*chromesegments)-1;
        Wgamsize = gamsize /numpops;
	if (Wgamsize % 2) Wgamsize++; /* add 1 if gamsize is odd */
	gamsize = Wgamsize * numpops;
	zsize = gamsize/2;
	Wzsize = Wgamsize/2;
	genomesize = sizeof(struct genome);
	tempheapptr = calloc(gamsize,genomesize);
	gpheap = tempheapptr;
	tempheapptr = calloc(gamsize,genomesize);
	zpheap = tempheapptr;
	for (j=0;j<=numpops;j++) Wzpheap[j] = j*Wgamsize;
        for (i=1,a1=0.0;i< gamsize;i++)    a1 += 1.0/(float) i;
        for (i=1,Wa1=0.0;i< Wgamsize;i++) Wa1 += 1.0/(float) i;
    /* set maxinum number of fitness classes to gamsize or MAXSGP */
	if (zsize < MAXSGP) sgp = zsize; else sgp = MAXSGP;
	allchromelength = SEGMENT_BITS*allchromesegments;
        sqrt_allchromelength = (unsigned int) sqrt((double) allchromelength);
	maxcrecs = allchromelength/ MINBASESPERREC;
	mutavail = allchromelength;
	mutavail_getting_low = 4.0*(pop_u_s + pop_u_n);
	minmutavail=allchromelength;
	meanmutavail=0.0;
	meanmutavail_counts = 0;
	emptyseg = 0;
	fullseg = 0;
	for ( i=0; i<= SEGMENT_BITS-1 ;i++ ){
		tempulong = ((SEGMENT_TYPE) 1) << i;
		fullseg |= tempulong;
	}
	for (i=0;i<= allchromesegments-1;i++){
		emptygenome.b[i] = emptyseg;
		fullgenome.b[i] = fullseg;
	}
	emptygenome.fit = 0.0;
	for (j=0; j <= gamsize - 1; j++) memcpy(zpheap+j,&emptygenome,genomesize);
	mean_numfitgps = 0.0;
	max_numfitgps = 0;
	countoutofroom = 0;
	availablespace = allchromelength;
	for (i=0;i<= allchromelength-1;i++) mutarray[i] =i+1;
#ifndef NOREC
	shuffle(&mutarray[0], allchromelength);
#endif
	for (i=0;i<= allchromesegments-1;i++) 
		for (j=0;j<=SEGMENT_BITS-1;j++) scoeff_array[i][j]=0.0;

/* initialize generation counting */
	countints = 0;
	numgen = timestartup + timebetween*(num_measurements-1);
	nowgen = 0;
	

/* set up selection parameters */
/* total_pop_u is the total mutation rate, selected and neutral, for the population */
/* s_prop is the fraction of all mutations that is selection */
/* del_prop is the fraction of selection mutations that are deleterious*/

	total_pop_u = pop_u_s + pop_u_n;
	total_u = total_pop_u/gamsize;
	s_prop = pop_u_s/ total_pop_u;
	scoeff = popscoeff/gamsize;
	if (fitmod==SADD){
		pscoeff = scoeff;
		nscoeff = -scoeff;
	}
	if (fitmod==SEPI){
		pscoeff = scoeff;
		nscoeff = -scoeff;
	}
	if (fitmod==SMULTI){
		pscoeff = 1.0+scoeff;
		nscoeff = 1.0/pscoeff;
		pscoeff_dom_adj = (1.0+dominance*scoeff)/(1.0 + scoeff);
		nscoeff_dom_adj = (1.0+ scoeff)/(1.0 + dominance*scoeff);
	}

/* initialize measurement accumulators */	
	if (gamsize >= 500) bins = 50;
	else {if (gamsize <= 20) bins = 2;
		else bins = gamsize/10;
	}
        if (Wgamsize >= 500) Wbins = 50;
	else {if (Wgamsize <= 20) Wbins = 2;
		else Wbins = Wgamsize/10;
	}
    lowbinbottom = gamsize - (bins * (gamsize/bins));
    Wlowbinbottom = Wgamsize - (Wbins * (Wgamsize/Wbins));
	if (numpops > 1) pcharlines = numpops * Wbins;
	else pcharlines = bins;
	regions = LDREGIONS;
	strcpy(polystr[POLYALL],"    ALL    ");
	strcpy(polystr[POLYNEG],"DELETERIOUS");
	strcpy(polystr[POLYPOS],"BENEFICIAL ");
	strcpy(polystr[POLYNEUT],"  NEUTRAL  ");
	for (polys = POLYALL; polys <= POLYNEUT; polys++) {
		polygenome[polys]=emptygenome;
		mutcount[polys] = 0;
		polycount[polys] = 0;
		lostcount[polys]  = 0;
		fixcount[polys] = 0;
		fixcountstart[polys] = 0;
		meanpolycount[polys] = 0;
	}
	totalmuts = 0;
	numcheckfix = 0;
	mutsincelastcheckfix = 0;
	fixcountgen = 0;
	
	polysbool[POLYALL] = _true;
	for (i=0;i<= bins-1;i++){
		absbins[i].fit = 0.0;
		for (polys = POLYALL; polys <= POLYNEUT; polys++) absbins[i].het[polys] = 0.0;
	}
	meanvarfit = 0.0;covarfit = 0.0;
	for (polys = POLYALL; polys <= POLYNEUT; polys++){
		meanvarhet[polys] = 0.0;
		covarhet[polys] = 0.0;
	    Fst[polys] =0.0; Nm[polys]=0.0;
		Hw[polys] = 0.0; Hb[polys] = 0.0;
	    numFst[polys] =0; numNm[polys] = 0;
		numHw[polys] = 0; numHb[polys] = 0;
		countstypes[SSHARED][polys] = 0;
		countstypes[SEXCLUSIVE][polys] = 0;
		countstypes[SFIXED][polys] = 0; 
		for (LDmeasure = LDR;LDmeasure <= LDDP; LDmeasure++) {
			LDall.n[LDmeasure][polys] = 0;
			LDall.LD[LDmeasure][polys] = 0.0;
			for (i=0;i< regions+1;i++){
				LD_by_regions[i].n[LDmeasure][polys] = 0;
				LD_by_regions[i].LD[LDmeasure][polys] = 0.0;
				for (j=0;j< bins;j++){
					LD_by_regions_by_fitness[i][j].n[LDmeasure][polys] = 0;
					LD_by_regions_by_fitness[i][j].LD[LDmeasure][polys] = 0.0;
				}
			}
			for (j=0;j< bins;j++){
				LD_mean_fitnesses[j].n[LDmeasure][polys] =0;
				LD_mean_fitnesses[j].LD[LDmeasure][polys] =0.0;
			}
		}
		tajd_sum[polys] = 0.0; tajd_n[polys] = 0;
		for (j=0;j<bins;j++){
			tajd_by_fitness_sum[j][polys]=0.0;
			tajd_by_fitness_n[j][polys]=0;
		}
		meanhet[polys]=0.0;
	        mean_s[polys] = 0.0;

	}
	
/* set up beginning of results file */	
    if ((rfile = fopen(rname,"a")) == NULL){
           printf("Error opening text file for writing\n"); exit(1);}
    FP"\n\n\n*************************************************************************************\n");
    FP"*    PARAMETER SET NUMBER   %d\n",datasetnum);
    FP"*    Parameter String: %s \n",p);
    FP"*************************************************************************************\n\n");
    if (polysbool[POLYPOS] && !((pop_u_s > 0) && (popscoeff > 0) && (del_prop != 1.0))){
	  FP"Parameters prevent BENEFICIAL mutations - analysis ignored \n");
	  polysbool[POLYPOS] = _false;
	}
    if (polysbool[POLYNEG] && !((pop_u_s > 0) && (popscoeff > 0) && (del_prop != 0.0))){
	  FP"Parameters prevent HARMFUL mutations - analysis ignored \n");
	  polysbool[POLYNEG] = _false;
	}
    if (polysbool[POLYNEUT] && !((pop_u_n > 0) || ((popscoeff = 0) && (pop_u_s > 0.0)))){
	  FP"Parameters Prevent NEUTRAL mutations - analysis ignored \n");
	  polysbool[POLYNEUT] = _false;
	}
	
	FP"\n   Reports of Generation Number and Population Mean Fitness \n");
	FP"   ----------------------------------------------------------------------------\n");
} /* start */


boolean mutroom(void){
float needu, urate,temp;
int  i;
boolean  c;
temp = 0.0;
for (i=1;i<=gamsize;i++) temp += 1.0/(float) i;
if (popscoeff == 0.0) urate = 2*(pop_u_s + pop_u_n); else urate = 2*pop_u_n;
needu = urate * temp * 1.1;
c = (needu < allchromelength); 
return(c);
}

void checkfix(struct genome *cptr){
struct genome notemptyc, fullc;
SEGMENT_TYPE  tempcomp;
unsigned long lost, fixed, poly,lostclass[POLYNEUT+1],fixclass[POLYNEUT+1],polyclass[POLYNEUT+1];
unsigned long spot, i, k,tempu;
static unsigned long lastpolycount[POLYNEUT+1] = {0}, lastpoly = 0;;

/*Deal with fixations and losses*/

if (resetfitness_absolute)  fitness_absolute = _false;
fixed = 0; lost = 0; poly = 0;
numcheckfix++;
fullc = fullgenome;
notemptyc = emptygenome;
for ( i=0;i<= gamsize -1 ;i++)
   for ( k=0;k<= allchromesegments-1 ;k++ ){
      fullc.b[k] &= cptr[i].b[k];
      notemptyc.b[k] |= cptr[i].b[k];
   }


/* count fixations, losses and polymorphisms in each class, 
   remove fixations from population 
   remove fixations and lost mutations from polygenome
   */
for (polys=POLYALL; polys <= POLYNEUT;polys++){
	   fixclass[polys]=0;
	   polyclass[polys] = 0;
	   lostclass[polys] = 0;
	   polycount[polys] = 0;
}

for ( k=0;k<=allchromesegments-1 ;k++ ){
	fixed += cardinality(fullc.b[k]);
	notemptyc.b[k] &= ~fullc.b[k];
	poly += cardinality(notemptyc.b[k]);
	if ( fullc.b[k] > 0 ) for ( i=0;i<= gamsize -1 ;i++) cptr[i].b[k] &= ~fullc.b[k];
	for (polys=POLYALL; polys <= POLYNEUT;polys++){
		tempu = cardinality(polygenome[polys].b[k] & fullc.b[k]);
		fixclass[polys]  += tempu;
		fixed += tempu;
		polygenome[polys].b[k] &= ~fullc.b[k];
		tempu = cardinality(polygenome[polys].b[k] & ~notemptyc.b[k]);
		lostclass[polys] += tempu;
		lost += tempu;
		polygenome[polys].b[k] &= notemptyc.b[k];
		tempu = cardinality(polygenome[polys].b[k]);
		polyclass[polys] += tempu;
		poly += tempu;
	}
}
for (polys=POLYALL; polys <= POLYNEUT;polys++) {
	fixcount[polys] += fixclass[polys];
	lostcount[polys] += lostclass[polys];
	lastpolycount[polys] = polycount[polys];
	polycount[polys] = polyclass[polys];
	meanpolycount[polys] += (float) polycount[polys];
    if (fixcountgen > 0) fixcountstart[polys] += fixclass[polys];
}

myassert (totalmuts == fixcount[POLYALL] + lostcount[POLYALL] + polycount[POLYALL]);
#ifdef  MYDEBUG
for (polys=POLYALL; polys <= POLYNEUT;polys++) {
  myassert (mutcount[polys] == fixcount[polys] + lostcount[polys] + polycount[polys]);
}
#endif

/* make unused sites available */
availablespace = 0;
for (k=0;k<= allchromesegments-1 ;k++ ){
    tempcomp = ~notemptyc.b[k];
    forall(i,tempcomp){
		spot = k*SEGMENT_BITS + i;
		mutarray[availablespace] = spot;
	  /* FOR DEBUGGING */ 
		myassert(mutarray[availablespace] > 0);
		scoeff_array[k][i-1] = 0.0;
		availablespace++;
	}
}
for (i=availablespace; i<= allchromelength-1; i++) mutarray[i] = 0;
mutavail = availablespace;
meanmutavail += mutavail;
meanmutavail_counts++;
mutsincelastcheckfix = 0;
/* now shuffle the mutation positions that have been put back into mutarray */
/* this shuffling randomizes the mutation positions that get used, but there will 
still be some tencency for mutation positions that were just used to be reused. 
I think this is OK*/
#ifndef NOREC
if (availablespace) shuffle(&mutarray[0], availablespace);
#endif
} /* checkfix*/

float assignfit(struct genome *zptr, int *numhom, int *numhet){
   /* determine the fitness of a z, also return the number of homozygous sites, and the number of heterozygous sites */
   /* for each mutation, look up the fitness in scoeff_array */
   /* if fitness coefficient of zero implies neutrality, whether fitness
      is multiplicative or additive*/
   double temp,tempfit,fit;
   int i,j, hom=0, het=0;
   int op,on,ep=0,en=0;
   short h;
   SEGMENT_TYPE am,af,ahomo,ahetero,a,b;

   if (fitness_absolute){
	   op = fixcount[POLYPOS];
 	   on = fixcount[POLYNEG];
	   switch (fitmod){
		   case SEPI   : fit = exp(scoeff * pow(op-on,epi));break;
		   case SMULTI : fit = pow(pscoeff,op) * pow(nscoeff,on); break;
		   case SADD   : fit = (op - on) * pscoeff; break;
	   }
	   if (fit > MYFLT_MAX || fit < MYFLT_MIN){
		   op = 0;
		   on = 0;
		   fit = 1.0;
		   resetfitness_absolute = _true;

	   }
   } else {
	   op = 0;
	   on = 0;
	   fit = 1.0;
   }
   if (fitmod==SEPI) fit = 1.0;
   for (j=0;j <= allchromesegments-1;j++) { 
/* the places where the two chromosome both share a mutant are in ahomo */
/* the places where just one has a mutant are in ahetero */
	  am = zptr->b[j];
	  af = (zptr+1)->b[j];
	  ahomo = am & af;
	  ahetero = (am ^ af);
	  for (h = 0;h<= 1; h++){
		  if (h) a = ahetero; else a = ahomo;
		  while ( a ){
		     b = a & -a;
			 a =  a ^ b;
			i=0;
			while ( b ){
				b >>=1;
				i++;
			}
			if (fitmod==SEPI){
				if (h){
					if (scoeff_array[j][i-1] != 0.0) {
						if (scoeff_array[j][i-1] > 0.0) ep++;
						else en++;
					}
				} else{
					if (scoeff_array[j][i-1] != 0.0) {
						if (scoeff_array[j][i-1] > 0.0) op++;
						else on++;
					}
				}
			}
            if (fitmod==SADD){
				if (h) {
						fit += dominance*scoeff_array[j][i-1];
						het++;
				}
				else {
					fit += scoeff_array[j][i-1];
					hom++;
				}
			}
			if (fitmod==SMULTI){
				if (h){
					if (scoeff_array[j][i-1] != 0.0) {
						if (scoeff_array[j][i-1] > 1.0)	fit *= pscoeff_dom_adj*scoeff_array[j][i-1]; 
							else fit *= nscoeff_dom_adj*scoeff_array[j][i-1]; 
						het++;
					}
				} else {
				  if (scoeff_array[j][i-1] != 0.0){
					fit *= scoeff_array[j][i-1]; 
					hom++;
				  }
			  }
			}
		}
	  }
    }
   if (fitmod==SEPI){
	   temp = (op - on)+ dominance*(ep - en);
	   if (temp < 0.0) tempfit = exp(-scoeff * pow(-temp,epi));
	   else tempfit = exp(scoeff * pow(temp,epi));
	   fit = tempfit;
   }
   if ( fit < 0 ) fit = 0;
   zptr->fit = fit;
   (zptr+1)->fit = fit;
   *numhom = hom;
   *numhet = het;
   return(fit);
} /*assignfit*/


void fitgroups(void){
  /* apply selection coefficients to each pair of chromosomes, and order them by groups */
 int i,j,hom, het;
 float fitemp;
 int group_places[MAXGAMSIZE]={0};
 int check_count[MAXSGP]={0}; /*use to debug */

 fitnesses[0] = 0.0;
 countfits[0] = 0.0;
 numfitgps = 0;
 for (i=0;i<= Wgamsize-1;i+= 2){
   fitemp =  assignfit(zpheap+Wzpheap[nowpop]+i,&hom,&het);
   j= -1;
   do { j++;}
   while (((zpheap+Wzpheap[nowpop]+i)->fit != fitnesses[j])&&(j< numfitgps));
   assert((j>-1)&&(j<=numfitgps));
   if (j==numfitgps){
      if (j==sgp){
 /* too many groups, put chromosome in with last group */
          toomanygps++;
          j--;
          countfits[j]++;
       } else {
         fitnesses[j] = (zpheap+Wzpheap[nowpop]+i)->fit;
         numfitgps++;
         countfits[j] = 1;
       }
      group_places[i] = j;
      }
    else {
         group_places[i]=j;
          countfits[j]++;
         }
   }
 begin_count[0] = 0;
 for (j = 1;j< numfitgps;j++) begin_count[j]  = begin_count[j-1] + countfits[j-1];
 for (j=numfitgps; j <= MAXSGP-1; j++) countfits[j] = 0;
 for (i=0;i<= numfitgps;i++) check_count[i] = 0;
/* at this point there are 'numfitgps'  distinct fitness classes
the fitenss values are in 'fitnesses' and the size of class j is in
'countfits[j]' */
 /* now mustt fill zptr_list, which is divided into numfitgps, but points
 to the zs  
     */
 for (i=0;i<= Wgamsize-1; i+=2){
         zptr_list[begin_count[group_places[i]]] = zpheap+Wzpheap[nowpop]+i;
	 begin_count[group_places[i]]++;
 }
 begin_count[0] = 0;
 for (j = 1;j< numfitgps;j++) begin_count[j]  = begin_count[j-1] + countfits[j-1];
 
 /*now put begin_counts back to beginnings of each group  */
 for (j = 1;j< numfitgps;j++) begin_count[j]  = begin_count[j-1] + countfits[j-1];
 mean_numfitgps += numfitgps;
if (numfitgps > max_numfitgps) max_numfitgps = numfitgps;
} /* fitgroups*/


double probability_mode(long mode, long gametes, double freq_good){
/* sqrt(2*pi) */
 #define sqrt2pi  2.50662827463
 #define e 2.718281828459
 double   temp1,temp2,temp3,temp4,temp5,temp6,temp7;
 double  rmode,rgametes,pmode;

 rmode = (double) mode;
 rgametes = (double) gametes;
 temp1 = ( rgametes + 1.0 ) / (( rmode + 1.0 )*( rgametes - rmode + 1.0 ));
 temp1 = sqrt( temp1 );
 temp1 = e * temp1 / sqrt2pi;
 temp2 = rgametes * log(( rgametes + 1.0 ) * ( 1.0 - freq_good ) /
          ( rgametes - rmode + 1.0 ));
 temp3 = rmode * log(( rgametes - rmode + 1.0 ) * freq_good /
          (( rmode + 1.0 ) * ( 1.0 - freq_good )));
 temp4 =( 1.0/(rgametes + 1.0) - 1.0/(rmode + 1.0) - 1.0/
           (rgametes - rmode + 1.0)) / 12.0;
 temp5 = 1.0/pow((rgametes + 1.0),3.0) -
          1.0/pow((rmode + 1.0),3.0) -
          1.0/pow((rgametes - rmode + 1.0),3.0);
 temp5 = temp5/360.0;
 temp6 = 1.0/pow((rgametes + 1.0),5.0) -
          1.0/pow((rmode + 1.0),5.0) -
          1.0/pow((rgametes - rmode + 1.0),5.0);
 temp6 = temp6/1260.0;
 temp7 = temp2 + temp3 + temp4 - temp5 + temp6;
 pmode = exp(log(temp1) + temp7);
 return(pmode);
 }

long binom_distr(double freq_good, double randnum, long gametes){
/*this is copied almost directly from Kemp, and I don't understand all of it*/
short int flip;
double   a,b;
long   j,x,mode;
double   freq_mode;

flip = 1;
if (freq_good > 0.5){
   freq_good = 1.0 - freq_good;
   flip = 0;
   }
mode = (long int) ((gametes + 1.0) * freq_good);
if (mode < 15) goto L10;
x = mode;
freq_mode = probability_mode(mode,gametes,freq_good);
randnum = randnum - freq_mode;
if (randnum <= 0.0) goto L20;
a = freq_mode;
b = freq_mode;
j = 1;
while (j <= mode) {
 a = a * (mode-j+1) * (1.0 - freq_good) /((gametes - mode + j) * freq_good);
 randnum = randnum-a;
 if (randnum <= 0.0){
         x = mode-j;
         goto L20;
         }
 b = b * (gametes-mode+1-j) * freq_good /((mode + j) * (1.0 - freq_good));
 randnum = randnum - b;
 if (randnum <= 0.0){
         x = mode + j;
         goto L20;
         }
 j = j + 1;
 }
if ((mode + mode) >= gametes) goto L20;
j = mode + mode + 1;
while (j <= gametes){
 b = b * (gametes + 1 - j) * freq_good /((1.0 - freq_good) * freq_good);
 randnum = randnum - b;
 if (randnum <= 0.0) {
         x = j;
         goto L20;
         }
 j = j + 1;
 }
goto L20;
L10: a = pow((1.0 - freq_good),gametes);
j = 0;
while (j <= gametes) {
 randnum = randnum - a;
 if (randnum <= 0.0){
         x = j;
         goto L20;
         }
 a = (a * (gametes - j) * freq_good) /((1.0 - freq_good) * (j + 1));
 j = j + 1;
 }
L20: if (flip == 0) x = gametes - x;
return(x);
}

void get_famsizes(void){
/* first determine the expected offspring # in each fitness class. This means
calculating the mean fitness (sum of products of freq * fitness). Then
for each class say i, find expected p_i = fitness_i*freq_i/mean fitness.

Design a loop that creates multinomial sampling out of binomial sampling:
 first for group 1, sample x_1 the number of offspring from this group.
 This is done with regular old binomial sampling using group 1 versus
 everything else.
 Then for group 2, the expected frequency is the value previously determined
 for group2 divided by 1-the expected value for group 1. Also the new
 population size for group 2 is the total population size minus the number
 of offspring already found for group 1.
In general for group i, get expected freq of group i, by dividing the
actual exprected freq by 1 - sum of all expected freqs of those groups
already done. Use a population size that is the real populatino size
minus the sum of the number of offspring found so far.

Probably a good idea to see that the largest group is the last one done */

int i;
int bigsize=0, biggp = 0, sizes[MAXSGP];
double meanfit, expfreqs[MAXSGP];
long sumtaken = 0, binompick;
float sumexpfreq = 0.0/*, checksum=0.0*/;
float randnum;

if (numfitgps > 1){
 for (i=0;i<numfitgps;i++){
	sizes[i] = countfits[i];
    if (sizes[i] > bigsize){
        biggp = i;
        bigsize = sizes[i];
        }
 }
 meanfit =0.0;
 for (i=0;i< numfitgps;i++)
    meanfit += sizes[i]*fitnesses[i]/Wzsize;
 for (i=0;i < numfitgps;i++){
    expfreqs[i] = (sizes[i]*fitnesses[i]/Wzsize)/meanfit;
 /*checksum += expfreqs[i];    */
}
 /* expfreqs contains the expected proportion of gametes to come out
 of each frequency class 
 gpsizes hasthe number of gametes to come out of each class */
 for (i=0;i< numfitgps ;i++)
 if (i != biggp){
  randnum = ran1(idum) ;
  binompick = binom_distr(expfreqs[i]/(1.0 - sumexpfreq),randnum,Wgamsize-sumtaken);
  gpsizes[i] = binompick;
  sumtaken += binompick;
  sumexpfreq += expfreqs[i];
  }
 }
if (sumtaken > Wgamsize){
	gpsizes[biggp] = 0;
	while (sumtaken > Wgamsize){
		/*i = rand() % numfitgps; */
		i = randint(numfitgps);
		if (gpsizes[i] > 0){
			gpsizes[i]--;
			sumtaken--;
		}
	}
}
else gpsizes[biggp] = Wgamsize-sumtaken;

} /*get_famsizes */


void pairs(void){
	/* shuffle the genomes in the gpheap */
struct genome gtempi, gtempj;
int i,j;
	for (i=Wgamsize-1; i> 0;i--){
		/*j  = rand() % (i+1); */
		j = randint(i);
		gtempi = *(gpheap+Wzpheap[nowpop]+i);
		gtempj = *(gpheap+Wzpheap[nowpop]+j);
		memcpy(gpheap+Wzpheap[nowpop]+i,&gtempj,genomesize);
		memcpy(gpheap+Wzpheap[nowpop]+j,&gtempi,genomesize);
	}
} /*pairs */


void add_scoeff(unsigned long mspot){
float temp;
int seg, inseg;
/* pick a random number and compare with s_prop to see if mutation is selected 
or neutral. If it is neutral do not add a selection coefficient to the array */
seg = segpos(mspot,&inseg);
add_inseg(&polygenome[POLYALL],mspot);
mutcount[POLYALL]++;
temp = ran1(idum);
if (temp <=  s_prop ){
  temp=ran1(idum);
  /* FOR DEBUGGING */ myassert(scoeff_array[seg][inseg-1] == 0.0);
   if (del_prop > temp){
      scoeff_array[seg][inseg-1] = nscoeff;
	  add_inseg(&polygenome[POLYNEG],mspot);
	  mutcount[POLYNEG]++;
/*         printf("n %7.5f  %7.5f    %7.5f   %7.5f    %7.5f\n",del_prop,temp,scoeff_array[seg][inseg-1]);  */

   }
   else {
      scoeff_array[seg][inseg-1] = pscoeff;
	  add_inseg(&polygenome[POLYPOS],mspot);
	  mutcount[POLYPOS]++;
/*         printf("p %7.5f  %7.5f    %7.5f   %7.5f    %7.5f\n",del_prop,temp,scoeff_array[seg][inseg-1]);  */

   }
} else {
	add_inseg(&polygenome[POLYNEUT],mspot);
	mutcount[POLYNEUT]++;
}
} /* add_scoeff */



void produce_gametes(void){
/* draw randomly 
	countfits has the numbers of parental pairs of each fitness class
	pgsizes has the number of offspring gametes from each
	randomly pick a pair
		generate an offspring w/ recombination
		add mutations
	repeat until gpheap is full */
unsigned  temprec;
struct genome  mask_genome, *m,*f,newgenome;
unsigned long newspot,recs,parent,maskspot,tospot,seg,smuts,mutspot; 
unsigned recspots[MAXALLCHROMELENGTH/MINBASESPERREC];
unsigned long i,j,k,c;
float randnum;
boolean  outofroom;
#ifdef MYDEBUG
#define MAXFAMSIZE  21
int  actualfamsizes[MAXGAMSIZE]={0};
int famsizedist[MAXFAMSIZE]={0};
#endif

for (i=0,newspot=0;i< numfitgps;i++)
	while (gpsizes[i]){
	  parent = randint(countfits[i]);  
		/* parent = rand() % countfits[i]; */
#ifdef MYDEBUG
       actualfamsizes[parent]++;
#endif
	  mask_genome = emptygenome;
	  randnum =  ran1(idum);
	  if (randnum < 0.5){
			m= zptr_list[begin_count[i] + parent];
			f= (zptr_list[begin_count[i] + parent]+1);
		} else{
			m= zptr_list[begin_count[i] + parent]+1;
			f= zptr_list[begin_count[i] + parent];
		}
	  for (c=0,j=0;c< num_chromosomes;c++){
		recs=pickpoisson(recrate);
		if (recs > maxcrecs)  recs=maxcrecs;
		k = 0;
		while ( k< recs ){
			if (c==0) temprec = randint(chromelength) + 1;
			else {
				temprec = ((chromosome_ends[c-1]+1)*SEGMENT_BITS)+randint(chromelength) + 1;
			}
			j=0;
			while (( j< k )&&(temprec != recspots[j])) j++;
			if (j == k ) {
				recspots[k] = temprec;
				k++;
			}
		}
		if (recs > 1) shell(&recspots[0],recs);
		if (recs){
			maskspot = ((chromosome_ends[c-1]+1)*SEGMENT_BITS)+1;
			for ( k=0;k<= recs ;k++){
				if ( k ==recs ) tospot = (chromosome_ends[c]+1)*SEGMENT_BITS+1;
					else tospot = recspots[k];
				for ( j=maskspot;j < tospot ; j++ )
					if ( (k%2) != 0 ) add_inseg(&mask_genome,j);
				maskspot = tospot;
			}
		}
	  }
	for (k=0;k < num_chromosomes;k++){
		randnum =  ran1(idum);
		if (randnum < 0.5){
			for (j = k*chromesegments;j<= chromosome_ends[k];j++) mask_genome.b[j] =  ~mask_genome.b[j];
		}
	}
	for ( seg=0;seg <= allchromesegments -1 ; seg++)
		newgenome.b[seg] = (m->b[seg] & mask_genome.b[seg]) | (f->b[seg] & ~(mask_genome.b[seg]));
	memcpy(gpheap+Wzpheap[nowpop]+newspot,&newgenome,genomesize);
	gpsizes[i]--;
	newspot++;
	}
#ifdef MYDEBUG
    for (i=0;i<(gamsize/2);i++){
		if (actualfamsizes[i] > MAXFAMSIZE-1) actualfamsizes[i] = MAXFAMSIZE-1;
		famsizedist[actualfamsizes[i]]++;
	}
#endif
/* now add mutations */
/* first check to see how much room there is
if there is not much, run checkfix() to clear up space*/
if (mutavail <= mutavail_getting_low) checkfix(gpheap);
outofroom = _false;
for (i=0;i<= Wgamsize-1; i++){
  smuts = pickpoisson(total_u);
/* for debugging */

#ifdef MYDEBUG
   meansmuts +=smuts;
   meansmutcalls++;
   mutscalled += smuts;
#endif

   if ( smuts > mutavail ) {
      smuts = mutavail;
	  outofroom = _true;
   }
  for (j=0 ;j< smuts ;j++ ){
     mutspot = mutarray[mutavail-1];
	 /* for debugging 
	 meanmutspot += mutspot;
	 sqrmeanmutspot += (double) mutspot*mutspot;
	 nummutspot++; */
     mutavail--;
     add_inseg(gpheap+Wzpheap[nowpop]+i,mutspot);
     add_scoeff(mutspot); 
	 totalmuts++;
	 mutsincelastcheckfix++;
#ifdef MYDEBUG
     mutsadded++;
#endif
  }
} 

if (mutavail < minmutavail) minmutavail = mutavail;
if (outofroom){
	countoutofroom++;
	FP"\n ** Insufficient genome space for mutations in generation  %d **\n",nowgen);
	printf("\n ** Insufficient genome space for mutations in generation  %d **\n",nowgen);
	if (countoutofroom >= MAXOUTOFROOM) {
		nowgen = numgen+1;
		mutroomok = _false;
	}
}
} /* produce gametes */

void migrate(void){
int i,j,m,from, to[MAXMIGRANTS];
    j=0;
	if (numpops >1) for (nowpop=0;nowpop<numpops;nowpop++){
		m = pickpoisson(mrate);
		for (i=0;i<m && j<MAXMIGRANTS; i++,j++){
			from = randint(Wgamsize);
			do{	to[j] = randint(gamsize);
			}while ((to[j]>=Wgamsize*nowpop)  && (to[j]<Wgamsize*(nowpop+1)));
			memcpy(&migrants[j],gpheap+Wzpheap[nowpop]+from,genomesize);
		}
	}
	for (i=0;i<j; i++){
		memcpy(gpheap+to[i],&migrants[i],genomesize);
	}
} /* migrate */

void generation(void){
/*
have just 2 heaps, that get alternated by swapping pointers gpheap and zpheap
For each population
 fitgroups()
	assign fitness to the pairs of gametes, pointed to by zpheap
	construct array zptr_list of groups of pointers, each group is a fitness class.
	group the array of pointers in zptr_list accordingly
 get_famsizes);
	figure out how many offspring each zygote has
 produce_gametes();
	drawing on zptr_list pick gametes
	do recombination and add new gametes to new heap, pointed to by gpheap
	add mutations within gametes pointed to by gpheap

For all populations
  migrate()  - copy some gametes from one population and send them to another.

For each population
    pairs() - scramble the order of gametes in each population

swap zpheap and gpheap
zpheap points to heap of new gametes.
  then (zpheap+0) and (zpheap+1) are the two gametes that make the first zygote, etc


*/
	for (nowpop=0;nowpop < numpops; nowpop++){
		fitgroups();
		get_famsizes();
		produce_gametes();
		}
	migrate(); 
    for (nowpop=0;nowpop < numpops; nowpop++){
            pairs();
        }
	tempheapptr = gpheap;
	gpheap = zpheap;
	zpheap = tempheapptr;

} /* generation */


/********  MEASUREMENT *************** */



void getmeanfit(void){
/* if pop < 0 return mean fit of the entire population */
int i,hom,het;
float temp;
    popmeanfit = 0;
	for ( i=0;i <= gamsize-1 ; i+= 2){
		temp = assignfit(zpheap+i, &hom, &het);
		popmeanfit +=temp;
		}
	popmeanfit /= zsize;
} /* getmeanfit */


void measureW(void){
static int p,pp,k,i,j,ii,jj, hom, het;
static struct genome polyW[MAXPOPS][POLYNEUT+1],fixedW[MAXPOPS][POLYNEUT+1];
static float temphw[POLYNEUT+1],temphb[POLYNEUT+1], pophw[MAXPOPS][POLYNEUT+1],popfit[MAXPOPS];
static float tempsum1, tempsum2, varfit, varhet[POLYNEUT+1];
static tpophw[MAXPOPS][POLYNEUT+1] = {{0}};

for ( p=0;p<numpops;p++ )popfit[p] =0.0;
for (polys = POLYALL; polys <= POLYNEUT; polys++){
    temphw[polys]=0.0;
    temphb[polys]=0.0;
	varhet[polys] = 0.0;
    for ( p=0;p<numpops;p++ ){
		fixedW[p][polys] = fullgenome;
		polyW[p][polys] = emptygenome;
		pophw[p][polys] = 0.0;
		
    }
}
for ( i=0,p=0,ii=0;i<gamsize;i++, ii++ ){
    if ( i>=Wzpheap[p+1] ){
		p++;
		ii = 0;
	}
    for ( k=0;k<= allchromesegments-1 ;k++ ){
		for (polys = POLYALL; polys <= POLYNEUT; polys++){
			fixedW[p][polys].b[k] &= ((zpheap + Wzpheap[p] + ii)->b[k] & polygenome[polys].b[k]);
			polyW[p][polys].b[k]  |= ((zpheap + Wzpheap[p] + ii)->b[k] & polygenome[polys].b[k]);
		}
    }
}
for ( p=0;p<numpops-1;p++ ){
  for ( k=0;k<= allchromesegments-1 ;k++ )
    for (polys = POLYALL; polys <= POLYNEUT; polys++) polyW[p][polys].b[k] &= ~fixedW[p][polys].b[k];
}

for ( p=0;p<numpops-1;p++ ){
    for ( pp=p+1;pp<numpops;pp++ ){
		for ( k=0;k<= allchromesegments-1 ;k++ ){
			for (polys = POLYALL; polys <= POLYNEUT; polys++){
				countstypes[SSHARED][polys] += cardinality(polyW[p][polys].b[k] & polyW[pp][polys].b[k]);
				countstypes[SEXCLUSIVE][polys] += cardinality(polyW[p][polys].b[k] ^ polyW[pp][polys].b[k]);
				countstypes[SFIXED][polys] += cardinality(fixedW[p][polys].b[k] ^ fixedW[pp][polys].b[k]);
			}
		}
	}
}

for ( i=0,p=0;i<gamsize-1;i++)
    for ( j=i+1,pp=0;j<gamsize;j++){
		if ( i>=Wzpheap[p+1] ){
			p++; 
		}
		while ( j>=Wzpheap[pp+1] ){
			pp++;
		}
		myassert(pp >= p);
		myassert(j >= i);
		for (polys = POLYALL; polys <= POLYNEUT; polys++){
			k = countdiffs(zpheap+i,zpheap+j,polys);
			if ( p==pp ) {
				temphw[polys] += k;
				pophw[p][polys] += k;
			}
			else temphb[polys] += k;
		}
	}
for ( i=0,p=0, ii=0;i<gamsize;i += 2, ii+=2 ){
	if ( i>=Wzpheap[p+1] ){
		p++;
		ii = 0;
	}
    popfit[p] +=  assignfit(zpheap+Wzpheap[p]+ii,&hom,&het);
}
for (tempsum1=0.0, tempsum2=0.0, p=0;p<numpops;p++ ) {
	tempsum1 += popfit[p]*popfit[p];
	tempsum2 += popfit[p];
}
covarfit += tempsum2/numpops;
varfit = (tempsum1 -  tempsum2*tempsum2/numpops)/numpops;
meanvarfit += varfit;

for (polys = POLYALL; polys <= POLYNEUT; polys++){
	for (tempsum1=0.0, tempsum2=0.0, p=0;p<numpops;p++ ) {
		pophw[p][polys] /= ((float) (Wgamsize *(Wgamsize-1))/ (float) 2);
		tpophw[p][polys] += pophw[p][polys];
		tempsum1 += pophw[p][polys]*pophw[p][polys];
		tempsum2 += pophw[p][polys];
	}
	varhet[polys] = (tempsum1 -  tempsum2*tempsum2/numpops)/numpops;
	covarhet[polys] += tempsum2/numpops;
	meanvarhet[polys] += varhet[polys];
}

/* for Wgamsize items in each of numpop populations
there are a total of 
numpops * (Wgamsize *(Wgamsize-1)/2) comparisons within populations
    and 
Wgamsize^2 * numpops*(numpops-1)/2 comparisons between 
    */
for (polys = POLYALL; polys <= POLYNEUT; polys++){
    temphw[polys] /= (numpops * (float) (Wgamsize *(Wgamsize-1))/ 2.0);
    temphb[polys] /= (Wgamsize*Wgamsize* (float) (numpops*(numpops-1))/ 2.0);
	Hb[polys] += temphb[polys];
	Hw[polys] += temphw[polys];
	numHb[polys]++;
	numHw[polys]++;
    if (temphb[polys] > 0.0){
		numFst[polys]++;
		Fst[polys] += 1.0 - temphw[polys]/temphb[polys];
		if (temphb[polys] > temphw[polys]){
			numNm[polys]++;
			Nm[polys] += temphw[polys]/(4.0 * (temphb[polys]-temphw[polys]));
		}
	}
}
} /* measureW */


void measure(void){
static struct binfithet tabsbins[MAXBINS];
SEGMENT_TYPE tempseg;
int i,k,ii, inc,hom,het;
unsigned long  c,j;
static float temphet[POLYNEUT+1], tempmeanfit;
static float temp_s[POLYNEUT+1];
int tajd_s, tajd_pi;
struct genome tajd_bins_s, tajd_bins_f;
int oneperbin[MAXBINS];
int pseudolist[MAXBINS*MAXPOPS], tempinlist;
boolean notinlist;
unsigned long  ilong;
int usebins;

	inc = zsize / bins;
	tempmeanfit = 0.0;
	for (polys = POLYALL; polys <= POLYNEUT; polys++){
		temphet[polys] = 0.0;
	    temp_s[polys] = 0.0;
	    
		
	    for ( k=0;k<=allchromesegments-1 ;k++ ){
		      temp_s[polys] += cardinality(polygenome[polys].b[k]);
		}
	}
	for ( i=0,j=0;i <gamsize-1 ; i++){
	    if (!(i % 2)){
		fh[j].fit = assignfit(zpheap+i, &hom, &het);
		tempmeanfit += fh[j].fit;
		fh[j].onum = i;
	        fh[j].zptr1 = (zpheap+i);
	        fh[j].zptr2 = (zpheap+i+1);
		j++;
		}
	     for (polys = POLYALL; polys <= POLYNEUT; polys++){
	       for (ii=i+1;ii<gamsize;ii++){
	            het = countdiffs(zpheap+i,zpheap+ii,polys);
		    if (ii == i+1) fh[j].het[polys] = het;
		    temphet[polys] += het;
	        }
	     }
	}
	tempmeanfit /= zsize;
	overallmeanfit += tempmeanfit;
	shellfithet(&fh[0],zsize);
	for (i=0,k=zsize-1;i<= bins-1;i++){
		oneperbin[i] = k;
		tabsbins[i].fit=0.0;
		for (j=0;j< inc;j++,k--) tabsbins[i].fit += fh[k].fit;
		tabsbins[i].fit /= inc;
		absbins[i].fit += tabsbins[i].fit;
	}
    myassert((k+1) * 2 == lowbinbottom);

	for (polys = POLYALL; polys <= POLYNEUT; polys++) {
		tajd_s = 0;
	    tajd_pi = 0.0;
		tajd_bins_s = emptygenome;
		tajd_bins_f = fullgenome;
		for(i=0;i< bins;i++){
			for ( k=0;k<=allchromesegments-1 ;k++ ){
				tajd_bins_s.b[k] |= fh[oneperbin[i]].zptr1->b[k]  & polygenome[polys].b[k];
				tajd_bins_f.b[k] &= fh[oneperbin[i]].zptr1->b[k]  & polygenome[polys].b[k];
			}
			for (ii =i+1; ii< bins;ii++){
					tajd_pi += countdiffs(fh[oneperbin[i]].zptr1,fh[oneperbin[ii]].zptr1,polys);
				}
		}
		for ( k=0;k<=allchromesegments-1 ;k++ ){
		    tajd_s += cardinality(tajd_bins_s.b[k] & ~tajd_bins_f.b[k]);
			}
		tajd_sum[polys] += dotajd(bins,tajd_s,tajd_pi);
		tajd_n[polys]++;
	}
	
	for (polys = POLYALL; polys <= POLYNEUT; polys++) {
		temphet[polys] /= (gamsize * (gamsize-1))/2.0;
		meanhet[polys] += temphet[polys];
		mean_s[polys] += temp_s[polys];
		for (i=0,k=zsize-1;i<= bins-1;i++){
			tabsbins[i].het[polys]=0.0;
			for (j=0;j< inc;j++,k--) {
				tabsbins[i].het[polys] += fh[k].het[polys];
			}
			tabsbins[i].het[polys] /= inc;
			absbins[i].het[polys] += tabsbins[i].het[polys];
		}
	}
	if ((countints==num_measurements) && (print_pseudo_data)){
		for (i=0;i< pcharlines;i++) pchars[i] = calloc(num_chromosomes * chromesegments*SEGMENT_BITS, sizeof(char));
		for (i=0;i< pcharlines;i++) {
			for (ilong = 0;ilong< num_chromosomes * chromesegments*SEGMENT_BITS;ilong++) {
				pchars[i][ilong] = 'A';
			}
		}
		if (numpops > 1) usebins = Wbins;
		else usebins = bins;
		for (i=0, k=0;i< numpops;i++)
			for (j=0;j< usebins;j++,k++) {
				do {
					notinlist = _true;
					tempinlist = randint(Wgamsize) + Wzpheap[i];
					for (ii=0;ii< k;ii++) notinlist = notinlist && (tempinlist != pseudolist[ii]);
				} while (!notinlist);
				pseudolist[k] = tempinlist;
			}
     	myassert(k==pcharlines);
		for (i=0;i< pcharlines;i++){
			for ( k=0;k<=allchromesegments-1 ;k++ ){
			   for (polys = POLYALL; polys <= POLYNEUT; polys++) {
				   forall(j,(zpheap + pseudolist[i])->b[k] & polygenome[polys].b[k]){
					   switch (polys){
							case POLYALL : pchars[i][k*SEGMENT_BITS + j] = 'T'; break;
							case POLYPOS : pchars[i][k*SEGMENT_BITS + j] = 'G'; break;
							case POLYNEG : pchars[i][k*SEGMENT_BITS + j] = 'C'; break;
							case POLYNEUT: pchars[i][k*SEGMENT_BITS + j] = 'T'; break;
						}
				   }
			   }
			}
		}
	}
} /* measure */

float dotajd(int n,int ss, float pi){
static float khat, sf, e1, e2, a1, c1,b2, b1, c2, a2, D, nf;
static int  i, s;
  
  nf = (float) n;
  for (i = 1,a1 = 0;i<= n-1;i++) a1 += 1/(float) i;
  a2 = 1.0;
  for (i = 2; i <= n - 1; i++) a2 += 1.0/(float) (i*i);
  b1 = (nf +1.0)/(3.0*(nf -1.0));
  b2 = 2.0*(nf*nf+nf+3.0)/(9.0*nf*(nf-1.0));
  c1 = b1- 1/a1;
  c2 = b2 - (nf+2.0)/(a1*nf)+a2/(a1*a1);
  e1 = c1/a1;
  e2 = c2/(a1*a1+a2);
  khat = pi/(nf*(nf-1.0)/2.0);
  s = ss;
  sf = (float) ss;
  if (sf > 0.0) D = (khat - (sf/a1))/sqrt(e1*sf+e2*sf*(sf-1));
         else D = 0.0;
  return(D);
  }/* dotajd */

void calcLD(unsigned int x[4], float *LD){
	float p1,q1,p2,q2;
	float xf[4];
	int  i;
	float n,D,denom,r,rsq, Dmax,Dprime;
/*enum {LDR,LDRS,LDD,LDABS,LDDP}; */
	for (i=0,n=0;i< 4;i++) n += x[i];
	if ((n==0)||(x[0] + x[1] == 0) ||(x[0] + x[2] == 0)||(x[2] + x[3] == 0)||(x[1] + x[3] == 0)) {
		LD[0] = Dfault;LD[1] = Dfault;LD[2] = Dfault;LD[3] = Dfault; LD[4]=Dfault;
	} else {
	for (i=0;i< 4;i++) xf[i] = (float) x[i] /(float) n;
	p1 = xf[0] + xf[1];
	p2 = xf[2] + xf[3];
	q1 = xf[0] + xf[2];
	q2 = xf[1] + xf[3];
	D = xf[0]*xf[3] - xf[1]*xf[2];
	if (D>0){
		if (p1*q2 < p2*q1) Dmax = p1*q2; else Dmax = p2*q1;
	} else{
		if (p1*q1 < p2*q2) Dmax = p1*q1; else Dmax = p2*q2;
	}
	Dprime = D/Dmax;
	LD[2] = D;
	LD[3] = fabs(D);
	LD[4] = Dprime;
    denom = p1*q1*p2*q2;
    if (n > 0.0 && denom > 0.0){
		denom = sqrt(denom);
		r = D/denom;
        rsq = r*r;
		LD[0] = r;
		LD[1] = rsq;
    }else {
		LD[0] = Dfault;
		LD[1] = Dfault;
	}
	}
} /* calcLD */


void measureLD(void){
int i, j,k,ii, jj, seg,inseg,ir,kr, ic;
unsigned long ilong;
static long numpoly[POLYNEUT+1];
static struct points{
	int seg;
	int inseg;
	int region;
	int chromosome;
} polypoints[MAXALLCHROMELENGTH][POLYNEUT+1];


static unsigned long LDregion_borders[LDREGIONS];
static int LDfitbin_borders[MAXBINS];
unsigned int x[4],sum1,sum2,fitbin,sum1f,sum2f,xfitbin[4];
float LD[LDDP+1];
boolean derived1,derived2,LDok[LDDP+1];
static int  finc, rinc;
int tajd_mark[POLYNEUT+1], tajd_s_sum[POLYNEUT+1], tajd_pi[POLYNEUT+1];
float tajd[POLYNEUT+1],tajd_pi_sum[POLYNEUT+1];
int tajd_by_fitness_s[MAXBINS][POLYNEUT+1],tajd_by_fitness_pi[POLYNEUT+1];
float tajd_by_fitness[MAXBINS][POLYNEUT+1],tajd_by_fitness_pisum[MAXBINS][POLYNEUT+1];

if (countints==1){
	finc = gamsize/bins;
	rinc = chromelength/regions;
	for (i=0;i< regions;i++) LDregion_borders[i] = (i+1)*rinc;
	for (i=0;i< bins;i++) LDfitbin_borders[i] = lowbinbottom + (i+1)*finc;
	for (polys = POLYALL; polys <= POLYNEUT; polys++) {
		for (i=0;i< bins;i++) tajd_by_fitness_s[i][polys] = 0;
		for (i=0;i< bins;i++) tajd_by_fitness_pisum[i][polys] = 0.0;
	}
}

for (polys = POLYALL; polys <= POLYNEUT; polys++) {
	numpoly[polys] = 0;
	tajd_by_fitness_pi[polys] = 0;
	for (ilong=1,ir = 0,ic=0;ilong<= allchromelength;ilong++){
		seg = segpos(ilong,&inseg);
		while (seg > chromosome_ends[ic]) {
			ic++;
			ir = 0;
		}
		if (element(inseg,polygenome[polys].b[seg])){
			while (((( signed long) ilong - (signed long) (ic*chromelength)) > (signed long) LDregion_borders[ir]) && (ir < regions)) ir++;
			polypoints[numpoly[polys]][polys].seg = seg;
			polypoints[numpoly[polys]][polys].inseg = inseg;
			polypoints[numpoly[polys]][polys].region = ir;
			polypoints[numpoly[polys]][polys].chromosome = ic;
			numpoly[polys]++;
		}
	}
	for (fitbin = 0,jj=lowbinbottom;jj <= gamsize-1 ; jj++){
	 while (jj >= LDfitbin_borders[fitbin]) fitbin++;
	 if (jj %2){
		 for (LDmeasure = LDR; LDmeasure <= LDDP; LDmeasure++) if (LDmeasurebool[LDmeasure]){
			LD_mean_fitnesses[fitbin].n[LDmeasure][polys]++;
			LD_mean_fitnesses[fitbin].LD[LDmeasure][polys] += fh[jj/2].fit;
		 } 
	 }	
	}
}
for (polys = POLYALL; polys <= POLYNEUT; polys++) {
	tajd_mark[polys] = 0;tajd_s_sum[polys] = 0;
	tajd_pi[polys] = 0; tajd_pi_sum[polys] = 0.0;
	for (i= 0;i<numpoly[polys]-1;i++){
	  tajd_mark[polys] = 1;
	  for (j=i+1;j< numpoly[polys];j++){
		  sum1=0;sum2=0;
		  x[0]=0;x[1]=0;x[2]=0;x[3]=0,
		  fitbin = 0;
		  sum1f=0;sum2f=0;
		  xfitbin[0]=0,xfitbin[1]=0,xfitbin[2]=0,xfitbin[3]=0; 
		  for (jj=0,k=1;jj <= gamsize-1 ; jj++,k++){
			if (jj % 2){
				derived1 = element(polypoints[i][polys].inseg,(fh[jj/2].zptr1->b[polypoints[i][polys].seg]));
				derived2 = element(polypoints[j][polys].inseg,(fh[jj/2].zptr1->b[polypoints[j][polys].seg]));
			} else {
				derived1 = element(polypoints[i][polys].inseg,(fh[jj/2].zptr2->b[polypoints[i][polys].seg]));
				derived2 = element(polypoints[j][polys].inseg,(fh[jj/2].zptr2->b[polypoints[j][polys].seg]));
			}
			sum1f += derived1;
			sum2f += derived2;
			if (derived1 && derived2)   xfitbin[0]++;
			if (derived1  && !derived2) xfitbin[1]++;
			if (!derived1 && derived2)  xfitbin[2]++;
			if (!derived1 && !derived2) xfitbin[3]++;
			if (tajd_mark[polys]) tajd_by_fitness_pi[polys] += derived1;
    		        if (k == finc){
				for (LDmeasure = LDR; LDmeasure <= LDDP; LDmeasure++) if (LDmeasurebool[LDmeasure]) LDok[LDmeasure]=_false;

   				if (sum1f && sum2f) {
					calcLD(xfitbin, &LD[0]);
					for (LDmeasure = LDR; LDmeasure <= LDDP; LDmeasure++) if (LDmeasurebool[LDmeasure]) LDok[LDmeasure] =  (LD[LDmeasure] < Dfault_check) ;
				} 
				for (LDmeasure = LDR; LDmeasure <= LDDP; LDmeasure++) if (LDmeasurebool[LDmeasure]) if (LDok[LDmeasure]){
					if (polypoints[i][polys].chromosome == polypoints[j][polys].chromosome){
						kr = abs(polypoints[i][polys].region-polypoints[j][polys].region);
						LD_by_regions_by_fitness[kr][fitbin].n[LDmeasure][polys]++;
						LD_by_regions_by_fitness[kr][fitbin].LD[LDmeasure][polys] +=LD[LDmeasure];
					} else {
						LD_by_regions_by_fitness[LDREGIONS][fitbin].n[LDmeasure][polys]++;
						LD_by_regions_by_fitness[LDREGIONS][fitbin].LD[LDmeasure][polys] += LD[LDmeasure];
					}
				}
				if (tajd_mark[polys]){
					if ((sum1f > 0)&&(sum1f < finc)){ 
						  tajd_by_fitness_s[fitbin][polys]++;
						  tajd_by_fitness_pisum[fitbin][polys] += tajd_by_fitness_pi[polys] * (finc - tajd_by_fitness_pi[polys]);
						  tajd_by_fitness_pi[polys] = 0;
					  } else {
						  tajd_by_fitness_pi[polys] = 0;
					  }
				}
				fitbin++;
				sum1f=0;sum2f=0;
				xfitbin[0]=0,xfitbin[1]=0,xfitbin[2]=0,xfitbin[3]=0; 
				k = 0;
				if (tajd_mark[polys]) tajd_pi[polys] += derived1;
				sum1 += derived1;
				sum2 += derived2;
				if (derived1 && derived2)   x[0]++;
				if (derived1  && !derived2) x[1]++;
				if (!derived1 && derived2) 	x[2]++;
				if (!derived1 && !derived2) x[3]++;
			}
		  }
		  for (LDmeasure = LDR; LDmeasure <= LDDP; LDmeasure++) if (LDmeasurebool[LDmeasure]) LDok[LDmeasure]=_false;
   		  if (sum1 && sum2) {
			  calcLD(x, &LD[0]);
			  for (LDmeasure = LDR; LDmeasure <= LDDP; LDmeasure++) if (LDmeasurebool[LDmeasure]) LDok[LDmeasure] =  (LD[LDmeasure] < Dfault_check) ;
		  }
		  for (LDmeasure = LDR; LDmeasure <= LDDP; LDmeasure++) if (LDmeasurebool[LDmeasure]) if (LDok[LDmeasure]){
				LDall.n[LDmeasure][polys]++;
				LDall.LD[LDmeasure][polys] += LD[LDmeasure];
				if (polypoints[i][polys].chromosome == polypoints[j][polys].chromosome){
					kr = abs(polypoints[i][polys].region-polypoints[j][polys].region);
					LD_by_regions[kr].n[LDmeasure][polys]++;
					LD_by_regions[kr].LD[LDmeasure][polys] +=LD[LDmeasure];
				} else {
					LD_by_regions[LDREGIONS].n[LDmeasure][polys]++;
					LD_by_regions[LDREGIONS].LD[LDmeasure][polys] += LD[LDmeasure];
				}
		  }
		  if (tajd_mark[polys]){
			  if ((sum1 > 0)&&(sum1 < bins)){ 
				  tajd_s_sum[polys]++;
				  tajd_pi_sum[polys] += tajd_pi[polys] * (bins - tajd_pi[polys]);
				  tajd_pi[polys] = 0;
				  tajd_mark[polys] = 0;
			  } else {
				  tajd_pi[polys] = 0.0;
				  tajd_mark[polys] = 0;
			  }
			  
		  }


	  }
	}
	for (fitbin = 0;fitbin < bins;fitbin++){
		if (tajd_by_fitness_s[fitbin][polys] > 0){
			tajd_by_fitness[fitbin][polys] = dotajd(finc,tajd_by_fitness_s[fitbin][polys],tajd_by_fitness_pisum[fitbin][polys]);
			tajd_by_fitness_sum[fitbin][polys] += tajd_by_fitness[fitbin][polys];
			tajd_by_fitness_n[fitbin][polys]++;
		}
	}
}

} /* measureLD */
/********  INPUT AND OUTPUT RELATED  FUNCTIONS *************** */


void process_options(void){
int i,j,k;
boolean b_gamsize=_false;
boolean b_chromesegments=_false;
boolean b_seed_for_ran1 = _false;
boolean b_pop_u_s = _false;
boolean b_popscoeff = _false;
boolean b_del_prop = _false; 
boolean b_pop_u_n = _false; 
boolean b_poprecrate = _false; 
boolean b_num_measurements = _false; 
boolean b_dominance = _false;
boolean b_bp = _false;
boolean b_cytotype = _false;
boolean b_pseudodata = _false;
boolean b_fitschema = _false;
boolean b_epistatic_coeff = _false;
boolean b_burn = _false;
boolean b_between = _false;
boolean b_polys = _false;
boolean b_numpops = _false;
boolean b_migraterate = _false;
boolean b_absolutefit = _false;
int pargc,pslength;
char *pargv[MAXOPTIONS], ch;

DoLD = _false;
if ( p[strlen(p)-1] == '\n' ) p[strlen(p)-1] = '\0';
pslength = strlen(p);
i=0;
pargc = 0;
while ( i <= pslength ){
   if ( p[i] == OPTION_CHAR ){
      pargc++;
      pargv[pargc] = &p[i];
   }
   i++;
}
for (i=1; i<=pargc;i++){
  j=1;
 if (*pargv[i] == OPTION_CHAR)
     switch(toupper(pargv[i][j])){
     case 'G' : b_gamsize = _true;
                gamsize = atoi(&pargv[i][j+1]);
                if (gamsize % 2) gamsize++; /* add 1 if gamsize is odd */
                break;
     case 'C' : b_chromesegments=_true;
                chromesegments = atof(&pargv[i][j+1]);
				if (chromesegments > MAXALLCHROMESEGMENTS){
					printf("ERROR :  chromesegments > MAXALLCHROMESEGMENTS \n\n\n\n");
					exit(0);
				}
                break;
     case 'V' : b_pop_u_s = _true;
                pop_u_s = atof(&pargv[i][j+1]);
				if (pop_u_s == 0.0) {
					b_popscoeff = _true;
					scoeff = 0.0;
					b_dominance = _true;
					dominance = 0.0;
				}
                break;
     case 'S' : b_popscoeff = _true;
                popscoeff = atof(&pargv[i][j+1]);
                scoeff = popscoeff / gamsize;
                break;
     case 'F' : b_del_prop = _true;
                del_prop = atof(&pargv[i][j+1]);
                break;
     case 'U' : b_pop_u_n = _true;
                pop_u_n = atof(&pargv[i][j+1]);
                break;
     case 'R' : b_poprecrate=_true;
                poprecrate = atof(&pargv[i][j+1]);
                recrate = poprecrate /gamsize;  
                break;
     case 'X' : b_cytotype=_true;
                num_chromosomes = atoi(&pargv[i][j+1]);
                break;
     case 'N' : b_num_measurements=_true;
                num_measurements = atoi(&pargv[i][j+1]);
                break;
     case 'A' : b_seed_for_ran1 = _true;
                seed_for_ran1 = atoi(&pargv[i][j+1]);
                if (!seed_for_ran1) seed_for_ran1 = 1;
                /*srand(seed_for_ran1); */
                break;
	 case 'H' : b_dominance = _true;
				dominance = atof(&pargv[i][j+1]);
				break;
	 case 'B' : b_bp = _true;
				bp = atof(&pargv[i][j+1]);
				break;
	 case 'T' : b_pseudodata = _true;
				print_pseudo_data = _true;
				break;
	 case 'L' : DoLD = _true;
				for (LDmeasure = LDR; LDmeasure <= LDDP; LDmeasure++) LDmeasurebool[LDmeasure] = _false;
				k=0;
				while(isalpha(pargv[i][k+2]) && k<5){
					LDruntime[k] = pargv[i][k+2];
					k++;
				}
				LDruntime[k] = '\0';
				k = 2;
				ch = pargv[i][k];
				while (ch != ' ' && ch !='\n'){
					switch (toupper(ch)){
					   case 'R':
							 LDmeasurebool[LDR] = _true;
							 break;
						case 'S':
							 LDmeasurebool[LDRS] = _true;
							 break;
						case 'P' :
							 LDmeasurebool[LDDP] = _true;
							 break;
						case 'D' :
							 LDmeasurebool[LDD] = _true;
							 break;
						case 'B' :
							 LDmeasurebool[LDABS] = _true;
							 break;
						 }
					k++;
					ch = pargv[i][k];
				}
				break;
     case 'O' : for (polys=POLYPOS;polys <= POLYNEUT;polys++) polysbool[polys] = _false;
		        polysbool[POLYALL] = _true;
				b_polys = _true;
				k=0;
				while(isalpha(pargv[i][k+2]) && k<5){
					polyruntime[k] = pargv[i][k+2];
					k++;
				}
				polyruntime[k] = '\0';
				k = 2;
				ch = pargv[i][k];
				while (ch != ' ' && ch !='\0'){
					switch (toupper(ch)){
						case 'B':  polysbool[POLYPOS] = _true; break;
						case 'H':  polysbool[POLYNEG] = _true; break;
						case 'N':  polysbool[POLYNEUT] = _true; break;
						 }
					k++;
					ch = pargv[i][k];
				}
				break;
	 case 'Y' : b_epistatic_coeff = _true;
				epi = atof(&pargv[i][j+1]);
				break;
	 case 'W' : b_fitschema = _true;
			    ch = pargv[i][j+1];
				switch (toupper(ch)){
					case 'A':  fitmod = SADD; break;
					case 'M':  fitmod = SMULTI; break;
					case 'E':  fitmod = SEPI; break;
				}
				break;
	case 'I'  : b_burn = _true;
		        startpopgens = atof(&pargv[i][j+1]);
		        timestartup = round(((float) gamsize) * startpopgens);
				if (timestartup < 1) timestartup = 1;
	
		        break;
	case 'J'  : b_between = _true;
		        measure_interval = atof(&pargv[i][j+1]);
		        timebetween = (int) (((float) gamsize)* measure_interval);
				if (timebetween < 1) timebetween = 1;
		        break;
    case 'M' :   b_migraterate = _true;
                 mrate = atof(&pargv[i][j+1]);
                 break;
    case 'K' :   b_numpops = _true;
    	         numpops = atoi(&pargv[i][j+1]);
    	         break;
	case 'Q' :   b_absolutefit = _true;
				 fitness_absolute = _true;
				 resetfitness_absolute = _false;
    	         break;
     default : printf("Bad Option : %c\n",pargv[i][j]);
     }
 }

if (!b_gamsize){
      printf("define P gamsize \n");
      exit(0);
} 
if ( !b_chromesegments ){
      printf("define C  chromesegments \n");
      exit(0);
}
if ( !b_pop_u_s ){
	if (!b_popscoeff || popscoeff ==0.0) {
		pop_u_s = 0.0;
		popscoeff = 0.0;
		b_popscoeff = _true;
	}
   else {
	   printf("define G pop_u_s \n");
	   exit(0);
   }
}
if ( !b_popscoeff ){
   if (!b_pop_u_s) pop_u_s = 0.0;
   else {
	   printf("define S  popscoeff    \n");
		exit(0);
   }
}
if ( !b_del_prop ){
	if (pop_u_s ==0.0 || !b_pop_u_s) del_prop = 0.0;
	else {
		printf("define F  del_prop   \n");
		exit(0);
	}
}
if ( !b_pop_u_n ){
   printf("define U  pop_u_n   \n");
   exit(0);
}
if ( !b_poprecrate ){
   printf("define R  poprecrate   \n");
   exit(0);
}
if (!b_num_measurements  ){
   printf("define N num_measurements  \n");
   exit(0);
}
if (!b_dominance){
	if (!b_popscoeff || popscoeff ==0.0) dominance = 0.0;
	else {
		printf("define dominance \n");
		exit(0);
	}
}
if (!b_migraterate){
  if (b_numpops && numpops==1) mrate = 0.0;
  if (b_numpops && numpops>1){
	printf("multiple subpopulations but no migration rate defined\n");
	exit(0);
  } 
  if (!b_numpops) mrate = 0.0;
}
if (!b_numpops)
  if (b_migraterate){
	printf("migration rate defined, but not number of subpopulations\n");
	exit(0);
  } else {
	  numpops = 1;
	  mrate = 0.0;
  }
if (b_numpops && b_migraterate){
	while (mrate > (gamsize/numpops)/4) {
		printf("Migration rate (%f) is too high given population sizes (%d)  \n",mrate,gamsize/numpops);
		printf("Enter migration rate per population per generation  : ");
        scanf("%f",&mrate);
	}
}
if (!b_epistatic_coeff){
	if (fitmod == SEPI){
		printf("define coefficient of epistasis \n");
		exit(0);
	}
}
if (!b_fitschema) fitmod = SMULTI;
if (!b_burn){
	startpopgens = 4.0;
    timestartup = round(((float) gamsize) * startpopgens);
	if (timestartup < 1) timestartup = 1;
}
if (!b_between){
	measure_interval = 1.0;
	timebetween = (int) (((float) gamsize)* measure_interval);
	if (timebetween < 1) timebetween = 1;
}
if (!b_polys) {
	for (polys=POLYALL;polys <= POLYNEUT;polys++) polysbool[polys] = _false;
	polysbool[POLYALL] = _true;
}
if (!b_cytotype)  num_chromosomes = 1;
if (!b_bp) bp= 1.0;	
if (!b_pseudodata) print_pseudo_data = _false;
if (!b_absolutefit){
	fitness_absolute = _false;
	resetfitness_absolute = _false;
}
if (!b_seed_for_ran1){
   seed_for_ran1 = time(NULL);
   /*srand(seed_for_ran1);  effectively random seed */
   printf("time variable selected to seed random number generator : %u \n",seed_for_ran1);
}


} /* process options */

void process_filename(int argc, char *argv[]){
/* gets input and output file names */
int i;
assert(argc ==3);
for (i=1; i<argc;i++){
     switch(toupper(argv[i][1])){
     case 'D' : strcpy(dname,&argv[i][2]);
				if ((dfile = fopen(dname,"r")) == NULL){
					printf("Error opening text file for reading\n"); exit(0);}
                break;

     case 'R' : strcpy(rname,&argv[i][2]);
				if (strlen(rname) > (_FNSIZE-4)) strdelete(rname, _FNSIZE-3, _FNSIZE);
				strcat(rname, ".fpg");
				if ((rfile = fopen(rname,"w")) == NULL){
						printf("Error opening text file for writing\n"); exit(1);}
                break;
    }
}
}

void setup(int argc, char *argv[]){
establish_signal_handler();
printf("**** FPG ****  FORWARD POPULATION GENETIC SIMULATION \n");
if (argc==3) process_filename(argc, argv);
if (!dfile) {
   printf("NAME OF INPUT DATA FILE : \n");
   scanf("%s",dname);
   if ((dfile = fopen(dname,"r")) == NULL){
	printf("Error opening text file for reading\n"); exit(0);}

}

if (!rfile) {
   printf("NAME OF RESULTS FILE : \n");
   scanf("%s",rname);
   if (strlen(rname) > (_FNSIZE-4)) strdelete(rname, _FNSIZE-3, _FNSIZE);
   strcat(rname, ".fpg");
   if ((rfile = fopen(rname,"w")) == NULL){
	printf("Error opening text file for writing\n"); exit(1);}
}


printf("  input file: %s    this output file: %s \n",dname,rname);
FP"**** FPG ****  FORWARD POPULATION GENETIC SIMULATION \n");
FP"  input file: %s    this output file: %s \n",dname,rname);
#ifdef NOREC
FP"   NO RECOMBINATION PERMITTED \n");
poprecrate = 0.0;
#endif

FPLINE;
FP"  Basic Input Parameters\n");
FP"  Label  InputFlag  Meaning \n");
FP"  -----  ---------  --------\n");
FP"  popsize      -G   the # of haploid genome copies (gametes) in the population \n");
FP"  polytype     -O   types of polymorphisms to analysis \n");
FP"  chromeseqs   -C   the # of segments per chromosome, %d bits per segment \n",isegment_bits);
FP"  chromenum    -X   number of chromosomes \n");
FP"  pop_u_s      -V   population selected mutation rate \n");
FP"  popscoeff    -S   population selection coefficient (gamsize * selection coefficient) \n");
FP"  del_prop     -F   proportion of selected mutations that are deleterious\n");
FP"  dominance    -H   dominance of selected mutations \n");
FP"  fit_schema   -W   fitness model: A - additive, M - multiplicatative, E - epistatic \n");
FP"  epistasis    -Y   degree of epistatis among selected mutations \n");
FP"  pop_u_n      -U   population neutral mutation rate \n");
FP"  poprecrate   -R   population  recombination rate (gamsize* rec.rate per chromosome)\n");
FP"  mrate        -M   migration rate out of each population per generation \n");
FP"  numpops      -K   number of subpopulations \n");
FP"  measurements -N   the number of measurements to make \n");
FP"  bp           -B   length of each chromosome in base pairs \n");
FP"  rand         -A   random number generator seed \n");
FP"  t_start      -I   length of burn-in period in units of population size generations\n");
FP"  t_between    -J   time between measurements after burnin, in units of population size generations\n");
FP"  LD_analysis  -L   invokes linkage disequilibrium analysis \n");
FP"  pseudo_data  -T   causes a pseudo data set, in SITES Program format to be included in output \n");
FP"  absolutefit  -Q   causes fitness measurements to include all those mutations that fixed \n");
FP"\n Result Values     Meaning\n");
FP" ------ ------     -------\n");
FP"  startmutavail    the number of positions for mutations =  chromesegs*chromenum*i_segment_bits\n");
FP"  meanmutavail     average number of available mutation spots found at check time  \n");
FP"                      includes those not used, those used and lost, and those used and fixed \n");
FP"  minmutavail      minimum number of available mutation spots \n");
FP"\n");
fclose(rfile);
p = malloc(200);
}

void write_results(char *p){
int i, j, inc;
unsigned long  ilong;
float mper,temp,temp1,temp2,temp3;
char LDstr[7];
FP"\n   ----------------------------------------------------------------------------\n");
if (fitmod==SADD) FP"\n   FITNESS IS ADDITIVE, VALUES LESS THAN 0 ARE TRUNCATED AT ZERO \n");
if (fitmod==SEPI) FP"\n   FITNESS IS EPISTATIC \n");
if (fitmod==SMULTI) FP"\n   FITNESS IS MULTIPLICATIVE \n");
if (fitmod==SEPI) FP"    Epistatic Fitness model - epistatis coefficent : %6.3f \n",epi);
FP"\n");
if (resetfitness_absolute == _true) {
	FP" AT START: FITNESS VALUES INCLUDE ALL FIXED MUTATIONS \n");
	FP"**  MEAN FITNESS OVERAN FLOATING POINT - FIXED MUTATIONS IGNORED THEREAFTER **\n");
} else {
	if (fitness_absolute) FP" FITNESS VALUES INCLUDE ALL FIXED MUTATIONS \n");
    else FP" FITNESS VALUES *DO NOT* INCLUDE FIXED MUTATIONS \n");
}


FP"Runtime Parameters \n");
FP"------- ----------\n");
FP"\tpopsize     :\t %-4u \n", gamsize); 
FP"\tpolytype    :\t %-s  \n", polyruntime); 
FP"\tchromesegs  :\t %-3u \n", chromesegments); 
FP"\tchromenum   :\t %-2d \n", num_chromosomes);
FP"\tpop_u_s     :\t %-7.2f \n",pop_u_s); 
FP"\tpopscoeff   :\t %-7.2f \n",popscoeff); 
FP"\tdel_prop    :\t %-7.2f \n",del_prop); 
FP"\tdominance   :\t %-6.3f \n",dominance); 
FP"\tfit schema  :\t ");
switch (fitmod){
	case SADD   : FP"A\n");break;
	case SMULTI : FP"M\n");break;
	case SEPI   : FP"E\n");
		          FP"\tepistasis  :\t %-7.3f \n",epi);
		          break;
}

FP"\tpop_u_n     :\t %-7.2f \n",pop_u_n); 
FP"\tpoprecrate  :\t %-7.2f \n",poprecrate); 
FP"\tmrate       :\t %-7.3f \n", mrate);
FP"\tnumpops     :\t %-u    \n",numpops);
FP"\tmeasurements:\t %-4u \n", num_measurements); 
FP"\tbp          :\t %-11.0f \n",bp); 
FP"\trand        :\t %-9u \n", seed_for_ran1); 
FP"\tt_start     :\t %-4u \n", timestartup); 
FP"\tt_between   :\t %-4u \n", timebetween); 
FP"\tLD_analysis :\t %-s  \n", LDruntime);

if (mutroomok){
    FP"\n\n");
	FP"-------------------------------------------\n");
	FP"** BASIC MEASUREMENTS  - FOR DATA SET %2d **\n",datasetnum);
	FP"-------------------------------------------\n\n");
	FP"\tstartmutavail :\t %12u \n", (unsigned) chromesegments*num_chromosomes*isegment_bits);
	FP"\tmeanmutavail  :\t %12.1f \n",meanmutavail/(float)meanmutavail_counts);
	FP"\tminmutavail   :\t %12u \n", minmutavail); 
	overallmeanfit /= countints;
	FP"\n Total Number of Generations                  : %-7lu  in popsize units: %-8.2f\n",numgen,(float)numgen/(float)gamsize);
	FP  " Burn-In Period Ended in Generation           : %-7d  in popsize units: %-8.2f\n",burnend,(float)burnend/(float)gamsize);
	FP  " Number of Generations After Burn-In Period   : %-7lu  in popsize units: %-8.2f\n",fixcountgen,(float)fixcountgen/(float)gamsize);
	FP  " Overall Mean Population Fitness After Burn-In: %-12.5f \n\n",overallmeanfit);

	FP" Total Mutation Counts, # Fixed and # Lost, Mean Polymorphism Counts \n");
	FP"   # Fixed and Fixation Rate (per gen) are also shown ('') after Burn-In period \n");
	FP" ---------------------------------------------------------------------------------\n\n");
 
	FP"    Poly Class  # Mutations    #Fixed   #Lost    #Poly   Mean#Poly  '#Fixed'  'FixRate'  'FixRate' per PopSize Generations\n");
	FP"    ----------  -----------    ------  -------  -------  ---------  --------  ---------  ---------------------------------\n");
	for (polys = POLYALL; polys <= POLYNEUT; polys++) if (polysbool[polys]){
		FP"    %11s   %8lu   %8lu %8lu %8lu   %9.1f %8lu  %10.8f %10.5f\n",polystr[polys],mutcount[polys],fixcount[polys],lostcount[polys],polycount[polys],meanpolycount[polys]/numcheckfix,fixcountstart[polys],fixcountstart[polys]/(float)fixcountgen, fixcountstart[polys]*(float)gamsize/(float) fixcountgen);
	}

	inc = zsize / bins;
	FP"\n Mean Heterozgyosity (H), #PolySites, Wattersons's Theta, Tajima's D (n=%d)\n",bins);
	FP" ----------------------------------------------------------------------------\n\n");
	FP"    Poly Class  Mean Obsv.H  Mean #Poly  Mean ThetaW  Tajima's D\n");
	FP"    ----------  -----------  ----------  -----------  ----------\n");
	for (polys = POLYALL; polys <= POLYNEUT; polys++) if (polysbool[polys]){
		meanhet[polys] /= countints;
		mean_s[polys] /= countints;
		FP"   %11s     %9.2f    %7.1f    %10.5f   % -7.3f\n",polystr[polys],meanhet[polys],mean_s[polys],mean_s[polys]/a1,tajd_sum[polys]/(float) tajd_n[polys]);
	}
	if (numpops > 1){
		FP"\n\n");
		FP"--------------------------------------------------------------------\n");
		FP"** MEASUREMENTS REGARDING MULTIPLE POPULATIONS -  FOR DATA SET %2d **\n",datasetnum);
		FP"--------------------------------------------------------------------\n\n");
		FP"\n Mean Variance in Fitness Among Populations: %12.6f  \n",meanvarfit /= countints);
		FP"               Coefficient of Variation (%%):   %10.6f \n\n", 100.0*sqrt(meanvarfit /= countints)/(covarfit/countints) );

		FP" Mean Variance in Heterozgyosity Among Subpopulations, by Polymorphism Class \n");
		FP" ---------------------------------------------------------------------------\n");
		FP"  Poly Class  Mean Variance    Coefficient of Variation (%%)\n");
		FP"  ----------  -------------    ---------------------------- \n");
		for (polys = POLYALL; polys <= POLYNEUT; polys++) if (polysbool[polys]){
		     FP"   %11s     %12.5f    %12.7f \n",polystr[polys], meanvarhet[polys]/= countints,100.0*sqrt(meanvarhet[polys] /= countints)/(covarhet[polys]/countints) );
		}


		FP"\n Mean # Shared and Exclusive Polymorphisms, Fixed Differences, between populations\n");
		FP  " -----------------------------------------------------------------------------------------------\n");
		FP"  Poly Class  Mean # Shared    Mean # Exclusive    Mean # Fixed \n");
		FP"  ----------  -------------    ----------------    ------------  \n");
		for (polys = POLYALL; polys <= POLYNEUT; polys++) if (polysbool[polys]){
			temp1=countstypes[SSHARED][polys]/(countints*numpops*(numpops-1.0)/2.0);
			temp2=countstypes[SEXCLUSIVE][polys]/(countints*numpops*(numpops-1.0)/2.0);
			temp3=countstypes[SFIXED][polys]/(countints*numpops*(numpops-1.0)/2.0);
			FP"   %11s     %9.2f        %9.2f           %9.2f \n",polystr[polys],temp1,temp2,temp3);
		}

		FP"\n Mean Heterozgyosity (H) Within and Between, Mean Fst, and Mean Nm estimates \n");
		FP  " ---------------------------------------------------------------------------\n");
		FP"  Poly Class    H Within    H Between      Fst         Nm\n");
		FP"  ----------    --------    ---------   --------    ------- \n");
		for (polys = POLYALL; polys <= POLYNEUT; polys++) if (polysbool[polys]){
		     FP"   %11s   %9.2f   %9.2f   %7.4f    %8.4f \n",polystr[polys],Hw[polys]/numHw[polys],Hb[polys]/numHb[polys],Fst[polys]/numFst[polys],Nm[polys]/numNm[polys]);
	     }
	}
	FP"\n\n");
	FP"------------------------------------------------------\n");
	FP"** MEASUREMENTS BY FITNESS CLASS -  FOR DATA SET %2d **\n",datasetnum);
	FP"------------------------------------------------------\n\n");
	FP"  Fitness (F) and Heterozygosity (H) Distributions  Averaged over Measurements\n");
	FP"        relative fitness is a percentage of the overall population mean\n");
	FP"        relative heterozygosity is a percentage of the overall observed value \n");
	FP"        Tajima's D values are given for each fitness class (n=%d) \n",gamsize/bins);

	for (polys = POLYALL; polys <= POLYNEUT; polys++) if (polysbool[polys]){
		FP"\n Polymorphism Class : %11s \n",polystr[polys]);
		FP  " ------------ -----   -----------\n");
		FP"   Avg %% F Rank   Avg F    %% Relative F   Avg Obs H   %% Relative Obs H     ");
		if (DoLD) FP"Mean Tajima's D\n"); else FP"\n");
		FP"   ------------   -----     ------------   ---------   ----------------    ");
		if (DoLD) FP"---------------\n"); else FP"\n");
		for (i=0;i<=bins-1;i++){
			mper = 100.0*( (float) (i * inc) + (float) ((i+1) * inc))/(2.0 * zsize);
			FP"     %7.2f   %10.6f  %10.4f     %10.6f   %10.4f     ",
				mper, absbins[i].fit/countints,100*absbins[i].fit/countints/overallmeanfit,absbins[i].het[polys]/(bp*countints),100*absbins[i].het[polys]/(bp*countints)/meanhet[polys]);
			if (DoLD) FP"      % -7.3f\n",tajd_by_fitness_sum[i][polys]/tajd_by_fitness_n[i][polys]); else FP"\n");
		}
		FP"\n");
	}
	if (DoLD){
		FP"\n\n");
		FP"--------------------------------------------------------\n");
		FP"** LINKAGE DISEQUILIBRIUM ANALYSES  - FOR DATA SET %2d **\n",datasetnum);
		FP"--------------------------------------------------------\n\n");

		for (polys = POLYALL; polys <= POLYNEUT; polys++) if (polysbool[polys]){
			FP"\n Polymorphism Class : %11s \n",polystr[polys]);
			FP  " ------------ -----   -----------\n");
			
			for (LDmeasure = LDR; LDmeasure <= LDDP; LDmeasure++) if (LDmeasurebool[LDmeasure]) { 
				if (LDall.n[LDmeasure][polys]) LDall.LD[LDmeasure][polys] /= (float) LDall.n[LDmeasure][polys];
				for (i=0;i< regions+1;i++){
					if (LD_by_regions[i].n[LDmeasure][polys]) LD_by_regions[i].LD[LDmeasure][polys] /= LD_by_regions[i].n[LDmeasure][polys];
					for (j=0;j< bins;j++){
						if (LD_by_regions_by_fitness[i][j].n[LDmeasure][polys]) LD_by_regions_by_fitness[i][j].LD[LDmeasure][polys] /= LD_by_regions_by_fitness[i][j].n[LDmeasure][polys];
					}
				}
				for (j=0;j< bins;j++) if (LD_mean_fitnesses[j].n[LDmeasure][polys]) LD_mean_fitnesses[j].LD[LDmeasure][polys] /= LD_mean_fitnesses[j].n[LDmeasure][polys];
				if (LDmeasure==LDR)  strcpy(LDstr," (r) ");
				if (LDmeasure==LDRS) strcpy(LDstr,"(r^2)");
				if (LDmeasure==LDD)  strcpy(LDstr," (D) ");
				if (LDmeasure==LDABS) strcpy(LDstr,"(|D|)");
				if (LDmeasure==LDDP) strcpy(LDstr," (D')");
				FP"\n  Linkage Disequilibrium Measure: %s  \n",LDstr);
				FP"\n  Overall mean of LD (in samples of size %d): % 8.4f \n\n",bins,LDall.LD[LDmeasure][polys]);
				FP"  Mean LD between sites separated by a number of chromosome regions \n");
				FP"     Each chromosome was divided into %d regions, each of size %d bits \n",regions,chromelength/regions);
				for (i=0;i< regions;i++) FP" %8d ",i);
				if (num_chromosomes > 1) FP" Diff Chromosome \n");
				else FP"\n");
				for (i=0;i< regions;i++) FP"  % 7.4f ",LD_by_regions[i].LD[LDmeasure][polys]);
				if (num_chromosomes > 1) FP"  % 7.4f ",LD_by_regions[LDREGIONS].LD[LDmeasure][polys]);
				FP"\n\n");
				temp  = (float) countints*regions*regions;
				FP"    Mean LD by fitness group (# zygotes/group: %d), and by distance (in #'s of regions) between points \n",zsize/bins);
				FP" Fitness  Rank  Mean Fitness");
				for (i=0;i< regions;i++)  FP"   %2d     ",i);
				if (num_chromosomes > 1) FP" Diff Chromosome \n");
				else FP"\n");
				for (j = bins-1 ;j >= 0;j--){
					FP" %8d        %10.8f",j,LD_mean_fitnesses[j].LD[LDmeasure][polys]);
					for (i=0;i< regions;i++) FP"  % 7.4f ",LD_by_regions_by_fitness[i][j].LD[LDmeasure][polys]);
					if (num_chromosomes > 1) FP"  % 7.4f ",LD_by_regions_by_fitness[LDREGIONS][j].LD[LDmeasure][polys]);
					FP"\n");
				}
				FP"\n\n");
			} 
		}

	}

	if (print_pseudo_data){
		FPLINE;
		FP"Pseudo Data Set - suitable for analysis by SITES Program - created by FPG from last measurement \n");
		FP"# Data Set # %d  FPG input file : %s  FPG results file : %s  \n",datasetnum,dname,rname);
		FP"# FPG command line: %s \n",p);
		FP"# Individuals were ranked by fitness and placed in equal size bins.  \n");
		FP"# One sequence was drawn at random from each fitness bin \n");
		FP"# 'A' for ancestral base \n");
		if (numpops > 1) FP"# %d  populations \n",numpops);
		else FP"# 1 population \n");
		if (polysbool[POLYPOS]) FP"# 'G' for derived BENEFICIAL mutations\n");
		if (polysbool[POLYNEG]) FP"# 'C' for derived HARMFUL mutations\n");
		if (polysbool[POLYNEUT])FP"# 'T' for derived NEUTRAL mutations\n");
		FP"%d   %d \n",pcharlines,num_chromosomes * chromesegments*SEGMENT_BITS);
		FP"1   1 \n");
		FP"1   %d\n",num_chromosomes * chromesegments*SEGMENT_BITS);
		FP"%d\n",numpops);
		if (numpops > 1) for (i=1;i<= numpops;i++) FP"pop%d    %d\n",i,Wbins);
		for (i=0;i< pcharlines;i++){
			FP"%-10d",i);
			for (ilong = 0;ilong< num_chromosomes * chromesegments*SEGMENT_BITS;ilong++) FP"%c",pchars[i][ilong]);
			FP"\n");
		}
		FPLINE;
		for (i=0;i< pcharlines;i++) 
            free(pchars[i]);
	}
	} else FP"\n**Insufficient Genome Size, given the neutral mutation rate** \n");
}


void run(void){
boolean donestart;

donestart=_false;
while (nowgen < numgen){
	  nowgen++;
	  if (donestart) {
		  fixcountgen++;
	  }
      generation();
/*	  printf("%d\n",nowgen); */
	  if ((countints==0) && (nowgen % timestartup ==0)) {
		donestart = _true;
	  }
	  if (nowgen % (gamsize ) == 0){ 
		  getmeanfit();
		  printf("%2d  %7d %16.8f %5d\n",datasetnum,nowgen,popmeanfit,availablespace); 
		  FP" %6d %16.8f  ",nowgen,popmeanfit);
	  }
	  if (nowgen % (4*gamsize) == 0) FP"\n");
      if (donestart && (countints>0)){
         if (checkint==timebetween){
			 countints++;
			 checkfix(zpheap); 
			 measure();
			 if (numpops > 1) measureW();
			 if (DoLD) measureLD();
             checkint = 1;
         }
         else {
            checkint++;
         }
      } 
      if (donestart && (countints==0)) {
		 burnend = nowgen;
		 printf("*********************** Measurements Begin *******\n");
		 getmeanfit();
		 countints++;
         checkfix(zpheap); 
		 measure();
		 if (numpops > 1) measureW();
		 if (print_pseudo_data || DoLD) measureLD();
         donestart = _true;
         checkint=1;
		 fixcountgen = 0;
      }
   }
} /* run */

/********  MAIN FUNCTION *************** */
int main(int argc, char *argv[]) {
#ifdef  __MWERKS__
  argc = ccommand(&argv); 
  SIOUXSettings.autocloseonquit = _true;
  SIOUXSettings.asktosaveonclose = _false;
#endif
setup(argc, argv);
while ( fgets(p,200,dfile) != NULL ){
   if ( strlen(p) < 10 ) exit(0);
   datasetnum++;
   process_options();
   start();
   printf("\n\nParameter Set %d  \n",datasetnum);
   mutroomok = mutroom();
   if (!mutroomok) {
	   printf("\n**Insufficient Genome Size, given the neutral mutation rate** \n");

   }
   else {
	   run();
   }
   write_results(p);
   FP"\n\n  - FPG  program by Jody Hey  updated and compiled June 8, 2001 \n");
   //fclose(rfile);
   //free(gpheap);
   //free(zpheap);
}
//fclose(rfile);
//fclose(dfile);
} /* main */

/******** END *************** */

/*
Miscellaneous notes - roughly in sequence  


notes made 1/5/01
history of this program is a bit obscure
used, debugged and updated extensively in summer of 00 for the study of the effects of many weakly selected mutations
on polymorphism levels

*/


