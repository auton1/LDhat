/******************************************************************/
/* Stat.c - routine to summarise output of Interval MCMC program */
/******************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "tools.h"
#include "ldhat.h"

long *idum;
void print_help(int argc, char *argv[]);

main(int argc, char *argv[]) {

  int nline, ndat, i, j, mid, l95, u95, k, nmiss=0;
  long seed = -setseed();
  double **dat, av, *arr;
  double drate[4];
  double *locs;
  int ns;
  double tlseq;
  char lc;
  int loc_file = 0;
  char fname[MAXNAME+1];
  char prefix[MAXNAME+1];
  FILE *ifp, *ofp;
  FILE *loc=NULL;

  char *in_str;
  int ask_questions = 1;

  print_help(argc, argv);
  strcpy(prefix, "");
  for(i = 0; i < argc; i++)
	{
		if(*argv[i] == '-')
		{ 
			ask_questions = 0;
			in_str = argv[i];
			if(strcmp(in_str, "-input") == 0) ifp = fopen(argv[i+1], "r");
			if(strcmp(in_str, "-burn") == 0) nmiss = atoi(argv[i+1]);
			if(strcmp(in_str, "-loc") == 0) {loc = fopen(argv[i+1], "r"); loc_file = 1;}
			if(strcmp(in_str, "-prefix") == 0) strcpy(prefix, argv[i+1]);
		}
	}
 
  if((ifp == NULL) || (ask_questions == 1)) 
  {
	  printf("\nCould not find input file in command line\n");
	  printf("\nInput file name:");
	  scanf("%s", fname);
	  ifp = fopen(fname, "r");
  }
  if (!ifp) nrerror("Cannot open file");
  idum = &seed;
  
  	if (loc != NULL)
  	{
		fscanf(loc, "%i %lf %c", &ns, &tlseq, &lc);
	  	locs = dvector(0, ns);
	  	locs[0] = -1.0;

		for (i=1; i<=ns; i++) 
		{
			fscanf(loc, "%lf", &locs[i]); 
			if (i>1 && locs[i]<=locs[i-1]) nrerror("Error in locs file: SNPs must be monotonically increasing");
		}
		fclose(loc);
	}

  if (ask_questions) 
  {
	  printf("\nNumber of burn-in samples:");
	  scanf("%i", &nmiss);
  }

  if (nmiss<0) nmiss=0;

  fscanf(ifp,"%i %i", &nline, &ndat);
  if (nline < 2) nrerror("Too few points for analysis");
  dat = dmatrix(1,nline,1,ndat);
  if (nmiss>nline) nrerror("Miss more than have data points!!");

  printf("\n\nReading data: %i lines (miss first %i) of %i points......",nline,nmiss,ndat);

  for (i=1;i<=nline;i++) for (j=1;j<=ndat;j++) {
    fscanf(ifp,"%lf", &dat[i][j]);
    if (feof(ifp)) nrerror("Reached end of file unexpectedly");
  }

  printf("...Data read successfully\n\n");
  fclose(ifp);
  
  if ((loc_file != 0) && (ndat != ns))
  {
  	printf("Loc file does not match input file. Ignoring\n\n");
  	loc_file = 0;
  }

  arr = dvector(1,nline-nmiss);
  strcpy(fname, prefix);
  ofp = fopen(strcat(fname, "res.txt"), "w");
  mid = (int) ((double) (nline-nmiss)/2);
  l95 = (int) ((double) (nline-nmiss)*0.025+1.0);
  u95 = (int) ((double) (nline-nmiss)*0.975);

  printf("\n\nMid = %i, L95=%i, U95=%i\n\n", mid, l95, u95);


  fprintf(ofp,"Loci\tMean_rho\tMedian\tL95\tU95");
  for (i=1;i<=ndat;i++) {
    printf("."); if (!i%50) printf("\n");
    for (j=1+nmiss,av=0;j<=nline;j++) {
      av+=dat[j][i];
      arr[j-nmiss]=dat[j][i];
    }
    av/=(double) (nline-nmiss);
    /*    printf("Initial  array: "); for (k=1;k<=nline-nmiss;k++) printf(" %8.5f",arr[k-nmiss]);*/
    sort_farray(arr,nline-nmiss);
    /*    printf("\nSorted array: "); for (k=1;k<=nline-nmiss;k++) printf(" %8.5f",arr[k-nmiss]);
    printf("\n\n");
    */

	if (loc_file == 0)
	    fprintf(ofp,"\n%6i\t%10.5f\t%10.5f\t%10.5f\t%10.5f",i-1,av,arr[mid],arr[l95],arr[u95]);
	else
	    fprintf(ofp,"\n%10.3f\t%10.5f\t%10.5f\t%10.5f\t%10.5f",locs[i-1],av,arr[mid],arr[l95],arr[u95]);
  }

  /*Count total rate jumped over region*/
  for (i=1+nmiss, drate[0]=drate[2]=0.0;i<=nline;i++) {
    for (j=2,drate[1]=drate[3]=0.0;j<ndat-1;j++) {
      drate[1]+=(double) fabs(dat[i][j+1]-dat[i][j]); 
      if (dat[i][j+1]!=dat[i][j]) drate[3]++;
    }
    drate[0]+=drate[1];
    drate[2]+=drate[3];
  }
  drate[0]/=(double) nline-nmiss;
  drate[2]/=(double) nline-nmiss;
  printf("\n\nAverage total change in rate = %.3f\nAverage total # changes = %.3f\n\n", drate[0],drate[2]);
  
  if (loc_file != 0)
  	free_dvector(locs, 0, ns);

  fclose(ofp);
}

sort_farray(double *arr, int n) {

  int pass, i;
  double tmp;

  for (pass=1;pass<=n;pass++) for (i=1;i<n;i++) {
    if (arr[i+1]<arr[i]) {
      tmp=arr[i+1];
      arr[i+1]=arr[i];
      arr[i]=tmp;
    }
  }
}


void print_help(int argc, char* argv[]) 
{
	int i;
	char *in_str;
	for(i = 0; i < argc; i++){
		if(*argv[i] == '-'){ 
			in_str = argv[i];
			if((strcmp(in_str, "-h") == 0) || (strcmp(in_str, "-?") == 0) || (strcmp(in_str, "-help") == 0)) {

				printf("\nstat\n");
				printf("Summarise output from interval\n\n");	
				printf("Required Options :\n");
				printf("-input <file>       Input file (rates.txt or bounds.txt)\n");
				printf("\n\n");	
				printf("Additional Options :\n");
				printf("-burn <int>         Number of burn-in samples (default 0)\n");
				printf("-loc <file>         SNP loci file\n");
				printf("-prefix <string>    Prefix of output files\n");
			
				printf("\n\n");	
				exit (0);
			}
		}
	}
}

