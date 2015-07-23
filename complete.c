#include <stdio.h>
#include <stdlib.h>

#include "tools.h"
#include "seqtools.h"
#include "ldhat.h"

int sizeofpset=100;
long *idum;

void print_help(int argc, char *argv[]);

main(int argc, char *argv[]) {
	int i,n_pts=0,K,nd, p1, p2, n11, hap[4], ns=0;
	int *data, hapname[2][4], ct=0;
	double mu[2]={0.5,0.5}, theta[2]={0.0,0.0}, *log_lik, rho_max, **P;
	int n_threads=0, thread=-1;
	int start, stop;
	char outputfilename[MAXNAME];
	char buffer[20];
	char *in_str;
	FILE *ofp;

	print_help(argc, argv);

	strcpy(outputfilename, "");
	for(i = 0; i < argc; i++)
	{
		if(*argv[i] == '-')
		{
			in_str = argv[i];
			if(strcmp(in_str, "-n") == 0) ns = atoi(argv[i + 1]);				/* no. of sequences */
			/*if(strcmp(in_str, "-lowtheta") == 0) theta[0] = atof(argv[i + 1]);		/* Low Theta */
			/*if(strcmp(in_str, "-hightheta") == 0) theta[1] = atof(argv[i + 1]);		/* High Theta */
			if(strcmp(in_str, "-theta") == 0) { theta[0] = atof(argv[i + 1]);
												theta[1] = atof(argv[i + 1]);}		/* Theta */
			if(strcmp(in_str, "-rhomax") == 0) rho_max = atof(argv[i + 1]);			/* max 4Ner */
			if(strcmp(in_str, "-n_pts") == 0) n_pts = atoi(argv[i + 1]);			/* Number of points on grid */
			if(strcmp(in_str, "-split") == 0) n_threads = atoi(argv[i + 1]);		/* Split the calculation into n threads */
			if(strcmp(in_str, "-element") == 0) thread = atoi(argv[i + 1]);			/* Current thread */
			if(strcmp(in_str, "-prefix") == 0) strcpy(outputfilename, argv[i+1]);					/* output file prefix */
		}
	}

	printf("\n\nEstimating coalescent likelihoods for complete set of configurations using Fearnhead and Donnelly (2001) method\n\n");
	if (ns == 0)
	{
		printf("\n\nInput number of chromosomes: ");
		scanf("%i", &ns);
	}
	if (theta[0] == 0.0)
	{
		printf("\n\nInput low theta per site: ");
		scanf("%lf", &theta[0]); 
	}
	if (theta[1] == 0.0)
	{
		printf("\n\nInput high theta per site: ");
		scanf("%lf", &theta[1]);
	}
	if (rho_max == 0.0)
	{
		printf("\n\nInput max value of 4Ner to calculate for (<=100): ");
		scanf("%lf", &rho_max);
	}
	if (n_pts == 0)
	{
		printf("\n\nInput number of points for grid: ");
		scanf("%i", &n_pts);
	}

	p1 = (int) ns/2;
	n11 = 1+p1+p1*(p1-1)*(p1+4)/6+(p1-1)*(p1+2)/2;

	if ((n_threads == 0) || (thread == -1))
	{
		start = 1;
		/*stop = (int) ns/2; */
		stop = n11+1;
		strcat( outputfilename, "new_lk.txt" );
	}
	else
	{
		start = n11 * ((double) thread / (double)n_threads)+1;
		stop = n11 * ((double) (thread + 1) / (double)n_threads)+1;
		sprintf(buffer, "%i", thread);
		strcat(outputfilename, "new_lk" );
		strcat(outputfilename, buffer);
		strcat(outputfilename, ".txt");
	}


	if (thread < 1)
	{
		ofp = fopen(outputfilename, "w");
		fprintf(ofp,"%i %i\n1 %.4f\n%i %.4f\n\n", ns, n11, theta[0], n_pts, rho_max);
		fprintf(ofp,"\n");
		fclose(ofp);
	}
	else
	{
		ofp = fopen(outputfilename, "w");
		fclose(ofp);
	}

	printf("\n\nFor %i chromosomes: total of %i pair configurations\n\n", ns, n11);

	K=2;
	log_lik=(double *)calloc(n_pts,sizeof(double));
	P=(double **)calloc(K,sizeof(double *));
	for(i=0;i<K;i++)
		P[i]=(double *)calloc(K,sizeof(double));

	P[0][0]=0.0;
	P[0][1]=1.0;
	P[1][0]=1.0;
	P[1][1]=0.0;
/*  P[1][2]=0.0;
  P[1][3]=0.0;
  P[2][0]=0.0;
  P[2][1]=0.0;
  P[2][2]=0.5;
  P[2][3]=0.5;
  P[3][0]=0.0;
  P[3][1]=0.0;
  P[3][2]=0.5;
  P[3][3]=0.5;
  */

	hapname[0][0]=hapname[0][1]=hapname[1][0]=hapname[1][2]=1;
	hapname[0][2]=hapname[0][3]=hapname[1][1]=hapname[1][3]=2;

	for (p1=1;p1<=(int) ns/2; p1++)
		for (p2=1;p2<=p1;p2++)
			for (n11=p2;n11>=0;n11--)
			{
				ct++;
				if ((ct >= start) & (ct < stop))
				{
					for (i=0;i<n_pts;i++) log_lik[i]=0.0;
					hap[3]=n11;
					hap[2]=p1-n11;
					hap[1]=p2-n11;
					hap[0]=ns-p1-p2+n11;
					for (i=0,nd=0;i<4;i++) if (hap[i]) nd++;
					data = (int *) malloc((size_t) 3*nd*sizeof(int));
					for (i=0,nd=0;i<4;i++) if (hap[i])
					{
						data[3*nd]=hapname[0][i];
						data[3*nd+1]=hapname[1][i];
						data[3*nd+2]=hap[i];
						nd++;
					}
					ofp = fopen(outputfilename, "a");
					printf("\nhap %3i: 00: %3i  01: %3i  10: %3i  11: %3i\n",ct,hap[0],hap[1],hap[2],hap[3]); 	
					pairs(nd,data,K,P,mu,theta,rho_max,n_pts,500000,log_lik);
					fprintf(ofp,"%i # %3i %3i %3i %3i : ", ct, hap[0], hap[2], hap[1], hap[3]);
					for (i=0;i<n_pts;i++) fprintf(ofp, "%10.3f ", log_lik[i]);
					fprintf(ofp,"\n");
					fclose(ofp);
					free(data);
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

				printf("\ncomplete\n");
				printf("Generate Likelihood lookup files for use with LDhat programs.\n\n");
				
				printf("Required Options :\n");
				printf("-n <int>                     Number of sequences\n");
				printf("-rhomax <float>              Maximum 4Ner\n");
				printf("-n_pts <int>                 Number of points on grid\n");
				/*printf("-lowtheta <float>            Low value of theta\n");*/
				/*printf("-hightheta <float>           High value of theta\n");*/
				printf("-theta <float>               Population Mutation Rate per SNP\n");
				printf("\n\n");	
				printf("Advanced Options :\n");
				printf("-split <int>                 Split the calculation into elements\n");
				printf("-element <int>               Element number to calculate (first element is 0)\n");
				printf("-prefix <string>             Prefix of output files\n");
				
				printf("\n\n");	
				exit (0);
			}
		}
	}
}
