#include "ldhat.h"
#include "tools.h"
#include "seqtools.h"

int sizeofpset=100;
long *idum;

void print_help(int argc, char* argv[]);

main (int argc, char *argv[]) 
{
	int i, npt, tcat=1, verb=0, npmax, nseq, ns=0, rcat;
	int p1, p2, n11, ct;
	char fname[MAXNAME+1];
	double theta, rmax;
	double **lkmat, **lkmatn, **tmp;
	struct site_type **pset;
	struct data_sum data;
	
	FILE *ifp=NULL;
	char *in_str;

	/* See if user is asking for help */
	print_help(argc, argv);
	strcpy(data.prefix, "");
	for(i = 0; i < argc; i++)
	{
		if(*argv[i] == '-')
		{ 
			in_str = argv[i];
			if(strcmp(in_str, "-lk") == 0) ifp = fopen(argv[i+1], "r");	
			if(strcmp(in_str, "-nseq") == 0) ns = atoi(argv[i + 1]);
			if(strcmp(in_str, "-concise") == 0) verb = 0;
			if(strcmp(in_str, "-prefix") == 0) strcpy(data.prefix, argv[i+1]);
		}
	}
    
	if (ifp == NULL)
	{
		printf("\nCould not find likelihood file in command line.");
		printf("\n\nInput name of likelihood file: ");
		scanf("%s", fname);
		ifp = fopen(fname, "r");
	}
	
	if (ifp==NULL) nrerror("Cannot open likelihood file");
	
	fscanf(ifp,"%i", &nseq);
	
	rewind(ifp);
	printf("\n\nInput data for n=%i\n\n", nseq);
	
	read_pars(ifp, &tcat, &theta, &rcat, &rmax);
	
	rewind(ifp);
	i = (int) (nseq/2);
        
	npmax = (int) 1+i+i*(i-1)*(i+4)/6+(i-1)*(i+2)/2;
	
	printf("\n\nMax number of types for n=%i is %i\n\n",nseq,npmax);
	
	lkmat = dmatrix(1,npmax,1,rcat);
	lkmatn = dmatrix(1,npmax,1,rcat);
	
	data.nseq = nseq; 
	data.hd=1; 
	data.exact=1;
	
	pset = init_pset(pset,1,ifp,&npt,&data);
	read_lk(ifp, lkmat, npt, tcat, rcat);
	check_exhaustive(pset,npt,nseq);

	if (ns <= 0)
	{
		printf("\n\nInput number of sequences to output likelihood file for:");
		scanf("%i", &ns);
		if (ns>nseq || ns<2) nrerror("Impossible output number");
	}
	i = (int) (ns/2);
	npmax = (int) 1+i+i*(i-1)*(i+4)/6+(i-1)*(i+2)/2;
	printf("\n\nMax number of output pair types = %i\n\n", npmax);

	for (p1=1,ct=0;p1<=i;p1++) 
		for (p2=1;p2<=p1;p2++) 
			for (n11=p2;n11>=0;n11--) 
			{
				ct++;
				pset[ct]->pt[3]=n11;
				pset[ct]->pt[1]=p1-n11;
				pset[ct]->pt[2]=p2-n11;
				pset[ct]->pt[0]=ns-p1-p2+n11;
			}
				  
	for (i=1;i<=nseq-ns;i++)
	{
		calc_nless1(nseq,nseq-i,lkmat,lkmatn, rcat);
		tmp=lkmat;
		lkmat=lkmatn;
		lkmatn=tmp;
	}
	
	print_lks(pset, ns, npmax, lkmat, tcat, theta, rcat, rmax, data.prefix);
	free_dmatrix(lkmat,1,npmax,1,rcat);
	free_dmatrix(lkmatn,1,npmax,1,rcat);

	printf("\n");
/*	system("PAUSE");*/
}


/*Calculate likelihood file for n-1*/

void calc_nless1(nseq0, nlessi, lkmat, lkmatn, rcat)
int nseq0, nlessi, rcat;
double **lkmat, **lkmatn;
{
	int i, p1, p2, n11, ptt, ptn, *hap, *hap_aug, nprev;
	int p1n, p2n;

	hap = (int *) malloc((size_t) 10*sizeof(int));
	hap_aug = (int *) malloc((size_t) 10*sizeof(int));
	nprev=nlessi+1;
	printf("\nn = %3i -> %3i",nprev, nlessi);

	for (i=0;i<10;i++) hap[i]=hap_aug[i]=0;
	for (p1=1;p1<=(int) nlessi/2;p1++) {
		for (p2=1;p2<=p1;p2++) {
			for (n11=p2;n11>=0;n11--) {
				hap[3]=n11;
				hap[1]=p1-n11;
				hap[2]=p2-n11;
				hap[0]=nlessi-p1-p2+n11;
				ptt=(int) p1*(p1-1)*(p1+4)/6+(p2-1)*(p2+2)/2+hap[2]+1;

/*Add first augmentation: 11*/
				for (i=0;i<4;i++) hap_aug[i]=hap[i];
				hap_aug[0]++; order_pt_hap(hap_aug,nprev);
				p1n = hap_aug[1]+hap_aug[3];
				p2n = hap_aug[2]+hap_aug[3];
				ptn = (int) p1n*(p1n-1)*(p1n+4)/6+(p2n-1)*(p2n+2)/2+hap_aug[2]+1;
				for (i=1;i<=rcat;i++) lkmatn[ptt][i]=lkmat[ptn][i];

/*Add second augmentation: 10*/
				for (i=0;i<4;i++) hap_aug[i]=hap[i];
				hap_aug[1]++; order_pt_hap(hap_aug,nprev);
				p1n = hap_aug[1]+hap_aug[3];
				p2n = hap_aug[2]+hap_aug[3];
				ptn = (int) p1n*(p1n-1)*(p1n+4)/6+(p2n-1)*(p2n+2)/2+hap_aug[2]+1;
				for (i=1;i<=rcat;i++) lkmatn[ptt][i] += (double) log(1+exp(lkmat[ptn][i]-lkmatn[ptt][i]));

/*Add third augmentation: 01*/
				for (i=0;i<4;i++) hap_aug[i]=hap[i];
				hap_aug[2]++; order_pt_hap(hap_aug,nprev);
				p1n = hap_aug[1]+hap_aug[3];
				p2n = hap_aug[2]+hap_aug[3];
				ptn = (int) p1n*(p1n-1)*(p1n+4)/6+(p2n-1)*(p2n+2)/2+hap_aug[2]+1;
				for (i=1;i<=rcat;i++) lkmatn[ptt][i] += (double) log(1+exp(lkmat[ptn][i]-lkmatn[ptt][i]));

/*Add fourth augmentation: 00*/
				for (i=0;i<4;i++) hap_aug[i]=hap[i];
				hap_aug[3]++; order_pt_hap(hap_aug,nprev);
				p1n = hap_aug[1]+hap_aug[3];
				p2n = hap_aug[2]+hap_aug[3];
				ptn = (int) p1n*(p1n-1)*(p1n+4)/6+(p2n-1)*(p2n+2)/2+hap_aug[2]+1;
				for (i=1;i<=rcat;i++) lkmatn[ptt][i] += (double) log(1+exp(lkmat[ptn][i]-lkmatn[ptt][i]));
			}
		}
	}

	free(hap);
	free(hap_aug);
}


/*Check that likelihood file is exhaustive for n*/
void check_exhaustive(pset,npt,nsamp)
	struct site_type **pset;
	int npt, nsamp;
{
	int p1, p2, i, ei;

	printf("\n\nChecking likelihood file is exhaustive:...");
	p1 = (int) nsamp/2;
  
	if (npt != (int) 1+p1+p1*(p1-1)*(p1+4)/6+(p1-1)*(p1+2)/2)
	{
		printf("\n\n!!npt = %i: E[npt] = %i\n\n",npt, 1+p1+p1*(p1-1)*(p1+4)/6+(p1-1)*(p1+2)/2);
		nrerror("npt and total do not agree: not exhaustive");
	}

  
	for (i=1;i<=npt;i++)
	{
		p1 = pset[i]->pt[1]+pset[i]->pt[3];
		p2 = pset[i]->pt[2]+pset[i]->pt[3];
		ei = (int) p1*(p1-1)*(p1+4)/6+(p2-1)*(p2+2)/2+pset[i]->pt[2]+1;
    
		if (i!=ei)
		{
			printf("\nError in pair-types: p1(%i), p2(%i), n11(%i): i(%i) ei(%i)\n\n",p1,p2,pset[i]->pt[3],i,ei);
			exit(1);
		}

	}
  
	printf("OK\n\n");
}


void print_lks(pset,nseq,npt,lkmat,tcat,theta,rcat,rmax, prefix) 
int nseq,npt,tcat,rcat;
double theta, rmax;
double **lkmat;
struct site_type **pset;
char *prefix;
{
	char fname[MAXNAME+1];
	int p, c1, c2, c3;
	FILE *ofp;
	
	strcpy(fname, prefix);
	
	ofp = fopen(strcat(fname, "new_lk.txt"), "w");
	fprintf(ofp, "%i %i\n%i ",nseq,npt,tcat);
	
	for (c1=1; c1<=tcat; c1++)
		fprintf(ofp, " %f ", theta);
	
	fprintf(ofp,"\n%i %f\n",rcat, rmax);
	fprintf(ofp, "\n\nType    00  01  10  11 Rho");
	
	for (c3=0; c3<rcat; c3++)
		fprintf(ofp, "%7.1f ",(double) c3*rmax/(rcat-1));
	
	fprintf(ofp,"\n\n");
	for (p=1; p<=npt; p++)
	{
		fprintf(ofp,"%4i # ", p);
		for (c1=0; c1<4; c1++)
			fprintf(ofp,"%3i ", pset[p]->pt[c1]);
		
		fprintf(ofp," :  ");
		
		for (c1=0; c1<tcat; c1++)
			for (c2=0; c2<tcat; c2++) 
				for (c3=1; c3<=rcat; c3++) 
					fprintf(ofp,"%7.2f ", lkmat[p][rcat*tcat*c1+rcat*c2+c3]);
		
		fprintf(ofp,"\n");
	}
	
	fclose(ofp);
}


void print_help(int argc, char* argv[]) 
{
	int i;
	char *in_str;
	for(i = 0; i < argc; i++){
		if(*argv[i] == '-'){ 
			in_str = argv[i];
			if((strcmp(in_str, "-h") == 0) || (strcmp(in_str, "-?") == 0) || (strcmp(in_str, "-help") == 0)) {
				printf("\nlkgen\n");
				printf("Create a lookup table from an existing, larger lookup table.\n\n");
				
				printf("Required Options :\n");
				printf("-lk <file>             Input lookup table.\n");
				printf("-nseq <int>            Number of sequences for output file.\n");
				printf("\n\n");	
				printf("Additional Options :\n");
				printf("-concise               Concise Output.\n");
				printf("-prefix <string>       Prefix of output files\n");
				printf("\n\n");	
				exit(0);
			}
		}
	}
}
