#include "ldhat.h"
#include "tools.h"
#include "seqtools.h"

void print_help(int argc, char* argv[]);


long *idum;
int sizeofpset;

main (int argc, char *argv[]) {
	int i, j, **seqs, **nall, l, u, fall=1, *fsnp, site, hd, nth, jmin, nout=0, *seq_out_list, out_ct;
	int nseq, nmin, lseq, fl=1, na, psite;

	long seed = -setseed();
	double fcut=0.0, pwd, sn, ns, ctaj[12], n, d, spwd[2], tlseq, *sloc, fgap=1.0;
	char fname[MAXNAME+1], bases[6]="0?TCAG", sc[6]="0?0120", c, **seqnames;
	char prefix[MAXNAME+1];

	int loc_file = 0;
	FILE *ifp=NULL, *ofp, *loc, *inloc=NULL;

	int ask_questions = 1;

	char *in_str;

	idum = &seed;
	print_help(argc, argv);
	strcpy(prefix, "");
	
	for(i = 0; i < argc; i++)
	{
		if(*argv[i] == '-')
		{ 
			in_str = argv[i];
			ask_questions = 0;
			if(strcmp(in_str, "-seq") == 0) 
			{
				ifp = fopen(argv[i+1], "r");
			}
			if(strcmp(in_str, "-loc") == 0) 
			{
				printf("\nUsing input location file: %s\n\n", argv[2]);
				inloc = fopen(argv[i+1], "r");
				if (!inloc) {printf("\n\nCannot open input locs\n\n"); exit(1);}
			}
			if(strcmp(in_str, "-freqcut") == 0) {fl = 2; fcut = (double) atof(argv[i+1]);}
			if(strcmp(in_str, "-missfreqcut") == 0) fgap = (double) atof(argv[i+1]);
			if(strcmp(in_str, "-sites") == 0) {fall = 0; l = atoi(argv[i+1]); u = atoi(argv[i+2]);}
			if(strcmp(in_str, "-nout") == 0) nout = atoi(argv[i+1]);
			if(strcmp(in_str, "-2only") == 0) fl=2;
			if(strcmp(in_str, "-prefix") == 0) strcpy(prefix, argv[i+1]);
		}
	}
	if (ifp == NULL)
	{
		printf("Could not find seq file in command line.\n");
		printf("\n\nInput filename for seqs:\n\n");
		scanf("%s", fname);
		ifp = fopen(fname, "r");
	}
	if (ifp == NULL) 
	{
		nrerror("Error in opening sequence file");
	}

    fscanf(ifp,"%i %i %i", &nseq, &lseq, &hd);
	if (hd<1 || hd>2) 
	{
		printf("\n\nThe first line of seqs file should have N_Seqs L_Seq Haploid(1)/Diploid(2)\n");
		nrerror("Error: have not read haploid/diploid status");
	}

	if (inloc != NULL) {
		fscanf(inloc,"%i %lf %c", &i, &tlseq, &c);
		if (i != lseq) {printf("\n\nError: loc file does not match sequence file\n\n"); exit(1);}
		sloc = dvector(1,lseq);
		for (i=1;i<=lseq;i++) {
			fscanf(inloc,"%lf", &sloc[i]);
			printf("\nSite %4i : position %10f", i, sloc[i]);
			if (i>1 && sloc[i]<=sloc[i-1]) {nrerror("Error in loc file: SNPs not monotonically increasing");}
		}
		fclose(inloc);
	}

	else if (ask_questions==1){
		printf("\n\nUse existing file with SNP locations?(0=no, 1=yes):");
		scanf("%i", &i);
		if (i) {
			printf("\nInput filename with SNP locations:");
			scanf("%s", fname);
			inloc = fopen(fname, "r");
			fscanf(inloc,"%i %lf %c", &i, &tlseq, &c);
			if (i != lseq) {printf("\n\nError: loc file does not match sequence file\n\n"); exit(1);}
			sloc = dvector(1,lseq);
			for (i=1;i<=lseq;i++) {
				fscanf(inloc,"%lf", &sloc[i]);
				printf("\nSite %4i : position %10f", i, sloc[i]);
				if (i>1 && sloc[i]<=sloc[i-1]) {nrerror("Error in loc file: SNPs not monotonically increasing");}
			}
			/*fclose(inloc);*/
		}
	}

	if (nseq > SEQ_MAX) {printf("\n\nMore than max no. sequences: Using first %i for analysis\n\n", SEQ_MAX); nseq=SEQ_MAX;}
	printf("\n\nReading %i sequences of length %i bases .........\n", nseq, lseq);

	seqs = (int **) imatrix(1, nseq+1, 1, lseq);
	seqnames = (char **) cmatrix(1, nseq+1, 1, MAXNAME+1);

	if (read_fasta(seqs, ifp, nseq, lseq, seqnames)) printf("\nSequences read succesfully\n\n");
	/*	fclose(ifp);*/

	nall = imatrix(1, lseq, 1, 6);
	fsnp = ivector(1,lseq);
	allele_count(seqs, nseq, lseq, nall, 1, hd, prefix);
	for (i=1;i<=lseq;i++) {fsnp[i]=1;}
	
	strcpy(fname, prefix);
	printf("\n\nSegregating sites written to file	: %ssites.txt\n", fname); 
	ofp = fopen(strcat(fname, "sites.txt"), "w");
	strcpy(fname, prefix);
	printf("Locations of segregating sites to file	: %slocs.txt\n\n",fname);
	loc = fopen(strcat(fname, "locs.txt"), "w");

	if (!ofp || !loc) {
		printf("\n\n***Error: cannot open output files ***\n\n"); exit(1);
	}

	if (ask_questions)
	{
		printf("\n\nAll segregating sites or those with 2 alleles? (1/2):"); scanf("%i", &fl);
		if (fl==2) 
		{
		  printf("\n\nFrequency cutoff? (0 for none):");
		  scanf("%lf", &fcut);
		}
		printf("\n\nFrequency cut-off for missing data:");
		scanf("%lf", &fgap);

		printf("\n\nPrint all sites? (1/0) :");
		scanf("%i", &fall);
		if (!fall) {
			printf("\nFrom site : ");
			scanf("%i", &l);
			printf("\nTo site   : ");
			scanf("%i", &u);
		}
		printf("\n\nNumber of sequences to output (0=all, otherwise random selection): ");
		scanf("%i", &nout);
	}
	if(fall) {l=1; u=lseq;}

/*	printf("\n\nPrint every nth SNP:");
	scanf("%i", &nth);*/

	nth=1;


	if (!nout) nout=nseq;
	else if (nout>nseq) nout=nseq;
	else if (nout<0) nrerror("Cannot output negative # seqs!!");
	seq_out_list = ivector(1,nseq);
	for (i=1;i<=nseq;i++) seq_out_list[i]=1;
	for (i=1,out_ct=0;out_ct<nseq-nout;) {
		j = (int) ((double) 1+ran2()*nseq);
		if (seq_out_list[j]) {seq_out_list[j]=0; out_ct++;}
	}

	if (l<1) l=1;
	if (u>lseq) u=lseq;
	if (l>lseq || u<1) nrerror("Impossible limits to output");
	for (i=l, psite=0;i<=u;i++) {
		        for (j=2,na=0,nmin=(int) nseq*hd, jmin=1;j<=5;j++) {
			  if (nall[i][j]) na++; 
			  if (nall[i][j] && nall[i][j]<nmin) {nmin=nall[i][j]; jmin=j;}
			}
			seqs[nseq+1][i]=jmin;/*Defines minor allele at each locus*/
			if ((fl==1)&&(na>1)) nall[i][6]=1;
			else if (fl==2 && na==2 && nmin!=nseq*hd && nmin>((double) nseq*hd*fcut)) nall[i][6]=1;
			else nall[i][6]=0;
			if ( nall[i][1] > ((double) nseq*hd*fgap)) nall[i][6]=0;
			if (nall[i][6] && na>2) fsnp[i]=0;  /*fsnp determines whether output is 0/1 or TCAG*/
			if ((i-l)%nth) nall[i][6]=0;
		if (nall[i][6]) psite++;
	}
	
	if (psite==0) {printf("\n\nNo data to output\n\n"); exit(1);}
	fprintf(ofp,"%i %i %i",nout,psite,hd);
	fprintf(loc,"%i %.3f %c", psite, (inloc ? sloc[u]/*-sloc[l]+2*/ : (double) u-l+1), (inloc ? c : 'L'));

	for (i=1;i<=nseq;i++) if (seq_out_list[i]) 
	{
		fprintf(ofp,"\n>%s\n",seqnames[i]);
		for (j=l,na=0;j<=u;j++) if (nall[j][6]) 
		{
			if (hd==2) 
				fprintf(ofp,"%c",sc[seqs[i][j]]);
			else if (!fsnp[j]) fprintf(ofp,"%c",  bases[seqs[i][j]]);
			else if (fsnp[j] && seqs[i][j]==1) fprintf(ofp,"%c",sc[1]);
			else fprintf(ofp,"%i",seqs[i][j]==seqs[nseq+1][j]); /*compares to minor allele*/
			na++; 
			if ((na%50)==0) fprintf(ofp,"\n");
		}
	}
	fclose(ofp);
	for (i=l;i<=u;i++) if (nall[i][6]) fprintf(loc,"\n%.3f", (inloc ? sloc[i]/*-sloc[l]+1*/ : (double) i-l+1));
	fclose(loc);

	for (i=1, sn=0, pwd=0, ns=0; i<=lseq; i++) if (nall[i][6]) {
	  sn++; 
	  for (j=2; j<=5; j++) if (nall[i][j]==1) ns++;
	  for (j=2;j<=5;j++) pwd += (double) nall[i][j]*nall[i][j]/(pow(nseq*hd-nall[i][1],2));
	}
	pwd = (double) sn-pwd;
	pwd *= (double) nseq*hd/(nseq*hd-1);

	if (hd==1) for (i=1, spwd[0]=spwd[1]=0.0;i<nseq;i++) for(j=i+1;j<=nseq;j++) {
		for(site=1,d=0; site<=lseq;site++) if (seqs[i][site]!=seqs[j][site]) d++;
		spwd[0]+=(double) d;
		spwd[1]+=(double) d*d;
	}

	for (i=1, ctaj[0]=ctaj[1]=0.0; i<nseq*hd; i++) {ctaj[0]+=(double) 1/i; ctaj[1]+= (double) 1/(i*i);}
	n = (double) nseq*hd;
	ctaj[2]=(double) (n+1)/(3*n-3);
	ctaj[3]=(double) 2*(n*n+n+3)/(9*n*(n-1));
	ctaj[4]=ctaj[2]-1/ctaj[0];
	ctaj[5]=(double) ctaj[3]-(n+2)/(n*ctaj[0])+ctaj[1]/(ctaj[0]*ctaj[0]);
	ctaj[6]=ctaj[4]/ctaj[0];
	ctaj[7]=ctaj[5]/(ctaj[0]*ctaj[0]+ctaj[1]);
	ctaj[8]=(double) 2*(n*ctaj[0]-2*(n-1))/((n-1)*(n-2));
        ctaj[9]=ctaj[8]+(n-2)/((n-1)*(n-1))+((2/(n-1))*(1.5-(2*(ctaj[0]+1/n)-3.0)/(n-2)-1/n));
	ctaj[10]=(n*n*ctaj[1]/((n-1)*(n-1))+ctaj[0]*ctaj[0]*ctaj[9]-2*n*ctaj[0]*(ctaj[0]+1)/((n-1)*(n-1)))/(ctaj[0]*ctaj[0]+ctaj[1]);
	ctaj[11]=n/(n-1)*(ctaj[0]-n/(n-1))-ctaj[10];

	if (DEBUG) {
	  printf("\n\nCoefficients for tests\n\n");
	  for (i=0;i<12;i++) printf("%8.4f\n", ctaj[i]);
	}

	spwd[0]*=(double) 2/(n*(n-1));
	spwd[1]*=(double) 2/(n*(n-1));
	spwd[1]-=(double) spwd[0]*spwd[0];
	printf("\n\nSummary of output data\n");
	printf("\nSegregating sites      = %8.0f", sn);
	printf("\nAverage PWD            = %8.3f", pwd);
	printf("\nWatterson theta        = %8.3f", (double) sn/ctaj[0]);
	printf("\nTajima D statistic     = %8.3f", (double) (pwd-sn/ctaj[0])/sqrt(ctaj[6]*sn+ctaj[7]*sn*(sn-1)));
	printf("\nFu and Li D* statistic = %8.3f", (double) (n/(n-1)*sn-ctaj[0]*ns)/sqrt(ctaj[11]*sn+ctaj[10]*sn*sn));
	if (hd==1) printf("\nVariance PWD           = %8.3f", (double) spwd[1]);
	printf("\n\n");

	/*system("PAUSE");*/

}


void print_help(int argc, char* argv[]) 
{
	int i;
	char *in_str;
	for(i = 0; i < argc; i++){
		if(*argv[i] == '-'){ 
			in_str = argv[i];
			if((strcmp(in_str, "-h") == 0) || (strcmp(in_str, "-?") == 0) || (strcmp(in_str, "-help") == 0)) {
				printf("\nconvert\n");
				printf("Convert FASTA-style file to LDhat format.\n\n");
				printf("Required Options :\n");
				printf("-seq <file>          Input FASTA-style format file.\n");
				printf("\n\n");	
				printf("Additional Options :\n");
				printf("-loc <file>          SNP positions in seq file. Assumed contiguous if absent\n");
				printf("-2only               Only output sites with exactly two alleles\n");
				printf("-freqcut <float>     Min Minor Allele Frequency (between 0 and 1: default=0.0)\n");
				printf("-missfreqcut <float> Max Missing data frequency (between 0 and 1: default=1.0)\n");
				printf("-sites <int> <int>   Only print sites between these two values: default=all\n");
				printf("-nout <int>          Number of sequences to output: default=all\n");
				printf("-prefix <string>     Prefix of output files\n");
				printf("\n\n");	
				exit(0);
			}
		}
	}
}

