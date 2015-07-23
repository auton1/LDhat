#include "ldhat.h"
#include "tools.h"
#include "seqtools.h"
#include "snp_sim.h"

long *idum;
double *lnfac_array;
int sizeofpset=100;

void print_help(int argc, char *argv[]);

main (int argc, char *argv[]) {

	int i, j, **seqs, **nall, ord=1, ns, **pij, lkf=0, npt=0, pnew=0, anc=0;
	int tcat=1, rcat=0, verb=1, miss=0, *flocs;

	int sw_flag=0, moment_flag=0, rmin_flag=0, sim_flag=0, test_flag=0;
	char fname[MAXNAME+1], **seqnames;
	long seed=-setseed();
	extern int sizeofpset;
	double *locs;

	double **lkmat, *lkres;
	FILE *ifp=NULL, *ifp2=NULL, *ifp3=NULL, *tfp;
	struct site_type **pset;
	struct data_sum *data;
	int ask_questions = 1;
	char *in_str;

	print_help(argc, argv);
	idum = &seed;
	data = malloc((size_t) sizeof(struct data_sum));
	data->exact = 0;
	strcpy(data->prefix, "");

	for(i = 0; i < argc; i++)
	{
		if(*argv[i] == '-')
		{ 
			in_str = argv[i];
			ask_questions = 0;
			if(strcmp(in_str, "-seq") == 0) ifp = fopen(argv[i+1], "r");		
			if(strcmp(in_str, "-loc") == 0) ifp2 = fopen(argv[i+1], "r");
			if(strcmp(in_str, "-lk") == 0) 
			{
				lkf = 1;
				ifp3 = fopen(argv[i+1], "r");
			}
			if(strcmp(in_str, "-exact") == 0) data->exact = 1;
			if(strcmp(in_str, "-concise") == 0) verb=0;
			if(strcmp(in_str, "-window") == 0) sw_flag=1;
			if(strcmp(in_str, "-moment") == 0) moment_flag=1;
			if(strcmp(in_str, "-simulate") == 0) sim_flag=1;
			if(strcmp(in_str, "-rmin_flag") == 0) rmin_flag=2;
			if(strcmp(in_str, "-test") == 0) test_flag=1;
			if(strcmp(in_str, "-prefix") == 0) strcpy(data->prefix, argv[i+1]);
		}
	}
	if (ifp == NULL) 
	{
		printf("\nCould not find seqs file in command line.\n");
		printf("\nInput filename for seqs:\n");
		scanf("%s", fname);
		ifp = fopen(fname, "r");
	}
	if (ifp == NULL) nrerror("Error in opening sequence file");

	
	fscanf(ifp,"%i%i%i", &data->nseq, &data->lseq, &data->hd);
	if ((data->nseq < 2) || (data->lseq < 2)) {printf("\n\nInsufficient data for analysis (n > 1, L > 1) \n\n"); exit(1);}
	if (data->nseq > SEQ_MAX) {printf("\n\nMore than max no. sequences: Using first %i for analysis\n\n", SEQ_MAX); data->nseq=SEQ_MAX;}
	printf("\nAnalysing %i (n=%i) sequences of length %i seg sites\n", data->nseq, data->hd, data->lseq);
	seqs = imatrix(1, data->nseq, 1, data->lseq);
    seqnames = cmatrix(1, data->nseq+11, 1, MAXNAME+11);
	if (read_fasta(seqs, ifp, data->nseq, data->lseq, seqnames)) printf("\nSequences read successfully\n");
    fclose(ifp);

	nall = imatrix(1, data->lseq, 1, 6);
	allele_count(seqs, data->nseq, data->lseq, nall,1, data->hd, data->prefix);

	/*Store lnfac values in array for speed of computation*/

	lnfac_array = (double *) malloc((size_t) ((int) (data->nseq+2)*(data->hd))*sizeof(double));

	lnfac_array[0]=lnfac_array[1]=0;

	for (j=2;j<=((int) data->nseq*(data->hd));j++) lnfac_array[j]=(double) lnfac_array[j-1]+log(j);


	/*Open file with location of seg sites and read in data*/	
	if (ifp2 == NULL) 
	{
		printf("\nCould not find locs file in command line.\n");
		printf("\nInput name of file containing location of seg sites\n\n");
		scanf("%s", fname);
		ifp2 = fopen(fname, "r");
	}

	if (ifp2 == NULL) nrerror("Cannot open loc file");
	fscanf(ifp2, "%i %lf %c", &ns, &data->tlseq, &data->lc);
	if (ns != data->lseq) nrerror("Lseq and Locs disagree");
	if ((data->lc != 'C')&&(data->lc != 'L')) nrerror("Must input linear(L)/conversion(C)");
	if (data->lc == 'C') {
	  data->avc=0;
	  while (data->avc <= 0) {
	    printf("\n\nInput average tract length for conversion model: ");scanf("%lf", &(data->avc));
	  }
	}

	locs = dvector(1, data->lseq);
	flocs = ivector(1, data->lseq); /*Array to use when simulating data*/


	for (i=1; i<=data->lseq; i++) {
		fscanf(ifp2, "%lf", &locs[i]); 
		if ((locs[i]==0)||(locs[i]>data->tlseq)) {printf("\n\nError in Loc file\n\n%f\n", data->tlseq); exit(1);}
		if (i>1 && locs[i]<=locs[i-1]) nrerror("Error in locs file: SNPs must be montonically increasing");
	}
	printf("\nLocation of seg sites\n\n");
	for (i=1; i<=data->lseq; i++) printf("%3i   %4.2f\n", i, locs[i]);
	fclose(ifp2);

	/*Read in likelihood file where needed*/
    if (ask_questions) 
	{
			printf("\n\nUse existing likelihood file? (yes=1, no=0):");
			scanf("%i", &lkf);  /*lkf is a flag: 1 means use existing likelihood file as starting point*/
			if (lkf) 
			{
				printf("\n\nInput name of likelihood file: ");
				scanf("%s", fname);
				ifp3 = fopen(fname, "r");
			}
			else 
				data->exact=0;

			if (lkf == 1)
			{
				printf("\n\nIs likelihood file an exact match to data?(no=0/yes=1): ");
				scanf("%i", &data->exact);
			}
	}

	if (lkf && !ifp3) nrerror("Cannot open likelihood file");
	if (!lkf && data->hd==2) nrerror("For diploid data need complete lookup table for sequences");

	/*Store pair-types in pij matrix - classify in pair_spectrum routine*/

	data->w	= data->lseq;  /*Note for this program use all data - pair_int restricts to a smaller window*/
	pij = imatrix((int) 1,(int) data->lseq,(int) 1,(int) data->w);

	for (i=1;i<=data->lseq;i++) for (j=1;j<=data->w;j++) pij[i][j]=0;

	pset = init_pset(pset, lkf, ifp3, &npt, data);  /*Reads in type configurations from likelihood file*/

	printf("\n\n*** Calculating distribution of pair types ***\n\n");
	pset = pair_spectrum(seqs, data, nall, pset, &npt, &pnew, &miss, anc, pij);
	printf("\n\n *** Completed classification of pair types ***\n\n");

	if (data->exact && (pnew || miss)) nrerror("Lookup table is not exact for sequences\n(possibly generated by interval)");
	printf("\n\nOld = %i: New = %i: Missing = %i\n\n", npt,pnew,miss);
	data->ptt = (int) npt+pnew+miss;  /*npt is number from likelihood file, pnew is number new with no missing data, miss is # new with missing data*/
	if (verb) {
		strcpy(fname, data->prefix);
		tfp = fopen(strcat(fname, "type_table.txt"), "w");
		if (!tfp) nrerror("Cannot open type file");
		type_print(pij, data->lseq, data->w,tfp);
		fclose(tfp);
	}
	if (verb) print_pairs(stdout, pset, npt+pnew, data->hd, data->nseq);

	/*Need a complete set for missing data or diploid data - check this*/
	if (!data->exact && (data->hd ==2 || miss)) {
		printf("\n\nMissing data or diploid: checking that likelihood table is exhaustive\n\n");
		check_exhaustive(pset,npt,(data->nseq)*((int) data->hd));
	}
	/*Read parameters and likelihoods from likelihood file - where appropriate*/
	if (lkf) {
		read_pars(ifp3, &tcat, &data->th, &data->rcat, &data->rmax);
		lkmat = dmatrix(1,npt+pnew+miss,1,data->rcat);
		if (lkf) read_lk(ifp3, lkmat, npt, tcat, data->rcat);
	}

	/*If haploid, but novel types, need to calculate new likelihoods and input parameter values*/
	if (data->hd ==1 && pnew) { /*Note can have pnew for diploid data, but this has been checked for already*/
		if (!lkf) {
			data->th=data->rmax=-1.0; data->rcat=0;
			printf("\n\nInput theta per site (suggest Watterson estimate of %.5f):",(double) data->lseq/(watterson(data->nseq*data->hd)*data->tlseq));
			while (data->th<0.0) scanf("%lf", &data->th);
			printf("\n\nMax 4Ner for grid (suggest 100):");
			while(data->rmax<0.0) scanf("%lf", &data->rmax);
			printf("\n\nNumber of points on grid (suggest 101, min=2):");
			while(data->rcat<2) scanf("%i", &data->rcat);
			lkmat = dmatrix(1,npt+pnew+miss,1,data->rcat);
		}
		lk_est(pset,npt,pnew,lkmat,data->th,data->rcat,data->rmax);
		data->exact=1;
	}

	/*Sum over missing data or resolve genotypes and sum over missing data+configurations*/
	else if (miss && data->hd==1) {  
		printf("\n\n*** Calculating likelihoods for missing data ***\n\n");
		for (i=1;i<=miss;i++) {
			lk_miss(pset[npt+i],lkmat[npt+i],lkmat,data);
			printf("\rType %i", i);
		}

		printf("  ...Done!\n\n");
	}


	/*Sum over resolutions for diploid data*/
	else if (data->hd==2 && !data->exact) {
	  printf("\n\n*** Resolving diploid data: %i ***\n\n",pnew+miss);
	  lkres = dvector(1,data->rcat);
	  for (i=1;i<=pnew+miss;i++) {
	    lk_resolve(lkres,pset[npt+i],lkmat[npt+i],lkmat,data);
	    printf("\rType %i", i); 
	  }
	  free_dvector(lkres,1,data->rcat); 

	  printf("  ...Done!\n\n");
	}

	/*If new likelihood generated can output likelihood file for future analyses*/
	if (verb) print_lks(pset, data, npt+pnew+miss, lkmat);


	/*Basic analysis - estimation of 4Ner asuming constant rate*/

	data->rme=data->rmax; data->rce=data->rcat;
	if (1) {
		printf("\n\nDo you wish to change grid over which to estimate likelihoods for (default = %i points, 4Ner 0 - %.1f) (1/0) :",data->rcat,data->rmax);
		scanf("%i", &lkf);
		if (lkf) {
			data->rme=-10; data->rce=0;
			printf("\n\nMax 4Ner for estimation           : ");
			while (data->rme < 0.0) scanf("%lf", &data->rme);  
       		printf("\n\nNumber of classes to estimate for: ");
       		while (data->rce < 1) scanf("%i", &data->rce);
		}
	}
	data->lksurf = dmatrix(1,data->rce,1,2);
	lk_surf(pset, pij, data, lkmat, data->th, locs, 1);


	/*Print marginal likelihood ratio test statistics for each pair of sites*/
	printf("\n\nCalculating fits\n\n");
	fit_pwlk(data,pij,locs,lkmat,verb);

	/*Sliding windows version*/
	if (1) {
		printf("\n\nDo you wish to carry out a sliding windows analysis? (yes=1/no=0):");
		scanf("%i", &sw_flag);
	}
	if (sw_flag) lk_win(pset,pij,data,lkmat,locs,nall);

	/*Nonparametric estimation of recombination rate*/
	if (1) {
		printf("\n\nPrint out table of Rmin values?\n(0=No, 1=Total only, 2=Full table):");
		scanf("%i", &rmin_flag);
	}

	if (rmin_flag) {
		rmin(data, pset, pij, locs, lkf-1);
		printf("\n\nLower bound on Rmin = %i\n\n",data->rmin);
	}

	/*Estimate 4Ner by Wakeley 1997 method*/
	if (1) {
		printf("\n\nEstimate 4Ner by moment method? (yes=1, no=0)");
		scanf("%i", &moment_flag);
	}

	if (moment_flag) wakeley_est(data, seqs, locs);

	/*Recombination tests - only available for haploid data!*/
	if (data->hd==1) {
		if (1) {
			printf("\n\nDo you wish to test for recombination? (yes=1, no=0): ");
			scanf("%i", &test_flag);
		}
		if (test_flag) {
			rec_test(data, pij, locs, lkmat, pset, npt+pnew+miss);
		}
	}

	/*Conditional simulation - only available for haploid data with a complete lk file*/
	if (data->hd==1 && !(data->exact)) {

		if (1) {
	  printf("\n\nDo you wish to test constant-rate model and estimate sampling distribution by simulation? (yes=1/no=0): ");
	  scanf("%i", &test_flag);
		}
	  if (test_flag) {
	    freq_min(locs, flocs, nall, data);
	    printf("\n\nHow many simulations? ");
	    scanf("%i", &lkf);
	    snp_sim(locs, flocs, pset, lkmat, lkf, data);
	  }
	}

	free_imatrix(pij,1,data->lseq,1,data->w);
	free_imatrix(seqs,1,data->nseq,1,data->lseq);
	free_imatrix(nall,1,data->lseq,1,5);
	for (i=1;i<sizeofpset;i++) free(pset[i]);
	free(pset);
	free(data);
	free_dvector(locs, 1, data->lseq);
	free_ivector(flocs, 1, data->lseq);

	/* system("PAUSE"); */
}


void freq_min(locs, flocs, nall, data)
     int *flocs, **nall;
     double *locs;
     struct data_sum *data;
{
  int i, j, fm, del=0;

  for (i=1;i<=data->lseq;i++) {
    for (j=2, fm=data->nseq;j<6;j++) if (nall[i][j] && nall[i][j]<fm) fm=nall[i][j];
	flocs[i]=fm;
  }
  for (i=1;i<=data->lseq;i++) {
	  if (flocs[i]==data->nseq) {
		  del++;
		  for (j=i;j<data->lseq;j++) {locs[j]=locs[j+1]; flocs[j]=flocs[j+1];}
		  i--;
	  }
  }
  data->lseq-=del;

  if (DEBUG) {
	  printf("\n\nLocs and flocs\n\n");
	for (i=1;i<=data->lseq;i++) printf("%10.2f %5i\n", locs[i], flocs[i]);
  }
}


/*Sliding windows version of same program*/
void lk_win(pset,pij,data,lkmat,locs,nall)
     struct site_type **pset;
     int **pij, **nall;
     double *locs;
	 double **lkmat;
     struct data_sum *data;
{
  int win, nwin, lwin, i, j, l, u, cat, nmin;
  double rhoe, rhoi, tc[8], avpwd;
  double lkmax, lkrun, lkg;
  char fname[MAXNAME+1];
  FILE *wfp;

  printf("\n\nInput number of windows: ");
  scanf("%i", &nwin);
  printf("\n\nInput length of windows in bp (total length = %.2f): ",data->tlseq);
  scanf("%i", &lwin);

  /*Calculate constants for Tajima D statistic: currently not used*/  
  for (i=2, tc[0]=tc[1]=1.0; i<data->nseq*data->hd; i++) {tc[0] += (double) 1/i; tc[1] += (double) 1/(i*i);}
  tc[2] = (double) (data->nseq*data->hd+1)/(3*(data->nseq*data->hd-1));
  tc[3] = (double) 2*((data->nseq*data->hd)*(data->nseq*data->hd)+(data->nseq*data->hd)+3)/(9*(data->nseq*data->hd)*(data->nseq*data->hd-1));
  tc[4] = (double) tc[2]-1/tc[0];
  tc[5] = (double) tc[3]-(data->nseq*data->hd+2)/(tc[0]*data->nseq*data->hd)+tc[1]/(tc[0]*tc[0]);
  tc[6] = (double) tc[4]/tc[0];
  tc[7] = (double) tc[5]/(tc[0]*tc[0]+tc[1]);

  strcpy(fname, data->prefix);
  wfp = fopen(strcat(fname,"window_out.txt"), "w");
  if (!wfp) nrerror("Cannot open windows outfile");
  fprintf(wfp,"\nSliding windows analysis - rho for total gene = %.5f per bp/kb\n\nSNP_L\t SNP_R\t #SNPs\t4Ner/bp/kb\t CLR\t  Tajima's D\n\n",\
	  (double) data->rho/data->tlseq);

  for (win=1;win<=nwin;win++) {
    for (l=1;locs[l]<(win-1)*(data->tlseq-lwin)/(nwin-1);l++);
    for (u=l;locs[u]<(win-1)*(data->tlseq-lwin)/(nwin-1)+lwin;u++) if (u>data->lseq) break; 
    u--;
    printf("\nWindow %3i: points: %i - %i",win,l,u);
    if (u!=l) {
      for (cat=1, lkmax=0.0, rhoi=0.0; cat<=data->rce;cat++) {
		rhoe = (double) (cat-1)*(data->rme)/(data->rce-1);  /*Input value of rho for whole window*/
		if (!lk_calc_win(pij,l,u,data,&lkrun,locs,rhoe,data->rho,lkmat)) nrerror("Error in calculating likelihoods");
		if (rhoe==data->rho) lkg = lkrun;
		if (cat==1) {lkmax=lkrun; rhoi=rhoe;}
		else if (lkrun>lkmax) {lkmax=lkrun; rhoi=rhoe;}
      }

 	  for (i=l,avpwd=0.0;i<=u;i++) {
		  for (j=2,nmin=((int) data->nseq*data->hd);j<=5;j++) if (nall[i][j] && nmin>nall[i][j]) nmin=nall[i][j];
		  avpwd += (double) nmin*(data->nseq*data->hd-nmin);
	  }

	  avpwd *= (double) 2/((data->nseq*data->hd)*(data->nseq*data->hd-1));


      fprintf(wfp,"%10.2f\t%10.2f\t%4i\t%10.5f\t%10.2f\t%10.3f\n", locs[l], locs[u], u-l+1, (double)rhoi/(locs[u]-locs[l]), lkmax-lkg, \
		  (double) (avpwd-(u-l+1)/tc[0])/sqrt(tc[6]*(u-l+1)+tc[7]*(u-l+1)*(u-l)));
    }
    else if (u==l) fprintf(wfp,"%10.2f - %10.2f : no data\n", locs[l], locs[u]);
    else nrerror("Error in window routine");
  }
  fclose(wfp);
   
}

/*Calculate likelihoods for diploid data*/
void lk_resolve(lkres,pset,lknew,lkmat,data)
     struct site_type *pset;
     struct data_sum *data;
     double *lkres,**lkmat, *lknew;
{
  int i, j, fl, pbase[9], p1, p2, hap;
  double mn;
  struct site_type tres;

  for (i=1;i<=data->rcat;i++) lkres[i]=lknew[i]=0.0;

  pbase[0]=2*pset->pt[5]+pset->pt[7]+pset->pt[13];
  pbase[1]=2*pset->pt[9]+pset->pt[11]+pset->pt[13];
  pbase[2]=2*pset->pt[6]+pset->pt[7]+pset->pt[14];
  pbase[3]=2*pset->pt[10]+pset->pt[11]+pset->pt[14];
  pbase[4]=2*pset->pt[4]+pset->pt[12];
  pbase[5]=2*pset->pt[8]+pset->pt[12];
  pbase[6]=2*pset->pt[1]+pset->pt[3];
  pbase[7]=2*pset->pt[2]+pset->pt[3];
  pbase[8]=2*pset->pt[0];
  fl = pbase[4]+pbase[5]+pbase[6]+pbase[7]+pbase[8];

  printf("\nResolving GT: "); for (i=0;i<16;i++) printf("%2i ",pset->pt[i]);

  for (i=0;i<=pset->pt[15];i++) {
    mn = (double) lnfac_array[pset->pt[15]]-lnfac_array[i]-lnfac_array[pset->pt[15]-i];
    for (j=0;j<9;j++) tres.pt[j]=pbase[j];
    tres.pt[0]+=i;
    tres.pt[3]+=i;
    tres.pt[1]+=pset->pt[15]-i;
    tres.pt[2]+=pset->pt[15]-i;
    order_pt_hap(tres.pt,data->nseq*2);
    if (DEBUG) {printf("\nRes : "); for (j=0;j<9;j++) printf("%2i ",tres.pt[j]);}
    if (fl) lk_miss(&tres,lkres,lkmat,data);
    else {
      p1 = (int) tres.pt[1]+tres.pt[3];
      p2 = (int) tres.pt[2]+tres.pt[3];
      hap = (int) p1*(p1-1)*(p1+4)/6+(p2-1)*(p2+2)/2+tres.pt[2]+1;
      for (j=1;j<=data->rcat;j++) lkres[j]=lkmat[hap][j];
    }
    if (!i) for (j=1;j<=data->rcat;j++) lknew[j]=lkres[j];
    else for (j=1;j<=data->rcat;j++) lknew[j]+= (double) log(1+exp((double) mn+lkres[j]-lknew[j]));
  }

  if (DEBUG) {printf("\nSummed likelihood: "); for (j=1;j<=data->rcat;j++) printf("%7.2f ",lknew[j]);}
}




/*Calculate likelihoods for missing data*/

void lk_miss(pset,lkmiss,lkmat,data)
		struct site_type *pset;
		struct data_sum *data;
		double **lkmat, *lkmiss;
{

  int j, a, b, c, d, e1, e2, e3, e4, pres[9], p1, p2, ct, k, ht;
  double cf, mn;

  if (DEBUG) {printf("\nHap: "); for (j=0;j<9;j++) printf("|%i|",pset->pt[j]);}
  cf = 0.0; for (j=1;j<=data->rcat;j++) lkmiss[j]=0.0;
  for (j=0;j<9;j++) pres[j]=0;
  for (k=1;k<=data->rcat;k++) lkmiss[k]=0.0;
  ct=0;
  for (a=0;a<=pset->pt[4];a++)
      for (b=0;b<=pset->pt[5];b++)
	for (c=0;c<=pset->pt[6];c++)
	  for (d=0;d<=pset->pt[7];d++)
	    for (e1=0;e1<=pset->pt[8];e1++)
	      for (e2=0;e2<=pset->pt[8]-e1;e2++)
		for (e3=0;e3<=pset->pt[8]-e1-e2;e3++) {
		  ct++;
		  e4 = pset->pt[8]-e1-e2-e3;
		  pres[0]=pset->pt[0]+a+c+e1;
		  pres[1]=pset->pt[1]+b+pset->pt[6]-c+e2;
		  pres[2]=pset->pt[2]+pset->pt[4]-a+d+e3;
		  pres[3]=pset->pt[3]+pset->pt[5]-b+pset->pt[7]-d+e4;
		  /*		  printf("\n%3i %3i %3i %3i: %i",pres[0], pres[1], pres[2], pres[3], data->nseq*((int) 1 + (data->hd==1?0:1)));*/
		  order_pt_hap(pres,(int) data->nseq*((int) 1 + (data->hd==1?0:1)));
		  p1 = (int) pres[1]+pres[3];
		  p2 = (int) pres[2]+pres[3];
		  ht = (int) p1*(p1-1)*(p1+4)/6+(p2-1)*(p2+2)/2+pres[2]+1;
		  /*		  printf("\n%3i %3i %3i %3i: %i %i %i",pres[0], pres[1], pres[2], pres[3], p1>data->nseq?1:0,p2>data->nseq?1:0,p2>p1?1:0 );
				  printf("\n ht = %i",ht);*/
		  if (DEBUG) {printf("\nRes %i: ht %i: ",ct,ht); for (k=0;k<4;k++) printf("|%i|",pres[k]);}
/*		  mn = (double) lognC2(pset->pt[4],a)+lognC2(pset->pt[5],b)+lognC2(pset->pt[6],c)+lognC2(pset->pt[7],d)+lognC4(pset->pt[8],e1,e2,e3,e4);*/
		  mn = (double) lnfac_array[pset->pt[4]]-lnfac_array[a]-lnfac_array[pset->pt[4]-a]+lnfac_array[pset->pt[5]]-lnfac_array[b]-lnfac_array[pset->pt[5]-b]+\
			  lnfac_array[pset->pt[6]]-lnfac_array[c]-lnfac_array[pset->pt[6]-c]+lnfac_array[pset->pt[7]]-lnfac_array[d]-lnfac_array[pset->pt[7]-d]+\
			  lnfac_array[pset->pt[8]]-lnfac_array[e1]-lnfac_array[e2]-lnfac_array[e3]-lnfac_array[e4];
		  if (DEBUG) printf("\nlog(mn coff) = %.2f\n",mn);

		  /*Better routine for summation*/
		  if (ct==1) {for (k=1;k<=data->rcat;k++) lkmiss[k] = lkmat[ht][k];}
		  else {for (k=1;k<=data->rcat;k++) lkmiss[k] += (double) log(1+exp(lkmat[ht][k]+mn-lkmiss[k]));}
		}

    if (DEBUG) {printf("\nLk : "); for (k=1;k<=data->rcat;k++) printf(" %.2f", lkmiss[k]);}
}



/*Check that likelihood file is exhaustive for n*/

void check_exhaustive(pset,npt,nsamp)
     struct site_type **pset;
     int npt, nsamp;
{
  int p1, p2, i, ei;

  printf("\n\nChecking likelihood file is exhaustive:...");

  p1 = (int) nsamp/2;
  if (npt != (int) 1+p1+p1*(p1-1)*(p1+4)/6+(p1-1)*(p1+2)/2) {
    printf("\n\n!!npt = %i: E[npt] = %i\n\n",npt, 1+p1+p1*(p1-1)*(p1+4)/6+(p1-1)*(p1+2)/2);
    nrerror("npt and total do not agree: not exhaustive: please use lkgen program to generate a complete lookup table for sample size");
  }

  for (i=1;i<=npt;i++) {
    p1 = pset[i]->pt[1]+pset[i]->pt[3];
    p2 = pset[i]->pt[2]+pset[i]->pt[3];
    ei = (int) p1*(p1-1)*(p1+4)/6+(p2-1)*(p2+2)/2+pset[i]->pt[2]+1;
    if (i!=ei) {printf("\nError in pair-types: p1(%i), p2(%i), n11(%i): i(%i) ei(%i)\n\n",p1,p2,pset[i]->pt[3],i,ei); exit(1);}
  }
  printf("OK\n\n");
}


/*Likelihood estimation for each pairwise comparison using the method of Fearnhead and Donnelly (2001) */

void lk_est(pset,npt,pnew,lkmat,stheta,rcat,rmax) 
int npt, pnew, rcat;
double  stheta, rmax;

double **lkmat;
struct site_type **pset;
{

	int i, j, k, p, nd, *data, K;
	double **P, *mu, theta[2],  *log_lk;

	printf("\n\nEstimating likelihood surface for data (%i pair types)\n", pnew);

	K = 2;
	mu = (double *) calloc(K,sizeof(double));
	for (i=0; i<K; i++) mu[i] = (double) 1/K;
	P = (double **) calloc(K, sizeof(double *));
	for (i=0; i<K; i++) P[i] = (double *) calloc(K, sizeof(double));
	for (i=0; i<K; i++) for (j=0; j<K; j++) {
		if (j!=i) P[i][j] = (double) 1/(K-1);
		else P[i][j] = 0.0;
	}
	log_lk = (double *) calloc(rcat, sizeof(double));

	for (i=npt+1; i<=npt+pnew; ) if (!(pset[i]->miss)) {
		for (j=0, nd=0; j<4; j++) if (pset[i]->pt[j]) nd++;
		data = (int *) malloc((size_t) 3*nd*sizeof(int)); 
		for (j=0,k=0; j<4; j++) 
		   if (pset[i]->pt[j]) {
			data[3*k] = (int) j/2 + 1;
			data[3*k+1] = (int) j%2 + 1;
			data[3*k+2] = (int) pset[i]->pt[j];
			k++;
		   }
		printf("\nHaplotype data: (%i of %i)\n",i,npt+pnew); for (j=0; j<nd; j++) printf("%i%i:%i\n",data[3*j], data[3*j+1],data[3*j+2]);
		theta[0]=stheta;
		theta[1]=stheta;
		printf("theta : %.3f %.3f\n", theta[0], theta[1]);
		pairs(nd, data, K, P, mu, theta, (double) rmax, rcat, NRUN, log_lk);
		for (p=1; p<=rcat; p++) lkmat[i][p] = (double) log_lk[p-1];
		free(data);
		i++;
	}

	free(mu);
	free(log_lk);
	for (i=0;i<K;i++) free(P[i]);
	free(P);
}	


void print_lks(pset,data,npt,lkmat) 
int npt;
double **lkmat;
struct site_type **pset;
struct data_sum *data;
{

	int p, i, nstate, ct=0;
	char fname[MAXNAME+1];
	FILE *ofp;

	strcpy(fname, data->prefix);
	ofp = fopen(strcat(fname, "new_lk.txt"), "w");
	if (data->hd == 2) nstate=16;
	else nstate=9;

	for (p=1;p<=npt;p++) if (pset[p]->nt>0) ct++;

	fprintf(ofp, "\n%i %i\n1 %.5f\n%i %f\n\n",(data->nseq*data->hd),ct,data->th,data->rcat,data->rmax);
	for (p=1; p<=npt; p++) if (pset[p]->nt>0) {
		fprintf(ofp,"\n%4i # ", p);
		for (i=0; i<nstate; i++) fprintf(ofp,"%3i ", pset[p]->pt[i]);
		fprintf(ofp," :  ");
		for (i=1; i<=data->rcat; i++) fprintf(ofp,"%8.3f ", lkmat[p][i]);
	}
	fclose(ofp);
}





/*Estimation of the likelihood surface for rho*/
void lk_surf(pset,pij,data,lkmat,theta,locs, ff) 
int **pij, ff;
double theta, *locs;
double **lkmat;
struct site_type **pset;
struct data_sum *data;
{

	int i, *pars, fl=0, j, rho_i;
	double **clist, lkmax;
	char fname[MAXNAME+1];
	FILE *ofp;

	clist = dmatrix(1,data->rce, 1, 5);
	pars = ivector(1,data->lseq); for (i=1;i<=data->lseq;i++) pars[i]=1;

	for (j=1; j<=data->rce; j++) {
	   	if (data->rce > 1) clist[j][1] = (double) (data->rme)*(j-1)/((data->rce)-1);
		if (ff) printf("\nLikelihood estimation with 4Ner = %6.2f: ", clist[j][1]);
	   	if ((fl=lk_calc(pij, 1, data->lseq, data, &lkmax, locs, clist[j][1], lkmat)) == 0) lkmax=0;
		clist[j][2]=lkmax;
		clist[j][3]=(double) fl;
		if (ff) printf(" %.3f", lkmax);
		if (j==1 || ((lkmax < 0) && (lkmax > data->lkmax))) {
			data->lkmax = lkmax;
			data->rho=(double) clist[j][1];
			rho_i = j;
		}
	}
	if (ff) {
		printf("\n\nResults of Estimation\n\n 4Ner(region)       Lkmax     Type\n\n");
		for (j=1; j<=data->rce; j++) {
			printf("%8.2f   ",clist[j][1]);
			if (clist[j][2] < 0.0) printf("%8.3f  %3.0f", clist[j][2], clist[j][3]);
			if (clist[j][1] > data->rmax) printf(" Beyond range of estimated likelihoods");
			printf("\n");
			data->lksurf[j][1]=clist[j][1];
			data->lksurf[j][2]=clist[j][2];
		}
		printf("\n\nMaximum at 4Ner(region) = %.3f : Lk = %.3f \n\n",data->rho, data->lkmax);
		strcpy(fname, data->prefix);
		ofp = fopen(strcat(fname, "outfile.txt"), "w");
		if (ofp == NULL) {printf("\n\nCannot open outfile\n\n"); exit(1);}
		fprintf(ofp,"Lk surface\n\nTheta = %.5f\n\n", theta);
		fprintf(ofp,"Maximum at 4Ner(region) = %.3f : Lk = %.3f\n\n", data->rho, data->lkmax);
		fprintf(ofp,"\n\n4Ner(region)   Pairwise Lk\n\n");
		for (j=1;j<=data->rce;j++) {
			fprintf(ofp,"%8.2f  ", clist[j][1]);
			if (clist[j][2] < 0.0) fprintf(ofp,"%8.3f", clist[j][2]);
			if (clist[j][1] > data->rmax) fprintf(ofp," Beyond range of likelihoods: Using 4Ner_max as estimate");
			fprintf(ofp,"\n");
		}
		fclose(ofp);
	}

	if (ff) {
		data->rho_i = rho_i;
	}
	else {
		data->clr = (double) clist[rho_i][2]-clist[data->rho_i][2];
	}
	free_ivector(pars,1,data->lseq);
    free_dmatrix(clist,1,data->rce,1,5);
}
		

/*Routine to calculate the pairwise likelihood for any given rho*/

int lk_calc(pij,l,u,data,lkrun,locs,ct,lkmat) 
int **pij, l,u;
double ct, *locs;
double **lkmat, *lkrun;
struct data_sum *data;
{

	int i, j, k, t, fl=1;
	double cij, d, lke=-100, dij;

	for (i=l, (*lkrun)=0.0; i<u; i++)
                for (j=i+1; j<=u; j++) {
                   t = pij[i][j-i];
		   if (t > 0) {
			if (data->lc == 'C') {
			  dij = (double) locs[j]-locs[i];
			  /*If wish to use circular model, have to comment out previous line, and use next*/
			  /*dij = mini(locs[j]-locs[i], locs[i]+(data->tlseq)-locs[j]);*/
			  cij = (double) 2*ct*(1-exp((double) -dij/(data->avc)));
			}
			else cij = (double) ct*(locs[j]-locs[i])/(data->tlseq);
			if (cij < 0) {printf(" C < 0!! "); return 0;}
			if (cij > data->rmax) {
				lke = (double) lkmat[t][data->rcat];
				fl=2;
			}
			else {
			   d = (double) cij*(data->rcat-1)/(data->rmax);
			   k = (int) d+1;
			   if (k<data->rcat) lke = (double) lkmat[t][k]+(d-k+1)*(lkmat[t][k+1]-lkmat[t][k]);
			   else lke = lkmat[t][data->rcat];
			}
			if (lke < 0.0) { (*lkrun) += (double) lke;}
			else {
				printf(" Lk >= 0!! "); return 0;
			}
		   }
                }
	return fl;
}




/*Routine to calculate the pairwise likelihood for sliding windows*/
/*Allows a window (from l to u) with a rate (rw) that is different from background(rb)*/

int lk_calc_win(pij,l,u,data,lkrun,locs,rw,rb,lkmat) 
int **pij, l,u;
double  rw, rb, *locs;
double *lkrun, **lkmat;
struct data_sum *data;
{

	int i, j, k, t, fl=1;
	double cij, d, lke=-100, dij;

	/*Input values as rho for the whole gene - rescale to rates*/
	rb /= (double) data->tlseq;
	rw /= (double) (locs[u]-locs[l]);

	for (i=1, (*lkrun)=0.0; i<u; i++) {
                for (j=maxi(i+1,l); j<=data->lseq; j++) {
                   t = pij[i][j-i];
		   if (t > 0) {
			if (data->lc == 'C') {
			  dij = (double) locs[j]-locs[i];
			  /*If wish to use circular model, have to comment out previous line, and use next*/
			  /*dij = mini(locs[j]-locs[i], locs[i]+(data->tlseq)-locs[j]);*/
			  cij = (double) 2*rw*(1-exp((double) -dij/(data->avc)));
			}
			else {
			  if (j<l || i>=u) cij = (double) (locs[j]-locs[i])*rb;
			  else if (i<l) {
			    if (j>u) cij = (double) (locs[u]-locs[l])*rw+(locs[l]-locs[i]+locs[j]-locs[u])*rb;
			    else cij = (double) (locs[j]-locs[l])*rw + (locs[l]-locs[i])*rb;
			  }
			  else {
			    if (j>u) cij = (double) (locs[u]-locs[i])*rw+(locs[j]-locs[u])*rb;
			    else cij = (double) (locs[j]-locs[i])*rw;
			  }
			}
			if (cij < 0) {printf(" C < 0!! "); return 0;}
			if (cij > data->rmax) {
				lke = (double) lkmat[t][data->rcat];
				fl=2;
			}
			else {
			   d = (double) cij*(data->rcat-1)/(data->rmax);
			   k = (int) d+1;
			   if (k<data->rcat) lke = (double) lkmat[t][k]+(d-k+1)*(lkmat[t][k+1]-lkmat[t][k]);
			   else lke = lkmat[t][data->rcat];
			}
			if (lke < 0.0) { (*lkrun) += (double) lke;}
			else {
				printf(" Lk >= 0!! "); return 0;
			}
		   }
                }
	}
	return fl;
}

void rec_test(data,pij,locs,lkmat,pset,npt) 
int **pij, npt;
double *locs;
double **lkmat;
struct data_sum *data;
struct site_type **pset;
{

	int i, j, **pijs, *ord, ngs[4], shuff, tmp, *pars, imax, l, u, t, k, *anal;
	double f[2], npairs;
	double lds[3], dav[6], lmax, lkrun;
	char fname[MAXNAME+1];
	FILE *ofp;
	
	printf("\n\nTesting for recombination and analysing patterns of linkage disequilibrium\n\n");

	for (i=1; i<=npt; i++){
		f[0]=(double) ((pset[i]->pt[1])+(pset[i]->pt[3]))/(data->nseq);
		f[1]=(double) ((pset[i]->pt[2])+(pset[i]->pt[3]))/(data->nseq);
		pset[i]->ld_stat[0]=(double) (pset[i]->pt[3])/(data->nseq)-f[0]*f[1];
		if (f[0]&&f[1]&&(1-f[0])&&(1-f[1])) {
		  pset[i]->ld_stat[1]=(double) pow((pset[i]->ld_stat[0]),2)/(f[0]*f[1]*(1-f[0])*(1-f[1]));
		  if (pset[i]->ld_stat[0] < 0) pset[i]->ld_stat[2]=(pset[i]->ld_stat[0])/(-1*f[0]*f[1]);
		  else pset[i]->ld_stat[2]=(pset[i]->ld_stat[0])/((double) f[0]<f[1] ? f[0]*(1-f[1]) : f[1]*(1-f[0]));
		}
		else pset[i]->ld_stat[1]=pset[i]->ld_stat[2]=0.0;
	}

	pijs = imatrix(1, data->lseq, 1, data->lseq);
	ord = ivector(1, data->lseq);
	pars = ivector(1, data->lseq);
	anal = ivector(1,data->lseq);

	for (i=1; i<=data->lseq; i++) ord[i]=i;
	for (i=1;i<data->lseq;i++) {
		for (j=i+1,t=0;j<=data->lseq;j++) if (pij[i][j-i]) t++;
		if (t) anal[i]=1; else anal[i]=0;
	}
	if (pij[(data->lseq)-1][1]) anal[data->lseq]=1; else anal[data->lseq]=0;

	for (i=0;i<6;i++) dav[i]=0.0;	
	for (i=1, npairs = 0;i<data->lseq;i++) for (j=i+1; j<=data->lseq;j++) 
		if (pij[i][j-i] != 0) {
			npairs++;
			if (data->lc=='L') {
				dav[0] += (double) locs[j]-locs[i];
				dav[3] += (double) (locs[j]-locs[i])*(locs[j]-locs[i]);
			}
			else {
				dav[0] += (double) minc(locs[i], locs[j], data->tlseq);
				dav[3] += (double) pow((double) minc(locs[i], locs[j], data->tlseq), 2);
			}
			dav[1] += (double) pset[pij[i][j-i]]->ld_stat[1];
			dav[2] += (double) pset[pij[i][j-i]]->ld_stat[2];
			dav[4] += (double) pow((pset[pij[i][j-i]]->ld_stat[1]), 2);
			dav[5] += (double) pow((pset[pij[i][j-i]]->ld_stat[2]), 2);
		}
	for(i=0; i<6; i++) dav[i] /= (double) npairs;
	for(i=3;i<6;i++) {dav[i]-=pow(dav[i-3],2); dav[i]=sqrt(dav[i]);}
	ld_calc(pset, pij, locs, data->ld, data);

	printf("\n\nObserved LD statistics\n\n");
	printf("	Npairs       = %8.0f\n", npairs);
	printf("    Lkmax        = %8.2f\n", data->lkmax);
	printf("	G4           = %8.0f\n", data->ld[0]);
	printf("	corr(r2,D)   = %8.5f\n", ((data->ld[1])/npairs-dav[0]*dav[1])/(dav[3]*dav[4]));
	if (dav[5]>0.0001) 
	printf("	corr(D',d)   = %8.5f\n", ((data->ld[2])/npairs-dav[0]*dav[2])/(dav[3]*dav[5]));
	else
	printf("	corr(D',d)   = na\n");	

	strcpy(fname, data->prefix);
	ofp = fopen(strcat(fname, "rdist.txt"), "w");
	if (!ofp) {printf("\nCannot open file for output\n\n"); exit(1);}
	fprintf(ofp,"Distribution of estimated Rho values\n\nShuffle   Rho    Lkmax      G4     corr(r2,d)   corr(D',d)\n\n");

	for (i=0;i<4;i++) ngs[i]=0;
	for (shuff=1; shuff<=NSHUFF;shuff++) {
			if ((shuff%100)==0) printf("\nShuffle %i", shuff);
			for (i=1; i<=data->lseq;i++) if (anal[i]) {
				k=0;
				while (k==0) {j=(int) ((double) ran2()*(data->lseq))+1; if (anal[j]) k=1;}
				tmp=ord[i]; ord[i]=ord[j]; ord[j]=tmp;
			}
			for (i=1; i<data->lseq; i++) 
				for (j=i+1;j<=data->lseq;j++) {
					l = mini(ord[i], ord[j]);
					u = maxi(ord[i], ord[j]);
					pijs[i][j-i] = pij[l][u-l];
				}
			for (i=1; i<=data->lseq;i++) pars[i]=1;
	       			for (i=1, lmax=-100000000;i<=data->rce;i++) {
				if(lk_calc(pijs, 1, data->lseq, data, &lkrun, locs, (double) (i-1)*(data->rme)/((data->rce)-1), lkmat))
				  if (lkrun > lmax) {imax = i; lmax=lkrun;}
			}
			
			ld_calc(pset, pijs, locs, lds, data);
			fprintf(ofp,"%4i  %8.3f  %8.3f  %.0f  %8.5f  %8.5f\n",shuff, (double) (imax-1)*(data->rme)/((data->rce)-1), lmax, \
				lds[0], (lds[1]/npairs-dav[0]*dav[1])/(dav[3]*dav[4]), \
				dav[5]>0.0001 ? (lds[2]/npairs-dav[0]*dav[2])/(dav[3]*dav[5]) : 0.00);
			if (lmax > data->lkmax) ngs[0]++;
			for (i=1; i<3; i++) if (lds[i]<=(data->ld[i])) ngs[i+1]++;
			if (lds[0]>=(data->ld[0])) ngs[1]++;
		}

	fclose(ofp);
	printf("\n\nResults of shuffling\n\n");
	printf("Proportion Lkmax <= estimated      = %.4f\n", (double) ngs[0]/NSHUFF);
	printf("Proportion G4 >= G4 estimated      = %.4f\n",(double) ngs[1]/NSHUFF);
	printf("Proportion corr(r2,d) <= estimated = %.4f\n", (double) ngs[2]/NSHUFF);
	printf("Proportion corr(D',d) <= estimated = %.4f\n", (double) ngs[3]/NSHUFF);

	free_imatrix(pijs,1,data->lseq, 1, data->lseq);
	free_ivector(ord, 1, data->lseq);
	free_ivector(pars, 1, data->lseq);
	free_ivector(anal,1,data->lseq);
}



void ld_calc(pset,pijs,locs,ldv, data) 
int **pijs;
double *locs;
double ldv[3];
struct site_type **pset;
struct data_sum *data;
{

	int i, j, pt;
	double dij;
	
	for (i=0; i<3; i++) ldv[i]=0.0;
	for (i=1; i<data->lseq; i++)
		for (j=i+1; j<=data->lseq;j++) {
			pt = pijs[i][j-i];
			if (pt) {
				if (data->lc=='L') dij = locs[j]-locs[i];
				else dij = minc(locs[i], locs[j], data->tlseq);
			 	if (pset[pt]->ld_stat[2] < 0.9999) ldv[0]+=(double) dij;
				ldv[1] += (double) pset[pt]->ld_stat[1]*dij;
				ldv[2] += (double) pset[pt]->ld_stat[2]*dij;
			}
		}
}


void fit_pwlk(data,pij,locs,lkmat, fl) 
int **pij,fl;
double *locs;
double **lkmat;
struct data_sum *data;
{
	int i, j, t, k;
	double cij, d, rmp;

	double lkmax, clk=0.0, lke, dij;
	char fname[MAXNAME+1];
	FILE *ofp;
	
	if (fl) {
		strcpy(fname, data->prefix);
		ofp = fopen(strcat(fname, "fit.txt"),"w");
		fprintf(ofp,"\nFit between data and ML model for each pair (+ve: Rho_max_pair > Rho_est, -ve Rho_max_pair < Rho_est)\n\n      ");
		for (i=1; i<data->lseq; i++) fprintf(ofp, "\t%6i", i+1);
	}
	for (i=1; i<data->lseq; i++) {
		if (fl) {
			fprintf(ofp,"\n%5i:", i);
			for (j=1; j<i; j++) fprintf(ofp,"\t");
		}
        for (j=i+1; j<=data->lseq; j++) {
                   t = pij[i][j-i];
                   if (t > 0) {
                        if (data->lc == 'C') {
							dij = (double) locs[j]-locs[i];
/*	Comment out previous line and use following for circular model*/
/*                                dij = mini(locs[j]-locs[i], locs[i]+(data->tlseq)-locs[j]);*/
                                cij = (double) 2*(data->rho)*(1-exp((double) -dij/(data->avc)));
                        }
                        else cij = (double) (data->rho)*(locs[j]-locs[i])/(data->tlseq);
                        if (cij > data->rmax) lke = (double) lkmat[t][data->rcat];
                        else {
                           d = (double) cij*(data->rcat-1)/data->rmax;
                           k = (int) d+1;
                           if (k<data->rcat) lke = (double) lkmat[t][k]+(d-k+1)*(lkmat[t][k+1]-lkmat[t][k]);
                           else lke = lkmat[t][data->rcat];
                        } /*printf(": extrapolated lk = %.3f", lke);*/
			for (k=1, lkmax=-100000000; k<=data->rcat;k++)
				if (lkmat[t][k]>lkmax) {lkmax = lkmat[t][k]; rmp=(double) (k-1)*(data->rmax)/(data->rcat-1);}
			if (fl) {
				if (rmp > cij) fprintf(ofp,"\t%6.2f",(double) lkmax - lke);
				else fprintf(ofp,"\t%6.2f", (double) lke-lkmax);
			}
			clk += lkmax-lke;
                   }
		   else if (fl) fprintf(ofp, "\t NA");
                }
	}

	if (fl) {
		fprintf(ofp, "\n\nTotal deviation = %.3f\n\n", clk);
		fclose(ofp);
	}
	data->fit = (double) clk;
}


void rmin(data, pset, pij, locs, print_flag)
struct data_sum *data;
struct site_type **pset;
int **pij, print_flag;
double *locs;
{
	int i, j, fl, k, rmin=0, **min_mat;
	char fname[MAXNAME+1];
	FILE *ofp;

	printf("\n\nCalculating Rmin (Hudson and Kaplan 1985)\n\n");
	for (i=1;i<=data->ptt;i++) if (pset[i]->nt) { /*Calculate Rmin for each pair type*/
		pset[i]->rm = rec_event(pset[i], data->hd);
	}

	min_mat = imatrix(1,data->lseq, 1, data->lseq);
	for (i=1;i<=data->lseq;i++) for (j=1;j<-data->lseq;j++) min_mat[i][j]=0;
	for (j=1;j<(print_flag?data->lseq:2);j++) {/*Fill in matrix of Rmin values*/
		min_mat[j][j+1]=pij[j][1] ? pset[pij[j][1]]->rm : 0;
		for (k=j+2;k<=data->lseq;k++) {
			for (i=j+1,fl=0;i<k;i++) if (pij[i][k-i] && min_mat[j][i]+pset[pij[i][k-i]]->rm > fl) fl=min_mat[j][i]+pset[pij[i][k-i]]->rm;
			min_mat[j][k]=fl;
		}
	}
	if (print_flag) {
		strcpy(fname, data->prefix);
		ofp = fopen(strcat(fname,"rmin.txt"), "w");
		fprintf(ofp,"\n\nHudson and Kaplan Rmin matrix\n\n");
		for (i=1;i<data->lseq;i++) {
			fprintf(ofp,"\n%4i: ",i);
			for (j=1;j<i;j++) fprintf(ofp,"\t%3f",(double) min_mat[j][i]/(locs[i]-locs[j]));
			fprintf(ofp,"\t"); j++;
			for (;j<=data->lseq;j++) fprintf(ofp,"\t%3i",min_mat[i][j]);
		}
		fclose(ofp);
	}
	data->rmin = (int) min_mat[1][data->lseq];
	free_imatrix(min_mat,1,data->lseq, 1, data->lseq);
}


/*Detects wether all 4 gametes present for haploid or diploid data*/
int rec_event(ptype, hd)
struct site_type *ptype;
int hd;
{
	int i, fl=1;

	if (hd==1) { /*Haploid - just look for all four haplotypes*/
		for (i=0;i<4;i++) if (!(ptype->pt[i])) {fl=0; break;}
	}
	else {/*Diploid - have to look at sets*/
		if (!(ptype->pt[5]+ptype->pt[7]+ptype->pt[13])) fl=0;
		else if (!(ptype->pt[9]+ptype->pt[11]+ptype->pt[13])) fl=0;
		else if (!(ptype->pt[6]+ptype->pt[7]+ptype->pt[14])) fl=0;
		else if (!(ptype->pt[10]+ptype->pt[11]+ptype->pt[14])) fl=0;
		else ; /*Has all four haplotypes*/
	}
	return fl;
}

/*Estimate 4Ner using Wakeley's (1997) moment estimator*/
void wakeley_est(data, seqs, locs)
struct data_sum *data;
int **seqs;
double *locs;
{
	int pos, s1, s2, itrn=0;
	double pwd, cons[3], x[3], y[3];
	printf("\n\nEstimatng 4Ner by Wakeley (1997): ");
	data->avpwd=0; data->varpwd=0;
	if (data->hd==1) {/*Is haploid*/
		for (s1=1;s1<data->nseq;s1++) 
			for (s2=s1+1;s2<=data->nseq;s2++) {
				for (pos=1,pwd=0.0;pos<=data->lseq;pos++)
					if (seqs[s1][pos]>1 && seqs[s2][pos]>1 && seqs[s1][pos]!=seqs[s2][pos]) pwd++;
				data->avpwd += (double) pwd;
				data->varpwd+= (double) pwd*pwd;
			}
		data->avpwd /= (double) data->nseq*(data->nseq-1)/2;
		data->varpwd /= (double) data->nseq*(data->nseq-1)/2;
		data->varpwd -= (double) (data->avpwd)*(data->avpwd);
	}
	else {/*Is diploid*/
		for (s1=1;s1<=data->nseq;s1++) {
			for (pos=1,pwd=0.0;pos<=data->lseq;pos++)
				if (seqs[s1][pos]==4) pwd++;
				data->avpwd += (double) pwd;
				data->varpwd+= (double) pwd*pwd;
			}
		data->avpwd /= (double) data->nseq;
		data->varpwd /= (double) data->nseq;
		data->varpwd -= (double) (data->avpwd)*(data->avpwd);
	}

	printf(" avPWD = %.3f, varPWD = %.3f :",data->avpwd, data->varpwd);
	x[0]=100;x[2]=0.01; /*Starting guess for limits to rho*/
	cons[2] = (double) data->avpwd;
	cons[1] = (double) data->nseq*data->hd;

	while ((y[0]=C_equation(x[0],cons)) > data->varpwd) {x[0]*=2;}
	while ((y[2]=C_equation(x[2],cons)) < data->varpwd) {x[2]*=0.5;}
	while (fabs(x[0]-x[2]) > 1e-4){/*Bisect routine to find rho*/
		x[1]=(x[0]+x[2])/2;
		y[1]=C_equation(x[1],cons);
		if (y[1] < data->varpwd) {x[0]=x[1]; y[0]=y[1];}
		else {x[2]=x[1]; y[2]=y[1];}
 /* If too many iterations give up */
		itrn++;
		if(itrn > 100000){ x[1] = -1; break; }
	}

	if (x[1]>0.0) printf(": 4Ner = %.3f\n\n",x[1]);
	else printf(": beyond range of estimation routine (>%.3f)\n\n",x[1]);
	data->rwak = x[1];
}


double C_equation(C, cons)
double C, cons[3];
{
  double Ia, Ib, gpi, f, S, s97 = (double) sqrt(97.0);

  /* Terms as given in Wakeley(1997) appendix */
  Ia = (double) log(C*C/18.0+13.0/18.0*C+1.0);
  Ib = (double) log((13.0-s97+2.0*C)*(13.0+s97)/((13.0+s97+2.0*C)*(13.0-s97)));

  gpi = (double) (cons[1]-2)/(cons[1]*(cons[1]-1)*C*C)*(-2*C*(cons[1]+1)-(cons[1]-7-C*(cons[1]+1))*Ia+\
							(49*cons[1]-55+C*(15*cons[1]-1))*Ib/s97);

  f = 2/(cons[1]*(cons[1]-1)*C*C)*(-2*C-(2*cons[1]*(cons[1]+1)-7-C)*Ia +(2*cons[1]*(cons[1]+1)*(13+2*C)-55-C)*Ib/s97);
  S = 2*(cons[1]-2)/(3*cons[1]-3)*cons[2]+gpi*(cons[2]*cons[2]-(cons[1]+1)/(3*cons[1]-3)*cons[2])/(f+1);

/* printf("\n %2.2f %2.2f %2.2f %2.2f  ", C, gpie, f, S);  */

 return S;

}


void print_help(int argc, char* argv[]) 
{
	int i;
	char *in_str;
	for(i = 0; i < argc; i++){
		if(*argv[i] == '-'){ 
			in_str = argv[i];
			if((strcmp(in_str, "-h") == 0) || (strcmp(in_str, "-?") == 0) || (strcmp(in_str, "-help") == 0)) {

				printf("\npairwise\n");
				printf("Estimate Constant Population Recombiantion Rate.\n\n");
				
				printf("Required Options :\n");
				printf("-seq <file>           Input sequence file.\n");
				printf("-loc <float>          Input locs file.\n");
				printf("\n\n");	
				printf("Additional Options :\n");
				printf("-lk <file>            Input existing likelihood file.\n");
				printf("-exact                Likelihood file is exact for this dataset.\n");
				/*
				printf("-window               Sliding windows analysis.\n");
				printf("-rmin                 Calculte minimum recombination events.\n");
				printf("-moment               Estimate 4Ner by moment method.\n");
				printf("-test                 Test hypothesis that 4Ner=0.\n");
				printf("-simulate             Simulate data to test model fit and estimate sampling distirbution.\n");
*/
				printf("-concise              Set to concise mode.\n");
				printf("-prefix <string>      Prefix of output files\n");
				printf("\n\n");	
				exit (0);
			}
		}
	}
}

