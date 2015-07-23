#include "ldhat.h"
#include "tools.h"
#include "seqtools.h"

long *idum;
int sizeofpset=100;
double *lnfac_array;


void print_help(int argc, char* argv[]);


int main (int argc, char *argv[])
{
	int i, j, **seqs, **nall, ord=1, ns, **pij, lkf=0, npt=0, pnew=0, anc=0, npmax;
	int tcat=1, rcat=0, verb=1, miss=0;
	char fname[MAXNAME+1], **seqnames;
	long seed=-setseed();
	extern int sizeofpset;
	double *locs;
	double **lkmat, *lkres;
	struct site_type **pset;
	struct data_sum *data;
	FILE *ifp=NULL, *ifp2=NULL, *ifp3=NULL, *tfp;
	char* in_str;

	print_help(argc, argv);

	data = malloc((size_t) sizeof(struct data_sum));
	data->bpen = -1;
	data->n_update = -1;
	data->r_update = -1;
	data->exact = 0;
	
	strcpy(data->prefix, "");
	
	for(i = 0; i < argc; i++)
	{
		if(*argv[i] == '-')
		{ 
			in_str = argv[i];
			if(strcmp(in_str, "-seed") == 0) {seed = atol(argv[i + 1]);						/* User defined random seed */
												if (seed > 0) seed = -seed;}
			if(strcmp(in_str, "-lk") == 0) {ifp3 = fopen(argv[i+1], "r"); lkf=1;}; 			/* Composite Likelihood filename */
			if(strcmp(in_str, "-seq") == 0) ifp = fopen(argv[i+1], "r");					/* Sequence file */
			if(strcmp(in_str, "-loc") == 0) ifp2 = fopen(argv[i+1], "r");					/* Loci file */
			if(strcmp(in_str, "-exact") == 0) data->exact=1;								/* Likelihood file is specifically for this dataset flag */
			if(strcmp(in_str, "-concise") == 0) verb=0;										/* Concise screen output */
			if(strcmp(in_str, "-bpen") == 0) data->bpen = atof(argv[i + 1]);				/* Block penalty */
			if(strcmp(in_str, "-its") == 0) data->n_update = atoi(argv[i + 1]);				/* Number of iterations */
			if(strcmp(in_str, "-samp") == 0) data->r_update = atoi(argv[i + 1]);			/* Sample every? */
			if(strcmp(in_str, "-prefix") == 0) strcpy(data->prefix, argv[i+1]);					/* Prefix for output filenames */
		}
	}

	/*	seed = 122343222; */
	idum = &seed;
	if (ifp == NULL)
	{
		printf("\n\nCould not find sequence file in command line\n");
		printf("\nInput filename for seqs\n\n");
		scanf("%s", fname);
		ifp = fopen(fname, "r");
	}
	if (ifp == NULL) nrerror("Error in opening sequence file");

/*Read #seqs, #SNPs and haploid/diploid*/
	fscanf(ifp,"%i %i %i\n", &data->nseq, &data->lseq, &data->hd);

	if ((data->nseq < 2) || (data->lseq < 2)) {printf("\n\nInsufficient data for analysis (n > 1, L > 1) \n\n"); exit(1);}
	if (data->nseq > SEQ_MAX) {printf("\n\nMore than max no. sequences: Using first %i for analysis\n\n", SEQ_MAX); data->nseq=SEQ_MAX;}
	if (data->hd==1) printf("\nAnalysing %i haploid sequences of length %i seg sites\n", data->nseq, data->lseq);
	else printf("\nAnalysing %i genotypes of length %i seg sites\n", data->nseq, data->lseq);

/*Read sequences */	
	seqs = imatrix(1, data->nseq, 1, data->lseq);
	seqnames = cmatrix(1, data->nseq+1, 1, MAXNAME+1);
	if (read_fasta(seqs, ifp, data->nseq, data->lseq, seqnames)) printf("\nSequences read succesfully\n");
	fclose(ifp);

	nall = imatrix(1, data->lseq, 1, 6);
	allele_count(seqs, data->nseq, data->lseq, nall, 0, data->hd, data->prefix);

	/*Store lnfac values in array for speed of computation*/
	lnfac_array = (double *) malloc((size_t) ((int) data->nseq*data->hd+2)*sizeof(double));
	lnfac_array[0]=lnfac_array[1]=0;
	for (j=2;j<=((int) data->nseq*data->hd);j++) lnfac_array[j]=lnfac_array[j-1]+log(j);

	if (ifp2 == NULL)
	{
		printf("\nInput name of file containing location of seg sites\n\n");
		scanf("%s", fname);
		ifp2 = fopen(fname, "r");
	}
	if (ifp2 == NULL) nrerror("Cannot open loc file");
	fscanf(ifp2, "%i %lf %c", &ns, &data->tlseq, &data->lc);
	printf("%f\n", data->tlseq);
	if (ns != data->lseq) nrerror("Lseq and Locs disagree");
	if ((data->lc != 'C')&&(data->lc != 'L')) nrerror("Must input linear(L)/circular(C)");
	if (data->lc == 'C') {
	  data->avc=0;
	  while (data->avc <= 0) {
	    printf("\n\nInput average tract length for conversion model: ");scanf("%lf", &(data->avc));
	  }}

	locs = dvector(1, data->lseq);

	for (i=1; i<=data->lseq; i++) {
		fscanf(ifp2, "%lf", &locs[i]); 
		/*printf("\n%i %f", i, locs[i]);*/
		if ((locs[i]<=0)||(locs[i]-locs[0]>data->tlseq)) 
		{
		   printf("%f\n", data->tlseq);
		   nrerror("Error in Loc file\n\n");
		}
		if (i>1 && locs[i]<=locs[i-1])
		{
			printf("Error at SNP %i and %i: %f %f\n", i, i-1, locs[i], locs[i-1]);
			nrerror("Error in locs file: SNPs must be monotonically increasing");
		}
	}
	/*printf("\nLocation of seg sites\n\n");
	for (i=1; i<=data->lseq; i++) printf("%3i   %9.3f\n", i, locs[i]); */
     	fclose(ifp2);

	if (ifp3 == NULL) 

	{
		printf("\n\nInput name of lookup likelihood file: ");
		scanf("%s", fname);

		ifp3 = fopen(fname, "r");
		printf("\n\nIs likelihood file exact (1/0) ?");

		scanf("%i", &data->exact);

		lkf=1;
	}
	if ((lkf) && (ifp3==NULL)) nrerror("Cannot open likelihood file");
	
	/*if (argc > 4) {printf("\n\nSet to concise output\n\n"); verb=0;}*/
	
	data->w = mini(data->lseq,MAXW);
	pij = imatrix(1,data->lseq,1,data->w);
	for (i=1;i<=data->lseq;i++) for (j=1;j<=data->w;j++) pij[i][j]=0;
	pset = init_pset(pset, lkf, ifp3, &npt, data);

	/*Check that all haplotypes are present in likelihood file*/
     if (!data->exact) check_exhaustive(pset,npt,(int) (data->nseq)*((int) data->hd));

	printf("\nCalculating distribution of pair types\n");
	pset = pair_spectrum(seqs, data, nall, pset, &npt, &pnew, &miss, anc, pij);
	printf("\nCompleted classification of pair types \n");
	if (lkf) printf("\nOld = %i: New = %i: Missing = %i\n", npt,pnew,miss);
	data->ptt = (int) npt+pnew+miss;
	if (data->hd==1 && pnew) nrerror("Cannot have haploid data and new types!!");

	if (verb)
	{
		strcpy(fname, data->prefix);
		tfp = fopen(strcat(fname, "type_table.txt"), "w");
		if (!tfp) nrerror("Cannot open type file");
		type_print(pij, data->lseq,data->w, tfp);
		fclose(tfp);
	}

	if (verb)
		print_pairs(stdout, pset, npt+pnew+miss, data->hd, data->nseq);
	i = (int) (data->hd==1? ((int) data->nseq/2): data->nseq);
	npmax = (int) 1+i+i*(i-1)*(i+4)/6+(i-1)*(i+2)/2;

	printf("\nMax number of haplotypes for n=%i is %i\n",data->nseq*((int) data->hd==1?1:2), npmax);

	if (lkf)
		read_pars(ifp3, &tcat, &data->th, &data->rcat, &data->rmax);
	else exit(1);

	lkmat = dmatrix(1,npt+pnew+miss,1,data->rcat);
	if (lkf)
		read_lk(ifp3, lkmat, npt, tcat, data->rcat);

	if (miss && data->hd==1)
	{
		printf("\nCalculating LKs for missing data (haplotypes: %i)\n",miss);
		for (i=1;i<=miss;i++)
		{
			lk_miss(pset[npt+i],lkmat[npt+i],lkmat,data);
		/*	printf("."); if (!(i%50)) printf(":%i\n",i); */
		}
	}
	else if (data->hd == 2 && !(data->exact))
	{
		printf("\nResolving diploid data: %i\n",pnew+miss);
		lkres = dvector(1,data->rcat);
		for (i=1;i<=pnew+miss;i++)
		{
			lk_resolve(lkres,pset[npt+i],lkmat[npt+i],lkmat,data);
			/* printf("."); if (!(i%50)) printf(": %i\n",i); */
		}
		free_dvector(lkres,1,data->rcat);
	}

	if (verb && !(data->exact))
		print_lks(pset, data, npt+pnew+miss, lkmat);
	data->rmap = (double *) malloc((size_t) (data->lseq+1)*sizeof(double));

	lk_block(pset, pij, data, lkmat, locs, verb);

	free_dvector(locs,1,data->lseq);
	free_imatrix(pij,1,data->lseq,1,data->w);
	free_imatrix(seqs,1,data->nseq,1,data->lseq);
	free_imatrix(nall,1,data->lseq,1,5);
	for (i=1;i<sizeofpset;i++) free(pset[i]);
	free(pset);
	free(data);

/*	system("PAUSE");*/
}





/*Check that likelihood file is exhaustive for n*/

void check_exhaustive(pset,npt,nsamp)
     struct site_type **pset;
     int npt, nsamp;
{
  int p1, p2, i, ei;

  printf("\nChecking likelihood file is exhaustive:...");

  p1 = (int) nsamp/2;
  if (npt != (int) 1+p1+p1*(p1-1)*(p1+4)/6+(p1-1)*(p1+2)/2) {
    printf("\n\n!!npt = %i: E[npt] = %i\n\n",npt, 1+p1+p1*(p1-1)*(p1+4)/6+(p1-1)*(p1+2)/2);
    nrerror("npt and total do not agree: not exhaustive. Do you need the -exact flag?");
  }

  for (i=1;i<=npt;i++) {
    p1 = pset[i]->pt[1]+pset[i]->pt[3];
    p2 = pset[i]->pt[2]+pset[i]->pt[3];
    ei = (int) p1*(p1-1)*(p1+4)/6+(p2-1)*(p2+2)/2+pset[i]->pt[2]+1;
    if (i!=ei) {printf("\nError in pair-types: p1(%i), p2(%i), n11(%i): i(%i) ei(%i)\n\n",p1,p2,pset[i]->pt[3],i,ei); exit(1);}
  }
  printf("OK\n");
}



void print_lks(pset,data,npt,lkmat) 
int npt;
double **lkmat;
struct site_type **pset;

struct data_sum *data;
{
	int p, i, nstate, ct=0;
	char filename[MAXNAME];
	FILE *ofp;

	strcpy(filename, data->prefix);
	ofp = fopen(strcat(filename, "new_lk.txt"), "w");
	if (data->hd == 2) nstate=16;
	else nstate=9;

	for (p=1;p<=npt;p++) if (pset[p]->nt>0) ct++;

	fprintf(ofp, "\n%i %i\n1 %.5f\n%i %f\n\n",(int) data->nseq*data->hd,ct,data->th,data->rcat,data->rmax);
	for (p=1; p<=npt; p++) if (pset[p]->nt>0) {
		fprintf(ofp,"\n%4i # ", p);
		for (i=0; i<nstate; i++) fprintf(ofp,"%3i ", pset[p]->pt[i]);
		fprintf(ofp," :  ");
		for (i=1; i<=data->rcat; i++) fprintf(ofp,"%7.2f ", lkmat[p][i]);
	}
	fclose(ofp);
}



/*
Routine to calculate the pairwise likelihood for any given rho
NB in previous version rmap was per interval, for the block routine it is cumulative

  lij is the matrix containing the likelihood for each pair of sites in current state
  pij is the matrix of pair types
  rl and ru and the lower and upper bounds of the regions being updated
  upate determines the update region - block and links for type 0, block, links and across for type 1
  dlk is a pointer to the difference in likelihood following update

*/

int lk_calc_site(lij, rl, ru, pij, data, dlk, lkmat, update) 
int **pij, rl, ru, update;
double *dlk, **lkmat, **lij;
struct data_sum *data;
{

	int i, j, k, t, fl=1, ct=0;
	double cij, d;

	/*To check MCMC routine - just return zero log likelihood
	*dlk = 0.0;
	return fl;
*/

	if (!update) {/*Only update blocks and associations to blocks*/
	  for (i=1, (*dlk)=0.0; i<=ru; i++){
		  for (j=(i<rl?rl:i+1); j<=(i<rl?mini(ru,i+data->w):mini(data->lseq,i+data->w)); j++) {/*Upper and lower bound depends on whether i<rl*/
			  t = pij[i][j-i];/*Type of pair*/
			  if (t > 0) {
				  cij = data->rmap[j]-data->rmap[i];/*4Ner between sites*/
				  if (cij > data->rmax) {
					  lij[i][j-i+data->w] = (double) lkmat[t][data->rcat];
					  fl=2;
				  }
				  else {
					  d = (double) cij*(data->rcat-1)/(data->rmax);
					  k = (int) d+1;
					  if (k<data->rcat) lij[i][j-i+data->w] = (double) lkmat[t][k]+(d-k+1)*(lkmat[t][k+1]-lkmat[t][k]);
					  else lij[i][j-i+data->w] = (double) lkmat[t][data->rcat];
				  }
				  if (lij[i][j-i+data->w] < 0.0) (*dlk) += (double) lij[i][j-i+data->w]-lij[i][j-i]; /*Add difference in LK to lkrun*/
				  else {
					  printf(" Lk >= 0!!  A t=%i cij=%f",t,cij);
					  exit(0); 
					  return 0;
				  }
			  }
		  }
	  } 
	  return fl;
	}
	else {/*Update block, associations to block and associations across blocks*/
		for (i=1, (*dlk)=0.0; i<=ru; i++) {
			for (j=(i<rl?rl:i+1); j<=mini(data->lseq, i+data->w); j++) {
				t = pij[i][j-i];
				if (t > 0) {
					cij = data->rmap[j]-data->rmap[i];/*4Ner between sites*/
					if (cij > data->rmax) {
						lij[i][j-i+data->w] = (double) lkmat[t][data->rcat];
						fl=2;
					}
					else {
						d = (double) cij*(data->rcat-1)/(data->rmax);
						k = (int) d+1;
						if (k<data->rcat) lij[i][j-i+data->w] = (double) lkmat[t][k]+(d-k+1)*(lkmat[t][k+1]-lkmat[t][k]);
						else lij[i][j-i+data->w]= (double) lkmat[t][data->rcat];
					}
					if (lij[i][j-i+data->w] < 0.0) { (*dlk) += (double) lij[i][j-i+data->w]-lij[i][j-i];}
					else {printf(" Lk >= 0!! B t=%i cij=%f",t,cij); exit(0); return 0;}
				}
			}
		}
		return fl;
	}
}


void lk_block(pset,pij,data,lkmat,loc, verb)
int **pij;
double *loc;
double **lkmat;
struct site_type **pset;
struct data_sum *data;
int verb;
{
	int i, j, run, nblock, nacc[4][2];
	double rt, bptrue;
	double pr[4];
	struct block *new_block, **blockp;
	char filename[MAXNAME];
	double **lij, lk[2];
	FILE *rfp, *bfp;

	printf("\nMCMC routine to define block structure\n");

	/*Set cumulative transition probabilities for events in proposal*/
	pr[0]=0.4;/*Change value within block*/
	pr[1]=0.6;/*Extend block by one*/
	pr[2]=0.8;/*Split blocks*/
	pr[3]=1.0;/*Merge blocks*/

	strcpy(filename, data->prefix);
	rfp = fopen(strcat(filename, "rates.txt"),"w");

	if (verb)
	{
		strcpy(filename, data->prefix);
		bfp = fopen(strcat(filename, "bounds.txt"),"w");
	}

	lij = dmatrix(1,data->lseq,1,2*(data->w));
	for (i=1;i<=data->lseq;i++) for (j=1;j<=2*(data->w);j++) lij[i][j]=0.0;

	rt = 100;
	rt/=(double) data->tlseq;
	if (rt<RMIN) rt=RMIN+0.000001;
	else if (rt>RMAX) rt=RMAX-0.000001;

	if (data->bpen < 0)
	{
		printf("\n\nInput block penalty:");
		scanf("%lf", &data->bpen);
	}

	if (data->n_update <= 0)
	{
		printf("\n\nHow many updates for MCMC? ");
		scanf("%i", &data->n_update);
	}

	if (data->n_update<BURNIN)
	{
		printf("\n\nWARNING: # updates < min allowed BURNIN: resetting updates to %i\n\n", (int) BURNIN*2); data->n_update=(int) BURNIN*2;
	}

	if (data->r_update < 0)
	{
		printf("\n\nNumber of updates between samples:");
		scanf("%i", &data->r_update);
	}

	bptrue=data->bpen; /*During burn-in, linear increase in penalty to allow algorithm to get to good starting point*/

	fprintf(rfp,"%i %i\n",(int) ((double) data->n_update/data->r_update), data->lseq);
	if (verb)
		fprintf(bfp,"%i %i\n",(int) ((double) data->n_update/data->r_update), data->lseq);

	/*Initialise block array and recombination map*/

	data->rmap[0]=data->rmap[1]=0.0;
	for (i=0;i<4;i++) nacc[i][0]=nacc[i][1]=0;
	nblock=data->lseq-1;
	blockp = (struct block **) malloc((size_t) nblock*sizeof(struct block *));
	for (i=0;i<nblock;i++)
	{
		new_block = (struct block *) malloc((size_t) sizeof(struct block));
		new_block->num=i;
		new_block->size=1;
		new_block->pos=i+1;
		new_block->rate=rt;

		if (i) new_block->bpl=blockp[i-1];
		else new_block->bpl=NULL;
		blockp[i]=new_block;
	}
	blockp[nblock-1]->bpr=NULL;
	for (i=nblock-2;i>=0;i--) blockp[i]->bpr=blockp[i+1];

	/*Convert blocks to recombination map*/
	block2map(blockp[0],data->rmap,loc,data->lseq);

	/*Calculate the likelihood of the data*/
	lk_calc_site(lij,1,data->lseq,pij,data,&lk[1],lkmat,1);
	for (i=1;i<=data->lseq;i++) for (j=i+1;j<=mini(data->lseq,i+data->w);j++) lij[i][j-i]=lij[i][j-i+data->w];
	lk[0]=lk[1];
	printf("\nInitial likelihood = %.3f\n",lk[0]);

	for (run=1;run<=data->n_update;run++)
	{
		if (run<BURNIN)
		{
			data->bpen=(double) bptrue*run/BURNIN; /*Allows algorithm to get to good starting point before penalty kicks in*/
			if (run<BURNIN/2) pr[0]=1.0;
			else pr[0]=0.4;
		}
		else
		{
			data->bpen=(double) bptrue;
			pr[0]=0.4;
		}

		blockp = update_blocks(pr,blockp,&nblock,loc,lk,nacc,data,lkmat,pij,lij);
		/*	  check_lk(lij,lk[0],data->lseq);
		check_blocks(blockp[0],data->rmap,loc);*/
		if (DEBUG)
		{
		  printf("\nRun %i\n", run);
		  print_block(blockp[0]);
		  print_rmap(data->rmap,data->lseq,loc);
		  check_blocks(blockp[0],data->rmap,loc);
		}
		if (!(run%(data->r_update)))
		{
			fprintf(rfp,"%8.6f\t",data->rmap[data->lseq]);
			if (verb)
			{
				fprintf(bfp,"%i\t",nblock);
				print_bounds(blockp[0],bfp);
			}
			print_rates(blockp[0],rfp);
			for (i=1,lk[0]=0.0;i<data->lseq;i++) for (j=i+1;j<=mini(data->lseq,i+data->w);j++) lk[0]+=lij[i][j-i];
			if (verb)
				printf("\nRun %8i: LK = %.3f : # blocks = %5i : Map length = %.3f",run,lk[0],nblock,data->rmap[data->lseq]);
		}
	}
	printf("\nAfter %i runs, likelihood = %f\n",data->n_update, lk[0]);
	/*print_block(blockp[0]);
	print_rmap(data->rmap,data->lseq,loc);*/
	printf("Final block length = %f\n", data->rmap[data->lseq]);

	printf("\n\nAcceptance rates: \nRate change = %.4f\nExtend = %.4f\nSplit = %.4f\nMerge = %.4f\n\n",\
	  (double) nacc[0][1]/nacc[0][0], (double) nacc[1][1]/nacc[1][0], (double) nacc[2][1]/nacc[2][0], (double) nacc[3][1]/nacc[3][0]);

	free_dmatrix(lij,1,data->lseq,1,2*(data->w)+2);
	fclose(rfp);
	if (verb)
	  fclose(bfp);
}


struct block ** update_blocks(pr,blockp,nblock,loc,lk,nacc,data,lkmat,pij,lij)
     struct block **blockp;
     int *nblock, nacc[4][2], **pij;
     double *loc;
     double pr[4];
	 double **lkmat, **lij, lk[2];
	 struct data_sum *data;
{

  int i, j, split, fl;
  double or, rt;
  double r, rl, ru, *rcheck, ll, lr, u;

  double hr, jcb;
  struct block *new_block, *obp;

/*  rcheck = malloc((size_t) (data->lseq+1)*sizeof(double)); */

  /*Choose block and event for updating*/
  i = (int) ((double) (*nblock)*ran2());
  if (DEBUG) printf("\nUpdating block %i\n", i+1);

  /*Changing value within block*/
  if ((r=ran2())<pr[0]) {
    or = blockp[i]->rate;
	nacc[0][0]++;
    blockp[i]->rate *= (double) exp(ran2()*2-1.0);
    if (blockp[i]->rate >= RMAX)  blockp[i]->rate = (double) RMAX-0.001;
	else if (blockp[i]->rate <= RMIN) blockp[i]->rate = (double) RMIN+0.000001;
    block2map(blockp[i],data->rmap,loc,data->lseq);

	/*Update lk following change*/
    lk_calc_site(lij,blockp[i]->pos,blockp[i]->pos+blockp[i]->size,pij,data,&lk[1],lkmat,1);/*Note lk[1] stores change in LK*/
    if (DEBUG) {
		printf("\n\nChanging value within block %i\n", i+1);
		printf("\n\nUpdated blocks and map\n");
		print_rmap(data->rmap,data->lseq,loc);
		printf("Rate %f ->%f | lk %f->%f ",or,blockp[i]->rate,lk[0],lk[1]);
	}

	hr = (double) blockp[i]->rate/or; /*HRatio for move*/

/*Rejection routine - includes difference in log lk for data, prior for change in rate & HR for scaling move*/
    if (log(ran2()) > lk[1]+(or-blockp[i]->rate)+log(hr)) {
      blockp[i]->rate = or; 
      block2map(blockp[i],data->rmap,loc,data->lseq);
      if (DEBUG) printf(": REJECT\n");
    }
    else {/*Accept*/
		lk[0]+=lk[1]; 
		nacc[0][1]++; 
		update_lij(lij,blockp[i]->pos,blockp[i]->pos+blockp[i]->size,data,1);
		if (DEBUG) printf(":ACCEPT\n");
	}
  }

  /*Extend block by one*/
  else if (r<pr[1]) {
    if (DEBUG) printf("\nExtending block - ");
	nacc[1][0]++;
    if (ran2()<0.5) {
		if (blockp[i]->bpr) {/*Extend block right - only if size of Rblock > 1*/
			if ((blockp[i]->bpr)->size>1) {
				if (DEBUG) printf("RIGHT\n");
				blockp[i]->size++;
				(blockp[i]->bpr)->size--;
				(blockp[i]->bpr)->pos++;

/*Change rate of block to keep average over interval constant - allows lk to be updated efficiently*/
/*NB no need to check for constraints as average must lie within bounds*/

/*Also, have to add randomness to move to keep reversibility*/

				if ((u=ran2())<0.5) {/*Update rate for block being extended*/
					or = blockp[i]->rate;
					blockp[i]->rate = (double) (data->rmap[(blockp[i]->bpr)->pos]-data->rmap[blockp[i]->pos])/(loc[(blockp[i]->bpr)->pos]-loc[blockp[i]->pos]);
					block2map(blockp[i],data->rmap,loc,(blockp[i]->bpr)->pos);/*Update recombination map - only up to next block*/
/*					for (j=0;j<=data->lseq;j++) rcheck[i]=rmap[i];
					printf("\nChecking maps agree - R\n");
					block2map(blockp[i],rmap,loc,data->lseq);
					for (j=1;j<=data->lseq;j++) if (fabs(rcheck[i]-rmap[i])>0.01) {printf("\n\nError in rmap\n\n"); exit(1);}
*/
					if (blockp[i]->rate<RMIN || blockp[i]->rate>RMAX) fl=1;
					else {lk_calc_site(lij,blockp[i]->pos,blockp[i]->pos+blockp[i]->size,pij,data,&lk[1],lkmat,0); fl=0;}

				}
				else {/*Update rate for block being shrunk*/
					or = (blockp[i]->bpr)->rate;
					(blockp[i]->bpr)->rate = (double) (data->rmap[(blockp[i]->bpr)->pos+(blockp[i]->bpr)->size]-data->rmap[(blockp[i]->bpr)->pos-1]-blockp[i]->rate*(loc[blockp[i]->bpr->pos]-loc[blockp[i]->bpr->pos-1]))/(loc[(blockp[i]->bpr)->pos+(blockp[i]->bpr)->size]-loc[(blockp[i]->bpr)->pos]);
					block2map(blockp[i],data->rmap,loc,(blockp[i]->bpr)->pos+(blockp[i]->bpr)->size);/*Update recombination map on RHS - only up to next block*/
					if ((blockp[i]->bpr)->rate<RMIN || (blockp[i]->bpr)->rate>RMAX) fl=1;
					else {lk_calc_site(lij,(blockp[i]->bpr)->pos-1,(blockp[i]->bpr)->pos+(blockp[i]->bpr)->size,pij,data,&lk[1],lkmat,0); fl=0;}
				}

				if (DEBUG) {
					printf("\n\nUpdated blocks and map\n");
					print_block(blockp[0]);
					print_rmap(data->rmap,data->lseq,loc);
					printf("block %i R: lk %f->%f ",i+1,lk[0], lk[1]);
				}

/*Rejection includes dLK and prior for rates - HR for move =1*/
				if (fl || log(ran2())>lk[1]+(or-(u<0.5?blockp[i]->rate:(blockp[i]->bpr)->rate))) {/*REJECT*/
					if (DEBUG) printf(": REJECT\n");
					blockp[i]->size--;
					(blockp[i]->bpr)->size++;
					(blockp[i]->bpr)->pos--;

					if (u<0.5) {
						blockp[i]->rate=or;
						block2map(blockp[i],data->rmap,loc,(blockp[i]->bpr)->pos+1);
					}
					else {
						(blockp[i]->bpr)->rate=or;
						block2map(blockp[i],data->rmap,loc,(blockp[i]->bpr)->pos+(blockp[i]->bpr)->size);
					}
				}
				else {/*ACCEPT*/
					if (DEBUG) printf(": ACCEPT\n");
					nacc[1][1]++;
					lk[0]+=lk[1];
					if (u<0.5) update_lij(lij,blockp[i]->pos,blockp[i]->pos+blockp[i]->size,data,0);
					else update_lij(lij,(blockp[i]->bpr)->pos-1,(blockp[i]->bpr)->pos+(blockp[i]->bpr)->size,data,0);
				}
			}
		}
	}

/*Else extend block left*/
	else if (blockp[i]->bpl) {/*Ignore if first block*/
		if ((blockp[i]->bpl)->size>1) {/*Only if size of previous block > 1*/
			if (DEBUG) printf("LEFT\n");
			blockp[i]->size++;
			blockp[i]->pos--;
			(blockp[i]->bpl)->size--;
			or = blockp[i]->rate;
/*As above*/
			if ((u=ran2())<0.5) {/*Update rate for block being extended*/
				or = blockp[i]->rate;
				blockp[i]->rate = (double) (data->rmap[blockp[i]->pos+blockp[i]->size]-data->rmap[(blockp[i]->pos)])/(loc[blockp[i]->pos+blockp[i]->size]-loc[blockp[i]->pos]);
/*				for (j=0;j<=data->lseq;j++) rcheck[i]=rmap[i];
				printf("\nChecking maps agree - L\n");
				block2map(blockp[i],rmap,loc,data->lseq);
				for (j=1;j<=data->lseq;j++) if (fabs(rcheck[i]-rmap[i])>0.01) {printf("\n\nError in rmap\n\n"); exit(1);}
*/
				if (blockp[i]->rate<RMIN || blockp[i]->rate>RMAX) fl=1;
				else {
					block2map(blockp[i]->bpl,data->rmap,loc,blockp[i]->pos+blockp[i]->size);
					lk_calc_site(lij,blockp[i]->pos,blockp[i]->pos+blockp[i]->size,pij,data,&lk[1],lkmat,0); 
					fl=0;
				}
			}
			else {/*Update rate for block being shrunk*/
				or = (blockp[i]->bpl)->rate;
				(blockp[i]->bpl)->rate= (double) (data->rmap[blockp[i]->pos+1]-data->rmap[(blockp[i]->bpl)->pos]-blockp[i]->rate*(loc[blockp[i]->pos+1]-loc[blockp[i]->pos]))/(loc[blockp[i]->pos]-loc[(blockp[i]->bpl)->pos]);
				if (blockp[i]->bpl->rate<RMIN || blockp[i]->bpl->rate>RMAX) fl=1;
				else {
					block2map(blockp[i]->bpl, data->rmap, loc, blockp[i]->pos+1);
					lk_calc_site(lij,(blockp[i]->bpl)->pos,blockp[i]->pos+1,pij,data,&lk[1],lkmat,0);
					fl=0;
				}
			}

			if (DEBUG) {
				printf("\n\nUpdated blocks and map\n");
				print_block(blockp[0]);
				print_rmap(data->rmap,data->lseq,loc);
				printf("block %i L: lk %f->%f ",i+1,lk[0],lk[1]);
			}

			if (fl || log(ran2())>lk[1]+(or-(u<0.5?blockp[i]->rate:blockp[i]->bpl->rate))) {/*REJECT*/
				if (DEBUG) printf(": REJECT\n");
				blockp[i]->size--;
				blockp[i]->pos++;
				(blockp[i]->bpl)->size++;
				if (u<0.5) {
					blockp[i]->rate=or;
					block2map(blockp[i]->bpl,data->rmap,loc,blockp[i]->pos+blockp[i]->size);
				}
				else {
					(blockp[i]->bpl)->rate=or;
					block2map(blockp[i]->bpl, data->rmap, loc, blockp[i]->pos+1);
				}
			}
			else {/*ACCEPT*/
				if (DEBUG) printf(": ACCEPT");
				nacc[1][1]++;
				lk[0]+=lk[1];
				if (u<0.5) update_lij(lij,blockp[i]->pos,blockp[i]->pos+blockp[i]->size,data,0);

				else update_lij(lij,(blockp[i]->bpl)->pos,blockp[i]->pos+1,data,0);
			}
		}
	}
  }


  /*Split block - generate new block: NB reversible jump step*/
  else if (r<pr[2]) {
	  nacc[2][0]++;
    if (DEBUG) printf("\nSplitting block");
	if (blockp[i]->size>1) {
    new_block = malloc((size_t) sizeof(struct block));
    (*nblock)++;
    j = (*nblock)-1;
    blockp = realloc(blockp,(*nblock)*sizeof(struct block *));
    blockp[j]=new_block;
    blockp[j]->num=j;
    split = (int) ((double) ran2()*(blockp[i]->size-1))+1;/*Choose size of LHS in terms of SNPs*/
    blockp[j]->pos=blockp[i]->pos+split;/*Locate start of RHS*/
    blockp[j]->size=blockp[i]->size-split;/*Define size of RHS*/
    blockp[i]->size=split;/*Define size of LHS in terms of SNPs*/
    blockp[j]->bpl=blockp[i];/*Define Lpointer of RHS*/
    blockp[j]->bpr=blockp[i]->bpr;/*Define Rpointer of RHS - can be NULL*/
    blockp[i]->bpr=blockp[j];/*Define Rpointer of LHS - points to existing one*/
    if(blockp[j]->bpr) (blockp[j]->bpr)->bpl=blockp[j];/*Update Lpointer of far RHS - only if it exists!*/
	or = blockp[i]->rate;


/*Choose which side of split to change rate of - rate chosen from prior (exp within bounds)*/

	rt = (double) or*(loc[blockp[j]->pos+blockp[j]->size]-loc[blockp[i]->pos]);/*In case split last have to be careful*/
/*Check that maps agree with rates*/
	if (fabs(rt-data->rmap[blockp[j]->pos+blockp[j]->size]+data->rmap[blockp[i]->pos])/rt >0.01) {
		printf("\n\nError in splitting: rt = %f, dmap=%f\n\n",rt,data->rmap[blockp[j]->pos+blockp[j]->size]-data->rmap[blockp[i]->pos]);
/*		exit(1);*/
	}
	lr = (double) loc[blockp[j]->pos+blockp[j]->size]-loc[blockp[j]->pos];
	ll = (double) loc[blockp[j]->pos]-loc[blockp[i]->pos];
	hr = (double) (*nblock-1)/(*nblock)*(blockp[i]->size+blockp[j]->size-1);

    /*Choose LHS rate first*/
	rl = (double) maxd((double) RMIN,(double) (rt-RMAX*lr)/ll);/*define bounds for new rate - approx 0*/
	ru = (double) mind((double) RMAX,(double) (rt-RMIN*lr)/ll);
	r = ran2();  /*Draw [U(0,1)] -> (1-U)/U for ratio of rates p1/p2*/
	blockp[i]->rate = or*(ll+lr)/(ll+r/(1-r)*lr);
	if (blockp[i]->rate > rl && blockp[i]->rate <ru) {
	  blockp[j]->rate = (or*(ll+lr)-blockp[i]->rate*ll)/lr; 
	  jcb = (double) pow(blockp[i]->rate+blockp[j]->rate,2)/or;  /*for geometric proposal*/
/*	  jcb = (double) (ll+lr)/lr*(1-exp(-or*(ll+lr)/ll))/exp(-blockp[i]->rate); /*For exponential prior update*/
/*	  jcb = (double) or*(ll+lr)*(ll+lr)/(ll*lr); /*Jacobian for linear choice*/
	  fl=1;
	}
	else fl=0;  /*fl checks whether new value is compatible - simple rejection if not*/

	if (fl) {
	  block2map(blockp[i],data->rmap,loc,blockp[j]->pos+blockp[j]->size);/*Update map*/
/*	  printf("\nChecking maps agree\n");
	  for (j=0;j<=data->lseq;j++) rcheck[i]=rmap[i];
	  block2map(blockp[i],rmap,loc,data->lseq);
	  for (j=1;j<=data->lseq;j++) if (fabs(rcheck[i]-rmap[i])>0.01) {printf("\n\nError in rmap\n\n"); exit(1);}
 */
	  if (DEBUG) {print_block(blockp[0]); print_rmap(data->rmap,data->lseq,loc);}
	  lk_calc_site(lij,blockp[i]->pos,blockp[j]->pos+blockp[j]->size,pij,data,&lk[1],lkmat,0);
	}

/*Reject on likelihood, rates and numbers of blocks - 
	NB # blocks is sufficient if we assume exponential for length distribution*/
    if (!fl || log(ran2())>lk[1]+(or-blockp[i]->rate-blockp[j]->rate)-(double) data->bpen+log(hr*jcb)) {/*REJECT*/
      if (DEBUG) printf(": REJECT\n");
      blockp[i]->size+=blockp[j]->size;/*Revert to original size*/
      blockp[i]->bpr=blockp[j]->bpr; /*Revert PR pointer*/
      if (blockp[j]->bpr) (blockp[i]->bpr)->bpl=blockp[i];/*Revert LH pointer of far RHS*/
	  blockp[i]->rate=or; /*Revert rate*/
      block2map(blockp[i],data->rmap,loc,blockp[i]->pos+blockp[i]->size);/*Update map*/
      free(blockp[j]);
      (*nblock)--;
      blockp = (struct block **) realloc(blockp,(*nblock)*sizeof(struct block *));
    }
    else {/*ACCEPT*/
      if (DEBUG) printf(": ACCEPT\n");
	  nacc[2][1]++;
	  lk[0]+=lk[1];
	  update_lij(lij,blockp[i]->pos,blockp[j]->pos+blockp[j]->size,data,0);
    }
	}
	else if (DEBUG) printf(" - cannot split block\n");
  }


  /*Join blocks - always merge with RHS: reversible jump step*/
  else {
    nacc[3][0]++;
    if (blockp[i]->bpr!=NULL) {/*Ignore if last block*/
      obp = blockp[i]->bpr;/*keep a pointer to the block being removed*/
      or = blockp[i]->rate;
      blockp[i]->size+=(blockp[i]->bpr)->size;/*Add size of RHS*/
      blockp[i]->bpr=(blockp[i]->bpr)->bpr;/*Update RH pointer - NB can be NULL*/
      if (obp->bpr) (blockp[i]->bpr)->bpl=blockp[i];/*Update LH pointer of far side - unless merging penultimate block*/

/*Update rates so that average is maintained*/
      blockp[i]->rate = (double) (data->rmap[blockp[i]->pos+blockp[i]->size]-data->rmap[blockp[i]->pos])/(loc[blockp[i]->pos+blockp[i]->size]-loc[blockp[i]->pos]);

      block2map(blockp[i],data->rmap,loc,blockp[i]->pos+blockp[i]->size);
/*	  printf("\nChecking maps agree\n");
	  for (j=0;j<=data->lseq;j++) rcheck[i]=rmap[i];
	  block2map(blockp[i],rmap,loc,data->lseq);
	  for (j=1;j<=data->lseq;j++) if (fabs(rcheck[i]-rmap[i])>0.01) {printf("\n\nError in rmap\n\n"); exit(1);}
 */
      if (DEBUG) {
		  printf("\nMerging block\n");
		  printf("\nJoining block %i to block %i\n", blockp[i]->num, obp->num);
		  printf("\nRates %f + %f ->%f",or,obp->rate,blockp[i]->rate);
		  printf("\n\nUpdated map and blocks\n");
		  print_block(blockp[0]); 
		  print_rmap(data->rmap,data->lseq,loc);
	  }
      lk_calc_site(lij,blockp[i]->pos,blockp[i]->pos+blockp[i]->size,pij,data,&lk[1],lkmat,0);

	  hr = (double) (*nblock)/((*nblock-1)*(blockp[i]->size-1));
	  lr = (double) loc[obp->pos+obp->size]-loc[obp->pos];
	  ll = (double) loc[obp->pos]-loc[blockp[i]->pos];
	  jcb = (double) pow(or+obp->rate,2)/blockp[i]->rate;
/*	  jcb = blockp[i]->rate*(ll+lr)*(ll+lr)/(ll*lr);/*For linear upwards proposal*/
/*	  jcb = (ll+lr)/(lr)*(1-exp(-blockp[i]->rate*(ll+lr)/ll))/exp(-or);/*For expontial prior upwards proposal*/


/*Acceptnce on likelihood, rates and block penalty*/
      if (log(ran2())>lk[1]+(or+obp->rate-blockp[i]->rate)+(double) data->bpen+log(hr/jcb)) {/*REJECT*/
	if (DEBUG) printf(": REJECT");
	blockp[i]->size-=obp->size;/*Subtract size of RHS*/
	blockp[i]->rate=or;/*Revert rate*/
	blockp[i]->bpr=obp;/*Revert RH pointer of Lblock*/
	if (obp->bpr) (obp->bpr)->bpl=obp;/*Revert LH pointer of far RHS*/

	block2map(blockp[i],data->rmap,loc,obp->pos+obp->size);
      }
      else {/*ACCEPT - have to rationalise block array*/
	if (DEBUG) printf(": ACCEPT\n");
	nacc[3][1]++;
	update_lij(lij,blockp[i]->pos,blockp[i]->pos+blockp[i]->size,data,0);
	for (j=obp->num;j<(*nblock)-1;j++) {blockp[j]=blockp[j+1]; blockp[j]->num--;}
	free(obp);
	(*nblock)--;
	blockp=(struct block **) realloc(blockp,(*nblock)*sizeof(struct block *));
	lk[0]+=lk[1];

      }
    }
}

/*  free(rcheck); */
  return blockp;
}


/*Routine to print blocks*/

void print_block(block0)
     struct block *block0;
{
  static int ct;

  if (!(block0->bpl)) {printf("\n\nPrinting blocks\n\n"); ct=1;}
  printf("\nBlock %3i(%3i): Rate = %8.3f, Start = %3i",ct++, block0->num,block0->rate,block0->pos);
  if (block0->bpr != NULL) print_block(block0->bpr);
  else printf("\n\nEnd of blocks\n\n");
}


/*Routine to convert blocks to recombination map*/

void block2map(block0,rmap,loc,u)
     struct block *block0; /*Pointer to block at which to start update*/
     double *rmap; /*Cumulative genetic map*/
     double *loc; /*List of SNP positions*/
     int u;/*SNP to finish updating at*/
{
  int i, usnp;

  if (!block0->bpl) rmap[block0->pos]=0.0;
  usnp = mini((int) u, block0->pos+block0->size);
  for (i=block0->pos+1;i<=usnp;i++) rmap[i] = (double) rmap[i-1]+(block0->rate)*(loc[i]-loc[i-1]);
  if ((--i)<u && block0->bpr) block2map(block0->bpr,rmap,loc,u);  /*Check to see if have reached end, otherwise call self function*/
}


/*Routine to print out recombination map*/

void print_rmap(rmap,nsnp,loc)
     double *loc;
     double *rmap;
     int nsnp;
{
  int i;

  printf("\n\nPrinting Recombination Map\n\n");
  for (i=1;i<=nsnp;i++) printf("SNP %3i: Phys pos = %9.3f : Map pos = %8.3f\n",i,loc[i], rmap[i]);
}



void check_blocks(block0,rmap,loc)
struct block *block0;
double *loc;
double *rmap;
{
	int i;

	if (block0->num==0) printf("\n\nChecking blocks........");
	else if (!block0->bpl) {printf("\nError: no BPL\n"); exit(1);}
	if (block0->rate<RMIN || block0->rate>RMAX) {printf("\nError in rate (%f)\n",block0->rate); exit(1);}

	for (i=block0->pos;i<=block0->pos+block0->size;i++) {
		if (loc[i]<=0 || loc[i]>10000000) {printf("\nError in loc (%i)\n",i); exit(1);}
		if (rmap[i]<rmap[i-1]) {printf("\nError in rmap %i(%f - %f)\n",i,rmap[i],rmap[i-1]); exit(1);}
	}

	if (!block0->bpr) printf("OK\n\n");
	else check_blocks(block0->bpr,rmap,loc);
}


void update_lij(lij,rl,ru,data,update)
double **lij;
int rl, ru, update;
struct data_sum *data;
{
	int i, j;

/*	printf("\nOld and New lk mat\n");
	for (i=1;i<=data->lseq;i++) {
		printf("\nPos %3i: ",i);
		for (j=1;j<=data->w;j++) printf(" %7.3f",lij[i][j-i]);
	}
*/
	if (update) {/*Have to update across block*/
		for (i=1;i<=ru;i++) for (j=(i<rl?rl:i+1);j<=mini(data->lseq,i+data->w); j++) lij[i][j-i]=lij[i][j-i+data->w];
/*		for (i=1;i<data->lseq;i++) for (j=i+1;j<=data->lseq;j++) lij[i][j]=lij[j][i]; */
	}
	else {/*Just update within block+associations*/
		for (i=1;i<=ru;i++) for (j=(i<rl?rl:i+1);j<=(i<rl?mini(ru,i+data->w):mini(data->lseq,i+data->w)); j++) lij[i][j-i]=lij[i][j-i+data->w];
	}

/*	printf("\nOld and New lk mat - updated from %i to %i (%i)\n",rl,ru,update);
	for (i=1;i<=data->lseq;i++) {
		printf("\nPos %3i: ");
		for (j=i+1;j<=mini(lseq,i+data->w);j++) printf(" %7.3f",lij[i][j-i]);
	}
	*/
}

void print_rates(block0,ofp)
struct block *block0;
FILE *ofp;
{
	int i;
	for (i=1;i<=block0->size;i++) fprintf(ofp,"%8.6f\t",block0->rate);
	if (block0->bpr) print_rates(block0->bpr,ofp);
	else fprintf(ofp,"\n");
}


void print_bounds(block0,ofp)
struct block *block0;
FILE *ofp;
{
	int i;

	fprintf(ofp,"\t1");
	for (i=2;i<=block0->size;i++) fprintf(ofp,"\t0");
	if (block0->bpr) print_bounds(block0->bpr,ofp);
	else fprintf(ofp,"\n");
}


void check_lk(lij,lk,lseq,w)
double **lij, lk;
int lseq, w;
{
	int i, j;
	double lkc=0.0;

	for (i=1;i<lseq;i++) for (j=i+1;j<=mini(i+w,lseq);j++) lkc+=lij[i][j-i];
	printf("\nlk = %f, expected = %f",lk,lkc);
	if (fabs(lkc-lk)>0.0001) nrerror("LK does not agree with update");
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

  if (!DEBUG) {
	  printf("\nResolving GT: "); 
	  for (i=0;i<16;i++) printf("%2i ",pset->pt[i]);
	  printf(" (%8i)",(int) (pbase[4]+1)*(pbase[5]+1)*(pbase[6]+1)*(pbase[7]+1)*(pbase[8]+1)*(pbase[8]+2)*(pbase[8]+3)*(pset->pt[15]+1)/6);
  }

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
    else for (j=1;j<=data->rcat;j++) lknew[j]+= (double) log(1+exp(mn+lkres[j]-lknew[j]));
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
  extern double *lnfac_array;
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



void print_help(int argc, char* argv[]) 
{
	int i;
	char *in_str;
	for(i = 0; i < argc; i++){
		if(*argv[i] == '-'){ 
			in_str = argv[i];
			if((strcmp(in_str, "-h") == 0) || (strcmp(in_str, "-?") == 0) || (strcmp(in_str, "-help") == 0)) {
				printf("\ninterval\n");
				printf(" Sample recombination rate profiles using composite likelihood and\npiecewise constant model.\n");
				printf(" Auton and McVean, Genome Research (2007)\n McVean et al, Science (2004)\n\n");
				printf(" Required Options :\n");
				printf(" -seq <file>               Sequence data file\n");
				printf(" -loc <file>               SNP positions data file\n");
				printf(" -lk <file>                Composite likelihood lookup data file\n");
				printf(" -its <int>                Number of MCMC iterations\n");
				printf(" -samp <int>               Number of MCMC iterations between samples\n");
				printf(" -bpen <float>             Background block penalty\n");
				printf("\n\n Additional Options :\n");
				printf(" -exact                    Likelihood file exact?\n");
				printf(" -seed <int>               User defined random seed\n");
				printf(" -prefix <string>          Prefix of output files\n");
				printf("\n\n");

				exit (0);
			}
		}
	}
}




