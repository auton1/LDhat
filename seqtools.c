#include "seqtools.h"
#include "tools.h"
#include "ldhat.h"


int read_fasta(seqs,ifp,nseq,lseq,seqnames) 
int **seqs, nseq, lseq;
char **seqnames;
FILE *ifp;
{
	int i, site, seq=0, cts[5];
	char line[MAXLINE+1], *c, bases[5]="TCAG-";
	
	printf("\nReading sequences in fasta format.\n");
	for (i=0;i<5;i++) cts[i]=0;

	while (!feof(ifp) && (seq<nseq)) 
	{
		fgets(line, MAXLINE, ifp);
		if ((c = (char *) strchr(line, '>')) != (char *)NULL) 
		{
			seq++;
			printf("Sequence :%3i ", seq);
			strncpy(seqnames[seq], (c+1), MAXNAME);
			for (i=1;i<=MAXNAME;i++) 
			{
				if ((seqnames[seq][i]=='\n') || (seqnames[seq][i]=='\r'))
					seqnames[seq][i]='\0';
			}
			printf("%s: ", seqnames[seq]);
			site=1;
			for (i=0;i<5;i++) cts[i]=0;
			while (site<=lseq) {
				fgets(line, MAXLINE, ifp);
				for (c=line; (*c)!='\0'; c++) 
				{
				   switch(*c) {
					case 'T': case 't': case '0': {
						seqs[seq][site] = 2;
						site++;
						cts[0]++;
						break;
					}
					case 'C': case 'c': case '1': {
						seqs[seq][site] = 3;
						site++;
						cts[1]++;
						break;
					}
					case 'A': case 'a': case '2': {
						seqs[seq][site] = 4;
						site++;
						cts[2]++;
						break;
					}
					case 'G': case 'g': case '3' :{
						seqs[seq][site] = 5;
						site++;
						cts[3]++;
						break;
					}
					case'-': case'N': case'n': case'?': case'R': case'Y': case'M': case'K': case'S': case'W': case'H': case'B': case'V': case'D' :{
						seqs[seq][site]=1;
						site++;
						cts[4]++;
						break;
					}
					case '>': {printf("\nError in sequence file(%i of %i bases read)\n",site,lseq); exit(1);}
					default: {
						break;
					}
				   }
				}
			}
			if (site-1 != lseq) {printf("\nSequences incorrect length (%i)\n\n",site-1); exit(1);}
			for (i=0;i<5;i++) printf("%c:%5i ",bases[i], cts[i]);
			printf("\n");
		}
	}
	if (seq!=nseq) {printf("\n\nDid not read %i sequences \n\n",nseq); exit(1);}
	return 1;
}



void allele_count(seqs,nseq,lseq,nall,fl,hd,prefix) 
int **seqs,nseq,lseq,**nall, fl, hd;
char *prefix;
{      
	int seq, site, i, ct;
	char filename[MAXNAME+1];
	FILE *ofp;
	for (site=1; site<=lseq; site++)
		for(i=1; i<=6; i++)
			nall[site][i]=0;
	for (site=1; site<=lseq; site++)
	{
		for (seq=1; seq<=nseq; seq++)
			nall[site][seqs[seq][site]]++;
		if (hd==2)
		{
			nall[site][2]=2*nall[site][2]+nall[site][4];
			nall[site][3]=2*nall[site][3]+nall[site][4];
			nall[site][4]=0;
			nall[site][5]=0;
			nall[site][1]*=2;
		}
	}
		
	if (fl)
	{
		strcpy(filename, prefix);
		ofp = fopen(strcat(filename,"freqs.txt"),"w");
		if (!ofp) nrerror("Cannot open file freqs.txt - OOM");
		fprintf(ofp,"\nAllele frequencies\n\n Site   -   T/0  C/1  A/2  G/3\n\n");
		for (site=1; site<=lseq; site++) {
			fprintf(ofp,"%4i ", site);
			for (i=1,ct=0; i<=5; i++) {fprintf(ofp,"%4i ", nall[site][i]);}

			fprintf(ofp,"\n");
		}
		fclose(ofp);
	}
}

double watterson(n) 
int n;
{

	int i;
	double cump=1.0;

	for (i=2;i<n;i++) cump+=(double) 1/i;
	return cump;
}

int check22(s1,s2,nall) 
int s1,s2,**nall;
{                      
        int i, na;

/*Commenting out this line allows the use of missing data*/
/*	if (nall[s1][1] || nall[s2][1]) return 0;*/

        for (na=0, i=2; i<=5; i++) if (nall[s1][i]>0) na++;
        if (na != 2) return 0;
        for (na=0, i=2; i<=5; i++) if (nall[s2][i]>0) na++;
        if (na != 2) return 0;
                           
        return 1; 
}



/* Routine to classify each pairwise comparison*/
struct site_type ** pair_spectrum(seqs,data,nall,pset,npt,pnew,miss,anc,pij) 
int **seqs,**nall,*npt,*pnew,anc,**pij,*miss;
struct data_sum *data;
struct site_type **pset;
{

/* pt and type (haploid)
00: 0
10: 1
01: 2
11: 3
0?: 4
1?: 5
?0: 6
?1: 7
??: 8

pt and type (diploid)

??_??:0
??_00:1
??_11:2
??_10:3
00_??:4
00_00:5
00_11:6
00_10:7
11_??:8
11_00:9
11_11:10
11_10:11
10_??:12
10_00:13
10_11:14
10_10:15
*/

	int i, j, seq, sites[2], states[2][2], *pt, nstate, ct;
	char bases[6]="n-TCAG";
	extern int sizeofpset;

	if (data->hd==1) nstate=9; else nstate=16;
	pt  = (int *) malloc((size_t) (nstate+1)*sizeof(int));

	printf("\n\n*** Classifying pairs of sites ***\n\n");

	for (sites[0]=1; sites[0]<data->lseq; sites[0]++)
	{
/*	  printf("\rSite %4i: ",sites[0]); */
		for (sites[1]=sites[0]+1; sites[1]<=mini(data->lseq,sites[0]+data->w); sites[1]++)
		{
			if (!check22(sites[0], sites[1], nall))
				pij[sites[0]][sites[1]-sites[0]]=0;
			else
			{
			   for (i=0; i<nstate; i++)
				   pt[i]=0;

			   if (data->hd==1)
			   {
			     for (j=0; j<2; j++)
			     {
			    	 for (i=2; i<=5; i++)
			    		 if (nall[sites[j]][i])
							 {states[j][0]=i;break;}
			    	 for (i++; i<=5; i++)
			    		 if (nall[sites[j]][i])
							 {states[j][1]=i;break;}
			     }
			     for (seq=1; seq<=data->nseq; seq++)
			     {
			    	 if (seqs[seq][sites[0]]==states[0][0])
			    	 {
			    		 if (seqs[seq][sites[1]]==states[1][0]) pt[0]++;
			    		 else if (seqs[seq][sites[1]]==states[1][1]) pt[2]++;
			    		 else pt[4]++;
			    	 }
			    	 else if (seqs[seq][sites[0]]==states[0][1])
			    	 {
			    		 if (seqs[seq][sites[1]]==states[1][0]) pt[1]++;
			    		 else if (seqs[seq][sites[1]]==states[1][1]) pt[3]++;
			    		 else pt[5]++;
			    	 }
			    	 else
			    	 {
			    		 if (seqs[seq][sites[1]]==states[1][0]) pt[6]++;
			    		 else if (seqs[seq][sites[1]]==states[1][1]) pt[7]++;
			    		 else pt[8]++;
			    	 }
			     }
			     if (anc == 0) pt =order_pt_hap(pt,data->nseq);
			   }

			   else
			   {
			     for (seq=1;seq<=data->nseq;seq++) pt[4*(seqs[seq][sites[0]]-1)+seqs[seq][sites[1]]-1]++; 
			     if (anc == 0) pt = order_pt_dip(pt, 2*data->nseq);
			     for (j=0,ct=0; j<nstate;j++) {ct+=pt[j];}
			     if (ct != data->nseq) nrerror("Error in calculating pair type");
			   }
			   
			   if ((*npt)+(*pnew)+(*miss)+10>sizeofpset) pset = add_pset(pset);
			   pij[sites[0]][sites[1]-sites[0]]=add_type(pset, pt, npt, pnew, miss, data);
			   if (DEBUG) printf("\nSites %i and %i: pij = %i",sites[0], sites[1], pij[sites[0]][sites[1]-sites[0]]);
			}
		}
	}
	free(pt);

	printf(" ....Done!\n\n");
	return pset;
}


/*Routine to add new pair type to existing set*/

int add_type(pset,cpt,ntc,pnew,miss,data) 
int *cpt,*ntc,*pnew, *miss;
struct site_type **pset;
struct data_sum *data;
{

	int t, fl, i, nstate=9, startp=1;
	extern int sizeofpset;
	
	if (data->hd==2) nstate=16;

/*Check that both sites are segregating*/
	if ((data->hd==1)&&((cpt[2]+cpt[3]+cpt[5]==0)||(cpt[2]+cpt[3]+cpt[7]==0))) return 0;
	else if ((data->hd==2)&&((cpt[8]+cpt[9]+cpt[10]+cpt[11]+cpt[12]+cpt[13]+cpt[14]+cpt[15]==0)||(cpt[2]+cpt[3]+cpt[6]+cpt[7]+cpt[10]+cpt[11]+cpt[14]+cpt[15]==0))) return 0;

/*If genotype data and likelihood file is not exact, start search at end of haplotype configurations*/
	if (data->hd==2 && !(data->exact)) startp=(*ntc);

	/*If next line not commented out - each new SNP pair given new type, only for GT data*/
	if (data->hd==2 && ! data->exact) startp = (*ntc)+(*pnew)+(*miss)+1;

	for (t=startp; t<=(*ntc)+(*pnew)+(*miss); t++) {
		for (i=0, fl=1; i<nstate; i++) {
			if (pset[t]->pt[i] != cpt[i]) {fl=0; break;}
		}
		if (fl==1) {
			pset[t]->nt++;
			return t;
		}
	}

	/*If gets here must be new type*/
	for (i=0; i<nstate; i++) pset[t]->pt[i]=cpt[i]; 
	if (data->hd==1) {for (i=4,fl=0;i<9;i++) if (cpt[i]) fl=1;}
	else if (cpt[0]+cpt[1]+cpt[2]+cpt[3]+cpt[4]+cpt[8]+cpt[12]) {
		fl=1;
	}
	pset[t]->nt=1;
	if (fl) {/*Missing data*/
		(*miss)++; pset[t]->miss=1;
		if (DEBUG) {
			printf("\nMissing data: ",*(miss));
			for (i=0;i<nstate;i++) printf("|%3i|",cpt[i]);
		}
	}
	else {(*pnew)++; pset[t]->miss=0; /*No missing data*/}
	return t;
}


void print_pairs(ofp,pset,nt,hd,nseq) 
int nt,hd,nseq;
FILE *ofp;
struct site_type **pset;
{

	int t, i, nstate, ct;

	char hap_types[9][3]={"00","10","01","11","0?","1?","?0","?1","??"};
	char dip_types[16][6]={"??_??","??_00","??_11","??_10","00_??","00_00","00_11","00_10","11_??","11_00","11_11","11_10","10_??","10_00","10_11","10_10"};

	if (hd==1) nstate=9;
	else nstate=16;

	fprintf(ofp,"\nPrinting pair types\n\nNum   "); 
	if (hd==1) for (i=0;i<nstate;i++) fprintf(ofp,"%5s ",hap_types[i]);
	else for (i=0;i<nstate;i++) fprintf(ofp,"%5s ",dip_types[i]);
	fprintf(ofp,"\n\n");
	for (t=1; t<=nt; t++) {
	   if (pset[t]->nt) {
		fprintf(ofp,"%5i ",t);
		for (i=0,ct=0; i<nstate; i++) {fprintf(ofp,"%4i ", pset[t]->pt[i]); ct+=pset[t]->pt[i];}
		fprintf(ofp,": %5i ",pset[t]->nt);
		fprintf(ofp, "\n");
		if (ct != nseq) {printf("\n%i %i",ct,nseq); nrerror("Error in pair types - do not match # sequences");}
	   }
	}
}

/* pt and type (haploid)
00: 0
10: 1
01: 2
11: 3
0?: 4
1?: 5
?0: 6
?1: 7
??: 8
*/

int * order_pt_hap(pt,nseq) 
int *pt, nseq;
{
	int fl=0;

	  if (2*(pt[1]+pt[3]+pt[5]) > nseq) fl += 2;
	  if (2*(pt[2]+pt[3]+pt[7]) > nseq) fl += 1;
	  switch(fl) {
		case 0 : {
			break;
		}
	  case 1 : { /*RHS only*/
			pswap(pt,0,2);
			pswap(pt,1,3);
			pswap(pt,6,7);
			break;
		}
	  case 2 : {/*LHS only*/
			pswap(pt,0,1);
			pswap(pt,2,3);
			pswap(pt,4,5);
			break;
		}
	  case 3 : {/*Both sides*/
			pswap(pt,0,3);
			pswap(pt,1,2);
			pswap(pt,4,5);
			pswap(pt,6,7);
			break;
		}
	  default : {
			printf("\n\nError in sorting\n\n");
			exit(1);
		}
	  }

/*Next Line must commented out if use multiple mutation rates*/

	  if (pt[2]+pt[7]>pt[1]+pt[5]) {
		pswap(pt,1,2);
		pswap(pt,5,7);
		pswap(pt,4,6);
	  }
	
	  return pt;
}

/*
pt and type (diploid)

??_??:0
??_00:1
??_11:2
??_10:3
00_??:4
00_00:5
00_11:6
00_10:7
11_??:8
11_00:9
11_11:10
11_10:11
10_??:12
10_00:13
10_11:14
10_10:15
*/

int * order_pt_dip(pt,nseq)
     int *pt, nseq;
{
  int fl=0;

  if (pt[12]+pt[13]+pt[14]+pt[15]+2*(pt[8]+pt[9]+pt[10]+pt[11]) > nseq) fl+=2;
  if (pt[3]+pt[7]+pt[11]+pt[15]+2*(pt[2]+pt[6]+pt[10]+pt[14]) > nseq) fl+=1;

  switch(fl) {
  case 0 : {
    break;
  }
  case 1 : {/*RHS only */
    pswap(pt,1,2);
    pswap(pt,5,6);
    pswap(pt,9,10);
    pswap(pt,13,14);
    break;
  }
  case 2 : {/*LHS only*/
    pswap(pt,4,8);
    pswap(pt,5,9);
    pswap(pt,6,10);
    pswap(pt,7,11);
    break;
  }
  case 3 : {/*Both*/
    pswap(pt,1,2);
    pswap(pt,5,10);
    pswap(pt,6,9);
    pswap(pt,7,11);
    pswap(pt,8,4);
    pswap(pt,13,14);
    break;
  }
  default : {
    printf("\n\nError in sorting\n\n");
    exit(1);
  }
  }
  if (pt[3]+pt[7]+pt[14]+2*(pt[2]+pt[6])>pt[11]+pt[12]+pt[13]+2*(pt[8]+pt[9])) {
       pswap(pt,1,4);
       pswap(pt,2,8);
       pswap(pt,3,12);
       pswap(pt,6,9);
       pswap(pt,7,13);
       pswap(pt,11,14);
  }
  return pt;
}




void type_print(pij,lseq,w,ofp) 
int **pij, lseq,w;
FILE *ofp;
{
	int i, j;
	if (!ofp) nrerror("No file to print to");

	fprintf(ofp,"Pair types\n\n        ");
        for (i=1; i<lseq; i++) fprintf(ofp, " %3i", i+1);
        for (i=1; i<lseq; i++) {
                fprintf(ofp,"\n%3i:", i);
                for (j=1; j<=i; j++) fprintf(ofp,"    ");
                for (j=i+1; j<=lseq && j-i <= w; j++) fprintf(ofp," %3i", pij[i][j-i]);
				for (;j<=lseq;j++) fprintf(ofp,"    ");
        }
}



void read_pt(ifp,pset,npt,data) 
int *npt;
FILE *ifp;
struct site_type **pset;
struct data_sum *data;
{

	int p=1, i;
	char c;

	printf("\n\n*** Reading pair types ***\n\n");

	printf("\n\nReading data for pair types\n\n");
	while((c=fgetc(ifp)) != EOF)
	if (c == '#') {
/*		printf("\rType %i",p); */
		if (!data->exact) {
		  for (i=0; i<4; i++) fscanf(ifp, "%i", &pset[p]->pt[i]);
		  for (i=4;i<16;i++) pset[p]->pt[i]=0;
		}
		else if (data->hd==1) {
		  for (i=0;i<9;i++) fscanf(ifp,  "%i", &pset[p]->pt[i]);
		  for (i=9;i<16;i++) pset[p]->pt[i]=0;
		}
		else {
		  for (i=0; i<16; i++) fscanf(ifp, "%i", &pset[p]->pt[i]);
		}
		p++;
	}
	if ((*npt) != p-1) {
		printf("\nWarning: No. entries in Likelihood file does not match total (%i)\n\n",p-1); 
		exit(1);
	}

	printf(" ....Done!\n\n");
	return;
}



struct site_type ** init_pset(pset,lkf,ifp,npt,data) 
int lkf, *npt;
FILE *ifp;
struct site_type **pset;
struct data_sum *data;
{
	int i, j, nsfile;
	struct site_type *new_pt;
	extern int sizeofpset;

	if (lkf) {
		fscanf(ifp,"%i %i", &nsfile, &(*npt));
		if (nsfile != data->nseq*data->hd) nrerror("Likelihood file for different no. seqs than data");
		sizeofpset = (*npt) + ADD;
	}

	pset = (struct site_type **) malloc((size_t) sizeofpset*sizeof(struct site_type *));
        for (i=1;i<sizeofpset;i++) {
                new_pt = (struct site_type *) malloc((size_t) sizeof(struct site_type));
                pset[i] = new_pt;
        }

        if (lkf) {
               	read_pt(ifp, pset, npt, data);
		rewind(ifp);
        }

	for (i=1;i<sizeofpset;i++) {
		pset[i]->nt=0;
		for (j=0;j<3;j++) pset[i]->ld_stat[j]=0.0;
		if ((!lkf) || (i>(*npt))) for (j=0;j<16;j++) pset[i]->pt[j]=0;
		pset[i]->miss=0;
	}
	return pset;
}



struct site_type ** add_pset(pset) 
struct site_type **pset;
{
	int i, j;
	extern int sizeofpset;
	struct site_type **npset, *new_pt;

	printf(" :New memory for pset: %i plus %i", sizeofpset, ADD); 
        npset = (struct site_type **) malloc((size_t) ((int) sizeofpset+ADD)*sizeof(struct site_type *));
        if (npset == NULL) {printf("\nError in reallocation\n"); exit(1);}
        for (i=1;i<sizeofpset;i++) npset[i]=pset[i];
        for (i=sizeofpset; i<(sizeofpset+ADD); i++) {
                        new_pt = (struct site_type *) malloc((size_t) sizeof(struct site_type));
			if (new_pt==NULL) {printf("\noom\n\n"); exit(1);}
                        for (j=0; j<16; j++) new_pt->pt[j]=0;
			for (j=0; j<3; j++) new_pt->ld_stat[j]=0.0;
                        new_pt->nt=0;
			npset[i] = new_pt;
        }
        sizeofpset += ADD; 
        free(pset);
	pset = npset;	

	return pset;
}



void read_pars(ifp,tcat,theta,rcat,rmax) 
int *tcat, *rcat;
double *theta, *rmax;
FILE *ifp;
{
	int ns, npt;

	fscanf(ifp, "%i %i", &ns, &npt);
	fscanf(ifp, "%i", &(*tcat));
	fscanf(ifp,"%lf", &(*theta));
	fscanf(ifp, "%i", &(*rcat));
	fscanf(ifp,"%lf", &(*rmax));
}


void read_lk(ifp,lkmat,npt,tcat,rcat) 
int npt, tcat, rcat;
double **lkmat;
FILE *ifp;
{
	int p=1, k;
	char c;

	printf("\n\n*** Reading likelihoods for pair types ***\n\n");
	while((c=fgetc(ifp)) != EOF)
		if (c == ':') {

	/*		printf("\rType %i", p); */
			for (k=1; k<=rcat; k++) fscanf(ifp,"%lf", &lkmat[p][k]);
			p++;
		}
	if (p-1 != npt) {printf("\n\nLikelihood file does not agree with header\n\n"); exit(1);}
	fclose(ifp);

	printf("  ...Done!\n\n");

	return;
}


