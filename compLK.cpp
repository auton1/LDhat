#include "compLK.h"

// Eventually would like to move all of this to likelihood.h

int check_exhaustive(struct site_type **pset, int npt, int nsamp)
{	//Check that likelihood file is exhaustive for n
	int p1, p2, i, ei;

	printf("Checking likelihood file is exhaustive:...");

	p1 = (int) nsamp/2;
	if (npt != (int) 1+p1+p1*(p1-1)*(p1+4)/6+(p1-1)*(p1+2)/2) 
	{
		printf("\n\n!!npt = %i: E[npt] = %i\n\n",npt, 1+p1+p1*(p1-1)*(p1+4)/6+(p1-1)*(p1+2)/2);
		printf("\nnpt and total do not agree: not exhaustive");
		return 0;
	}

	for (i=1;i<=npt;i++) 
	{
		p1 = pset[i]->pt[1]+pset[i]->pt[3];
		p2 = pset[i]->pt[2]+pset[i]->pt[3];
		ei = (int) p1*(p1-1)*(p1+4)/6+(p2-1)*(p2+2)/2+pset[i]->pt[2]+1;
		if (i!=ei) {printf("\nError in pair-types: p1(%i), p2(%i), n11(%i): i(%i) ei(%i)\n\n",p1,p2,pset[i]->pt[3],i,ei); my_exit("Not Exhaustive", 2);}
	}
	printf("OK\n");
	return 1;
}




/*Calculate likelihoods for diploid data*/
void lk_resolve(double *lkres, struct site_type *pset, double *lknew, double **lkmat,data *mydata, double *lnfac_array, int rcat)
{
  int i, j, fl, pbase[9], p1, p2, hap;
  double mn;
  struct site_type tres;

  for (i=1;i<=rcat;i++) lkres[i]=lknew[i]=0.0;

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

  for (i=0;i<=pset->pt[15];i++) {
    mn = lnfac_array[pset->pt[15]]-lnfac_array[i]-lnfac_array[pset->pt[15]-i];
    for (j=0;j<9;j++) 
		tres.pt[j]=pbase[j];
    tres.pt[0]+=i;
    tres.pt[3]+=i;
    tres.pt[1]+=pset->pt[15]-i;
    tres.pt[2]+=pset->pt[15]-i;
    order_pt_hap(tres.pt,mydata->nseq*2);

    if (fl) lk_miss(&tres,lkres,lkmat,mydata, lnfac_array, rcat);
    else {
	 p1 = (int) tres.pt[1]+tres.pt[3];
	 p2 = (int) tres.pt[2]+tres.pt[3];
	 hap = (int) p1*(p1-1)*(p1+4)/6+(p2-1)*(p2+2)/2+tres.pt[2]+1;
	 for (j=1;j<=rcat;j++) lkres[j]=lkmat[hap][j];
    }
    if (!i) for (j=1;j<=rcat;j++) lknew[j]=lkres[j];
    else for (j=1;j<=rcat;j++) lknew[j]+= (double) log(1+exp(mn+lkres[j]-lknew[j]));
  }
}



//Calculate likelihoods for missing data
void lk_miss(struct site_type *pset,double *lkmiss,double **lkmat,data *mydata, double *lnfac_array, int rcat)
{

	int j, a, b, c, d, e1, e2, e3, e4, pres[9], p1, p2, ct, k, ht;
	double cf, mn;

	cf = 0.0; for (j=1;j<=rcat;j++) lkmiss[j]=0.0;
	for (j=0;j<9;j++) pres[j]=0;
	for (k=1;k<=rcat;k++) lkmiss[k]=0.0;
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
		  order_pt_hap(pres,(int) mydata->nseq*((int) 1 + (mydata->hd==1?0:1)));
		  p1 = (int) pres[1]+pres[3];
		  p2 = (int) pres[2]+pres[3];
		  ht = (int) p1*(p1-1)*(p1+4)/6+(p2-1)*(p2+2)/2+pres[2]+1;
		  mn = lnfac_array[pset->pt[4]]-lnfac_array[a]-lnfac_array[pset->pt[4]-a]+lnfac_array[pset->pt[5]]-lnfac_array[b]-lnfac_array[pset->pt[5]-b]+\
			  lnfac_array[pset->pt[6]]-lnfac_array[c]-lnfac_array[pset->pt[6]-c]+lnfac_array[pset->pt[7]]-lnfac_array[d]-lnfac_array[pset->pt[7]-d]+\
			  lnfac_array[pset->pt[8]]-lnfac_array[e1]-lnfac_array[e2]-lnfac_array[e3]-lnfac_array[e4];

		  /*Better routine for summation*/
		  if (ct==1) {for (k=1;k<=rcat;k++) lkmiss[k] = lkmat[ht][k];}
		  else {for (k=1;k<=rcat;k++) lkmiss[k] += (double) log(1+exp(lkmat[ht][k]+mn-lkmiss[k]));}
		}
}

int check22(int s1, int s2, int **nall) 
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




/*Routine to add new pair type to existing set*/

int add_type(struct site_type **pset, int *cpt, int *ntc, int *pnew, int *miss, data *mydata) 
{

	int t, fl=1, i, nstate=9, startp=1;
	extern int sizeofpset;
	
	if (mydata->hd==2) nstate=16;

/*Check that both sites are segregating*/
	if ((mydata->hd==1)&&((cpt[2]+cpt[3]+cpt[5]==0)||(cpt[2]+cpt[3]+cpt[7]==0))) return 0;
	else if ((mydata->hd==2)&&((cpt[8]+cpt[9]+cpt[10]+cpt[11]+cpt[12]+cpt[13]+cpt[14]+cpt[15]==0)||(cpt[2]+cpt[3]+cpt[6]+cpt[7]+cpt[10]+cpt[11]+cpt[14]+cpt[15]==0))) return 0;

/*If genotype data and likelihood file is not exact, start search at end of haplotype configurations*/
	if (mydata->hd==2 && !(mydata->lk_exact)) startp=(*ntc);

	for (t=1; t<=(*ntc)+(*pnew)+(*miss); t++) {
		for (i=0, fl=1; i<nstate; i++) {
			if (pset[t]->pt[i] != cpt[i]) {fl=0; break;}
		}
		if (fl==1) {
			pset[t]->nt++;
			return t;
		}
	}
	for (i=0; i<nstate; i++) pset[t]->pt[i]=cpt[i]; 
	if (mydata->hd==1) {for (i=4,fl=0;i<9;i++) if (cpt[i]) fl=1;}
	else if (cpt[0]+cpt[1]+cpt[2]+cpt[3]+cpt[4]+cpt[8]+cpt[12]) {
		fl=1;
	}
	pset[t]->nt=1;
	if (fl) {
		(*miss)++; //pset[t]->miss=1;
	}
	else {(*pnew)++; 
	//pset[t]->miss=0;
	}
	return t;
}


void print_pairs(FILE *ofp, struct site_type **pset, int nt, int hd, int nseq) 
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



int * order_pt_hap(int *pt, int nseq) 
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
			my_exit("Error in sorting", 5);
		}
	  }

/*Next Line must commented out if use multiple mutation rates*/

	  if (pt[2]+pt[7]>pt[1]+pt[5]) 
	  {
		pswap(pt,1,2);
		pswap(pt,5,7);
		pswap(pt,4,6);
	  }
	
	  return pt;
}

int * order_pt_dip(int *pt, int nseq)
{
  int fl=0;

  if (pt[12]+pt[13]+pt[14]+pt[15]+2*(pt[8]+pt[9]+pt[10]+pt[11]) > nseq) fl+=2;
  if (pt[3]+pt[7]+pt[11]+pt[15]+2*(pt[2]+pt[6]+pt[10]+pt[14]) > nseq) fl+=1;

  switch(fl) 
  {
  case 0 : 	break;
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
    my_exit("Error in sorting", 6);
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




void type_print(int **pij, int lseq, int w, FILE *ofp) 
{
	int i, j;
	if (!ofp) nrerror("No file to print to");

	fprintf(ofp,"Pair types\n\n        ");
        for (i=1; i<lseq; i++) fprintf(ofp, " %3i", i+1);
        for (i=1; i<lseq; i++) 
		{
                fprintf(ofp,"\n%3i:", i);
                for (j=1; j<=i; j++) fprintf(ofp,"    ");
                for (j=i+1; j<=lseq && j-i <= w; j++) fprintf(ofp," %3i", pij[i][j-i]);
				for (;j<=lseq;j++) fprintf(ofp,"    ");
        }
}



void read_pt(FILE *ifp, struct site_type **pset, int *npt, data *mydata) 
{

	int p=1, i;
	int c;

	printf("Reading data for pair types\n\n");
	while((c=fgetc(ifp)) != EOF)
	if (c == '#') {
		//if ((p%50) == 0) printf("\n");
		//printf(".");
		if (!mydata->lk_exact) {
		  for (i=0; i<4; i++) fscanf(ifp, "%i", &pset[p]->pt[i]);
		  for (i=4;i<16;i++) pset[p]->pt[i]=0;
		}
		else if (mydata->hd==1) {
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
		my_exit("Error in read_pt", 7);
	}
}



struct site_type ** init_pset(struct site_type **pset, int lkf, FILE *ifp, int *npt, data *mydata) 
{

	int i, j, nsfile;
	struct site_type *new_pt;
	extern int sizeofpset;

	if (lkf) 
	{
		fscanf(ifp,"%i %i", &nsfile, &(*npt));
		//if (nsfile != (mydata->nseq)*((int)(mydata->hd==1?1:2)))	// Possible Bug Here?
		if (nsfile != (mydata->nseq)*((int)(mydata->hd)))	// Possible Bug Here?
			nrerror("Likelihood file for different no. seqs than data");
		sizeofpset = (*npt) + ADD;
	}

	//pset = (struct site_type **) malloc((size_t) sizeofpset*sizeof(struct site_type *));
	pset = new site_type*[sizeofpset];
	if (pset == NULL) nrerror("pset");
        for (i=1;i<sizeofpset;i++) 
		{
                //new_pt = (struct site_type *) malloc((size_t) sizeof(struct site_type));
				new_pt = new site_type;
				if (new_pt == NULL) nrerror("new_pt");
                pset[i] = new_pt;
        }

        if (lkf) 
		{
            read_pt(ifp, pset, npt, mydata);
			rewind(ifp);
        }

	for (i=1;i<sizeofpset;i++) 
	{
		pset[i]->nt=0;
//		for (j=0;j<3;j++) pset[i]->ld_stat[j]=0.0;
		if ((!lkf) || (i>(*npt))) for (j=0;j<16;j++) pset[i]->pt[j]=0;
//		pset[i]->miss=0;
	}
	return pset;
}



struct site_type ** add_pset(struct site_type **pset) 
{

	int i, j;
	extern int sizeofpset;
	struct site_type **npset, *new_pt;

	printf("New memory for pset: %i plus %i\n", sizeofpset, ADD); 
        //npset = (struct site_type **) malloc((size_t) ((int) sizeofpset+ADD)*sizeof(struct site_type *));
		npset = new site_type*[sizeofpset+ADD];
		if (npset == NULL) nrerror("npset");
        for (i=1;i<sizeofpset;i++) npset[i]=pset[i];
        for (i=sizeofpset; i<(sizeofpset+ADD); i++) 
		{
			//new_pt = (struct site_type *) malloc((size_t) sizeof(struct site_type));
			new_pt = new site_type;
			if (new_pt==NULL) nrerror("new_pt");
			for (j=0; j<16; j++) new_pt->pt[j]=0;
			new_pt->nt=0;
			npset[i] = new_pt;
        }
        sizeofpset += ADD; 
        //free(pset);
		delete pset;
	pset = npset;	

	return pset;
}



void read_pars(FILE *ifp, int *tcat, double *theta, int *rcat, double *rmax) 
{
	int ns, npt;
	double temp;

	fscanf(ifp, "%i %i", &ns, &npt);
	fscanf(ifp, "%i", &(*tcat));
	fscanf(ifp,"%lf", &(temp));
	*theta = (double)temp;
	fscanf(ifp, "%i", &(*rcat));
	fscanf(ifp,"%lf", &(temp));
	*rmax = (double)temp;
}


