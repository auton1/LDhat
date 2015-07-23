#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "tools.h"
#include "gamma.h"
#include "fin.h"

long int *idum;
int tree_size;

#define DEBUG 0
#define NPANEL 8
#define NCP 1.1
#define SMC 0
#define TREEPRINT 0
#define PHASE 0

main(int argc, char *argv[])
{
	int nmut, nrec, run, i, **segl, site, **seqs;
	int na_c=0, nac=0, fl=0;
	double rtot, stot;
	extern int tree_size;
	struct node_tree ***tree_ptr, *new_node;
	struct results res;
	struct control con;
	int tmp;

	read_input(&con, argc, argv);

	printf("\n\nCoalescent simulations with mutations and recombination: After Hudson (1991)\n\n");
	if (con.asc>0) con.nsamp += (int) con.asc;

	/*Define tree size and nodes*/
	tree_size = 2*con.nsamp-1;
	tree_ptr = (struct node_tree ***) malloc((size_t) (con.len+1)*sizeof(struct node_tree **));
	for (site=1;site<=con.len;site++) {
		tree_ptr[site] = (struct node_tree **) malloc((size_t) (tree_size+1)*sizeof(struct node_tree *));
		for (i=1;i<=tree_size;i++) {
			new_node = (struct node_tree *) malloc((size_t) sizeof(struct node_tree));
			new_node->d[0]=NULL; new_node->d[1]=NULL;
			new_node->a=NULL; 
			new_node->site=site;
			new_node->time_above=new_node->time_below=0.0;
			if (i<=con.nsamp) {new_node->time=0.0; new_node->node_num=i; new_node->ndesc=1;}
			tree_ptr[site][i] = new_node;
		}
	}

	segl = imatrix(1,con.len,1,2);
	/*Keeps track of how many coalescent events for each site: can also use to tell which node to use next*/

	res.fdist=dvector(1,con.nsamp);
	set_res(&res, &con);
	seqs = imatrix(1,con.nsamp,1,con.len);

	for (run=0, rtot=0, stot=0; run<con.nrun; run++) {
		if ((run%1000)==0)
			printf("\nRun %i\n", run);
		if(DEBUG) printf("\nrun %4i\n", run);
		nmut=nrec;
		tree_ptr = make_tree(tree_ptr, &con, &nrec, segl);
		res.nr+=(double) nrec;
		tree_summary(tree_ptr, &con, &res, seqs);
	}

	printf("\nFinished simulations\n");
	print_res(res, con);

	for (site=1;site<=con.len;site++) {
		for (i=1; i<=tree_size; i++) free(tree_ptr[site][i]);
		free(tree_ptr[site]);
	}
	free(tree_ptr);
	
	free_imatrix(segl,1,con.len,1,2);
	free_imatrix(seqs,1,con.nsamp,1,con.len);
}




/*****************************/
/* Initialise results        */
/*****************************/

void set_res(res, con)
struct results *res;
struct control *con;
{
	int i;
	res->nm = res->nr = res->sn = res->pwd = res->mhm = 0;
	for (i=1;i<=con->nsamp;i++) res->fdist[i]=0.0;
	res->cf=res->cf2=0.0;
}






/****************************************/
/*The key routine: creates the genealogy*/
/****************************************/

struct node_tree ***  make_tree(tree_ptr, con, revts, segl) 
struct node_tree ***tree_ptr;
struct control *con;
int *revts, **segl;
{
	int k, i, j, nrec=0, nco=0, ncb, k_eff, fl, pos;
	extern int tree_size;
	double rr, avR, cump, t=0, told;
	struct node_list **list, *new_node_list;
	char fname[128];
	FILE *ofp;

	list = (struct node_list **) malloc((size_t)((con->nsamp+1)*sizeof(struct node_list *)));
	for (i=1;i<=con->nsamp;i++) {
		new_node_list = (struct node_list *) malloc((size_t) sizeof(struct node_list));
		new_node_list->asite = (struct node_tree **) malloc((size_t) (con->len+1)*sizeof(struct node_tree *));
		for (j=1;j<=con->len;j++) new_node_list->asite[j]=tree_ptr[j][i]; /*Start pointers at base of trees at each site*/
		new_node_list->nanc=1;
		new_node_list->rlen = (double) con->R;
		new_node_list->time=0.0;
		list[i]=new_node_list;
	}
	for (i=1;i<=con->len;i++)  {segl[i][1]=con->nsamp;}

	if (DEBUG>2) {
	  printf("\n\nInitial state of sample\n\n");
	  print_lin(list, con->len, con->nsamp);
	}

	k=con->nsamp;	/*Number of lineages active in list*/
	if (EVENT_PRINT)
	{
		strcpy(fname, con->prefix);
		ofp = fopen(strcat(fname, "events.txt"), "w");
	}
	
	while (k>1) {
	
/*		check_lin(list, &k, con);
 		printf("\n%i lineages currently active\n", k); */


	/*Count lineages which can contribute to rec events*/
		for(i=1, rr=0; i<=k; i++) {
			rr += list[i]->rlen;
		}			
		avR = (double) rr/k;
		if (DEBUG) printf("\nTotal recombination length = %.0f , Re = %.3f", rr, avR);

	/*Count pairs of lineages that can coalesce under the SMC - only needed if printing*/
		if (EVENT_PRINT) {
			if (SMC) {
				for (i=1,k_eff=0;i<k;i++) for (j=i+1;j<=k;j++) {
					for (pos=1,fl=0;pos<=con->len;pos++) if (list[i]->asite[pos] && list[j]->asite[pos]) {fl=1;break;}
					if (fl) k_eff++;
				}
			}
			else k_eff = k*(k-1)/2;
			fprintf(ofp,"%i\t%i\t%.3f\t%.3f\t",k, k_eff, (double) avR/2, -log((double) k_eff+avR/2));
		}

	/*Choose time for next event*/
		told=t;
		choose_time(&t, k, avR, con);
		if (DEBUG>2) printf("\nNext event at time %.3f\n", t);

	/*Choose type of next event: update list and tree accordingly*/
		if ((con->bneck)&&(told<con->tb)&&(t>con->tb)) { /*Bottleneck*/
			t=con->tb;
			told = con->tb+0.0001;
			for(cump=0,ncb=0;cump<con->strb;) {
				cump+=(double) -log(ran2())*2/(k*(k-1));
				if (cump<con->strb) {
					k--;nco++; ncb++;
					coalesce(&k,con,list,tree_ptr,segl,&t);
				}
			}
			if (DEBUG) printf("\n\nA total of %i coalescent events occurred during the bottlneck leaving %i\n\n",ncb,k);
		}


		else if (ran2() < (double) avR/(avR+(k-1)*exp((con->lambda)*t))) {  /*Event is a recombination*/
			k++; nrec++; 
			list = recombine(&k, &rr, con, list, &t, ofp);
		}

		else { /*Event is a coalescent*/
			nco++; k--; 
			list = coalesce(&k,con,list,tree_ptr,segl,&t, ofp);
		}

		/*		if (DEBUG) print_lin(list, con->len, k);*/
	}

/*	printf("\nEnd of run: %.3lf Nrec: %i Nco: %i\n", t, nrec, nco); */

	if (EVENT_PRINT) fclose(ofp);

	if (DEBUG > 0)	print_nodes(tree_ptr[1], con->len);
	free(list);
	(*revts) = nrec;

	return tree_ptr;
}




/*******************************/
/* Generate coalescent         */
/*******************************/

struct node_list ** coalesce(k,con,list,tree_ptr,segl,t,ofp)
int *k, **segl;
double *t;
struct control *con;
struct node_list **list;
struct node_tree ***tree_ptr;
FILE *ofp;
{
	int i, j, pos, node_pos, swap, fl;
	static int tally=0;
	extern int tree_size;

	i=1+((*k)+1)*ran2(); /*Choose first lineage to coalesce*/
	j=i;
	while(j==i)
		j = 1 +((*k)+1)*ran2(); /*Choose second lineage to coalesce*/
	if (DEBUG)
		printf("Lineages %i and %i coalesce", i,j);

	if (SMC)
	{
		for (pos=1,fl=0;pos<=con->len;pos++)
			if (list[i]->asite[pos] && list[j]->asite[pos])
			{
				fl=1;break;
			}
		if (!fl)
		{
			(*k)++;
			if (EVENT_PRINT)
				fprintf(ofp,"0\n");
			return list;
		}
	}

	if (EVENT_PRINT) fprintf(ofp,"1\n");

	/*Find all sites for which coalescent event has occured and update tree accordingly*/
	for (pos=1;pos<=con->len;pos++)
	{
		if (list[i]->asite[pos] && list[j]->asite[pos])
		{/*Coalescent event at site*/
			segl[pos][1]--;
			node_pos = 2*con->nsamp-segl[pos][1]; /*Number in tree_ptr of next unused pointer for pos*/
			tree_ptr[pos][node_pos]->d[0]=list[i]->asite[pos];
			tree_ptr[pos][node_pos]->d[1]=list[j]->asite[pos];
			(tree_ptr[pos][node_pos]->d[0])->a=tree_ptr[pos][node_pos];
			(tree_ptr[pos][node_pos]->d[0])->time_above=(*t)-(tree_ptr[pos][node_pos]->d[0])->time;
			(tree_ptr[pos][node_pos]->d[1])->a=tree_ptr[pos][node_pos];
			(tree_ptr[pos][node_pos]->d[1])->time_above=(*t)-(tree_ptr[pos][node_pos]->d[1])->time;
			tree_ptr[pos][node_pos]->time = (*t);
			tree_ptr[pos][node_pos]->time_below = (double) (tree_ptr[pos][node_pos]->d[0])->time_below+\
				(tree_ptr[pos][node_pos]->d[0])->time_above+(tree_ptr[pos][node_pos]->d[1])->time_below+\
				(tree_ptr[pos][node_pos]->d[1])->time_above;
			tree_ptr[pos][node_pos]->ndesc = (tree_ptr[pos][node_pos]->d[0])->ndesc+(tree_ptr[pos][node_pos]->d[1])->ndesc;
			tree_ptr[pos][node_pos]->node_num = node_pos;
			if (segl[pos][1]>1)
				list[j]->asite[pos]=tree_ptr[pos][node_pos];
			else
			{
				list[j]->asite[pos]=NULL;
				/*printf("\rSite %5i MRCA : tally = %i",pos, ++tally);*/
			}
		}
		else
			list[j]->asite[pos] = (struct node_tree *) (list[i]->asite[pos]?list[i]->asite[pos]:list[j]->asite[pos]); /*Just copy existing*/
	}
	count_rlen(list[j], con);
	list[j]->time=(*t);

	/*Reallocate list - check to see if can remove new node because all sites found MRCA*/
	if (!(list[j]->nanc)) {
		free(list[i]->asite);
		free(list[i]);
		free(list[j]->asite);
		free(list[j]);
		if (j==(*k)) {
			list[i]=list[(*k)+1];
			list[j]=list[(*k)];
		}
		else {
			list[i]=list[(*k)];
			list[j]=list[(*k)+1];
		}
		(*k)--;
	}
	else {
		free(list[i]->asite);
		free(list[i]);
		list[i]=list[(*k)+1];/*Replace i with last in (previous) list*/
	}
	list = (struct node_list **) realloc(list, (size_t) (*k+1)*sizeof(struct node_list *));
	return list;
}


/***********************/
/*Genereate recombinant*/
/***********************/

struct node_list ** recombine(k, rr, con, list, t, ofp) 
int *k;
double *t, *rr;
struct control *con;
struct node_list **list;
FILE *ofp;
{

	int i,j,j0, upper=con->len;
	double cump, ri;
	struct node_list *new_node_list;

	list = (struct node_list **) realloc(list, (size_t) (*k+1)*sizeof(struct node_list *)); /*Resize: use k+1th pos to fill in*/
	new_node_list = (struct node_list *) malloc((size_t) sizeof(struct node_list));
	new_node_list->asite = (struct node_tree **) malloc((size_t) (con->len+1)*sizeof(struct node_tree *));
	list[(*k)]=new_node_list;
	
	/*Choose lineage to recombine*/
	ri = (double) ran2()*(*rr);
	for (i=1, cump=0; cump<ri; i++)
		cump += list[i]->rlen;
	i--;
	if (DEBUG > 1) printf("\nLineage %i recombines at position ", i);
	ri = (double) ran2()*(list[i]->rlen);  /*Choose point along sequence to recombine*/
	for (j0=1;(list[i]->asite[j0])==NULL;j0++); /*Find beginning of ancestral material*/

	for (cump=0,j=j0;cump<ri;j++)
		cump=con->rmap[j]-con->rmap[j0]; /*Choose last point on LHS to come from parent 1*/
	j0=j-2;

	if (DEBUG > 1) printf("%i\n", j0);
	for (j=1; j<=j0; j++) {list[(*k)]->asite[j]=NULL;}/*Copy LHS to parent*/

	if (EVENT_PRINT) fprintf(ofp,"%.3f\n",(double) log((con->rmap[j0+1]-con->rmap[j0])/2));

	/*Choose if crossing-over or gene conversion event and define upper point to copy until*/

	if (con->conv && ran2()<con->cratio/(1+con->cratio)) upper = (int) ((double) j0 - log(ran2())*con->clen);
	upper = mini(upper, con->len);
	for (;j<=upper; j++)
	{
		list[(*k)]->asite[j]=list[i]->asite[j];
		list[i]->asite[j]=NULL;
	}
	for (;j<=con->len;j++)
	{
		list[(*k)]->asite[j]=NULL;
	}

	list[i]->time=list[(*k)]->time=(*t);

	count_rlen(list[i],con);
	count_rlen(list[(*k)],con);			

	return list;
}



/*********************************************************/
/*Prints a list of lineages still active in the genealogy*/
/*********************************************************/

void print_lin(list, len, k) 
struct node_list **list;
int len, k;
{
	int i, j;

	printf("\n\nLineages active (%i)\n\n", k);
	for (i=1; i<=k; i++) {
		printf("lineage %3i: time %.3f: rlen %6.2f :ANC ",\
			i, list[i]->time, list[i]->rlen);
		for (j=1; j<=len; j++) printf("%i", (list[i]->asite[j]!=NULL?1:0));
		printf("\n");
	}
	printf("\n\n");
}


/******************************************************/
/*Print nodes in tree: useful for checking simulations*/
/******************************************************/

void print_nodes(tree, len) 
struct node_tree **tree;
int len;
{
	int n;
	extern int tree_size;

	printf("\n\nPrinting nodes in tree\n\n");

	for (n=1; n<=tree_size; n++) {
		printf("node %3i: time %.3f ", tree[n]->node_num, tree[n]->time);
		if (tree[n]->d[0] != NULL) printf(" : DESC %i",(tree[n]->d[0])->node_num);
		if (tree[n]->d[1] != NULL) printf(" %i", (tree[n]->d[1])->node_num);
		if (tree[n]->a != NULL) printf(" : ANC %i", (tree[n]->a)->node_num);
		printf("\n");
	}
	printf("\n\n");
}


/***************************************/
/*Calculate potential for recombination*/
/***************************************/

void count_rlen(nodel, con) 
struct node_list *nodel;
struct control *con;
{
	int i;

	nodel->rlen=0.0;
	nodel->nanc=1;
	for (i=1; i<=con->len; i++)
		if (nodel->asite[i])
			break;
	if (i>con->len)
	{ /*If no longer any ancestral material*/
		nodel->rlen = 0;
		nodel->nanc=0;
		return;
	}
	nodel->rlen -= (double) con->rmap[i];
	for (i=con->len; i>0; i--)
		if (nodel->asite[i])
			break;
	nodel->rlen += (double) con->rmap[i];
}


/******************************************************************/
/*Add mutations to tree and calculate summary statistics of sample*/
/******************************************************************/

void tree_summary(tree, con, res, seqs) 
int **seqs;
struct node_tree ***tree;
struct control *con;
struct results *res;
{

	int site, mrca, nmuts, i, j, fi, sn;
	int mhm, fl=1, fl_asc;
	char bases[5]="NTCAG";
	double tree_len=0, r1, ft, pwd=0, *f, **mutmat, *mm, cump, tt;
	double cf;
	struct node_tree *node_mut;
	char fname[128];
	FILE *ifp, *ofp;

/*	printf("\n\nTree summary");*/ 

	for (i=1; i<=con->nsamp; i++) for (j=1; j<=con->len; j++) seqs[i][j]=0;

	if (TREEPRINT) {
		printf("\n\nInput position to print tree:");
		scanf("%i", &site);
		strcpy(fname, con->prefix);
		ofp = fopen(strcat(fname, "tree.txt"), "w");
		fprintf(ofp,"Node\tTime\tLpos\tRpos\tD1\tD2\tNmut");
		for (i=1;i<=tree_size;i++) {
			fprintf(ofp,"\n%i\t%.3f\t1\t%i",tree[site][i]->node_num, tree[site][i]->time, con->len);
			if (tree[site][i]->d[0]) fprintf(ofp,"\t%i\t%i\t0",tree[site][i]->d[0]->node_num, tree[site][i]->d[1]->node_num);
			else fprintf(ofp,"\t0\t0\t0");
		}
		fclose(ofp);
	}

/*	Read user-defined mutation model*/

	if (con->mut > 0) {
	  mm=dvector(1,4);
	  f = dvector(1,4);
	  mutmat=dmatrix(1,4,1,4);
	  printf("\nReading input data for mutation model\n\n");
	  ifp = fopen(con->mut_file, "r");
	  if (!ifp) {printf("\nCannot open mutation file\n\n"); exit(1);}
	  for (i=1, ft=0.0;i<=4;i++) {fscanf(ifp,"%lf", &f[i]); ft+=f[i];}
	  if (fabs(ft-1.0)>0.001) {printf("\nAllele freqs do not sum to 1\n"); exit(1);}
	  for (i=1, ft=0.0; i<=4; i++) for (j=1, mm[i]=0.0;j<=4;j++){
		fscanf(ifp,"%lf", &mutmat[i][j]);
		mm[i]+=mutmat[i][j];
	  }
	  for (i=1, ft=0.0; i<=4; i++) ft += f[i]*mm[i];
	  for (i=1; i<=4; i++) {
		for (j=1;j<=4;j++) mutmat[i][j]/=ft; 
		mm[i]/=ft;
	  }
	  printf("\nMutation matrix\n\n     ");
	  for (i=1;i<=4;i++) printf(" %5c",bases[i]);
	  for (i=1; i<=4;i++) {
	    printf("\n%4c:", bases[i]);
	    for (j=1;j<=4;j++) {printf(" %.3f", mutmat[i][j]); mutmat[i][j]/=mm[i];}
	  }
	  fclose(ifp);
	}	


/*	Add mutations to tree	*/

	for (site=1, cf=0.0; site<=con->len; site++) {
		nmuts= (con->cond ? 1 : 0);			/*Default: if no flocs file, assume all segregate*/
		mrca = tree_size;
		if (DEBUG>1) printf("\nSite %3i : mrca = %3i", site, mrca);
		if (con->mut) {/*Routines for user defined mutation matrix*/
			r1=ran2(); 
			for (i=1, cump=0.0;cump<r1;i++) cump+=f[i]; i--;
			if (DEBUG>1) printf("\nAncestral base = %c\n", bases[i]);
			tree[site][mrca]->nuc=i;
			seq_mut(tree[site][mrca], seqs, site, i);
			evolve(tree[site][mrca], seqs, con, mm, mutmat, &nmuts, site);
		}
		else {
			tree_len = tree[site][mrca]->time_below;
			if (con->cond) {
			        nmuts=con->slocs[site][1];
				if (nmuts>1) nmuts=1;		/*Cannot condition on more than one segregating site*/
				if (!nmuts) cf += (double) -(tree_len/2-(con->w0))*(con->theta[site]);
				else {
				  if (!(con->slocs[site][2])) {
					cf+=(double) log(tree_len/(2*con->w0))-(tree_len/2-(con->w0))*(con->theta[site]);
					tt = (double) tree_len*ran2();
					node_mut = tree[site][add_mut(tree[site], &tt)];
					if (DEBUG>1) printf(" : added to node %i\n", node_mut->node_num);
				  }
				  else node_mut = add_mut_f(con->slocs[site][2], con, &cf, tree[site]);
				  if (DEBUG>1) printf(" : mutation added to node %i\n", node_mut->node_num);
				  seq_mut(node_mut, seqs, site, 1);
				}
			}
			else {
				nmuts=rpoiss((con->theta[site])*tree_len/2);
				if ((con->inf)&&(nmuts))
					nmuts=1;
				for (i=0; i<nmuts; i++) {
				  tt = (double) ((double) tree_len*ran2());
				  j = add_mut(tree[site], &tt);
				  node_mut = tree[site][j];
				  if (DEBUG>1)
					  printf(" : added to node %i: time interval = %.3f - %.3f\n", node_mut->node_num, node_mut->time, (node_mut->a)->time);
				  /*printf("site=%i\n", site);*/
				  seq_mut(node_mut, seqs, site, i+1);
				}
				cf = 0.0;
			}
		}
	(*res).nm += (double) nmuts;
	}

	/*Ascertainment routine - wipes out any segregating sites not in ascertainment set*/

	if (con->asc > 0) {/*Look to see if segregating in first n sequences: if not, wipes out site*/
	  for (site=1;site<=con->len;site++) {
	    for (i=2,fl_asc=0;i<=con->asc;i++) {if (seqs[i][site]!=seqs[1][site]) fl_asc=1;}
	    if (!fl_asc) for (i=1;i<=con->nsamp;i++) seqs[i][site]=seqs[1][site];
	  }
	}
	
	if (con->p_genotype_error > 0)
	{
		printf("Adding genotype error (genoErr=%f)\n", con->p_genotype_error);
		add_genotype_error(seqs, con);
	}


	if (con->p_switch_error > 0)
	{
		printf("Adding switch error (switchErr=%f)\n", con->p_switch_error);
		add_switch_error(seqs, con);
	}

	if (con->p_bad_site > 0)
	{
		printf("Adding site error (siteErr=%f)\n", con->p_bad_site);
		add_bad_sites(seqs, con);
	}

	if (con->fmin > 0)
	{
		printf("Filtering sites by frequency (f=%i)\n", con->fmin);
		remove_sites_by_frequency(seqs, con);
	}

	if (con->print) print_seqs(seqs, con);

/*Summary statistics for polymorphism*/
	for (site=1, sn=0, mhm=0; site<=con->len; site++)
	{
		for (i=(int) (con->asc==2 ? NPANEL+3 : 2), fi=0; i<=con->nsamp; i++) 
			if (seqs[i][site]!=seqs[(int) (con->asc==2 ? NPANEL+2 : 1)][site])
				fi++;
		if (fi && (fi!=con->nsamp))
		{
			sn++; 
			if ((con->theta[site])>(con->theta[1]))
				mhm++;
			if (fl)
				res->fdist[fi]++;
		}
		pwd += (double) fi*((con->nsamp-((int) con->asc==2?NPANEL+1:0))-fi);
	}
	pwd *= (double) 2/((con->nsamp-((int) con->asc==2?NPANEL+1:0))*((con->nsamp-((int) con->asc==2?NPANEL+1:0))-1));

	if (fl)
	{
	  (*res).sn += (double) sn*exp(cf);
	  (*res).pwd += (double) pwd*exp(cf);
	  if (sn) (*res).mhm += (double) mhm/sn*exp(cf);
	  (*res).cf+=(double) exp(cf);
	  (*res).cf2+=(double) exp(2*cf);
	  if (DEBUG) printf("\nWeighting factor for sample = %.5f\n", exp(cf));
	}

	if (con->mut)
	{
	  free_dvector(f,1,4);
	  free_dvector(mm,1,4);
	  free_dmatrix(mutmat,1,4,1,4);
	}
}

int add_genotype_error(int **seqs, struct control *con)
{
	int i, j, site, fi;
	int *seg;
	double u;
	int count=0;

	seg = ivector(1,con->len);
	for (i=1; i<=con->len; i++)
	{
		fi = 0;
		for (j=2; j<=con->nsamp; j++)
			if (seqs[j][i]!= seqs[1][i])
				fi++;
		if (fi != 0)
			seg[i]=1;
		else
			seg[i]=0;
	}

	for (site=1;site<=con->len;site++)
	{
		if (seg[site] == 1)
		{
			for (i=1;i<=con->nsamp;i++)
			{
				u = ran2();
				if (u < con->p_genotype_error)
				{
					seqs[i][site] = !seqs[i][site];
					count++;
				}
			}
		}
	}

	free_ivector(seg,1,con->len);
	return count;
}

int add_bad_sites(int **seqs, struct control *con)
{
	int i, site, j;
	double u;
	int tmp;
	int count=0;
	for (site=1;site<=con->len;site++)
	{
		u = ran2();
		if (u < con->p_bad_site)
		{
			for (i=con->nsamp-1;i>0;i--)
			{
				j = (int)((i * ran2()) + 1.5);
				tmp = seqs[i][site];
				seqs[i][site] = seqs[j][site];
				seqs[j][site] = tmp;
			}
			count++;
		}
	}
	return count;
}

int add_switch_error(int **seqs, struct control *con)
{
	int i, j, site, site2, fi;
	double u;
	int tmp;
	int *seg;
	int count=0;

	seg = ivector(1,con->len);
	for (i=1; i<=con->len; i++)
	{
		fi = 0;
		for (j=2; j<=con->nsamp; j++)
			if (seqs[j][i]!= seqs[1][i])
				fi++;
		if (fi != 0)
			seg[i]=1;
		else
			seg[i]=0;
	}

	for (i=1;i<=con->nsamp;i+=2)
	{
		for (site=1;site<con->len;site++)
		{
			if (seg[site] == 1)
			{
				u = ran2();
				if (u < con->p_switch_error)
				{
					for (site2=site+1;site2<=con->len;site2++)
					{
						tmp = seqs[i][site2];
						seqs[i][site2] = seqs[i+1][site2];
						seqs[i+1][site2] = tmp;
					}
					count++;
				}
			}
		}
	}


	free_ivector(seg,1,con->len);
	return count;
}

int remove_sites_by_frequency(int **seqs, struct control *con)
{
	int i, site;
	int count=0, freq;
	for (site=1;site<=con->len;site++)
	{
		freq = 0;
		for (i=1;i<=con->nsamp;i++)
			freq += seqs[i][site];

		if ((freq < con->fmin) || (con->nsamp - freq < con->fmin))
		{
			for (i=2;i<=con->nsamp;i++)
				seqs[i][site] = seqs[1][site];
			count++;
		}

	}
	return count;
}



/*************************************************/
/*Routine to find total tree length for each site*/
/*************************************************/

double tree_time(node) 
struct node_tree *node;
{

	double ttime=0;

	if (node->d[0] == NULL) ttime = node->time; /*Is a terminal node*/
	else {
		ttime += (node->time) - (node->d[0])->time;
		ttime += tree_time(node->d[0]);
		ttime += (node->time) - (node->d[1])->time;
		ttime += tree_time(node->d[1]);
	}
	return ttime;
}

/*********************************************/
/*Routine to place mutations on the genealogy*/
/*********************************************/

int add_mut(tree, fl) 
struct node_tree **tree;
double *fl;
{
	int i;
	static int nn=0;
	extern int tree_size;

	for (i=1;i<=tree_size &  *fl > 0.0;i++) *fl -= tree[i]->time_above;
	return (--i);

	/*

	if (!node->d[0]) return nn;
	else {
		if ((*fl)<0) return nn;
		(*fl) -= node->time-(node->d[0])->time;
		if ((*fl)<0) nn=(node->d[0])->node_num;
		else nn = add_mut(node->d[0], fl);

		if ((*fl)<0) return nn;
		(*fl) -= node->time-(node->d[1])->time;
		if ((*fl)<0) nn=(node->d[1])->node_num;
		else nn = add_mut(node->d[1], fl);

	}
	*/
}

/******************************************************************************/
/*Routine to place mutations on the genealogy conditioning on allele frequency*/
/*fsim is an integer with the MAF                                             */
/******************************************************************************/

struct node_tree * add_mut_f(fsim, con, cf, tree_site) 
struct node_tree **tree_site;
struct control *con;
double *cf;
int fsim;
{

	int i, maf;
	double *flist, cump, r1, part[3], mx=0.0;
	static double cons=0.0;
	extern int tree_size;

	fsim = mini(fsim, con->nsamp-fsim);
	if (!cons) cons = (double) log(2)+2*lnfac(con->nsamp)-lnfac(2*con->nsamp)-lnfac(fsim-1)-lnfac(con->nsamp-fsim-1);

	flist = dvector(1,tree_size);
	/*Weight branches by marginal Pr(j | i)*/
	for (i=1;i<tree_size;i++) {
		maf = mini(tree_site[i]->ndesc, con->nsamp-tree_site[i]->ndesc);
		part[0] = (double) lnfac(con->nsamp+fsim-maf-1)+lnfac(con->nsamp-fsim+maf-1);
		part[1] = (double) lnfac(fsim+maf-1)+lnfac(2*con->nsamp-fsim-maf-1);
		part[2] = (double) lnfac(maf)+lnfac(con->nsamp-maf);
		flist[i]=(double) tree_site[i]->time_above*(exp(cons+part[0]-part[2])+exp(cons+part[1]-part[2]));
		if (flist[i]>mx) mx=flist[i];
	}

	/*Normalise to workable numbers*/
	for (i=1,cump=0.0;i<tree_size;i++) {
		flist[i]/=(double) mx;
		cump += (double) flist[i];
	}
	r1 = (double) ran2()*cump;
	for (i=1,cump=0.0;i<tree_size && cump<r1;i++) cump+= (double) flist[i];
	i--;

	free_dvector(flist,1,tree_size);
	return tree_site[i];
}






/**************************************************/
/*Routine to mutate sequences at tips of genealogy*/
/**************************************************/

void seq_mut(nm, seqs, site, base) 
struct node_tree *nm;
int **seqs, base, site;
{

	if ((nm->d[0]==NULL) && (nm->d[1]==NULL)) { /*terminal*/
		seqs[nm->node_num][site] = base;
	}
	else  { 
		seq_mut(nm->d[0], seqs, site, base);
		seq_mut(nm->d[1], seqs, site, base);
	}

}


/******************************************************************/
/*Evolve sequences down genealogy with user-defined mutation model*/
/*NB: default mutation model is infinite sites.                   */
/******************************************************************/

void evolve(np, seqs, con, mm, mutmat, muts, site) 
struct node_tree *np;
int **seqs, *muts, site;
struct control *con;
double *mm, **mutmat;
{

	int nb, i, j, base;
	double t, dt, u;

	if ((np->d[0]==NULL) && (np->d[1]==NULL)) return; /*Note in current system only 2 or zero daughter nodes*/

	base = np->nuc;
	dt = (np->time) - ((np->d[0])->time);
	for (t=0, u=(con->theta[site])*mm[base]/2; t<dt; ) {
		t += (double) -log(ran2())/u;
		if (t<dt) {
			select_base(&nb, base, mutmat);
			base = nb;
			u = (con->theta[site])*mm[base]/2;
			seq_mut(np->d[0], seqs, site, base);
			(np->d[0])->nuc=base;
			(*muts)++;
		}
		else {/*printf(": No mut\n");*/ (np->d[0])->nuc=base;}
	}
	evolve(np->d[0], seqs, con, mm, mutmat, muts, site);

	base = np->nuc;
	dt = (np->time) - ((np->d[1])->time);
	for (t=0, u=(con->theta[site])*mm[base]/2; t<dt; ) {
		t += (double) -log(ran2())/u;
		if (t<dt) {
			select_base(&nb, base, mutmat);
			base = nb;
			u = (con->theta[site])*mm[base]/2;
			seq_mut(np->d[1], seqs, site, base);
			(np->d[1])->nuc=base;
			(*muts)++;
		}
		else {/*printf(": No mut\n");*/ (np->d[1])->nuc=base;}
	}
	evolve(np->d[1], seqs, con, mm, mutmat, muts, site);
}


/*****************************************************************/
/*Print sequences to stdout, segregating sites to file "seq" and */
/*locations of segregating sites to file "loc".                  */
/*****************************************************************/

void print_seqs(seqs, con) 
int **seqs;
struct control *con;
{
	int i, j, *seg, site, fi, ns;
	char c;
	char fname[128];
	FILE *seq, *loc;

	/*	
	printf("\n\nPrintout of sequences\n");
	for (i=1; i<=con->nsamp; i++) {
		printf("\nInd %3i: ", i);
		for (site=1; site<=con->len; site++) printf("%i", seqs[i][site]);
	}
	printf("\n\n");
	*/ 

	if (PHASE) {
		seg = ivector(1,con->len);
		strcpy(fname, con->prefix);
		seq = fopen(strcat(fname, "seq.phase"), "w");
		for (i=1;i<=con->len;i++) {
			for (fi=0,j=2; j<=con->nsamp; j++) if (seqs[j][i]!= seqs[1][i]) fi++;
			if (fi) {seg[i]=1;ns++;}
			else seg[i]=0;
		}
		for (i=1;i<=con->len;i++) if (seg[i]) fprintf(seq,"1 "); fprintf(seq,"\n");
		for (i=1;i<=con->len;i++) if (seg[i]) fprintf(seq,"%i ",i); fprintf(seq,"\n");
		for (j=1;j<=con->nsamp;j++) {
			for (i=1;i<=con->len;i++) if (seg[i])  fprintf(seq,"%i ",seqs[j][i]);
			fprintf(seq,"\n");
		}
		fclose(seq);
		exit(0);

	}

	if (con->gt && !(con->inf)) nrerror("Cannot have GT model without infinite sites!!");

	seg = ivector(1,con->len);
	for (i=1, ns=0; i<=con->len; i++) {
		for (fi=0,j=2; j<=con->nsamp; j++) if (seqs[j][i]!= seqs[1][i]) fi++;
		if (fi) {seg[i]=1; ns++;}
		else seg[i]=0;
	}
	strcpy(fname, con->prefix);
	seq=fopen(strcat(fname, "sim.seq"),"w");

	if (con->asc>0) fprintf(seq,"%i %i  %i",(int) con->nsamp-con->asc, ns,con->gt+1);
	else fprintf(seq,"%i %i  %i",con->nsamp,ns,con->gt+1);

	if (!con->gt) {
	  for (i=(int) (con->asc>0 ? con->asc+1 : 1); i<=con->nsamp;i++) {
		fprintf(seq,"\n>Seq%i\n",i);
		for (j=1,fi=0;j<=con->len;j++) {
			if (seg[j]) {
			  if (con->mut) {
				c = num_to_nuc(seqs[i][j]);
				fprintf(seq,"%c",c); 
			  }
			  else fprintf(seq,"%i",seqs[i][j]%2);
			  fi++;
			  if ((fi%50)==0) fprintf(seq,"\n");
			}
		}
	  }
	}

	else {
	  for (i=(int) (con->asc>0 ? con->asc+1 : 1); i<=con->nsamp;i+=2) {
	    fprintf(seq,"\n>GT%i\n",(int) i/2+1);
	    for (j=1,fi=0;j<=con->len;j++) {
	      if (seg[j]) {
		fprintf(seq,"%i",seqs[i][j]==seqs[i+1][j]?seqs[i][j]:2);
		fi++;
		if (!(fi%50)) fprintf(seq,"\n"); 
	      }
	    }
	  }
	}

	fclose(seq);
	strcpy(fname, con->prefix);
	loc=fopen(strcat(fname, "sim.loc"),"w");
	fprintf(loc,"%i %.3f L", ns, (double) con->rescale);
	for (i=1;i<=con->len;i++) 
		if (seg[i]) 
			fprintf(loc,"\n%.3f", ((double)i)/((double)con->len)*con->rescale);
	fclose(loc);
	free_ivector(seg,1,con->len);

}


/********************************************/
/*Print results to stdout: summaries of data*/
/********************************************/

void print_res(res, con) 
struct results res;
struct control con;
{

	int i;
	double an, em, es;

	for (i=2, an=1; i<con.nsamp; i++)
		an += (double) 1/i;
	for (i=1, em=0.0, es=0.0;i<=con.len; i++)
	{
		em+= (double) con.theta[i];
		es += (double) 1-exp(-con.theta[i]*an);
	}
	for (i=1;i<=con.nsamp;i++)
		res.fdist[i]/=res.cf;

	printf("\n\nResults of %i simulations\n\n", con.nrun);
	printf("Average no. mutations  = %.3f\n", (double) res.nm/res.cf);
	printf("Expected               = %.3f\n", (double) em*an);
	printf("Average Sn             = %.3f\n", (double) res.sn/res.cf);
	printf("Expected               = %.3f\n", (double) es);
	printf("Average Pairwise Diff. = %.3f\n", (double) res.pwd/res.cf);
	printf("Expected               = %.3f\n", (double) em);
	printf("Average no. rec events = %.3f\n", (double) res.nr/res.cf);
	printf("Expected               = %.3f\n", (double) con.R*an);
	if (con.hyp)
		printf("Proportion muts HM     = %.3f\n", (double) res.mhm/res.cf);
	printf("ESS                    = %.3f\n", (double) (res.cf)*(res.cf)/(res.cf2));
	printf("\n\n");

	printf("Frequency distribution\n\n");
	for (i=1;i<con.nsamp-((int) con.asc==2?NPANEL+1:0);i++)
		printf("%3i %8.3f\n",i,res.fdist[i]);
	printf("\n\n");
}


/*****************************************/
/*Convert numerical coding to nucleotides*/
/*****************************************/

char num_to_nuc(i)
int i;
{
	char c;
	switch (i) {
		case 0:
			c = 'N';
			break;
		case 1 : 
			c = 'T';
			break;
		case 2 : 
			c = 'C';
			break;
		case 3 : 
			c = 'A';
			break;
		case 4:
			c = 'G';
			break;
		default :
			printf("\n\nError in nuc_to_num\n\n");
			exit(1);
	}
	return c;
}


/***********************************************************/
/*Read input parameters from command line or by user prompt*/
/***********************************************************/


void read_input(con, argc, argv) 
struct control *con;
int argc;
char *argv[];
{
	int i, ss, tl, nmut, freq;	
	extern long int *idum;
	char c;
	FILE *ifp;
	char *in_str;

	con->a=0; con->hyp=0; con->mut=0, con->seed=0; con->print=0; con->growth=0; con->lambda=0.0;
	con->inf=0; con->fmin=0; con->cond=0; con->rm=0; con->asc=0;
	con->bneck=0; con->tb=con->strb=0.0;
	con->conv=0;con->clen=0; con->cratio=0;
	con->gt=0;
	con->nsamp = 0; con->len = 0; con->R = -1.0; con->nrun = 1;
	con->p_genotype_error = 0.0;
	con->p_bad_site = 0.0;
	con->p_switch_error = 0.0;
	strcpy(con->prefix, "fin-sim");

	print_help(argc, argv);

	for(i = 0; i < argc; i++)
	{
		if(*argv[i] == '-')
		{ 
			in_str = argv[i];
			if(strcmp(in_str, "-nsamp") == 0) con->nsamp = atoi(argv[i + 1]);		/* no. of samples */
			if(strcmp(in_str, "-len") == 0) 
			{
				con->len = atoi(argv[i + 1]);			/* Sequence Length */
			}
			if(strcmp(in_str, "-R") == 0) con->R = atof(argv[i + 1]);			/* max 4Ner */
			if(strcmp(in_str, "-nruns") == 0) con->nrun = atof(argv[i + 1]);			/* Number of runs */
			if(strcmp(in_str, "-prefix") == 0) strcpy(con->prefix, argv[i+1]);
			if(strcmp(in_str, "-genoErr") == 0) con->p_genotype_error = atof(argv[i+1]);	/* Genotype error probability */
			if(strcmp(in_str, "-siteErr") == 0) con->p_bad_site = atof(argv[i+1]);	/* Bad site probability */
			if(strcmp(in_str, "-switchErr") == 0) con->p_switch_error = atof(argv[i+1]);	/* Switch Error probability */
		}
	}

	if (con->nsamp == 0)
	{
		printf("\nSample size?\n");
		scanf("%i", &(con->nsamp));
	}

	if (con->len == 0.0)
	{
		printf("\nLength?\n");
		scanf("%i", &(con->len));
	}

	con->theta = dvector(1,con->len);
	con->rmap = dvector(1,con->len);
	con->theta[1] = -1.0;

	for(i = 0; i < argc; i++)
	{
		if(*argv[i] == '-')
		{ 
			in_str = argv[i];
			if(strcmp(in_str, "-theta") == 0) con->theta[1] = atof(argv[i+1]);		/* Theta */
		}
	}

	if (con->theta[1] < 0)
	{
		printf("\nTheta?\n");
		scanf("%lf", &(con->theta[1]));
	}
	con->rescale=(double)con->len;
	read_flags(con, argc, argv);

	for (i=1, con->w0=0.0;i<con->nsamp;i++) con->w0 += (double) 1/i;
	printf("\nWatterson's const for n=%i = %.3f\n", con->nsamp, con->w0);

	if (con->seed == 0)
		con->seed = -setseed();
	idum = &(con->seed);

	printf("Random Seed: %li\n", con->seed);

	for (i=2;i<=con->len;i++) con->theta[i]=con->theta[1];

	if (con->a>0) 
		for (i=2;i<=con->len;i++) con->theta[i] *= (double) rgamm(con->a);

	if (con->hyp) {
		con->hyp = (int) ((double) (con->len)*(con->phm));
		printf("\nTotal of %i hypermutable sites: ", con->hyp);
		for (i=1;i<=con->hyp;i++) {
			nmut = (int) ((double) ((con->len)-1)*ran2())+2; 
			con->theta[nmut]=con->theta[1]*(con->rhm); 
			printf(" %i",nmut);
		}
	}

	if ((con->a)||(con->hyp)) {
	  printf("\n\nThetas for sites\n\n");
	  for (i=1;i<=con->len;i++) printf("%4i : %.5f\n", i, con->theta[i]);
	}

	if (con->cond) {
	  con->slocs = imatrix(1,con->len,1,2);
	  ifp = fopen(con->flocs_file, "r");
	  for (i=1; i<=con->len; i++) con->slocs[i][1] = con->slocs[i][2]=0;
	  if (ifp) {
	    printf("\n\nReading seg sites from file: %s\n\n", con->flocs_file);
	    fscanf(ifp,"%i %i %c", &ss, &tl, &c);
	    printf("\nCheck: %i %i %c \n\n", ss, tl, c);
	    if (tl!=con->len) {printf("\n\nSimulated and conditioning length not the same\n\n"); exit(1);}
	    for (i=1;i<=ss;i++) {
	      if (feof(ifp)) {printf("\n\nError in flocs file\n\n"); exit(1);}
	      fscanf(ifp,"%i %i", &tl, &freq);
	      con->slocs[tl][1]=1; con->slocs[tl][2]=freq;
	    }
	    fclose(ifp);
	  }
	  else {
	    printf("\n\nAssuming all sites are segregating\n\n");
	    for (i=1;i<=con->len;i++) {con->slocs[i][1]=1; con->slocs[i][2]=0;}
	  }
	  printf("\n\nConditioning on sites segregating: Freq(0=not conditioned)\n\n");
	  for(i=1;i<=con->len; i++) if (con->slocs[i][1]) printf("%5i %3i\n", i, con->slocs[i][2]);
	}

	if (con->rm)
	{
		ifp = fopen(con->rmap_file,"r");
		if (!ifp) {printf("\n\nCannot open file: %s\n\n", con->rmap_file); exit(1);}
		fscanf(ifp,"%i",&tl);
		if (tl!=con->len) {
			printf("\nError: recombination for incorrect sequence length\n");
			printf("tl=%i con->len=%i\n\n",tl,con->len);
			exit(1);}
		for (i=1;i<=con->len;i++) {
		  fscanf(ifp,"%i %lf", &ss, &(con->rmap[i]));
		  if ((i>1)&&(con->rmap[i]<con->rmap[i-1])) {printf("\n\nError in Map\n\n"); exit(1);}
		}
		if (DEBUG)
		{
			printf("\n\nRecombination map:\n\n");
			for (i=1;i<=con->len;i++) printf(" %4i %5.2f\n", i, con->rmap[i]);
		}
		con->R = con->rmap[con->len];
		printf("Recombination map length: %f\n", con->R);
	}
	else
	{
		if (con->R < 0)
		{
			printf("\nR?\n");
			scanf("%lf", &(con->R));
		}
		for (i=1;i<=con->len;i++)
			con->rmap[i]=(double) (i-1)*(con->R)/((con->len)-1);
	}

	printf("\n\nParameters: sample %i: Length %i: theta %.4f: R %.4f: Nrun %i : seed = %li\n", \
		con->nsamp, con->len, con->theta[1], con->R, con->nrun, con->seed);
	if (con->hyp) printf("HM: phm = %.3f rhm = %.3f\n", con->phm, con->rhm);

}


/******************************/
/*Read flags from command line*/
/******************************/

void read_flags(con, argc, argv) 
struct control *con;
int argc;
char *argv[];
{
	int i;
	char *in_str;
	
	for(i = 0; i < argc; i++)
	{
		if(*argv[i] == '-')
		{ 
			in_str = argv[i];
			if ((strcmp(in_str, "-s") == 0) || (strcmp(in_str, "-seed") == 0))
			{	/*User defined seed*/
				if ((i+1)<argc) con->seed = atol(argv[i+1]);
				else {printf("\nIncorrect flag use: -s seed\n"); exit(1);}
				if (con->seed > 0)
					con->seed = -con->seed;
			}
			if(strcmp(in_str, "-v") == 0)
			{	/*Include hypermutable sites: phm=Proportion, rhm=Relative rate */
				con->hyp = 1;
				if ((i+2)<argc) {
					con->phm = atof(argv[i+1]);
					con->rhm = atof(argv[i+2]);
				}
				else {printf("\n\nIncorrect flag use: -h pHM rHM\n"); exit(1);}
				if ((con->phm<0.0)||(con->phm>1.0)) {printf("\n0<=pHM<=1!\n"); exit(1);}
				if (con->rhm<0.0) {printf("\nrHM>0!\n"); exit(1);}
				printf("\nIncluding hypermutable sites: pHM = %.4f, rHM = %.1f\n",con->phm, con->rhm);
			}
			if(strcmp(in_str, "-m") == 0)
			{	/*User-defined mutation model*/
				printf("\nIncluding non-standard mutation model: in file mut.in\n");
				con->mut=1;
				strcpy(con->mut_file, argv[i+1]); 
			}
			if(strcmp(in_str, "-a") == 0)
			{	/*Gamma distribution of mutation rates - mean = a*/
				if ((i+1)<argc) con->a = atof(argv[i+1]);
				else {printf("\nIncorrect flag usage: -a alpha value\n"); exit(1);}
				if (con->a>=0) printf("\nRate heterogeneity: alpha = %.2f\n", con->a);
				else {printf("\na>=0!\n"); exit(1);}
			}
			if(strcmp(in_str, "-p") == 0)
			{	/*Print sequences to file: sites and locs of seg sites to file: locs*/
				printf("\nPrinting sequences (segregating sites only)\n");
				con->print = 1;
			}
			if(strcmp(in_str, "-g") == 0)
			{	/*Population growth*/
				printf("\nModel with population growth: ");
				con->growth=1;
				if ((i+1)<argc) con->lambda = atof(argv[i+1]);
				else {printf("\n\nMust enter growth rate for population\n\n"); exit(1);}
				printf(" rate = %.3f\n\n", con->lambda);
			}
			if(strcmp(in_str, "-b") == 0)
			{	/*Population bottleneck*/
				printf("\nModel with population bottleneck: ");
				con->bneck=1;
				if ((i+2)<argc) {con->tb = atof(argv[i+1]); con->strb = atof(argv[i+2]);}
				else {printf("\n\nMust enter time and strength of bottleneck\n\n"); exit(1);}
				printf(" time = %.3f  strength = %.3f\n\n", con->tb, con->strb);
			}
			if(strcmp(in_str, "-f") == 0)
			{	/*Min freq for mutations*/
                if ((i+1)<argc) con->fmin = atoi(argv[i+1]);
                if ((con->fmin>0)&&(con->fmin)<(int)con->nsamp/2) 
                	printf("\n\nRestricting analysis to derived mutations represented > %i in sample\n\n", con->fmin);
                else {printf("\n\nInappropriate cut-off frequency\n\n"); exit(1);}
			}
			if(strcmp(in_str, "-i") == 0)
			{	/*Infinite-sites model*/
				printf("\nUsing infinite site model\n");
				con->inf = 1;
			}
			if(strcmp(in_str, "-c") == 0)
			{	/*Condition on mutations occuring at all sites*/
                printf("\nConditioning on mutations occuring at sites\n");
                con->cond=1;
                strcpy(con->flocs_file, argv[i+1]);
			}
			if(strcmp(in_str, "-r") == 0)
			{	/*User defined recombination map*/				
				printf("\n\nUser defined recombination map: read from file\n\n");
				strcpy(con->rmap_file, argv[i+1]);
				con->rm=1;
			}
			if(strcmp(in_str, "-n") == 0)
			{	/*SNP ascertainment strategy - simulate an additional n chromosomes*/
				con->asc = atoi(argv[i+1]);
				if (con->asc<2) nrerror("*** Error: ascertainment sample must be at least 2 chromosomes ***");	      
			}
			if(strcmp(in_str, "-x") == 0)
			{	/*Gene conversion model - input ratio of gene conversion to recombination and average tract length*/
			      con->conv=1;
				if ((i+2)<argc) {
				  con->cratio = atof(argv[i+1]);
				  con->clen = atof(argv[i+2]);
				}
				else nrerror("Insufficient parameters for conversion model");
				printf("\n\nGene conversion model: ratio = %.3f: geometric distribution of tract lengths with mean = %.3f\n\n",con->cratio, con->clen);
			}
			if(strcmp(in_str, "-d") == 0)
			{	/*Print out sequences as genotypes in 0=00, 1=11, 2=01 format*/
			    con->gt=1;
			}
			if(strcmp(in_str, "-L") == 0)
			{	/* Rescale output loci file so that loci between 0 and L */
				printf("\nRescaling loci in output locs file\n");
				con->rescale = atof(argv[i+1]);
			}
		}
	}
}


/**********************************************************/
/*Choose mutation according to user-defined mutation model*/
/**********************************************************/

void select_base(nb, base, mut_mat) 
int *nb, base;
double **mut_mat;
{
	int i;
	double cump=0, r1=ran2();

	for (i=1; cump<r1; i++) cump+=mut_mat[base][i];
	(*nb)=i-1;
}



/********************************************************/
/*Routine to find number of descendants for a given node*/
/********************************************************/

int count_desc(node)
struct node_tree *node;
{
	int ndesc=0;
	if (node->d[0] == NULL) ndesc=1;
	else {
		ndesc += count_desc(node->d[0]);
		ndesc += count_desc(node->d[1]);
	}

	return ndesc;
}





/**************************************************/
/*Routine to choose time for next coalescent event*/
/**************************************************/

void choose_time(t, k, rho, con) 
int k;
double *t, rho;
struct control *con;
{
  double u=-log(ran2()), cons[3];

  if (!(con->growth)) { /*Constant population size*/
    (*t) += (double) 2*u/(k*rho+k*(k-1));
  }
  else { /*Exponential growth*/
    if (con->lambda == 0) {printf("\n\nCannot have zero growth\n\n"); exit(1);}
    cons[0]=(double) rho;
    cons[1] = (double) (k-1)/(con->lambda)*exp((*t)*(con->lambda));
    cons[2]= (double) (con->lambda);
    (*t) += (double) bisect(tgrowth, (double) 2*u/k, cons);
  }
}

/***************************************************/
/*Routine to solve monotonic functions by bisection*/
/***************************************************/

double bisect(bfunc, val, cons) 
double val, *cons, (*bfunc)(double *, double **);
{

  double x[3], y[3];

  x[0]=-1.0; x[2]=1.0;
  while ((y[0]=bfunc(&x[0],&cons))>val) x[0]*=2.0;
  while ((y[2]=bfunc(&x[2],&cons))<val) x[2]*=2.0;

  while (fabs(x[0]-x[2])>BACC) {
    x[1]=(double) (x[0]+x[2])/2; 
	y[1]=bfunc(&x[1],&cons);
    if (y[1]<val) {x[0]=x[1]; y[0]=y[1];}
    else {x[2]=x[1]; y[2]=y[1];}
  }

  return x[1];
}

/**************************/
/*Population growth module*/
/**************************/

double tgrowth(var, cons) 
double *var, **cons;
{
	return (double) (*cons)[0]*(*var)+(*cons)[1]*(exp((*cons)[2]*(*var))-1);;
}


void check_lin(list, k, con)
int *k;
struct node_list **list;
struct control *con;
{
	int i;

	for (i=1;i<=*k;i++) if (list[i]->rlen<0) {
		printf("\n\nError in rlen at %i (%f)!!",i, (double) list[i]->rlen);
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
				printf("\nfin\n");
				printf("Simulate LDhat formatted data using Coalescent with recombination\n\n");
				
				printf("Required Options :\n");
				printf("-nsamp <int>       Number of sampled chromosomes\n");
				printf("-len <int>         Sequence Length (sites)\n");
				printf("-theta <float>     Mutation Rate (4Nu) per site\n");
				printf("-R <float>         Total Recombination Map Length (4Nr)\n");
				printf("\n\n");	
				printf("Additional Options :\n");
				
				printf("-i                 Infinite-sites model for sequences\n");
				printf("-p                 Prints segregating sites and locations to files \n");
				printf("-c <file>          Condition on mutations at sites\n");
				printf("                     (if file not found, condition on all sites)\n");
				printf("-f <int>           Restrict Analysis to mutations segregating at minimum\n");
				printf("                     frequency (specified: e.g. -f 2 to remove singletons)\n");
				printf("-g <float>         Population growth with exponential growth parameter\n");
				printf("                     (e.g. -g 0.1)\n");
				printf("-b <float> <float> Population bottleneck time and strength (e.g. -b 0.1 10)\n");
				printf("-v <float> <int>   Hypermutable sites model: proportion sites : ratio rates\n");
				printf("                     (e.g. -v 0.01 100)\n");
				printf("-a <float>         Gamma distributed rates across sites : shape parameter\n");
				printf("                      (e.g. -a 1.0)\n");
				printf("-m <file>          User-defined mutation model file\n");
				printf("-r <file>          User-specifed recombination map file\n");
				printf("-s <int>           User-defined seed (e.g. -s 193429012)\n");
				printf("-n 2               n=2 SNP ascertainment model\n");
				printf("-x <float> <float> Gene conversion model: input ratio conversion to x-over\n");
				printf("                      and average tract length (e.g. -x 2 100)\n");
				printf("-d                 Print out sequences as genotypes\n");
				printf("-L <float>         Rescale output loci file so that loci between 0 and L\n");
				printf("                      (e.g. -L 6.5)\n");
				printf("-nruns <int>       Number of simulations (will overwrite seq and loc files)\n");
				printf("-prefix <string>   Prefix of output files\n");
				printf("-genoErr <float>   Probability of a genotype error\n");
				printf("-siteErr <float>   Probability of bad site\n");
				printf("-switchErr <float> Probability of switch error\n");
				
				printf("\n\n");	
				exit(0);
			}
		}
	}
}

