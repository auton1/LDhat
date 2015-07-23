#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>


#include "tools.h"

#include "seqtools.h"
#include "snp_sim.h"
#include "ldhat.h"

int n_node=0;
int tree_size;

#define SDEBUG 0
#define FSAMP 10


void snp_sim(locs,flocs,pset,lkmat,nrun,data)
int *flocs, nrun;
double *locs;
double **lkmat;
struct site_type **pset;
struct data_sum *data;
{
	int nmut, nrec, run, i, **segl, n5;
	double rtot, stot, *rho_sim, *clr_sim, *vmin, *vmax;
	extern int n_node;
	extern int tree_size;
	struct node **tree_ptr, *new_node;
	struct control con;
	char fname[MAXNAME+1];
	FILE *ofp, *sfp;

	strcpy(con.prefix, data->prefix);
	printf("\n\nCoalescent simulations with mutations and recombination: After Hudson (1991)\n\n");
	strcpy(fname, data->prefix);
	ofp = fopen(strcat(fname, "sim_out.txt"), "w");
	if (!ofp) nrerror("Cannot open outfile for simulations");
	fprintf(ofp,"\nRho     Fit     CLR     CF\n\n");

	  con.nsamp = data->nseq;
	  con.len = data->lseq;
	  con.rm = 1;
	  con.rmap = dvector(1,con.len);
	  for (i=1;i<=data->lseq;i++) con.rmap[i] = (double) (locs[i]-1)*(data->rho)/(data->tlseq-1);
	  con.cond=1;
	  con.growth=con.print=con.fmin=con.bneck=con.lambda=0;
	  con.slocs = imatrix(1,data->lseq,1,2);
	  for (i=1;i<=data->lseq;i++) {con.slocs[i][1]=1; con.slocs[i][2] = flocs[i];}
	  data->rho_drive = data->rho;
	  data->fit_obs = data->fit;
	  con.cl_fac=0.0;
	  printf("\nObserved fit = %.3f (rho = %.3f)\n", data->fit_obs, data->rho_drive);

	rho_sim = dvector(1,nrun);
	clr_sim = dvector(1,nrun);


	tree_size = 2*(con.nsamp);
	tree_ptr = (struct node **) malloc((size_t) tree_size*sizeof(struct node *));
	segl = imatrix(1,2,1,con.len);
	data->ng[0]=data->ng[1]=0;
	for (i=1; i<tree_size; i++) {
		new_node = (struct node *) malloc((size_t) sizeof(struct node));
		new_node->asite = (int *) malloc((size_t) ((con.len)+1)*sizeof(int));
		new_node->d[0]=NULL; new_node->d[1]=NULL;
		new_node->a[0]=NULL; new_node->a[1]=NULL;
		tree_ptr[i] = new_node;
	}

	for (run=0, rtot=0, stot=0; run<nrun; run++) {
		if ((run%1000)==0) printf("\nRun %i", run);
		nmut=nrec=n_node=0;
		tree_ptr = make_tree(tree_ptr, &con, &nrec, segl);
		tree_summary(tree_ptr, n_node, &con, segl, data, pset, lkmat, locs, ofp);
		rho_sim[run+1]=data->rho;
		clr_sim[run+1]=data->clr;
	}
	
	for (i=1; i<tree_size; i++) {
		free(tree_ptr[i]->asite);
		free(tree_ptr[i]);
	}
	free(tree_ptr);
	free_imatrix(segl,1,2,1,con.len);
	free_dvector(con.rmap,1,con.len);
	free_imatrix(con.slocs,1,data->lseq,1,2);

	fprintf(ofp, "\n\nNumber of simulated fits > obs = %i (of %i)\n", data->ng[0], nrun);
	fprintf(ofp, "Number of simulated rho  > obs = %i\n\n", data->ng[1]);
	fclose(ofp);
	
	printf("\n\nNumber of simulated fits > obs = %i (of %i)\n", data->ng[0], nrun);
	printf("\n\nRho = %.3f: ", data->rho_drive);

	n5 = (int) nrun/40+1;
	vmin = dvector(1,n5);
	vmax = dvector(1,n5);
	if (nrun<5) nrerror("\n\nToo few runs to calculate 95%CIs\n\n");

	for (run=1;run<=n5;run++) vmin[run]=vmax[run]=rho_sim[run];
	sort(vmin,n5); sort(vmax,n5);
	for (run=n5+1;run<=nrun;run++) {
	  if (rho_sim[run]<vmin[n5]) {vmin[n5]=rho_sim[run]; sort(vmin,n5);}
	  else if (rho_sim[run]>vmax[1]) {vmax[1]=rho_sim[run]; sort(vmax,n5);}
	}	
	printf("\nPBootstrap Rho CI: %.2f - %.2f\n\n",vmin[n5],vmax[1]);

	free_dvector(rho_sim,1,nrun);
	free_dvector(clr_sim, 1, nrun);
	free_dvector(vmin,1,n5);
	free_dvector(vmax,1,n5);
}




/****************************************/
/*The key routine: creates the genealogy*/
/****************************************/

struct node **  make_tree(tree_ptr, con, revts, segl) 
struct node **tree_ptr;
struct control *con;
int *revts, **segl;
{

	int k, i, j, nrec=0, nco=0, ncb;
	extern int n_node;
	extern int tree_size;
	double rr, t=0, avR, told, cump;
	struct node **list;

	list = (struct node **) malloc((size_t)((con->nsamp+1)*sizeof(struct node *)));
	for (i=1; i<=con->len; i++) segl[1][i] = 1;

	for (i=1; i<=con->nsamp; i++) {
		tree_ptr[i]->time = (double) 0;
		tree_ptr[i]->node_num = i;
		tree_ptr[i]->d[0] = tree_ptr[i]->d[1] = NULL;
		tree_ptr[i]->a[0] = tree_ptr[i]->a[1] = NULL;
		for (j=1; j<=con->len; j++) tree_ptr[i]->asite[j] = 1;
		tree_ptr[i]->rlen=(double) con->rmap[con->len]-con->rmap[1];
		list[i] = tree_ptr[i];
	}
	for (i=(con->nsamp)+1; i<tree_size; i++) {
		tree_ptr[i]->d[0] = tree_ptr[i]->d[1] = NULL;
		tree_ptr[i]->a[0] = tree_ptr[i]->a[1] = NULL;
	}

	if (SDEBUG) {
	  printf("\n\nInitial state of sample\n\n");
	  print_lin(list, con->len, con->nsamp);
	}

	k=con->nsamp;	
	n_node = con->nsamp;
	
	while (k>1) {
	
	/*Count lineages which can contribute to rec events*/

		for(i=1, rr=0; i<=k; i++) {
			rr += list[i]->rlen;
		}			
		avR = (double) rr/k;
		if (SDEBUG) printf("\nTotal recombination length = %.0f , Re = %.3f", rr, avR);

	/*Choose time for next event*/
		
		told=t;
		choose_time(&t, k, avR, con);
		if (n_node+2>= tree_size) tree_ptr = add_tree(tree_ptr,  n_node-2,  con->len);
		if (SDEBUG>2) printf("\nNext event at time %.3f\n", t);

	/*Choose type of next event: update list and tree accordingly*/

		if ((con->bneck)&&(told<con->tb)&&(t>con->tb)) { /*Bottleneck*/
			t=con->tb;
			told = con->tb+0.0001;
			for(cump=0,ncb=0;cump<con->strb;) {
				cump+=(double) -log(ran2())*2/(k*(k-1));
				if (cump<con->strb) {
					k--;nco++;n_node+=1; ncb++;
					coalesce(&k,con,list,tree_ptr,segl,&t);
					list = (struct node **) realloc(list, (size_t) (k+1)*sizeof(struct node *));
					if (n_node+2>= tree_size) tree_ptr = add_tree(tree_ptr,  n_node-2,  con->len);
				}
			}
			if (SDEBUG) printf("\n\nA total of %i coalescent events occurred during the bottlneck leaving %i\n\n",ncb,k);
		}


		else if (ran2() < (double) avR/(avR+(k-1)*exp((con->lambda)*t))) {  /*Event is a recombination*/
			k++; nrec++; n_node+=2; 
			list = (struct node **) realloc(list, (size_t) (k+1)*sizeof(struct node *));
			recombine(&k, &rr, con, list, tree_ptr, &t);
		}

		else { /*Event is a coalescent*/
			nco++; k--; n_node += 1;
			coalesce(&k,con,list,tree_ptr,segl,&t);
			list = (struct node **) realloc(list, (size_t) (k+1)*sizeof(struct node *));
		}

		if (SDEBUG >1) print_lin(list, con->len, k);
	}

	if (SDEBUG) printf("\nEnd of run: %.3f Nrec: %i Nco: %i\n", t, nrec, nco);
	if (SDEBUG) printf("\n\nMRCAs for each site\n\nSite  MRCA\n\n");
	for (i=1; i<=con->len; i++) {
		if (SDEBUG) printf("%4i  %4i\n", i, segl[2][i]); 
		tree_ptr[segl[2][i]]->asite[i]=1;
	}
	if (SDEBUG >1)	print_nodes(tree_ptr, n_node, con->len);
	if (n_node != (con->nsamp)+nco+2*nrec) {printf("\nN_Node does not tally !\n\n"); exit(1);}
	free(list);
	(*revts) = nrec;

	return tree_ptr;
}



void coalesce(k,con,list,tree_ptr,segl,t)
int *k, **segl;
double *t;
struct control *con;
struct node **list, **tree_ptr;
{
	int i, j, fs;

	extern int n_node;
	extern int tree_size;

			i=(int) ((double) 1+((*k)+1)*ran2()); 
			if (DEBUG ==3) printf(" %i", i);
			tree_ptr[n_node]->node_num = n_node;
			tree_ptr[n_node]->d[0] = list[i];
			for (j=1; j<=con->len; j++) tree_ptr[n_node]->asite[j]=list[i]->asite[j];
			tree_ptr[n_node]->time = (*t);
			list[i]->a[0] = tree_ptr[n_node];
			list[i]=list[(*k)+1];
			i=(int) ((double) 1+(*k)*ran2());
			if (SDEBUG ==3) printf(" %i", i);
			tree_ptr[n_node]->d[1] = list[i];
			for (j=1; j<=con->len; j++) if (list[i]->asite[j]) tree_ptr[n_node]->asite[j] = 1;
			count_rlen(tree_ptr[n_node]->asite, con->len, &(tree_ptr[n_node]->rlen),con);
			list[i]->a[0] = tree_ptr[n_node];
			list[i] = tree_ptr[n_node];

		/*Remove sites -> lineages for which MRCA reached*/

			for (i=1; i<=con->len; i++) {
				if (segl[1][i]) {
					for (j=1, fs=0; j<=(*k); j++) if (list[j]->asite[i]) fs++;
					if (fs < 2) {
						if (SDEBUG > 1) printf("\n\nMRCA reached for site %i\n\n", i); 
						segl[1][i] = 0;
						for (j=1; j<=(*k); j++) list[j]->asite[i]=0;
						segl[2][i] = n_node;
						for (j=1; j<=(*k); j++) count_rlen(list[j]->asite, con->len, &(list[j]->rlen),con);
					}
				}
			}
			for (i=1; i<=(*k); i++) {
				for (j=1, fs=0; j<=con->len; j++) if (segl[1][j]) fs += list[i]->asite[j];
				if (fs == 0) {
					if (SDEBUG >1) printf("\n\nRemoving lineage %i\n\n", i);
					list[i] = list[*k];
					(*k)--;
				}
			}
}


/***********************/
/*Genereate recombinant*/
/***********************/

void recombine(k, rr, con, list, tree_ptr, t) 
int *k;
double *t, *rr;
struct control *con;
struct node **list, **tree_ptr;
{


	extern int n_node;
	extern int tree_size;

	int i,j,j0;
	double cump, ri;
	
		/*Choose lineage to recombine*/

			ri = ran2()*(*rr);
			for (i=1, cump=0; cump<ri; i++) cump += list[i]->rlen;
			i--;
			if (SDEBUG > 1) printf("\nLineage %i recombines at position ", i);
			tree_ptr[n_node-1]->node_num = n_node-1;
			tree_ptr[n_node]->node_num   = n_node;
			tree_ptr[n_node-1]->d[0] = list[i];
			tree_ptr[n_node]->d[0] =   list[i];
			list[i]->a[0] = tree_ptr[n_node-1];
			list[i]->a[1] = tree_ptr[n_node];

			ri = (double) ran2()*(list[i]->rlen);
			for (j0=1;list[i]->asite[j0]<1;j0++);
			for (cump=0,j=j0;cump<ri;j++) cump=con->rmap[j]-con->rmap[j0];
			j0=j-2;
			if (SDEBUG > 1) printf("%i\n", j0);
			for (j=1; j<=j0; j++) {tree_ptr[n_node-1]->asite[j]=list[i]->asite[j]; tree_ptr[n_node]->asite[j]=0;}
			for (;j<=con->len; j++) {tree_ptr[n_node-1]->asite[j]=0; tree_ptr[n_node]->asite[j] = list[i]->asite[j];}
			count_rlen(tree_ptr[n_node-1]->asite, con->len, &(tree_ptr[n_node-1]->rlen),con);
			count_rlen(tree_ptr[n_node]->asite, con->len, &(tree_ptr[n_node]->rlen),con);			
	
			list[i] = tree_ptr[n_node-1];
			list[*k] = tree_ptr[n_node];
			list[*k]->time=(*t); list[i]->time=(*t);
}



/*********************************************************/
/*Prints a list of lineages still active in the genealogy*/
/*********************************************************/

void print_lin(list, len, k) 
struct node **list;
int len, k;
{

	int i, j;

	printf("\n\nLineages active (%i)\n\n", k);
	for (i=1; i<=k; i++) {
		printf("lineage %3i: node %3i: time %.3f: rlen %6.2f :ANC ",\
			i, list[i]->node_num, list[i]->time, list[i]->rlen);
		for (j=1; j<=len; j++) printf("%i", list[i]->asite[j]);
		printf("\n");
	}
	printf("\n\n");
}


/******************************************************/
/*Print nodes in tree: useful for checking simulations*/
/******************************************************/

void print_nodes(tree, nn, len) 
struct node **tree;
int nn, len;
{

	int n, i;

	printf("\n\nPrinting nodes in tree\n\n");

	for (n=1; n<=nn; n++) {
		printf("node %3i: time %.3f  ANC ", tree[n]->node_num, tree[n]->time);
		for (i=1; i<=len; i++) printf("%i", tree[n]->asite[i]);
		if (tree[n]->d[0] != NULL) printf(" : DESC %i",(tree[n]->d[0])->node_num);
		if (tree[n]->d[1] != NULL) printf(" %i", (tree[n]->d[1])->node_num);
		if (tree[n]->a[0] != NULL) printf(" : ANC %i", (tree[n]->a[0])->node_num);
		if (tree[n]->a[1] != NULL) printf(" %i", (tree[n]->a[1])->node_num);
		printf("\n");
	}
	printf("\n\n");
}


/***************************************/
/*Calculate potential for recombination*/
/***************************************/

void count_rlen(asite, len, rlen, con) 
int *asite, len;
double *rlen;
struct control *con;
{

	int i;

	(*rlen) = 0;
	for (i=1; i<=len; i++) if (asite[i]) break;
	if (i>len) {(*rlen) = 0; return;}
	(*rlen) -= (double) con->rmap[i];
	for (i=len; i>0; i--) if (asite[i]) break;
	(*rlen) += (double) con->rmap[i];
}


/******************************************************************/
/*Add mutations to tree and calculate summary statistics of sample*/
/******************************************************************/

void tree_summary(tree, nn, con, segl, data, pset, lkmat, locs, ofp) 
struct node **tree;
int nn, **segl, *locs;
double **lkmat;
struct control *con;
struct site_type **pset;
struct data_sum *data;
FILE *ofp;
{

	int site, mrca, nmuts, i, j;
	int **seqs, node_mut, fl=1, pnew=0, miss=0, **pij, **nall;
	char bases[5]="NTCAG";
	extern int n_node;
	double tree_len=0, pwd=0, tt;
	double cf;

	seqs = imatrix(1, con->nsamp, 1, con->len);
	for (i=1; i<=con->nsamp; i++) for (j=1; j<=con->len; j++) seqs[i][j]=2;

/*	Add mutations to tree	*/

	for (site=1, cf=0.0; site<=con->len; site++) {
		mrca = segl[2][site];
		if (SDEBUG>1) printf("\nSite %3i : mrca = %3i", site, mrca);
		tree_len = tree_time(tree[mrca], site);
		if (con->cond) {
			        nmuts=con->slocs[site][1];
				if (!nmuts) exit(0);
				else {
				  if (!(con->slocs[site][2])) {
					cf+=(double) log(tree_len/(2*con->w0))-(tree_len/2-(con->w0))*(con->theta[site]);
					tt = (double) tree_len*ran2();
					node_mut = add_mut(tree[mrca], site, &tt);
					if (SDEBUG>1) printf(" : added to node %i\n", node_mut);
				  }
				  else node_mut = add_mut_f(con->slocs[site][2], nn, con->nsamp, &cf, tree, site, data->th);
				  if (SDEBUG>1) printf(" : mutation added to node %i\n", node_mut);
				  if (node_mut>0) seq_mut(tree[node_mut], seqs, con->nsamp, site, 3);
				}
		}
		else exit(0);
      	}
	
/*	print_seqs(seqs, con);*/

	pij = imatrix(1,con->len, 1, con->len);
	nall = imatrix(1,con->len,1,6);
	allele_count(seqs,con->nsamp,con->len,nall,0,1);
	pair_spectrum(seqs,data,nall,pset,&(data->ptt),&pnew,&miss,0,pij);


/*Would like to do this in iterative fashion - need C++ probably for memory handling*/	

	if (pnew+miss) nrerror("Should be no new or missing data in simulated");
	lk_surf(pset,pij,data,lkmat,data->th,locs,0);
	fit_pwlk(data,pij,locs,lkmat,0);
	free_imatrix(pij,1,con->len,1,con->len);
	free_imatrix(nall,1,con->len,1,6);
	if (data->fit > data->fit_obs) data->ng[0]++;
	if (data->rho > data->rho_drive) data->ng[1]++;
	printf("rho=%8.2f fit=%8.2f clr=%8.3f cf=%8.3f\n", data->rho, data->fit, data->clr, cf);
	fprintf(ofp,"%8.2f %8.2f %8.3f %8.3f\n", data->rho, data->fit, data->clr, cf);
	con->cl_fac += (double) data->clr;

	free_imatrix(seqs, 1,con->nsamp,1,con->len);
}



/*********************************/
/*Routine to double nodes in tree*/
/*********************************/

struct node ** add_tree(tree, n_node, len) 
struct node **tree;
int n_node, len;
{

	int i;
	struct node *new_node;
	struct node **new_tree;
	extern int tree_size;
	
	printf("\n.....Increasing size of memory for tree.....\n");

	new_tree = (struct node **) malloc((size_t) 2*tree_size*sizeof(struct node *));
	for (i=1; i<tree_size; i++) new_tree[i] = tree[i];
	for (i=tree_size; i<2*tree_size; i++) {
		new_node = (struct node *) malloc((size_t) sizeof(struct node));
		new_node->asite = (int *) malloc((size_t) (len+1)*sizeof(int));
		new_node->d[0]=NULL; new_node->d[1]=NULL;
		new_node->a[0]=NULL; new_node->a[1]=NULL;
		new_tree[i] = new_node;
	}
	free(tree);
	tree = new_tree;
	tree_size  *=2;
	return tree;
}


/*************************************************/
/*Routine to find total tree length for each site*/
/*************************************************/

double tree_time(node, site) 
struct node *node;
int site;
{

	double ttime=0;
	
	if ((node->d[0] == NULL) && (node->d[1] == NULL)) ttime = node->time;
	else if (node -> asite[site]) {
		if ((node->d[0] != NULL) && ((node->d[0])->asite[site])) {
			ttime+=(node->time) - ((node->d[0])->time);
			ttime+=tree_time(node->d[0], site);
		}
		if ((node->d[1] != NULL) && ((node->d[1])->asite[site])) {
			ttime+=node->time - (node->d[1])->time;
			ttime+=tree_time(node->d[1], site);
		}
	}
	else ttime+=0;
	return ttime;
}

/*********************************************/
/*Routine to place mutations on the genealogy*/
/*********************************************/

int add_mut(node, site, fl) 
struct node *node;
double *fl;
int site;
{

	static int nn=0;

	/*	printf("\nNode = %i, cumt=%.4f, nn=%i", node->node_num, *fl, nn);*/

	if (!(node->asite[site])) {printf("\n\nError in add_mut\n\n"); exit(1);}
	
	if (node->d[1]) {
	  if ((node->d[1])->asite[site]) {
	    if ((*fl)<0) return nn;
	      (*fl) -= node->time-(node->d[1])->time;
	      if ((*fl) < 0.0) {
		/*		printf(": mutation occurs on node %i ", (node->d[1])->node_num);*/
		nn = (node->d[1])->node_num;
	      }
	      else nn = add_mut(node->d[1], site, fl);
	  }
	}
	if (node->d[0]) {
	  if ((node->d[0])->asite[site]) {
	     if ((*fl)<0) return nn;
	      (*fl) -= node->time-(node->d[0])->time;
	      if ((*fl)<0.0) {
		/*		printf(": mutation occurs on node %i ", (node->d[0])->node_num);*/
		nn = (node->d[0])->node_num;
	      }
	      else nn = add_mut(node->d[0], site, fl);
	  }
	}	  

	return nn;

}

/******************************************************************************/
/*Routine to place mutations on the genealogy conditioning on allele frequency*/
/******************************************************************************/

int add_mut_f(fsim, nnode, nsamp, cf, tree, site, theta) 
struct node **tree;
double theta;
double *cf;
int nsamp, site, fsim, nnode;
{

	int i, desc, j=0, fl;
	double r1, len=0.0, ct=0, **lnode;

	lnode = dmatrix(1,nnode,1,3);
	for (i=1;i<=nnode;i++) lnode[i][1]=lnode[i][2]=lnode[i][3]=0;

	if (SDEBUG) {printf("\nSite %3i, fsim = %3i\n", site, fsim);}

	/*Create list of numbers of descendants at each site*/ 

	for (i=1;i<=nnode; i++) {
	  if (tree[i]->asite[site]) {
	    desc = count_desc(tree[i], site);
	    if (SDEBUG>1) printf("Node %3i (%i): #desc = %3i:", tree[i]->node_num, i, desc);
	    if (tree[i]->a[0]) { /*Either coalecsent or recombinant*/ 
	      if (SDEBUG>1) printf(" a0:");
	      if ((tree[i]->a[0])->asite[site]) { /*Add if ancestral*/ 
		lnode[i][1]=(double) tree[i]->node_num; 
		lnode[i][2]=((tree[i]->a[0])->time-tree[i]->time);
		lnode[i][3]=(double) desc;
		if (SDEBUG>1) printf("a");
	      }
	      else if (tree[i]->a[1]) { /*Recombinant*/ 
		if (SDEBUG>1) printf(" n a1:");
		if ((tree[i]->a[1])->asite[site]) { /*Add time if ancestral*/ 
		  lnode[i][1]=(double) tree[i]->node_num; 
		  lnode[i][2]=((tree[i]->a[1])->time-tree[i]->time);
		  lnode[i][3]=(double) desc;
		  if (SDEBUG>1) printf("a");
		}
		else { /*Is Site MRCA*/
		  lnode[i][1]=(double) tree[i]->node_num;
		  lnode[i][2]=0.0;
		  lnode[i][3]=0;
		  if (SDEBUG>1) printf("mrca");
		}
	      }
	    }
	    else { /*Is GMRCA*/
	      lnode[i][1]=(double) tree[i]->node_num;
	      lnode[i][2]=0.0;
	      lnode[i][3]=0;
	      if (SDEBUG>1) printf(" GMRCA ");
	    }
	    if (SDEBUG>1) printf(" : dtime = %.3f\n", lnode[i][2]);
	  }
	  else {  /*site not ancestral*/
	    lnode[i][1]= (double) tree[i]->node_num;
	    lnode[i][2]= 0.0;
	    lnode[i][3]=0;
	  }
	}



	/*Choose node to add mutation to: start with exact fsim and increase tolerance*/ 

	for (desc=0, fl=0; fl<1; desc++) {
	  for (i=1, len=0.0; i<=nnode; i++) 
	    if (lnode[i][3]) 
	      if ((fabs(lnode[i][3]-fsim)==desc)||(fabs(nsamp-lnode[i][3]-fsim)==desc)) len+=lnode[i][2];
	  if (len) {
	    r1 = ran2()*len;
	    for (i=1; r1>0.0; i++) 
	      if ((fabs(lnode[i][3]-fsim)==desc)||(fabs(nsamp-lnode[i][3]-fsim)==desc))
		r1 -= lnode[i][2];
	    i--;
	    fl=1;
	  }
	  else if (SDEBUG) printf("Cannot generate ndesc = %i or %i\n",fsim+desc, fsim-desc); 
	}
	if (SDEBUG) {printf(": node %.0f chosen with %.0f/%.0f  descs (fsim = %i)", lnode[i][1], lnode[i][3],nsamp-lnode[i][3], fsim);}
	
	(*cf) += (double) log((double) len*fsim*(nsamp-fsim)/(2*nsamp))-theta*((double) len/2-nsamp/(fsim*(nsamp-fsim)));

	j = (int) lnode[i][1];
	free_dmatrix(lnode,1,nnode,1,3); 
	return (int) j;

}






/**************************************************/
/*Routine to mutate sequences at tips of genealogy*/
/**************************************************/

void seq_mut(nm, seqs, nsamp, site, base) 
struct node *nm;
int **seqs, nsamp, site, base;
{

	if ((nm->d[0]==NULL) && (nm->d[1]==NULL)) {
		seqs[nm->node_num][site] = base;
	}
	else if (nm->asite[site]) {
		if ((nm->d[0]!=NULL)&&((nm->d[0])->asite[site])) 
			seq_mut(nm->d[0], seqs, nsamp, site, base);
		if ((nm->d[1]!=NULL)&&((nm->d[1])->asite[site])) 
			seq_mut(nm->d[1], seqs, nsamp, site, base);
	}
	else {printf("\n\nError: node %i not ancestral to site %i\n\n", nm->node_num, site); exit(1);}
}



/*****************************************************************/
/*Print sequences to stdout, segregating sites to file "seq" and */
/*locations of segregating sites to file "loc".                  */
/*****************************************************************/

void print_seqs(seqs, con) 
int **seqs;
struct control *con;
{
	char fname[MAXNAME+1];
	int i, j, *seg, site, fi, ns;
	FILE *seq, *loc;

	if (SDEBUG) {	
	printf("\n\nPrintout of sequences\n");
	for (i=1; i<=con->nsamp; i++) {
		printf("\nInd %3i: ", i);
		for (site=1; site<=con->len; site++) printf("%i", seqs[i][site]);
	}
	printf("\n\n");
	}

	seg = ivector(1,con->len);
	for (i=1, ns=0; i<=con->len; i++) {
		for (fi=0,j=2; j<=con->nsamp; j++) if (seqs[j][i]!= seqs[1][i]) fi++;
		if (fi) {seg[i]=1; ns++;}
		else seg[i]=0;
	}
	
	strcpy(fname, con->prefix);
	seq=fopen(strcat(fname, "seq"),"w");

	fprintf(seq,"%i %i",con->nsamp,ns);
	for (i=1;i<=con->nsamp;i++) {
		fprintf(seq,"\n>Seq%i\n",i);
		for (j=1,fi=0;j<=con->len;j++) {
			if (seg[j]) {
			  fprintf(seq,"%i",seqs[i][j]);
			  fi++;
			  if ((fi%50)==0) fprintf(seq,"\n");
			}
		}
	}
	fclose(seq);
	strcpy(fname, con->prefix);
	loc=fopen(strcat(fname, "loc"),"w");
	fprintf(loc,"%i %i L",ns,con->len);
	for (i=1;i<=con->len;i++) if (seg[i]) fprintf(loc,"\n%i",i);
	fclose(loc);
	free_ivector(seg,1,con->len);

}



/********************************************************/
/*Routine to find number of descendants for a given node*/
/********************************************************/

int count_desc(node, site)
struct node *node;
int site;
{
	int ndesc=0;
	if ((node->d[0] == NULL) && (node->d[1] == NULL)) ndesc=1;
	else if (node -> asite[site]) {
		if ((node->d[0] != NULL) && ((node->d[0])->asite[site])) 
			ndesc += count_desc(node->d[0],site);
		if ((node->d[1] != NULL) && ((node->d[1])->asite[site])) 
			ndesc += count_desc(node->d[1],site);
	}
	else ndesc+=0;
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

