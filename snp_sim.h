#if !defined SNPSIM_H
#define SNPSIM_H

struct node {
	int node_num;
	double time;
	struct node *d[2];
	struct node *a[2];
	int *asite;
	double rlen;
	int nuc;
};

struct results {
	int nrun;
	double nm;
	double pwd;
	double sn;
	double nr;
	double cf;
	double cf2;
};

struct control{
	int nsamp;
	int len;
	long int seed;
	double *theta;
	double *rmap;
	double R;
	int hyp;
	double phm;
	double rhm;
	double a;
	int mut;
	int inf;
	int nrun;
	int print;
	int fmin;
	int cond;
    int **slocs;
    int growth;
    double lambda;
	int bneck;
	double tb;
	double strb;
    double w0;
	int rm;

	double cl_fac;
	char prefix[127];
};

#define BACC 1e-6




void snp_sim(double *locs, int *flocs, struct site_type **pset, double **lkmat, int nrun, struct data_sum *data);
struct node ** make_tree(struct node **tree_ptr, struct control *con, int *revts, int **segl);
void print_lin(struct node **list, int len, int k);
void  count_rlen(int *asite, int len, double *rlen, struct control *con);
void print_nodes(struct node **tree, int nn, int len);
void tree_summary(struct node **tree, int nn, struct control *con, int **segl, struct data_sum *data, struct site_type **pset, double **lkmat, double *locs, FILE *ofp);
struct node ** add_tree(struct node **tree, int n_node, int len);
double tree_time(struct node *node, int site);
int add_mut(struct node *node, int site, double *fl);
int add_mut_f(int fsim, int nnode, int nsamp, double *cf, struct node **tree, int site, double theta);
void seq_mut(struct node *nm, int **seqs, int nsamp, int site, int base);
void print_seqs(int **seqs, struct control *con);
char num_to_nuc(); /* declared but never defined or called in this program (only in fin.c under fin.h) - left as-is */
int count_desc(struct node *node, int site);
void choose_time(double *t, int k, double rho, struct control *con);
double bisect(double (*bfunc)(double *, double **), double val, double *cons);
double tgrowth(double *var, double **cons);

void recombine(int *k, double *rr, struct control *con, struct node **list, struct node **tree_ptr, double *t);
void coalesce(int *k, struct control *con, struct node **list, struct node **tree_ptr, int **segl, double *t);

#endif


