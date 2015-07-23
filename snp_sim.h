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




void snp_sim();
struct node ** make_tree();
void print_lin();
void  count_rlen();
void print_nodes();
void tree_summary();
struct node ** add_tree();
double tree_time();
int add_mut();
int add_mut_f();
void seq_mut();
void print_seqs();
char num_to_nuc();
int count_desc();
void choose_time();
double bisect();
double tgrowth();

void recombine();
void coalesce();

#endif


