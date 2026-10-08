#if !defined FIN_H
#define FIN_H


#define EVENT_PRINT 0

void print_help(int argc, char* argv[]);

struct node_tree
{
	int node_num;
	int site;
	float time;
	float time_above;
	float time_below;	/* Floats save memory */
	int ndesc;
	struct node_tree *d[2];
	struct node_tree *a;
	unsigned short nuc;	/* Short saves memory */
};



struct node_list
{
	int node_num;
	struct node_tree **asite;
	short nanc;	/* Is ancestral material? */
	double rlen;
	double time;
};


struct results {
	int nrun;
	double nm;
	double pwd;
	double sn;
	double nr;
	double mhm;
	double cf;
	double cf2;
	double *fdist;
	double covG[10];   /*For genealogical covariances*/
	double r2[3];
    double Dp;
	double G4;
	double cfld;
    double cfld2;
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
	int conv;
	double clen;
	double cratio;
	double p_genotype_error;
	double p_bad_site;
	double p_switch_error;
	int asc;
  	int gt;
	double rescale;	/* Rescale output loci file */
  	char prefix[127];
  	char rmap_file[127];
  	char flocs_file[127];
  	char mut_file[127];
};

#define BACC 1e-6

struct node_tree *** make_tree(struct node_tree ***tree_ptr, struct control *con, int *revts, int **segl);

void set_res(struct results *res, struct control *con);
void print_lin(struct node_list **list, int len, int k);
void  count_rlen(struct node_list *nodel, struct control *con);
void print_nodes(struct node_tree **tree, int len);
void tree_summary(struct node_tree ***tree, struct control *con, struct results *res, int **seqs);
double tree_time(struct node_tree *node);
int add_mut(struct node_tree **tree, double *fl);
struct node_tree * add_mut_f(int fsim, struct control *con, double *cf, struct node_tree **tree_site);
void seq_mut(struct node_tree *nm, int **seqs, int site, int base);
void print_seqs(int **seqs, struct control *con);
void print_res(struct results res, struct control con);
void read_input(struct control *con, int argc, char *argv[]);
void read_flags(struct control *con, int argc, char *argv[]);
void select_base(int *nb, int base, double **mut_mat);
void evolve(struct node_tree *np, int **seqs, struct control *con, double *mm, double **mutmat, int *muts, int site);
char num_to_nuc(int i);
int count_desc(struct node_tree *node);

int add_genotype_error(int **seqs, struct control *con);
int add_bad_sites(int **seqs, struct control *con);
int add_switch_error(int **seqs, struct control *con);
int remove_sites_by_frequency(int **seqs, struct control *con);

void choose_time(double *t, int k, double rho, struct control *con);
double bisect(double (*bfunc)(double *, double **), double val, double *cons);
double tgrowth(double *var, double **cons);

struct node_list ** recombine(int *k, double *rr, struct control *con, struct node_list **list, double *t, FILE *ofp);
struct node_list ** coalesce(int *k, struct control *con, struct node_list **list, struct node_tree ***tree_ptr, int **segl, double *t, FILE *ofp);

void check_lin(struct node_list **list, int *k, struct control *con);


#endif


