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

struct node_tree *** make_tree();

void set_res();
void print_lin();
void  count_rlen();
void print_nodes();
void tree_summary();
double tree_time();
int add_mut();
struct node_tree * add_mut_f();
void seq_mut();
void print_seqs();
void print_res();
void read_input();
void read_flags();
void select_base();
void evolve();
char num_to_nuc();
int count_desc();

int add_genotype_error();
int add_bad_sites();
int add_switch_error();
int remove_sites_by_frequency();

void choose_time();
double bisect();
double tgrowth();

struct node_list ** recombine();
struct node_list ** coalesce();

void check_lin();


#endif


