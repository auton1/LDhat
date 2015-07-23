#if !defined COMPLK_H
#define COMPLK_H
#pragma warning(disable:4786) 
#include <math.h>

#include "rhomap_tools.h"
#include "data.h"

const int ADD=10000; //Size of extra pait types to be added when more memory needed

struct site_type 
{
  int pt[16];	// Haplotype pair type
  int nt;		// Number of such type in data
};

int check_exhaustive(struct site_type **pset, int npt, int nsamp);
int check22(int s1, int s2, int **nall) ;
int add_type(struct site_type **pset, int *cpt, int *ntc, int *pnew, int *miss, data *mydata) ;
void print_pairs(FILE *ofp, struct site_type **pset, int nt, int hd, int nseq) ;
int * order_pt_hap(int *pt, int nseq) ;
int * order_pt_dip(int *pt, int nseq);
void type_print(int **pij, int lseq, int w, FILE *ofp) ;
void read_pt(FILE *ifp, struct site_type **pset, int *npt, data *mydata) ;
struct site_type ** init_pset(struct site_type **pset, int lkf, FILE *ifp, int *npt, data *mydata) ;
struct site_type ** add_pset(struct site_type **pset) ;
void read_pars(FILE *ifp, int *tcat, double *theta, int *rcat, double *rmax) ;
void lk_resolve(double *lkres, struct site_type *pset, double *lknew, double **lkmat,data *mydata, double *lnfac_array, int rcat);
void lk_miss(struct site_type *pset,double *lkmiss,double **lkmat,data *mydata, double *lnfac_array, int rcat);

#endif
