#if !defined SEQTOOLS_H
#define SEQTOOLS_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "ldhat.h"
#include "tools.h"
#include "seqtools.h"

int read_fasta(int **seqs, FILE *ifp, int nseq, int lseq, char **seqnames);
void allele_count(int **seqs, int nseq, int lseq, int **nall, int fl, int hd, char *prefix);
double watterson(int n);
int check22(int s1, int s2, int **nall);
struct site_type ** pair_spectrum(int **seqs, struct data_sum *data, int **nall, struct site_type **pset, int *npt, int *pnew, int *miss, int anc, int **pij);
int add_type(struct site_type **pset, int *cpt, int *ntc, int *pnew, int *miss, struct data_sum *data);
void print_pairs(FILE *ofp, struct site_type **pset, int nt, int hd, int nseq);
int * order_pt_hap(int *pt, int nseq);
int * order_pt_dip(int *pt, int nseq);
void type_print(int **pij, int lseq, int w, FILE *ofp);
void read_pt(FILE *ifp, struct site_type **pset, int *npt, struct data_sum *data);
struct site_type ** init_pset(struct site_type **pset, int lkf, FILE *ifp, int *npt, struct data_sum *data);
struct site_type ** add_pset(struct site_type **pset);
void read_pars(FILE *ifp, int *tcat, double *theta, int *rcat, double *rmax);
void read_lk(FILE *ifp, double **lkmat, int npt, int tcat, int rcat);

#endif


