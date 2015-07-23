#if !defined LDHAT_H
#define LDHAT_H

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>


#define DEBUG 0 /*Debug option*/


/*Constants for composite likelihood part*/
#define NSHUFF 1000 /*Number of permutations in tests for recombination*/
#define NRUN 1000000 /*Number of proposals in IS estimation of coalescent likelihoods*/
#define ADD 10000 /*Size of extra PTs to be added when more memory needed*/
#define SEQ_MAX 1000/*Max number of sequences*/
#define MAXNAME 65535 /*Max length of sequences names*/
#define MAXLINE 65535 /*Max length of line*/
#define MAXW 50 /*MAXW*2 = Max number of SNPs to consider for likelihood - i.e. ignore SNPS > MAXW apart*/
#define BURNIN 100000


/*Constants for block routines*/
#define RMIN 0.00001 /*Min per unit rho*/
#define RMAX 250 /*Max per unit rho*/

struct site_type {

  int pt[16];/*Haplotype pair type*/
  int nt;/*Number of such type in data*/
  double ld_stat[3];/*LD statistics for pair type*/
  int miss;/*Integer indicating whether PT contains missing data*/
  double lkptmx;/*Max Likelihood for PT*/
  double rmpt;/*Rho_max for PT*/

  int rm; /*Min number of recombination events for pair - 0 or 1*/
};

struct data_sum {

  int nseq;/*Number of sequences*/
  int lseq;/*No. segregating sites*/
  double tlseq;/*Total length of sequence (may be in kb)*/
  int w;/*Size of window to be used in analysis*/
  int hd;/*Haploid (1) or diploid (2) data*/
  char lc;/*Crossing-over (L) or gene conversion (C) model*/
  int ptt;/*Number of pair types*/
  double avpwd;/*Average pairwise differences*/
  double varpwd;/*Sample variance in pairwise differences*/
  int rmin;/*Lower bound on minimum number of recombination events*/
  double rwak;/*4Ner estimated by Wakeley 1997*/
  double avc;/*Average conversion tract length*/
  double th;/*Theta per site*/
  double rho;/*Rho for whole gene (or gamma for conversion model)*/
  int rho_i;/*position of maximum rho*/
  double rho_drive;/*Rho to be used in driving simulations*/
  double fit_obs;/*Observed fit*/
  double *rmap;/*Array for recombination map*/
  double lkmax;/*Maximum composite likelihood*/
  double **lksurf;/*Array for likelihood surface*/
  double ld[4];/*LD statistics*/
  double rme;/*Number of points for rho in coalescent likelihood estimation*/
  double rmax;/*Max rho in coalescent likelihood estimation*/
  int rce;/*Maximum rho for estimation - can be >> rmax*/
  int rcat;/*Number of categories for estimating rho - can be >> rme*/
  double fit;/*Fit for simulations*/
  double clr;/*Composite likelihood ratio - from simulations*/
  int ng[2];/*Counters for P values in simulations*/
  int n_update;/*Number of updates in MCMC*/
  int r_update;/*Number of updates between samples from MCMC*/
  double bpen;/*Block penalty for MCMC*/
  int exact;/*Switch to speed up events if exact set of likelihoods will be inputed*/
  
  char prefix[MAXNAME+1];	/* Prefix for output filenames */
};

struct block {
  int num;/*position in array of active blocks*/
  double rate;/*Recombination rate (per kb) in blocks*/
  int pos;/*Which SNP the block starts at*/
  int size;/*Length of block in SNPs*/
  struct block *bpr;/*Pointer to RH block*/
  struct block *bpl;/*Pointer to LH block*/
};

void rec_test();
void ld_calc();
void lk_est();
void print_lks();
void lk_surf();
int lk_calc();

int lk_calc_win();
void print_par();
void print_lkres();
void ld_test();
void ld_calc2();
void check_exhaustive();
void lk_miss();
void lk_resolve();
void lk_win();
int lk_calc_site();
void freq_min();
void calc_nless1();


/*Routine for coalescent estimation of likelihoods - from Paul Fearnhead*/
void pairs();


/*Routines in MCMC block part*/
void lk_block();
void print_block();
void block2map();
void print_rmap();
struct block **update_blocks();
void check_blocks();
void update_lij();
void print_rates();
void print_bounds();
void check_lk();

/*Extra routines in pairdip*/
void fit_pwlk();
void rec_test();
void rmin();
int rec_event();
void wakeley_est();
double C_equation();



#endif

