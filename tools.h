#if !defined TOOLS_H
#define TOOLS_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double *dvector();
double **dmatrix();
int *ivector();
int **imatrix();
char **cmatrix();
void free_dvector();
void free_ivector();
void free_dmatrix();
void free_imatrix();
void free_cmatrix();
void nrerror(const char error_text[]);
int mini();
int maxi();
double minc();
double mind();
double maxd();
double lnfac();
void pswap();
double lognC2();
double lognC4();
void sort();
long setseed();
double ran2();

int rpoiss();

#endif


