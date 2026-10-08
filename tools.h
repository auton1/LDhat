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

double *dvector(int nl, int nh);
double **dmatrix(int nrl, int nrh, int ncl, int nch);
int *ivector(int nl, int nh);
int **imatrix(int nrl, int nrh, int ncl, int nch);
char **cmatrix(int nrl, int nrh, int ncl, int nch);
void free_dvector(double *v, int nl, int nh);
void free_ivector(int *v, int nl, int nh);
void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch);
void free_imatrix(int **m, int nrl, int nrh, int ncl, int nch);
void free_cmatrix(char **m, int nrl, int nrh, int ncl, int nch);
void nrerror(const char error_text[]);
int mini(int i, int j);
int maxi(int i, int j);
double minc(double l1, double l2, double ls);
double mind(double f1, double f2);
double maxd(double f1, double f2);
double lnfac(int i);
void pswap(int *pt, int s1, int s2);
double lognC2(int n, int a);
double lognC4(int n, int a, int b, int c, int d);
void sort(double *array, int ne);
long setseed(void);
double ran2(void);

int rpoiss(double x);

#endif


