#include "tools.h"

void nrerror(const char error_text[])
{
	fprintf(stderr,"\n\nRun-time error\n\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"\n..now exiting to system...\n");
	exit(1);
}



int *ivector(int nl, int nh)
{
	int *v;

	v=(int *)malloc((unsigned) (nh-nl+1)*sizeof(int));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl;
}

double *dvector(int nl, int nh)
{
	double *v;

	v=(double *)malloc((unsigned) (nh-nl+1)*sizeof(double));
	if (!v) nrerror("allocation failure in dvector()");
	return v-nl;
}


double **dmatrix(int nrl, int nrh, int ncl, int nch)
{
	int i;
	double **m;

	m=(double **) malloc((unsigned) (nrh-nrl+1)*sizeof(double*));
	if (!m) nrerror("allocation failure 1 in dmatrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(double *) malloc((unsigned) (nch-ncl+1)*sizeof(double));
		if (!m[i]) nrerror("allocation failure 2 in dmatrix()");
		m[i] -= ncl;
	}
	return m;
}

int **imatrix(int nrl, int nrh, int ncl, int nch)
{
	int i,**m;

	m=(int **)malloc((unsigned) (nrh-nrl+1)*sizeof(int*));
	if (!m) nrerror("allocation failure 1 in imatrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(int *)malloc((unsigned) (nch-ncl+1)*sizeof(int));
		if (!m[i]) nrerror("allocation failure 2 in imatrix()");
		m[i] -= ncl;
	}
	return m;
}

char **cmatrix(int nrl, int nrh, int ncl, int nch)
{
    int i;
	char **m;

	m=(char **)malloc((unsigned) (nrh-nrl+1)*sizeof(char*));
    if (!m) nrerror("allocation failure 1 in cmatrix()");
    m -= nrl;

    for(i=nrl;i<=nrh;i++) {
    	m[i]=(char *)malloc((unsigned) (nch-ncl+1)*sizeof(char));
        if (!m[i]) nrerror("allocation failure 2 in cmatrix()");
        m[i] -= ncl;
    }
    return m;
}

void free_ivector(int *v, int nl, int nh)
{
	free((char*) (v+nl));
}

void free_dvector(double *v, int nl, int nh)
{
	free((char*) (v+nl));
}


void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch)
{
	int i;

	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}

void free_imatrix(int **m, int nrl, int nrh, int ncl, int nch)
{
	int i;

	for(i=nrh;i>=nrl;i--) if ((m[i]+ncl) != NULL) free((char*) (m[i]+ncl));
	if ((m+nrl) != NULL) free((char*) (m+nrl));
}


void free_cmatrix(char **m, int nrl, int nrh, int ncl, int nch)
{
        int i;

        for(i=nrh;i>=nrl;i--) if ((m[i]+ncl) != NULL) free((char*) (m[i]+ncl));
        if ((m+nrl) != NULL) free((char*) (m+nrl));
}



int mini(int i, int j)
{
        if (i<j) return i;
        else return j;
}

int maxi(int i, int j)
{
	if (i>j) return i;
	else return j;
}

double mind(double f1, double f2)
{
  if (f1<f2) return f1;
  else return f2;
}

double maxd(double f1, double f2)
{
  if (f1>f2) return f1;
  else return f2;
}

double lnfac(int i)
{
  int j;
  double cp=0.0;

  for (j=2;j<=i;j++) cp += log((double) j);
  return cp;
}


double minc(double l1, double l2, double ls)
{
        double d;
        if ((d=(l2-l1)) >  (ls*0.5)) d = (double) ls-l2+l1;
        return d;
}


void pswap(int *pt, int s1, int s2)
{
        int tmp;

        tmp = pt[s2];
        pt[s2]=pt[s1];
        pt[s1]=tmp;
}


double lognC2(int n, int a)
{
	int i;
	double x=0;

	if (a>n) nrerror("Error in lognC2 (1)");
	if ((n<2)||(a==0)||(a==n)) return 0.0;

	for (i=a+1;i<=n;i++) x += (double) log(i);
	for (i=2;i<=n-a;i++) x -= (double) log(i);


	if (x>0) return x;
	else nrerror("Error in logCn2 (2)");
	return 0;
}

double lognC4(int n, int a, int b, int c, int d)
{
	int i;
	double x=0;

	if (a+b+c+d != n) nrerror("Error in logCn4 (1)");
	if ((n<2)||(a==n)||(b==n)||(c==n)||(d==n)) return 0.0;

	for (i=a+1;i<=n;i++) x += (double) log(i);
	for (i=2;i<=b;i++) x -= (double) log(i);
	for (i=2;i<=c;i++) x -= (double) log(i);
	for (i=2;i<=d;i++) x -= (double) log(i);

	if (x>0) return x;
	else nrerror("Error in lognC4 (2)");
	return 0;
}


void sort(double *array, int ne)
{
	int pass, i;
	double tmp;

	for (pass=1;pass<=ne;pass++)
		for (i=1;i<ne;i++)
			if (array[i+1]<array[i]) {
				tmp = array[i];
				array[i]=array[i+1];
				array[i+1]=tmp;
			}
}

/*Routine to set seed from clock*/
long setseed(void) {
	time_t lt;
	lt=time(NULL);
	return lt;
}

/*C routine random number generator from
Numerical Recipes in C: Press et al*/

double ran2(void) {
	int j;
	long k;
	extern long *idum;
	static long idum2=123456789;
	static long iy=0;
	static long iv[NTAB];
	double temp;

	if (*idum <= 0) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		idum2=(*idum);
		for (j=NTAB+7; j>=0; j--) {
			k=(*idum)/IQ1;
			*idum=IA1*(*idum-k*IQ1)-k*IR1;
			if (*idum < 0) *idum += IM1;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ1;
	*idum=IA1*(*idum-k*IQ1)-k*IR1;
	if (*idum < 0) *idum += IM1;
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;
	if (idum2 < 0) idum2 += IM2;
	j=iy/NDIV;
	iy=iv[j]-idum2;
	iv[j]=*idum;
	if (iy < 1) iy += IMM1;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}

int rpoiss(double x)
{

	int i=0;
	double r1=ran2(), cump=0.0, fn;

	fn = exp(-x);
	while (cump<r1){
		cump+=fn;
		i++;
		fn*=(double) x/i;
	}
	i--;
	return i;
}


