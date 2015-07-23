#include "rhomap_tools.h"

void nrerror(const char error_text[])
{
	fprintf(stderr,"Run-time error...\n\n");
	fprintf(stderr,"%s\n\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	my_exit("Bye!", -1);
}

double lnfac(int i) 
{

  int j;
  double cp=0.0;

  for (j=2;j<=i;j++) cp += log((double) j);
  return cp;
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

	for (i=a+1;i<=n;i++) x += log(i);
	for (i=2;i<=n-a;i++) x -= log(i);

	if (x>0) return x;
	else nrerror("Error in logCn2 (2)");
	return -1.0;
}

double lognC4(int n, int a, int b, int c, int d)
{
	int i;
	double x=0;

	if (a+b+c+d != n) nrerror("Error in logCn4 (1)");
	if ((n<2)||(a==n)||(b==n)||(c==n)||(d==n)) return 0.0;

	for (i=a+1;i<=n;i++) x += log(i);
	for (i=2;i<=b;i++) x -= log(i);
	for (i=2;i<=c;i++) x -= log(i);	
	for (i=2;i<=d;i++) x -= log(i);

	if (x>0) return x;
	else nrerror("Error in lognC4 (2)");
	return -1.0;
}

__inline double Add_Log(double Summand1, double Summand2)
{
	if (Summand1 < Summand2)
		return Summand2 + log(1.0 + exp(Summand1 - Summand2));
	else
		return Summand1 + log(1.0 + exp(Summand2 - Summand1));
}


// This funciton takes in two logged numbers and outputs the log of the difference of the numbers.
__inline double Subtract_Log(double Summand, double subtrahend)
{
	if (Summand <= subtrahend)
	{	
		return HUGE_VAL;
	}
	return Summand + log(1.0 - exp(subtrahend - Summand));
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
	return (long)lt;
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
	double EPS = 1.2e-7;
	double RNMX=(1.0-EPS);

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

// Returns LnGam(xx)
// From Numerical Recipes in C
double gammln(double xx)
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
	24.01409824083091,-1.231739572450155,
	0.1208650973866179e-2,-0.5395239384953e-5};
	int j;
	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}

// Returns a factorial
double factorial(int n)
{
	static int ntop=4;
	static double a[33]={1.0,1.0,2.0,6.0,24.0};
	int j;
	if (n < 0) nrerror("Negative factorial in routine factrl");
	if (n > 32) return exp(gammln(n+1.0));
	while (ntop<n) { 
		j=ntop++;
		a[ntop]=a[j]*ntop;
	}
	return a[n];
}

double gamma_function(int xx)
{
	return factorial(xx - 1);
}

// Normal Distrubted random number generator
// From Numerical Recipes in C, Press et al. 1988
// Returns a normally distributed random number with zero mean and unit variance
double normrnd(void)
{
	extern long *idum;
	static int iset=0;
	static double gset;
	double fac, r, v1, v2;

	if (iset == 0)
	{
		do
		{
			v1 = 2.0 * ran2() - 1.0;
			v2 = 2.0 * ran2() - 1.0;
			r = v1 * v1 + v2 * v2;
		} while (r >= 1.0);
		
		fac = sqrt(-2.0*log(r) / r);
		gset = v1 * fac;
		iset = 1;
		return v2*fac;
	} 
	else
	{
		iset = 0;
		return gset;
	}
}

// Returns random numbers distributed as a gamma distribution
// From Devroye 1986
double gamrnd(double ia, double ib)
{
	double b, c, d, u, v, w, x, y, z;
	int accept;
	double ret = 0;

	if (ia == 1.0)
	{
		// gamma is exponetial (Devroye pg 405)
		ret = -ib * log(ran2());
	}
	else if ((ia < 1) && (ia > 0))
	{
		c = 1 / ia;
		d = 1 / (1 - ia);
		accept = 0;
		do 
		{
			u = ran2();
			v = ran2();
			x = pow(u, c);
			y = pow(v, d);
			z = x + y;
			if (z <= 1.0)
				accept = 1;
		} while (accept != 1);

		ret = -ib * log(ran2()) * x / z;
	}
	else if (ia > 1)
	{
		b = ia - 1;
		c = 3.0 * ia - 0.75;
		accept = 0;
		do
		{
			u = ran2();
			v = ran2();
			w = u * (1 - u);
			y = (u - 0.5) * sqrt(c / w);
			x = b + y;
			if (x >= 0.0)
			{
				z = 64.0 * pow(w, 3) * pow(v, 2);
				if (z <= (1 - 2 * y * y / x))
				{ accept = 1; }
				else
				{	
					if (log(z) <= (2 * (b * log(x / b) - y)))
					{ accept = 1; }
				}
			}
		} while (accept != 1);

			ret = ib * x;
	}
	return ret;
}

// Returns the density of a gamma pdf at x
double log_gamma_pdf(double alpha, double beta, double x)
{
	if (x < 0)
		return 0;
	double out;
	
	out = (alpha - 1.0) * log(x) - alpha*log(beta) - gammln(alpha) - (x / beta);
	return out;
}


// Poisson Distrubted random number generator
// For Numerical Recipes in C, Press et al. 1988
// Returns a poisson distributed random number with mean xm

int poissrnd(double xm)
{
	static double sq, alxm, g, oldm=(-1.0);
	double em,t,y;
	
	if (xm < 12.0) { 
		if (xm != oldm) {
			oldm=xm;
			g=exp(-xm); 
		}
		em = -1;
		t=1.0;
		do {
			++em;
			t *= ran2();
		} while (t > g);
	} else { 
		if (xm != oldm) { 
			oldm=xm;
			sq=sqrt(2.0*xm);
			alxm=log(xm);
			g = xm*alxm - gammln(xm + 1.0);
		}
		do {
			do { 
				y=tan(PI * ran2());
				em=sq*y+xm; 
			} while (em < 0.0); 
			em=floor(em); 
			t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
		} while (ran2() > t);
	}
	return (int)em;
}

double exprnd(double mean)
{
	return -mean * log(ran2());
}


double betai(double a, double b, double x)
{
double bt;
if (x < 0.0 || x > 1.0) nrerror("Bad x in routine betai");
if (x == 0.0 || x == 1.0) 
	bt=0.0;
else 
	bt=exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.0-x));
if (x < (a+1.0)/(a+b+2.0)) 
	return bt*betacf(a,b,x)/a;
else 
	return 1.0-bt*betacf(b,a,1.0-x)/b;
}

double betacf(double a, double b, double x)
{
	double EPS = 3.0e-7;
	int m,m2;
	double aa,c,d,del,h,qab,qam,qap;
	qab=a+b; 
	qap=a+1.0;
	qam=a-1.0;
	c=1.0; 
	d=1.0-qab*x/qap;
	if (fabs(d) < FPMIN) d=FPMIN;
	d=1.0/d;
	h=d;
	for (m=1;m<=MAXIT;m++) {
	m2=2*m;
	aa=m*(b-m)*x/((qam+m2)*(a+m2));
	d=1.0+aa*d; 
	if (fabs(d) < FPMIN) d=FPMIN;
	c=1.0+aa/c;
	if (fabs(c) < FPMIN) c=FPMIN;
	d=1.0/d;
	h *= d*c;
	aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
	d=1.0+aa*d; 
	if (fabs(d) < FPMIN) d=FPMIN;
	c=1.0+aa/c;
	if (fabs(c) < FPMIN) c=FPMIN;
	d=1.0/d;
	del=d*c;
	h *= del;
	if (fabs(del-1.0) < EPS) break; 
	}
	if (m > MAXIT) nrerror("a or b too big, or MAXIT too small in betacf");
	return h;
}

double cumulative_binomial_prob(int success, int trials, double p)
{
	return betai(success, trials, p);
}


double *dvector(int nh)
{
	int i;
	double *v = new double[nh];

	for (i=0;i<nh; i++)
		v[i] = 0.0;

	return v;
}

void free_dvector(double *v)
{
	delete [] v;
}

int *ivector(int nh)
{
	int i;
	int *v = new int[nh];

	for (i=0;i<nh; i++)
		v[i] = 0;
	return v;
}

void free_ivector(int *v)
{
	delete [] v;
}

double **dmatrix(int nrh, int nch)
{
	int i,j;
	double **m;

	m = new double*[nrh];
	for (i=0;i<nrh; i++)
		m[i] = new double [nch];

	for (i=0;i<nrh; i++)
		for (j=0;j<nch; j++)
			m[i][j] = 0.0;

	return m;
}

void free_dmatrix(double **m, int nrh)
{
	int i;

	for(i=0;i<nrh;i++) 
		delete [] m[i];
	delete [] m;
}

int **imatrix(int nrh, int nch)
{
	int i,j;
	int **m;

	m = new int*[nrh];
	for (i=0;i<nrh; i++)
		m[i] = new int [nch];

	for (i=0;i<nrh; i++)
		for (j=0;j<nch; j++)
			m[i][j] = 0;
	return m;
}

void free_imatrix(int **m, int nrh)
{
	int i;

	for(i=0;i<nrh;i++) 
		delete [] m[i];
	delete [] m;
}


int my_finite(double x)
{
	#ifdef _MSC_VER
		return _finite(x);
	#else
		return finite(x);
	#endif
}

void my_exit(const char message[], int code)
{
	printf("\n%s\n", message);
	exit(code);
}



int compare (const void * a, const void * b)
{
  return ( *(int*)a - *(int*)b );
}


