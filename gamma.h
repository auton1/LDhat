#if !defined GAMMA_H
#define GAMMA_H

#define ACC 0.001
#define EPSg 3.0e-7
#define ITMAX 100

double gammln(double xx);
void gcf(double *gammcf, double a, double x, double *gln);
void gser(double *gamser, double a, double x, double *gln);
double gammp(double a, double x);
double gammq(double a, double x);
double rgamm(double a);

#endif
