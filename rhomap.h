#if !defined RHOMAP_H
#define RHOMAP_H
#pragma warning(disable:4786)
#include <cstring>
#include <string>
#include <math.h>
#include <time.h>

using namespace std;

#include "MCMC.h"
#include "data.h"
#include "likelihood.h"

const double VERSION=1.0;

int sizeofpset = 100;
long *idum;

void print_help(int argc, char* argv[]);
void read_arguments(int argc, char* argv[], data *mydata, MCMC *myMCMC, likelihood *myLikelihood);
void output_parameters_to_file(int argc, char *argv[], const data *mydata, const MCMC *myMCMC, const likelihood *myLikelihood);

#endif

