#if !defined SEQTOOLS_H
#define SEQTOOLS_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "ldhat.h"
#include "tools.h"
#include "seqtools.h"

int read_fasta();
void allele_count();
double watterson();
int check22();
struct site_type ** pair_spectrum();
int add_type();
void print_pairs();
int * order_pt_hap();
int * order_pt_dip();
void type_print();
void read_pt();
struct site_type ** init_pset();
struct site_type ** add_pset();
void read_pars();
void read_lk();

#endif


