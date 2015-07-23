
#include "rhomap.h"

int main(int argc, char* argv[])
{
	print_help(argc, argv);		// Output help to user if required

	clock_t start, finish;
	start = clock();

	cout << "Starting rate estimation" << endl;

	data mydata;				// Object to hold data
	MCMC myMCMC;				// Object to perform MCMC
	likelihood myLikelihood;	// Object that holds likelihood tables etc.

	read_arguments(argc, argv, &mydata, &myMCMC, &myLikelihood);	// Read in command line arguments.

	srand((unsigned int) myMCMC.seed);

	mydata.readseqs();	// Read in sequences
	mydata.readlocs();	// Read in SNP loci

	mydata.allele_count();	// Calculate some statistics about the sequences

	myLikelihood.init(&mydata);			// Initalise the likelihood object
	myLikelihood.read_lk_file(&mydata);	// Read in the likelihood file

	//output_parameters_to_file(argc, argv, &mydata, &myMCMC, &myLikelihood);	// Output parameters (useful for debugging)

	myMCMC.rjMCMC(&mydata, &myLikelihood);	// Run the rjMCMC

	cout << "\n\nTotal Program Run Time" << endl;
	finish = clock();
	double duration = (double)(finish - start) / CLOCKS_PER_SEC;
	cout << duration << " seconds" << endl;

	// Stop exit when running Windows...
	#ifdef _MSC_VER
		printf("\n\nPress Return to exit\n\n");
		int ch = getchar();
	#endif

	cout << "\n" << endl;

	return 0;
}


void print_help(int argc, char* argv[])
{	// See if user is asking for help and output in needed
	int i;
	char *in_str;
	for(i = 0; i < argc; i++)
	{
		if(*argv[i] == '-')
		{ 
			in_str = argv[i];
			if((strcmp(in_str, "-h") == 0) || (strcmp(in_str, "-?") == 0) || (strcmp(in_str, "-help") == 0)) 
			{
				printf("\nrhomap %0.1f\n", VERSION);
				printf(" Sample recombination rate profiles using flattened composite likelihood\n");
				printf(" and a hotspot prior model.\n");
				printf(" Auton and McVean, Genome Research (2007)\n\n");
	
				printf(" Required Options :\n");
				printf(" -seq <file>               Sequence data file\n");
				printf(" -loc <file>               SNP positions data file\n");
				printf(" -lk <file>                Composite likelihood lookup data file\n");
				printf(" -its <int>                Number of MCMC iterations\n");
				printf(" -samp <int>               Number of MCMC iterations between samples\n");
				printf(" -burn <int>               Number of MCMC burn-in iterations\n");
				printf("\n\n Additional Options :\n");
				printf(" -exact                    Likelihood file exact?\n");
				printf(" -seed <int>               User defined random seed\n");
				printf(" -no_files                 Do not output MCMC samples\n");
				printf(" -prefix <string>          Prefix for ALL output files\n");

				printf("\n\n");	
				printf(" Prior Parameters (defaults as in Auton and McVean 2007):\n");
				printf(" -bpen <float>             Background block penalty (default 0)\n");
				printf(" -hpen <float>             Hotspot penalty (default 0)\n");
				printf(" -bgAlpha <float>          Background rate prior alpha (Gamma Distribution)\n");
				printf(" -bgBeta <float>           Background rate prior beta (Gamma Distribution)\n");
				printf(" -HeatAlpha <float>        Hotspot heat prior alpha (Gamma Distribution)\n");
				printf(" -HeatBeta <float>         Hotspot heat prior beta (Gamma Distribution)\n");
				printf(" -ScaleAlpha <float>       Hotspot scale prior alpha (Gamma Distribution)\n");
				printf(" -ScaleBeta <float>        Hotspot scale prior beta (Gamma Distribution)\n");
				printf(" -T <float>                Expected Hotspot seperation\n");
				printf(" -m <float>                Maximum extent of hotspots (default 5kb)\n");

				// Commented out options are unsupported.
//				printf("\n\n");	
//				printf(" Advanced Options :\n");
//				printf(" -verbose                  Output additional information to file\n");
//				printf(" -nonverbose               Output less information to screen\n");
//				printf(" -freqs <file>                Output allele frequencies to file\n");
//				printf(" -paramfile <file>            Output parameter file\n");
//				printf(" -ratesfile <file>            Output rates file\n");
//				printf(" -hotspotfile <file>          Output hotspot file\n");
//				printf(" -blockfile <file>            Output block file\n");
//				printf(" -hotMoveSigma <float>     Standard deviation of hotspot move\n");
//				printf(" -print_blocks             Print Background blocks\n");
//				printf(" -print_hotspots           Print Hotspots blocks\n");
//				printf(" -w <int>                  Composite Likelihood SNP window size\n");
//				printf(" -oldlk                    Use old uncorrected composite likelihood\n");
//				printf(" -no_hotspots              No hotspots\n");
//				printf(" -nolk                     Do not use likelihood. Used to test MCMC\n");
//				printf(" -seqsubset <int>          Choose random sequence subset\n");
//				printf(" -init_rate <float>        Initial rate. Sets MCMC starting point\n");

				my_exit(" ", 0);
			}
		}
	}
}

void read_arguments(int argc, char* argv[], data *mydata, MCMC *myMCMC, likelihood *myLikelihood)
{	// Read in arguments from command line
	int i;
	char* in_str;
	for(i = 0; i < argc; i++)
	{
		if(*argv[i] == '-')
		{ 
			in_str = argv[i];
			if(strcmp(in_str, "-seed") == 0) {myMCMC->seed = atol(argv[i + 1]);						// User defined random seed
											  if (myMCMC->seed > 0) myMCMC->seed = -myMCMC->seed;}
			if(strcmp(in_str, "-lk") == 0) myLikelihood->lk_filename = argv[i + 1]; 				// Composite Likelihood filename
			if(strcmp(in_str, "-seq") == 0) mydata->seq_filename = argv[i + 1];						// Sequence filename
			if(strcmp(in_str, "-loc") == 0) mydata->loc_filename = argv[i + 1];						// Locus filename
			if(strcmp(in_str, "-freq") == 0) {mydata->output_freqs = true; 
												mydata->freq_filename = argv[i + 1];}				// Freqs filename
			if(strcmp(in_str, "-paramfile") == 0) mydata->params_filename = argv[i + 1];			// Parameter filename
			if(strcmp(in_str, "-exact") == 0) mydata->lk_exact = true;								// Likelihood file exact? flag
			if(strcmp(in_str, "-verbose") == 0) mydata->verbose = true;								// Verbose mode
			if(strcmp(in_str, "-nonverbose") == 0) myMCMC->nonverbose = true;						// Non-verbose mode.
			if(strcmp(in_str, "-bpen") == 0) myMCMC->bpen = atof(argv[i + 1]);						// Block penalty
			if(strcmp(in_str, "-hpen") == 0) myMCMC->hpen = atof(argv[i + 1]);						// Hotspot penalty
			if(strcmp(in_str, "-its") == 0) myMCMC->n_update = atoi(argv[i + 1]);					// Number of iterations
			if(strcmp(in_str, "-samp") == 0) myMCMC->r_update = atoi(argv[i + 1]);					// Sample every?
			if(strcmp(in_str, "-burn") == 0) myMCMC->burn_in = atoi(argv[i + 1]);					// Burn In
			if(strcmp(in_str, "-no_hotspots") == 0) myMCMC->no_hotspots = true;						// No hotspots
			if(strcmp(in_str, "-no_blocks") == 0) myMCMC->no_blocks = true;						// Single background block
			if(strcmp(in_str, "-w") == 0) myLikelihood->MAXW = atoi(argv[i+1]);					// Composite Likelihood window size
			if(strcmp(in_str, "-print_blocks") == 0) myMCMC->print_blocks = true;
			if(strcmp(in_str, "-print_hotspots") == 0) myMCMC->print_hotspots = true;

			if(strcmp(in_str, "-ratesfile") == 0) myMCMC->rates_file = argv[i + 1];					// Rates filename
			if(strcmp(in_str, "-hotspotfile") == 0) myMCMC->hotspots_file = argv[i + 1];			// Hotspot filename
			if(strcmp(in_str, "-blockfile") == 0) {myMCMC->blocks_file = argv[i + 1]; 
														myMCMC->print_blocks = true;}		// Blocks filename

			if(strcmp(in_str, "-bgBeta") == 0) myMCMC->blockparam.RATE_BETA = atof(argv[i + 1]);			// Mean prior BG rate
			if(strcmp(in_str, "-bgAlpha") == 0) myMCMC->blockparam.RATE_ALPHA = atof(argv[i + 1]);			// Mean prior BG rate
			if(strcmp(in_str, "-HeatAlpha") == 0) myMCMC->hotspotparam.HEAT_ALPHA = atof(argv[i + 1]);		// Hotspot heat prior alpha
			if(strcmp(in_str, "-HeatBeta") == 0) myMCMC->hotspotparam.HEAT_BETA = atof(argv[i + 1]);		// Hotspot heat prior beta
			if(strcmp(in_str, "-ScaleAlpha") == 0) myMCMC->hotspotparam.MU_ALPHA = atof(argv[i + 1]);		// Hotspot scale prior alpha
			if(strcmp(in_str, "-ScaleBeta") == 0) myMCMC->hotspotparam.MU_BETA = atof(argv[i + 1]);			// Hotspot scale prior beta
			if(strcmp(in_str, "-T") == 0) myMCMC->hotspotparam.AVG_DIST_BETWEEN_HOTSPOTS = atof(argv[i + 1]);	// Expected distance between hotspots
			if(strcmp(in_str, "-m") == 0) myMCMC->hotspotparam.MAX_WIDTH = atof(argv[i + 1]);	// Expected distance between hotspots

			if(strcmp(in_str, "-oldlk") == 0) myLikelihood->old_lk = true;								// Use uncorrected version of likelihood.
			if(strcmp(in_str, "-nolk") == 0) myLikelihood->no_lk = true;								// Don't calculate likelihood - used to test MCMC
			if(strcmp(in_str, "-seqsubset") == 0) mydata->random_seq_subset_size = atoi(argv[i+1]);		// Set size of sequence subset
			
			if(strcmp(in_str, "-hotMoveSigma") == 0) myMCMC->hotspotparam.MOVE_SIGMA = atof(argv[i + 1]);	// Hotspot move sigma
			if(strcmp(in_str, "-init_rate") == 0) myMCMC->init_rate = atof(argv[i + 1]);						// MCMC initialisation rate
			if(strcmp(in_str, "-no_files") == 0) myMCMC->no_files = true;								// Don't output MCMC samples
		}
	}

	idum = &(myMCMC->seed);
	cout << "Seed = " << myMCMC->seed << endl;

	char tempchar[255];
	if (mydata->seq_filename == "") { printf("\n\nInput sequence filename:"); scanf("%s", tempchar); mydata->seq_filename=tempchar; }
	if (mydata->loc_filename == "") { printf("\n\nInput locus filename:"); scanf("%s", tempchar); mydata->loc_filename=tempchar; }
	if (myLikelihood->lk_filename == "") { printf("\n\nInput likelihood filename:"); scanf("%s", tempchar); myLikelihood->lk_filename=tempchar; }
	if (myMCMC->bpen == -1 ) { printf("\n\nInput block penalty:"); scanf("%lf", &myMCMC->bpen); }
	if (myMCMC->hpen == -1 ) { printf("\n\nInput hotspot penalty:"); scanf("%lf", &myMCMC->hpen); }
	if (myMCMC->n_update == -1 ) { printf("\n\nHow many updates for MCMC? ");	scanf("%i", &myMCMC->n_update); }
	if (myMCMC->r_update == -1 ) { printf("\n\nNumber of updates between samples:"); scanf("%i", &myMCMC->r_update); }
	if (myMCMC->burn_in == -1 ) { printf("\n\nInput number of Burn In iterations:"); scanf("%i", &myMCMC->burn_in); }

	for(i = 0; i < argc; i++)
	{	// Prefix filenames - do this after everything else to allow user defined filenames.
		if(*argv[i] == '-')
		{ 
			in_str = argv[i];
			if ((strcmp(in_str, "-fileprefix") == 0)  || (strcmp(in_str, "-prefix") == 0))
			{
				string prefix = argv[i+1];
				mydata->freq_filename = prefix + mydata->freq_filename;
				mydata->params_filename = prefix + mydata->params_filename;
				myMCMC->rates_file = prefix + myMCMC->rates_file;
				myMCMC->hotspots_file = prefix + myMCMC->hotspots_file;
				myMCMC->blocks_file = prefix + myMCMC->blocks_file;
				myMCMC->acceptance_rates_file = prefix + myMCMC->acceptance_rates_file;
				myMCMC->summary_file = prefix + myMCMC->summary_file;
				myLikelihood->lk_out_filename = prefix + myLikelihood->lk_out_filename;
			}
		}
	}
}

void output_parameters_to_file(int argc, char *argv[], const data *mydata, const MCMC *myMCMC, const likelihood *myLikelihood)
{	// Output certain parameters to file - useful for debugging
	int i;
	FILE *params;
	params = fopen(mydata->params_filename.c_str(), "w");
	fprintf(params, "Input string: ");
	for (i=0; i<argc; i++)
		fprintf(params, "%s ", argv[i]);
	fprintf(params, "\n");
	fprintf(params, "Sequence file: %s\n", mydata->seq_filename.c_str());
	fprintf(params, "Locus file: %s\n", mydata->loc_filename.c_str());
	fprintf(params, "LK file: %s\n", myLikelihood->lk_filename.c_str());
	fprintf(params, "Seed: %li\n", myMCMC->seed);
	fprintf(params, "Block Penalty: %f\n", myMCMC->bpen);
	fprintf(params, "Hotspot Penalty: %f\n", myMCMC->hpen);
	fprintf(params,"Total Updates: %i\n", myMCMC->n_update);
	fprintf(params,"Updates Between Samples: %i\n", myMCMC->r_update);
	fprintf(params,"Burn In: %i\n", myMCMC->burn_in);

	if (mydata->random_seq_subset_size > 0)
	{
		fprintf(params, "\n\nSequence Subset Numbers\n");
		for (i=1; i<mydata->nseq+1; i++)
			fprintf(params, "%i\n", mydata->order[i]);
	}
	
	fclose(params);
}



