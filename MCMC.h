#if !defined MCMC_H
#define MCMC_H
#pragma warning(disable:4786) 
#include <math.h>
#include <vector>
#include <string>

using namespace std;
#include "params.h"
#include "block.h"
#include "rhomap_tools.h"
#include "likelihood.h"
#include "data.h"


class MCMC
{
public:
	long seed;						// Random seed.
	double bpen;					// Block penalty
	double hpen;					// Hotspot penalty
	int n_update;					// Number of iterations
	int r_update;					// Sample every?
	int burn_in;					// Burn In
	int run;						// Current iteration number

	struct block_params blockparam;		// Structure holding block prior parameters
	struct hotspot_params hotspotparam;	// Structure holding hotspot prior parameters

	string rates_file;				// Output filename
	string hotspots_file;			// Output filename
	string blocks_file;				// Output filename
	string acceptance_rates_file;	// Output filename
	string summary_file;			// Output filename

	bool rate_move, move_block_move, split_merge_move;
	bool heat_move, move_hot_move, insert_delete_hot_move, scale_move;
	double lk[2];					// Store the change in likelihood
	int nacc[9][2];					// Acceptance rates;

	double *rmap;					// Array for recombination map
	double *rmap_cache;				// Cache Array for recombination map
	double rmax;

	double *rmap_sum;				// Rmap sum (used for summary)
	double *hot_sum;				// Hotspot sum (used for summary)

	double log_theta;				// Pre-calculated hotspot prior constant (for performance)
	double log_gamma;				// Pre-calculated block prior mean (for performance)
	double theta;					// Pre-calculated hotspot prior constant (for performance)
	double gamma;					// Pre-calculated block prior mean (for performance)

	int numchangepoints;			// Number of blocks
	int numhotspots;				// Number of hotspots
	vector< block > blocks;			// Blocks vector

	bool no_hotspots;				// Flag to indicate that no hotspots allowed
	bool no_blocks;					// Flag to indicate only one background block
	bool nonverbose;				// Flag to indicate non-verbose mode.
	bool one_hotspot_per_block;			// Flag to only allow one hotspot per block

	int max_changepoints;			// Max number of changepoints
	int max_hotspots;				// Max number of hotspots

	double init_rate;				// Block initialisation rate

	double split_merge_constant;
	double hot_birth_death_constant;

	double log_block_prior_const;
	double log_hotspot_prior_const;

	bool print_blocks;				// Output blocks from chain to file?
	bool print_hotspots;			// Output hotspots from chain to file?

	bool no_files;					// Output chains to files?


	MCMC()
	{
		int i;
		run = 0;
		seed = -setseed();

		no_hotspots = false;
		no_blocks = false;

		print_blocks = false;
		no_files = false;

		one_hotspot_per_block = false;

		for (i=0;i<9;i++) nacc[i][0]=nacc[i][1]=0;	// Set Acceptance counts to zero

		rates_file = "rates.txt";
		hotspots_file = "hotspots.txt";
		blocks_file = "bgblocks.txt";
		summary_file = "summary.txt";

		acceptance_rates_file = "acceptance_rates.txt";

		numchangepoints = 0;
		numhotspots = 0;

		set_default_parameters(blockparam, hotspotparam); // Set default parameters
		log_theta = 0.0;
		
		bpen = 0; hpen = 0; n_update = -1; r_update = -1;	burn_in = -1;

		nonverbose = false;
		init_rate = 0.0;
	}

	~MCMC()
	{
		free_dvector(rmap);
		free_dvector(rmap_cache);
		free_dvector(rmap_sum);
		free_dvector(hot_sum);
	}

	void set_moves(bool rate, bool move_block, bool split_merge, 
							  bool heat, bool move_hot, bool insert_delete_hot, bool scale)
	{
		rate_move = rate;
		move_block_move = move_block;
		split_merge_move = split_merge;
		heat_move = heat;
		move_hot_move = move_hot;
		insert_delete_hot_move = insert_delete_hot;
		scale_move = scale;

		if (no_hotspots)
			heat_move = move_hot_move = insert_delete_hot_move = scale_move = false;
		if (no_blocks)
			move_block_move = split_merge_move = heat_move = move_hot_move = insert_delete_hot_move = scale_move = false;
	}


	void rjMCMC(data *mydata, likelihood *myLikelihood)
	{
		int i;
		FILE *rfp, *hotspotfile=NULL, *blockfile=NULL;

		rmap = dvector(mydata->lseq+1);
		rmap_cache = dvector(mydata->lseq+1);
		rmap_sum = dvector(mydata->lseq+1);
		hot_sum = dvector(mydata->lseq+1);
		rmax = myLikelihood->rmax;

		// Set constants for performance.
		log_theta = log(mind(mydata->tlseq / hotspotparam.AVG_DIST_BETWEEN_HOTSPOTS, mydata->lseq-2.0)) - hpen;
		theta = (mind(mydata->tlseq / hotspotparam.AVG_DIST_BETWEEN_HOTSPOTS, mydata->lseq-2.0))*exp(-hpen);
		// kb prior
		log_gamma = log(mind(mydata->tlseq, mydata->lseq - 2.0)) - bpen;
		gamma = mind(mydata->tlseq, mydata->lseq - 2.0) * exp(-bpen);
		// SNP prior
		//log_gamma = log(mydata->lseq-2) - bpen;
		//gamma = (mydata->lseq-2) * exp(-bpen);

		max_changepoints = mydata->lseq;	// Arbitary
		max_hotspots = mydata->lseq;

		double a1 = blockparam.RATE_ALPHA;
		double b1 = blockparam.RATE_BETA;
		log_block_prior_const = -a1*log(b1)-gammln(a1);

		a1 = hotspotparam.HEAT_ALPHA;
		b1 = hotspotparam.HEAT_BETA;
		log_hotspot_prior_const = -a1*log(b1)-gammln(a1);

		split_merge_constant = calc_split_merge_constant();
		hot_birth_death_constant = calc_hot_birth_death_constant();

		if (!no_files)
		{
			rfp = fopen(rates_file.c_str(),"w");
			if (print_hotspots) hotspotfile = fopen(hotspots_file.c_str(), "w");
			if (print_blocks) blockfile = fopen(blocks_file.c_str(), "w");
			fprintf(rfp,"%i %i\n",(int) ((double) ((n_update - burn_in)/r_update)), mydata->lseq);
		}

		if (no_blocks)
			numchangepoints = 0;
		else
			numchangepoints = mydata->lseq-1;
		blocks.resize(numchangepoints+1);

		set_moves(1, 1, 1, 0, 0, 0, 0);
//		set_moves(1, 0, 0, 1, 0, 1, 1);

		if (no_blocks)
			if (init_rate > 0.0)
			{
				if (init_rate <= blockparam.RMIN)
					init_rate = blockparam.RMIN+(blockparam.RMIN*0.1);
				blocks[0].init(init_rate, mydata->locs[1], mydata->locs[mydata->lseq], 1, mydata->lseq);
			}
			else
				blocks[0].init(&blockparam, mydata->locs[1], mydata->locs[mydata->lseq], 1, mydata->lseq);
			
		else
		{
			for(i=0; i<=numchangepoints; i++)
				if (init_rate > 0.0)
				{
					if (init_rate <= blockparam.RMIN)
						init_rate = blockparam.RMIN+(blockparam.RMIN*0.1);
					blocks[i].init(init_rate, mydata->locs[i+1], mydata->locs[i+2], i+1, i+1);
				}
				else
					blocks[i].init(&blockparam, mydata->locs[i+1], mydata->locs[i+2], i+1, i+1);
		}
		calc_rmap(mydata, 1, mydata->lseq);

		lk[0] = myLikelihood->lk_calc(1, mydata->lseq, mydata, rmap, 1);
		lk[1] = lk[0];
		myLikelihood->update_lij(1, mydata->lseq, mydata, 1);

		//myLikelihood->print_lij(mydata);

		if (!nonverbose) printf("\nInitial likelihood = %.3f\n",lk[0]);

		if (burn_in > 0) printf("\n\nRunning Burn-In\n");
		else 
		{	// No burn in, so allow hotspots from the start
			set_moves(1,1,1,1,1,1,1);
		//	set_moves(1, 0, 0, 1, 0, 1, 1);
		}

		// Begin main MCMC sequence
		for (run=1;run<=n_update;run++) 
		{
			do_move(mydata, myLikelihood);	// DO RJMCMC MOVES HERE!
			
			if (run == burn_in)
			{	// Completed burn-in. Output burn-in details to screen.
				lk[0] = myLikelihood->update_lk0(mydata);
				printf("\nFinal Burn In LK %.3f: ChangePoints %4i: Hotspots %4i: MapLen %.3f",lk[0],numchangepoints,numhotspots, rmap[mydata->lseq]);

				printf("\n\nBurn In Acceptance rates:\n");
				printf("Change Block Rate =   \t %.4f\t Move Changepoint =\t %.4f\n", (double) nacc[0][1]/nacc[0][0], (double) nacc[1][1]/nacc[1][0]);
				printf("Split Block =         \t %.4f\t Merge blocks =    \t %.4f\n", (double) nacc[2][1]/nacc[2][0], (double) nacc[3][1]/nacc[3][0]);
				printf("Change Hotspot Heat = \t %.4f\t Move Hotspot =    \t %.4f\n", (double) nacc[4][1]/nacc[4][0], (double) nacc[5][1]/nacc[5][0]);
				printf("Insert Hotspot =      \t %.4f\t Delete Hotspot =  \t %.4f\n", (double) nacc[6][1]/nacc[6][0], (double) nacc[7][1]/nacc[7][0]);
				printf("Change Hotspot Scale =\t %.4f\n",(double) nacc[8][1]/nacc[8][0]);

				for (i = 0; i < 9; i++)
				{	// Reset the move counters
					nacc[i][1] = 0;
					nacc[i][0] = 0;
				}
			}

			if (run == (int)burn_in * 0.25)
			{	// Time to allow hotspots
				set_moves(1,1,1,1,1,1,1);
				//set_moves(1, 0, 0, 1, 0, 1, 1);
			}

			if (!(run%(r_update)))
			{	// Output a sample
				lk[0] = myLikelihood->update_lk0(mydata);
				if (!nonverbose) printf("\nRun %8i: LK %.3f: ChangePoints %4i: Hotspots %4i: MapLen %.3f",run,lk[0],numchangepoints,numhotspots, rmap[mydata->lseq]);

				if (run >= burn_in)
				{
					if (!no_files)
					{
						print_rates(mydata, rfp);	// Output current progress
						print_blocks_and_hotspots(blockfile, hotspotfile);
					}
					store_sums(mydata);
				}
			}
		}

		printf("\nAfter %i runs - likelihood = %f\n\n",n_update, lk[0]);
		printf("Final block length = %f\n", rmap[mydata->lseq]);

		printf("\n\nFinal Acceptance rates:\n");
		printf("Change Block Rate =   \t %.4f\t Move Changepoint =\t %.4f\n", (double) nacc[0][1]/nacc[0][0], (double) nacc[1][1]/nacc[1][0]);
		printf("Split Block =         \t %.4f\t Merge blocks =    \t %.4f\n", (double) nacc[2][1]/nacc[2][0], (double) nacc[3][1]/nacc[3][0]);
		printf("Change Hotspot Heat = \t %.4f\t Move Hotspot =    \t %.4f\n", (double) nacc[4][1]/nacc[4][0], (double) nacc[5][1]/nacc[5][0]);
		printf("Insert Hotspot =      \t %.4f\t Delete Hotspot =  \t %.4f\n", (double) nacc[6][1]/nacc[6][0], (double) nacc[7][1]/nacc[7][0]);
		printf("Change Hotspot Scale =\t %.4f\n",(double) nacc[8][1]/nacc[8][0]);

		if (!no_files)
		{
			if (print_hotspots) fclose(hotspotfile);
			if (print_blocks) fclose(blockfile);
			fclose(rfp);
		}

		final_sum(mydata);

		output_acceptance_rates();
	}

	void store_sums(data *mydata)
	{
		int i, j, k;
		double pos;
		for(i=1; i<mydata->lseq; i++)
			rmap_sum[i] += rmap[i];

		for(j=0; j<numchangepoints; j++)
		{
			for(k=0; k<blocks[j].hotspots.size(); k++)
			{
				pos = blocks[j].hotspots[k].position;
				for(i=0; i<mydata->lseq; i++)
				{
					if((pos >= mydata->locs[i+1]) && (pos < mydata->locs[i+2]))
					{
						hot_sum[i+1]++;
						break;
					}
				}
			}
		}
	}

	void final_sum(data *mydata)
	{
		int i;
		for(i=1; i<mydata->lseq; i++)
			rmap_sum[i] = rmap_sum[i] / (((n_update - burn_in)/r_update)+1.0);

		for(i=1; i<mydata->lseq; i++)
			hot_sum[i] = hot_sum[i] / (mydata->locs[i+1] - mydata->locs[i]) / (((n_update - burn_in)/r_update) + 1.0);


		FILE *res;
		double rate;
		res = fopen(summary_file.c_str(), "w");
		fprintf(res, "Position(kb)\t4Ner\t4Ner/kb\tHotspot_Density\n");
		for (i=1; i<mydata->lseq; i++)
		{
			rate = (rmap_sum[i] - rmap_sum[i-1]) / (mydata->locs[i+1] - mydata->locs[i]);
			fprintf(res, "%10.4f\t%f\t%f\t%f\n", mydata->locs[i], rmap_sum[i-1], rate, hot_sum[i-1]);
		}
		fclose(res);
	}

	void output_acceptance_rates()
	{
		FILE *output = fopen(acceptance_rates_file.c_str(), "w");

		fprintf(output, "Change Block Rate =   \t %.4f\nMove Changepoint =\t %.4f\n", (double) nacc[0][1]/nacc[0][0], (double) nacc[1][1]/nacc[1][0]);
		fprintf(output, "Split Block =         \t %.4f\nMerge blocks =    \t %.4f\n", (double) nacc[2][1]/nacc[2][0], (double) nacc[3][1]/nacc[3][0]);
		fprintf(output, "Change Hotspot Heat = \t %.4f\nMove Hotspot =    \t %.4f\n", (double) nacc[4][1]/nacc[4][0], (double) nacc[5][1]/nacc[5][0]);
		fprintf(output, "Insert Hotspot =      \t %.4f\nDelete Hotspot =  \t %.4f\n", (double) nacc[6][1]/nacc[6][0], (double) nacc[7][1]/nacc[7][0]);
		fprintf(output, "Change Hotspot Scale =\t %.4f\n",(double) nacc[8][1]/nacc[8][0]);

		fclose(output);
	}

	void restore_rmap_from_cache(int start_loci, int end_loci)
	{
		int i;
		for (i=start_loci; i<=end_loci; i++)
			rmap[i] = rmap_cache[i];
	}

	void calc_rmap(data *mydata, int start_loci, int end_loci)
	{
		unsigned int i, j;
		int start_block;

		// Cache current rmap
		for (i=start_loci; i<=end_loci; i++)
			rmap_cache[i] = rmap[i];

		//start_loci = 1;
		//end_loci = mydata->lseq;

		// Find first block with contribution to this SNP
		start_block = 0;
		for (i=0; i<blocks.size(); i++)
		{
			if ((blocks[i].position >= mydata->locs[start_loci]) || (blocks[i].endposition > mydata->locs[start_loci]))
			{
				start_block = i;
				break;
			}
		}
		// Calculate rmap between limits.
		i = start_loci;
		rmap[i] = rmap[i-1];
		j = start_block;
		double pos1 = 0.0, pos2 = 0.0;
		do
		{ 
			pos1 = maxd(mydata->locs[i], blocks[j].position);
			pos2 = mind(mydata->locs[i+1], blocks[j].endposition);

			if (pos1 != pos2)
			{
				rmap[i] += blocks[j].rate * (pos2 - pos1);		// ****************
				// contribution from hotspots
				rmap[i] += blocks[j].calc_contribution_to_rmap_from_hotspots(pos1, pos2);
			}

			if (blocks[j].endposition == pos2)
			{
				j++;	// Move to next block
			}
			if (mydata->locs[i+1] == pos2)
			{
				i++;	// Move to next SNP
				if ((pos1 < mydata->locs[end_loci+1]) && (pos2 < mydata->locs[end_loci+1]))
					rmap[i] = rmap[i-1];	// If moving onto the next SNP, reset the rmap
			}
		} while ((pos1 < mydata->locs[end_loci+1]) && (pos2 < mydata->locs[end_loci+1]));

		rmap[mydata->lseq] = rmap[mydata->lseq-1];
		//check_rmap(mydata);		// For debugging

	}


	bool check_rmap(data *mydata)
	{
		int i;
		for (i=1; i<=mydata->lseq; i++)
			if ((rmap[i-1] - rmap[i]) > 0.0000001)
			{
				cout << "Error in rmap!" << endl;
				FILE *file1, *file2;
				file1 = fopen("errorblocks.txt", "w");
				file2 = fopen("errorhotspots.txt", "w");
				print_blocks_and_hotspots(file1, file2);
				fclose(file1); fclose(file2);
				file1 = fopen("errorrates.txt", "w");
				print_rates(mydata, file1);
				fclose(file1);
				return false;
			}
		return true;
	}

	void print_blocks_and_hotspots(FILE *blockfile, FILE *hotspotfile)
	{
		unsigned int i;
		for (i=0; i<blocks.size(); i++)
		{
			if (print_blocks) blocks[i].print_block(blockfile, run);
			if (print_hotspots) blocks[i].print_hotspots(hotspotfile, run);
		}
	}

	void print_rates(data *mydata, FILE *outfile)
	{
		int i;
		double rate;
		fprintf(outfile, "%.6f ", rmap[mydata->lseq]);
		for(i=1; i<mydata->lseq; i++)
		{
			rate = (rmap[i] - rmap[i-1]) / (mydata->locs[i+1] - mydata->locs[i]);
			fprintf(outfile, "%.6f ", rate);
		}
		fprintf(outfile, "\n");
	}

	void find_ith_hotspot(unsigned int hot_num, int &container_block, int &hot_num_in_block)
	{
		unsigned int i, j, k=0;
		hot_num_in_block = -1;
		container_block=-1;
		
		for (i=0; i<blocks.size(); i++)	// Find the selected hotspot
		{
			container_block = i;
			for (j = 0; j < blocks[i].hotspots.size(); j++)
			{
				if (k == hot_num)
				{
					hot_num_in_block = j;
					i = (unsigned int)blocks.size();
					break;
				}
				k++;
			}
		}
		if ((container_block == -1) || (hot_num_in_block == -1))
		{
			cout << "\nBlock Positions" << endl;
			for (i=0; i<blocks.size(); i++)
				cout << blocks[i].position << "\t" << blocks[i].endposition << endl;

			cout << "\nHotspot Positions" << endl;
			for (i=0; i<blocks.size(); i++)
				for (j = 0; j < blocks[i].hotspots.size(); j++)
					cout << blocks[i].hotspots[j].position << endl;

			nrerror("Cannot find hotspot");
		}
	}

	// Calculate move probablities
	double calc_ck(int k)
	{
		if (k == max_changepoints)
			return 0.0;
		else
			return mind(1.0, gamma / (k + 1.0));
	}

	double calc_dk(int k)
	{
		if (k == 0)
			return 0.0;
		else
            return mind(1.0, (k) / gamma);
	}

	double calc_h_ck(int k)
	{
		if (k == max_hotspots)
			return 0.0;
		else
			return mind(1.0, theta / (k + 1.0));
	}

	double calc_h_dk(int k)
	{
		if (k == 0.0)
			return 0.0;
		else
			return mind(1.0, k / theta);
	}

	double calc_split_merge_constant()
	{
		int i;
		double c = mind(1.0, mind(1.0/mind(1.0, gamma), 1.0 / mind(1.0, max_changepoints / gamma)));
		for(i=1; i<= max_changepoints-1; i++)
			c = mind(c, 1.0/(mind(1.0, gamma / (i+1.0))+mind(1.0, i / gamma)));
		c = 0.9*c;
		return c;
	}

	double calc_hot_birth_death_constant()
	{
		int i;
		double c = mind(1.0, mind(1.0/mind(1.0, theta), 1.0 / mind(1.0, max_hotspots / theta)));
		for(i=1; i<= max_hotspots-1; i++)
			c = mind(c, 1.0/(mind(1.0, theta / (i+1.0))+mind(1.0, i / theta)));
		c = 0.9*c;
		return c;
	}


	void do_move(data *mydata, likelihood *myLikelihood)
	{
		double r;
		double dk, ck, h_dk, h_ck;

		dk = split_merge_constant * calc_dk(numchangepoints);
		ck = split_merge_constant * calc_ck(numchangepoints);

		h_ck = hot_birth_death_constant * calc_h_ck(numhotspots);
		h_dk = hot_birth_death_constant * calc_h_dk(numhotspots);

		r = ran2();

		if (split_merge_move)
		{
			if (r < dk)
				merge_blocks(mydata, myLikelihood);
			else if (r < (ck + dk))
				split_block(mydata, myLikelihood);
		}

		if (rate_move) change_block_rate(mydata, myLikelihood);
		if (move_block_move) move_block_boundary(mydata, myLikelihood);

		r = ran2();
		if (insert_delete_hot_move)
		{
			if  (r < h_dk)
				delete_hotspot(mydata, myLikelihood);
			else if (r < (h_ck + h_dk))
				insert_hotspot(mydata, myLikelihood);
		}


		if (heat_move) change_hotspot_heat(mydata, myLikelihood);
		if (move_hot_move) move_hotspot(mydata, myLikelihood);
		if (scale_move) change_hotspot_scale(mydata, myLikelihood);

	}

	int change_block_rate(data *mydata, likelihood *myLikelihood)
	{
		nacc[0][0]++;
		int i;
		double orig_rate;

		i=(int)(((double)numchangepoints+1)*ran2());
		orig_rate = blocks[i].rate;
		// Change block rate
		double u = ran2() - 0.5;
		blocks[i].rate = blocks[i].rate * exp(u);

		if (blocks[i].rate < blockparam.RMIN)
		{
			blocks[i].rate = orig_rate;
			return 0;
		}
		// Calc recombination map
		calc_rmap(mydata, blocks[i].start_loci, mydata->lseq);
		
		// Calc lk following change to this block (need to calc until the end of the data)
		lk[1] = myLikelihood->lk_calc(blocks[i].start_loci, blocks[i].end_loci, mydata, rmap, 1);
		
		double log_alpha = (orig_rate - blocks[i].rate)/blockparam.RATE_BETA;
		log_alpha += blockparam.RATE_ALPHA*log(blocks[i].rate / orig_rate);		// Proposal density

		if ((log(ran2()) > (log_alpha + lk[1])))
		{	// Reject
			blocks[i].rate = orig_rate;
			restore_rmap_from_cache(blocks[i].start_loci, mydata->lseq);
		}
		else 
		{	// Accept
			lk[0] += lk[1]; 
			myLikelihood->update_lij(blocks[i].start_loci, blocks[i].end_loci, mydata, 1);
			nacc[0][1]++; 
		}
		return 1;
	}

	int move_block_boundary(data *mydata, likelihood *myLikelihood)
	{
		if (numchangepoints == 0)
			return 0;
		nacc[1][0]++;
		int i = ((int) (((double) numchangepoints)*ran2())+1);	//Choose block for updating

		unsigned int j,k;
		
		double orpos = blocks[i].position;	// Store the original position
			
		double leftlim = blocks[i-1].position;
		double rightlim = blocks[i].endposition;

		// choose a new position
		blocks[i].position = leftlim + ran2() * (rightlim - leftlim);
		blocks[i-1].endposition = blocks[i].position;

		if ((blocks[i].position <= leftlim + 0.001) || (blocks[i].position >= rightlim - 0.001))
		{	// Can't move boundary outside limits...
			blocks[i].position = orpos;
			blocks[i-1].endposition = blocks[i].position;
			return 0;
		}

		int start_loci = blocks[i-1].start_loci;

		// See if move has moved across SNP boundary...
		int orig_loci = blocks[i].start_loci;
		for (j = start_loci; j <= blocks[i].end_loci; j++)
			if ((mydata->locs[j] <= blocks[i].position) &&
				(mydata->locs[j+1] > blocks[i].position))
			{
				blocks[i].start_loci = j;
				blocks[i-1].end_loci = j;
				break;
			}
		
		// Calculate total recomb over region
		double orig_rate[2];
		orig_rate[0] = blocks[i-1].rate;
		orig_rate[1] = blocks[i].rate;			// Cache old rates
		
		double totalrecomb = (orig_rate[0] * (orpos - leftlim)) + (orig_rate[1] * (rightlim - orpos));
		
		vector <hotspot> local_hotspots;		// Store hotspots contained in blocks...
		for(k=0; k<2; k++)
			for(j=0; j<blocks[i-1+k].hotspots.size(); j++)
			{
				local_hotspots.push_back(blocks[i-1+k].hotspots[j]);
			}

		// Calculate new rate so that total recombination over region is unchanged.
		// AT THIS POINT WE NEED TO UPDATE THE LEFT OR RIGHT BLOCK WITH EQUAL PROB!
		double u = ran2();
		double temp, temp2;
		if (u > 0.5)
		{  //Update rh side
			k = 1;
			temp = totalrecomb - (orig_rate[0] * (blocks[i].position - leftlim));
			temp2 = rightlim - blocks[i].position;
		}
		else	// Update lh side
		{
			k = 0;
			temp = totalrecomb - (orig_rate[1] * (rightlim - blocks[i].position));
			temp2 = blocks[i].position - leftlim;
		}

		blocks[i-1+k].rate = temp / temp2;		// Set new rate
		
		if ((blocks[i].rate < blockparam.RMIN) || (blocks[i-1].rate < blockparam.RMIN))
		{	// Can't have a negative rate, orig_rate one over the max allowed, so undo the move
			blocks[i].position = orpos;
			blocks[i].rate = orig_rate[1]; 
			blocks[i-1].rate = orig_rate[0]; 
			blocks[i-1].endposition = blocks[i].position;
			blocks[i].start_loci = orig_loci;
			blocks[i-1].end_loci = orig_loci;
			return 0;
		}

		double log_alpha = (orig_rate[k] - blocks[i-1+k].rate) / blockparam.RATE_BETA;
		log_alpha += (blockparam.RATE_ALPHA-1.0)*log(blocks[i].rate / orig_rate[k]);		// Proposal density
		double s_star;
		s_star = blocks[i].position;
		log_alpha += log((rightlim - s_star) * (s_star - leftlim) / ((rightlim - orpos) * (orpos - leftlim)));

		blocks[i-1].find_hotspots_in_block(local_hotspots);
		blocks[i].find_hotspots_in_block(local_hotspots);

		//double totalrecombafter = (blocks[i-1].rate * (blocks[i-1].endposition - leftlim)) + (blocks[i].rate * (rightlim - blocks[i].position));

		// Only need to calc between the two block boundaries
		calc_rmap(mydata, blocks[i-1].start_loci, blocks[i].end_loci);
		lk[1] = myLikelihood->lk_calc(blocks[i-1].start_loci, blocks[i].end_loci, mydata, rmap, 0);

		bool reject = false;
		if (one_hotspot_per_block && ((blocks[i-1].hotspots.size() > 1) || (blocks[i].hotspots.size() > 1)))
			reject = true;
        
		if (reject || (log(ran2()) > (log_alpha + lk[1])))
		{	// Reject
			blocks[i].position = orpos; 
			blocks[i].rate = orig_rate[1]; 
			blocks[i-1].rate = orig_rate[0]; 
			blocks[i-1].endposition = blocks[i].position;
			blocks[i-1].find_hotspots_in_block(local_hotspots);	blocks[i].find_hotspots_in_block(local_hotspots);
			blocks[i].start_loci = orig_loci;
			blocks[i-1].end_loci = orig_loci;
			restore_rmap_from_cache(blocks[i-1].start_loci, blocks[i].end_loci);
		}
		else 
		{	// Accept
			lk[0] += lk[1]; 
			myLikelihood->update_lij(blocks[i-1].start_loci, blocks[i].end_loci, mydata, 0);
			nacc[1][1]++; 
		}
		return 1;
	}

	int split_block(data *mydata, likelihood *myLikelihood)
	{
		unsigned int j;
		bool reject = false;
		nacc[2][0]++;

		int i;
		//Choose block for updating
		double new_position = mydata->locs[1] + ran2() * (mydata->tlseq);	// New changepoint location
		for (i=0; i<=numchangepoints; i++)
			if ((new_position >= blocks[i].position) && (new_position < blocks[i].endposition))
				break;

		// Now do the sums before we change anything
		double w1 = new_position - blocks[i].position;
		double w2 = blocks[i].endposition - new_position;
		double h = blocks[i].rate;
		double u = ran2();
		double h1 = (u * h * (w1 + w2)) / (u*w1 + (1-u)*w2);
		double h2 = h1 * ((1.0 - u)/u);
		double a1 = blockparam.RATE_ALPHA;
		double leng = mydata->tlseq;
		double b1 = blockparam.RATE_BETA;
		int k = numchangepoints;

		double log_ratio = 0.0;
		log_ratio += log(gamma / (k+1.0));	//log(pk[k+1] / pk[k]);
		log_ratio += log((2.0*k+2.0)*(2.0*k+3.0)/leng/leng);
		log_ratio += log(w1*w2/(w1+w2));
		log_ratio += log_block_prior_const + (a1 - 1.0)*log(h1*h2/h);
		log_ratio -= (h1+h2-h) / b1;
		// Proposal
		log_ratio += log(calc_dk(k+1)*leng / (calc_ck(k) * (k+1.0)));
		// Jacobian
		log_ratio += log(h1*h2/(h*u*(1.0-u)));

		int start_loci = blocks[i].start_loci;

		// Create new block
		for (j = start_loci; j <= blocks[i].end_loci; j++)
			if ((mydata->locs[j] <= new_position) &&
				(mydata->locs[j+1] > new_position))
			{
				start_loci = j;
				break;
			}

		block new_block;
		new_block.init(&blockparam, new_position, blocks[i].endposition, start_loci, blocks[i].end_loci);
		blocks.insert(blocks.begin() + i + 1, new_block);	// Insert the block into the block list
		blocks[i].endposition = new_block.position;
		blocks[i].end_loci = new_block.start_loci;

		// Calculations required to figure out the new rates....
		vector<hotspot> local_hotspots;
		for (j=0; j<blocks[i].hotspots.size(); j++)
		{
			local_hotspots.push_back(blocks[i].hotspots[j]);
		}

		blocks[i+1].find_hotspots_in_block(local_hotspots);
		blocks[i].find_hotspots_in_block(local_hotspots);

		blocks[i].rate = h1;
		blocks[i+1].rate = h2;

		// Reject the move if it gives us a negative rate, orig_rate the rate is above the max
		double alpha = 0.0;
		calc_rmap(mydata, blocks[i].start_loci, blocks[i+1].end_loci);
		if ((blocks[i].rate < blockparam.RMIN) || (blocks[i+1].rate < blockparam.RMIN))
		{
			reject = true;
		}
		else
		{
			lk[1] = myLikelihood->lk_calc(blocks[i].start_loci, blocks[i+1].end_loci, mydata, rmap, 0);

			alpha = exp(maxd(-30.0, mind(0.0, log_ratio + lk[1])));
		}

		if (one_hotspot_per_block && ((blocks[i].hotspots.size() > 1) || (blocks[i+1].hotspots.size() > 1)))
			reject = true;
	
		if (reject || (ran2() > alpha))
		{	// REJECT
			blocks[i].rate = h;
			blocks[i].endposition = new_block.endposition;
			blocks[i].end_loci = new_block.end_loci;
			blocks[i].find_hotspots_in_block(local_hotspots);		//Reset hotspot pointers
			blocks.erase(blocks.begin() + i + 1);
			restore_rmap_from_cache(blocks[i].start_loci,  blocks[i].end_loci);
		}
		else 
		{ //ACCEPT
			numchangepoints++;
			nacc[2][1]++;
			lk[0] += lk[1]; 
			myLikelihood->update_lij(blocks[i].start_loci, blocks[i+1].end_loci, mydata, 0);
		}


		return 1;
	}

	int merge_blocks(data *mydata, likelihood *myLikelihood)
	{
		if (numchangepoints == 0)
			return 0;	// Can't merge if no changepoints
		unsigned int j;
		nacc[3][0]++;

		int i = (int) (((double) numchangepoints)*ran2());	// Choose changepoint to remove

		double w1 = blocks[i].endposition - blocks[i].position;
		double w2 = blocks[i+1].endposition - blocks[i+1].position;
		double h1 = blocks[i].rate;
		double h2 = blocks[i+1].rate;
		double h = (h1*w1 + h2*w2) / (w1 + w2);
		double u = 1.0/(1.0+(h2/h1));

		double a1 = blockparam.RATE_ALPHA;
		double leng = mydata->tlseq;
		double b1 = blockparam.RATE_BETA;
		int k = numchangepoints;

		double log_ratio = 0.0;
		// Prior
		log_ratio += log( gamma / double(k)); //log(pk[k] / pk[k-1]);
		log_ratio += log((2.0*k)*(2.0*k+1.0)/leng/leng);
		log_ratio += log(w1*w2/(w1+w2));
		log_ratio += log_block_prior_const + (a1 - 1.0)*log(h1*h2/h);
		log_ratio -= (h1+h2-h) / b1;
		// Proposal
		log_ratio += log(calc_dk(k)*leng / (calc_ck(k-1) * k));
		// Jacobian
		log_ratio += log(h1*h2/(h*u*(1-u)));

		log_ratio = -log_ratio;

		block blockcache = blocks[i+1]; // Keep a copy of the block being removed
		int orig_end_loci = blocks[i].end_loci;

		// Set hotspot pointers
		vector<hotspot> local_hotspots;
		for(j=0; j<blocks[i].hotspots.size(); j++)
			local_hotspots.push_back(blocks[i].hotspots[j]);
		for(j=0; j<blockcache.hotspots.size(); j++)
			local_hotspots.push_back(blockcache.hotspots[j]);

		blocks[i].endposition = blockcache.endposition;
		blocks[i].end_loci = blockcache.end_loci;

		blocks[i].rate = h;

		double log_prior = -HUGE_VAL, log_hr = -HUGE_VAL, log_jcb = -HUGE_VAL;
		double alpha = 0.0;
		bool reject = false; 
			
		blocks[i].hotspots = local_hotspots;	// The block now contains all the hotspots.
		blocks.erase(blocks.begin() + i + 1);	// Remove the next block

		calc_rmap(mydata, blocks[i].start_loci, blocks[i].end_loci);
		lk[1] = myLikelihood->lk_calc(blocks[i].start_loci, blocks[i].end_loci, mydata, rmap, 0);
		alpha = exp(maxd(-30.0, mind(0.0, log_ratio + lk[1])));
		
		if (blocks[i].rate < blockparam.RMIN)
			reject = true;	// Should never happen

		//Acceptance on likelihood, rates and block penalty
		if (reject || (ran2() > alpha))
		{	// REJECT
			blocks[i].rate=h1;
			blocks[i].endposition = blockcache.position;
			blocks[i].end_loci = orig_end_loci;
			blocks[i].find_hotspots_in_block(local_hotspots);
			blocks.insert(blocks.begin() + i + 1, blockcache);
			restore_rmap_from_cache(blocks[i].start_loci, blocks[i+1].end_loci);
		} 
		else 
		{	//ACCEPT
			numchangepoints--;
			nacc[3][1]++;
			lk[0] += lk[1];
			myLikelihood->update_lij(blocks[i].start_loci, blocks[i].end_loci, mydata, 0);
		}
		return 1;
	}

	int change_hotspot_heat(data *mydata, likelihood *myLikelihood)
	{
		if (numhotspots == 0)
			return 0;	// Can't change the heat if there are no hotspots!

		nacc[4][0]++;

		int hot_num = (int)(((double)numhotspots) * ran2());
		int block_hot_num = -1;
		int container = -1;

		find_ith_hotspot(hot_num, container, block_hot_num);
		double orig_lamda = blocks[container].hotspots[block_hot_num].lamda;
		double orig_loglamda = blocks[container].hotspots[block_hot_num].loglamda;

		double u = ran2() - 0.5;
		blocks[container].hotspots[block_hot_num].loglamda += u;
		blocks[container].hotspots[block_hot_num].lamda = exp(blocks[container].hotspots[block_hot_num].loglamda);

		if ((blocks[container].hotspots[block_hot_num].lamda <= 0.0) || (blocks[container].hotspots[block_hot_num].lamda >= HUGE_VAL))
		{
			blocks[container].hotspots[block_hot_num].lamda = orig_lamda;
			blocks[container].hotspots[block_hot_num].loglamda = orig_loglamda;
			return 0;
		}
		
		calc_rmap(mydata, blocks[container].start_loci, mydata->lseq);
	
		// Calc lk following change to this hotspot
		lk[1] = myLikelihood->lk_calc(blocks[container].start_loci, blocks[container].end_loci, mydata, rmap, 1);
	
		// Prior on hotspot change
		double log_alpha = (hotspotparam.HEAT_ALPHA) * log(blocks[container].hotspots[block_hot_num].lamda / orig_lamda);
		log_alpha += ((orig_lamda - blocks[container].hotspots[block_hot_num].lamda) / hotspotparam.HEAT_BETA);

		if ((log(ran2()) > (log_alpha + lk[1])))
		{	// Reject
			blocks[container].hotspots[block_hot_num].lamda = orig_lamda;
			blocks[container].hotspots[block_hot_num].loglamda = orig_loglamda;
			restore_rmap_from_cache(blocks[container].start_loci, mydata->lseq);
		}
		else 
		{	// Accept
			lk[0] += lk[1];
			myLikelihood->update_lij(blocks[container].start_loci, blocks[container].end_loci, mydata, 1);
			nacc[4][1]++; 
		}	
		return 1;
	}

	int move_hotspot(data *mydata, likelihood *myLikelihood)
	{
		if (numhotspots == 0)
			return 0;
		nacc[5][0]++;

		bool reject=false;
		int i;
		int hot_num = (int)(((double)numhotspots) * ran2());
		int block_hot_num = -1;
		int container = -1;

		find_ith_hotspot(hot_num, container, block_hot_num);

		hotspot orig_hotspot = blocks[container].hotspots[block_hot_num];				// Store the original hotspot
		double orig_pos = blocks[container].hotspots[block_hot_num].position;			// Store the original position

		double movedistance = hotspotparam.MOVE_SIGMA * normrnd();
		blocks[container].hotspots[block_hot_num].position += movedistance;				// Set a new position

		// Cannot move the hotspot out of the data range
		if ((blocks[container].hotspots[block_hot_num].position <= mydata->locs[1]) || 
			(blocks[container].hotspots[block_hot_num].position >= mydata->locs[mydata->lseq]))
		{
			blocks[container].hotspots[block_hot_num].position = orig_pos;
			return 0;	// reject move
		}

		blocks[container].hotspots[block_hot_num].leftlim += movedistance;
		blocks[container].hotspots[block_hot_num].rightlim += movedistance;

		int new_container = container, new_hot_num = block_hot_num;
		// Has the hotspot moved into a new block?
		if (blocks[container].hotspots[block_hot_num].position < blocks[container].position)
		{	// Moved left of current block
			for (i=0; i < container; i++)
			{
				if ((blocks[container].hotspots[block_hot_num].position >= blocks[i].position) && 
					(blocks[container].hotspots[block_hot_num].position < blocks[i].endposition))
				{
					hotspot new_hotspot = blocks[container].hotspots[block_hot_num];	// Copy the hotspot
					if (one_hotspot_per_block && (blocks[i].hotspots.size() > 0))
						reject = true;	// Reject
					blocks[i].hotspots.push_back(new_hotspot);							// Insert hotspot into new block
					new_container = i;
					new_hot_num = (int)blocks[i].hotspots.size() - 1;
					blocks[container].hotspots.erase(blocks[container].hotspots.begin() + block_hot_num);	// Erase the original hotspot
					break;
				}
			}
		}
		else if (blocks[container].hotspots[block_hot_num].position >= blocks[container].endposition)
		{	// Moved Right of current block
			for (i = container; i < (int)blocks.size(); i++)
			{
				if ((blocks[container].hotspots[block_hot_num].position >= blocks[i].position) && 
					(blocks[container].hotspots[block_hot_num].position < blocks[i].endposition))
				{
					hotspot new_hotspot = blocks[container].hotspots[block_hot_num];	// Copy the hotspot
					if (one_hotspot_per_block && (blocks[i].hotspots.size() > 0))
						reject = true;	// Reject
					blocks[i].hotspots.push_back(new_hotspot);							// Insert hotspot into new block
					new_container = i;
					new_hot_num = (int)blocks[i].hotspots.size() - 1;
					blocks[container].hotspots.erase(blocks[container].hotspots.begin() + block_hot_num);	// Erase the original hotspot
					break;
				}
			}
		}
		
		int left_loci, right_loci; 
		if (new_container == container)
		{
			left_loci = blocks[container].start_loci; right_loci = blocks[container].end_loci;
		}
		else if (new_container < container)
		{
			left_loci = blocks[new_container].start_loci; right_loci = blocks[container].end_loci;
		}
		else
		{
			left_loci = blocks[container].start_loci; right_loci = blocks[new_container].end_loci;
		}

		calc_rmap(mydata, left_loci, right_loci);
		if (!reject)
		{
			lk[1] = myLikelihood->lk_calc(left_loci, right_loci, mydata, rmap, 0);
		}

		// No change in prior probabilities!
		if (reject || (log(ran2()) > lk[1]))
		{	// Reject
			blocks[new_container].hotspots.erase(blocks[new_container].hotspots.begin() + new_hot_num);;
			blocks[container].hotspots.push_back(orig_hotspot);
			restore_rmap_from_cache(left_loci, right_loci);
		}
		else 
		{	// Accept
			lk[0] += lk[1];
			myLikelihood->update_lij(left_loci, right_loci, mydata, 0);
			nacc[5][1]++; 
		}
		return 1;
	}

	int change_hotspot_scale(data *mydata, likelihood *myLikelihood)
	{
		if (numhotspots == 0)
			return 0;
		nacc[8][0]++;

		int hot_num = (int)(((double)numhotspots) * ran2());
		int block_hot_num = -1;
		int container = -1;

		find_ith_hotspot(hot_num, container, block_hot_num);

		double u = ran2()-0.5;

		double orig_mu = blocks[container].hotspots[block_hot_num].mu;
		blocks[container].hotspots[block_hot_num].mu *= exp(u);
		
		if (blocks[container].hotspots[block_hot_num].mu <= 0.0)
		{
			blocks[container].hotspots[block_hot_num].mu = orig_mu;
			return 0;
		}
		
		calc_rmap(mydata, blocks[container].start_loci, blocks[container].end_loci);
		// Calc lk following change to this hotspot
		lk[1] = myLikelihood->lk_calc(blocks[container].start_loci, blocks[container].end_loci, mydata, rmap, 0);

		double log_alpha = (hotspotparam.MU_ALPHA) * log(blocks[container].hotspots[block_hot_num].mu / orig_mu);
		log_alpha += ((orig_mu - blocks[container].hotspots[block_hot_num].mu) / hotspotparam.MU_BETA);
		
		if ((log(ran2()) > (log_alpha + lk[1])))
		{	// Reject
			blocks[container].hotspots[block_hot_num].mu = orig_mu;
			restore_rmap_from_cache(blocks[container].start_loci, blocks[container].end_loci);
		}
		else 
		{	// Accept
			lk[0] += lk[1];
			myLikelihood->update_lij(blocks[container].start_loci, blocks[container].end_loci, mydata, 0);
			nacc[8][1]++;
		}	
		return 1;
	}

	int insert_hotspot(data *mydata, likelihood *myLikelihood)
	{
		nacc[6][0]++;
		unsigned int i;
		int container_block=-1;
		bool reject = false;
	
		hotspot new_hotspot;
	
		double position = mydata->locs[1] + ran2() * (mydata->tlseq);
		new_hotspot.init(&hotspotparam, position);

		// Find block that contains hotspot.
		for (i=0; i<blocks.size(); i++)
		{
			if ((position >= blocks[i].position) && (position < blocks[i].endposition))
			{
				if (one_hotspot_per_block && (blocks[i].hotspots.size() > 0))
					reject = true;	// Reject
				blocks[i].hotspots.push_back(new_hotspot);
				container_block = i;
				break;
			}
		}

		numhotspots++;

		calc_rmap(mydata, blocks[container_block].start_loci, mydata->lseq);

		if (!reject)
		{
			lk[1] = myLikelihood->lk_calc(blocks[container_block].start_loci, blocks[container_block].end_loci, mydata, rmap, 1);
		}
		
		// Log version (more accurate)
		double log_alpha=0.0;
		log_alpha = (log_theta  - log((double)numhotspots));
		log_alpha += log(calc_h_dk(numhotspots + 1)/calc_h_ck(numhotspots));

		if (reject || (log(ran2()) > (log_alpha + lk[1])))
		{	// Reject
			blocks[container_block].hotspots.pop_back();
			numhotspots--;
			restore_rmap_from_cache(blocks[container_block].start_loci, mydata->lseq);
		}
		else 
		{	// Accept
			lk[0] += lk[1];
			myLikelihood->update_lij(blocks[container_block].start_loci, blocks[container_block].end_loci, mydata, 1);
			nacc[6][1]++; 
		}
		return 1;
	}

	int delete_hotspot(data *mydata, likelihood *myLikelihood)
	{
		if (numhotspots == 0)
			return 0;
		nacc[7][0]++;

		int hot_num = (int)(((double)numhotspots) * ran2());
		int deleted_hotspot=-1;
		int container_block=-1;

		find_ith_hotspot(hot_num, container_block, deleted_hotspot);

		double orig_lamda = blocks[container_block].hotspots[deleted_hotspot].lamda;

		// Just set the lamda of the hotspot to zero, and free it if we accept the move....
		blocks[container_block].hotspots[deleted_hotspot].lamda = 0.0;

		calc_rmap(mydata, blocks[container_block].start_loci, mydata->lseq);
		lk[1] = myLikelihood->lk_calc(blocks[container_block].start_loci, blocks[container_block].end_loci, mydata, rmap, 1);

		// Log version

		double log_alpha=0.0;
		log_alpha = (log((double)numhotspots) - log_theta);

		log_alpha += log(calc_h_ck(numhotspots-1)/calc_h_dk(numhotspots));

		//log_alpha += log((numhotspots));


		if ((log(ran2()) > (log_alpha + lk[1])))
		{	// Reject
			blocks[container_block].hotspots[deleted_hotspot].lamda = orig_lamda;
			restore_rmap_from_cache(blocks[container_block].start_loci, mydata->lseq);
		}
		else 
		{	// Accept
			numhotspots--;
			blocks[container_block].hotspots.erase(blocks[container_block].hotspots.begin() + deleted_hotspot);
			lk[0] += lk[1];
			myLikelihood->update_lij(blocks[container_block].start_loci, blocks[container_block].end_loci, mydata, 1);
			nacc[7][1]++; 
		}
		return 1;
	}

};

#endif

