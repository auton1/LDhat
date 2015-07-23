#if !defined BLOCK_H
#define BLOCK_H
#pragma warning(disable:4786) 

#include <iostream>
#include <vector>
#include <math.h>

using namespace std;
#include "params.h"
#include "rhomap_tools.h"
#include "hotspot.h"

class block
{
public:
	double rate;

	unsigned int start_loci;
	unsigned int end_loci;
	double position;
	double endposition;

	vector<hotspot> hotspots;	// Vector of hotspots contained in block

	block()
	{
	}

	~block()
	{
	}

	block copy()
	{
		unsigned int i;
		block new_block;
		new_block.init(rate, position, endposition, start_loci, end_loci);
		for(i=0; i<hotspots.size(); i++)
			new_block.hotspots.push_back(hotspots[i].copy());

		return new_block;
	}

	void init(const struct block_params *myblockparams, double pos, double endpos, int startloci, int endloci)
	{
		rate = gamrnd(myblockparams->RATE_ALPHA, myblockparams->RATE_BETA);
		position = pos;
		endposition = endpos;
		start_loci = startloci;
		end_loci = endloci;
	}

	void init(double Rate, double pos, double endpos, int startloci, int endloci)
	{
		rate = Rate;
		position = pos;
		endposition = endpos;
		start_loci = startloci;
		end_loci = endloci;
	}

	double calc_contribution_to_rmap_from_hotspots(double leftlim, double rightlim)
	{
		unsigned int i;
		double recomb=0.0;
		for(i=0; i<hotspots.size(); i++)
			recomb += hotspots[i].calc_contribution_to_rmap(position, endposition, leftlim, rightlim);
		return recomb;
	}

	void print_block(FILE *outfile, int run)
	{
		fprintf(outfile, "%i %f %f\n", run, position, rate);
	}

	void print_hotspots(FILE *outfile, int run)
	{
		unsigned int i;
		for(i=0; i<hotspots.size(); i++)
			hotspots[i].print_hotspot(outfile, run, rate);
	}

	void print_hotspots_to_screen(int run)
	{
		unsigned int i;
		for(i=0; i<hotspots.size(); i++)
			hotspots[i].print_hotspot_to_screen(run, rate);
	}

	void find_hotspots_in_block(vector<hotspot> hotspot_list)
	{
		hotspots.clear();
		unsigned int i;
		for(i=0; i<hotspot_list.size(); i++)
			if ((hotspot_list[i].position < endposition) && (hotspot_list[i].position >= position))
				hotspots.push_back(hotspot_list[i]);
	}

};

#endif
