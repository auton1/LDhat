#if !defined HOTSPOT_H
#define HOTSPOT_H
#pragma warning(disable:4786) 
#include <vector>
#include <math.h>
using namespace std;

#include "rhomap_tools.h"

class hotspot
{
public:
	double lamda;		// Hotspot mass
	double loglamda;	// Log of hotspot mass
	double mu;			// Hotspot scale;
	double position;	// Hotspot position

	double leftlim;		// Absolute Mimimum limit of hotspot
	double rightlim;	// Absolute Maximum limit of hotspot

	hotspot()
	{
	}
	
	~hotspot()
	{
	}


	void init(const struct hotspot_params *myhotspotparams, double pos)
	{
		position = pos;
		leftlim = pos - myhotspotparams->MAX_WIDTH;
		rightlim = pos + myhotspotparams->MAX_WIDTH;

		mu = gamrnd(myhotspotparams->MU_ALPHA, myhotspotparams->MU_BETA);

		//lamda = myhotspotparams->minLamda + gamrnd(myhotspotparams->HEAT_ALPHA, myhotspotparams->HEAT_BETA);
		lamda = 0.0 + gamrnd(myhotspotparams->HEAT_ALPHA, myhotspotparams->HEAT_BETA);
		loglamda = log(lamda);
	}

	void init(double Lamda, double Mu, double pos, double left_lim, double right_lim)
	{
		position = pos;
		leftlim = left_lim;
		rightlim = right_lim;
		mu = Mu;
		lamda = Lamda;
		loglamda = log(lamda);
	}

	hotspot copy()
	{
		hotspot new_hotspot;
		new_hotspot.init(lamda, mu, position, leftlim, rightlim);
		return new_hotspot;
	}

	double calc_contribution_to_rmap(double block_pos, double block_end, double loci_pos, double loci_end)
	{
		double recomb = 0;
		double mass = lamda;

		double min1 = maxd(leftlim, block_pos);		// Find minimum limit of hotspot
		double max1 = mind(rightlim, block_end);	// Find maximum limit of hotspot

		// Case 0: Loci block has no contribution from hotspot
		if ((min1 >= loci_end) || (max1 <= loci_pos))
			return recomb;
		
		// Case 1: Loci block entirely contains hotspot
		if ((min1 >= loci_pos) && (max1 <= loci_end))
		{
			recomb = mass;
			return recomb;
		}

		double t = position;
		// Calculate increase in rate due to trunctated tails
		double left_tail_mass = 0.5 * mass * exp((min1 - t) / mu);
		double left_tail_rate = left_tail_mass / (t - min1);
		double right_tail_mass = mass * 0.5 * exp((t - max1) / mu);
		double right_tail_rate = right_tail_mass / (max1 - t);

		// Case 2: Loci block contains minimum limit and peak, but not maximum limit
		if ((min1 >= loci_pos) && (t <= loci_end) && (max1 >= loci_end))
		{	// Contribution from left tail contained in body
			recomb = right_tail_rate * (loci_end - t);					// Contribution from right tail
			recomb += mass * (1.0 - (0.5 * exp((t - loci_end) / mu)));	// Contribution from body
			return recomb;
		}

		// Case 3: Block contains peak and maximium limit, but not minimum limit
		if ((min1 <= loci_pos) && (t >= loci_pos) && (t <= loci_end) && (max1 <= loci_end))
		{
			recomb = left_tail_rate * (t - loci_pos);					// Contribution from left tail
			recomb += mass * (1.0 - (0.5 * exp((t - max1) / mu)));		// Contribution from body
			recomb -= 0.5 * mass * exp((loci_pos - t) / mu);			// Contribution from body
			recomb += right_tail_mass;									// Contribution from right tail
			return recomb;
		}

		// Case 4: Block contains mimimum limit but not peak
		if ((min1 >= loci_pos) && (min1 <= loci_end) && (t >= loci_pos) && (t >= loci_end) && (max1 >= loci_end))
		{	
			recomb = 0.5 * mass * exp((loci_end - t) / mu);				// Contribution from body
			recomb -= left_tail_mass;
			recomb += left_tail_rate * (loci_end - min1);
			return recomb;
		}

		// Case 5: Block left of peak, and contained within limits
		if ((min1 <= loci_pos) && (t >= loci_end) && (max1 >= loci_end))
		{
			recomb = left_tail_rate * (loci_end - loci_pos);			// Contribution from left tail
			recomb += (0.5 * mass * exp((loci_end - t) / mu));			// Contribution from body
			recomb -= (0.5 * mass * exp((loci_pos - t) / mu));			// Contribution from body
			return recomb;
		}

		// Case 6: Block contains peak, and contained within limits
		if ((min1 <= loci_pos) && (t >= loci_pos) && (t <= loci_end) && (max1 >= loci_end))
		{
			recomb = left_tail_rate * (t - loci_pos);					// Contribution from left tail
			recomb += right_tail_rate * (loci_end - t);					// Contribution from right tail
			recomb += mass * (1.0 - (0.5 * exp((t - loci_end) / mu)));	// Contribution from body
			recomb -= (0.5 * mass * exp((loci_pos - t) / mu));			// Contribution from body
			return recomb;
		}

		// Case 7: Block right of peak, and contained with limits
		if ((min1 <= loci_pos) && (t <= loci_pos) && (max1 >= loci_end))
		{
			recomb += right_tail_rate * (loci_end - loci_pos);			// Contribution from right tail
			recomb += mass * (1.0 - (0.5 * exp((t - loci_end) / mu)));	// Contribution from body
			recomb -= mass * (1.0 - (0.5 * exp((t - loci_pos) / mu)));	// Contribution from body
			return recomb;
		}

		// Case 8; Block right of peak, and contains maximum limit
		if ((min1 <= loci_pos) && (t <= loci_pos) && (max1 <= loci_end))
		{
			recomb = right_tail_rate * (max1 - loci_pos);				// Contribution from right tail
			recomb += mass * (1.0 - (0.5 * exp((t - max1) / mu)));		// Contribution from body
			recomb -= mass * (1.0 - (0.5 * exp((t - loci_pos) / mu)));	// Contribution from body
			return recomb;
		}

		cout << "problem here?" << endl;
		return recomb;
	}

	void print_hotspot(FILE *outfile, int run, double bg_rate)
	{
		fprintf(outfile, "%i %f %f %f %f\n", run, position, lamda, mu, bg_rate);
	}

	void print_hotspot_to_screen(int run, double bg_rate)
	{
		printf("%i %f %f %f %f\n", run, position, lamda, mu, bg_rate);
	}

};

#endif
