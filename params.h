#if !defined PARAMS_H
#define PARAMS_H

struct block_params
{
	double RMIN;				// Minimum recombination rate
	double RATE_ALPHA;			// Background rate hyperprior parameter
	double RATE_BETA;			// Background rate hyperprior parameter
};

struct hotspot_params
{
	double MAX_WIDTH;			// Maximum allowed width of hotspot
	double MOVE_SIGMA;			// Hotspot move sigma
	double MU_ALPHA;			// Scale prior hyperparameter
	double MU_BETA;				// Scale prior hyperparameter
	double HEAT_ALPHA;			// Heat prior hyperparameter
	double HEAT_BETA;			// Heat prior hyperparameter
	double minLamda;			// Minimum allowed value of lamda

	double AVG_DIST_BETWEEN_HOTSPOTS;	// Mean distance between hotspots
};

void set_default_parameters(block_params &blockparam, hotspot_params &hotspotparam)
{
		blockparam.RATE_ALPHA = 1.0;
		blockparam.RATE_BETA = 0.05;
		blockparam.RMIN = 0.00001;

		hotspotparam.MOVE_SIGMA = 0.5;

		hotspotparam.MU_ALPHA = 23.98429474371825;
		hotspotparam.MU_BETA = 0.01048804333547;	
		
		hotspotparam.MAX_WIDTH = 5.0;

		hotspotparam.HEAT_ALPHA = 0.61247605196969;
		hotspotparam.HEAT_BETA = 52.47027753716133;
		hotspotparam.minLamda = 0.00001;

		hotspotparam.AVG_DIST_BETWEEN_HOTSPOTS = 50.0;
}

#endif

