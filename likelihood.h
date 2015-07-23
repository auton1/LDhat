#if !defined LIKELIHOOD_H
#define LIKELIHOOD_H
#pragma warning(disable:4786) 
#include <string>
#include <math.h>

using namespace std;

#include "data.h"
#include "compLK.h"
#include "rhomap_tools.h"

class likelihood
{
public:
	string lk_filename;				// Likelihood filename

	string lk_out_filename;			// Output file of likelihood curve
	bool output_lk_curve;			// Output likelihood curve

	double *lnfac_array;			// Precomputed table of log factorials

	int **pij;						// Matrix of pair types at each loci
	int size_pij;					// Size of above matrix
	int npt;						// Number of pair types in data
	int pnew;
	int npmax;						// Maximum number of site types		
	int tcat;						// Number of mutation rate categories (should be 1?)
	int miss;
	double **lkmat;					// Precalculated likelihood values [site type][rho]
	struct site_type **pset;		// Array containing the site type at each position
	double **lij;					// Composite likelihood of data [position][distance from position] (contains a 'temp' and 'non-temp' version)
	int MAXW;						// MAXW*2 = Max number of SNPs to consider for likelihood - i.e. ignore SNPS > MAXW apart
	bool init_called;				// Boolean to record if initialisation routine called.
	double *snp_weights;			// Attempt to compensate for non-independence by weighting likelihood between SNPs

	bool old_lk;					// Use original composite likelihood
	bool no_lk;						// Use original composite likelihood
	double rmax;					// Max rho in coalescent likelihood estimation
	int rcat;						// Number of categories for estimating rho - can be >> rme
	double theta_per_site;			// Theta per site
	int window;						// Number of SNPs to be summed when calculating likelihood

	likelihood()
	{
		lk_filename = "";
		lk_out_filename = "likelihood.txt";
		npt=0; pnew=0;
		tcat = 1;
		miss = 0;
		MAXW = 50;
		init_called = false;
		output_lk_curve = false;
		old_lk = false;
		no_lk = false;
	}

	~likelihood()
	{
		if (init_called)
		{
			free_dvector(lnfac_array);
			free_imatrix(pij, size_pij);
			free_dmatrix(lij, size_pij);
			free_dvector(snp_weights);
			free_dmatrix(lkmat, npt+pnew+miss+1);
			delete pset;
		}
	}

	void init(data *mydata)
	{
		int i, j;
		init_called = true;

		lnfac_array = dvector(mydata->nseq*mydata->hd+3);	//Store lnfac values in array for speed of computation
		lnfac_array[0]=lnfac_array[1]=0;
		for (j=2;j<=((int) mydata->nseq*mydata->hd);j++) lnfac_array[j]=lnfac_array[j-1]+log((double)j);

		window = mini(mydata->lseq,MAXW);
		pij = imatrix(mydata->lseq+1,window+1);
		size_pij = mydata->lseq+1;
		for (i=1;i<=mydata->lseq;i++) for (j=1;j<=window;j++) pij[i][j]=0;

		lij = dmatrix(mydata->lseq+1,(2*window)+1);	// define likelihood matrix

		snp_weights = dvector(mydata->lseq+1);
		calc_snp_weights(mydata);
	}

	void calc_snp_weights(data *mydata)
	{
		int i;
		if (old_lk)
		{
			for (i=0; i<=mydata->lseq; i++)
				snp_weights[i] = 1.0;
		}
		else
		{
			int a,b,d;
			for (i=1; i<=mydata->lseq; i++)
			{
				//d = (2.0*window - 1);
				//d = mini(mydata->lseq - i + 1, (2*window)-1);



				a = maxi(i - window, 1);

				b = mini(i + window, mydata->lseq-1);

				d = mini(b-a, 2*window);
				snp_weights[i] = 1.0/((double)d);
			}
		}
	}

	void read_lk_file(data *mydata)
	{
		int i;
		FILE *ifp;
		
		ifp = fopen(lk_filename.c_str(), "r");
		if (ifp==NULL) 
			nrerror("Cannot open likelihood file");

		pset = NULL;
		pset = init_pset(pset, 1, ifp, &npt, mydata);
		
		//Check that all haplotypes are present in likelihood file
		if (!mydata->lk_exact) 
		{
			int check = check_exhaustive(pset,npt,(int) (mydata->nseq)*((int) mydata->hd));;
			if (check != 1)
			{
				int tempint;
				cout << "\n\nLikelihood file does not seem to be exhaustive." << endl;
				cout << "Is the likelihood pre-calculated for this dataset? (0/1)" << endl;
				scanf("%i", &tempint);
				if (tempint == 1)
				{
					mydata->lk_exact = true;
					pset = init_pset(pset, 1, ifp, &npt, mydata);
				}
				else
				{
					my_exit("Exiting Program", 5);
				}
			}
		}


		printf("Calculating distribution of pair types\n");

		pset = pair_spectrum(mydata, mydata->nall);
		printf("Completed classification of pair types\n");
		printf("Old = %i: New = %i: Missing = %i\n", npt,pnew,miss);
		if (mydata->hd==1 && pnew) nrerror("Cannot have haploid data and new types!!");

		i = (int) (mydata->hd==1? ((int) mydata->nseq/2): mydata->nseq);
		npmax = (int) 1+i+i*(i-1)*(i+4)/6+(i-1)*(i+2)/2;
		printf("Max number of haplotypes for n=%i is %i\n",mydata->nseq*((int) mydata->hd==1?1:2), npmax);

		if (mydata->verbose) 
		{
			FILE *tfp;
			tfp = fopen("type_table", "w");
			if (!tfp) nrerror("Cannot open type file");
			type_print(pij, mydata->lseq, window, tfp);
			fclose(tfp);
		}

		read_pars(ifp, &tcat, &theta_per_site, &rcat, &rmax);

		lkmat = dmatrix(npt+pnew+miss+1,rcat+1);
		read_lk(ifp);
		fclose(ifp);

		if (miss && mydata->hd==1)  
		{
			printf("Calculating LKs for missing data (haplotypes: %i)\n",miss);
			for (i=1;i<=miss;i++) 
			{
				lk_miss(pset[npt+i],lkmat[npt+i],lkmat,mydata, lnfac_array, rcat);
			}
		}
		else if (mydata->hd == 2 && !(mydata->lk_exact)) 
		{
			printf("Resolving diploid data: %i\n",pnew+miss);
			double *lkres;
			lkres = dvector(rcat+1);
			for (i=1;i<=pnew+miss;i++) 
			{
				lk_resolve(lkres,pset[npt+i],lkmat[npt+i],lkmat,mydata, lnfac_array, rcat);
			}
			free_dvector(lkres);
		}

		if (mydata->verbose && !(mydata->lk_exact))  print_lks(mydata, npt+pnew+miss);
	}

	void read_lk(FILE *ifp) 
	{
		int p=1, k;
		int c;

		printf("Reading likelihoods for pair types\n");
		while((c=fgetc(ifp)) != EOF)
			if (c == ':') 
			{
				for (k=1; k<=rcat; k++) 
				{
					fscanf(ifp,"%lf", &lkmat[p][k]);
				}
				p++;
			}
		if (p-1 != npt) {printf("Likelihood file does not agree with header\n"); 
			my_exit("Error in likelihood file", 8);
		}
	}

// Routine to calculate the pairwise likelihood for any given rho
// Returns the change in likelihood due to the rjMCMC move

// Updates the temp section of lij
// lij is the matrix containing the likelihood for each pair of sites in current state
// pij is the matrix of pair types
// rl and ru and the lower and upper bounds of the regions being updated
// update determines the update region - block and links for type 0, block, links and across for type 1
// dlk is a pointer to the difference in likelihood following update
// lkmat is the matrix  containing the precalculated likelihoods from the likelihood file
	double lk_calc(int rl, int ru, data *mydata, double *rmap, int update) 
	{
		int i, j, k, m, t;
		double d;
		double cij;
		double dlk=0.0;
		double rcatDIVrmax;
		int im, jm;

	// To check MCMC routine - just return no change log likelihood
		if (no_lk == true)
			return 0.0;
		int a, b;

		rcatDIVrmax = ((double)(rcat-1))/rmax;

		if (!update) 
		{//Only update blocks and associations to blocks
			for (i=1; i<=ru; i++)
			{
				if (i<rl) a = rl;
				else a = i+1;
				if (i<rl) b = mini(ru, i+window);
				else b = mini(mydata->lseq, i+window);
				im = i - 1;
				for (j = a; j <= b; j++) 
				{//Upper and lower bound depends on whether i<rl
					t = pij[i][j-i];	//Type of pair
					if (t > 0) 
					{
						jm = j - 1;
						cij = rmap[jm]-rmap[im];	//4Ner between sites
						m = j-i+window;
						if (cij >= rmax) 
						{
							lij[i][m] = lkmat[t][rcat];
						}
						else 
						{
							d = cij*rcatDIVrmax;	//d = (double) cij*(mydata->rcat-1)/(mydata->rmax);
							k = (int) d+1;
							lij[i][m] = lkmat[t][k]+(d-k+1)*(lkmat[t][k+1]-lkmat[t][k]);
						}
						if (lij[i][m] < 0.0) 
							dlk += (lij[i][m] - lij[i][j-i])*snp_weights[i];		// Difference in log-likelihood 
						else 
						{
							printf(" Log Likelihood >= 0!!  Code Location:A t=%i cij=%f\n",t,cij);
							my_exit("Error in Likelihood", 3); return 0;
						}
					}
				}
			} 
		}
		else 
		{//Update block, associations to block and associations across blocks
			for (i=1; i<=ru; i++) 
			{
				if (i<rl) a = rl;
				else a = i+1;
				b = mini(mydata->lseq, i+window);
				im = i - 1;
				for (j = a; j <= b; j++) 
				{
					t = pij[i][j-i];
					if (t > 0) 
					{
						m = j-i+window;
						jm = j - 1;
						cij = rmap[jm]-rmap[im];	//4Ner between sites
						if (cij >= rmax) 
						{
							lij[i][m] = lkmat[t][rcat];
						}
						else 
						{
							d = cij*rcatDIVrmax;	//d = (double) cij*(mydata->rcat-1)/(mydata->rmax);
							k = (int) (d+1.0);
							lij[i][m] = lkmat[t][k]+(d-k+1)*(lkmat[t][k+1]-lkmat[t][k]);
						}
						if (lij[i][m] < 0.0) 
						{
							dlk += (lij[i][m] - lij[i][j-i])*snp_weights[i];
						}
						else 
						{
							printf(" Log Likelihood >= 0!! Code Location:B t=%i cij=%f\n",t,cij);
							my_exit("Error in Likelihood", 4); return 0;
						}
					}
				}
			}
		}
		return dlk;
	}


	void print_lks(data *mydata, int NPT) 
	{
		int p, i, nstate, ct=0;
		FILE *ofp;

		ofp = fopen("new_lk", "w");
		if (mydata->hd == 2) nstate=16;
		else nstate=9;

		for (p=1;p<=NPT;p++) if (pset[p]->nt>0) ct++;

		fprintf(ofp, "\n%i %i\n1 %.5f\n%i %f\n\n",mydata->nseq*((int)(mydata->hd==1?1:2)),ct,theta_per_site,rcat,rmax);
		for (p=1; p<=NPT; p++) if (pset[p]->nt>0) 
		{
			fprintf(ofp,"\n%4i # ", p);
			for (i=0; i<nstate; i++) fprintf(ofp,"%3i ", pset[p]->pt[i]);
			fprintf(ofp," :  ");
			for (i=1; i<=rcat; i++) fprintf(ofp,"%7.2f ", lkmat[p][i]);
		}
		fclose(ofp);
	}


	// Copies the temp section of lij to the permanent section of lij
	void update_lij(int rl, int ru, data *mydata, int update)
	{
		int i, j;
		int a, b;

		if (update) 
		{	// Have to update across block
			for (i=1;i<=ru;i++) 
			{
				if (i<rl) a = rl;
				else a = i+1;
				b = mini(mydata->lseq, i+window);
				for (j=a;j<=b; j++) 
					lij[i][j-i]=lij[i][j-i+window];
			}
		}
		else 
		{	//Just update within block+associations
			for (i=1;i<=ru;i++) 
			{
				if (i<rl) a = rl;
				else a = i+1;
				if (i<rl) b = mini(ru, i+window);
				else b = mini(mydata->lseq, i+window);
				for (j=a;j<=b; j++) 
					lij[i][j-i]=lij[i][j-i+window];
			}
		}
	}

	// Routine recalculates likelihood of the data to prevent loss of accuracy
	double update_lk0(data *mydata)
	{
		int i, j;
		double dlk = 0.0;
		for (i=1;i<mydata->lseq;i++) 
			for (j=i+1;j<=mini(mydata->lseq,i+window);j++)
				dlk += lij[i][j-i];
		return dlk;
	}

	void print_lij(data *mydata)
	{
		int i, j;
		FILE *out;
		out = fopen("lij.txt", "w");
		for (i=1;i<mydata->lseq;i++)
		{
			for (j=i+1;j<=mini(mydata->lseq,i+window);j++)
				fprintf(out, "%f\t", lij[i][j-i]);
			fprintf(out, "\n");
		}
		fclose(out);
	}


	// Routine to classify each pairwise comparison
	struct site_type ** pair_spectrum(data *mydata, int **nall) 
	{

	/* pt and type (haploid)
	00: 0
	10: 1
	01: 2
	11: 3
	0?: 4
	1?: 5
	?0: 6
	?1: 7
	??: 8

	pt and type (diploid)

	??_??:0
	??_00:1
	??_11:2
	??_10:3
	00_??:4
	00_00:5
	00_11:6
	00_10:7
	11_??:8
	11_00:9
	11_11:10
	11_10:11
	10_??:12
	10_00:13
	10_11:14
	10_10:15
	*/

		int i, j, seq, sites[2], states[2][2], *pt, nstate, ct;
		//char bases[6]="n-TCAG";
		extern int sizeofpset;

		if (mydata->hd==1) nstate=9; 
		else nstate=16;
		//pt  = (int *) malloc((size_t) (nstate+1)*sizeof(int));
		pt = ivector(nstate+1);
		if (pt == NULL) nrerror("pt");

		for (sites[0]=1; sites[0]<mydata->lseq; sites[0]++) 
		{
	/*	  printf("\nSite %4i: ",sites[0]);*/
			for (sites[1]=sites[0]+1; sites[1]<=mini(mydata->lseq,sites[0]+window); sites[1]++) 
			{
				if (!check22(sites[0], sites[1], nall)) 
					pij[sites[0]][sites[1]-sites[0]]=0; 
				else 
				{
				   for (i=0; i<nstate; i++) 
					   pt[i]=0;

				   if (mydata->hd==1) 
				   {
						for (j=0; j<2; j++) 
						{
							for (i=2; i<=5; i++) if (nall[sites[j]][i]) {states[j][0]=i;break;}
							for (i++; i<=5; i++) if (nall[sites[j]][i]) {states[j][1]=i;break;}
						}
						for (seq=1; seq<=mydata->nseq; seq++) 
						{
							if (mydata->seqs[seq][sites[0]]==states[0][0]) 
							{
								if (mydata->seqs[seq][sites[1]]==states[1][0]) pt[0]++;
								else if (mydata->seqs[seq][sites[1]]==states[1][1]) pt[2]++;
								else pt[4]++;
							}
							else if (mydata->seqs[seq][sites[0]]==states[0][1]) 
							{
								if (mydata->seqs[seq][sites[1]]==states[1][0]) pt[1]++;
								else if (mydata->seqs[seq][sites[1]]==states[1][1]) pt[3]++;
								else pt[5]++;
							}
							else 
							{
								if (mydata->seqs[seq][sites[1]]==states[1][0]) pt[6]++;
								else if (mydata->seqs[seq][sites[1]]==states[1][1]) pt[7]++;
								else pt[8]++;
							}
						}
						pt =order_pt_hap(pt,mydata->nseq);
					}

					else 
					{
						for (seq=1;seq<=mydata->nseq;seq++) 
							pt[4*(mydata->seqs[seq][sites[0]]-1)+mydata->seqs[seq][sites[1]]-1]++; 
						pt = order_pt_dip(pt, 2*mydata->nseq);
						for (j=0,ct=0; j<nstate;j++) 
							{ct+=pt[j];}
						if (ct != mydata->nseq) nrerror("Error in calculating pair type");
					}
				   
					if ((npt)+(pnew)+(miss)+10>sizeofpset) 
						pset = add_pset(pset);
					pij[sites[0]][sites[1]-sites[0]]=add_type(pset, pt, &npt, &pnew, &miss, mydata);
				}
			}
		}
		free_ivector(pt);
		return pset;
	}

};

#endif


