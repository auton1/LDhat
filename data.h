#if !defined DATA_H
#define DATA_H
#pragma warning(disable:4786) 
#include <vector>
#include <string>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <algorithm>

using namespace std;

#include "rhomap_tools.h"


#define MAXNAME 255						// Max length of sequences names
#define MAXLINE 2048					// Max length of line


class data
{
public:
	
	string seq_filename;				// Sequence filename
	string loc_filename;				// Locus filename
	string freq_filename;				// Allele frequency file
	string params_filename;				// Paramater output file (for debugging)
	int **seqs;							// Sequence data
	//vector<  string  > seqnames;		// sequence header data
	double *locs;						// Loci data
	int lseq;							// length of seq (in snps)
	int nseq;							// number of seqs
	int hd;								// haplotype/genotype data
	int ns;								// Number of sites
	double tlseq;						// Total length of sequence in kb?
	char lc;							// Crossing-over (L) or gene conversion (C) model (no longer supported)
	bool verbose;						// Verbose output mode
	bool output_freqs;					// Output allele frequencies to file

	double theta;						// Watterson estimate of theta

	bool lk_exact;						// Likelihood file is precalculated for this dataset
	int random_seq_subset_size;			// Use a random subset of sequences to match likelihood file
	vector<int> order;					// Used to define subset ordering

	int **nall;			// Matrix of allele frequencies

	data()
	{
		verbose=false;
		output_freqs=false;
		lk_exact=false;
		random_seq_subset_size = 0;

		seq_filename = "";
		loc_filename = "";

		params_filename = "params.txt";
	}

	~data()
	{
		free_dvector(locs);
		free_imatrix(seqs,nseq+1);
		free_imatrix(nall,lseq+1);
	}

	void readseqs()
	{
		int i, site, seq=1, cts[5];
		char line[MAXLINE], line2[MAXLINE], bases[6]="TCAG-";
		const char* filename = seq_filename.c_str();
		FILE *ifp;
		ifp = fopen(filename, "r"); 
		if (ifp == NULL) nrerror("Error in opening sequence file");
		fscanf(ifp,"%i %i %i", &nseq, &lseq, &hd);

		if ((nseq < 2) || (lseq < 2)) {my_exit("Insufficient data for analysis (n > 1, L > 1)", 10);}

		if (hd==1) printf("Analysing %i haploid sequences of length %i seg sites\n", nseq, lseq);
		else printf("Analysing %i genotypes of length %i seg sites\n", nseq, lseq);

		//seqnames.resize(nseq);
		seqs = imatrix(nseq+1, lseq+1);
	
		printf("Reading sequences in fasta format\n");
		for (i=0;i<5;i++) cts[i]=0;

		while (!feof(ifp) && (seq<nseq+1)) 
		{
			fgets(line, MAXLINE, ifp);
			if (line[0] == '>') 
			{
				if (verbose) printf("Sequence :%3i ", seq);
				//seqnames[seq-1] = line+1;
				//seqnames[seq-1].replace(seqnames[seq-1].find("\n"), 2, "\0");
				//if (verbose) cout << seqnames[seq-1] << endl;
				site=1;
				for (i=0;i<5;i++) cts[i]=0;
				while (site<lseq+1) 
				{
					fgets(line2, MAXLINE, ifp);
					i = 0;
					while(line2[i] != '\0')
					{
					   switch(line2[i]) 
					   {
							case 'T': case 't': case '0':	// Homozygote 1
							{
								seqs[seq][site] = 2;
								site++;
								cts[0]++;
								break;
							}
							case 'C': case 'c': case '1': 	// Homozygote 2
							{
								seqs[seq][site] = 3;
								site++;
								cts[1]++;
								break;
							}
							case 'A': case 'a': case '2': 	// Heterozygote
							{
								seqs[seq][site] = 4;
								site++;
								cts[2]++;
								break;
							}
							case 'G': case 'g': case '3':
							{
								seqs[seq][site] = 5;
								site++;
								cts[3]++;
								break;
							}
							case '-': case'N': case'n': case'?': case'R': case'Y': case'M': case'K': case'S': case'W': case'H': case'B': case'V': case'D' :
							{
								seqs[seq][site]=1;
								site++;
								cts[4]++;
								break;
							}
							case '>': 
							{
								printf("%s\n\n", line2);
								printf("\nError in sequence file(%i of %i bases read)\n",site,lseq); my_exit("Exiting", 11);
							}
							default: 
							{
								break;
							}
						}
						i++;
					}
				}
				seq++;
				if (site != lseq+1) 
				{
					nrerror("Sequences incorrect length");
				}
				if (verbose) 
				{	
					for (i=0;i<5;i++) printf("%c:%5i ",bases[i], cts[i]);
					printf("\n");
				}
			}
		}
		if (seq!=nseq+1) {
			printf("\n\nDid not read %i sequences \n\n",nseq); my_exit("Exiting", 12);}
		fclose(ifp);

		if (random_seq_subset_size > 0)
			take_random_sequence_subset();		// Take a subset of sequences to match LK file

		theta = watterson(nseq);
	}

	void take_random_sequence_subset()
	{
		// Function to reduce the number of sequences by taking a random subset.
		int i, j;
		int **new_seqs;
		new_seqs = imatrix(random_seq_subset_size+1, lseq+1);

		order.resize(random_seq_subset_size+1);
		for (i=0; i<random_seq_subset_size+1; i++)
			order[i] = i;

		random_shuffle(order.begin()+1, order.end());

		for (i=1; i<random_seq_subset_size+1; i++)
			for (j=0; j<lseq+1; j++)
				new_seqs[i][j] = seqs[order[i]][j];

		free_imatrix(seqs,nseq+1);
		seqs = new_seqs;
		nseq = random_seq_subset_size;

	}

	void readlocs()
	{
		int i;
		const char* filename = loc_filename.c_str();
		FILE *ifp;
		ifp = fopen(filename, "r");
		if (ifp == NULL) nrerror("Error in opening locus file");
		fscanf(ifp, "%i %lf %c", &ns, &tlseq, &lc);
		locs = dvector(lseq+2);
		for (i=1; i<=lseq; i++) 
			fscanf(ifp, "%lf", &locs[i]);
		fclose(ifp);
		locs[lseq+1] = locs[lseq];
		tlseq = locs[lseq] - locs[1];
	}

	void allele_count()
	{
		int seq, site, i, ct;
		int temp;
		FILE *ofp;

		nall = imatrix(lseq+1, 7);
        for (site=1; site<=lseq; site++) 
		{
			for (seq=1; seq<=nseq; seq++) 
				nall[site][seqs[seq][site]]++;
			if (hd==2) 
			{
				nall[site][2]=2*nall[site][2]+nall[site][4];	// Allele 1
				nall[site][3]=2*nall[site][3]+nall[site][4];	// Allele 2
				nall[site][4]=0;
				nall[site][5]=0;
				nall[site][1]*=2;								// Missing data
			}
		}
		if (output_freqs) 
		{
			ofp = fopen(freq_filename.c_str(),"w");
			if (ofp == NULL) nrerror("Error in opening freqs file");
			temp = nseq*hd;
			fprintf(ofp,"\nAllele frequencies\n\n Site   -   T/0  C/1  A/2  G/3\n\n");
			for (site=1; site<=lseq; site++) 
			{
			fprintf(ofp,"%4i ", site);
			for (i=1,ct=0; i<=5; i++) {fprintf(ofp,"%4i ", nall[site][i]); ct+=nall[site][i];}
			if (ct != temp) nrerror("Error in counting allele frequencies");
			fprintf(ofp,"\n");
			}
		fclose(ofp);
		}
	}

	double watterson(int n) 
	{
		// Wattersons estimator of theta
		int i;
		double cump=1.0;

		for (i=2;i<n;i++) cump+=(double) 1/i;
		return 1.0 / cump;
	}

};

#endif
