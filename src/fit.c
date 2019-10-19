///////////////////////////////////////////////////////////////////////////////////////////////////
// Emerge code - File fit.c                                                                      //
///////////////////////////////////////////////////////////////////////////////////////////////////
///
/// \file fit.c
/// \brief Contains functions that fit the model parameters
///
/// This file contains the functions that fit the parameters of the model given the data.
/// The implemented methods include an affine invariant ensemble sampler MCMC in #aies_mcmc()
/// and the hybrid optimisation method in #optimize_hybrid(). There are also functions that
/// initialise the parameters on each task and that compute the probability on each master task
/// for a given set of parameters. This can be done on multiple universes in parallel. To do this
/// there are several functions that compute \f$\chi^2\f$ for a given statistics.
///
///////////////////////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "allvars.h"
#include "proto.h"


/*! \brief This function runs the affine invariant ensemble sampler MCMC
 *
 *  The MCMC can be either started (restart==0) or restarted from a stored walker file (restart==1).
 *  If newly started it first initialises all walkers with parameters according to the specified range
 *  in the parameter file. Then the probability for each walker is computed for all universes in
 *  parallel. If restarted all walkers including their parameters and probabilities are read from a
 *  walker file. In the main MCMC loop the walkers are split into two sets with equal numbers of walkers.
 *  For each walker a walker from the other set is randomly selected. The walker parameters are then
 *  'stretched' along the direction to the other walker and the probability is computed. This is done
 *  in parallel for all universes. A random number [0,1) is then computed and compared to the acceptance
 *  probability. If the new walker is accepted it replaces the old one. During each chain step all walkers
 *  are printed to a file. After each step a walker file is written containing the parameters and
 *  probabilities for all walkers. The chain is terminated if the predicted ending time for the next
 *  step is longer then the time limit specified in the parameter file. This algorithm is an implementation
 *  of the method presented by <a href="http://msp.org/camcos/2010/5-1/p04.xhtml"> Goodman & Weare (2010)</a>.
 *
 *  \param restart Flag indicating if the MCMC chain is started (0) or restarted from a walker file (1)
 */
void aies_mcmc(int restart)
{
	const int nline = 500;
	char outfname[NSTRING], infname[NSTRING], line[nline], *tmp;
	int i,j,k,iwalker,jwalker,istart,itask,istep,isave,isplit,doloop,accept,nrun,nfilepar,nfilewalk;
	float px,py;
	double q,lastlooptime;
	double acctot = 0.0;
	FILE *ifp,*ofp,*wfp;
	MPI_Status status;
	gsl_rng *rng_uniform_mcmc;

	//Define the scale parameter for the MCMC
	double a = sqrt(All.mcmca);

	//Define the walker array x, the proposed walker array y, and the probability array prob
	double *x;     //Walkers, Length: Nwalkers*Nparam
	double *y;     //Updated walker
	double *z;     //Scale of the proposed step
	float  *prob;  //Current probability for each walker

	//The function can only be called with a restart of 0 (start new chain) or 1 (restart from walker file)
	if (restart < 0 || restart > 1)
	{
		if (ThisTask == 0) printf("Affine invariant ensemble sampler MCMC called with wrong restart flag");
		return;
	}

	//Initialise the uniform generator
	rng_uniform_mcmc = gsl_rng_alloc(gsl_rng_ranlxd2);
	gsl_rng_set(rng_uniform_mcmc, All.Seed);

	//Allocate all arrays
	prob       = emalloc_movable(&prob, "Prob", All.Nwalkers * sizeof(float));
	x          = emalloc_movable(&x,    "X",    All.Nwalkers * All.Nparam * sizeof(double));
	y          = emalloc_movable(&y,    "Y",    All.Nparam * All.NUniverses * sizeof(double));
	z          = emalloc_movable(&z,    "Z",    All.NUniverses * sizeof(double));

	//This is the first run
	nrun = 0;

	//If this is the first run of the MCMC then set up the run
	if (restart == 0)
	{
		//Only the main task opens and writes to the output file and prints to screen
		if (ThisTask == 0)
		{
			//Open the MCMC output file
			sprintf(outfname,"%s/mcmc%d.%03d.out",All.OutputDir,All.Seed,nrun);
			ofp = fopen(outfname,"w");

			//Print what is done...
			printf("%s\n%s Starting affine invariant ensemble sampler MCMC...\n",All.fullline,All.startline);
			printf("%s Setting up %d Walkers...\n",All.startline,All.Nwalkers);

			//Print the header for the MCMC steps
			print_mcmc_header();
		}

		//Initialization of walkers
		init_param(x);

		//Calculate the starting probability for each walker
		for (j = 0; j < All.Nwalkers/All.NUniverses; j++)
		{
			//Compute the probability for each universe
			px   = lnprob(&x[j*All.NUniverses*All.Nparam]);
			//Set the current time
			All.timenow = second() - All.starttime;
			//Go through all master tasks
			for (itask = 0; itask < All.NUniverses; itask++)
			{
				//If this is a master task but not task 0
				if (ThisTask == All.NTaskPerUniverse * itask && itask > 0)
				{
					//Send the probability to task 0
					MPI_Ssend(&px, 1, MPI_FLOAT, 0, TAG_PROB, MPI_COMM_WORLD);
				}
				//If this is task 0 and it is not the first loop
				if (ThisTask == 0 && itask > 0)
				{
					//Receive the probability from the other master task
					MPI_Recv(&px, 1, MPI_FLOAT, All.NTaskPerUniverse * itask, TAG_PROB, MPI_COMM_WORLD, &status);
				}
				//Write the probability to the correct index (only valid for task 0)
				prob[j*All.NUniverses+itask] = px;
				//If this is task 0 write all parameters and the probability to screen and output file
				if (ThisTask == 0)
				{
					printf("%8d %8.1f",j*All.NUniverses+itask,All.timenow);
					for (i = 0; i < All.Nparam; i++) printf(" %11.8f",x[(j*All.NUniverses+itask)*All.Nparam+i]);
					printf(" %f\n",prob[j*All.NUniverses+itask]);

					fprintf(ofp,"%f",prob[j*All.NUniverses+itask]);
					for (i = 0; i < All.Nparam; i++) fprintf(ofp," %f",x[(j*All.NUniverses+itask)*All.Nparam+i]);
					fprintf(ofp," 0 0 0 %d\n",j*All.NUniverses+itask);
					fflush(ofp);
				}
			}//End loop through master tasks
		}//End loop through all walkers

		//Set the start index to zero
		istart = 0;

	}//End restart 0

	//If the chain is restarted from a file with walkers
	if (restart == 1)
	{
		//Only the main task opens and writes to the output file and prints to screen
		if (ThisTask == 0)
		{
			//Open the walker file
			sprintf(infname,"%s/walkers.mcmc.%d.out",All.OutputDir,All.Seed);
			if(!(ifp=fopen(infname,"r")))
			{
				sprintf(line, "Cannot open file `%s' for reading the walkers.", infname);
				endrun(line);
			}

			//Get first line of walker file
			fgets(line,nline,ifp);
			sscanf(line,"#Walkers/Parameters/seed/i %d %d %d %d\n", &nfilewalk, &nfilepar, &nrun, &istart);

			//Sanity Checks
			if (nfilewalk != All.Nwalkers)
			{
				sprintf(line, "The walker file %s contains %d walkers (while the code needs %d)",infname,nfilewalk,All.Nwalkers);
				endrun(line);
			}
			if (nfilepar != All.Nparam)
			{
				sprintf(line, "The walker file %s contains %d parameters (while the code needs %d)",infname,nfilepar,All.Nparam);
				endrun(line);
			}
			if (nrun < 0)
			{
				sprintf(line, "The walker file %s has a negative run number (%02d)",infname,nrun);
				endrun(line);
			}

			//Increment the Run Number
			nrun++;

			//Increment the starting step
			istart++;

			//Print what is done...
			printf("%s\n%s Restarting affine invariant ensemble sampler MCMC...\n",All.fullline,All.startline);
			printf("%s Reading %d Walkers from %s...\n",All.startline,All.Nwalkers,infname);

			//Set output file and open file pointer
			sprintf(outfname,"%s/mcmc%d.%03d.out",All.OutputDir,All.Seed,nrun);
			ofp = fopen(outfname,"w");

			//Get each line in the file
			for (k = 0; k < All.Nwalkers; k++)
			{
				//Get the new line
				fgets(line,nline,ifp);
				tmp = strtok(line, " ");
				//Get the probability
				prob[k] = atof(tmp);
				//Get all the parameters for the walker
				for (i = 0; i < All.Nparam; i++)
				{
					tmp = strtok(NULL, " ");
					x[k*All.Nparam+i] = atof(tmp);
				}
			} //End reading walker file
			fclose(ifp);
		}

		//Re-initialize random numbers
		gsl_rng_set(rng_uniform_mcmc, All.Seed + 100 * nrun);

		//Broadcast the walker array x from task 0
		MPI_Bcast(x, All.Nwalkers * All.Nparam, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}

	//Broadcast the probability array from task 0
	MPI_Bcast(prob, All.Nwalkers, MPI_FLOAT, 0, MPI_COMM_WORLD);

	//Define loop iteration variable
	istep=istart;
	isave=istart*All.Nwalkers;
	doloop=1;

	//Inform everyone we are starting with the main MCMC loop
	if (ThisTask == 0) printf("%s\n%s Starting main MCMC loop...\n",All.fullline,All.startline);

	//Start main MCMC loop, stop loop after some specified time is reached
	while (doloop)
	{
		//Initialise Loop start time
		lastlooptime = second() - All.starttime;

		//Print the header for the MCMC steps
		if (ThisTask == 0) print_mcmc_header();

		//Split the walkers in two sets with the same number of walkers
		for (isplit = 0; isplit < 2; isplit++)
		{
			//Loop through half of the walkers
			for (j = 0; j < All.Nwalkers/All.NUniverses/2; j++)
			{
				//Set the current time
				All.timenow = second() - All.starttime;
				//Set up new trial walkers for all universes

				for (itask = 0; itask < All.NUniverses; itask++)
				{
					//Current walker
					iwalker = isplit*All.Nwalkers/2 + j*All.NUniverses + itask;
					//For walker iwalker, select random walker jwalker from other set
					jwalker = (int)((gsl_rng_uniform(rng_uniform_mcmc)+1.0-(double)(isplit))*(double)(All.Nwalkers)*0.5);
					//Draw random number from g(z) = (z)^(-1/2) within [1/a,a]
					z[itask] = pow((a-1.0)*gsl_rng_uniform(rng_uniform_mcmc)+1.0,2.0)/a;
					//Stretch walker iwalker along jwalker to obtain new propasal walker y
					for (i = 0; i < All.Nparam; i++)
						y[itask*All.Nparam+i] = x[jwalker*All.Nparam+i] + z[itask] * (x[iwalker*All.Nparam+i] - x[jwalker*All.Nparam+i]);
				}

				//Broadcast the walker array y from task 0
				MPI_Bcast(y, All.NUniverses * All.Nparam, MPI_DOUBLE, 0, MPI_COMM_WORLD);

				//Compute the probability for each universe
				py = lnprob(y);
				MPI_Barrier(MPI_COMM_WORLD);

				//Go through all master tasks
				for (itask = 0; itask < All.NUniverses; itask++)
				{
					//Set index of the current walker
					iwalker = isplit*All.Nwalkers/2 + j*All.NUniverses + itask;
					//Get this walkers previous probability
					px   = prob[iwalker];
					//If this is a master task but not task 0
					if (ThisTask == All.NTaskPerUniverse * itask && itask > 0)
					{
						//Send the probability to task 0
						MPI_Ssend(&py, 1, MPI_FLOAT, 0, TAG_PROB, MPI_COMM_WORLD);
					}
					//If this is task 0 and it is not the first loop
					if (ThisTask == 0 && itask > 0)
					{
						//Receive the probability from the other master task
						MPI_Recv(&py, 1, MPI_FLOAT, All.NTaskPerUniverse * itask, TAG_PROB, MPI_COMM_WORLD, &status);
					}

					//Sanity check if the computed probability is in a sensible range
					if (!(py < 0 && py > -1.e10) ) py = -2.e10;

					//Get probability to accept walker and decide whether to accept y
					q = ((double)(All.Nparam)-1.0)*log(z[itask])+py-px;
					accept = (exp(q) > gsl_rng_uniform(rng_uniform_mcmc));

					//If the probability is beyond the edge do not accept
					if (py < -1.e10) accept = 0;

					//If this is task 0 write all parameters and the probability to screen and output file
					if (ThisTask == 0)
					{
						//If y is accepted set walker k to y and store probability and update accepted fraction
						if (accept)
						{
							for (i = 0; i < All.Nparam; i++) x[iwalker*All.Nparam+i] = y[itask*All.Nparam+i];
							prob[iwalker] = py;
							acctot+=1.0;
						}

						//Print the walkers to screen
						printf("%8d %8.1f",isave,All.timenow);
						for (i = 0; i < All.Nparam; i++) printf(" %11.8f",x[iwalker*All.Nparam+i]);
						printf(" %f\n",prob[iwalker]);

						//Output to file
						fprintf(ofp,"%f",prob[iwalker]);
						for (i = 0; i < All.Nparam; i++) fprintf(ofp," %f",x[iwalker*All.Nparam+i]);
						fprintf(ofp," %f %d %d %d\n",acctot/(isave+1),isave,istep,iwalker);
						fflush(ofp);

					}//End task 0

					//Update output index
					isave++;

				}//loop through master tasks

				//Broadcast the walker array x from task 0
				MPI_Bcast(x, All.Nwalkers * All.Nparam, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				//Broadcast the probability array from task 0
				MPI_Bcast(prob, All.Nwalkers, MPI_FLOAT, 0, MPI_COMM_WORLD);

			}//loop through half walkers
		}//isplit

		//Only task 0 writes the walker file
		if (ThisTask == 0)
		{
			//Write all current walkers to a file
			sprintf(outfname,"%s/walkers.mcmc.%d.out",All.OutputDir,All.Seed);
			wfp = fopen(outfname,"w");
			fprintf(wfp,"#Walkers/Parameters/seed/i %d %d %d %d\n",All.Nwalkers,All.Nparam,nrun,istep);
			for (iwalker = 0; iwalker < All.Nwalkers; iwalker++)
			{
				fprintf(wfp,"%f",prob[iwalker]);
				for (i = 0; i < All.Nparam; i++) fprintf(wfp," %f",x[iwalker*All.Nparam+i]);
				fprintf(wfp,"\n");
			}
			fclose(wfp);
		}

		//Update step
		istep++;

		//Get time the initialisation loop took
		lastlooptime = second() - All.starttime -lastlooptime;
		All.timenow  = second() - All.starttime;

		//End loop if projected finish time of next loop would be larger than the time limit
		if (All.timenow+1.5*lastlooptime>All.timelimit) doloop = 0;

	}//loop

	//Inform everyone we are done with the main MCMC loop
	if (ThisTask == 0) printf("%s Done with MCMC...\n",All.startline);

	//Free all arrays and the output file pointer
	efree_movable(z);
	efree_movable(y);
	efree_movable(x);
	efree_movable(prob);
	if (ThisTask == 0) fclose(ofp);
	gsl_rng_free(rng_uniform_mcmc);

}


/*! \brief This function runs the HYBRID optimization method
 *
 *  The method can be either started (restart==0) or restarted from a stored walker file (restart==1).
 *  If newly started it first initialises all walkers with parameters according to the specified range
 *  in the parameter file. Then \f$\chi^2\f$ for each walker is computed for all universes in
 *  parallel. Also the mean starting \f$\chi^2\f$ is computed and stored. If restarted all walkers including
 *  their parameters and probabilities are read from a walker file. In the main HYBRID loop first the
 *  mean current \f$\chi^2\f$ over all walkers is computed and broadcast from task 0. The value q (mean current
 *  \f$\chi^2\f$ over mean starting \f$\chi^2\f$) and the resulting g value are computed. For each walker the p
 *  value (\f$\chi^2\f$ over mean current \f$\chi^2\f$) and the resulting f value are computed. Then a new proposal
 *  walker is computed from the current walker its width sigma and the g and f values and the probability
 *  is computed. This is done in parallel for all universes. A random number [0,1) is then computed and
 *  compared to the acceptance probability given the current temperature. If the new walker is accepted
 *  it replaces the old one. During each chain step all walkers are printed to a file. After each step
 *  a walker file is written containing the parameters and probabilities for all walkers. The chain is
 *  terminated if the predicted ending time for the next step is longer then the time limit specified
 *  in the parameter file. This algorithm is an implementation of the method presented by
 *  <a href="https://arxiv.org/abs/astro-ph/0602338"> Elson et al. (2007)</a>.
 *
 *  \param restart Flag indicating if the HYBRID chain is started (0) or restarted from a walker file (1)
 */
void optimize_hybrid(int restart)
{

//Define slopes
#ifndef HYBRID_ALPHA
#define HYBRID_ALPHA 0.4
#endif

#ifndef HYBRID_BETA
#define HYBRID_BETA 1.0
#endif

#ifndef HYBRID_GAMMA
#define HYBRID_GAMMA 1.0
#endif

	const int nline = 500;
	char outfname[NSTRING], infname[NSTRING], line[nline], *tmp;
	int i,j,k,iwalker,istart,itask,istep,isave,doloop,accept,nrun,nfilepar,nfilewalk;
	float cx,cy,mean_starting_chi2,mean_current_chi2,min_current_chi2,max_current_chi2;
	double alpha,q,p,f,g,lastlooptime;
	FILE *ifp,*ofp,*wfp;
	MPI_Status status;
	gsl_rng *rng_uniform_hybrid, *rng_gaussian_hybrid;

	int    *ac;    //Acceptance counter for each walker
	double *x;     //Walkers, Length: Nwalkers*Nparam
	double *y;     //Updated walker
	double *sigma; //Width for the chosen parameter space
	float  *chi2;  //Current chi^2 for each walker

	//The function can only be called with a restart of 0 (start new chain) or 1 (restart from walker file)
	if (restart < 0 || restart > 1)
	{
		if (ThisTask == 0) printf("HYBRID optimization method called with wrong restart flag");
		return;
	}

	//Allocate all arrays
	ac    = emalloc_movable(&ac,    "AC",    All.Nwalkers * sizeof(int));
	chi2  = emalloc_movable(&chi2,  "Chi2",  All.Nwalkers * sizeof(float));
	sigma = emalloc_movable(&sigma, "Sigma", All.Nparam * sizeof(double));
	x     = emalloc_movable(&x,     "X",     All.Nwalkers * All.Nparam * sizeof(double));
	y     = emalloc_movable(&y,     "Y",     All.Nparam * All.NUniverses * sizeof(double));

	//This is the first run
	nrun = 0;

	// Get a sigma (width) for the chosen parameter space
	init_sigma(sigma);

	//Initialise the uniform generator
	rng_uniform_hybrid = gsl_rng_alloc(gsl_rng_ranlxd2);
	gsl_rng_set(rng_uniform_hybrid, All.Seed);

	//Initialise the gaussian generator
	rng_gaussian_hybrid = gsl_rng_alloc(gsl_rng_ranlxd2);
	gsl_rng_set(rng_gaussian_hybrid, All.Seed);

	//If this is the first run of the MCMC then set up the run
	if (restart == 0)
	{
		//Only the main task opens and writes to the output file and prints to screen
		if (ThisTask == 0)
		{
			//Open the MCMC output file
			sprintf(outfname,"%s/hybrid%d.%03d.out",All.OutputDir,All.Seed,nrun);
			ofp = fopen(outfname,"w");

			//Print what is done...
			printf("%s\n%s Starting HYBRID optimization method...\n",All.fullline,All.startline);
			printf("%s Setting up %d Walkers...\n",All.startline,All.Nwalkers);

			//Print the header for the chain steps
			print_mcmc_header();
		}

		//Initialization of walkers
		init_param(x);

		//Calculate the starting chi^2 for each walker
		for (j = 0; j < All.Nwalkers/All.NUniverses; j++)
		{
			//Compute the probability for each universe
			cx   = -2.0 * lnprob(&x[j*All.NUniverses*All.Nparam]);
			//Set the current time
			All.timenow = second() - All.starttime;
			//Go through all master tasks
			for (itask = 0; itask < All.NUniverses; itask++)
			{
				//If this is a master task but not task 0
				if (ThisTask == All.NTaskPerUniverse * itask && itask > 0)
				{
					//Send the chi^2 to task 0
					MPI_Ssend(&cx, 1, MPI_FLOAT, 0, TAG_PROB, MPI_COMM_WORLD);
				}
				//If this is task 0 and it is not the first loop
				if (ThisTask == 0 && itask > 0)
				{
					//Receive the chi^2 from the other master task
					MPI_Recv(&cx, 1, MPI_FLOAT, All.NTaskPerUniverse * itask, TAG_PROB, MPI_COMM_WORLD, &status);
				}
				//Write the chi^2 to the correct index (only valid for task 0)
				chi2[j*All.NUniverses+itask] = cx;
				//If this is task 0 write all parameters and the chi^2 to screen and output file
				if (ThisTask == 0)
				{
					printf("%8d %8.1f",j*All.NUniverses+itask,All.timenow);
					for (i = 0; i < All.Nparam; i++) printf(" %11.8f",x[(j*All.NUniverses+itask)*All.Nparam+i]);
					printf(" %f\n",-0.5*chi2[j*All.NUniverses+itask]);

					fprintf(ofp,"%f",chi2[j*All.NUniverses+itask]);
					for (i = 0; i < All.Nparam; i++) fprintf(ofp," %f",x[(j*All.NUniverses+itask)*All.Nparam+i]);
					fprintf(ofp," 0 0 0 0 0 0 0 0 0 %d\n",j*All.NUniverses+itask);
					fflush(ofp);
				}
			}//End loop through master tasks
		}//End loop through all walkers

		// Compute mean starting chi2
		mean_starting_chi2 = 0;
		for (k = 0; k < All.Nwalkers; k++) mean_starting_chi2 += chi2[k];
		mean_starting_chi2 /= All.Nwalkers;

		// Initialise acceptance counters and starting step
		for (k = 0; k < All.Nwalkers; k++) ac[k] = 0;
		istart = 0;

	}//End restart 0

	//If the chain is restarted from a file with walkers
	if (restart == 1)
	{
		//Only the main task opens and writes to the output file and prints to screen
		if (ThisTask == 0)
		{
			//Open the walker file
			sprintf(infname,"%s/walkers.hybrid.%d.out",All.OutputDir,All.Seed);
			if(!(ifp=fopen(infname,"r")))
			{
				sprintf(line, "Cannot open file `%s' for reading the walkers.", infname);
				endrun(line);
			}

			//Get first line of walker file
			fgets(line,nline,ifp);
			sscanf(line,"#Walkers/Parameters/seed/i/mean_starting_chi2/Temperature %d %d %d %d %f %lf\n", &nfilewalk, &nfilepar, &nrun, &istart, &mean_starting_chi2, &All.temperature);

			//Sanity Checks
			if (nfilewalk != All.Nwalkers)
			{
				sprintf(line, "The walker file %s contains %d walkers (while the code needs %d)",infname,nfilewalk,All.Nwalkers);
				endrun(line);
			}
			if (nfilepar != All.Nparam)
			{
				sprintf(line, "The walker file %s contains %d parameters (while the code needs %d)",infname,nfilepar,All.Nparam);
				endrun(line);
			}
			if (nrun < 0)
			{
				sprintf(line, "The walker file %s has a negative run number (%02d)",infname,nrun);
				endrun(line);
			}
			if (istart < 0)
			{
				sprintf(line, "The walker file %s has a negative starting step number (%02d)",infname,istart);
				endrun(line);
			}

			//Increment the Run Number
			nrun++;

			//Increment the starting step
			istart++;

			//Print what is done...
			printf("%s\n%s Restarting HYBRID optimization method...\n",All.fullline,All.startline);
			printf("%s Reading %d Walkers from %s...\n",All.startline,All.Nwalkers,infname);

			//Set output file and open file pointer
			sprintf(outfname,"%s/hybrid%d.%03d.out",All.OutputDir,All.Seed,nrun);
			ofp = fopen(outfname,"w");

			//Get each line in the file
			for (k = 0; k < All.Nwalkers; k++)
			{
				//Get the new line
				fgets(line,nline,ifp);
				tmp = strtok(line, " ");
				//Get the chi^2
				chi2[k] = atof(tmp);
				//Get all the parameters for the walker
				for (i = 0; i < All.Nparam; i++)
				{
					tmp = strtok(NULL, " ");
					x[k*All.Nparam+i] = atof(tmp);
				}
				//Get acceptance counters for each walker
				tmp = strtok(NULL, " ");
				ac[k] = atoi(tmp);
			} //End reading walker file
			fclose(ifp);
		}

		//Broadcast the nrun, istart and the temperature from task 0
		MPI_Bcast(&nrun,   1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&istart, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&All.temperature, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		//Re-initialize random numbers
		gsl_rng_set(rng_uniform_hybrid,  All.Seed + 100 * nrun);
		gsl_rng_set(rng_gaussian_hybrid, All.Seed + 100 * nrun);

		//Broadcast the walker array x from task 0
		MPI_Bcast(x, All.Nwalkers * All.Nparam, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		//Broadcast the acceptance counter array from task 0
		MPI_Bcast(ac, All.Nwalkers, MPI_INT, 0, MPI_COMM_WORLD);
	}

	//Broadcast the chi2 array from task 0
	MPI_Bcast(chi2, All.Nwalkers, MPI_FLOAT, 0, MPI_COMM_WORLD);
	//Broadcast starting state of the chains
	MPI_Bcast(&mean_starting_chi2, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

	//Define loop iteration variable
	istep=istart;
	isave=istart*All.Nwalkers;
	doloop=1;

	//Inform everyone we are starting with the main MCMC loop
	if (ThisTask == 0) printf("%s\n%s Starting main loop for HYBRID optimization method at istart = %d...\n",All.fullline,All.startline,istart);

	//Start main MCMC loop, stop loop after some specified time is reached
	while (doloop)
	{
		//Initialise Loop start time
		lastlooptime = second() - All.starttime;

		if (ThisTask == 0)
		{
			//Print the header for the chain steps
			print_mcmc_header();

			// Compute mean/min/max current chi2
			mean_current_chi2 = 0;
			min_current_chi2 = max_current_chi2 = chi2[0];
			for (k = 0; k < All.Nwalkers; k++)
			{
				mean_current_chi2 += chi2[k];
				if (min_current_chi2 > chi2[k]) min_current_chi2 = chi2[k];
				if (max_current_chi2 < chi2[k]) max_current_chi2 = chi2[k];
			}
			mean_current_chi2 /= All.Nwalkers;
		}

		//Broadcast mean, min and max chi^2 to other tasks
		MPI_Bcast(&mean_current_chi2, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&min_current_chi2,  1, MPI_FLOAT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&max_current_chi2,  1, MPI_FLOAT, 0, MPI_COMM_WORLD);

		// Check how well current sample is doing wrt to initial sample
		q = mean_current_chi2 / mean_starting_chi2;
		g = pow(q,HYBRID_ALPHA);

		// Sets the likelihood of accepting based on the current temperature
		if (All.temperature > 0.0) alpha = log(((double)(istep))+2.71828182845905) / All.temperature;
		else if (All.temperature == 0.0) alpha = 1.e20;
		else alpha = -1.0/All.temperature;

		//Loop through the walkers
		for (j = 0; j < All.Nwalkers/All.NUniverses; j++)
		{
			//Set the current time
			All.timenow = second() - All.starttime;
			//Set up new trial walkers for all universes

			for (itask = 0; itask < All.NUniverses; itask++)
			{
				//Current walker
				iwalker = j*All.NUniverses + itask;
				//Get this walkers previous probability
				cx = chi2[iwalker];
				// Project chi2 squares to p around [-1:1] and compute f
				// If chi2 is larger than mean, p < 0, use linear slope
				if (cx > mean_current_chi2)
				{
					p = ((mean_current_chi2-cx)/(max_current_chi2-mean_current_chi2));
					f = 1.0 - HYBRID_BETA * p;
				} else
				// Now chi2 is smaller than mean, p > 0, and we use power law decay
				{
					p = ((mean_current_chi2-cx)/(mean_current_chi2-min_current_chi2));
					f = pow(p + 1.0, -1.0 * HYBRID_GAMMA);
				}

				//Stretch walker iwalker along jwalker to obtain new propasal walker y
				for (i = 0; i < All.Nparam; i++)
					y[itask*All.Nparam+i] = x[iwalker*All.Nparam+i] + sigma[i] * f * g * gsl_ran_gaussian(rng_gaussian_hybrid,1.0);
			}

			//Broadcast the walker array y from task 0
			MPI_Bcast(y, All.NUniverses * All.Nparam, MPI_DOUBLE, 0, MPI_COMM_WORLD);

			//Compute the probability for each universe
			cy   = -2.0*lnprob(y);
			MPI_Barrier(MPI_COMM_WORLD);

			//Go through all master tasks
			for (itask = 0; itask < All.NUniverses; itask++)
			{
				//Set index of the current walker
				iwalker = j*All.NUniverses + itask;
				//Get this walkers previous probability
				cx = chi2[iwalker];
				//Compute p and f again for this task
				if (cx > mean_current_chi2)
				{
					p = ((mean_current_chi2-cx)/(max_current_chi2-mean_current_chi2));
					f = 1.0 - HYBRID_BETA * p;
				} else
				// Now chi2 is smaller than mean, p > 0, and we use power law decay
				{
					p = ((mean_current_chi2-cx)/(mean_current_chi2-min_current_chi2));
					f = pow(p + 1.0, -1.0 * HYBRID_GAMMA);
				}
				//If this is a master task but not task 0
				if (ThisTask == All.NTaskPerUniverse * itask && itask > 0)
				{
					//Send the chi^2 to task 0
					MPI_Ssend(&cy, 1, MPI_FLOAT, 0, TAG_PROB, MPI_COMM_WORLD);
				}
				//If this is task 0 and it is not the first loop
				if (ThisTask == 0 && itask > 0)
				{
					//Receive the chi^2 from the other master task
					MPI_Recv(&cy, 1, MPI_FLOAT, All.NTaskPerUniverse * itask, TAG_PROB, MPI_COMM_WORLD, &status);
				}

				//Sanity check if the computed chi^2 is in a sensible range
				if (!(cy > 0 && cy < 1.e10) ) cy = 2.e10;

				// Check if new step is accepted
				accept = (alpha * 0.5 * (cx-cy) > log(gsl_rng_uniform(rng_uniform_hybrid)));

				//If the probability is beyond the edge do not accept
				if (cy > 1.e10) accept = 0;

				//If this is task 0 write all parameters and the chi^2 to screen and output file
				if (ThisTask == 0)
				{
					//If y is accepted set walker k to y and store chi^2 and update accepted fraction
					if (accept)
					{
						for (i = 0; i < All.Nparam; i++) x[iwalker*All.Nparam+i] = y[itask*All.Nparam+i];
						cx = cy;
						chi2[iwalker] = cx;
						ac[iwalker]++;
					}

					//Print the walkers to screen
					printf("%8d %8.1f",isave,All.timenow);
					for (i = 0; i < All.Nparam; i++) printf(" %11.8f",x[iwalker*All.Nparam+i]);
					printf(" %f\n",-0.5*chi2[iwalker]);

					//Output to file
					fprintf(ofp,"%f",chi2[iwalker]);
					for (i = 0; i < All.Nparam; i++) fprintf(ofp," %f",x[iwalker*All.Nparam+i]);
					fprintf(ofp," %f %f %f %f %f %f %f %d %d %d\n",mean_current_chi2,p,q,f,g,alpha,((double)(ac[iwalker]))/((double)(istep+1)),isave,istep,iwalker);
					fflush(ofp);

				}//End task 0

				//Update output index
				isave++;

			}//loop through master tasks

			//Broadcast the walker array x from task 0
			MPI_Bcast(x, All.Nwalkers * All.Nparam, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			//Broadcast the chi2 array from task 0
			MPI_Bcast(chi2, All.Nwalkers, MPI_FLOAT, 0, MPI_COMM_WORLD);
			//Broadcast the acceptance counter array from task 0
			MPI_Bcast(ac, All.Nwalkers, MPI_INT, 0, MPI_COMM_WORLD);

		}//loop through walkers

		//Only task 0 writes the walker file
		if (ThisTask == 0)
		{
			//Write all current walkers to a file
			sprintf(outfname,"%s/walkers.hybrid.%d.out",All.OutputDir,All.Seed);
			wfp = fopen(outfname,"w");
			fprintf(wfp,"#Walkers/Parameters/seed/i/mean_starting_chi2/Temperature %d %d %d %d %f %f\n",All.Nwalkers,All.Nparam,nrun,istep,mean_starting_chi2,All.temperature);
			for (iwalker = 0; iwalker < All.Nwalkers; iwalker++)
			{
				fprintf(wfp,"%f",chi2[iwalker]);
				for (i = 0; i < All.Nparam; i++) fprintf(wfp," %f",x[iwalker*All.Nparam+i]);
				fprintf(wfp," %d\n",ac[iwalker]);
			}
			fclose(wfp);
		}

		//Update step
		istep++;

		//Get time the initialisation loop took
		lastlooptime = second() - All.starttime - lastlooptime;
		All.timenow  = second() - All.starttime;

		//End loop if projected finish time of next loop would be larger than the time limit
		if (All.timenow+1.5*lastlooptime>All.timelimit) doloop = 0;

	}//loop

	//Inform everyone we are done with the main MCMC loop
	if (ThisTask == 0) printf("%s Done with HYBRID optimization method...\n",All.startline);

	//Free all arrays and the output file pointer
	efree_movable(y);
	efree_movable(x);
	efree_movable(sigma);
	efree_movable(chi2);
	efree_movable(ac);
	if (ThisTask == 0) fclose(ofp);
	gsl_rng_free(rng_gaussian_hybrid);
	gsl_rng_free(rng_uniform_hybrid);

}


/*! \brief This function runs the parallel tempering method
 *
 *  The method can be either started (restart==0) or restarted from a stored walker file (restart==1).
 *  If newly started it first initialises all walkers with parameters according to the specified range
 *  in the parameter file. Then the probability for each walker is computed for all universes in
 *  parallel. The fraction of cold walkers with a temperature of T=1 is given by PT_COLD_FRACTION. The
 *  temperature for the remaining walkers is taken from a log uniform distribution up to the temperature
 *  given in the parameter file. If restarted all walkers including their parameters, probabilities,
 *  and temperatures are read from a walker file. In the main PT loop first each walker can swap the
 *  position and probability with a random other walker of different temperature. Then each walker is
 *  evolved using the Metropolis-Hastings algorithm, i.e. a new proposal walker is computed given
 *  the current walker and a proposal sigma. The probability of the new walker is then compared to the
 *  current one and the acceptance probability at the given temperature is computed. If accepted the
 *  new position and probability replaces the current one. This is done in parallel for all universes.
 *  During each chain step all walkers are printed to a file. After each step a walker file is written
 *  containing the parameters and probabilities for all walkers. The chain is terminated if the predicted
 *  ending time for the next step is longer then the time limit specified in the parameter file.
 *
 *  \param restart Flag indicating if the PT chain is started (0) or restarted from a walker file (1)
 */
void parallel_tempering(int restart)
{

//First define the fraction of cold walkers with T=1
#ifndef PT_COLD_FRACTION
#define PT_COLD_FRACTION 0.25
#endif

//Now define the target acceptance rate for each walker
#ifndef PT_TARGET_ACCEPTANCE
#define PT_TARGET_ACCEPTANCE 0.30
#endif

//And define the number of chain steps after which the proposal width sigma is adjusted
#ifndef PT_STEPS_SCALE_ADJUST
#define PT_STEPS_SCALE_ADJUST 10
#endif

	const int nline = 500;
	char outfname[NSTRING], infname[NSTRING], line[nline], *tmp;
	int i,j,k,iwalker,istart,itask,istep,isave,doloop,accept,nrun,nfilepar,nfilewalk,ncold;
	float px,py;
	double xswap, lastlooptime, facc;
	FILE *ifp,*ofp,*wfp;
	MPI_Status status;
	gsl_rng *rng_uniform_pt, *rng_gaussian_pt;

	int    *ac;    //Acceptance counter for each walker
	double *x;     //Walkers, Length: Nwalkers*Nparam
	double *y;     //Updated walker
	double *sigma; //Width for the chosen parameter space for each walker
	double *scale; //Scale parameter to adjust sigma based on acceptance rate for each walker
	float  *prob;  //Current chi^2 for each walker
	float  *beta;  //Inverse temperature for each walker

	//The function can only be called with a restart of 0 (start new chain) or 1 (restart from walker file)
	if (restart < 0 || restart > 1)
	{
		if (ThisTask == 0) printf("Parallel tempering method called with wrong restart flag");
		return;
	}

	//Allocate all arrays
	ac    = emalloc_movable(&ac,    "AC",    All.Nwalkers * sizeof(int));
	prob  = emalloc_movable(&prob,  "Prob",  All.Nwalkers * sizeof(float));
	beta  = emalloc_movable(&beta,  "BetaT", All.Nwalkers * sizeof(float));
	scale = emalloc_movable(&beta,  "Scale", All.Nwalkers * sizeof(double));
	sigma = emalloc_movable(&sigma, "Sigma", All.Nparam   * sizeof(double));
	x     = emalloc_movable(&x,     "X",     All.Nwalkers   * All.Nparam * sizeof(double));
	y     = emalloc_movable(&y,     "Y",     All.NUniverses * All.Nparam * sizeof(double));

	//This is the first run
	nrun = 0;

	// Get a sigma (width) for the chosen parameter space
	init_sigma(sigma);

	//Determine the number of cold walkers
	ncold = (int)(PT_COLD_FRACTION * (double)(All.Nwalkers));

	//Initialise the uniform generator
	rng_uniform_pt = gsl_rng_alloc(gsl_rng_ranlxd2);
	gsl_rng_set(rng_uniform_pt, All.Seed);

	//Initialise the gaussian generator
	rng_gaussian_pt = gsl_rng_alloc(gsl_rng_ranlxd2);
	gsl_rng_set(rng_gaussian_pt, All.Seed);

	//If this is the first run of the MCMC then set up the run
	if (restart == 0)
	{
		//Only the main task opens and writes to the output file and prints to screen
		if (ThisTask == 0)
		{
			//Open the MCMC output file
			sprintf(outfname,"%s/pt%d.%03d.out",All.OutputDir,All.Seed,nrun);
			ofp = fopen(outfname,"w");

			//Print what is done...
			printf("%s\n%s Starting parallel tempering method...\n",All.fullline,All.startline);
			printf("%s Setting up %d Walkers...\n",All.startline,All.Nwalkers);

			//Print the header for the chain steps
			print_mcmc_header();
		}

		//Set the temperature for each walker
		for (k = 0; k < All.Nwalkers; k++)
		{
			//First the cold walkers
			if (k < ncold) beta[k] = 1.0;
			//And the temperature ladder
			else beta[k] = 1.0 / exp((double)(k-ncold+1) / (double)(All.Nwalkers-ncold) * log(All.temperature));
		}

		//Initialization of walkers
		init_param(x);

		//Calculate the starting probability for each walker
		for (j = 0; j < All.Nwalkers/All.NUniverses; j++)
		{
			//Compute the probability for each universe
			px   = lnprob(&x[j*All.NUniverses*All.Nparam]);
			//Set the current time
			All.timenow = second() - All.starttime;
			//Go through all master tasks
			for (itask = 0; itask < All.NUniverses; itask++)
			{
				//If this is a master task but not task 0
				if (ThisTask == All.NTaskPerUniverse * itask && itask > 0)
				{
					//Send the probability to task 0
					MPI_Ssend(&px, 1, MPI_FLOAT, 0, TAG_PROB, MPI_COMM_WORLD);
				}
				//If this is task 0 and it is not the first loop
				if (ThisTask == 0 && itask > 0)
				{
					//Receive the probability from the other master task
					MPI_Recv(&px, 1, MPI_FLOAT, All.NTaskPerUniverse * itask, TAG_PROB, MPI_COMM_WORLD, &status);
				}
				//Write the probability to the correct index (only valid for task 0)
				prob[j*All.NUniverses+itask] = px;
				//If this is task 0 write all parameters and the probability to screen and output file
				if (ThisTask == 0)
				{
					printf("%8d %8.1f",j*All.NUniverses+itask,All.timenow);
					for (i = 0; i < All.Nparam; i++) printf(" %11.8f",x[(j*All.NUniverses+itask)*All.Nparam+i]);
					printf(" %f %6.3f 1.0\n",prob[j*All.NUniverses+itask],1./beta[j*All.NUniverses+itask]);

					fprintf(ofp,"%f",prob[j*All.NUniverses+itask]);
					for (i = 0; i < All.Nparam; i++) fprintf(ofp," %f",x[(j*All.NUniverses+itask)*All.Nparam+i]);
					fprintf(ofp," %6.3f 1.0 0 0 0 %d\n",1./beta[j*All.NUniverses+itask],j*All.NUniverses+itask);
					fflush(ofp);
				}
			}//End loop through master tasks
		}//End loop through all walkers

		// Initialise acceptance counters, scale adjustments and starting step
		for (k = 0; k < All.Nwalkers; k++)
		{
			ac[k]    = 0;
			scale[k] = 1.0;
		}
		istart = 0;

	}//End restart 0

	//If the chain is restarted from a file with walkers
	if (restart == 1)
	{
		//Only the main task opens and writes to the output file and prints to screen
		if (ThisTask == 0)
		{
			//Open the walker file
			sprintf(infname,"%s/walkers.pt.%d.out",All.OutputDir,All.Seed);
			if(!(ifp=fopen(infname,"r")))
			{
				sprintf(line, "Cannot open file `%s' for reading the walkers.", infname);
				endrun(line);
			}

			//Get first line of walker file
			fgets(line,nline,ifp);
			sscanf(line,"#Walkers/Parameters/seed/i %d %d %d %d\n", &nfilewalk, &nfilepar, &nrun, &istart);

			//Sanity Checks
			if (nfilewalk != All.Nwalkers)
			{
				sprintf(line, "The walker file %s contains %d walkers (while the code needs %d)",infname,nfilewalk,All.Nwalkers);
				endrun(line);
			}
			if (nfilepar != All.Nparam)
			{
				sprintf(line, "The walker file %s contains %d parameters (while the code needs %d)",infname,nfilepar,All.Nparam);
				endrun(line);
			}
			if (nrun < 0)
			{
				sprintf(line, "The walker file %s has a negative run number (%02d)",infname,nrun);
				endrun(line);
			}
			if (istart < 0)
			{
				sprintf(line, "The walker file %s has a negative starting step number (%02d)",infname,istart);
				endrun(line);
			}

			//Increment the Run Number
			nrun++;

			//Increment the starting step
			istart++;

			//Print what is done...
			printf("%s\n%s Restarting parallel tempering method...\n",All.fullline,All.startline);
			printf("%s Reading %d Walkers from %s...\n",All.startline,All.Nwalkers,infname);

			//Set output file and open file pointer
			sprintf(outfname,"%s/pt%d.%03d.out",All.OutputDir,All.Seed,nrun);
			ofp = fopen(outfname,"w");

			//Get each line in the file
			for (k = 0; k < All.Nwalkers; k++)
			{
				//Get the new line
				fgets(line,nline,ifp);
				tmp = strtok(line, " ");
				//Get the probability
				prob[k] = atof(tmp);
				//Get beta
				tmp = strtok(NULL, " ");
				beta[k] = 1./atof(tmp);
				//Get scale
				tmp = strtok(NULL, " ");
				scale[k] = atof(tmp);
				//Get acceptance counters for each walker
				tmp = strtok(NULL, " ");
				ac[k] = atoi(tmp);
				//Get all the parameters for the walker
				for (i = 0; i < All.Nparam; i++)
				{
					tmp = strtok(NULL, " ");
					x[k*All.Nparam+i] = atof(tmp);
				}
			} //End reading walker file
			fclose(ifp);
		}

		//Broadcast the nrun and istart from task 0
		MPI_Bcast(&nrun,   1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&istart, 1, MPI_INT, 0, MPI_COMM_WORLD);

		//Re-initialize random numbers
		gsl_rng_set(rng_uniform_pt,  All.Seed + 100 * nrun);
		gsl_rng_set(rng_gaussian_pt, All.Seed + 100 * nrun);

		//Broadcast the walker array x from task 0
		MPI_Bcast(x, All.Nwalkers * All.Nparam, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		//Broadcast the scale array from task 0
		MPI_Bcast(scale, All.Nwalkers, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		//Broadcast the scale array from task 0
		MPI_Bcast(beta, All.Nwalkers, MPI_FLOAT, 0, MPI_COMM_WORLD);
		//Broadcast the acceptance counter array from task 0
		MPI_Bcast(ac, All.Nwalkers, MPI_INT, 0, MPI_COMM_WORLD);
	}

	//Broadcast the probability array from task 0
	MPI_Bcast(prob, All.Nwalkers, MPI_FLOAT, 0, MPI_COMM_WORLD);

	//Define loop iteration variable
	istep=istart;
	isave=istart*All.Nwalkers;
	doloop=1;

	//Inform everyone we are starting with the main MCMC loop
	if (ThisTask == 0) printf("%s\n%s Starting main loop for parallel tempering at istart = %d...\n",All.fullline,All.startline,istart);

	//Start main MCMC loop, stop loop after some specified time is reached
	while (doloop)
	{
		//Initialise Loop start time
		lastlooptime = second() - All.starttime;

		//After every PT_STEPS_SIGMA_ADJUST chain step we adjust the scale for each walker
		if (istep % PT_STEPS_SCALE_ADJUST == 0 && istep > 0)
		{
			for (k = 0; k < All.Nwalkers; k++)
			{
				facc = ((double)(ac[k]))/((double)(istep));
				//If the acceptance fraction is smaller than the target rate make the proposal step smaller
				if ( facc < PT_TARGET_ACCEPTANCE)
					scale[k] /= 1.0 + (PT_TARGET_ACCEPTANCE-facc) / PT_TARGET_ACCEPTANCE;
				//If the acceptance fraction is higher than the target rate make the proposal step larger
				if ( facc > PT_TARGET_ACCEPTANCE)
					scale[k] *= 1.0 + (facc-PT_TARGET_ACCEPTANCE) / (1.0-PT_TARGET_ACCEPTANCE);
			}
		}

		//Loop over all walkers and check if they are swapped for another walker
		for (k = 0; k < All.Nwalkers; k++)
		{
			//Randomly choose walker for swapping
			do j = (int)(gsl_rng_uniform(rng_uniform_pt)*(double)(All.Nwalkers));
			while (j==k);

			// Check if the two walkers are swapped (positions and probability)
			accept = ((beta[k]-beta[j]) * (prob[j]-prob[k]) > log(gsl_rng_uniform(rng_uniform_pt)));

			//Do not swap if the temperatures are equal
			if (beta[k] == beta[j]) accept = 0;

			//If the swap is accepted
			if (accept)
			{
				//Print the swapping to screen
				if (ThisTask==0) printf("%s Swapping walker %d (T = %f, lnprob=%f) and %d (T = %f, lnprob=%f)\n",All.startline,k,1./beta[k],prob[k],j,1./beta[j],prob[j]);
				//Swap walker k and walker j
				for (i = 0; i < All.Nparam; i++)
				{
					//parameters
					xswap                 = x[k*All.Nparam+i];
					x[k*All.Nparam+i]     = x[j*All.Nparam+i];
					x[j*All.Nparam+i]     = xswap;
				}
				//probability
				px      = prob[k];
				prob[k] = prob[j];
				prob[j] = px;
			}//End accept
		}//End loop over walkers

		//Print the header for the chain steps
		if (ThisTask == 0) print_mcmc_header();

		//Loop through the walkers
		for (j = 0; j < All.Nwalkers/All.NUniverses; j++)
		{
			//Set the current time
			All.timenow = second() - All.starttime;

			//Set up new trial walkers for all universes
			for (itask = 0; itask < All.NUniverses; itask++)
			{
				//Current walker
				iwalker = j*All.NUniverses + itask;
				//Get this walkers previous probability
				px = prob[iwalker];
				//Propose new walker
				for (i = 0; i < All.Nparam; i++)
					y[itask*All.Nparam+i] = x[iwalker*All.Nparam+i] + sigma[i] * scale[iwalker] * gsl_ran_gaussian(rng_gaussian_pt,1.0);
			}

			//Broadcast the walker array y from task 0
			MPI_Bcast(y, All.NUniverses * All.Nparam, MPI_DOUBLE, 0, MPI_COMM_WORLD);

			//Compute the probability for each universe
			py = lnprob(y);
			MPI_Barrier(MPI_COMM_WORLD);

			//Go through all master tasks
			for (itask = 0; itask < All.NUniverses; itask++)
			{
				//Set index of the current walker
				iwalker = j*All.NUniverses + itask;
				//Get this walkers previous probability
				px = prob[iwalker];
				//If this is a master task but not task 0
				if (ThisTask == All.NTaskPerUniverse * itask && itask > 0)
				{
					//Send the probability to task 0
					MPI_Ssend(&py, 1, MPI_FLOAT, 0, TAG_PROB, MPI_COMM_WORLD);
				}
				//If this is task 0 and it is not the first loop
				if (ThisTask == 0 && itask > 0)
				{
					//Receive the probability from the other master task
					MPI_Recv(&py, 1, MPI_FLOAT, All.NTaskPerUniverse * itask, TAG_PROB, MPI_COMM_WORLD, &status);
				}

				//Sanity check if the computed probability is in a sensible range
				if (!(py < 0 && py > -1.e10) ) py = -2.e10;

				//Check if new step is accepted
				accept = (beta[iwalker] * (py-px) > log(gsl_rng_uniform(rng_uniform_pt)));

				//If the probability is beyond the edge do not accept
				if (py < -1.e10) accept = 0;

				//If this is task 0 write all parameters and the probability to screen and output file
				if (ThisTask == 0)
				{
					//If y is accepted set walker k to y and store probability and update accepted fraction
					if (accept)
					{
						for (i = 0; i < All.Nparam; i++) x[iwalker*All.Nparam+i] = y[itask*All.Nparam+i];
						px = py;
						prob[iwalker] = px;
						ac[iwalker]++;
					}

					//Print the walkers to screen
					printf("%8d %8.1f",isave,All.timenow);
					for (i = 0; i < All.Nparam; i++) printf(" %11.8f",x[iwalker*All.Nparam+i]);
					printf(" %f %6.3f %6.3f\n",prob[iwalker],1./beta[iwalker],scale[iwalker]);

					//Output to file
					fprintf(ofp,"%f",prob[iwalker]);
					for (i = 0; i < All.Nparam; i++) fprintf(ofp," %f",x[iwalker*All.Nparam+i]);
					fprintf(ofp," %f %f %f %d %d %d\n",1./beta[iwalker],scale[iwalker],((double)(ac[iwalker]))/((double)(istep+1)),isave,istep,iwalker);
					fflush(ofp);

				}//End task 0

				//Update output index
				isave++;

			}//loop through master tasks

			//Broadcast the walker array x from task 0
			MPI_Bcast(x, All.Nwalkers * All.Nparam, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			//Broadcast the probability array from task 0
			MPI_Bcast(prob, All.Nwalkers, MPI_FLOAT, 0, MPI_COMM_WORLD);
			//Broadcast the acceptance counter array from task 0
			MPI_Bcast(ac, All.Nwalkers, MPI_INT, 0, MPI_COMM_WORLD);

		}//loop through walkers

		//Only task 0 writes the walker file
		if (ThisTask == 0)
		{
			//Write all current walkers to a file
			sprintf(outfname,"%s/walkers.pt.%d.out",All.OutputDir,All.Seed);
			wfp = fopen(outfname,"w");
			fprintf(wfp,"#Walkers/Parameters/seed/i %d %d %d %d\n",All.Nwalkers,All.Nparam,nrun,istep);
			for (iwalker = 0; iwalker < All.Nwalkers; iwalker++)
			{
				fprintf(wfp,"%f %f %f %d",prob[iwalker],1./beta[iwalker],scale[iwalker],ac[iwalker]);
				for (i = 0; i < All.Nparam; i++) fprintf(wfp," %f",x[iwalker*All.Nparam+i]);
				fprintf(wfp,"\n");
			}
			fclose(wfp);
		}

		//Update step
		istep++;

		//Get time the initialisation loop took
		lastlooptime = second() - All.starttime - lastlooptime;
		All.timenow  = second() - All.starttime;

		//End loop if projected finish time of next loop would be larger than the time limit
		if (All.timenow+1.5*lastlooptime>All.timelimit) doloop = 0;

	}//loop


	//Inform everyone we are done with the main MCMC loop
	if (ThisTask == 0) printf("%s Done with parallel tempering method...\n",All.startline);

	//Free all arrays and the output file pointer
	efree_movable(y);
	efree_movable(x);
	efree_movable(sigma);
	efree_movable(scale);
	efree_movable(beta);
	efree_movable(prob);
	efree_movable(ac);
	if (ThisTask == 0) fclose(ofp);
	gsl_rng_free(rng_gaussian_pt);
	gsl_rng_free(rng_uniform_pt);

}


/*! \brief This function prints a header line for the MCMC
 *
 *  For each parameter the corresponding parameter name is printed with the correct indention.
 */
void print_mcmc_header(void)
{
	if (ThisTask > 0) return;
	int i = 1;
	printf("%s Step Time     ",All.startline);
	if (All.DeltaM0 >= 0.0) printf("(%02d)M0      ",i);
	if (All.DeltaM0 >= 0.0) i++;
	if (All.DeltaEpsilon0 >= 0.0) printf("(%02d)E0      ",i);
	if (All.DeltaEpsilon0 >= 0.0) i++;
	if (All.DeltaBeta0 >= 0.0) printf("(%02d)B0      ",i);
	if (All.DeltaBeta0 >= 0.0) i++;
#ifndef SFE_SAME_SLOPE
	if (All.DeltaGamma0 >= 0.0) printf("(%02d)G0      ",i);
	if (All.DeltaGamma0 >= 0.0) i++;
#endif
#ifdef SFE_MPEAK_ZEVOLV
	if (All.DeltaMZ >= 0.0) printf("(%02d)MZ      ",i);
	if (All.DeltaMZ >= 0.0) i++;
#endif
#ifdef SFE_NORM_ZEVOLV
	if (All.DeltaEpsilonZ >= 0.0) printf("(%02d)EZ      ",i);
	if (All.DeltaEpsilonZ >= 0.0) i++;
#endif
#ifdef SFE_BETA_ZEVOLV
	if (All.DeltaBetaZ >= 0.0) printf("(%02d)BZ      ",i);
	if (All.DeltaBetaZ >= 0.0) i++;
#endif
#ifdef SFE_GAMMA_ZEVOLV
#ifndef SFE_SAME_SLOPE
	if (All.DeltaGammaZ >= 0.0) printf("(%02d)GZ      ",i);
	if (All.DeltaGammaZ >= 0.0) i++;
#endif
#endif
	if (All.DeltaFesc >= 0.0) printf("(%02d)Fesc    ",i);
	if (All.DeltaFesc >= 0.0) i++;
	if (All.DeltaFstrip >= 0.0) printf("(%02d)Fstrip  ",i);
	if (All.DeltaFstrip >= 0.0) i++;
	if (All.DeltaTau0 >= 0.0) printf("(%02d)Tau0    ",i);
	if (All.DeltaTau0 >= 0.0) i++;
#ifdef SAT_SFR_EXP_DECAY
	if (All.DeltaTauD >= 0.0) printf("(%02d)TauD    ",i);
	if (All.DeltaTauD >= 0.0) i++;
#endif
#ifdef SAT_QUENCH_MASS_DEPENDENT
	if (All.DeltaTauS >= 0.0) printf("(%02d)TauS    ",i);
	if (All.DeltaTauS >= 0.0) i++;
#endif
	printf("(%02d)lnprob",i);
	i++;
	printf("\n");
}


/*! \brief This function initialises a walker with random parameters within the specified range
 *
 *  This function s initialises a walker p with parameters within a range that is specified in
 *  the parameter file. Within this range the parameters are uniformly sampled.
 *
 *  \param *p Walker array that will be initialised
 */
void init_param(double *p)
{
	int i,k;

	//Initialise the uniform generator
	gsl_rng *rng_uniform_init = gsl_rng_alloc(gsl_rng_ranlxd2);
	gsl_rng_set(rng_uniform_init, All.Seed + 661894);

	for (k = 0; k < All.Nwalkers; k++)
	{
		i = 0;

		//M0
		if (All.DeltaM0 >= 0.0) p[k*All.Nparam + i] = (P.M0-All.DeltaM0/2.0) + gsl_rng_uniform(rng_uniform_init)*All.DeltaM0;
		if (All.DeltaM0 >= 0.0) i++;

		//E0
		if (All.DeltaEpsilon0 >= 0.0) p[k*All.Nparam + i] = (log10(P.Epsilon0)-All.DeltaEpsilon0/2.0) + gsl_rng_uniform(rng_uniform_init)*All.DeltaEpsilon0;
		if (All.DeltaEpsilon0 >= 0.0) i++;

		//B0
		if (All.DeltaBeta0 >= 0.0) p[k*All.Nparam + i] = (P.Beta0-All.DeltaBeta0/2.0) + gsl_rng_uniform(rng_uniform_init)*All.DeltaBeta0;
		if (All.DeltaBeta0 >= 0.0) i++;

		//G0
#ifndef SFE_SAME_SLOPE
		if (All.DeltaGamma0 >= 0.0) p[k*All.Nparam + i] = (P.Gamma0-All.DeltaGamma0/2.0) + gsl_rng_uniform(rng_uniform_init)*All.DeltaGamma0;
		if (All.DeltaGamma0 >= 0.0) i++;
#endif

		//MZ
#ifdef SFE_MPEAK_ZEVOLV
		if (All.DeltaMZ >= 0.0) p[k*All.Nparam + i] = (P.MZ+P.M0-p[k*All.Nparam]-All.DeltaMZ/2.0) + gsl_rng_uniform(rng_uniform_init)*All.DeltaMZ;
		if (All.DeltaMZ >= 0.0) i++;
#endif

		//EZ
#ifdef SFE_NORM_ZEVOLV
		if (All.DeltaEpsilonZ >= 0.0) p[k*All.Nparam + i] = (P.EpsilonZ-All.DeltaEpsilonZ/2.0) + gsl_rng_uniform(rng_uniform_init)*All.DeltaEpsilonZ;
		if (All.DeltaEpsilonZ >= 0.0) i++;
#endif

		//BZ
#ifdef SFE_BETA_ZEVOLV
		if (All.DeltaBetaZ >= 0.0) p[k*All.Nparam + i] = (P.BetaZ+P.Beta0-p[k*All.Nparam+2]-All.DeltaBetaZ/2.0) + gsl_rng_uniform(rng_uniform_init)*All.DeltaBetaZ;
		if (All.DeltaBetaZ >= 0.0) i++;
#endif

		//GZ
#ifdef SFE_GAMMA_ZEVOLV
#ifndef SFE_SAME_SLOPE
		if (All.DeltaGammaZ >= 0.0) p[k*All.Nparam + i] = (P.GammaZ-All.DeltaGammaZ/2.0) + gsl_rng_uniform(rng_uniform_init)*All.DeltaGammaZ;
		if (All.DeltaGammaZ >= 0.0) i++;
#endif
#endif

		//Fesc
		if (All.DeltaFesc >= 0.0) p[k*All.Nparam + i] = (log10(P.Fesc)-All.DeltaFesc/2.0) + gsl_rng_uniform(rng_uniform_init)*All.DeltaFesc;
		if (All.DeltaFesc >= 0.0) i++;

		//Fstrip
		if (All.DeltaFstrip >= 0.0) p[k*All.Nparam + i] = (log10(P.Fstrip)-All.DeltaFstrip/2.0) + gsl_rng_uniform(rng_uniform_init)*All.DeltaFstrip;
		if (All.DeltaFstrip >= 0.0) i++;

		//Tau0
		if (All.DeltaTau0 >= 0.0) p[k*All.Nparam + i] = (log10(P.Tau0)-All.DeltaTau0/2.0) + gsl_rng_uniform(rng_uniform_init)*All.DeltaTau0;
		if (All.DeltaTau0 >= 0.0) i++;

		//TauD
#ifdef SAT_SFR_EXP_DECAY
		if (All.DeltaTauD >= 0.0) p[k*All.Nparam + i] = (P.TauD-All.DeltaTauD/2.0) + gsl_rng_uniform(rng_uniform_init)*All.DeltaTauD;
		if (All.DeltaTauD >= 0.0) i++;
#endif

		//TauS
#ifdef SAT_QUENCH_MASS_DEPENDENT
		if (All.DeltaTauS >= 0.0) p[k*All.Nparam + i] = (P.TauS-All.DeltaTauS/2.0) + gsl_rng_uniform(rng_uniform_init)*All.DeltaTauS;
		if (All.DeltaTauS >= 0.0) i++;
#endif

	}

	gsl_rng_free(rng_uniform_init);

}


/*! \brief This function initialises the standard deviation with the specified range
 *
 *  This function s initialises the standard deviation s with values as specified in
 *  the parameter file. The range is modified with the MCMC scale parameter.
 *
 *  \param *s Standard deviation array that will be initialised
 */
void init_sigma(double *s)
{
	int i = 0;

	//M0
	if (All.DeltaM0 >= 0.0) s[i] = All.DeltaM0/pow(2.0,All.mcmca);
	if (All.DeltaM0 >= 0.0) i++;

	//E0
	if (All.DeltaEpsilon0 >= 0.0) s[i] = All.DeltaEpsilon0/pow(2.0,All.mcmca);
	if (All.DeltaEpsilon0 >= 0.0) i++;

	//B0
	if (All.DeltaBeta0 >= 0.0) s[i] = All.DeltaBeta0/pow(2.0,All.mcmca);
	if (All.DeltaBeta0 >= 0.0) i++;

	//G0
#ifndef SFE_SAME_SLOPE
	if (All.DeltaGamma0 >= 0.0) s[i] = All.DeltaGamma0/pow(2.0,All.mcmca);
	if (All.DeltaGamma0 >= 0.0) i++;
#endif

	//MZ
#ifdef SFE_MPEAK_ZEVOLV
	if (All.DeltaMZ >= 0.0) s[i] = All.DeltaMZ/pow(2.0,All.mcmca);
	if (All.DeltaMZ >= 0.0) i++;
#endif

	//EZ
#ifdef SFE_NORM_ZEVOLV
	if (All.DeltaEpsilonZ >= 0.0) s[i] = All.DeltaEpsilonZ/pow(2.0,All.mcmca);
	if (All.DeltaEpsilonZ >= 0.0) i++;
#endif

	//BZ
#ifdef SFE_BETA_ZEVOLV
	if (All.DeltaBetaZ >= 0.0) s[i] = All.DeltaBetaZ/pow(2.0,All.mcmca);
	if (All.DeltaBetaZ >= 0.0) i++;
#endif

	//GZ
#ifdef SFE_GAMMA_ZEVOLV
#ifndef SFE_SAME_SLOPE
	if (All.DeltaGammaZ >= 0.0) s[i] = All.DeltaGammaZ/pow(2.0,All.mcmca);
	if (All.DeltaGammaZ >= 0.0) i++;
#endif
#endif

	//Fesc
	if (All.DeltaFesc >= 0.0) s[i] = All.DeltaFesc/pow(2.0,All.mcmca);
	if (All.DeltaFesc >= 0.0) i++;

	//Fstrip
	if (All.DeltaFstrip >= 0.0) s[i] = All.DeltaFstrip/pow(2.0,All.mcmca);
	if (All.DeltaFstrip >= 0.0) i++;

	//Tau0
	if (All.DeltaTau0 >= 0.0) s[i] = All.DeltaTau0/pow(2.0,All.mcmca);
	if (All.DeltaTau0 >= 0.0) i++;

	//TauD
#ifdef SAT_SFR_EXP_DECAY
	if (All.DeltaTauD >= 0.0) s[i] = All.DeltaTauD/pow(2.0,All.mcmca);
	if (All.DeltaTauD >= 0.0) i++;
#endif

	//TauS
#ifdef SAT_QUENCH_MASS_DEPENDENT
	if (All.DeltaTauS >= 0.0) s[i] = All.DeltaTauS/pow(2.0,All.mcmca);
	if (All.DeltaTauS >= 0.0) i++;
#endif

}


/*! \brief This function computes the log probability for the walkers on each master task
 *
 *  This function first sets the values of the local parameter structure P to the specified values
 *  in the array p that consists of NUniverse walkers (one for each universe that is computed in parallel).
 *  Then a few basic physical boundary conditions are check and if not met a very large value is returned.
 *  Otherwise the galaxies are assembled for each universe and \f$\chi^2\f$ for each universe is computed.
 *  Each master task returns its respective \f$-\chi^2/2\f$ (all other tasks return 0).
 *
 *  \param *p Array with NUniverse walkers that will be used to compute the probability
 */
float lnprob(double *p)
{

	int i = 0;
	int thisUniverse = MasterTask / All.NTaskPerUniverse;
	float chi2 = 0.0;

	// Initialize the total DoF
	DoF      = 0;
	DoFsmf   = 0;
	DoFfq    = 0;
	DoFcsfrd = 0;
	DoFssfr  = 0;
	DoFwp    = 0;

	// Update all the parameters for each task / universe

	//M0
	if (All.DeltaM0 >= 0.0) P.M0 = p[i + All.Nparam * thisUniverse];
	if (All.DeltaM0 >= 0.0) i++;

	//E0
	if (All.DeltaEpsilon0 >= 0.0) P.Epsilon0 = pow(10.,p[i + All.Nparam * thisUniverse]);
	if (All.DeltaEpsilon0 >= 0.0) i++;

	//B0
	if (All.DeltaBeta0 >= 0.0) P.Beta0 = p[i + All.Nparam * thisUniverse];
	if (All.DeltaBeta0 >= 0.0) i++;

	//G0
#ifndef SFE_SAME_SLOPE
	if (All.DeltaGamma0 >= 0.0) P.Gamma0 = p[i + All.Nparam * thisUniverse];
	if (All.DeltaGamma0 >= 0.0) i++;
#else
	P.Gamma0 = P.Beta0;
#endif

	//MZ
#ifdef SFE_MPEAK_ZEVOLV
	if (All.DeltaMZ >= 0.0) P.MZ = p[i + All.Nparam * thisUniverse];
	if (All.DeltaMZ >= 0.0) i++;
#else
	P.MZ = 0.0;
#endif

	//EZ
#ifdef SFE_NORM_ZEVOLV
	if (All.DeltaEpsilonZ >= 0.0) P.EpsilonZ = p[i + All.Nparam * thisUniverse];
	if (All.DeltaEpsilonZ >= 0.0) i++;
#else
	P.EpsilonZ = 0.0;
#endif

	//BZ
#ifdef SFE_BETA_ZEVOLV
	if (All.DeltaBetaZ >= 0.0) P.BetaZ = p[i + All.Nparam * thisUniverse];
	if (All.DeltaBetaZ >= 0.0) i++;
#else
	P.BetaZ = 0.0;
#endif

	//GZ
#ifdef SFE_GAMMA_ZEVOLV
#ifndef SFE_SAME_SLOPE
	if (All.DeltaGamma0 >= 0.0) P.GammaZ = p[i + All.Nparam * thisUniverse];
	if (All.DeltaGamma0 >= 0.0) i++;
#else
	P.GammaZ = P.BetaZ;
#endif
#else
	P.GammaZ = 0.0;
#endif

	//Fesc
	if (All.DeltaFesc >= 0.0) P.Fesc = pow(10.,p[i + All.Nparam * thisUniverse]);
	if (All.DeltaFesc >= 0.0) i++;

	//Fstrip
	if (All.DeltaFstrip >= 0.0) P.Fstrip = pow(10.,p[i + All.Nparam * thisUniverse]);
	if (All.DeltaFstrip >= 0.0) i++;

	//Tau0
	if (All.DeltaTau0 >= 0.0) P.Tau0 = pow(10.,p[i + All.Nparam * thisUniverse]);
	if (All.DeltaTau0 >= 0.0) i++;

	//TauD
#ifdef SAT_SFR_EXP_DECAY
	if (All.DeltaTauD >= 0.0) P.TauD = p[i + All.Nparam * thisUniverse];
	if (All.DeltaTauD >= 0.0) i++;
#else
	P.TauD = 0.0;
#endif

	//TauS
#ifdef SAT_QUENCH_MASS_DEPENDENT
	if (All.DeltaTauS >= 0.0) P.TauS = p[i + All.Nparam * thisUniverse];
	if (All.DeltaTauS >= 0.0) i++;
#else
	P.TauS = 0.0;
#endif

	//Physical boundary conditions
	if (P.Fesc > 1.0) return -1.e10;
#ifndef SAT_STRIP_USE_MSTAR
	if (P.Fstrip > 1.0) return -1.e10;
#endif
	if (P.Epsilon0+P.EpsilonZ <= 0.0) return -1.e10;
	if (P.M0+P.MZ <= 0.0) return -1.e10;

	//Populate haloes with galaxies
	make_galaxies();

	//Get statistics
	get_statistics();

	//Calculate Chi^2 for this master task
#ifdef READ_SMF
	chi2 += chi2_smf();
#endif

#ifdef READ_FQ
	chi2 += chi2_fq();
#endif

#ifdef READ_CSFRD
	chi2 += chi2_csfrd();
#endif

#ifdef READ_SSFR
	chi2 += chi2_ssfr();
#endif

#ifdef READ_WP
	chi2 += chi2_wp();
#endif

	return -chi2/2.0;

}


/*! \brief This function gets the total \f$\chi^2\f$ and prints it to screen.
 *
 *  This function gets the contribution to \f$\chi^2\f$ from each observable. All contributions and the total are
 *  printed to screen.
 */
float print_chi2(void)
{
	int i;
	float *chi2;

	//Allocate the chi2 array for All.Nobs values plus the total value
	chi2 = emalloc("CHI2", (All.Nobs + 1) * sizeof(float));
	//Compute Chi2
	get_chi2(chi2);
	//Init
	i = 0;

	//Only the main task gets the chi^2 and prints to screen
	if (ThisTask == 0)
	{
		//Print to screen
		printf("%s\n",All.fullline);

		//Print each contribution to the total chi^2
#ifdef READ_SMF
		printf("%s Chi2 from SMF                 = %12.5f   -> Chi2/DoF = %12.5f\n",All.startline,chi2[i],chi2[i]/DoFsmf);
		i++;
#endif

#ifdef READ_FQ
		printf("%s Chi2 from Quenched Fractions  = %12.5f   -> Chi2/DoF = %12.5f\n",All.startline,chi2[i],chi2[i]/DoFfq);
		i++;
#endif

#ifdef READ_CSFRD
		printf("%s Chi2 from Cosmic SFR Density  = %12.5f   -> Chi2/DoF = %12.5f\n",All.startline,chi2[i],chi2[i]/DoFcsfrd);
		i++;
#endif

#ifdef READ_SSFR
		printf("%s Chi2 from Specific SFRs       = %12.5f   -> Chi2/DoF = %12.5f\n",All.startline,chi2[i],chi2[i]/DoFssfr);
		i++;
#endif

#ifdef READ_WP
		printf("%s Chi2 from Clustering          = %12.5f   -> Chi2/DoF = %12.5f\n",All.startline,chi2[i],chi2[i]/DoFwp);
		i++;
#endif

		//Print total chi^2
		printf("%s Chi2 (Total)                  = %12.5f   -> Chi2/DoF = %12.5f\n",All.startline,chi2[i],chi2[i]/DoF);

	}

	//Free the chi2 array
	efree(chi2);

	return chi2[i];

}


/*! \brief This function computes the total \f$\chi^2\f$ and returns it
 *
 *  This function initialises the total \f$\chi^2\f$ and the total and partial degrees of freedom to zero
 *  and computes the contribution to \f$\chi^2\f$ from each observable. All contributions and the total are
 *  returned as a pointer.
 */
void get_chi2(float *chi2)
{
	int i;
	float chi2now,chi2tot;

	//Initialise everything to zero
	chi2tot  = 0.0;
	DoF      = 0;
	DoFsmf   = 0;
	DoFfq    = 0;
	DoFcsfrd = 0;
	DoFssfr  = 0;
	DoFwp    = 0;
	i        = 0;

	//Only the master task computes chi^2
	if (ThisTask == MasterTask)
	{
#ifdef READ_SMF
		chi2now      = chi2_smf();
		chi2tot     += chi2now;
		chi2[i]      = chi2now;
		i++;
#endif

#ifdef READ_FQ
		chi2now      = chi2_fq();
		chi2tot     += chi2now;
		chi2[i]      = chi2now;
		i++;
#endif

#ifdef READ_CSFRD
		chi2now      = chi2_csfrd();
		chi2tot     += chi2now;
		chi2[i]      = chi2now;
		i++;
#endif

#ifdef READ_SSFR
		chi2now      = chi2_ssfr();
		chi2tot     += chi2now;
		chi2[i]      = chi2now;
		i++;
#endif

#ifdef READ_WP
		chi2now      = chi2_wp();
		chi2tot     += chi2now;
		chi2[i]      = chi2now;
		i++;
#endif

		//Get total chi^2
		chi2[i] = chi2tot;

	}

}


/*! \brief This function computes the contribution to chi^2 from the SMF
 *
 *  This function loops over all data points of the SMFs in the interpolation range and sums
 *  their contribution. The total squared error used is the sum of the squared error of the observations
 *  and the model (Poisson). Each data point is weighted by its nfit value. The total and partial
 *  number of degrees of freedom is incremented as well.
 *  The returned \f$\chi^2\f$ is given by \f[\chi^2 = \sum_i n_\mathrm{fit,i}
 *  \frac{(\Phi_\mathrm{obs,i}-\Phi_\mathrm{mod,i})^2}{\sigma_\mathrm{obs,i}^2+\sigma_\mathrm{mod,i}^2} .\f]
 */
float chi2_smf(void)
{
	int i;
	float mod, obs, sigmamod, sigmaobs;
	float chi2 = 0.0;

	//If this is not a master task return 0
	if (ThisTask != MasterTask) return 0.0;

	//Go through all data points of the SMF
	for (i = 0; i < All.Nsmf; i++)
	{
		//If the data point is within the interpolation range add it to the chi^2
		if (Smf[i].mod_y > 0.0 && Smf[i].nfit > 0)
		{
			obs      = Smf[i].obs_y;
			mod      = log10(Smf[i].mod_y);
			sigmaobs = Smf[i].obs_sigma;
			sigmamod = log10(Smf[i].mod_sigma);
			chi2    += (obs-mod) * (obs-mod) / (sigmaobs*sigmaobs + sigmamod*sigmamod) * (float)(Smf[i].nfit);
			DoF     += Smf[i].nfit;
			DoFsmf  += Smf[i].nfit;
		}
	}

	return chi2;
}


/*! \brief This function computes the contribution to \f$\chi^2\f$ from the quenched fractions
 *
 *  This function loops over all data points of the quenched fractions in the interpolation range and sums
 *  their contribution. The total squared error used is the sum of the squared error of the observations
 *  and the model (Poisson). Each data point is weighted by its nfit value. The total and partial
 *  number of degrees of freedom is incremented as well.
 *  The returned \f$\chi^2\f$ is given by \f[\chi^2 = \sum_i n_\mathrm{fit,i}
 *  \frac{(f_\mathrm{q,obs,i}-f_\mathrm{q,mod,i})^2}{\sigma_\mathrm{obs,i}^2+\sigma_\mathrm{mod,i}^2} .\f]
 */
float chi2_fq(void)
{
	int i;
	float mod, obs, sigmamod, sigmaobs;
	float chi2 = 0.0;

	//If this is not a master task return 0
	if (ThisTask != MasterTask) return 0.0;

	//Go through all data points of the FQs
	for (i = 0; i < All.Nfq; i++)
	{
		//If the data point is within the interpolation range add it to the chi^2
		if (Fq[i].mod_y >= 0.0 && Fq[i].nfit > 0)
		{
			obs      = Fq[i].obs_y;
			mod      = Fq[i].mod_y;
			sigmaobs = Fq[i].obs_sigma;
			sigmamod = Fq[i].mod_sigma;
			chi2    += (obs-mod) * (obs-mod) / (sigmaobs*sigmaobs + sigmamod*sigmamod) * (float)(Fq[i].nfit);
			DoF     += Fq[i].nfit;
			DoFfq   += Fq[i].nfit;
		}
	}

	return chi2;
}


/*! \brief This function computes the contribution to chi^2 from the cosmic SFR density
 *
 *  This function loops over all data points of the CSFRD in the interpolation range and sums
 *  their contribution. The total error used is only the observational error. Each data point
 *  is weighted by its nfit value. The total and partial number of degrees of freedom is
 *  incremented as well.
 *  The returned \f$\chi^2\f$ is given by \f[\chi^2 = \sum_i n_\mathrm{fit,i}
 *  \frac{(\dot\rho^*_\mathrm{obs,i}-\dot\rho^*_\mathrm{mod,i})^2}{\sigma_\mathrm{obs,i}^2+\sigma_\mathrm{mod,i}^2} .\f]
 */
float chi2_csfrd(void)
{
	int i;
	float mod, obs, sigmaobs;
	float chi2 = 0.0;

	//If this is not a master task return 0
	if (ThisTask != MasterTask) return 0.0;

	//Go through all data points of the CSFRD
	for (i = 0; i < All.Ncsfrd; i++)
	{
		//If the data point is within the interpolation range add it to the chi^2
		if (Csfrd[i].mod_y > 0 && Csfrd[i].nfit > 0)
		{
			obs       = Csfrd[i].obs_y;
			mod       = log10(Csfrd[i].mod_y);
			sigmaobs  = Csfrd[i].obs_sigma;
			chi2    += (obs-mod) * (obs-mod) / (sigmaobs*sigmaobs) * (float)(Csfrd[i].nfit);
			DoF      += Csfrd[i].nfit;
			DoFcsfrd += Csfrd[i].nfit;
		}
	}

	return chi2;
}


/*! \brief This function computes the contribution to chi^2 from the specific SFRs
 *
 *  This function loops over all data points of the SSFRs in the interpolation range and sums
 *  their contribution. The total squared error used is the sum of the squared error of the observations
 *  and the model (Poisson). Each data point is weighted by its nfit value. The total and partial
 *  number of degrees of freedom is incremented as well.
 *  The returned \f$\chi^2\f$ is given by \f[\chi^2 = \sum_i n_\mathrm{fit,i}
 *  \frac{(\Psi_\mathrm{obs,i}-\Psi_\mathrm{mod,i})^2}{\sigma_\mathrm{obs,i}^2+\sigma_\mathrm{mod,i}^2} .\f]
 */
float chi2_ssfr(void)
{
	int i;
	float mod, obs, sigmamod, sigmaobs;
	float chi2 = 0.0;

	//If this is not a master task return 0
	if (ThisTask != MasterTask) return 0.0;

	//Go through all data points of the SSFR
	for (i = 0; i < All.Nssfr; i++)
	{
		//If the data point is within the interpolation range add it to the chi^2
		if (Ssfr[i].mod_y > 0.0 && Ssfr[i].nfit > 0)
		{
			obs      = Ssfr[i].obs_y;
			mod      = log10(Ssfr[i].mod_y / All.t_unit);
			sigmaobs = Ssfr[i].obs_sigma;
			sigmamod = log10(Ssfr[i].mod_sigma);
			chi2    += (obs-mod) * (obs-mod) / (sigmaobs*sigmaobs + sigmamod*sigmamod) * (float)(Ssfr[i].nfit);
			DoF     += Ssfr[i].nfit;
			DoFssfr += Ssfr[i].nfit;
		}
	}

	return chi2;
}


/*! \brief This function computes the contribution to chi^2 from the projected correlation functions
 *
 *  This function loops over all data points of the projected correltion functions up to the maximum distance
 *  given by #WP_RMAX times the box length and sums their contribution. The total squared error used is the
 *  sum of the squared error of the observations and the model (Poisson). Each data point is weighted by its
 *  nfit value. The total and partial number of degrees of freedom is incremented as well.
 *  The returned \f$\chi^2\f$ is given by \f[\chi^2 = \sum_i n_\mathrm{fit,i}
 *  \frac{(w_\mathrm{p,obs,i}-w_\mathrm{p,mod,i})^2}{\sigma_\mathrm{obs,i}^2+\sigma_\mathrm{mod,i}^2} .\f]
 */
float chi2_wp(void)
{
	int i;
	float mod, obs, sigmamod, sigmaobs;
	float chi2 = 0.0;

	//If this is not a master task return 0
	if (ThisTask != MasterTask) return 0.0;

	//Go through all data points of the WP
	for (i = 0; i < All.Nwp; i++)
	{
		//If the data point is within the interpolation range add it to the chi^2
		if (Wp[i].mod_y > 0.0 && Wp[i].nfit > 0 && Wp[i].obs_x < WP_RMAX * All.Lbox)
		{
			obs      = Wp[i].obs_y;
			mod      = Wp[i].mod_y;
			sigmaobs = Wp[i].obs_sigma;
			sigmamod = Wp[i].mod_sigma;
			chi2    += (obs-mod) * (obs-mod) / (sigmaobs*sigmaobs + sigmamod*sigmamod) * (float)(Wp[i].nfit);
			DoF     += Wp[i].nfit;
			DoFwp   += Wp[i].nfit;
		}
	}

	return chi2;
}
