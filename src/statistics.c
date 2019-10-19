///////////////////////////////////////////////////////////////////////////////////////////////////
// Emerge code - File statistics.c                                                               //
///////////////////////////////////////////////////////////////////////////////////////////////////
///
/// \file statistics.c
/// \brief Contains functions that compute galaxy statistics
///
/// This file contains all functions that compute galaxy statistics. This includes the function
/// that adds a galaxy to the local statistics on a given task, and a function that combines
/// the local statistics on each task to a global statistics. It also includes a function that
/// writes all statistics of each universe to a file.
///
///////////////////////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <sys/stat.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include "allvars.h"
#include "proto.h"

#ifdef HDF5_SUPPORT
#include <hdf5.h>
#endif


/*! \brief This function adds a galaxy to all relevant statistics
 *
 *  First the observed stellar mass and star formation rate is calculated. For this a gaussian
 *  scatter scaled for the current redshift is used. The scatter is either taken from the random
 *  number table, or if not available directly computed. Then the galaxy is added to the SMF, the
 *  SMF for quenched galaxies, the SSFRs, and the CSFRD.
 *
 *  \param ihalo The index of the halo that is added to the statistics
 *  \param thisthread Number of the OpenMP thread that is responsible
 */
void add_galaxy_to_statistics(int ihalo, int thisthread)
{
	int imass, iscale;
	float sigmaobs, mstarobs, sfrobs, scatter, logmass;

	// Add observational error to stellar masses
	scatter = get_gaussian_random_number(ihalo);
	sigmaobs = All.obssigma0+(1./H[ihalo].a-1.)*All.obssigmaz;
	if (1./H[ihalo].a-1. > ZMAX_SMFERR) sigmaobs = All.obssigma0+(ZMAX_SMFERR)*All.obssigmaz;
	mstarobs = H[ihalo].mstar * pow(10.,sigmaobs*scatter);

	// Add observational error to SFRs
	scatter = get_gaussian_random_number(ihalo + RANDOM_NUMBER_TABLE/2);
	sfrobs = H[ihalo].sfr * pow(10.,sigmaobs*scatter);

	if (mstarobs > 0) logmass = log10(mstarobs*All.m_unit);
	else logmass = All.mstarmin - 1.0;
	imass  = (int)((logmass-All.mstarmin)/All.dmstar);
	iscale = H[ihalo].iscale;

	if (logmass >= All.mstarmin && imass < All.Nmstar)
	{
		//Add galaxy to SMF
		Modelsmf[iscale*All.Nmstar+imass+thisthread*All.NStatistics]
		  += 1.0/All.Lbox/All.Lbox/All.Lbox/All.dmstar;
		//Add galaxy to quenched SMF
		if (sfrobs <= 0.0) Modelsmfred[iscale*All.Nmstar+imass+thisthread*All.NStatistics]
			  += 1.0/All.Lbox/All.Lbox/All.Lbox/All.dmstar;
		else if (log10(sfrobs)-logmass < log10(SSFRTHRESH/CosmicTime[iscale]/All.t_unit))
			Modelsmfred[iscale*All.Nmstar+imass+thisthread*All.NStatistics] += 1.0/All.Lbox/All.Lbox/All.Lbox/All.dmstar;
		//Add galaxy to SSFR
		Modelssfr[iscale*All.Nmstar+imass+thisthread*All.NStatistics] += sfrobs*All.m_unit/All.t_unit/mstarobs;
	}
	//Add galaxy to CSFRD
	Modelcsfrd[iscale+thisthread*All.NTimesteps] += sfrobs * All.m_unit / All.t_unit / All.Lbox / All.Lbox / All.Lbox;

}


/*! \brief This function computes the final statistics from all processors on each Master Task
 *
 *  The stellar mass function in each slave task is copied to a local array smfComm, and then added
 *  to the master task's stellar mass function. The model stellar mass function for each observed
 *  data point is then interpolated from a 2D grid (mstar/a) using bicubic splines. The quenched
 *  fraction of galaxies is set as SMF/SMFred and interpolated as well.
 */
void get_statistics(void)
{
	int i,j,task;
	float fquenched;
	float  *smfComm, *csfrdComm;
	double *xa, *ya, *yb, *za;
	const gsl_interp2d_type *T = gsl_interp2d_bicubic;
	size_t nx,ny;
	gsl_spline *spline1D;
	gsl_spline2d *spline;
	gsl_interp_accel *xacc, *yacc;
	MPI_Status status;

	//Allocate the buffer arrays
	smfComm   = emalloc("SMFComm",   All.NTimesteps * All.Nmstar * sizeof(float));
	csfrdComm = emalloc("CSFRDComm", All.NTimesteps * sizeof(float));

	//If OpenMP is set we need to first collect all statistics using the sub-array for the first thread
#if defined(_OPENMP) && (OPENMPTHREADS > 1)
	for (i = 1; i < NThread; i++)
	{
		for (j = 0; j < All.NStatistics; j++)
		{
			Modelsmf[j]    += Modelsmf[j+i*All.NStatistics];
			Modelsmfred[j] += Modelsmfred[j+i*All.NStatistics];
			Modelssfr[j]   += Modelssfr[j+i*All.NStatistics];
		}
		for (j = 0; j < All.NTimesteps; j++) Modelcsfrd[j] += Modelcsfrd[j+i*All.NTimesteps];
	}
#endif

	//Go through each slave task
	for (task = MasterTask + 1; task < MasterTask + All.NTaskPerUniverse; task++)
	{
		//if this is a slave task...
		if (ThisTask == task)
		{
			//...send the model SMF to the master task
			MPI_Ssend(Modelsmf, All.NTimesteps * All.Nmstar, MPI_FLOAT, MasterTask, TAG_SMF, MPI_COMM_WORLD);
			MPI_Ssend(Modelsmfred, All.NTimesteps * All.Nmstar, MPI_FLOAT, MasterTask, TAG_SMFRED, MPI_COMM_WORLD);
			MPI_Ssend(Modelcsfrd, All.NTimesteps, MPI_FLOAT, MasterTask, TAG_CSFRD, MPI_COMM_WORLD);
			MPI_Ssend(Modelssfr, All.NTimesteps * All.Nmstar, MPI_FLOAT, MasterTask, TAG_SSFR, MPI_COMM_WORLD);
		}
		//if this is the master task...
		if (ThisTask == MasterTask)
		{
			//...receive the SMF from the slave task, write it to the buffer and add it to the local SMF
			MPI_Recv(smfComm, All.NTimesteps * All.Nmstar, MPI_FLOAT, task, TAG_SMF, MPI_COMM_WORLD, &status);
			for (i = 0; i < All.NTimesteps * All.Nmstar; i++) Modelsmf[i] += smfComm[i];
			MPI_Recv(smfComm, All.NTimesteps * All.Nmstar, MPI_FLOAT, task, TAG_SMFRED, MPI_COMM_WORLD, &status);
			for (i = 0; i < All.NTimesteps * All.Nmstar; i++) Modelsmfred[i] += smfComm[i];
			MPI_Recv(csfrdComm, All.NTimesteps, MPI_FLOAT, task, TAG_CSFRD, MPI_COMM_WORLD, &status);
			for (i = 0; i < All.NTimesteps; i++) Modelcsfrd[i] += csfrdComm[i];
			MPI_Recv(smfComm, All.NTimesteps * All.Nmstar, MPI_FLOAT, task, TAG_SSFR, MPI_COMM_WORLD, &status);
			for (i = 0; i < All.NTimesteps * All.Nmstar; i++) Modelssfr[i] += smfComm[i];
		}
	}
	//Free the buffers
	efree(csfrdComm);
	efree(smfComm);

	//The master task now does the interpolation
	if (ThisTask == MasterTask)
	{
		//Define the number of 2D grid points for GSL
		nx = (size_t)(All.NTimesteps);
		ny = (size_t)(All.Nmstar);
		//Allocate double arrays that will be passed to the splines
		xa = emalloc("XA", nx * sizeof(double));
		ya = emalloc("YA", ny * sizeof(double));
		yb = emalloc("YB", nx * sizeof(double));
		za = emalloc("ZA", nx * ny * sizeof(double));
		//Allocate the spline and the acceleration arrays for GSL
		spline = gsl_spline2d_alloc(T, nx, ny);
		spline1D = gsl_spline_alloc(gsl_interp_cspline,nx);
		xacc = gsl_interp_accel_alloc();
		yacc = gsl_interp_accel_alloc();
		//Write the scale factors and stellar masses to the arrays
		for (i = 0; i < nx; i++) xa[i] = (double)(ScaleFactor[i]);
		for (j = 0; j < ny; j++) ya[j] = (double)(Mstar[j]);

		//Now write the SMF to the za array
		for (i = 0; i < nx; i++)
		{
			//For each grid point set the spline
			for (j = 0; j < ny; j++) gsl_spline2d_set(spline, za, i, j, (double)(Modelsmf[i*ny+j]));
		}
		//Initialise the spline
		gsl_spline2d_init(spline, xa, ya, za, nx, ny);
#ifdef READ_SMF
		//Go through all observed data points
		for (i = 0; i < All.Nsmf; i++)
		{
			//If the observed stellar mass is in range do a bicubic spline interpolation
			if (Smf[i].obs_x > Mstar[0] && Smf[i].obs_x < Mstar[All.Nmstar-1] && Smf[i].bin > ScaleFactor[0] && Smf[i].bin <= ScaleFactor[All.NTimesteps-1])
				Smf[i].mod_y = gsl_spline2d_eval(spline, Smf[i].bin, Smf[i].obs_x, xacc, yacc);
			//Otherwise set the model value to -1
			else Smf[i].mod_y = -1.0;
			//And the poisson error
			if (Smf[i].mod_y > 0) Smf[i].mod_sigma = 1.+1./sqrt(Smf[i].mod_y*All.Lbox*All.Lbox*All.Lbox*All.dmstar);
			else Smf[i].mod_sigma = 1.0;
		}
#endif
#ifdef READ_FQ
		//Go through all observed data points
		for (i = 0; i < All.Nfq; i++)
		{
			//If the observed stellar mass is in range do a bicubic spline interpolation
			if (Fq[i].obs_x > Mstar[0] && Fq[i].obs_x < Mstar[All.Nmstar-1] && Fq[i].bin > ScaleFactor[0] && Fq[i].bin <= ScaleFactor[All.NTimesteps-1])
				Fq[i].mod_sigma = gsl_spline2d_eval(spline, Fq[i].bin, Fq[i].obs_x, xacc, yacc);
			//Otherwise set the model value to -1
			else Fq[i].mod_sigma = -1.0;
			if (Fq[i].mod_sigma > 0) Fq[i].mod_sigma = 1./sqrt(Fq[i].mod_sigma*All.Lbox*All.Lbox*All.Lbox*All.dmstar);
			else Fq[i].mod_sigma = 1.0;
		}
#endif
#ifdef READ_SSFR
		//Go through all observed data points
		for (i = 0; i < All.Nssfr; i++)
		{
			if (Ssfr[i].bin > Mstar[0] && Ssfr[i].bin < Mstar[All.Nmstar-1] && 1./(Ssfr[i].obs_x+1.0) > ScaleFactor[0] && 1./(Ssfr[i].obs_x+1.0) <= ScaleFactor[All.NTimesteps-1])
				Ssfr[i].mod_sigma = gsl_spline2d_eval(spline, 1./(Ssfr[i].obs_x+1.0), Ssfr[i].bin, xacc, yacc);
			//Otherwise set the model value to -1
			else Ssfr[i].mod_sigma = -1.0;
			if (Ssfr[i].mod_sigma > 0) Ssfr[i].mod_sigma = 1.+1./sqrt(Ssfr[i].mod_sigma*All.Lbox*All.Lbox*All.Lbox*All.dmstar);
			else Ssfr[i].mod_sigma = 1.0;
		}
#endif
#ifdef READ_FQ
		//Now write the SMF for quenched galaxies to the za array
		for (i = 0; i < nx; i++)
		{
			//For each grid point set the spline
			for (j = 0; j < ny; j++) {
				//Define the Quenched Fraction
				if (Modelsmf[i*ny+j] > 0) fquenched = (double)(Modelsmfred[i*ny+j]) / (double)(Modelsmf[i*ny+j]);
				else fquenched = 0.0;
				gsl_spline2d_set(spline, za, i, j, fquenched);
			}
		}
		//Initialise the spline
		gsl_spline2d_init(spline, xa, ya, za, nx, ny);
		//Go through all observed data points
		for (i = 0; i < All.Nfq; i++)
		{
			//If the observed stellar mass is in range do a bicubic spline interpolation
			if (Fq[i].obs_x > Mstar[0] && Fq[i].obs_x < Mstar[All.Nmstar-1] && Fq[i].bin > ScaleFactor[0] && Fq[i].bin <= ScaleFactor[All.NTimesteps-1])
				Fq[i].mod_y = gsl_spline2d_eval(spline, Fq[i].bin, Fq[i].obs_x, xacc, yacc);
			//Otherwise set the model value to -1
			else Fq[i].mod_y = -1.0;
			Fq[i].mod_sigma *= Fq[i].mod_y;
		}
#endif
#ifdef READ_CSFRD
		//Initialise the spline
		for (i = 0; i < nx; i++) yb[i] = (double)(Modelcsfrd[i]);
		gsl_spline_init(spline1D, xa, yb, nx);
		//Go through all observed data points and do a bicubic spline interpolation
		for (i = 0; i < All.Ncsfrd; i++)
		{
			if (1./(Csfrd[i].obs_x+1.0) > ScaleFactor[0] && 1./(Csfrd[i].obs_x+1.0) <= ScaleFactor[All.NTimesteps-1])
				Csfrd[i].mod_y = gsl_spline_eval(spline1D, 1./(Csfrd[i].obs_x+1.0), xacc);
			else Csfrd[i].mod_y = -1.0;
		}
#endif
		//Normalise the SSFRs
		for (i = 0; i < All.NTimesteps * All.Nmstar; i++)
		{
			if (Modelsmf[i] > 0) Modelssfr[i] /= (Modelsmf[i] * All.Lbox*All.Lbox*All.Lbox*All.dmstar);
			else Modelssfr[i] = 0.0;
		}
#ifdef READ_SSFR
		//Now write the SMF to the za array
		for (i = 0; i < nx; i++)
		{
			//For each grid point set the spline
			for (j = 0; j < ny; j++) gsl_spline2d_set(spline, za, i, j, (double)(Modelssfr[i*ny+j]));
		}
		//Initialise the spline
		gsl_spline2d_init(spline, xa, ya, za, nx, ny);
		//Go through all observed data points
		for (i = 0; i < All.Nssfr; i++)
		{
			//If the observed stellar mass is in range do a bicubic spline interpolation
			if (Ssfr[i].bin > Mstar[0] && Ssfr[i].bin < Mstar[All.Nmstar-1] &&  1./(Ssfr[i].obs_x+1.0) > ScaleFactor[0] &&  1./(Ssfr[i].obs_x+1.0) <= ScaleFactor[All.NTimesteps-1])
				Ssfr[i].mod_y = gsl_spline2d_eval(spline, 1./(Ssfr[i].obs_x+1.0), Ssfr[i].bin, xacc, yacc);
			//Otherwise set the model value to -1.0
			else Ssfr[i].mod_y = -1.0;
		}
#endif
		//Free all used arrays
		efree(za);
		efree(yb);
		efree(ya);
		efree(xa);
		gsl_interp_accel_free(yacc);
		gsl_interp_accel_free(xacc);
		gsl_spline_free(spline1D);
		gsl_spline2d_free(spline);
	}

#ifdef READ_WP
	compute_wp();
#endif

}


/*! \brief This function prints all statistics to either ascii or hdf5 output files
 *
 */
void write_statistics(void)
{
	if (All.OutputFormat == 2)
#ifdef HDF5_SUPPORT
		write_statistics_init_hdf5();
#else
		write_statistics_init_ascii();
#endif
	else if (All.OutputFormat == 1) write_statistics_init_ascii();
	else write_statistics_init_ascii();

}


/*! \brief This function generates the folder for the statistiscs ascii files if needed and set the directory path
 *
 */
void write_statistics_init_ascii(void)
{
	char outdir[NSTRING], outfname[NSTRING];
	FILE *ofp;

	//Get path
	sprintf(outdir,"%s/statistics",All.OutputDir);
	//Create directory if needed
	if (ThisTask == 0) {
		printf("%s\n%s Writing statistics to folder %s\n",All.fullline,All.startline,outdir);
		if (mkdir(outdir, 02755) == 0)
		{
			if (All.verbose >= VERBOSE_MIN) printf("%s\n%s Creating statistics directory '%s'\n",All.fullline,All.startline,outdir);
		}
	}
	//Wait for task 0
	MPI_Barrier(MPI_COMM_WORLD);

	//Task 0 creates the chi2 file
	if (ThisTask ==  0)
	{ //Create Chi2^2 file
		sprintf(outfname,"%s/chi2.out",outdir);
		ofp = fopen(outfname,"w");
		fclose(ofp);
	}

	//Write all chi2 values to the output file
	write_statistics_chi2(outdir,0);
	//Write statistics files
	write_statistics_ascii(outdir,(int)(MasterTask/All.NTaskPerUniverse));

}


#ifdef HDF5_SUPPORT
/*! \brief This function generates the HDF5 file for the statistiscs and set the file name
 *
 */
void write_statistics_init_hdf5(void)
{
	char outfname[NSTRING];
	int itask;
	hid_t stats_file;
	herr_t status;

	//Set the file name
	sprintf(outfname,"%s/statistics.h5",All.OutputDir);

	//Create file
	if (ThisTask == 0)
	{ //Print what's done to the screen
		printf("%s\n%s Writing statistics to file %s\n",All.fullline,All.startline,outfname);
		//Create statistics file by opening and closing
		stats_file = H5Fcreate(outfname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
		status = H5Fclose(stats_file);
		//Warn if there are problems
		if (status != 0) printf("%s Task 0 failed to create file %s\n",All.startline,outfname);
	}

	//Go through all tasks
	for (itask = 0; itask < NTask; itask += All.NTaskPerUniverse)
	{ //If this is a master task write the statistics to a file
		if (itask == ThisTask && ThisTask == MasterTask)
		{ //Each master task writes its own universe to the statistics file
			write_statistics_hdf5(outfname, (int)(MasterTask/All.NTaskPerUniverse));
		}
		//Block all tasks from continuing (so that universes are written sequentially
		MPI_Barrier(MPI_COMM_WORLD);
	}//End Master Task
}
#endif


/*! \brief This function prints all statistics to ascii files
 *
 */
void write_statistics_ascii(char *outdir, int iuniverse)
{
	char outfname[NSTRING];
	int i,j;
	FILE *ofp;

	//The master task now does the interpolation and output
	if (ThisTask == MasterTask)
	{
		//Write the model SMF to a file for each master task
		sprintf(outfname,"%s/smfmod.%d.out",outdir,iuniverse);
		ofp = fopen(outfname,"w");
		//The first column is the stellar mass followed by a SMF value for that mass at each scale factor
		for (i = 0; i < All.Nmstar; i++)
		{
			fprintf(ofp,"%f",Mstar[i]);
			for (j = 0; j < All.NTimesteps; j++)
				fprintf(ofp," %e",Modelsmf[j*All.Nmstar+i]);
			fprintf(ofp,"\n");
		}
		fclose(ofp);
		//Write the model quenched fractions to a file for each master task
		sprintf(outfname,"%s/fqmod.%d.out",outdir,iuniverse);
		ofp = fopen(outfname,"w");
		//The first column is the stellar mass followed by a FQ value for that mass at each scale factor
		for (i = 0; i < All.Nmstar; i++)
		{
			fprintf(ofp,"%f",Mstar[i]);
			for (j = 0; j < All.NTimesteps; j++)
				fprintf(ofp," %e",Modelsmfred[j*All.Nmstar+i]/Modelsmf[j*All.Nmstar+i]);
			fprintf(ofp,"\n");
		}
		fclose(ofp);
		//Write the model cosmic SFR density to a file for each master task
		sprintf(outfname,"%s/csfrdmod.%d.out",outdir,iuniverse);
		ofp = fopen(outfname,"w");
		for (i = 0; i < All.NTimesteps; i++) fprintf(ofp,"%f %e\n",1.0/ScaleFactor[i]-1.0,Modelcsfrd[i]);
		fclose(ofp);
		//Write the model specific SFRs to a file for each master task
		sprintf(outfname,"%s/ssfrmod.%d.out",outdir,iuniverse);
		ofp = fopen(outfname,"w");
		//The first column is the stellar mass followed by a SSFR value for that mass at each scale factor
		for (i = 0; i < All.Nmstar; i++)
		{
			fprintf(ofp,"%f",Mstar[i]);
			for (j = 0; j < All.NTimesteps; j++)
				fprintf(ofp," %e",Modelssfr[j*All.Nmstar+i]);
			fprintf(ofp,"\n");
		}
		fclose(ofp);
#ifdef READ_SMF
		//Now write the observed data points and the computed model SMF values to a file
		sprintf(outfname,"%s/smfobs.%d.out",outdir,iuniverse);
		ofp = fopen(outfname,"w");
		for (i = 0; i < All.Nsmfset; i++)
		{
			fprintf(ofp,"# %d %f %f %s\n",SmfSet[i].ndata,SmfSet[i].min,SmfSet[i].max,SmfSet[i].tag);
			for (j = SmfSet[i].offset; j < SmfSet[i].offset + SmfSet[i].ndata; j++)
				fprintf(ofp,"%f %f %f %f\n",Smf[j].obs_x,Smf[j].obs_y,Smf[j].obs_sigma,log10(Smf[j].mod_y));
			fprintf(ofp,"\n");
		}
		fclose(ofp);
#endif
#ifdef READ_FQ
		//Now write the observed data points and the computed model SMF values to a file
		sprintf(outfname,"%s/fqobs.%d.out",outdir,iuniverse);
		ofp = fopen(outfname,"w");
		for (i = 0; i < All.Nfqset; i++)
		{
			fprintf(ofp,"# %d %f %f %s\n",FqSet[i].ndata,FqSet[i].min,FqSet[i].max,FqSet[i].tag);
			for (j = FqSet[i].offset; j < FqSet[i].offset + FqSet[i].ndata; j++)
				fprintf(ofp,"%f %f %f %f\n",Fq[j].obs_x,Fq[j].obs_y,Fq[j].obs_sigma,Fq[j].mod_y);
			fprintf(ofp,"\n");
		}
		fclose(ofp);
#endif
#ifdef READ_CSFRD
		//Now write the observed data points and the computed model CSFRD values to a file
		sprintf(outfname,"%s/csfrdobs.%d.out",outdir,iuniverse);
		ofp = fopen(outfname,"w");
		for (i = 0; i < All.Ncsfrdset; i++)
		{
			fprintf(ofp,"# %d %s\n",CsfrdSet[i].ndata,CsfrdSet[i].tag);
			for (j = CsfrdSet[i].offset; j < CsfrdSet[i].offset + CsfrdSet[i].ndata; j++)
				fprintf(ofp,"%f %f %f %f\n",Csfrd[j].obs_x,Csfrd[j].obs_y,Csfrd[j].obs_sigma,log10(Csfrd[j].mod_y*All.m_unit/All.t_unit/pow(All.x_unit,3.0)));
			fprintf(ofp,"\n");
		}
		fclose(ofp);
#endif
#ifdef READ_SSFR
		//Now write the observed data points and the computed model SSFR values to a file
		sprintf(outfname,"%s/ssfrobs.%d.out",outdir,iuniverse);
		ofp = fopen(outfname,"w");
		for (i = 0; i < All.Nssfrset; i++)
		{
			fprintf(ofp,"# %d %s\n",SsfrSet[i].ndata,SsfrSet[i].tag);
			for (j = SsfrSet[i].offset; j < SsfrSet[i].offset + SsfrSet[i].ndata; j++)
				fprintf(ofp,"%f %f %f %f %f\n",Ssfr[j].obs_x,Ssfr[j].bin,Ssfr[j].obs_y,Ssfr[j].obs_sigma,log10(Ssfr[j].mod_y/All.t_unit));
			fprintf(ofp,"\n");
		}
		fclose(ofp);
#endif
#ifdef READ_WP
		//Write the model 3d correlation functions to a file for each master task
		sprintf(outfname,"%s/ximod.%d.out",outdir,iuniverse);
		ofp = fopen(outfname,"w");
		//The three columns are the radius, xi and xi's error - repeated for each stellar mass bin
		fprintf(ofp,"#");
		for (j = 0; j < All.Nwpset; j++) fprintf(ofp," [%f %f]",WpSet[j].min,WpSet[j].max);
		fprintf(ofp,"\n");
		for (i = 0; i < WP_RBINS_INT; i++)
		{
			for (j = 0; j < All.Nwpset; j++) fprintf(ofp,"%e %e %e ",Radius[WP_RBINS_INT * j + i], Modelxi[WP_RBINS_INT * j + i], Modelxierr[WP_RBINS_INT * j + i]);
			fprintf(ofp,"\n");
		}
		fclose(ofp);
		//Now write the observed data points and the computed model wp values to a file
		sprintf(outfname,"%s/wpobs.%d.out",outdir,iuniverse);
		ofp = fopen(outfname,"w");
		for (i = 0; i < All.Nwpset; i++)
		{
			fprintf(ofp,"# %d %f %f %f %s\n",WpSet[i].ndata,WpSet[i].min,WpSet[i].max,WpSet[i].cut,WpSet[i].tag);
			for (j = WpSet[i].offset; j < WpSet[i].ndata + WpSet[i].offset; j++)
				fprintf(ofp,"%e %e %e %e %e\n",Wp[j].obs_x,Wp[j].obs_y,Wp[j].obs_sigma,Wp[j].mod_y,Wp[j].mod_sigma);
			fprintf(ofp,"\n");
		}
		fclose(ofp);
#endif
	}

}


#ifdef HDF5_SUPPORT
/*! \brief This function prints all statistics of one universe to a HDF5 file
 *
 */
void write_statistics_hdf5(char *outfname, int iuniverse)
{
	int i, j;
	float   *data, *chi2;;
	char path[NSTRING], buffer[SSTRING];
	hid_t stats_file, universe_group, set_group, data_group, dspace, dset, dtype, ftype, strtype;
	hsize_t dims2[2], dims1[1];
	herr_t status;

	struct comp_smfset {int i1, i2; float f1, f2; char tag[SSTRING];} *compsmfset;
	struct comp_fqset {int i1, i2; float f1, f2; char tag[SSTRING];} *compfqset;
	struct comp_csfrdset {int i1, i2; char tag[SSTRING];} *compcsfrdset;
	struct comp_ssfrset {int i1, i2; char tag[SSTRING];} *compssfrset;
	struct comp_wpset {int i1, i2; float f1, f2, f3; char tag[SSTRING];} *compwpset;
	struct comp_smf {float f1, f2, f3, f4, f5, f6;} *compsmf;
	struct comp_fq {float f1, f2, f3, f4, f5, f6;} *compfq;
	struct comp_csfrd {float f1, f2, f3, f4;} *compcsfrd;
	struct comp_ssfr {float f1, f2, f3, f4, f5, f6;} *compssfr;
	struct comp_xi {float f1, f2, f3;} *compxi;
	struct comp_wp {float f1, f2, f3, f4, f5;} *compwp;

	//If this is the first master task open a new file
	stats_file = H5Fopen(outfname, H5F_ACC_RDWR, H5P_DEFAULT);

	//Set the path name for this universe
	sprintf(path,"/Universe_%06d",iuniverse);
	//Create new group for this universe
	universe_group = H5Gcreate(stats_file, path, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	//////////////////////////////////////////////////////////////////////////////////////////////////////
	//    Stellar mass functions                                                                        //
	//////////////////////////////////////////////////////////////////////////////////////////////////////
	//Set the path name for the SMFs
	sprintf(path,"/Universe_%06d/SMF",iuniverse);
	//Create new group for the SMFs
	set_group = H5Gcreate(stats_file, path, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	//First we write the model SMFs
	//Allocate data array
	data = emalloc("DATA", (All.NTimesteps+1) * (All.Nmstar+1) * sizeof(float));
	//Fill data array
	data[0] = 0.0;
	for (i = 0; i < All.NTimesteps; i++) data[i+1] = ScaleFactor[i];
	for (i = 0; i < All.Nmstar; i++)
	{
		data[(All.NTimesteps+1)*(i+1)] = Mstar[i];
		for (j = 0; j < All.NTimesteps; j++) data[(All.NTimesteps+1)*(i+1)+j+1] = log10(Modelsmf[j*All.Nmstar+i]);
	}
	//Set dimensions array
	dims2[0] = All.Nmstar + 1;
	dims2[1] = All.NTimesteps + 1;
	//Set the path name for the model SMFs
	sprintf(path,"/Universe_%06d/SMF/Model",iuniverse);
	//Create dataspace and dataset, and write the dataset
	dspace = H5Screate_simple(2, dims2, NULL);
	dset   = H5Dcreate(stats_file, path, H5T_NATIVE_FLOAT, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dwrite(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
	// Close the dataset
	status = H5Dclose(dset);
	// Close the data space
	status = H5Sclose(dspace);
	//Free data array
	efree(data);

#ifdef READ_SMF
	//The next block is the information on the data sets
	//Allocate compound
	compsmfset = emalloc("COMPOUND_SMF", All.Nsmfset * sizeof(struct comp_smfset));
	//Fill compound
	for (i = 0; i < All.Nsmfset; i++)
	{
		compsmfset[i].i1 = SmfSet[i].ndata;
		compsmfset[i].i2 = SmfSet[i].offset;
		compsmfset[i].f1 = SmfSet[i].min;
		compsmfset[i].f2 = SmfSet[i].max;
		strcpy(compsmfset[i].tag,SmfSet[i].tag);
	}
	//Create variable-length string datatype.
	strtype= H5Tcopy(H5T_C_S1);
	status = H5Tset_size(strtype, SSTRING);
	//Specify data type
	dtype  = H5Tcreate(H5T_COMPOUND, sizeof(struct comp_smfset));
	status = H5Tinsert(dtype, "Ndata",        HOFFSET(struct comp_smfset, i1), H5T_NATIVE_INT);
	status = H5Tinsert(dtype, "Offset",       HOFFSET(struct comp_smfset, i2), H5T_NATIVE_INT);
	status = H5Tinsert(dtype, "Redshift_min", HOFFSET(struct comp_smfset, f1), H5T_NATIVE_FLOAT);
	status = H5Tinsert(dtype, "Redshift_max", HOFFSET(struct comp_smfset, f2), H5T_NATIVE_FLOAT);
	status = H5Tinsert(dtype, "Tag",          HOFFSET(struct comp_smfset, tag), strtype);
	//Specify file type
	ftype  = H5Tcreate(H5T_COMPOUND, 4 + 4 + 4 + 4 + SSTRING);
	status = H5Tinsert(ftype, "Ndata",        0,       H5T_STD_I32LE);
	status = H5Tinsert(ftype, "Offset",       4,       H5T_STD_I32LE);
	status = H5Tinsert(ftype, "Redshift_min", 4+4,     H5T_IEEE_F32LE);
	status = H5Tinsert(ftype, "Redshift_max", 4+4+4,   H5T_IEEE_F32LE);
	status = H5Tinsert(ftype, "Tag",          4+4+4+4, strtype);
	//Set dimensions array
	dims1[0] = All.Nsmfset;
	//Set the path name for the model SMFs
	sprintf(path,"/Universe_%06d/SMF/Sets",iuniverse);
	//Create dataspace and dataset, and write the compound
	dspace = H5Screate_simple(1, dims1, NULL);
	dset   = H5Dcreate(stats_file, path, ftype, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dwrite(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, compsmfset);
	// Close the dataset
	status = H5Dclose(dset);
	// Close the data space
	status = H5Sclose(dspace);
	// Close the file type
	status = H5Tclose(ftype);
	// Close the data type
	status = H5Tclose(dtype);
	// Close the string type
	status = H5Tclose(strtype);
	//Free compound
	efree(compsmfset);

	//Last entry are all observed SMFs and the corresponding model values
	//Set the path name for the SMFs
	sprintf(path,"/Universe_%06d/SMF/Data",iuniverse);
	//Create new group for the SMF data
	data_group = H5Gcreate(stats_file, path, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	//Specify data type
	dtype  = H5Tcreate(H5T_COMPOUND, sizeof(struct comp_smf));
	status = H5Tinsert(dtype, "Stellar_mass",     HOFFSET(struct comp_smf, f1), H5T_NATIVE_FLOAT);
	status = H5Tinsert(dtype, "Phi_observed",     HOFFSET(struct comp_smf, f2), H5T_NATIVE_FLOAT);
	status = H5Tinsert(dtype, "Sigma_observed",   HOFFSET(struct comp_smf, f3), H5T_NATIVE_FLOAT);
	status = H5Tinsert(dtype, "Mean_ScaleFactor", HOFFSET(struct comp_smf, f4), H5T_NATIVE_FLOAT);
	status = H5Tinsert(dtype, "Phi_model",        HOFFSET(struct comp_smf, f5), H5T_NATIVE_FLOAT);
	status = H5Tinsert(dtype, "Sigma_model",      HOFFSET(struct comp_smf, f6), H5T_NATIVE_FLOAT);
	//Specify file type
	ftype  = H5Tcreate(H5T_COMPOUND, 4 + 4 + 4 + 4 + 4 + 4);
	status = H5Tinsert(ftype, "Stellar_mass",     0,         H5T_IEEE_F32LE);
	status = H5Tinsert(ftype, "Phi_observed",     4,         H5T_IEEE_F32LE);
	status = H5Tinsert(ftype, "Sigma_observed",   4+4,       H5T_IEEE_F32LE);
	status = H5Tinsert(ftype, "Mean_ScaleFactor", 4+4+4,     H5T_IEEE_F32LE);
	status = H5Tinsert(ftype, "Phi_model",        4+4+4+4,   H5T_IEEE_F32LE);
	status = H5Tinsert(ftype, "Sigma_model",      4+4+4+4+4, H5T_IEEE_F32LE);
	//Go through all data sets
	for (i = 0; i < All.Nsmfset; i++)
	{ //Allocate compound
		compsmf = emalloc("COMPOUND_SMF", SmfSet[i].ndata * sizeof(struct comp_smf));
		//Fill compound
		for (j = 0; j < SmfSet[i].ndata; j++)
		{
			compsmf[j].f1 = Smf[SmfSet[i].offset + j].obs_x;
			compsmf[j].f2 = Smf[SmfSet[i].offset + j].obs_y;
			compsmf[j].f3 = Smf[SmfSet[i].offset + j].obs_sigma;
			compsmf[j].f4 = Smf[SmfSet[i].offset + j].bin;
			compsmf[j].f5 = log10(Smf[SmfSet[i].offset + j].mod_y);
			compsmf[j].f6 = log10(Smf[SmfSet[i].offset + j].mod_sigma);
		}
		//Set dimensions array
		dims1[0] = SmfSet[i].ndata;
		//Write tag to buffer and remove slashes so that the path is defined correctly
		strcpy(buffer,SmfSet[i].tag);
		for (j = 0; j < SSTRING; j++) if (buffer[j]=='/') buffer[j] = ' ';
		//Set the path name for the model SMFs
		sprintf(path,"/Universe_%06d/SMF/Data/%03d %s (z = %.2lf - %.2lf)",iuniverse,i,buffer,SmfSet[i].min,SmfSet[i].max);
		//Create dataspace and dataset, and write the compound
		dspace = H5Screate_simple(1, dims1, NULL);
		dset   = H5Dcreate(stats_file, path, ftype, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		status = H5Dwrite(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, compsmf);
		// Close the dataset
		status = H5Dclose(dset);
		// Close the data space
		status = H5Sclose(dspace);
		//Free compound
		efree(compsmf);
	}
	// Close the file type
	status = H5Tclose(ftype);
	// Close the data type
	status = H5Tclose(dtype);
	//Close data group
	status = H5Gclose(data_group);
#endif
	//Close SMFs group
	status = H5Gclose(set_group);

	//////////////////////////////////////////////////////////////////////////////////////////////////////
	//    Quenched Fractions                                                                            //
	//////////////////////////////////////////////////////////////////////////////////////////////////////
	//Set the path name for the FQs
	sprintf(path,"/Universe_%06d/FQ",iuniverse);
	//Create new group for the FQs
	set_group = H5Gcreate(stats_file, path, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	//First we write the model FQs
	//Allocate data array
	data = emalloc("DATA", (All.NTimesteps+1) * (All.Nmstar+1) * sizeof(float));
	//Fill data array
	data[0] = 0.0;
	for (i = 0; i < All.NTimesteps; i++) data[i+1] = ScaleFactor[i];
	for (i = 0; i < All.Nmstar; i++)
	{
		data[(All.NTimesteps+1)*(i+1)] = Mstar[i];
		for (j = 0; j < All.NTimesteps; j++) data[(All.NTimesteps+1)*(i+1)+j+1] = Modelsmfred[j*All.Nmstar+i]/Modelsmf[j*All.Nmstar+i];
	}
	//Set dimensions array
	dims2[0] = All.Nmstar + 1;
	dims2[1] = All.NTimesteps + 1;
	//Set the path name for the model FQs
	sprintf(path,"/Universe_%06d/FQ/Model",iuniverse);
	//Create dataspace and dataset, and write the dataset
	dspace = H5Screate_simple(2, dims2, NULL);
	dset   = H5Dcreate(stats_file, path, H5T_NATIVE_FLOAT, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dwrite(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
	// Close the dataset
	status = H5Dclose(dset);
	// Close the data space
	status = H5Sclose(dspace);
	//Free data array
	efree(data);

#ifdef READ_FQ
	//The next block is the information on the data sets
	//Allocate compound
	compfqset = emalloc("COMPOUND_FQ", All.Nfqset * sizeof(struct comp_fqset));
	//Fill compound
	for (i = 0; i < All.Nfqset; i++)
	{
		compfqset[i].i1 = FqSet[i].ndata;
		compfqset[i].i2 = FqSet[i].offset;
		compfqset[i].f1 = FqSet[i].min;
		compfqset[i].f2 = FqSet[i].max;
		strcpy(compfqset[i].tag,FqSet[i].tag);
	}
	//Create variable-length string datatype.
	strtype= H5Tcopy(H5T_C_S1);
	status = H5Tset_size(strtype, SSTRING);
	//Specify data type
	dtype  = H5Tcreate(H5T_COMPOUND, sizeof(struct comp_fqset));
	status = H5Tinsert(dtype, "Ndata",        HOFFSET(struct comp_fqset, i1), H5T_NATIVE_INT);
	status = H5Tinsert(dtype, "Offset",       HOFFSET(struct comp_fqset, i2), H5T_NATIVE_INT);
	status = H5Tinsert(dtype, "Redshift_min", HOFFSET(struct comp_fqset, f1), H5T_NATIVE_FLOAT);
	status = H5Tinsert(dtype, "Redshift_max", HOFFSET(struct comp_fqset, f2), H5T_NATIVE_FLOAT);
	status = H5Tinsert(dtype, "Tag",          HOFFSET(struct comp_fqset, tag), strtype);
	//Specify file type
	ftype  = H5Tcreate(H5T_COMPOUND, 4 + 4 + 4 + 4 + SSTRING);
	status = H5Tinsert(ftype, "Ndata",        0,       H5T_STD_I32LE);
	status = H5Tinsert(ftype, "Offset",       4,       H5T_STD_I32LE);
	status = H5Tinsert(ftype, "Redshift_min", 4+4,     H5T_IEEE_F32LE);
	status = H5Tinsert(ftype, "Redshift_max", 4+4+4,   H5T_IEEE_F32LE);
	status = H5Tinsert(ftype, "Tag",          4+4+4+4, strtype);
	//Set dimensions array
	dims1[0] = All.Nfqset;
	//Set the path name for the model FQs
	sprintf(path,"/Universe_%06d/FQ/Sets",iuniverse);
	//Create dataspace and dataset, and write the compound
	dspace = H5Screate_simple(1, dims1, NULL);
	dset   = H5Dcreate(stats_file, path, ftype, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dwrite(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, compfqset);
	// Close the dataset
	status = H5Dclose(dset);
	// Close the data space
	status = H5Sclose(dspace);
	// Close the file type
	status = H5Tclose(ftype);
	// Close the data type
	status = H5Tclose(dtype);
	// Close the string type
	status = H5Tclose(strtype);
	//Free compound
	efree(compfqset);

	//Last entry are all observed FQs and the corresponding model values
	//Set the path name for the FQs
	sprintf(path,"/Universe_%06d/FQ/Data",iuniverse);
	//Create new group for the FQ data
	data_group = H5Gcreate(stats_file, path, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	//Specify data type
	dtype  = H5Tcreate(H5T_COMPOUND, sizeof(struct comp_fq));
	status = H5Tinsert(dtype, "Stellar_mass",     HOFFSET(struct comp_fq, f1), H5T_NATIVE_FLOAT);
	status = H5Tinsert(dtype, "Fq_observed",      HOFFSET(struct comp_fq, f2), H5T_NATIVE_FLOAT);
	status = H5Tinsert(dtype, "Sigma_observed",   HOFFSET(struct comp_fq, f3), H5T_NATIVE_FLOAT);
	status = H5Tinsert(dtype, "Mean_ScaleFactor", HOFFSET(struct comp_fq, f4), H5T_NATIVE_FLOAT);
	status = H5Tinsert(dtype, "Fq_model",         HOFFSET(struct comp_fq, f5), H5T_NATIVE_FLOAT);
	status = H5Tinsert(dtype, "Sigma_model",      HOFFSET(struct comp_fq, f6), H5T_NATIVE_FLOAT);
	//Specify file type
	ftype  = H5Tcreate(H5T_COMPOUND, 4 + 4 + 4 + 4 + 4 + 4);
	status = H5Tinsert(ftype, "Stellar_mass",     0,         H5T_IEEE_F32LE);
	status = H5Tinsert(ftype, "Fq_observed",      4,         H5T_IEEE_F32LE);
	status = H5Tinsert(ftype, "Sigma_observed",   4+4,       H5T_IEEE_F32LE);
	status = H5Tinsert(ftype, "Mean_ScaleFactor", 4+4+4,     H5T_IEEE_F32LE);
	status = H5Tinsert(ftype, "Fq_model",         4+4+4+4,   H5T_IEEE_F32LE);
	status = H5Tinsert(ftype, "Sigma_model",      4+4+4+4+4, H5T_IEEE_F32LE);
	//Go through all data sets
	for (i = 0; i < All.Nfqset; i++)
	{ //Allocate compound
		compfq = emalloc("COMPOUND_FQ", FqSet[i].ndata * sizeof(struct comp_fq));
		//Fill compound
		for (j = 0; j < FqSet[i].ndata; j++)
		{
			compfq[j].f1 = Fq[FqSet[i].offset + j].obs_x;
			compfq[j].f2 = Fq[FqSet[i].offset + j].obs_y;
			compfq[j].f3 = Fq[FqSet[i].offset + j].obs_sigma;
			compfq[j].f4 = Fq[FqSet[i].offset + j].bin;
			compfq[j].f5 = Fq[FqSet[i].offset + j].mod_y;
			compfq[j].f6 = Fq[FqSet[i].offset + j].mod_sigma;
		}
		//Set dimensions array
		dims1[0] = FqSet[i].ndata;
		//Write tag to buffer and remove slashes so that the path is defined correctly
		strcpy(buffer,FqSet[i].tag);
		for (j = 0; j < SSTRING; j++) if (buffer[j]=='/') buffer[j] = ' ';
		//Set the path name for the model FQs
		sprintf(path,"/Universe_%06d/FQ/Data/%03d %s (z = %.2lf - %.2lf)",iuniverse,i,buffer,FqSet[i].min,FqSet[i].max);
		//Create dataspace and dataset, and write the compound
		dspace = H5Screate_simple(1, dims1, NULL);
		dset   = H5Dcreate(stats_file, path, ftype, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		status = H5Dwrite(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, compfq);
		// Close the dataset
		status = H5Dclose(dset);
		// Close the data space
		status = H5Sclose(dspace);
		//Free compound
		efree(compfq);
	}
	// Close the file type
	status = H5Tclose(ftype);
	// Close the data type
	status = H5Tclose(dtype);
	//Close data group
	status = H5Gclose(data_group);
#endif
	//Close FQs group
	status = H5Gclose(set_group);

	//////////////////////////////////////////////////////////////////////////////////////////////////////
	//    Specific star formation rate density                                                          //
	//////////////////////////////////////////////////////////////////////////////////////////////////////
	//Set the path name for the CSFRD
	sprintf(path,"/Universe_%06d/CSFRD",iuniverse);
	//Create new group for the CSFRD
	set_group = H5Gcreate(stats_file, path, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	//First we write the model CSFRD
	//Allocate data array
	data = emalloc("DATA", All.NTimesteps * 2 * sizeof(float));
	//Fill data array
	for (i = 0; i < All.NTimesteps; i++)
	{
		data[i*2]   = 1.0 / ScaleFactor[i] - 1.0;
		data[i*2+1] = log10(Modelcsfrd[i]);
	}
	//Set dimensions array
	dims2[0] = All.NTimesteps;
	dims2[1] = 2;
	//Set the path name for the model CSFRD
	sprintf(path,"/Universe_%06d/CSFRD/Model",iuniverse);
	//Create dataspace and dataset, and write the dataset
	dspace = H5Screate_simple(2, dims2, NULL);
	dset   = H5Dcreate(stats_file, path, H5T_NATIVE_FLOAT, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dwrite(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
	// Close the dataset
	status = H5Dclose(dset);
	// Close the data space
	status = H5Sclose(dspace);
	//Free data array
	efree(data);

#ifdef READ_CSFRD
	//The next block is the information on the data sets
	//Allocate compound
	compcsfrdset = emalloc("COMPOUND_CSFRD", All.Ncsfrdset * sizeof(struct comp_csfrdset));
	//Fill compound
	for (i = 0; i < All.Ncsfrdset; i++)
	{
		compcsfrdset[i].i1 = CsfrdSet[i].ndata;
		compcsfrdset[i].i2 = CsfrdSet[i].offset;
		strcpy(compcsfrdset[i].tag,CsfrdSet[i].tag);
	}
	//Create variable-length string datatype.
	strtype= H5Tcopy(H5T_C_S1);
	status = H5Tset_size(strtype, SSTRING);
	//Specify data type
	dtype  = H5Tcreate(H5T_COMPOUND, sizeof(struct comp_csfrdset));
	status = H5Tinsert(dtype, "Ndata",        HOFFSET(struct comp_csfrdset, i1), H5T_NATIVE_INT);
	status = H5Tinsert(dtype, "Offset",       HOFFSET(struct comp_csfrdset, i2), H5T_NATIVE_INT);
	status = H5Tinsert(dtype, "Tag",          HOFFSET(struct comp_csfrdset, tag), strtype);
	//Specify file type
	ftype  = H5Tcreate(H5T_COMPOUND, 4 + 4 + SSTRING);
	status = H5Tinsert(ftype, "Ndata",        0,   H5T_STD_I32LE);
	status = H5Tinsert(ftype, "Offset",       4,   H5T_STD_I32LE);
	status = H5Tinsert(ftype, "Tag",          4+4, strtype);
	//Set dimensions array
	dims1[0] = All.Ncsfrdset;
	//Set the path name for the model FQs
	sprintf(path,"/Universe_%06d/CSFRD/Sets",iuniverse);
	//Create dataspace and dataset, and write the compound
	dspace = H5Screate_simple(1, dims1, NULL);
	dset   = H5Dcreate(stats_file, path, ftype, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dwrite(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, compfqset);
	// Close the dataset
	status = H5Dclose(dset);
	// Close the data space
	status = H5Sclose(dspace);
	// Close the file type
	status = H5Tclose(ftype);
	// Close the data type
	status = H5Tclose(dtype);
	// Close the string type
	status = H5Tclose(strtype);
	//Free compound
	efree(compcsfrdset);

	//Last entry are all observed CSFRDs and the corresponding model values
	//Set the path name for the CSFRD
	sprintf(path,"/Universe_%06d/CSFRD/Data",iuniverse);
	//Create new group for the CSFRD data
	data_group = H5Gcreate(stats_file, path, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	//Specify data type
	dtype  = H5Tcreate(H5T_COMPOUND, sizeof(struct comp_csfrd));
	status = H5Tinsert(dtype, "Redshift",       HOFFSET(struct comp_csfrd, f1), H5T_NATIVE_FLOAT);
	status = H5Tinsert(dtype, "Csfrd_observed", HOFFSET(struct comp_csfrd, f2), H5T_NATIVE_FLOAT);
	status = H5Tinsert(dtype, "Sigma_observed", HOFFSET(struct comp_csfrd, f3), H5T_NATIVE_FLOAT);
	status = H5Tinsert(dtype, "Csfrd_model",    HOFFSET(struct comp_csfrd, f4), H5T_NATIVE_FLOAT);
	//Specify file type
	ftype  = H5Tcreate(H5T_COMPOUND, 4 + 4 + 4 + 4);
	status = H5Tinsert(ftype, "Redshift",       0,       H5T_IEEE_F32LE);
	status = H5Tinsert(ftype, "Csfrd_observed", 4,       H5T_IEEE_F32LE);
	status = H5Tinsert(ftype, "Sigma_observed", 4+4,     H5T_IEEE_F32LE);
	status = H5Tinsert(ftype, "Csfrd_model",    4+4+4,   H5T_IEEE_F32LE);
	//Go through all data sets
	for (i = 0; i < All.Ncsfrdset; i++)
	{ //Allocate compound
		compcsfrd = emalloc("COMPOUND_CSFRD", CsfrdSet[i].ndata * sizeof(struct comp_csfrd));
		//Fill compound
		for (j = 0; j < CsfrdSet[i].ndata; j++)
		{
			compcsfrd[j].f1 = Csfrd[CsfrdSet[i].offset + j].obs_x;
			compcsfrd[j].f2 = Csfrd[CsfrdSet[i].offset + j].obs_y;
			compcsfrd[j].f3 = Csfrd[CsfrdSet[i].offset + j].obs_sigma;
			compcsfrd[j].f4 = log10(Csfrd[CsfrdSet[i].offset + j].mod_y * All.m_unit / All.t_unit * pow(All.x_unit,-3.0));
		}
		//Set dimensions array
		dims1[0] = CsfrdSet[i].ndata;
		//Write tag to buffer and remove slashes so that the path is defined correctly
		strcpy(buffer,CsfrdSet[i].tag);
		for (j = 0; j < SSTRING; j++) if (buffer[j]=='/') buffer[j] = ' ';
		//Set the path name for the model CSFRD
		sprintf(path,"/Universe_%06d/CSFRD/Data/%03d %s",iuniverse,i,buffer);
		//Create dataspace and dataset, and write the compound
		dspace = H5Screate_simple(1, dims1, NULL);
		dset   = H5Dcreate(stats_file, path, ftype, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		status = H5Dwrite(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, compcsfrd);
		// Close the dataset
		status = H5Dclose(dset);
		// Close the data space
		status = H5Sclose(dspace);
		//Free compound
		efree(compcsfrd);
	}
	// Close the file type
	status = H5Tclose(ftype);
	// Close the data type
	status = H5Tclose(dtype);
	//Close data group
	status = H5Gclose(data_group);

#endif
	//Close SFRDS group
	status = H5Gclose(set_group);

	//////////////////////////////////////////////////////////////////////////////////////////////////////
	//    Specific Star Formation Rates                                                                 //
	//////////////////////////////////////////////////////////////////////////////////////////////////////
	//Set the path name for the SSFRs
	sprintf(path,"/Universe_%06d/SSFR",iuniverse);
	//Create new group for the SSFRs
	set_group = H5Gcreate(stats_file, path, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	//First we write the model SSFRs
	//Allocate data array
	data = emalloc("DATA", (All.NTimesteps+1) * (All.Nmstar+1) * sizeof(float));
	//Fill data array
	data[0] = 0.0;
	for (i = 0; i < All.Nmstar; i++) data[i+1] = Mstar[i];
	for (i = 0; i < All.NTimesteps; i++)
	{
		data[(All.Nmstar+1)*(i+1)] = 1.0 / ScaleFactor[i] - 1.0;
		for (j = 0; j < All.Nmstar; j++) data[(All.Nmstar+1)*(i+1)+j+1] = log10(Modelssfr[i*All.Nmstar+j]/All.t_unit);
	}
	//Set dimensions array
	dims2[0] = All.NTimesteps + 1;
	dims2[1] = All.Nmstar + 1;
	//Set the path name for the model SSFRs
	sprintf(path,"/Universe_%06d/SSFR/Model",iuniverse);
	//Create dataspace and dataset, and write the dataset
	dspace = H5Screate_simple(2, dims2, NULL);
	dset   = H5Dcreate(stats_file, path, H5T_NATIVE_FLOAT, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dwrite(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
	// Close the dataset
	status = H5Dclose(dset);
	// Close the data space
	status = H5Sclose(dspace);
	//Free data array
	efree(data);

#ifdef READ_SSFR
	//The next block is the information on the data sets
	//Allocate compound
	compssfrset = emalloc("COMPOUND_SSFR", All.Nssfrset * sizeof(struct comp_ssfrset));
	//Fill compound
	for (i = 0; i < All.Nssfrset; i++)
	{
		compssfrset[i].i1 = SsfrSet[i].ndata;
		compssfrset[i].i2 = SsfrSet[i].offset;
		strcpy(compssfrset[i].tag,SsfrSet[i].tag);
	}
	//Create variable-length string datatype.
	strtype= H5Tcopy(H5T_C_S1);
	status = H5Tset_size(strtype, SSTRING);
	//Specify data type
	dtype  = H5Tcreate(H5T_COMPOUND, sizeof(struct comp_ssfrset));
	status = H5Tinsert(dtype, "Ndata",  HOFFSET(struct comp_ssfrset, i1), H5T_NATIVE_INT);
	status = H5Tinsert(dtype, "Offset", HOFFSET(struct comp_ssfrset, i2), H5T_NATIVE_INT);
	status = H5Tinsert(dtype, "Tag",    HOFFSET(struct comp_ssfrset, tag), strtype);
	//Specify file type
	ftype  = H5Tcreate(H5T_COMPOUND, 4 + 4 + SSTRING);
	status = H5Tinsert(ftype, "Ndata",        0,   H5T_STD_I32LE);
	status = H5Tinsert(ftype, "Offset",       4,   H5T_STD_I32LE);
	status = H5Tinsert(ftype, "Tag",          4+4, strtype);
	//Set dimensions array
	dims1[0] = All.Nssfrset;
	//Set the path name for the model SSFRs
	sprintf(path,"/Universe_%06d/SSFR/Sets",iuniverse);
	//Create dataspace and dataset, and write the compound
	dspace = H5Screate_simple(1, dims1, NULL);
	dset   = H5Dcreate(stats_file, path, ftype, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dwrite(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, compssfrset);
	// Close the dataset
	status = H5Dclose(dset);
	// Close the data space
	status = H5Sclose(dspace);
	// Close the file type
	status = H5Tclose(ftype);
	// Close the data type
	status = H5Tclose(dtype);
	// Close the string type
	status = H5Tclose(strtype);
	//Free compound
	efree(compssfrset);

	//Last entry are all observed SSFRs and the corresponding model values
	//Set the path name for the SSFRs
	sprintf(path,"/Universe_%06d/SSFR/Data",iuniverse);
	//Create new group for the SSFR data
	data_group = H5Gcreate(stats_file, path, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	//Specify data type
	dtype  = H5Tcreate(H5T_COMPOUND, sizeof(struct comp_ssfr));
	status = H5Tinsert(dtype, "Redshift",       HOFFSET(struct comp_ssfr, f1), H5T_NATIVE_FLOAT);
	status = H5Tinsert(dtype, "Ssfr_observed",  HOFFSET(struct comp_ssfr, f2), H5T_NATIVE_FLOAT);
	status = H5Tinsert(dtype, "Sigma_observed", HOFFSET(struct comp_ssfr, f3), H5T_NATIVE_FLOAT);
	status = H5Tinsert(dtype, "Stellar_mass",   HOFFSET(struct comp_ssfr, f4), H5T_NATIVE_FLOAT);
	status = H5Tinsert(dtype, "Ssfr_model",     HOFFSET(struct comp_ssfr, f5), H5T_NATIVE_FLOAT);
	status = H5Tinsert(dtype, "Sigma_model",    HOFFSET(struct comp_ssfr, f6), H5T_NATIVE_FLOAT);
	//Specify file type
	ftype  = H5Tcreate(H5T_COMPOUND, 4 + 4 + 4 + 4 + 4 + 4);
	status = H5Tinsert(ftype, "Redshift",       0,         H5T_IEEE_F32LE);
	status = H5Tinsert(ftype, "Ssfr_observed",  4,         H5T_IEEE_F32LE);
	status = H5Tinsert(ftype, "Sigma_observed", 4+4,       H5T_IEEE_F32LE);
	status = H5Tinsert(ftype, "Stellar_mass",   4+4+4,     H5T_IEEE_F32LE);
	status = H5Tinsert(ftype, "Ssfr_model",     4+4+4+4,   H5T_IEEE_F32LE);
	status = H5Tinsert(ftype, "Sigma_model",    4+4+4+4+4, H5T_IEEE_F32LE);
	//Go through all data sets
	for (i = 0; i < All.Nssfrset; i++)
	{ //Allocate compound
		compssfr = emalloc("COMPOUND_SSFR", SsfrSet[i].ndata * sizeof(struct comp_ssfr));
		//Fill compound
		for (j = 0; j < SsfrSet[i].ndata; j++)
		{
			compssfr[j].f1 = Ssfr[SsfrSet[i].offset + j].obs_x;
			compssfr[j].f2 = Ssfr[SsfrSet[i].offset + j].obs_y;
			compssfr[j].f3 = Ssfr[SsfrSet[i].offset + j].obs_sigma;
			compssfr[j].f4 = Ssfr[SsfrSet[i].offset + j].bin;
			compssfr[j].f5 = log10(Ssfr[SsfrSet[i].offset + j].mod_y / All.t_unit);
			compssfr[j].f6 = log10(Ssfr[SsfrSet[i].offset + j].mod_sigma);
		}
		//Set dimensions array
		dims1[0] = SsfrSet[i].ndata;
		//Write tag to buffer and remove slashes so that the path is defined correctly
		strcpy(buffer,SsfrSet[i].tag);
		for (j = 0; j < SSTRING; j++) if (buffer[j]=='/') buffer[j] = ' ';
		//Set the path name for the model FQs
		sprintf(path,"/Universe_%06d/SSFR/Data/%03d %s",iuniverse,i,buffer);
		//Create dataspace and dataset, and write the compound
		dspace = H5Screate_simple(1, dims1, NULL);
		dset   = H5Dcreate(stats_file, path, ftype, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		status = H5Dwrite(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, compssfr);
		// Close the dataset
		status = H5Dclose(dset);
		// Close the data space
		status = H5Sclose(dspace);
		//Free compound
		efree(compssfr);
	}
	// Close the file type
	status = H5Tclose(ftype);
	// Close the data type
	status = H5Tclose(dtype);
	//Close data group
	status = H5Gclose(data_group);
#endif
	//Close SSFRs group
	status = H5Gclose(set_group);

	//////////////////////////////////////////////////////////////////////////////////////////////////////
	//    Clustering                                                                                    //
	//////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef READ_WP
	//Set the path name for the clustering
	sprintf(path,"/Universe_%06d/Clustering",iuniverse);
	//Create new group for the clustering
	set_group = H5Gcreate(stats_file, path, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	//First entry are all 3D model correlation functions
	//Set the path name for the 3D correlation functions
	sprintf(path,"/Universe_%06d/Clustering/Model_3D",iuniverse);
	//Create new group for the 3D correlation functions
	data_group = H5Gcreate(stats_file, path, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	//Specify data type
	dtype  = H5Tcreate(H5T_COMPOUND, sizeof(struct comp_xi));
	status = H5Tinsert(dtype, "Radius",   HOFFSET(struct comp_xi, f1), H5T_NATIVE_FLOAT);
	status = H5Tinsert(dtype, "Xi",       HOFFSET(struct comp_xi, f2), H5T_NATIVE_FLOAT);
	status = H5Tinsert(dtype, "Sigma_xi", HOFFSET(struct comp_xi, f3), H5T_NATIVE_FLOAT);
	//Specify file type
	ftype  = H5Tcreate(H5T_COMPOUND, 4 + 4 + 4);
	status = H5Tinsert(ftype, "Radius",   0,   H5T_IEEE_F32LE);
	status = H5Tinsert(ftype, "Xi",       4,   H5T_IEEE_F32LE);
	status = H5Tinsert(ftype, "Sigma_xi", 4+4, H5T_IEEE_F32LE);
	//Go through all data sets
	for (i = 0; i < All.Nwpset; i++)
	{ //Allocate compound
		compxi = emalloc("COMPOUND_XI", WP_RBINS_INT * sizeof(struct comp_xi));
		//Fill compound
		for (j = 0; j < WP_RBINS_INT; j++)
		{
			compxi[j].f1 = Radius[WP_RBINS_INT * i + j];
			compxi[j].f2 = Modelxi[WP_RBINS_INT * i + j];
			compxi[j].f3 = Modelxierr[WP_RBINS_INT * i + j];
		}
		//Set dimensions array
		dims1[0] = WP_RBINS_INT;
		//Set the path name for the model correlation function
		sprintf(path,"/Universe_%06d/Clustering/Model_3D/%03d - m = %.2f-%.2f - z = %.2f",iuniverse,i,log10(WpSet[i].min*All.m_unit),log10(WpSet[i].max*All.m_unit),All.wpredshift);
		//Create dataspace and dataset, and write the compound
		dspace = H5Screate_simple(1, dims1, NULL);
		dset   = H5Dcreate(stats_file, path, ftype, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		status = H5Dwrite(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, compxi);
		// Close the dataset
		status = H5Dclose(dset);
		// Close the data space
		status = H5Sclose(dspace);
		//Free compound
		efree(compxi);
	}
	// Close the file type
	status = H5Tclose(ftype);
	// Close the data type
	status = H5Tclose(dtype);
	//Close data group
	status = H5Gclose(data_group);

	//The next block is the information on the data sets
	//Allocate compound
	compwpset = emalloc("COMPOUND_WP", All.Nwpset * sizeof(struct comp_wpset));
	//Fill compound
	for (i = 0; i < All.Nwpset; i++)
	{
		compwpset[i].i1 = WpSet[i].ndata;
		compwpset[i].i2 = WpSet[i].offset;
		compwpset[i].f1 = log10(WpSet[i].min*All.m_unit);
		compwpset[i].f2 = log10(WpSet[i].max*All.m_unit);
		compwpset[i].f3 = WpSet[i].cut;
		strcpy(compwpset[i].tag,WpSet[i].tag);
	}
	//Create variable-length string datatype.
	strtype= H5Tcopy(H5T_C_S1);
	status = H5Tset_size(strtype, SSTRING);
	//Specify data type
	dtype  = H5Tcreate(H5T_COMPOUND, sizeof(struct comp_wpset));
	status = H5Tinsert(dtype, "Ndata",        HOFFSET(struct comp_wpset, i1), H5T_NATIVE_INT);
	status = H5Tinsert(dtype, "Offset",       HOFFSET(struct comp_wpset, i2), H5T_NATIVE_INT);
	status = H5Tinsert(dtype, "Minimum_Mass", HOFFSET(struct comp_wpset, f1), H5T_NATIVE_FLOAT);
	status = H5Tinsert(dtype, "Maximum_Mass", HOFFSET(struct comp_wpset, f2), H5T_NATIVE_FLOAT);
	status = H5Tinsert(dtype, "Pi_max",       HOFFSET(struct comp_wpset, f3), H5T_NATIVE_FLOAT);
	status = H5Tinsert(dtype, "Tag",          HOFFSET(struct comp_wpset, tag), strtype);
	//Specify file type
	ftype  = H5Tcreate(H5T_COMPOUND, 4 + 4 + 4 + 4 + 4 + SSTRING);
	status = H5Tinsert(ftype, "Ndata",        0,         H5T_STD_I32LE);
	status = H5Tinsert(ftype, "Offset",       4,         H5T_STD_I32LE);
	status = H5Tinsert(ftype, "Minimum_Mass", 4+4,       H5T_IEEE_F32LE);
	status = H5Tinsert(ftype, "Maximum_Mass", 4+4+4,     H5T_IEEE_F32LE);
	status = H5Tinsert(ftype, "Pi_max",       4+4+4+4,   H5T_IEEE_F32LE);
	status = H5Tinsert(ftype, "Tag",          4+4+4+4+4, strtype);
	//Set dimensions array
	dims1[0] = All.Nwpset;
	//Set the path name for the model projected correlation functions
	sprintf(path,"/Universe_%06d/Clustering/Sets",iuniverse);
	//Create dataspace and dataset, and write the compound
	dspace = H5Screate_simple(1, dims1, NULL);
	dset   = H5Dcreate(stats_file, path, ftype, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dwrite(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, compwpset);
	// Close the dataset
	status = H5Dclose(dset);
	// Close the data space
	status = H5Sclose(dspace);
	// Close the file type
	status = H5Tclose(ftype);
	// Close the data type
	status = H5Tclose(dtype);
	// Close the string type
	status = H5Tclose(strtype);
	//Free compound
	efree(compwpset);

	//Last entry are all observed projected correlation functions and the corresponding model values
	//Set the path name for the WPs
	sprintf(path,"/Universe_%06d/Clustering/Data",iuniverse);
	//Create new group for the WP data
	data_group = H5Gcreate(stats_file, path, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	//Specify data type
	dtype  = H5Tcreate(H5T_COMPOUND, sizeof(struct comp_wp));
	status = H5Tinsert(dtype, "Radius",         HOFFSET(struct comp_wp, f1), H5T_NATIVE_FLOAT);
	status = H5Tinsert(dtype, "Wp_observed",    HOFFSET(struct comp_wp, f2), H5T_NATIVE_FLOAT);
	status = H5Tinsert(dtype, "Sigma_observed", HOFFSET(struct comp_wp, f3), H5T_NATIVE_FLOAT);
	status = H5Tinsert(dtype, "Wp_model",       HOFFSET(struct comp_wp, f4), H5T_NATIVE_FLOAT);
	status = H5Tinsert(dtype, "Sigma_model",    HOFFSET(struct comp_wp, f5), H5T_NATIVE_FLOAT);
	//Specify file type
	ftype  = H5Tcreate(H5T_COMPOUND, 4 + 4 + 4 + 4 + 4);
	status = H5Tinsert(ftype, "Radius",         0,       H5T_IEEE_F32LE);
	status = H5Tinsert(ftype, "Wp_observed",    4,       H5T_IEEE_F32LE);
	status = H5Tinsert(ftype, "Sigma_observed", 4+4,     H5T_IEEE_F32LE);
	status = H5Tinsert(ftype, "Wp_model",       4+4+4,   H5T_IEEE_F32LE);
	status = H5Tinsert(ftype, "Sigma_model",    4+4+4+4, H5T_IEEE_F32LE);
	//Go through all data sets
	for (i = 0; i < All.Nwpset; i++)
	{ //Allocate compound
		compwp = emalloc("COMPOUND_WP", WpSet[i].ndata * sizeof(struct comp_wp));
		//Fill compound
		for (j = 0; j < WpSet[i].ndata; j++)
		{
			compwp[j].f1 = Wp[WpSet[i].offset + j].obs_x;
			compwp[j].f2 = Wp[WpSet[i].offset + j].obs_y;
			compwp[j].f3 = Wp[WpSet[i].offset + j].obs_sigma;
			compwp[j].f4 = Wp[WpSet[i].offset + j].mod_y;
			compwp[j].f5 = Wp[WpSet[i].offset + j].mod_sigma;
		}
		//Set dimensions array
		dims1[0] = WpSet[i].ndata;
		//Write tag to buffer and remove slashes so that the path is defined correctly
		strcpy(buffer,WpSet[i].tag);
		for (j = 0; j < SSTRING; j++) if (buffer[j]=='/') buffer[j] = ' ';
		//Set the path name for the model WPs
		sprintf(path,"/Universe_%06d/Clustering/Data/%03d %s - m = %.2f-%.2f - z = %.2f",iuniverse,i,buffer,log10(WpSet[i].min*All.m_unit),log10(WpSet[i].max*All.m_unit),All.wpredshift);
		//Create dataspace and dataset, and write the compound
		dspace = H5Screate_simple(1, dims1, NULL);
		dset   = H5Dcreate(stats_file, path, ftype, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		status = H5Dwrite(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, compwp);
		// Close the dataset
		status = H5Dclose(dset);
		// Close the data space
		status = H5Sclose(dspace);
		//Free compound
		efree(compwp);
	}
	// Close the file type
	status = H5Tclose(ftype);
	// Close the data type
	status = H5Tclose(dtype);
	//Close data group
	status = H5Gclose(data_group);

	//Close clustering group
	status = H5Gclose(set_group);
#endif

	//////////////////////////////////////////////////////////////////////////////////////////////////////
	//    Model parameters                                                                              //
	//////////////////////////////////////////////////////////////////////////////////////////////////////
	//Specify data type
	dtype  = H5Tcreate(H5T_COMPOUND, sizeof(struct parameters));
	status = H5Tinsert(dtype, "M0",       HOFFSET(struct parameters, M0),       H5T_NATIVE_DOUBLE);
	status = H5Tinsert(dtype, "Epsilon0", HOFFSET(struct parameters, Epsilon0), H5T_NATIVE_DOUBLE);
	status = H5Tinsert(dtype, "Beta0",    HOFFSET(struct parameters, Beta0),    H5T_NATIVE_DOUBLE);
	status = H5Tinsert(dtype, "Gamma0",   HOFFSET(struct parameters, Gamma0),   H5T_NATIVE_DOUBLE);
	status = H5Tinsert(dtype, "MZ",       HOFFSET(struct parameters, MZ),       H5T_NATIVE_DOUBLE);
	status = H5Tinsert(dtype, "EpsilonZ", HOFFSET(struct parameters, EpsilonZ), H5T_NATIVE_DOUBLE);
	status = H5Tinsert(dtype, "BetaZ",    HOFFSET(struct parameters, BetaZ),    H5T_NATIVE_DOUBLE);
	status = H5Tinsert(dtype, "GammaZ",   HOFFSET(struct parameters, GammaZ),   H5T_NATIVE_DOUBLE);
	status = H5Tinsert(dtype, "Fesc",     HOFFSET(struct parameters, Fesc),     H5T_NATIVE_DOUBLE);
	status = H5Tinsert(dtype, "Fstrip",   HOFFSET(struct parameters, Fstrip),   H5T_NATIVE_DOUBLE);
	status = H5Tinsert(dtype, "Tau0",     HOFFSET(struct parameters, Tau0),     H5T_NATIVE_DOUBLE);
	status = H5Tinsert(dtype, "TauS",     HOFFSET(struct parameters, TauS),     H5T_NATIVE_DOUBLE);
	status = H5Tinsert(dtype, "TauD",     HOFFSET(struct parameters, TauD),     H5T_NATIVE_DOUBLE);
	//Specify file type
	ftype  = H5Tcreate(H5T_COMPOUND, 8 * 13);
	status = H5Tinsert(ftype, "M0",       0,  H5T_IEEE_F64LE);
	status = H5Tinsert(ftype, "Epsilon0", 8,  H5T_IEEE_F64LE);
	status = H5Tinsert(ftype, "Beta0",    16, H5T_IEEE_F64LE);
	status = H5Tinsert(ftype, "Gamma0",   24, H5T_IEEE_F64LE);
	status = H5Tinsert(ftype, "MZ",       32, H5T_IEEE_F64LE);
	status = H5Tinsert(ftype, "EpsilonZ", 40, H5T_IEEE_F64LE);
	status = H5Tinsert(ftype, "BetaZ",    48, H5T_IEEE_F64LE);
	status = H5Tinsert(ftype, "GammaZ",   56, H5T_IEEE_F64LE);
	status = H5Tinsert(ftype, "Fesc",     64, H5T_IEEE_F64LE);
	status = H5Tinsert(ftype, "Fstrip",   72, H5T_IEEE_F64LE);
	status = H5Tinsert(ftype, "Tau0",     80, H5T_IEEE_F64LE);
	status = H5Tinsert(ftype, "TauS",     88, H5T_IEEE_F64LE);
	status = H5Tinsert(ftype, "TauD",     96, H5T_IEEE_F64LE);
	//Set dimensions array
	dims1[0] = 1;
	//Set the path name for the parameters
	sprintf(path,"/Universe_%06d/Model_Parameters",iuniverse);
	//Create dataspace and dataset, and write the compound
	dspace = H5Screate_simple(1, dims1, NULL);
	dset   = H5Dcreate(stats_file, path, ftype, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dwrite(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &P);
	// Close the dataset
	status = H5Dclose(dset);
	// Close the data space
	status = H5Sclose(dspace);
	// Close the file type
	status = H5Tclose(ftype);
	// Close the data type
	status = H5Tclose(dtype);

	//////////////////////////////////////////////////////////////////////////////////////////////////////
	//    Chi2 values                                                                                   //
	//////////////////////////////////////////////////////////////////////////////////////////////////////

	//Allocate the chi2 array for All.Nobs values plus the total value
	chi2 = emalloc("CHI2", (All.Nobs + 1) * sizeof(float));
	//Compute Chi2
	get_chi2(chi2);
	//Set dimensions array
	dims1[0] = All.Nobs+1;
	//Set the path name for the chi2 output
	sprintf(path,"/Universe_%06d/Chi2",iuniverse);
	//Create dataspace and dataset, and write the dataset
	dspace = H5Screate_simple(1, dims1, NULL);
	dset   = H5Dcreate(stats_file, path, H5T_NATIVE_FLOAT, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dwrite(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, chi2);
	// Close the dataset
	status = H5Dclose(dset);
	// Close the data space
	status = H5Sclose(dspace);
	//Free chi2 array
	efree(chi2);

	//////////////////////////////////////////////////////////////////////////////////////////////////////
	//    Done with this Universe                                                                       //
	//////////////////////////////////////////////////////////////////////////////////////////////////////

	//Close universe group
	status = H5Gclose(universe_group);
	//Close the file and print to screen if writing was unsuccessful
	status = H5Fclose(stats_file);
	if (status != 0) printf("%s Task %5d failed to write to file %s\n",All.startline,ThisTask,outfname);

}
#endif


/*! \brief This function writes all chi2 values of all universes
 *
 */
void write_statistics_chi2(char *outdir, int step)
{
	if (All.OutputFormat == 2)
#ifdef HDF5_SUPPORT
		//For HDF5 the chi2 values are written directly into each universe
		return;
#else
		write_statistics_chi2_ascii(outdir, step);
#endif
	else if (All.OutputFormat == 1) write_statistics_chi2_ascii(outdir, step);
	else write_statistics_chi2_ascii(outdir, step);
}


/*! \brief This function writes all chi2 values of all universes into an ascii file
 *
 */
void write_statistics_chi2_ascii(char *outdir, int step)
{
	char outfname[NSTRING];
	int i,itask;
	float   *chi2;
	FILE    *ofp;
	MPI_Status status;

	//The master task now do the work
	if (ThisTask == MasterTask)
	{
		//Allocate the chi2 array for All.Nobs values plus the total value
		chi2 = emalloc("CHI2", (All.Nobs + 1) * sizeof(float));
		//Compute Chi2
		get_chi2(chi2);
		//Get the file name
		sprintf(outfname,"%s/chi2.out",outdir);

		//Task 0 creates the file and writes its chi2 values
		if (ThisTask ==  0)
		{ //Create Chi2^2 file
			ofp = fopen(outfname,"a");
			fprintf(ofp,"# TOTAL");
#ifdef READ_SMF
			fprintf(ofp," SMF");
#endif
#ifdef READ_FQ
			fprintf(ofp," FQ");
#endif
#ifdef READ_CSFRD
			fprintf(ofp," CSFRD");
#endif
#ifdef READ_SSFR
			fprintf(ofp," SSFR");
#endif
#ifdef READ_WP
			fprintf(ofp," WP");
#endif
			fprintf(ofp,"\n");
			fprintf(ofp,"%d %f",step,chi2[All.Nobs]);
			for (i = 0; i < All.Nobs; i++) fprintf(ofp," %f",chi2[i]);
			fprintf(ofp,"\n");
			fclose(ofp);
		}

		//Loop over all other master tasks
		for (itask = All.NTaskPerUniverse; itask < NTask; itask += All.NTaskPerUniverse)
		{ //Task 0 revceives the chi2 values and prints them to the file
			if (ThisTask == 0)
			{
				MPI_Recv(chi2, All.Nobs + 1, MPI_FLOAT, itask, TAG_PROB, MPI_COMM_WORLD, &status);
				ofp = fopen(outfname,"a");
				fprintf(ofp,"%d %f",step+itask/All.NTaskPerUniverse,chi2[All.Nobs]);
				for (i = 0; i < All.Nobs; i++) fprintf(ofp," %f",chi2[i]);
				fprintf(ofp,"\n");
				fclose(ofp);
			}
			//All other master tasks send their chi2 values to task 0
			else MPI_Ssend(chi2, All.Nobs + 1, MPI_FLOAT, 0, TAG_PROB, MPI_COMM_WORLD);
		}//Loop over master tasks done

		//Free the chi2 array
		efree(chi2);
	}

}


/*! \brief This function prints all walker statistics to either ascii or hdf5 output files
 *
 */
void write_walker_statistics(void)
{
	const int nline = 500;
	char outfname[NSTRING], infname[NSTRING], line[nline], *tmp;
	int i,j,k,itask,Nfilepar,Nfilewalk,Ncoldwalkers,Nsteps,Nleft;
	double temp;
	double *x;     //Walkers, Length: Nwalkers*Nparam
	double *y;     //Walkers, Length: NUniverses*Nparam
	FILE *ifp;

	//Allocate Arrays
	x = emalloc_movable(&x, "X",    All.Nwalkers * All.Nparam * sizeof(double));
	y = emalloc_movable(&y, "Y",    All.NUniverses * All.Nparam * sizeof(double));

	//Only the main task opens the walker file and prints to screen
	if (ThisTask == 0)
	{
		//Open the walker file
		if (All.Mode == 12)
		{
			sprintf(infname,"%s/walkers.mcmc.%d.out",All.OutputDir,All.Seed);
		}
		else if (All.Mode == 22)
		{
			sprintf(infname,"%s/walkers.hybrid.%d.out",All.OutputDir,All.Seed);
		}
		else if (All.Mode == 32)
		{
			sprintf(infname,"%s/walkers.pt.%d.out",All.OutputDir,All.Seed);
		}
		else
		{
			sprintf(line, "Mode %d not available...",All.Mode);
			endrun(line);
		}

		if(!(ifp=fopen(infname,"r")))
		{
			sprintf(line, "Cannot open file `%s' for reading the walkers.", infname);
			endrun(line);
		}

		//Get first line of walker file
		fgets(line,nline,ifp);
		sscanf(line,"#Walkers/Parameters/seed/i %d %d\n", &Nfilewalk, &Nfilepar);

		//Sanity Checks
		if (Nfilewalk != All.Nwalkers)
		{
			sprintf(line, "The walker file %s contains %d walkers (while the code needs %d)",infname,Nfilewalk,All.Nwalkers);
			endrun(line);
		}
		if (Nfilepar != All.Nparam)
		{
			sprintf(line, "The walker file %s contains %d parameters (while the code needs %d)",infname,Nfilepar,All.Nparam);
			endrun(line);
		}

		//Print what is done...
		printf("%s\n%s Computing statistics for each walker...\n",All.fullline,All.startline);
		printf("%s Reading %d Walkers from %s...\n",All.startline,All.Nwalkers,infname);

		//Initialise the number of cold walkers
		Ncoldwalkers = 0;

		//Get each line in the file
		for (k = 0; k < All.Nwalkers; k++)
		{ //Get the new line
			fgets(line,nline,ifp);
			tmp = strtok(line, " ");
			//Check for parallel tempering
			if (All.Mode == 32)
			{ //Get the temperature
				tmp = strtok(NULL, " ");
				temp = atof(tmp);
				//Get the scaling
				tmp = strtok(NULL, " ");
				//Get the acceptance count
				tmp = strtok(NULL, " ");
			} //If not parallel tempering...
			else
			{ //... set temperature to 1
				temp = 1.0;
			}
			//Skip if not cold walker
			if (temp > 1.0) continue;
			Ncoldwalkers++;
			//Get all the parameters for the walker
			for (i = 0; i < All.Nparam; i++)
			{
				tmp = strtok(NULL, " ");
				x[k*All.Nparam+i] = atof(tmp);
			}
		} //End reading walker file
		fclose(ifp);
	}

	//Broadcast the number of cold walkers from task 0
	MPI_Bcast(&Ncoldwalkers, 1, MPI_INT, 0, MPI_COMM_WORLD);

	//Reallocate walker array
	x = (double *) erealloc_movable(x, Ncoldwalkers * All.Nparam * sizeof(double));

	//Broadcast the walker array x from task 0
	MPI_Bcast(x, Ncoldwalkers * All.Nparam, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	//Calculate the number of steps that have to be taken, and the number of walkers left afterwards
	Nsteps = (int)((double)(Ncoldwalkers) / (double)(All.NUniverses));
	Nleft  = Ncoldwalkers % All.NUniverses;

	//Check format to see which files to create
	if (All.OutputFormat == 2)
#ifdef HDF5_SUPPORT
		write_walker_statistics_init_hdf5(outfname);
#else
		write_walker_statistics_init_ascii(outfname);
#endif
	else if (All.OutputFormat == 1) write_walker_statistics_init_ascii(outfname);
	else write_walker_statistics_init_ascii(outfname);

	//Go through all steps
	for (k = 0; k < Nsteps; k++)
	{ //For each step go through all universes
		for (j = 0; j < All.NUniverses; j++)
		{ //Go through each task in the universe
			for (itask = MasterTask; itask < MasterTask + All.NTaskPerUniverse; itask++)
			{ //If this task is in the current universe
				if (itask == ThisTask && j == MasterTask / All.NTaskPerUniverse)
				{//Go through all parameters
					for (i = 0; i < All.Nparam; i++)
					{ //Write the corresponding walker parameters to y
						y[j*All.Nparam+i] = x[(k*All.NUniverses+j)*All.Nparam+i];
					}//Params
				}//ThisTask
			}//TaskPerUniverse
		}//Universes
		//Calculate statistics for these universes
		lnprob(y);
		//Write all chi2 values to the output file
		write_statistics_chi2(outfname, k * All.NUniverses);
		//Go through all universes again
		for (j = 0; j < All.NUniverses; j++)
		{ //If this task is the master task of the current universe
			if (ThisTask == MasterTask && j == MasterTask / All.NTaskPerUniverse)
			{ //Check which format to use
				if (All.OutputFormat == 2)
#ifdef HDF5_SUPPORT
					write_walker_statistics_hdf5(outfname, k*All.NUniverses+j);
#else
					write_walker_statistics_ascii(outfname, k*All.NUniverses+j);
#endif
				else if (All.OutputFormat == 1) write_walker_statistics_ascii(outfname, k*All.NUniverses+j);
				else write_walker_statistics_ascii(outfname, k*All.NUniverses+j);
			}//ThisTask
			//Make the other tasks wait until writing is completed
			MPI_Barrier(MPI_COMM_WORLD);
		}//Universes
	}//Steps

	//Reallocate walker array
	y = (double *) erealloc_movable(y, Nleft * All.Nparam * sizeof(double));

	//Go through remaning universes
	for (j = 0; j < Nleft; j++)
	{ //Go through each task in the universe
		for (itask = MasterTask; itask < MasterTask + All.NTaskPerUniverse; itask++)
		{ //If this task is in the current universe
			if (itask == ThisTask && j == MasterTask / All.NTaskPerUniverse)
			{//Go through all parameters
				for (i = 0; i < All.Nparam; i++)
				{ //Write the corresponding walker parameters to y
					y[j*All.Nparam+i] = x[(Nsteps*All.NUniverses+j)*All.Nparam+i];
				}//Params
			}//ThisTask
		}//TaskPerUniverse
	}//Universes
	//Calculate statistics for these universes
	lnprob(y);
	//Write all chi2 values to the output file
	write_statistics_chi2(outfname, Nsteps * All.NUniverses);
	//Go through all universes again
	for (j = 0; j < Nleft; j++)
	{ //If this task is the master task of the current universe
		if (ThisTask == MasterTask && j == MasterTask / All.NTaskPerUniverse)
		{ //Check which format to use
			if (All.OutputFormat == 2)
#ifdef HDF5_SUPPORT
				write_walker_statistics_hdf5(outfname, k*All.NUniverses+j);
#else
				write_walker_statistics_ascii(outfname, k*All.NUniverses+j);
#endif
			else if (All.OutputFormat == 1) write_walker_statistics_ascii(outfname, k*All.NUniverses+j);
			else write_walker_statistics_ascii(outfname, k*All.NUniverses+j);
		}//ThisTask
		//Make the other tasks wait until writing is completed
		MPI_Barrier(MPI_COMM_WORLD);
	}//Universes

	//Free walker arrays
	efree_movable(y);
	efree_movable(x);

}


/*! \brief This function generates the folder for the walker statistiscs ascii files if needed and set the directory path
 *
 */
void write_walker_statistics_init_ascii(char *outfname)
{
	FILE *ofp;

	//Set the directory name for all tasks
	sprintf(outfname,"%s/walkerstats/%d",All.OutputDir,All.Seed);

	//Task 0 creates new directory if needed
	if (ThisTask == 0)
	{ //Print what's done to the screen
		printf("%s\n%s Writing walker statistics in folder %s\n",All.fullline,All.startline,outfname);
		//Get name of directory
		sprintf(outfname,"%s/walkerstats",All.OutputDir);
		//Create walkerstats directory if needed
		if (mkdir(outfname, 02755) == 0)
		{
			if (All.verbose >= VERBOSE_MIN) printf("%s\n%s Creating walker statistics directory '%s'\n",All.fullline,All.startline,outfname);
		}
		//Get name of directory
		sprintf(outfname,"%s/walkerstats/%d",All.OutputDir,All.Seed);
		//Create seed name directory if needed
		if (mkdir(outfname, 02755) == 0)
		{
			if (All.verbose >= VERBOSE_MIN) printf("%s\n%s Creating walker statistics directory '%s'\n",All.fullline,All.startline,outfname);
		}

		//Create the chi2 file
		sprintf(outfname,"%s/walkerstats/%d/chi2.out",All.OutputDir,All.Seed);
		ofp = fopen(outfname,"w");
		fclose(ofp);

		//Reset the directory name
		sprintf(outfname,"%s/walkerstats/%d",All.OutputDir,All.Seed);

	}//End task 0

}


#ifdef HDF5_SUPPORT
/*! \brief This function generates the HDF5 file for the walker statistiscs and set the file name
 *
 */
void write_walker_statistics_init_hdf5(char *outfname)
{
	hid_t stats_file;
	herr_t status;

	//Set the HDF5 file name
	if (All.Mode == 12) sprintf(outfname,"%s/stats.mcmc.%d.h5",  All.OutputDir,All.Seed);
	else if (All.Mode == 22) sprintf(outfname,"%s/stats.hybrid.%d.h5",All.OutputDir,All.Seed);
	else if (All.Mode == 32) sprintf(outfname,"%s/stats.pt.%d.h5",    All.OutputDir,All.Seed);
	else sprintf(outfname,"%s/stats.%d.h5",       All.OutputDir,All.Seed);
	//Create file
	if (ThisTask == 0)
	{ //Print what's done to the screen
		printf("%s\n%s Writing walker statistics to file %s\n",All.fullline,All.startline,outfname);
		//Create statistics file by opening and closing
		stats_file = H5Fcreate(outfname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
		status = H5Fclose(stats_file);
		//Warn if there are problems
		if (status != 0) printf("%s Task 0 failed to create file %s\n",All.startline,outfname);
	}

}
#endif


/*! \brief This function calls write_statistics_ascii() and prints to screen
 *
 */
void write_walker_statistics_ascii(char *outfname, int iuniverse)
{ //Print to screen what is done (use ascii)
	printf("%s Writing Universe %3d to folder %s\n",All.startline,iuniverse,outfname);
	//Write this universe to the HDF5 file
	write_statistics_ascii(outfname, iuniverse);
}


#ifdef HDF5_SUPPORT
/*! \brief This function calls write_statistics_hdf5() and prints to screen
 *
 */
void write_walker_statistics_hdf5(char *outfname, int iuniverse)
{ //Print to screen what is done (use HDF5)
	printf("%s Writing Universe %3d to file %s\n",All.startline,iuniverse,outfname);
	//Write this universe to the HDF5 file
	write_statistics_hdf5(outfname, iuniverse);
}
#endif
