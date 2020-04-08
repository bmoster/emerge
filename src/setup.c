///////////////////////////////////////////////////////////////////////////////////////////////////
// Emerge code - File setup.c                                                                    //
// Parts of these functions have been adapted from the GADGET code developed by Volker Springel  //
///////////////////////////////////////////////////////////////////////////////////////////////////
///
/// \file setup.c
/// \brief Contains functions that initialise the emerge code and read the input parameters
///
/// This file contains all functions that are needed to read the parameters and initialuse all
/// variables that are needed to run emerge. It also contains all functions needed to finalise
/// the code.
///
///////////////////////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <sys/stat.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "allvars.h"
#include "proto.h"


/*! \brief This function parses the parameter file.
 *
 *  Each parameter is defined by a keyword (`tag'), and can be either
 *  of type douple, int, or character string. Three arrays containing the name,
 *  type and address of the parameter are filled first. The routine then parses
 *  the parameter file and fills the referenced variables. The routine makes sure that
 *  each parameter appears exactly once in the parameter file, otherwise
 *  error messages are produced that complain about the missing parameters.
 *  Basic checks are performed on the supplied parameters in the end.
 *
 *  \param fname The file name of the parameter file
 */
void read_parameterfile(char *fname)
{

#define FLOAT 1        ///< Tag for floating point numbers
#define STRING 2       ///< Tag for strings
#define INT 3          ///< Tag for integers
#define MAXTAGS 300    ///< Maximum number of entries the parameter file can contain

	FILE *fd;
	char buf[NSTRING], buf1[NSTRING], buf2[NSTRING], buf3[NSTRING];
	int i, j, nt;
	int id[MAXTAGS];
	void *addr[MAXTAGS];
	char tag[MAXTAGS][50];
	int errorFlag = 0;

	// Check sizes of types on this machine
	if(sizeof(long long) != 8) endrun("Type `long long' is not 64 bit on this platform.");
	if(sizeof(int) != 4) endrun("Type `int' is not 32 bit on this platform.");
	if(sizeof(float) != 4) endrun("Type `float' is not 32 bit on this platform.");
	if(sizeof(double) != 8) endrun("Type `double' is not 64 bit on this platform.");

	// Read parameter file on process 0
	if(ThisTask == 0)
	{

		nt = 0;

		strcpy(tag[nt], "TreefileName");
		addr[nt] = All.treefile_name;
		id[nt++] = STRING;

		strcpy(tag[nt], "MaxMemSize");
		addr[nt] = &All.MaxMemSize;
		id[nt++] = INT;

		strcpy(tag[nt], "BufferSize");
		addr[nt] = &All.BufferSize;
		id[nt++] = INT;

		strcpy(tag[nt], "UniversesInParallel");
		addr[nt] = &All.NUniverses;
		id[nt++] = INT;

		strcpy(tag[nt], "NumOutputFiles");
		addr[nt] = &All.NumOutputFiles;
		id[nt++] = INT;

		strcpy(tag[nt], "NumFilesInParallel");
		addr[nt] = &All.NumFilesInParallel;
		id[nt++] = INT;

		strcpy(tag[nt], "OutputFormat");
		addr[nt] = &All.OutputFormat;
		id[nt++] = INT;

		strcpy(tag[nt], "MCMCseed");
		addr[nt] = &All.Seed;
		id[nt++] = INT;

		strcpy(tag[nt], "NumberOfMCMCWalkers");
		addr[nt] = &All.Nwalkers;
		id[nt++] = INT;

		strcpy(tag[nt], "ModelName");
		addr[nt] = All.model_name;
		id[nt++] = STRING;

		strcpy(tag[nt], "MCMCScaleParameter");
		addr[nt] = &All.mcmca;
		id[nt++] = FLOAT;

		strcpy(tag[nt], "ChainTemperature");
		addr[nt] = &All.temperature;
		id[nt++] = FLOAT;

		strcpy(tag[nt], "TimeLimit");
		addr[nt] = &All.timelimit;
		id[nt++] = FLOAT;

		strcpy(tag[nt], "OutputRedshifts");
		addr[nt] = All.output_redshifts;
		id[nt++] = STRING;

		strcpy(tag[nt], "OutputMassThreshold");
		addr[nt] = &All.minmass;
		id[nt++] = FLOAT;

		strcpy(tag[nt], "Verbose");
		addr[nt] = &All.verbose;
		id[nt++] = INT;

		strcpy(tag[nt], "HubbleParam");
		addr[nt] = &All.h_100;
		id[nt++] = FLOAT;

		strcpy(tag[nt], "Omega0");
		addr[nt] = &All.Omega_0;
		id[nt++] = FLOAT;

		strcpy(tag[nt], "OmegaLambda");
		addr[nt] = &All.Omega_Lambda_0;
		id[nt++] = FLOAT;

		strcpy(tag[nt], "OmegaBaryon");
		addr[nt] = &All.Omega_Baryon_0;
		id[nt++] = FLOAT;

		strcpy(tag[nt], "BoxSize");
		addr[nt] = &All.Lbox;
		id[nt++] = FLOAT;

		strcpy(tag[nt], "Mstar_min");
		addr[nt] = &All.mstarmin;
		id[nt++] = FLOAT;

		strcpy(tag[nt], "Mstar_max");
		addr[nt] = &All.mstarmax;
		id[nt++] = FLOAT;

		strcpy(tag[nt], "Delta_Mstar");
		addr[nt] = &All.dmstar;
		id[nt++] = FLOAT;

		strcpy(tag[nt], "Observation_Error_0");
		addr[nt] = &All.obssigma0;
		id[nt++] = FLOAT;

		strcpy(tag[nt], "Observation_Error_z");
		addr[nt] = &All.obssigmaz;
		id[nt++] = FLOAT;

		strcpy(tag[nt], "UnitLength_in_Mpc");
		addr[nt] = &All.x_unit;
		id[nt++] = FLOAT;

		strcpy(tag[nt], "UnitTime_in_yr");
		addr[nt] = &All.t_unit;
		id[nt++] = FLOAT;

		strcpy(tag[nt], "UnitMass_in_Msun");
		addr[nt] = &All.m_unit;
		id[nt++] = FLOAT;

		strcpy(tag[nt], "Eff_MassPeak");
		addr[nt] = &All.M0;
		id[nt++] = FLOAT;

		strcpy(tag[nt], "Eff_Normalisation");
		addr[nt] = &All.Epsilon0;
		id[nt++] = FLOAT;

#ifdef SFE_MPEAK_ZEVOLV
		strcpy(tag[nt], "Eff_MassPeak_Z");
		addr[nt] = &All.MZ;
		id[nt++] = FLOAT;
#endif

#ifdef SFE_NORM_ZEVOLV
		strcpy(tag[nt], "Eff_Normalisation_Z");
		addr[nt] = &All.EpsilonZ;
		id[nt++] = FLOAT;
#endif

		strcpy(tag[nt], "Eff_LowMassSlope");
		addr[nt] = &All.Beta0;
		id[nt++] = FLOAT;

#ifndef SFE_SAME_SLOPE
		strcpy(tag[nt], "Eff_HighMassSlope");
		addr[nt] = &All.Gamma0;
		id[nt++] = FLOAT;
#endif

#ifdef SFE_BETA_ZEVOLV
		strcpy(tag[nt], "Eff_LowMassSlope_Z");
		addr[nt] = &All.BetaZ;
		id[nt++] = FLOAT;
#endif

#ifdef SFE_GAMMA_ZEVOLV
#ifndef SFE_SAME_SLOPE
		strcpy(tag[nt], "Eff_HighMassSlope_Z");
		addr[nt] = &All.GammaZ;
		id[nt++] = FLOAT;
#endif
#endif

		strcpy(tag[nt], "Fraction_Escape_ICM");
		addr[nt] = &All.Fesc;
		id[nt++] = FLOAT;

		strcpy(tag[nt], "Fraction_Stripping");
		addr[nt] = &All.Fstrip;
		id[nt++] = FLOAT;

		strcpy(tag[nt], "Timescale_Quenching");
		addr[nt] = &All.Tau0;
		id[nt++] = FLOAT;

#ifdef SAT_QUENCH_MASS_DEPENDENT
		strcpy(tag[nt], "Slope_Quenching");
		addr[nt] = &All.TauS;
		id[nt++] = FLOAT;
#endif

#ifdef SAT_SFR_EXP_DECAY
		strcpy(tag[nt], "Decay_Quenching");
		addr[nt] = &All.TauD;
		id[nt++] = FLOAT;
#endif

		strcpy(tag[nt], "Eff_MassPeak_Range");
		addr[nt] = &All.DeltaM0;
		id[nt++] = FLOAT;

		strcpy(tag[nt], "Eff_Normalisation_Range");
		addr[nt] = &All.DeltaEpsilon0;
		id[nt++] = FLOAT;

#ifdef SFE_MPEAK_ZEVOLV
		strcpy(tag[nt], "Eff_MassPeak_Z_Range");
		addr[nt] = &All.DeltaMZ;
		id[nt++] = FLOAT;
#endif

#ifdef SFE_NORM_ZEVOLV
		strcpy(tag[nt], "Eff_Normalisation_Z_Range");
		addr[nt] = &All.DeltaEpsilonZ;
		id[nt++] = FLOAT;
#endif

		strcpy(tag[nt], "Eff_LowMassSlope_Range");
		addr[nt] = &All.DeltaBeta0;
		id[nt++] = FLOAT;

#ifndef SFE_SAME_SLOPE
		strcpy(tag[nt], "Eff_HighMassSlope_Range");
		addr[nt] = &All.DeltaGamma0;
		id[nt++] = FLOAT;
#endif

#ifdef SFE_BETA_ZEVOLV
		strcpy(tag[nt], "Eff_LowMassSlope_Z_Range");
		addr[nt] = &All.DeltaBetaZ;
		id[nt++] = FLOAT;
#endif

#ifdef SFE_GAMMA_ZEVOLV
#ifndef SFE_SAME_SLOPE
		strcpy(tag[nt], "Eff_HighMassSlope_Z_Range");
		addr[nt] = &All.DeltaGammaZ;
		id[nt++] = FLOAT;
#endif
#endif

		strcpy(tag[nt], "Fraction_Escape_ICM_Range");
		addr[nt] = &All.DeltaFesc;
		id[nt++] = FLOAT;

		strcpy(tag[nt], "Fraction_Stripping_Range");
		addr[nt] = &All.DeltaFstrip;
		id[nt++] = FLOAT;

		strcpy(tag[nt], "Timescale_Quenching_Range");
		addr[nt] = &All.DeltaTau0;
		id[nt++] = FLOAT;

#ifdef SAT_QUENCH_MASS_DEPENDENT
		strcpy(tag[nt], "Slope_Quenching_Range");
		addr[nt] = &All.DeltaTauS;
		id[nt++] = FLOAT;
#endif

#ifdef SAT_SFR_EXP_DECAY
		strcpy(tag[nt], "Decay_Quenching_Range");
		addr[nt] = &All.DeltaTauD;
		id[nt++] = FLOAT;
#endif

#ifdef READ_SMF
		strcpy(tag[nt], "SMFfileName");
		addr[nt] = All.smffile_name;
		id[nt++] = STRING;
#endif

#ifdef READ_FQ
		strcpy(tag[nt], "FQfileName");
		addr[nt] = All.fqfile_name;
		id[nt++] = STRING;
#endif

#ifdef READ_CSFRD
		strcpy(tag[nt], "CSFRDfileName");
		addr[nt] = All.csfrdfile_name;
		id[nt++] = STRING;
#endif

#ifdef READ_SSFR
		strcpy(tag[nt], "SSFRfileName");
		addr[nt] = All.ssfrfile_name;
		id[nt++] = STRING;
#endif

#ifdef READ_WP
		strcpy(tag[nt], "WPfileName");
		addr[nt] = All.wpfile_name;
		id[nt++] = STRING;
#endif

#ifdef GLOBAL_SIGMA_SMF
		strcpy(tag[nt], "Global_Sigma_SMF_LZ");
		addr[nt] = &All.GlobalSigmaSmfLz;
		id[nt++] = FLOAT;

		strcpy(tag[nt], "Global_Sigma_SMF_HZ");
		addr[nt] = &All.GlobalSigmaSmfHz;
		id[nt++] = FLOAT;
#endif

#ifdef GLOBAL_SIGMA_FQ
		strcpy(tag[nt], "Global_Sigma_FQ");
		addr[nt] = &All.GlobalSigmaFq;
		id[nt++] = FLOAT;
#endif

#ifdef GLOBAL_SIGMA_CSFRD
		strcpy(tag[nt], "Global_Sigma_CSFRD");
		addr[nt] = &All.GlobalSigmaCsfrd;
		id[nt++] = FLOAT;
#endif

#ifdef GLOBAL_SIGMA_SSFR
		strcpy(tag[nt], "Global_Sigma_SSFR");
		addr[nt] = &All.GlobalSigmaSsfr;
		id[nt++] = FLOAT;
#endif

#ifdef GLOBAL_SIGMA_WP
		strcpy(tag[nt], "Global_Sigma_WP");
		addr[nt] = &All.GlobalSigmaWp;
		id[nt++] = FLOAT;
#endif

#ifdef WRITE_MAINBRANCH
		strcpy(tag[nt], "MainBranchMasses");
		addr[nt] = All.output_mass_mb;
		id[nt++] = STRING;

		strcpy(tag[nt], "MainBranchBinSize");
		addr[nt] = All.output_mass_mb_bin;
		id[nt++] = STRING;

		strcpy(tag[nt], "MainBranchMassType");
		addr[nt] = &All.MainBranchMassType;
		id[nt++] = INT;

		strcpy(tag[nt], "MainBranchRedshift");
		addr[nt] = &All.mainBranchRedshift;
		id[nt++] = FLOAT;
#endif

		if((fd = fopen(fname, "r")))
		{
			while(!feof(fd))
			{
				buf[0] = 0;
				fgets(buf, 200, fd);

				if(sscanf(buf, "%s%s%s", buf1, buf2, buf3) < 2)
					continue;

				if(buf1[0] == '%')
					continue;

				for(i = 0, j = -1; i < nt; i++)
					if(strcmp(buf1, tag[i]) == 0)
					{
						j = i;
						tag[i][0] = 0;
						break;
					}

				if(j >= 0)
				{
					switch (id[j])
					{
					case FLOAT:
						*((double *) addr[j]) = atof(buf2);
						break;
					case STRING:
						strcpy(addr[j], buf2);
						break;
					case INT:
						*((int *) addr[j]) = atoi(buf2);
						break;
					}
					printf("%s %-30s %s\n",All.startline,buf1,buf2);
				}
				else
				{
					fprintf(stdout, "%s Error in file %s:   Tag '%s' not allowed or multiple defined.\n", All.startline, fname, buf1);
					errorFlag = 1;
				}
			}
			fclose(fd);

		}
		else
		{
			fprintf(stdout, "%s Parameter file %s not found.\n", All.startline, fname);
			errorFlag = 1;
		}


		for(i = 0; i < nt; i++)
		{
			if(*tag[i])
			{
				fprintf(stdout, "%s Error. I miss a value for tag '%s' in parameter file '%s'.\n", All.startline, tag[i], fname);
				errorFlag = 1;
			}
		}

		//Write name of output directory
		sprintf(All.OutputDir,"output/%s",All.model_name);

	}

	//Communicate ErrorFlag
	MPI_Bcast(&errorFlag, 1, MPI_INT, 0, MPI_COMM_WORLD);

	if(errorFlag)
	{
		MPI_Finalize();
		exit(0);
	}

	//Communicate all global parameters to the other processes
	MPI_Bcast(&All, sizeof(struct global_data), MPI_BYTE, 0, MPI_COMM_WORLD);

	//Make sure the number of tasks is a multiple of the number of universes
	if(NTask % All.NUniverses != 0)
	{
		if(ThisTask == 0) printf("%s UniversesInParallel = %d, Number of processors = %d\n",All.startline,All.NUniverses,NTask);
		endrun("Number of processors MUST be a multiple of number of UniversesInParallel");
	}

	//Make sure the number of walkers is a multiple of the number of universes
	if(All.Nwalkers % (All.NUniverses*2) != 0 && (All.Mode==1 || All.Mode==11))
	{
		if(ThisTask == 0) printf("%s UniversesInParallel = %d, Number of walkers = %d\n",All.startline,All.NUniverses,All.Nwalkers);
		endrun("Number of walkers MUST be a multiple of twice the number of UniversesInParallel");
	}
	if(All.Nwalkers % (All.NUniverses) != 0 && (All.Mode==2 || All.Mode==21 || All.Mode==3 || All.Mode==31))
	{
		if(ThisTask == 0) printf("%s UniversesInParallel = %d, Number of walkers = %d\n",All.startline,All.NUniverses,All.Nwalkers);
		endrun("Number of walkers MUST be a multiple of the number of UniversesInParallel");
	}

#ifndef HDF5_SUPPORT
	if (All.OutputFormat == 2)
	{
		All.OutputFormat = 1;
		if (ThisTask == 0) printf("%s\n%s The output format cannot be set to 2 as HDF5 support is not enabled.\n",All.fullline,All.startline);
		if (ThisTask == 0) printf("%s If HDF5 output files are required the code must be compiled with the option HDF5_SUPPORT.\n",All.startline);
		if (ThisTask == 0) printf("%s Setting the output format to 1 (normal ascii output).\n",All.startline);
	}
#endif

	if (All.OutputFormat < 1 || All.OutputFormat > 2)
	{
		All.OutputFormat = 1;
		if (ThisTask == 0) printf("%s\n%s The output format is out of the defined range.\n",All.fullline,All.startline);
		if (ThisTask == 0) printf("%s Setting the output format to 1 (normal ascii output).\n",All.startline);
	}

#undef FLOAT
#undef STRING
#undef INT
#undef MAXTAGS
}


/*! \brief This function initialises most global variables
 *
 *  The number of tasks per universe are set and the Master task is stored. The memory manager is initialised.
 *  The starting time, the universal baryon fraction, the box length, the hubble time, and the minimum mass for
 *  the output are set. The starting parameters are set and the number of parameters is counted. All starting
 *  parameter values are copied to the #parameters structure #P. If not yet present, the output folder is created.
 *  The log files are opened and the random number generator or table is set. The stellar mass array for the
 *  interpolation #Mstar is computed. Finally, all global parameters that are the same on all tasks #All are
 *  broadcast from task 0.
 */

void setup(void){

	int i;
	char buf[NSTRING], tmp1[NSTRING], *tmp2;

	//Print a line
	if (ThisTask == 0) printf("%s\n",All.fullline);

	//Set the number of tasks per universe
	All.NTaskPerUniverse = NTask / All.NUniverses;

	//Identify master task which is the main task for each universe
	MasterTask = (ThisTask/All.NTaskPerUniverse) * All.NTaskPerUniverse;

	//Initialise the memory manager
	malloc_init();

	//Initialise the time
	All.starttime  = second();
	All.timenow    = second() - All.starttime;

	//Define some variables
	All.f_baryon = All.Omega_Baryon_0/All.Omega_0;
	All.Lbox /= All.h_100;
	All.Hubbletime = 9.7776E09/All.h_100;

	// Minimum mass for galaxy out put...
	All.minmass = pow(10.,All.minmass)/All.m_unit;

	// Set number of parameters based on Initial Range (not a parameter to fit when Delta < 0)
	All.Nparam = 0;
	if (All.DeltaM0       >= 0.0) All.Nparam++;
	if (All.DeltaEpsilon0 >= 0.0) All.Nparam++;
	if (All.DeltaBeta0    >= 0.0) All.Nparam++;
	if (All.DeltaGamma0   >= 0.0) All.Nparam++;
	if (All.DeltaFesc     >= 0.0) All.Nparam++;
	if (All.DeltaFstrip   >= 0.0) All.Nparam++;
	if (All.DeltaTau0     >= 0.0) All.Nparam++;

#ifdef SFE_MPEAK_ZEVOLV
	if (All.DeltaMZ       >= 0.0) All.Nparam++;
#else
	All.MZ                 = 0.0;
	All.DeltaMZ            = -1.0;
#endif

#ifdef SFE_NORM_ZEVOLV
	if (All.DeltaEpsilonZ >= 0.0) All.Nparam++;
#else
	All.EpsilonZ           = 0.0;
	All.DeltaEpsilonZ      = -1.0;
#endif

#ifdef SFE_BETA_ZEVOLV
	if (All.DeltaBetaZ    >= 0.0) All.Nparam++;
#else
	All.BetaZ              = 0.0;
	All.DeltaBetaZ         = -1.0;
#endif

#ifdef SFE_GAMMA_ZEVOLV
	if (All.DeltaGammaZ   >= 0.0) All.Nparam++;
#else
	All.GammaZ             = 0.0;
	All.DeltaGammaZ        = -1.0;
#endif

#ifdef SFE_SAME_SLOPE
#ifdef SFE_GAMMA_ZEVOLV
	if (All.DeltaGamma0   >= 0.0) All.Nparam--;
	if (All.DeltaGammaZ   >= 0.0) All.Nparam--;
	All.Gamma0             = All.Beta0;
	All.GammaZ             = All.BetaZ;
	All.DeltaGamma0        = -1.0;
	All.DeltaGammaZ        = -1.0;
#else
	if (All.DeltaGamma0   >= 0.0) All.Nparam--;
	All.Gamma0             = All.Beta0;
	All.DeltaGamma0        = -1.0;
#endif
#endif

#ifdef SAT_SFR_EXP_DECAY
	if (All.DeltaTauD     >= 0.0) All.Nparam++;
#else
	All.TauD               = 0.0;
	All.DeltaTauD          = -1.0;
#endif

#ifdef SAT_QUENCH_MASS_DEPENDENT
	if (All.DeltaTauS     >= 0.0) All.Nparam++;
#else
	All.TauS               = 0.0;
	All.DeltaTauS          = -1.0;
#endif

	//Initialise local parameters with value read from parameter file
	P.M0       = All.M0;
	P.Epsilon0 = All.Epsilon0;
	P.Beta0    = All.Beta0;
	P.Gamma0   = All.Gamma0;
	P.MZ       = All.MZ;
	P.EpsilonZ = All.EpsilonZ;
	P.BetaZ    = All.BetaZ;
	P.GammaZ   = All.GammaZ;
	P.Fesc     = All.Fesc;
	P.Fstrip   = All.Fstrip;
	P.Tau0     = All.Tau0;
	P.TauS     = All.TauS;
	P.TauD     = All.TauD;

	//Print number of free parameters
	if (ThisTask == 0) printf("%s The total number of free parameters is %d\n",All.startline,All.Nparam);

	// Create Output directory
	if (ThisTask == 0)
	{
		mkdir("output", 02755);
		if (mkdir(All.OutputDir, 02755) == 0)
		{
			if (All.verbose >= VERBOSE_MIN) printf("%s Creating output directory '%s'\n",All.startline,All.OutputDir);
		}
#ifdef WRITE_GALAXY_CATALOG
		sprintf(buf, "%s/galaxies", All.OutputDir);
		if (mkdir(buf, 02755) == 0)
		{
			if (All.verbose >= VERBOSE_MIN) printf("%s Creating galaxy catalogue directory '%s'\n",All.startline,buf);
		}
#endif
#ifdef WRITE_HALO_CATALOG
		sprintf(buf, "%s/haloes", All.OutputDir);
		if (mkdir(buf, 02755) == 0)
		{
			if (All.verbose >= VERBOSE_MIN) printf("%s Creating halo catalogue directory '%s'\n",All.startline,buf);
		}
#endif
#ifdef WRITE_MAINBRANCH
		sprintf(buf, "%s/mainbranches", All.OutputDir);
		if (mkdir(buf, 02755) == 0)
		{
			if (All.verbose >= VERBOSE_MIN) printf("%s Creating mainbranch directory '%s'\n",All.startline,buf);
		}
#endif
	}

	//Open all log files
	open_logfiles();

	//Initialise random number generator or table
	init_random_numbers();

	//Initialise the number of observational statistics
	All.Nobs = 0;

	//Initialise the number of time steps to zero
	//If at some point it is larger we know the scale factors have been computed
	All.NTimesteps = 0;

	//Create the stellar mass array
	All.Nmstar = (int)((All.mstarmax - All.mstarmin) / All.dmstar) + 1;
	//Check that the number of mass bins is not too low (at least 4 for the splines - though in practice much higher)
	if (All.Nmstar < 4)
	{
		sprintf(buf, "The number of mass bins is %d which is lower than the minimim allowed number of 4.", All.Nmstar);
		endrun(buf);
	}
	Mstar = emalloc("MStar", All.Nmstar * sizeof(float));
	for (i = 0; i < All.Nmstar; i++) Mstar[i] = (float)(All.mstarmin) + ((float)(i)+0.5) * (float)(All.dmstar);

	//Initialise the number of output reshifts
	All.Noutputredshifts = i = 1;
	//Parse the output redshift string for the first time
	strcpy(tmp1,All.output_redshifts);
	tmp2 = strtok(tmp1, ",");
	//Parse the output redshift strings until all values have been parsed and increment the number of output redshifts
	while ((tmp2 = strtok(NULL, ","))) All.Noutputredshifts++;
	//Allocate the arrays that store the output redshifts and the corresponding snapshot numbers
	OutputRedshifts = emalloc("OutputRedshifts", All.Noutputredshifts * sizeof(float));
	Output_iscale   = emalloc("Output_iscale",   All.Noutputredshifts * sizeof(int));
	//Parse the output redshift string for the first time (again)
	strcpy(tmp1,All.output_redshifts);
	tmp2 = strtok(tmp1, ",");
	//Read the first output redshift
	OutputRedshifts[0] = atof(tmp2);
	if (OutputRedshifts[0] < 0.0) endrun("The output redshifts cannot be negative!");
	//Parse the output redshift strings until all values have been parsed and increment the index
	while ((tmp2 = strtok(NULL, ",")))
	{
		OutputRedshifts[i] = atof(tmp2);
		if (OutputRedshifts[i] < 0.0) endrun("The output redshifts cannot be negative!");
		i++;
	}

#ifdef WRITE_MAINBRANCH
	//Initialise the number of masses for the main branch histories
	All.Noutputbranch = i = 1;
	//Parse the mass string for the first time
	strcpy(tmp1,All.output_mass_mb);
	tmp2 = strtok(tmp1, ",");
	//Parse the mass string until all values have been parsed and increment the number of output masses
	while ((tmp2 = strtok(NULL, ","))) All.Noutputbranch++;
	//Allocate the arrays that store the output redshifts and the corresponding snapshot numbers
	MainBranchMasses = emalloc("MainBranchMasses", All.Noutputbranch * sizeof(float));
	//Parse the mass string for the first time (again)
	strcpy(tmp1,All.output_mass_mb);
	tmp2 = strtok(tmp1, ",");
	//Read the first mass
	MainBranchMasses[0] = atof(tmp2);
	if (MainBranchMasses[0] < 0.0) endrun("The main branch mass cannot be negative!");
	//Parse the mass string until all values have been parsed and increment the index
	while ((tmp2 = strtok(NULL, ",")))
	{
		MainBranchMasses[i] = atof(tmp2);
		if (MainBranchMasses[i] < 0.0) endrun("The main branch mass cannot be negative!");
		i++;
	}

	//Initialise the number of mass bin sizes for the main branch histories
	i = 1;
	//Parse the mass string for the first time
	strcpy(tmp1,All.output_mass_mb_bin);
	tmp2 = strtok(tmp1, ",");
	//Parse the mass string until all values have been parsed and increment the number of output masses
	while ((tmp2 = strtok(NULL, ","))) i++;
	//Check if the number of masses is the same as the number of bin sizes and abort if needed
	if (i != All.Noutputbranch) endrun("The number of mass bin sizes needs to be equal to the number of masses for the main branch output!");
	//Reset i
	i = 1;
	//Allocate the arrays that store the output redshifts and the corresponding snapshot numbers
	MainBranchBinSize = emalloc("MainBranchBinSize", All.Noutputbranch * sizeof(float));
	//Parse the mass string for the first time (again)
	strcpy(tmp1,All.output_mass_mb_bin);
	tmp2 = strtok(tmp1, ",");
	//Read the first mass
	MainBranchBinSize[0] = atof(tmp2);
	if (MainBranchBinSize[0] < 0.0) endrun("The main branch mass bin size cannot be negative!");
	//Parse the mass string until all values have been parsed and increment the index
	while ((tmp2 = strtok(NULL, ",")))
	{
		MainBranchBinSize[i] = atof(tmp2);
		if (MainBranchBinSize[i] < 0.0) endrun("The main branch mass bin size cannot be negative!");
		i++;
	}

	//If main branch mass type has a value that is not pre-defined set it to the default of 0, i.e. halo mass
	if (All.MainBranchMassType < 0 && All.MainBranchMassType > 1)
	{
		All.MainBranchMassType = 0;
		if (ThisTask == 0) printf("%s Main Branch Mass Type is out of bounds. Selecting default (halo mass).\n",All.startline);
	}

	//If main branch redshift is negative, set it to redshift 0
	if (All.mainBranchRedshift < 0.0)
	{
		All.mainBranchRedshift = 0.0;
		if (ThisTask == 0) printf("%s Main Branch Redshift cannot be negative. Setting it to z = 0.\n",All.startline);
	}
#endif

	//Communicate all global parameters to the other processes
	MPI_Bcast(&All, sizeof(struct global_data), MPI_BYTE, 0, MPI_COMM_WORLD);

}


/*! \brief Frees all memory at the end of the run.
 *
 *   This function frees all memory that is still allocated and closes the log files.
 */
void finalize (void)
{
	efree_movable(MassLeft);
	efree_movable(PopAge);
	efree_movable(DynTime);
	efree_movable(Timestep);
	efree_movable(CosmicTime);
	efree_movable(ScaleFactor);
	efree_movable(H);
#ifdef COMPUTE_ICM
	efree_movable(NhalosInForest);
	efree_movable(NtreesInForest);
#endif
	efree_movable(OffsetHalos);
	efree_movable(NhalosInTree);
#ifdef READ_WP
	efree(Wp);
	efree(WpSet);
#endif
#ifdef READ_SSFR
	efree(Ssfr);
	efree(SsfrSet);
#endif
#ifdef READ_CSFRD
	efree(Csfrd);
	efree(CsfrdSet);
#endif
#ifdef READ_FQ
	efree(Fq);
	efree(FqSet);
#endif
#ifdef READ_SMF
	efree(Smf);
	efree(SmfSet);
#endif
#ifdef WRITE_MAINBRANCH
	efree(MainBranchBinSize);
	efree(MainBranchMasses);
#endif
	efree(Output_iscale);
	efree(OutputRedshifts);
	efree(Mstar);
#if (RANDOM_NUMBER_TABLE > 0)
	efree(RndTableUniform);
	efree(RndTableGaussian);
#endif
	close_logfiles();
	gsl_rng_free(rng_gaussian);
}


/*! \brief Open files for logging.
 *
 *   This function opens various log-files that report on the status and performance. Upon restart, the code
 *  will append to these files.
 */
void open_logfiles(void)
{
	char mode[2], buf[NSTRING], msg[NSTRING];

	//if(RestartFlag == 0)
	strcpy(mode, "w");
	//else
	//strcpy(mode, "a");

	if(ThisTask != 0)             // only the root processors writes to the log files
		return;

	sprintf(buf, "%s/%s", All.OutputDir, "memory.txt");
	if(!(FpMemory = fopen(buf, mode)))
	{
		sprintf(msg, "error in opening file '%s'\n", buf);
		endrun(msg);
	}

}


/*! \brief  Close the global log-files.
 */
void close_logfiles(void)
{
	if(ThisTask != 0)             // only the root processors writes to the log files
		return;

	fclose(FpMemory);

}


/*! \brief This function initialises the random numbers.
 *
 *  This function initialises the random numbers. If #RANDOM_NUMBER_TABLE is selected the random number
 *  table is allocated and filled with gaussian random numbers on each task
 */
void init_random_numbers(void)
{
#if (RANDOM_NUMBER_TABLE > 0)
	int i;
#endif

	//Initialise the generators and set the seed
	rng_gaussian = gsl_rng_alloc(gsl_rng_ranlxs2);
	rng_uniform  = gsl_rng_alloc(gsl_rng_ranlxs2);
	gsl_rng_set(rng_gaussian, 140314081082 + (long)(ThisTask % All.NTaskPerUniverse));
	gsl_rng_set(rng_uniform,  160982140314 + (long)(ThisTask % All.NTaskPerUniverse));

	//If the random number table option is selected we store random numbers in a table
#if (RANDOM_NUMBER_TABLE > 0)
	if (RANDOM_NUMBER_TABLE<NRNDTABLEMIN && ThisTask == 0) printf("%s\n%s The random number tables have fewer than %d entries. Check if this is suitable!\n", All.fullline,All.startline,NRNDTABLEMIN);
	//Allocate the random number table and put numbers in it
	RndTableGaussian = emalloc("RndTableGaussian", RANDOM_NUMBER_TABLE * sizeof(float));
	RndTableUniform  = emalloc("RndTableUniform",  RANDOM_NUMBER_TABLE * sizeof(float));
	for (i = 0; i < RANDOM_NUMBER_TABLE; i++)
	{
		RndTableGaussian[i] = gsl_ran_gaussian(rng_gaussian,1.0);
		RndTableUniform[i]  = gsl_rng_uniform(rng_uniform);
	}
#endif

}


/*! \brief This function prints a nice banner at startup
 */
void print_banner(void)
{
	if (ThisTask==0)
	{
		printf("%s\n",All.fullline);
		printf("%s\n",All.startline);
		printf("%s                        ######   \n",All.startline);
		printf("%s                  ######         \n",All.startline);
		printf("%s               #####             \n",All.startline);
		printf("%s            #####  #########     \n",All.startline);
		printf("%s         ###### ##############      #########                                                \n",All.startline);
		printf("%s        ##### ####    #########     ##     ##                                                \n",All.startline);
		printf("%s      ###### ###        #######     ##        ##     ## ######## ########   ######   ########\n",All.startline);
		printf("%s     ####### ##       # #######     ##        ###   ### ##       ##     ## ##    ##  ##      \n",All.startline);
		printf("%s     ####### #       ## #######     #######   #### #### ##       ##     ## ##        ##      \n",All.startline);
		printf("%s     #######        ### ######      ##        ## ### ## ######   ########  ##   #### ######  \n",All.startline);
		printf("%s     #########    #### #####        ##        ##     ## ##       ##   ##   ##    ##  ##      \n",All.startline);
		printf("%s      ############## ######         ##     ## ##     ## ##       ##    ##  ##    ##  ##      \n",All.startline);
		printf("%s        #########  #####            ######### ##     ## ######## ##     ##  ######   ########\n",All.startline);
		printf("%s                #####            \n",All.startline);
		printf("%s            ######               \n",All.startline);
		printf("%s      ######                     \n",All.startline);
		printf("%s\n",All.startline);
	}
}
