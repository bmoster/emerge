///////////////////////////////////////////////////////////////////////////////////////////////////
// Emerge code - File allvars.h                                                                  //
///////////////////////////////////////////////////////////////////////////////////////////////////
///
/// \file allvars.h
/// \brief Contains global variables and declarations
///
/// This file declares all global variables. Further variables should be added here, and declared
/// as 'extern'. The actual existence of these variables is provided by the file 'allvars.c'.
/// To produce 'allvars.c' from 'allvars.h', do the following:
///
///   - Erase all #define statements
///   - add #include "allvars.h"
///   - delete all keywords 'extern'
///   - delete all struct definitions enclosed in {...},
///      e.g. "extern struct global_data {....} All;"
///      becomes "struct global_data All;"
///
///////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef ALLVARS_H
#define ALLVARS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <math.h>
#include <gsl/gsl_rng.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "codeoptions.h"

#define CODE      "EMERGE"             ///< Code name
#define BRANCH    "master"             ///< Branch name
#define VERSION   "1.0.2"              ///< Code version
#define SYMBOL    '#'                  ///< Symbol used to start output lines

#define NSTRING              200       ///< Number of elements used for a generic string
#define SSTRING               50       ///< Number of elements used for a short string
#define LINESIZE             100       ///< Number of characters used for the screen output
#define ALLOC_TOLERANCE      0.2       ///< Tolerance in the memory allocation
#define NINTQUAD              16       ///< Number of integration nodes for gaussian quadrature

#define ZMAX_SMFERR          4.0       ///< Maximum redshift up to which the observed stellar mass error increases
#define SSFRTHRESH           0.3       ///< Fraction of the inverse Hubble time used for the SSFR threshold
#define SIGMASMFZTHRESH      0.1       ///< Redshift that divides the SMF for the global sigma assignment

#define MAXTIMESTEPS         10000     ///< Maximum number of timesteps for the simulation
#define NRNDTABLEMIN         1000      ///< Minimum number of entries for the random number table
#define NALLOCATION          1000      ///< Number of elements that are used for dynamic allocation

#define TREEDIM              3         ///< Number of dimensions of the kd-tree used for the clustering

#if !defined(OPENMPTHREADS)
#define OPENMPTHREADS 1                ///< Default number of OpenMP threads (1: OpenMP disabled)
#endif

#if !defined(RANDOM_NUMBER_TABLE)
#define RANDOM_NUMBER_TABLE 0          ///< Default number entries for the random number tables (0: disabled)
#endif

#if !defined(WP_RBINS)
#define WP_RBINS 20                    ///< Number of radial bins for the 3D correlation functions
#endif

#if !defined(WP_RBINS_INT)
#define WP_RBINS_INT 1000              ///< Number of radial bins for the interpolation of the 3D correlation function
#endif

#if !defined(WP_RMAX)
#define WP_RMAX 0.1                    ///< Maximum radius (given as a fraction of the box size) of the 3D correlation functions
#endif

#if !defined(WP_NLEAF_MIN)
#define WP_NLEAF_MIN 4                 ///< Minimum number of galaxies to split a node of the kd-tree
#endif

#if !defined(WP_NODE_WIDTH_MIN)
#define WP_NODE_WIDTH_MIN 0.01         ///< Minimum width of each node in units of the box size
#endif

#ifndef LONGIDS
typedef unsigned int IDType;           ///< Integer type used for IDs
#else
typedef unsigned long long IDType;     ///< Integer type used for IDs if LONGIDS is selected
#endif

//
// Macro functions \cond
#define  endrun(...) {if(EndFlag==0) {char termbuf1[1000], termbuf2[1000]; sprintf(termbuf1, "%s Code was stopped on task=%d, function %s(), file %s, line %d", All.startline, ThisTask, __FUNCTION__, __FILE__, __LINE__); sprintf(termbuf2, __VA_ARGS__); printf("%s: %s\n", termbuf1, termbuf2); fflush(stdout); EndFlag=1; MPI_Abort(MPI_COMM_WORLD, 1);} exit(0);}

//Minimum and maximum of two values
#define min(X, Y)  ((X) < (Y) ? (X) : (Y))
#define max(X, Y)  ((X) > (Y) ? (X) : (Y))

//Square and cube of a value
#define sqr(X)  ((X)*(X))
#define cub(X)  ((X)*(X)*(X))

//Nearest distance to a different position given periodic boundary conditions
#define  NEAREST(x) (((x)>0.5*All.Lbox) ? ((x)-All.Lbox) : (((x)<-0.5*All.Lbox) ? ((x)+All.Lbox) : (x)))

//Memory manager Macro functions
#ifndef DISABLE_MEMORY_MANAGER
#define  emalloc(x, y)            malloc_fullinfo(x, y, __FUNCTION__, __FILE__, __LINE__)
#define  emalloc_movable(x, y, z) malloc_movable_fullinfo(x, y, z, __FUNCTION__, __FILE__, __LINE__)

#define  erealloc(x, y)           realloc_fullinfo(x, y, __FUNCTION__, __FILE__, __LINE__)
#define  erealloc_movable(x, y)   realloc_movable_fullinfo(x, y, __FUNCTION__, __FILE__, __LINE__)

#define  efree(x)                 free_fullinfo(x, __FUNCTION__, __FILE__, __LINE__)
#define  efree_movable(x)         free_movable_fullinfo(x, __FUNCTION__, __FILE__, __LINE__)

#define  report_memory_usage(x, y) report_detailed_memory_usage_of_largest_task(x, y, __FUNCTION__, __FILE__, __LINE__)

//Standard memory allocation if the Memory manager is disabled
#else
#define  emalloc(x, y)            malloc(y)
#define  emalloc_movable(x, y, z) malloc(z)

#define  erealloc(x, y)           realloc(x, y)
#define  erealloc_movable(x, y)   realloc(x, y)

#define  efree(x)                 free(x)
#define  efree_movable(x)         free(x)

#define  report_memory_usage(x, y) printf("Memory manager disabled.\n")
#endif
// \endcond


///////////////////////////////////////////////////////////////////////////////////////////////////
// Structures                                                                                    //
///////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Structure containing the same data for all tasks
///
/// This structure contains data which is the SAME for all tasks (mostly code parameters read from
/// the parameter file).  Holding this data in a structure is convenient for writing/reading the
/// restart file, and it allows the introduction of new global variables in a simple way. The only
/// thing to do is to introduce them into this structure.
///////////////////////////////////////////////////////////////////////////////////////////////////
extern struct global_data
{
	char treefile_name[NSTRING];         ///< Base name of the tree file
	char model_name[NSTRING];            ///< Name of the model
	char OutputDir[NSTRING];             ///< Output directory (will be generated if not existent)
	char output_redshifts[NSTRING];      ///< Comma separated string of output redshifts
	char output_mass_mb[NSTRING];        ///< Comma separated string of masses for main branch output
	char output_mass_mb_bin[NSTRING];    ///< Comma separated string of mass bin sizes for main branch output
	char fullline[NSTRING];              ///< String with #LINESIZE elements given by #SYMBOL
	char startline[NSTRING];             ///< String with 3 elements given by #SYMBOL

	int Mode;                            ///< Model mode stating what will be done (MCMC, etc)
	int verbose;                         ///< States how verbose the output is
	int NUniverses;                      ///< Number of universes that will be computed in parallel
	int NTaskPerUniverse;                ///< Number of tasks per universe
	int NumOutputFiles;                  ///< Number of output files for the galaxy catalogue
	int NumFilesInParallel;              ///< Number of files read/written in parallel
	int OutputFormat;                    ///< Determines which file format is used (1: ascii, 2:hdf5)
	int MaxTrees;                        ///< Maximum number of trees that can be held on one task
	int MaxHalos;                        ///< Maximum number of haloes that can be held on one task
	int MaxOrphans;                      ///< Maximum number of orphans that can be held on one task
	int MaxMemSize;                      ///< Maximum memory size allocated at the start in MB
	int BufferSize;                      ///< Size of the Communication Buffer in MB
	int NTimesteps;                      ///< Number of timesteps in the simulation
	int Nparam;                          ///< Number of free parameters that will be fitted
	int Nobs;                            ///< Number of observed statistics that will be read
	int MaxNLeaves;                      ///< Maximum number of leaves per tree
	int NMassFormed;                     ///< Number of entries per thread for MassFormed (MaxNLeaves*NTimesteps)
	int NStatistics;                     ///< Number of entries per thread for Statistics (Nmstar*NTimesteps)
	int Seed;                            ///< Initial seed for the random number generator
	int Nwalkers;                        ///< Number of walkers for probing parameter space
	int Noutputredshifts;                ///< Number of redshifts for the output catalogues
	int Noutputbranch;                   ///< Number of mass bins for main branch history output
	int MainBranchMassType;              ///< Mass type for main branch output (0: halo, 1: stellar)
	int MainBranch_iscale;               ///< Simulation snapshot numbers for main branch selection redshift

	unsigned long long TotNtrees;        ///< Total number of trees
	unsigned long long TotNhalos;        ///< Total number of haloes

	double h_100;                        ///< Hubble parameter in units of 100km/s/Mpc
	double Omega_0;                      ///< Matter density at redshift 0
	double Omega_Lambda_0;               ///< Dark energy density at redshift 0
	double Omega_Baryon_0;               ///< Baryon density at redshift 0
	double f_baryon;                     ///< Universal baryon fraction
	double Hubbletime;                   ///< Hubble time in internal time units
	double Lbox;                         ///< Side length of the simulation box in internal length units
	double obssigma0;                    ///< Normalisation at redshift 0 of the observational mass uncertainty
	double obssigmaz;                    ///< Redshift slope of the observational mass uncertainty
	double mcmca;                        ///< Scale parameter for the MCMC
	double temperature;                  ///< Starting temperature for the chains
	double timelimit;                    ///< Time limit in seconds after which the run is stopped
	double timenow;                      ///< Current time
	double starttime;                    ///< Time at the start of the run
	double minmass;                      ///< Minimum stellar mass of the galaxy catalogue
	double mainBranchRedshift;           ///< Redshift at which masses are selected for main branch output

	// Model parameters read from the parameter file
	double M0;                           ///< Characteristic mass of instantaneous conversion efficiency
	double Epsilon0;                     ///< Normalisation of instantaneous conversion efficiency
	double Beta0;                        ///< Low mass slope of instantaneous conversion efficiency
	double Gamma0;                       ///< High mass slope of instantaneous conversion efficiency
	double MZ;                           ///< Evolution of characteristic mass of instantaneous conversion efficiency
	double EpsilonZ;                     ///< Evolution of normalisation of instantaneous conversion efficiency
	double BetaZ;                        ///< Evolution of low mass slope of instantaneous conversion efficiency
	double GammaZ;                       ///< Evolution of high mass slope of instantaneous conversion efficiency
	double Fesc;                         ///< Fracion of stellar mass going to ICM in merger
	double Fstrip;                       ///< Fracion of peak halo mass for which stars get stripped
	double Tau0;                         ///< Fracion of dynamical time after peak halo mass when galaxy is quenched
	double TauS;                         ///< Slope for quenching as function of stellar mass
	double TauD;                         ///< Decay timescale for quenching as fracion of dynamical time

	// Range of model parameters read from the parameter file (prior)
	double DeltaM0;                      ///< Characteristic mass of instantaneous conversion efficiency (Range)
	double DeltaEpsilon0;                ///< Normalisation of instantaneous conversion efficiency (Range)
	double DeltaBeta0;                   ///< Low mass slope of instantaneous conversion efficiency (Range)
	double DeltaGamma0;                  ///< High mass slope of instantaneous conversion efficiency (Range)
	double DeltaMZ;                      ///< Evolution of characteristic mass of instantaneous conversion efficiency (Range)
	double DeltaEpsilonZ;                ///< Evolution of normalisation of instantaneous conversion efficiency (Range)
	double DeltaBetaZ;                   ///< Evolution of low mass slope of instantaneous conversion efficiency (Range)
	double DeltaGammaZ;                  ///< Evolution of high mass slope of instantaneous conversion efficiency (Range)
	double DeltaFesc;                    ///< Fracion of stellar mass going to ICM in merger (Range)
	double DeltaFstrip;                  ///< Fracion of peak halo mass for which stars get stripped (Range)
	double DeltaTau0;                    ///< Fracion of dynamical time after peak halo mass when galaxy is quenched (Range)
	double DeltaTauS;                    ///< Slope for quenching as function of stellar mass (Range)
	double DeltaTauD;                    ///< Decay timescale for quenching as fracion of dynamical time (Range)

	double x_unit;                       ///< Internal length unit
	double t_unit;                       ///< Internal time unit
	double m_unit;                       ///< Internal mass unit

	//Data
	char smffile_name[NSTRING];          ///< Name of the file containing the stellar mass functions
	char fqfile_name[NSTRING];           ///< Name of the file containing the quenched fractions
	char csfrdfile_name[NSTRING];        ///< Name of the file containing the cosmic star formation rate density
	char ssfrfile_name[NSTRING];         ///< Name of the file containing the specific star formation rates
	char wpfile_name[NSTRING];           ///< Name of the file containing the projected correlation functions
	int Nsmf;                            ///< Number of data points for the stellar mass functions
	int Nfq;                             ///< Number of data points for the quenched fractions
	int Ncsfrd;                          ///< Number of data points for the cosmic star formation rate density
	int Nssfr;                           ///< Number of data points for the specific star formation rates
	int Nwp;                             ///< Number of data points for the projected correlation functions
	int Nsmfset;                         ///< Number of data sets for the stellar mass functions
	int Nfqset;                          ///< Number of data sets for the quenched fractions
	int Ncsfrdset;                       ///< Number of data sets for the cosmic star formation rate density
	int Nssfrset;                        ///< Number of data sets for the specific star formation rates
	int Nwpset;                          ///< Number of mass bins for the projected correlation functions
	int Wpscale;                         ///< Simulation snapshot number for the projected correlation functions
	double wpredshift;                   ///< Global redshift of the projected correlation functions
	double wpmmin;                       ///< Minimum stellar mass for the projected correlation functions
	double wpmmax;                       ///< Maximum stellar mass for the projected correlation functions
	double GlobalSigmaSmfLz;             ///< Global error for the low redshift stellar mass functions
	double GlobalSigmaSmfHz;             ///< Global error for the high redshift stellar mass functions
	double GlobalSigmaFq;                ///< Global error for the quenched fractions
	double GlobalSigmaCsfrd;             ///< Global error for the cosmic star formation rate density
	double GlobalSigmaSsfr;              ///< Global error for the specific star formation rates
	double GlobalSigmaWp;                ///< Global error for the projected correlation functions

	//Statistics
	int Nmstar;                          ///< Number of stellar mass bins for the interpolation
	double mstarmin;                     ///< Minimum stellar mass for the interpolation
	double mstarmax;                     ///< Maximum stellar mass for the interpolation
	double dmstar;                       ///< Stellar mass bin width for the interpolation

} All;                                 ///< Stores all global variables that are the same on each task


///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Halo structure that holds all information on a halo
///
/// This structure contains all information for each dark matter halo and the galaxy at its centre.
/// The array H is allocated on each task so that it containts all haloes on a given task.
///////////////////////////////////////////////////////////////////////////////////////////////////
extern struct halo
{
	IDType haloid;              ///< ID of the halo
	IDType descid;              ///< ID of the halo's descendant
	IDType upid;                ///< ID of the most massive host halo (upid)

	int iscale;                 ///< Index of the scale factor array corresponding to this halo
	int imain;                  ///< Index of the main co-progenitor
	int idesc;                  ///< Index of the descandant
	int iprog;                  ///< Index of the main progenitor
	int impeak;                 ///< Index of the halo for which the mass reached its peak

	//Galaxy
	int icoprog;                ///< Index of the next co-progenitor
	int idescgal;               ///< Index of the descandant
	int imass;                  ///< Index of the 0th element of the mass array MassFormed
	int itdf;                   ///< Index of when dynamical friction time was last computed

	unsigned short np;          ///< Number of progenitors
	unsigned short mmp;         ///< Flag for most massive progenitor (1: yes, 0: no)
	unsigned short type;        ///< Halo type (0: main halo, 1: subhalo, 2: orphan)
	unsigned short gone;        ///< Flag indicating if this halo is present or gone

	float pos[3];               ///< 3D position
	float vel[3];               ///< 3D velocity
	float a;                    ///< Scale factor
	float mvir;                 ///< Virial mass
	float rvir;                 ///< Virial radius
	float c;                    ///< Concentration \f$(r_\mathrm{vir}/r_\mathrm{s})\f$
	float lambda;               ///< Spin parameter
	float mdotbary;             ///< Baryonic accretion rate averaged over one dynamical time

	//Galaxy
	float sfr;                  ///< Star formation rate
	float mstar;                ///< Stellar mass
	float icm;                  ///< Intra cluster mass
	float tdf;                  ///< Dynamical friction time

#ifdef COMPUTE_ICM
	IDType forestid;            ///< ID of the forest
	int ihost;                  ///< Index of the host halo (possibly on another task)
#endif

} *H;                         ///< Stores all haloes and their properties on a given task


///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Halo structure that is is used to read in the merger trees in binary form
///
/// For any merger tree to read the haloes MUST be provided in this structure. Merger trees
/// created with the consistent-trees code can be converted to this format with the
/// convert_CT_to_emerge.c code in the tools/ folder.
///////////////////////////////////////////////////////////////////////////////////////////////////
extern struct haloread
{
	IDType haloid;        ///< ID of the halo
	IDType descid;        ///< ID of the halo's descendant
	IDType upid;          ///< ID of the most massive host halo (upid)

	unsigned short np;    ///< Number of progenitors
	unsigned short mmp;   ///< Flag for most massive progenitor (1: yes, 0: no)

	float a;              ///< Scale factor
	float mvir;           ///< Virial mass
	float rvir;           ///< Virial radius
	float c;              ///< Concentration (rvir/rs)
	float lambda;         ///< Spin parameter

	float pos[3];         ///< 3D position
	float vel[3];         ///< 3D velocity

} *HR;                  ///< Stores all haloes and properties during reading


///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Parameter structure that holds all parameters
///
/// On each task all current parameter values are stored in the structure P. They can be modified
/// during a chain loop but need to the same for each universe.
///////////////////////////////////////////////////////////////////////////////////////////////////
extern struct parameters
{
	double M0;                  ///< Characteristic mass of instantaneous conversion efficiency
	double Epsilon0;            ///< Normalisation of instantaneous conversion efficiency
	double Beta0;               ///< Low mass slope of instantaneous conversion efficiency
	double Gamma0;              ///< High mass slope of instantaneous conversion efficiency
	double MZ;                  ///< Evolution of characteristic mass of instantaneous conversion efficiency
	double EpsilonZ;            ///< Evolution of normalisation of instantaneous conversion efficiency
	double BetaZ;               ///< Evolution of low mass slope of instantaneous conversion efficiency
	double GammaZ;              ///< Evolution of high mass slope of instantaneous conversion efficiency
	double Fesc;                ///< Fracion of stellar mass going to ICM in merger
	double Fstrip;              ///< Fracion of peak halo mass for which stars get stripped
	double Tau0;                ///< Fracion of dynamical time after peak halo mass when galaxy is quenched
	double TauS;                ///< Slope for quenching as function of stellar mass
	double TauD;                ///< Decay timescale for quenching as fracion of dynamical time

} P;                          ///< Stores all parameters on a given task (all the same within a universe)


///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Galaxy data structure that holds all information for each observational data point
///
/// For each observation the structure member can have a different meaning. For the stellar mass
/// function and the quenched fractions obs_x is the stellar mass and obs_bin is the scale factor.
/// For the specific star formation rate obs_x is the scale factor and obs_bin is the stellar mass.
///////////////////////////////////////////////////////////////////////////////////////////////////
struct galaxy_data
{
	int nfit;          ///< Weighting in the fit. Should only be different from 1 for testing.
	float obs_x;       ///< X-value of the observation (SMF/FQ: mstar, CSFRD/SSFR: a)
	float obs_y;       ///< Y-value of the observation
	float obs_sigma;   ///< Uncertainty of the observation
	float bin;         ///< Bin for which the observation has been measured (SMF/FQ: a, CSFRD/SSFR: mstar)
	float mod_y;       ///< Y-value predicted by the model
	float mod_sigma;   ///< Uncertainty predicted by the model
};


///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Data set structure that holds all information for each specific data set
///
/// A structure member contains all information for any given data set.
///
///////////////////////////////////////////////////////////////////////////////////////////////////
struct data_set
{
	float min;            ///< Minimum value of each bin (SMF, FQ: redshift; SSFR, WP: stellar mass)
	float max;            ///< Maximum value of each bin (SMF, FQ: redshift; SSFR, WP: stellar mass)
	float cut;            ///< WP: Maximum radius along the line of sight for the projection \f$\pi_\mathrm{max}\f$
	int ndata;            ///< Number of data points in each data set
	int offset;           ///< Offset for each data set (points to the first data point for each data set)
	char tag[SSTRING];    ///< Name tag to identify the data set
};


///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Simple structure to create a searchable ID list
///
/// Each element has an ID and an original position. The list can then be sorted according to the
/// ID and the original position be retrieved by the function binary_search_id.
///////////////////////////////////////////////////////////////////////////////////////////////////
struct search_list_IDs
{
	IDType ID;        ///< ID of the element
	int position;     ///< Original position of the element
};


///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Simple structure to sort the load
///
/// Each element stores the number of haloes and the load on its task as well as the task number.
/// The list can then be sorted according to the load and the task number be retrieved with
/// TaskSorted.
///////////////////////////////////////////////////////////////////////////////////////////////////
struct task_list
{
	int NhalosInTask;     ///< Number of haloes on this task
	double LoadInTask;    ///< Load in this task
	int TaskSorted;       ///< Task number
};


///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief kd-tree structure that is used to search for pairs to compute the correlation function
///
/// A kd-tree consists of node_count nodes, the dimension, the number of objects it contains.
/// Additionally the starting level containing nstartnodes nodes is included.
///////////////////////////////////////////////////////////////////////////////////////////////////
struct kd_tree
{
	struct kd_node *nodes;       ///< Array that stores all nodes of the tree
	int node_count;              ///< Number of nodes in the tree
	int count;                   ///< Number of objects in the tree
	int dim;                     ///< Dimension of the tree
	int lmin;                    ///< Minimum level for the parallel computation of pairs
	int lmax;                    ///< Maximum level for the parallel computation of pairs
	int midlevel;                ///< Level of the tree for which all nodes will be computed in parallel
	int nstartnodes;             ///< Number of nodes that will be computed in parallel
};


///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Node structure for the kd-tree that stores the children nodes the count and the bounding box
///
/// A kd-tree node can be either a leaf (leaf = 1) or a parent node (leaf = 0). If it is a leaf the
/// integers i0 and i1 store the indices of the first and the last object in the leaf. If it is a
/// parent node the integers i0 and i1 store the indices of the children nodes. The number of objects
/// in the node is stored in cnt. The bounding box that is defined by the minimum volume hyper-rectangle
/// around all objects in the node is spanned by the two positions min and max.
///////////////////////////////////////////////////////////////////////////////////////////////////
struct kd_node
{
	int i0;                 ///< Index pointing to the left child node (if not leaf) or to the first object in the node (if leaf).
	int i1;                 ///< Index pointing to the right child node (if not leaf) or to the last object in the node (if leaf).
	int leaf;               ///< Flag that defines if the node is a leaf (1) or a parent node (0).
	int cnt;                ///< Number of objects in the node.
	float min[TREEDIM];     ///< Minimum position of the bounding box in all dimensions.
	float max[TREEDIM];     ///< Maximum position of the bounding box in all dimensions.
};


///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Position structure for the kd-tree that stores the position of all objects in the tree
///
/// Each object has a position defined by a #TREEDIM-dimensional array x.
///////////////////////////////////////////////////////////////////////////////////////////////////
struct kd_pos
{
	float x[TREEDIM];       ///< Position array with dimension #TREEDIM.
};


///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Galaxy structure holding all galaxies that are used in the correlation functions
///
/// Each object has a position defined by a #TREEDIM-dimensional array x, a mass and a sfr.
///////////////////////////////////////////////////////////////////////////////////////////////////
struct galaxy
{
	float x[TREEDIM];       ///< Position array with dimension #TREEDIM.
	float mass;             ///< Stellar mass of the galaxy
	float sfr;              ///< Star formation rate of the galaxy
};


///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Helper structure to build the kd-tree
///
/// This structure is used to build the kd-tree and stores the index of the node and all needed
/// information about the objects that are contained in the current node.
///////////////////////////////////////////////////////////////////////////////////////////////////
struct kd_build_task
{
	int node;               ///< Index of the current node.
	int first;              ///< Index of the first object in the current node.
	int last;               ///< Index of the last object in the current node.
};


///////////////////////////////////////////////////////////////////////////////////////////////////
// Global variables                                                                              //
///////////////////////////////////////////////////////////////////////////////////////////////////

extern int Ntrees;                 ///< Number of trees on this task
extern int Nhalos;                 ///< Number of haloes on this task
extern int Nforests;               ///< Number of forests on this task
extern int *NhalosInTree;          ///< Number of haloes in each tree on this task
extern int *NhalosInForest;        ///< Number of haloes in each forest on this task
extern int *NtreesInForest;        ///< Number of trees in each forest on this task
extern int *OffsetHalos;           ///< Offset of the first halo in a tree/forest
extern int Norphans;               ///< Number of orphans on this task
extern IDType *TreeIDInTree;       ///< ID of each tree on this task

extern int NThread;                ///< OMP number of threads
extern int NTask;                  ///< MPI number of tasks
extern int ThisTask;               ///< MPI number of this tasks
extern int MasterTask;             ///< MPI number of master tasks for each universe

extern int EndFlag;                ///< Flag for the endrun Macro

extern float *ScaleFactor;         ///< List of all scale factors in the merger trees
extern float *CosmicTime;          ///< Cosmic time corresponding to the above scale factors
extern float *Timestep;            ///< Timesteps corresponding the above scale factors
extern float *DynTime;             ///< Dynamical halo times corresponding the above scale factors

extern float *PopAge;              ///< Age of a stellar population
extern float *MassLeft;            ///< Stores the fraction of the mass left between each scale factor
extern float *MassFormed;          ///< Stores the mass of all populations ever formed in a branch
extern float *ICMFormed;           ///< Stores the mass of all populations in the ICM

extern int   *Output_iscale;       ///< Simulation snapshot numbers for the output catalogue
extern float *OutputRedshifts;     ///< Output redshifts of the galaxy catalogue

//Observational data
extern struct galaxy_data *Smf;    ///< Stellar Mass Functions
extern struct galaxy_data *Fq;     ///< Quenched Fractions
extern struct galaxy_data *Csfrd;  ///< Cosmic star formation rate density
extern struct galaxy_data *Ssfr;   ///< Specific star formation rates
extern struct galaxy_data *Wp;     ///< Projected correlation functions
extern struct data_set *SmfSet;    ///< Information about each data set of the Stellar Mass Functions
extern struct data_set *FqSet;     ///< Information about each data set of the Quenched Fractions
extern struct data_set *CsfrdSet;  ///< Information about each data set of the Cosmic star formation rate density
extern struct data_set *SsfrSet;   ///< Information about each data set of the Specific star formation rates
extern struct data_set *WpSet;     ///< Information about each data set of the projected correlation function

//Model statistics
extern float *Mstar;               ///< Array of stellar masses used for model statistics
extern float *Modelsmf;            ///< Model stellar mass function (NTimestep * Nmstar)
extern float *Modelsmfred;         ///< Model stellar mass function for red galaxies (NTimestep * Nmstar)
extern float *Modelcsfrd;          ///< Model cosmic SFR density (NTimestep)
extern float *Modelssfr;           ///< Model specific SFRs (NTimestep * Nmstar)
extern float *Radius;              ///< Array of radii used for the correlation function
extern float *Modelxi;             ///< Model 3d auto-correlation function
extern float *Modelxierr;          ///< Model 3d auto-correlation function error
extern struct kd_pos *Galpos;      ///< Positions of all galaxies in a stellar mass bin of the correlation function

//Fitting
extern int DoF;                    ///< Total degrees of freedom
extern int DoFsmf;                 ///< Degrees of freedom for the SMF
extern int DoFfq;                  ///< Degrees of freedom for the FQ
extern int DoFcsfrd;               ///< Degrees of freedom for the CSFRD
extern int DoFssfr;                ///< Degrees of freedom for the SSFR
extern int DoFwp;                  ///< Degrees of freedom for the WP

//Memory
extern size_t AllocatedBytes;      ///< Number of allocated Bytes
extern size_t HighMarkBytes;       ///< Highest mark of allocated Bytes
extern size_t FreeBytes;           ///< Number of free Bytes
extern size_t HighMark;            ///< High Mark

extern gsl_rng *rng_gaussian;      ///< Gaussian random number generator
extern gsl_rng *rng_uniform;       ///< Uniform random number generator

extern void *CommBuffer;           ///< Points to communication buffer, which is used at a few places

extern FILE *FpInfo,               ///< File pointer for info.txt log-file.
            *FpMemory;             ///< File pointer for memory.txt log-file.

#if (RANDOM_NUMBER_TABLE > 0)
extern float *RndTableGaussian;    ///< Random number table filled with gaussian random numbers
extern float *RndTableUniform;     ///< Random number table filled with uniform random numbers
#endif

#ifdef WRITE_MAINBRANCH
extern float *MainBranchMasses;    ///< Masses for the main branch history outputs
extern float *MainBranchBinSize;   ///< Mass bin sizes for the main branch history outputs
#endif


///////////////////////////////////////////////////////////////////////////////////////////////////
// MPI Tags                                                                                      //
///////////////////////////////////////////////////////////////////////////////////////////////////

#define TAG_NTREES          10      ///< MPI tag for the number of trees
#define TAG_NHALOS          11      ///< MPI tag for the number of haloes
#define TAG_NFORESTS        12      ///< MPI tag for the number of forests
#define TAG_NORPHANS        13      ///< MPI tag for the number of orphans
#define TAG_TREEID          14      ///< MPI tag for the tree ID
#define TAG_NHALOSINTREE    15      ///< MPI tag for the number of haloes in a tree
#define TAG_NTREESINFOREST  16      ///< MPI tag for the number of trees in a forest
#define TAG_NHALOSINFOREST  17      ///< MPI tag for the number of haloes in a forest
#define TAG_OFFSETHALOS     18      ///< MPI tag for the offset of haloes in a tree or forest
#define TAG_N               19      ///< MPI tag for the memory allocation
#define TAG_HDATA           20      ///< MPI tag for the halo array
#define TAG_IDS             21      ///< MPI tag for the IDs
#define TAG_SMF             22      ///< MPI tag for the stellar mass function
#define TAG_SMFRED          23      ///< MPI tag for the stellar mass function of quenched galaxies
#define TAG_CSFRD           24      ///< MPI tag for the cosmic star formation rate density
#define TAG_SSFR            25      ///< MPI tag for the specific star formation rate
#define TAG_PROB            26      ///< MPI tag for the probability
#define TAG_REQUEST         27      ///< MPI tag for the pair count requests
#define TAG_COUNT           28      ///< MPI tag for the pair count addition


///////////////////////////////////////////////////////////////////////////////////////////////////
// Verbose Tags                                                                                  //
///////////////////////////////////////////////////////////////////////////////////////////////////
#define VERBOSE_MIN 1               ///< Tag for mimimal screen output
#define VERBOSE_ALL 2               ///< Tag for all screen output

#endif
