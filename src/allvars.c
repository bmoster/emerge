///////////////////////////////////////////////////////////////////////////////////////////////////
// Emerge code - File allvars.c                                                                  //
///////////////////////////////////////////////////////////////////////////////////////////////////
///
/// \file allvars.c
/// \brief Contains global variables and declarations
///
/// This file declares all global variables. Further variables should be added here and in the
/// file allvarsh.h where they should be declared as 'extern'.
///
///////////////////////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "allvars.h"


///////////////////////////////////////////////////////////////////////////////////////////////////
// Structures                                                                                    //
///////////////////////////////////////////////////////////////////////////////////////////////////

struct global_data All;  ///< Stores all global variables that are the same on each task
struct halo *H;          ///< Stores all haloes and their properties on a given task
struct haloread *HR;     ///< Stores all haloes and properties during reading
struct parameters P;     ///< Stores all parameters on a given task (all the same within a universe)


///////////////////////////////////////////////////////////////////////////////////////////////////
// Global variables                                                                              //
///////////////////////////////////////////////////////////////////////////////////////////////////

int Ntrees;                 // Number of trees on this task
int Nhalos;                 // Number of haloes on this task
int Nforests;               // Number of forests on this task
int *NhalosInTree;          // Number of haloes in each tree on this task
int *NhalosInForest;        // Number of haloes in each forest on this task
int *NtreesInForest;        // Number of trees in each forest on this task
int *OffsetHalos;           // Offset of the first halo in a tree/forest
int Norphans;               // Number of orphans on this task
IDType *TreeIDInTree;       // ID of each tree on this task

int NThread;                // OMP number of threads
int NTask;                  // MPI number of tasks
int ThisTask;               // MPI number of this tasks
int MasterTask;             // MPI number of master tasks for each universe

float *ScaleFactor;         // List of all scale factors in the merger trees
float *CosmicTime;          // Cosmic time corresponding to the above scale factors
float *Timestep;            // Timesteps corresponding the above scale factors
float *DynTime;             // Dynamical halo times corresponding the above scale factors

float *PopAge;              // Age of a stellar population
float *MassLeft;            // Stores the fraction of the mass left between each scale factor
float *MassFormed;          // Stores the mass of all populations ever formed in a branch
float *ICMFormed;           // Stores the mass of all populations in the ICM

int   *Output_iscale;       // Simulation snapshot numbers for the output catalogue
float *OutputRedshifts;     // Output redshifts of the galaxy catalogue

//Observational data
struct galaxy_data *Smf;    // Stellar Mass Functions
struct galaxy_data *Fq;     // Quenched Fractions
struct galaxy_data *Csfrd;  // Cosmic star formation rate density
struct galaxy_data *Ssfr;   // Specific star formation rates
struct galaxy_data *Wp;     // Projected correlation functions
struct data_set *SmfSet;    // Information about each data set of the Stellar Mass Functions
struct data_set *FqSet;     // Information about each data set of the Quenched Fractions
struct data_set *CsfrdSet;  // Information about each data set of the Cosmic star formation rate density
struct data_set *SsfrSet;   // Information about each data set of the Specific star formation rates
struct data_set *WpSet;     // Information about each bin of the projected correlation function

//Model statistics
float *Mstar;               // Array of stellar masses used for model statistics
float *Modelsmf;            // Model stellar mass function (NTimestep * Nmstar);
float *Modelsmfred;         // Model stellar mass function for red galaxies (NTimestep * Nmstar);
float *Modelcsfrd;          // Model cosmic SFR density (NTimestep)
float *Modelssfr;           // Model specific SFRs (NTimestep * Nmstar)
float *Radius;              // Array of radii used for the correlation function
float *Modelxi;             // Model 3d auto-correlation function
float *Modelxierr;          // Model 3d auto-correlation function error
struct kd_pos *Galpos;      // Positions of all galaxies in a stellar mass bin of the correlation function

//Fitting
int DoF;                    // Total degrees of freedom
int DoFsmf;                 // Degrees of freedom for the SMF
int DoFfq;                  // Degrees of freedom for the FQ
int DoFcsfrd;               // Degrees of freedom for the CSFRD
int DoFssfr;                // Degrees of freedom for the SSFR
int DoFwp;                  // Degrees of freedom for the WP

//Memory
int EndFlag;                // Flag for the endrun Macro
size_t AllocatedBytes;      // Number of allocated Bytes
size_t HighMarkBytes;       // Highest mark of allocated Bytes
size_t FreeBytes;           // Number of free Bytes
size_t HighMark;            // High Mark

gsl_rng *rng_gaussian;      // Gaussian random number generator
gsl_rng *rng_uniform;       // Uniform random number generator

void *CommBuffer;           // Points to communication buffer, which is used at a few places */

FILE *FpInfo,               // File pointer for info.txt log-file.
     *FpMemory;             // File pointer for memory.txt log-file.

#if (RANDOM_NUMBER_TABLE > 0)
float *RndTableGaussian;    // Random number table filled with gaussian random numbers
float *RndTableUniform;     // Random number table filled with uniform random numbers
#endif

#ifdef WRITE_MAINBRANCH
float *MainBranchMasses;    // Masses for the main branch history outputs
float *MainBranchBinSize;   // Mass bin sizes for the main branch history outputs
#endif
