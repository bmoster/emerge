///////////////////////////////////////////////////////////////////////////////////////////////////
// Emerge code - File galaxies.c                                                                 //
///////////////////////////////////////////////////////////////////////////////////////////////////
///
/// \file galaxies.c
/// \brief Contains functions that populate haloes in a merger tree with galaxies
///
/// This file contains all functions that are needed to compute galaxy related quantities.
/// This includes a function to initialise all arrays, a function that generates all galaxies
/// by calling a function that resets all galaxy quantities followed by calling a function that
/// assembles all galaxies, and a function that finalises the galaxy computation by freeing all
/// memory. Also this file contains some functions that are needed for the galaxy computation,
/// such as functions that compute the current stellar mass and ICM by summing over the stellar
/// populations formed at each scale factor, and a function to compute the dynamical friction time.
///
///////////////////////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_rng.h>
#include "allvars.h"
#include "proto.h"

#ifdef HDF5_SUPPORT
#include <hdf5.h>
#endif


/*! \brief This function allocates and initialises all arrays needed for galaxy formation
 *
 *  This function should be called once before the make_galaxies is called to set up all arrays
 *  that are needed. This includes the array #MassFormed that stores the stellar populations in the galaxy
 *  that formed in each branch and the array #ICMFormed that stores the same for the ICM. Also the imass
 *  index is set which links each halo to a mass array.
 */
void init_galaxies(void)
{
	int ihalo,imass;

#ifndef COMPUTE_ICM
	int itree;
#else
	int iforest;
#endif

	//First find the maximum number of leaves
#ifndef COMPUTE_ICM
	//Loop through all trees
	for (itree = 0; itree < Ntrees; itree++)
	{
		//Go through all haloes in this tree
		for (ihalo = OffsetHalos[itree], imass = 0; ihalo < OffsetHalos[itree] + NhalosInTree[itree]; ihalo++)
		{
#else
	//Loop through all forests
	for (iforest = 0; iforest < Nforests; iforest++)
	{
		//Go through all haloes in this tree
		for (ihalo = OffsetHalos[iforest], imass = 0; ihalo < OffsetHalos[iforest] + NhalosInForest[iforest]; ihalo++)
		{
#endif
			//If this halo does not have any progenitors
			if (H[ihalo].np < 1)
			{
				//Set the branch number for this halo (equal to the number of its leaf)
				H[ihalo].imass = imass;
				//Increment the number of leaves
				imass++;
			}
			//If the halo has a descendant set its imass to this halo's imass
			if (H[ihalo].idesc >= 0) H[H[ihalo].idesc].imass = H[ihalo].imass;
		}
		//If imass in this tree/forest is larger than the global All.MaxNLeaves increase the global one
		if (imass > All.MaxNLeaves) All.MaxNLeaves = imass;
	}

	//For formating purposes add the two braces
#ifdef NEVERDEF
}
}
#endif

	//Get the maximum number of haloes per scale factor in any tree throughout all tasks
	MPI_Allreduce(MPI_IN_PLACE, &All.MaxNLeaves, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

	//Allocate the array MassFormed which stores the mass of all populations in a branch
	MassFormed = emalloc_movable(&MassFormed,"MassFormed", All.MaxNLeaves * All.NTimesteps * NThread * sizeof(float));
#ifdef COMPUTE_ICM
	ICMFormed  = emalloc_movable(&ICMFormed, "ICMFormed",  All.MaxNLeaves * All.NTimesteps * NThread * sizeof(float));
#endif

	//Allocate arrays for statistics
	Modelsmf    = emalloc_movable(&Modelsmf,    "ModelSMF",     All.NTimesteps * All.Nmstar * NThread * sizeof(float));
	Modelsmfred = emalloc_movable(&Modelsmfred, "ModelSMFRed",  All.NTimesteps * All.Nmstar * NThread * sizeof(float));
	Modelcsfrd  = emalloc_movable(&Modelcsfrd,  "ModelCSFRD",   All.NTimesteps * NThread * sizeof(float));
	Modelssfr   = emalloc_movable(&Modelssfr,   "ModelSSFR",    All.NTimesteps * All.Nmstar * NThread * sizeof(float));

#ifdef READ_WP
	//Allocate arrays for the correlation functions
	Radius      = emalloc_movable(&Radius,     "RADIUS",     All.Nwpset * WP_RBINS_INT * sizeof(float));
	Modelxi     = emalloc_movable(&Modelxi,    "MODELXI",    All.Nwpset * WP_RBINS_INT * sizeof(float));
	Modelxierr  = emalloc_movable(&Modelxierr, "MODELXIERR", All.Nwpset * WP_RBINS_INT * sizeof(float));
#endif

	//Set number of entries per thread for MassFormed
	All.NMassFormed = All.MaxNLeaves * All.NTimesteps;
	//Set number of entries per thread for Statistics
	All.NStatistics = All.Nmstar * All.NTimesteps;

	//Print the memory usage
	HighMark = 0;
	report_memory_usage(&HighMark, "Init_Galaxies");

}


/*! \brief This function frees all arrays needed for galaxy formation
 *
 *  This function should be called once after make_galaxies has been called (as often as needed) to free
 *  all arrays that were needed.
 */
void finish_galaxies(void)
{
#ifdef READ_WP
	efree(Modelxierr);
	efree(Modelxi);
	efree(Radius);
#endif
	efree(Modelssfr);
	efree(Modelcsfrd);
	efree(Modelsmfred);
	efree(Modelsmf);
#ifdef COMPUTE_ICM
	efree_movable(ICMFormed);
#endif
	efree_movable(MassFormed);

	//And print the memory usage
	report_memory_usage(&HighMark, "Finish_Galaxies");

}


/*! \brief This function creates galaxies in the dark matter halo trees
 *
 *  This function should be called each time the haloes shall be populated with galaxies. It needs to have
 *  #init_galaxies be called once at some point before.
 */
void make_galaxies(void)
{
	//Print what is done...
	if (ThisTask == 0 && All.verbose >= VERBOSE_ALL)
		printf("%s\n%s Populating trees with galaxies...\n", All.fullline,All.startline);

	reset_galaxies();
	assemble_galaxies();

	//Print that we are done here...
	if (ThisTask == 0 && All.verbose >= VERBOSE_ALL)
		printf("%s ... done populating the trees with galaxies.\n",All.startline);

}


/*! \brief This function resets all galaxy related variables
 *
 *  This function should be called before each call of #assemble_galaxies. It resets all variables to
 *  the default values without containing any stars. All mmp are set to 1 and the coprog indices are set
 *  to -1 so that there are no mergers (yet). Also all haloes are potentially still present so gone is set
 *  to zero. All SFRs and stellar masses are set to zero and the dynamical friction time is set to the current
 *  time. All orphan virial masses are set to the mass at the time they became an orphan.
 *  The mass arrays that hold the stellar masses formed in each branch are set to zero. Finally the seed for
 *  the gaussian random number generator is reset if the random number table is not used.
 */
void reset_galaxies(void)
{
	int i;

	//Reset all galaxy properties in the haloes to zero and the dynamical friction time to the current cosmic time
#if defined(_OPENMP) && (OPENMPTHREADS > 1)
#pragma omp parallel for schedule(static)
#endif
	for (i = 0; i < Nhalos; i++)
	{
		H[i].mmp      =  1;
		H[i].icoprog  = -1;
		H[i].itdf     = -1;
		H[i].idescgal = H[i].idesc;
		H[i].gone     = 0;

		if (H[i].np > 1) H[i].np = 1;

		H[i].sfr      = 0.0;
		H[i].mstar    = 0.0;
		H[i].icm      = 0.0;
		H[i].tdf      = CosmicTime[H[i].iscale];

		if (H[i].type == 2) H[i].mvir = H[H[i].iprog].mvir;
	}

	//Reset the arrays for the stellar populations formed at any time in the galaxy and ICM to zero
	for (i = 0; i < All.MaxNLeaves * All.NTimesteps * NThread; i++)
	{
		MassFormed[i] = 0.0;
#ifdef COMPUTE_ICM
		ICMFormed[i]  = 0.0;
#endif
	}

#if (RANDOM_NUMBER_TABLE == 0)
	// start-up seed for gaussian random numbers
	gsl_rng_set(rng_gaussian, 140314081082 + (long)(ThisTask % All.NTaskPerUniverse));
#endif

	//Reset the model statistics to zero
	for (i = 0; i < All.NTimesteps * All.Nmstar * NThread; i++)
	{
		Modelsmf[i]    = 0.0;
		Modelsmfred[i] = 0.0;
		Modelssfr[i]   = 0.0;
	}
	for (i = 0; i < All.NTimesteps * NThread; i++) Modelcsfrd[i]  = 0.0;

}


/*! \brief This function creates a galaxy in each dark matter halo
 *
 *  This function goes through all trees/forests and haloes and forms a galaxy. It is a wrapper around
 *  assemble_tree, which populates haloes in individual trees/forests with galaxies. Depending on whether
 *  #COMPUTE_ICM is set, the function loops over trees and calls assemble_tree with itrees or loops over
 *  forests and calls assemble_tree with iforest.
 */
void assemble_galaxies(void)
{
#ifndef COMPUTE_ICM
	int itree;
#else
	int iforest;
#endif
	int thisthread = 0;

	//If we do not need the ICM we can go tree-by-tree
#ifndef COMPUTE_ICM
	//Loop through all trees
#if defined(_OPENMP) && (OPENMPTHREADS > 1)
#pragma omp parallel private(thisthread)
	{
		thisthread = omp_get_thread_num();
#pragma omp for schedule(dynamic,100)
#endif
	for (itree = 0; itree < Ntrees; itree++) assemble_tree(itree,thisthread);
#if defined(_OPENMP) && (OPENMPTHREADS > 1)
}
#endif
	//If we do need the ICM we must go forest-by-forest
#else
	//Loop through all forests
#if defined(_OPENMP) && (OPENMPTHREADS > 1)
#pragma omp parallel private(thisthread)
	{
		thisthread = omp_get_thread_num();
#pragma omp for schedule(dynamic,100)
#endif
	for (iforest = 0; iforest < Nforests; iforest++) assemble_tree(iforest,thisthread);
#if defined(_OPENMP) && (OPENMPTHREADS > 1)
}
#endif
#endif

}


/*! \brief This function populates haloes in a tree or forest with galaxies
 *
 *  This function goes through a specific tree/forest and haloes and forms galaxies. It is first checked if the halo
 *  is still present or if it merged or got stripped already. Then the dynamical friction time is set or
 *  possibly reset for orphans. The SFR of each halo is calculated as the product of the efficiency and
 *  the baryonic accretion rate. If the halo mass has fallen below its peak mass the SFR is kept constant.
 *  If the quenching time (tau * t_decay) has elapsed the SFR is set to zero. Given the SFR the stellar mass
 *  that has formed in this timestep in this branch is added to the #MassFormed array at the correct position.
 *  If a halo becomes a subhalo all its ICM is moved to the host halo. If the virial mass has dropped below
 *  a fraction of its peak mass (Fstrip * mpeak) the stellar mass is moved to the ICM of the host halo and
 *  the halo is marked as gone. If the dynamical friction time of an orphan has elapsed we merge it back to
 *  its main branch and add all stellar mass and ICM to the main galaxy (type 0 or 1 halo). The orphan is
 *  then marked as gone. Finally the stellar mass and ICM at each timestep is computed from the arrays
 *  #MassFormed and #ICMFormed given the stellar loss rate that has been stored in a table. At the end of each
 *  tree/forest the #MassFormed and #ICMFormed arrays are reset to zero.
 *
 *  \param itree Index of the tree or forest that is populated with galaxies
 *  \param thisthread Number of the OpenMP thread that calculates this tree
 */
void assemble_tree(int itree, int thisthread)
{
	int i, j, ihalo, iprog, imass, itdf, itdfmain, hstart, hend;
	float t_now, t_since_peak, t_decay, t_quench, radius, theta, phi, dx, dy, dz;

#ifndef COMPUTE_ICM
	hstart = OffsetHalos[itree];
	hend   = OffsetHalos[itree] + NhalosInTree[itree];
#else
	hstart = OffsetHalos[itree];
	hend   = OffsetHalos[itree] + NhalosInForest[itree];
#endif

	//Go through all haloes in this tree/forest
	for (ihalo = hstart, imass = 0; ihalo < hend; ihalo++)
	{

//////////////////////////////////////////////////////////////////////////////////////////////////////
//    Check if galaxy in halo still exists - otherwise skip                                         //
//////////////////////////////////////////////////////////////////////////////////////////////////////
		//Check if the halo is still there - otherwise mark its descendant as gone and skip
		if (H[ihalo].gone)
		{
			if (H[ihalo].idesc >= 0) H[H[ihalo].idesc].gone = 1;
			continue;
		}

//////////////////////////////////////////////////////////////////////////////////////////////////////
//    Compute the current cosmic time                                                               //
//////////////////////////////////////////////////////////////////////////////////////////////////////
		t_now = CosmicTime[H[ihalo].iscale];

//////////////////////////////////////////////////////////////////////////////////////////////////////
//    Check if halo is a leaf                                                                       //
//////////////////////////////////////////////////////////////////////////////////////////////////////
		//If this is a leaf increase the leaf index and set iprog to 0 (in case we ever need to call it)
		if (H[ihalo].np == 0)
		{
			imass++;
			iprog = 0;
		}
		//If the halo has at least one progenitor
		else
		{
			//Set its actual iprog index
			iprog = H[ihalo].iprog;
//////////////////////////////////////////////////////////////////////////////////////////////////////
//    Check if halo is an orphan and possibly reset the dynamical friction time                     //
//////////////////////////////////////////////////////////////////////////////////////////////////////
			if (H[ihalo].type == 2)
			{
				//If the halo has become an orphan this timestep
				//Set its dynamical friction time using iprog as satellite and H[imain].iprog as central
				if (H[iprog].type != 2)
				{
					H[ihalo].tdf = CosmicTime[H[iprog].iscale] + get_tdf(iprog,H[H[ihalo].imain].iprog);
					H[ihalo].itdf = iprog;
				}
				//If the halo is not the first orphan - check if merger clock needs reset
				//If the progenitor's main's descendant is an orphan - if so reset with H[iprog].imain as sat and H[imain].iprog as cen
				else if (H[H[H[iprog].imain].idesc].type == 2)
				{
					//H[ihalo].tdf = CosmicTime[H[iprog].iscale] + get_tdf(H[iprog].imain,H[H[ihalo].imain].iprog);
					H[ihalo].tdf = CosmicTime[H[iprog].iscale] + get_tdf(iprog,H[H[ihalo].imain].iprog);
					H[ihalo].itdf = iprog;
				} //End recomputation of tdf
			} //End Orphan
		} //End leaf/no leaf

//////////////////////////////////////////////////////////////////////////////////////////////////////
//    Check if halo is an orphan, compute position and inherit the dynamical friction time          //
//////////////////////////////////////////////////////////////////////////////////////////////////////
		//If the halo is an orphan set its descandants dynamical friction time
		if (H[ihalo].type == 2)
		{ //Get index of progenitors when tdf was computed
			itdf     = H[ihalo].itdf;
			itdfmain = H[H[H[itdf].idesc].imain].iprog;
#ifdef ORPHAN_MASSLOSS
			//Let the halo mass of the orphan decay on a timescale tdf-tpeak
			if (H[itdf].iscale > H[H[ihalo].impeak].iscale)
			{
				H[ihalo].mvir = H[iprog].mvir * pow(H[H[ihalo].impeak].mvir/H[itdf].mvir, -Timestep[H[ihalo].iscale]/(CosmicTime[H[itdf].iscale] - CosmicTime[H[H[ihalo].impeak].iscale]));
			}
			else H[ihalo].mvir = H[iprog].mvir;
#endif
			//Compute the initial radius scaled with respect to the dynamical friction time
			dx = NEAREST(H[itdf].pos[0] - H[itdfmain].pos[0]);
			dy = NEAREST(H[itdf].pos[1] - H[itdfmain].pos[1]);
			dz = NEAREST(H[itdf].pos[2] - H[itdfmain].pos[2]);
			radius = sqrt(dx*dx+dy*dy+dz*dz) * sqrt((H[ihalo].tdf-t_now)/(H[ihalo].tdf-CosmicTime[H[itdf].iscale]));
#ifdef ORPHAN_NONRADIAL_INFALL
			theta = get_uniform_random_number(ihalo) * M_PI;
			phi   = get_uniform_random_number(ihalo + RANDOM_NUMBER_TABLE/2) * 2.0 * M_PI;
#else
			theta = get_uniform_random_number(itdf) * M_PI;
			phi   = get_uniform_random_number(itdf + RANDOM_NUMBER_TABLE/2) * 2.0 * M_PI;
#endif
			//Position the orphan on a sphere with r=radius around imain
			H[ihalo].pos[0]   = H[H[ihalo].imain].pos[0] + radius * cos(phi) * sin(theta);
			H[ihalo].pos[1]   = H[H[ihalo].imain].pos[1] + radius * sin(phi) * sin(theta);
			H[ihalo].pos[2]   = H[H[ihalo].imain].pos[2] + radius * cos(theta);
			if (H[ihalo].pos[0] >= All.Lbox) H[ihalo].pos[0] -= All.Lbox;
			if (H[ihalo].pos[1] >= All.Lbox) H[ihalo].pos[1] -= All.Lbox;
			if (H[ihalo].pos[2] >= All.Lbox) H[ihalo].pos[2] -= All.Lbox;
			if (H[ihalo].pos[0] < 0.0) H[ihalo].pos[0] += All.Lbox;
			if (H[ihalo].pos[1] < 0.0) H[ihalo].pos[1] += All.Lbox;
			if (H[ihalo].pos[2] < 0.0) H[ihalo].pos[2] += All.Lbox;

			//If there is a descendent let it inherit the dynamical friction time and the index when it was computed
			if (H[ihalo].idesc >= 0)
			{
				H[H[ihalo].idesc].tdf  = H[ihalo].tdf;
				H[H[ihalo].idesc].itdf = H[ihalo].itdf;
			}
		}

//////////////////////////////////////////////////////////////////////////////////////////////////////
//    Compute relevant time scales                                                                  //
//////////////////////////////////////////////////////////////////////////////////////////////////////
		// Compute time that has elapsed since peak mass was reached
		t_since_peak = t_now - CosmicTime[H[H[ihalo].impeak].iscale];
		// Compute timescale on which sfr declines, here: dynamical time of the halo when peak mass was reached
		t_decay = DynTime[H[H[ihalo].impeak].iscale];

#ifdef SAT_QUENCH_MASS_DEPENDENT
		// Get the slope for the quenching timescale
		if (H[H[ihalo].impeak].mstar > 0) t_quench = t_decay * P.Tau0 * pow(H[H[ihalo].impeak].mstar*All.m_unit/1.0e10,-P.TauS);
		else t_quench = 2.0 * t_since_peak;
		if (t_quench < t_decay * P.Tau0) t_quench = t_decay * P.Tau0;
#else
		t_quench = P.Tau0 * t_decay;
#endif

//////////////////////////////////////////////////////////////////////////////////////////////////////
//    Compute star formation rate from efficiency and growth rate                                   //
//////////////////////////////////////////////////////////////////////////////////////////////////////
		//Set the SFR equal to the efficiency times the baryonic accretion rate
		H[ihalo].sfr = sfe(H[ihalo].mvir,H[ihalo].a) * H[ihalo].mdotbary;

//////////////////////////////////////////////////////////////////////////////////////////////////////
//    Satellite keeps forming stars at constant rate                                                //
//////////////////////////////////////////////////////////////////////////////////////////////////////
		// If Halo does not gain mass (passed peak mass) keep SFR constant
		if (H[ihalo].mvir < H[H[ihalo].impeak].mvir) H[ihalo].sfr = H[iprog].sfr;

//////////////////////////////////////////////////////////////////////////////////////////////////////
//    Satellite keeps forming stars at declining rate                                               //
//////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef SAT_SFR_EXP_DECAY
		// If there is an exponential decline in the sfr after peak mass reduce current sfr
		if (H[ihalo].mvir < H[H[ihalo].impeak].mvir) H[ihalo].sfr = H[H[ihalo].impeak].sfr * exp(-P.TauD * t_since_peak / t_decay);
#endif

//////////////////////////////////////////////////////////////////////////////////////////////////////
//    Satellite quenching                                                                           //
//////////////////////////////////////////////////////////////////////////////////////////////////////
		// If certain amount of time after tpeak has elapsed, quench galaxy
		if (t_since_peak > t_quench) H[ihalo].sfr = 0.0;

//////////////////////////////////////////////////////////////////////////////////////////////////////
//    Form stellar population for this scale factor and store in MassFormed array                   //
//////////////////////////////////////////////////////////////////////////////////////////////////////
		//Add the mass that has formed this timestep (sfr*dt) to the mass array
		MassFormed[H[ihalo].imass*All.NTimesteps+H[ihalo].iscale+thisthread*All.NMassFormed]
		  += H[ihalo].sfr * Timestep[H[ihalo].iscale];

//////////////////////////////////////////////////////////////////////////////////////////////////////
//    Add ICM of halo that becomes a subhalo to the main halo                                       //
//////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef COMPUTE_ICM
		// Add ICM from satellites in a host halo to the host's ICM if a galaxy becomes a satellite
		if (H[ihalo].type != 0 && H[iprog].type == 0 && H[ihalo].np > 0)
		{
			for (j = 0; j <= H[ihalo].iscale; j++)
			{
				//Add the ICM of the satellite to the main halo
				ICMFormed[H[H[ihalo].ihost].imass*All.NTimesteps+j+thisthread*All.NMassFormed]
				  += ICMFormed[H[ihalo].imass*All.NTimesteps+j+thisthread*All.NMassFormed];
				//Set the ICM array to zero so there is no ICM left in descendant subhaloes
				ICMFormed[H[ihalo].imass*All.NTimesteps+j+thisthread*All.NMassFormed] = 0.0;
			}
		}
#endif

//////////////////////////////////////////////////////////////////////////////////////////////////////
//    Satellite stripping                                                                           //
//////////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef SAT_STRIP_USE_MSTAR
		// STRIPPING: If halo has lost mass down to some fraction fstrip * mpeak assume that the stellar mass is stripped and goes to the hosts icm
		if (H[ihalo].mvir < P.Fstrip * H[H[ihalo].impeak].mvir)
#else
		//Alternative: If halo has lost mass down to some fraction fstrip * mstar assume that the stellar mass is stripped and goes to the hosts icm
		if (H[ihalo].mvir < P.Fstrip * H[iprog].mstar && H[ihalo].np > 0)
#endif
		{
#ifdef COMPUTE_ICM
			//If the host is on this task
			for (j = 0; j <= H[ihalo].iscale; j++)
			{
				//Add the stellar mass of the satellite to the main halo's ICM
				ICMFormed[H[H[ihalo].ihost].imass*All.NTimesteps+j+thisthread*All.NMassFormed]
				  += MassFormed[H[ihalo].imass*All.NTimesteps+j+thisthread*All.NMassFormed]
				     + ICMFormed[H[ihalo].imass*All.NTimesteps+j+thisthread*All.NMassFormed];
			}
			// descid points to where icm is distributed
			H[iprog].descid   = H[H[ihalo].ihost].haloid;
#endif
			//If stripped remove galaxy from run
			H[ihalo].gone = 1;
			if (H[ihalo].idesc >= 0) H[H[ihalo].idesc].gone = 1;
			H[iprog].idescgal = -1;
			continue;
		}

//////////////////////////////////////////////////////////////////////////////////////////////////////
//    Satellite merging                                                                             //
//////////////////////////////////////////////////////////////////////////////////////////////////////
		//If the dynamical friction time has run out merge the orphan
		if (H[ihalo].tdf < t_now && H[ihalo].type == 2)
		{
			//Check if the main galaxy is gone
			if (H[H[ihalo].imain].gone == 1)
			{
#ifdef COMPUTE_ICM
				//Go through all formed populations up to now
				for (j = 0; j <= H[ihalo].iscale; j++)
				{
					//Add the stellar mass of the satellite to the main halo's ICM
					ICMFormed[H[H[ihalo].ihost].imass*All.NTimesteps+j+thisthread*All.NMassFormed]
					  += MassFormed[H[ihalo].imass*All.NTimesteps+j+thisthread*All.NMassFormed]
					     + ICMFormed[H[ihalo].imass*All.NTimesteps+j+thisthread*All.NMassFormed];
				}
#endif
				//Remove satellite galaxy from run
				H[ihalo].gone = 1;
				if (H[ihalo].idesc >= 0) H[H[ihalo].idesc].gone = 1;
				H[iprog].idescgal = -1;
				continue;
			}
			//Otherwise add all formed populations up to now to main galaxy and its ICM
			for (j = 0; j <= H[ihalo].iscale; j++)
			{
				//Add mass to stellar mass of main descendant
				MassFormed[H[H[ihalo].imain].imass*All.NTimesteps+j+thisthread*All.NMassFormed]
				  += (1.0 - P.Fesc) * MassFormed[H[ihalo].imass*All.NTimesteps+j+thisthread*All.NMassFormed];
#ifdef COMPUTE_ICM
				//Add mass to ICM of main descendant
				ICMFormed[H[H[H[ihalo].imain].ihost].imass*All.NTimesteps+j+thisthread*All.NMassFormed]
				  += P.Fesc * MassFormed[H[ihalo].imass*All.NTimesteps+j+thisthread*All.NMassFormed];
#endif
			}
			//Set the galaxy descandant ID of the progenitor to the first galaxy's index (orphans always have a progenitor)
			iprog = H[ihalo].iprog;
			H[iprog].idescgal = H[ihalo].imain;
			//After merging set the descendant ID to that of main halo(galaxy)
			H[iprog].descid   = H[H[ihalo].imain].haloid;
			//Increase the number of progenitors for the galaxy it merges with
			H[H[ihalo].imain].np++;
			//Set the progenitor's mmp to zero
			H[iprog].mmp = 0;
			//Go to the first progenitor of the halo's descendant
			iprog = H[H[ihalo].imain].iprog;
			//While the icoprog is already used set the iprog to the next co-progenitor
			while (H[iprog].icoprog >= 0) iprog = H[iprog].icoprog;
			//As soon as icoprog is free write the halo index to the icoprog
			H[iprog].icoprog = H[ihalo].iprog;
			//Mark the galaxy as gone and if has a descendant mark it as gone as well
			H[ihalo].gone = 1;
			if (H[ihalo].idesc >= 0) H[H[ihalo].idesc].gone = 1;
			continue;
		}

//////////////////////////////////////////////////////////////////////////////////////////////////////
//    Compute total stellar mass and ICM from populations                                           //
//////////////////////////////////////////////////////////////////////////////////////////////////////
		//Compute the stellar mass from the formed stellar populations
		H[ihalo].mstar = stellar_mass(H[ihalo].imass,H[ihalo].iscale,thisthread);
#ifdef COMPUTE_ICM
		//Compute the intra cluster mass from the formed stellar populations (only for main haloes)
		if (H[ihalo].type == 0) H[ihalo].icm = intra_cluster_mass(H[ihalo].imass,H[ihalo].iscale,thisthread);
		else H[ihalo].icm = 0.0;
#endif

//////////////////////////////////////////////////////////////////////////////////////////////////////
//    Add this galaxy to the global statistics                                                      //
//////////////////////////////////////////////////////////////////////////////////////////////////////
		add_galaxy_to_statistics(ihalo,thisthread);

	}//End loop through all haloes - increment ihalo

//////////////////////////////////////////////////////////////////////////////////////////////////////
//    Reset MassFormed and ICMFormed to zero                                                        //
//////////////////////////////////////////////////////////////////////////////////////////////////////
	//After the tree/forest is done reset the MassFormed array to 0
	for (i = 0; i < imass * All.NTimesteps; i++)
	{
		MassFormed[i+thisthread*All.NMassFormed] = 0.0;
#ifdef COMPUTE_ICM
		ICMFormed[i+thisthread*All.NMassFormed]  = 0.0;
#endif
	}

}


/*! \brief This function is the heart of the model
 *
 *  This function returns the instantaneous baryon conversion efficiency for a given halo mass and scale
 *  factor. The parameters are those stored on each task so it should be made sure that in a given
 *  universe all tasks have the same values for the parameters.
 *  The instantaneous baryon conversion efficiency is defined as:
 *  \f[ \epsilon(M,z) = 2 \;\epsilon_\mathrm{N} \left[ \left(\frac{M}{M_1}\right)^{-\beta} +
 *  \left(\frac{M}{M_1}\right)^{\gamma}\right]^{-1} \; ,\f]
 *  where \f$M_1\f$,\f$\epsilon_\mathrm{N}\f$, \f$\beta\f$, and \f$\gamma\f$ depend on redshift and are
 *  given by:
 *  \f{eqnarray*}{
 *  \log_{10} M_1(z) &=& M_0 + M_\mathrm{z}\frac{z}{z+1}\\
 *  \epsilon_\mathrm{N}(z) &=& \epsilon_0 + \epsilon_\mathrm{z}\frac{z}{z+1}\\
 *  \beta(z) &=& \beta_0 + \beta_\mathrm{z}\frac{z}{z+1}\\
 *  \gamma(z) &=& \gamma_0 + \gamma_\mathrm{z}\frac{z}{z+1} \; .
 *  \f}
 *
 *  \param hmass Virial mass of the dark matter halo in code units
 *  \param a Scale factor at which the dark matter halo is located
 */
float sfe(float hmass, float a)
{
	//Possibly make static inline?
	double sfe_e,sfe_m,sfe_b,sfe_g;
	sfe_e = P.Epsilon0    + P.EpsilonZ  * (1.-a);
	sfe_m = pow(10., P.M0 + P.MZ        * (1.-a))/All.m_unit;
	sfe_b = P.Beta0       + P.BetaZ     * (1.-a);
	sfe_g = P.Gamma0      + P.GammaZ    * (1.-a);
	return (float)(2.0*sfe_e/(pow(hmass/sfe_m,-sfe_b)+pow(hmass/sfe_m,sfe_g)));
}


/*! \brief This function calculates the stellar mass from the #MassFormed array
 *
 *  This function returns the stellar mass in a given branch imass. The elements in the #MassFormed array
 *  are summed up from 0 to the current time specified by the index of the scale factor array scalemax.
 *  Each element is multiplied by the fraction of stars still left at this time which is stored in #MassLeft.
 *
 *  \param imass The branch of the tree for which the stellar mass is calculated
 *  \param scalemax The index of the scale factor array up to which the stellar mass is calculated.
 *  \param thisthread Number of the OpenMP thread that is responsible
 */
float stellar_mass(int imass, int scalemax, int thisthread)
{
	int i;
	float mass = 0.0;
	//Set the thread number
	for (i = 0; i <= scalemax; i++)
		mass += MassFormed[imass*All.NTimesteps+i+thisthread*All.NMassFormed] * MassLeft[scalemax*All.NTimesteps+i];
	return mass;
}


/*! \brief This function calculates the ICM from the #ICMFormed array
 *
 *  This function returns the ICM in a given branch imass. The elements in the #ICMFormed array
 *  are summed up from 0 to the current time specified by the index of the scale factor array scalemax.
 *  Each element is multiplied by the fraction of stars still left at this time which is stored in #MassLeft.
 *
 *  \param imass The branch of the tree for which the ICM is calculated
 *  \param scalemax The index of the scale factor array up to which the ICM is calculated.
 *  \param thisthread Number of the OpenMP thread that is responsible
 */
float intra_cluster_mass(int imass, int scalemax, int thisthread)
{
	int i;
	float mass = 0.0;
	for (i = 0; i <= scalemax; i++)
		mass += ICMFormed[imass*All.NTimesteps+i+thisthread*All.NMassFormed] * MassLeft[scalemax*All.NTimesteps+i];
	return mass;
}


/*! \brief This function calculates the dynamical friction time of an orphan
 *
 *  This function returns the dynamical friction time of an orphan ihalo merging with a main galaxy imain.
 *  It either uses the default dynamical friction formula by
 *  <a href="http://adsabs.harvard.edu/abs/1987gady.book.....B"> Binney & Tremain (1987)</a> or one based on
 *  B<a href="https://arxiv.org/abs/0707.2960"> Boylan-Kolchin et al. (2008)</a>, which can be enabled by the
 *  option #DF_USE_BK08. For BT97 a factor of 2 is used to account for the longer merging times
 *  found by others. For BK08 the power d of the satellite radius is set to 2 since local satellite mass is used
 *  (not the one at infall). For both, the stellar mass of both galaxies is taken into account as well when
 *  calculating the dynamical frition time.

 *  \param ihalo The index of the orphan halo that is merging
 *  \param imain The index of the main (non-orphan) halo the orphan is merging with
 */
float get_tdf(int ihalo, int imain)
{
#define G         4.30071e-9   ///< Gravitational constant in (km/s)^2 * Mpc / Msol
#define MPC_IN_KM 3.08568e19   ///< Megaparsecs in kilometers
#define YR_IN_SEC 31556926.0   ///< Year in seconds

	double vvir, rsat, mu, eta, tau, clog, dx, dy, dz;

	//Compute virial velocity, current radius, satellite mass and Coulomb logarithm
	vvir = sqrt(G*H[imain].mvir*All.m_unit/H[imain].rvir/All.x_unit);
	dx   = NEAREST(H[ihalo].pos[0]-H[imain].pos[0]);
	dy   = NEAREST(H[ihalo].pos[1]-H[imain].pos[1]);
	dz   = NEAREST(H[ihalo].pos[2]-H[imain].pos[2]);
	rsat = sqrt(dx*dx+dy*dy+dz*dz)*H[imain].a;
	clog = log(1.0+(H[imain].mvir+H[imain].mstar)/(H[ihalo].mvir+H[ihalo].mstar));
	mu   = (H[ihalo].mvir+H[ihalo].mstar)/(H[imain].mvir+H[imain].mstar);
	tau  = (H[imain].rvir*All.x_unit/vvir) * MPC_IN_KM / YR_IN_SEC / All.t_unit;

#ifdef DF_USE_BK08
	eta  = 0.5 + 0.214 * get_gaussian_random_number(ihalo+imain); //Zentner et al. (2005)
	if (eta < 0.0) eta = 0.0;
	if (eta > 1.0) eta = 1.0;
	//Use a modified version of the Boylan-Kolchin formula (see section 3.3.2).
	//d=2 is used because we do not take the satellite mass when it was at the virial radius but the current one
	const float a = 0.216;
	const float b = 1.3;
	const float c = 1.9;
	const float d = 2.0;
#else
	eta = 1.0;
	//Use formulat from Binney & Tremain with a factor of 2
	const float a = 2.0 * 1.17;
	const float b = 1.0;
	const float c = 0.0;
	const float d = 2.0;
#endif

	return a / clog * pow(mu,-b) * exp(c*eta) * pow(rsat/H[imain].rvir,d) * tau;

#undef YR_IN_SEC
#undef MPC_IN_KM
#undef G
}
