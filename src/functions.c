///////////////////////////////////////////////////////////////////////////////////////////////////
// Emerge code - File functions.c                                                                //
// Contains general functions                                                                    //
///////////////////////////////////////////////////////////////////////////////////////////////////
///
/// \file functions.c
/// \brief Contains general functions
///
/// This file contains general functions that are used to compute various quantities.
/// There are several functions that compare two values and can be used in the qsort algorithm.
/// The file also contains functions to search numbers in a sorted list.
/// Also included are cosmology functions to calculate the cosmic time.
/// The random number generator or table is accessed through a function in this file.
/// A function to compute the current CPU time is included as well.
///
///////////////////////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "allvars.h"
#include "proto.h"


/*! \brief This function compares two floating point numbers
 *
 *  \param a First floating point number
 *  \param b Second floating point number
 */
int compare_float(const void *a, const void *b)
{
	float c = *(const float*) a;
	float d = *(const float*) b;
	return (c > d) - (c < d);
}

/*! \brief This function compares the IDs of two search_list_IDs structures
 *
 *  \param a First search_list_IDs structure
 *  \param b Second search_list_IDs structure
 */
int compare_id(const void *a, const void *b)
{
	const struct search_list_IDs *c = a;
	const struct search_list_IDs *d = b;
	return (c->ID > d->ID) - (c->ID < d->ID);
}


/*! \brief This function compares the loads of two task_load structures
 *
 *  \param a First task_load structure
 *  \param b Second task_load structure
 */
int compare_load(const void *a, const void *b)
{
	const struct task_list *c = a;
	const struct task_list *d = b;
	return (c->LoadInTask < d->LoadInTask) - (c->LoadInTask > d->LoadInTask);
}


/*! \brief This function compares the scale factors of two haloes
 *
 *  If they are the same compare the haloid of the haloes.
 *  This is done to make the sorting stable.
 *
 *  \param a First halo structure
 *  \param b Second halo structure
 */
int compare_scale(const void *a, const void *b)
{
	const struct halo *c = a;
	const struct halo *d = b;
	int scale = (c->a > d->a) - (c->a < d->a);
	return scale != 0 ? scale : (d->haloid > c->haloid) - (d->haloid < c->haloid);
}


/*! \brief This function compares the imain of two haloes
 *
 *  This is done to sort the haloes by the task they need to go.
 *  If they are the same compare the iscale of the haloes which gives the original index.
 *  This is done to make the sorting stable.
 *
 *  \param a First halo structure
 *  \param b Second halo structure
 */
int compare_totask(const void *a, const void *b)
{
	const struct halo *c = a;
	const struct halo *d = b;
	int totask = c->imain - d->imain;
	return totask != 0 ? totask : (c->iscale > d->iscale) - (c->iscale < d->iscale);
}


/*! \brief This function compares the forestid of two haloes
 *
 *  This is done to sort the haloes by the forest ID
 *  If they are the same compare the iscale of the haloes which gives the original index.
 *  This is done to make the sorting stable.
 *
 *  \param a First halo structure
 *  \param b Second halo structure
 */
#ifdef COMPUTE_ICM
int compare_forest(const void *a, const void *b)
{
	const struct halo *c = a;
	const struct halo *d = b;
	int forest = (c->forestid > d->forestid) - (c->forestid < d->forestid);
	return forest != 0 ? forest : (c->iscale > d->iscale) - (c->iscale < d->iscale);
}
#endif


/*! \brief This function compares the tree or forest length of two haloes
 *
 *  This is done to sort the haloes by their tree/forest length.
 *  First compare the idescgal of the haloes which gives the tree or forest length.
 *  If they are the same compare the icoprog of the haloes which gives the original index.
 *  This is done to make the sorting stable.
 *
 *  \param a First halo structure
 *  \param b Second halo structure
 */
int compare_size(const void *a, const void *b)
{
	const struct halo *c = a;
	const struct halo *d = b;
	int size = (c->idescgal < d->idescgal) - (c->idescgal > d->idescgal);
	return size != 0 ? size : (c->icoprog > d->icoprog) - (c->icoprog < d->icoprog);
}


/*! \brief This function searches an array of search_list_IDs structures
 *
 *  The search_list_IDs structure array a with length n is searched for the ID key.
 *  If it is found the position of the key in the unsorted list is returned.
 *  If it is found a value of -99 is returned.
 *
 *  \param a search_list_IDs structure that will be scanned
 *  \param n Size of the search_list_IDs structure array
 *  \param key ID that is searched
 */
int binary_search_id(struct search_list_IDs *a, int n, IDType key)
{
	int imin = 0;
	int imax = n-1;
	int imid = (imin+imax)/2;
	// continually narrow search until just one element remains
	while (imin <= imax)
	{
		// determine which subarray to search
		if (a[imid].ID < key)
			// change min index to search upper subarray
			imin = imid + 1;
		else if ((imin==imax)&&(a[imid].ID!=key))
			return -99;
		else if (a[imid].ID == key)
			// key found at index imid
			return a[imid].position;
		else
			// change max index to search lower subarray
			imax = imid - 1;
		//calculate the midpoint for roughly equal partition
		imid = (imin+imax)/2;
	}
	// key not found
	return -99;
}


/*! \brief This function searches an array of floats
 *
 *  The float array a with length n is searched for the key.
 *  If it is found the index of the key is returned.
 *  If it is found a value of -99 is returned.
 *
 *  \param a float array that will be scanned
 *  \param n Size of the float array
 *  \param key that is searched
 */
int binary_search_float(float *a, int n, float key)
{
	int imin = 0;
	int imax = n-1;
	int imid = (imin+imax)/2;
	// continually narrow search until just one element remains
	while (imin <= imax)
	{
		// determine which subarray to search
		if (a[imid] < key)
			// change min index to search upper subarray
			imin = imid + 1;
		else if ((imin==imax)&&(a[imid]!=key))
			return -99;
		else if (a[imid] == key)
			// key found at index imid
			return imid;
		else
			// change max index to search lower subarray
			imax = imid - 1;
		//calculate the midpoint for roughly equal partition
		imid = (imin+imax)/2;
	}
	// key not found
	return -99;
}


/*! \brief This function computes the cosmic time as a function of the scale factor
 *
 *  The age of the universe at any scale factor (cosmic time) is computed with the approximation
 *  given in <a href="http://adsabs.harvard.edu/abs/1992ARA%26A..30..499C"> Carroll, Press & Turner (1992)</a>,
 *  equations 17/18.
 *
 *  \param a The scale factor at which the cosmic time is computed
 */
float cosmictime(float a)
{
	float z,omega_L,f,absf,arg,t;
	z = 1.0/a-1.0;
	if (z < 0.0) return -1;
	if (All.Omega_0+All.Omega_Lambda_0==1.0)
	{
		omega_L = a*a*a*All.Omega_Lambda_0/(a+All.Omega_0*(1.0-a)+All.Omega_Lambda_0*(a*a*a-a));
		f = 0.7*Omega(z) - 0.3*omega_L + 0.3;
		if (1.0-f < 0) absf = f-1.0;
		else absf = 1.0-f;
		arg = sqrt(absf/f);
		if (f < 1) t = (2.0*All.Hubbletime/(3.0*Epeebles(z)*sqrt(absf)))*asinh(arg);
		else t = (2.0*All.Hubbletime/(3.0*Epeebles(z)*sqrt(absf)))*asin(arg);
		return t/All.t_unit;
	}
	else if (All.Omega_0+All.Omega_Lambda_0<1.0)
	{
		return All.Hubbletime*(sqrt(1+All.Omega_0*z)/((1-All.Omega_0)*(1+z)) - (All.Omega_0/(2.*pow(1-All.Omega_0, 1.5))) * acosh(1+2.*(1-All.Omega_0)/(All.Omega_0*(1+z))))/All.t_unit;
	}
	else return -1;
}


/*! \brief This function computes the matter density at a given redshift
 *
 *  \param z The redshift at which the matter density is computed
 */
float Omega(float z)
{
	//omega as a function of redshift
	if (z < 0.005)
		return All.Omega_0;
	else if (All.Omega_0+All.Omega_Lambda_0==1.0)
		return All.Omega_0*pow(1+z, 3.0)/(1-All.Omega_0 + pow(1+z, 3.0)*All.Omega_0);
	else if (All.Omega_0+All.Omega_Lambda_0<1.0)
		return All.Omega_0*(1+z)/(1 + All.Omega_0*z);
	return 0.0;
}


/*! \brief This function computes the peebles function \f$E(z) = H(z)/\mathrm{H}_0\f$ at a given redshift
 *  (<a href="http://press.princeton.edu/titles/5263.html"> Peebles 2010</a>, p. 312)
 *
 *  \param z The redshift at which the matter density is computed
 */
float Epeebles(float z)
{
	float omega_R, omega_Lambda;
	omega_R = 0.0;
	omega_Lambda = 0.0;
	if (All.Omega_0+All.Omega_Lambda_0==1.0) omega_Lambda = 1.0-All.Omega_0;
	else if (All.Omega_0+All.Omega_Lambda_0<1.0) omega_R = 1.0-All.Omega_0;
	return sqrt(All.Omega_0*(1.0+z)*(1.0+z)*(1.0+z)+omega_R*(1.0+z)*(1.0+z)+omega_Lambda);
}


/*! \brief This function computes the dynamical time using the virial overdensity by
 *  <a href="https://arxiv.org/abs/astro-ph/9710107"> Bryan & Norman (1998)</a>
 *
 *  \param a The scale factor at which the cosmic time is computed
 */
float tdyn(float a)
{
	float x,virdens,z;
	z = 1.0/a-1.0;
	x = (Omega(z)-1.0);
	virdens = 18.0*M_PI*M_PI+82.0*x-39.0*x*x;
	return All.Hubbletime/All.t_unit/sqrt(virdens/2.0)/Epeebles(z);
}


/*! \brief This function returns a gaussian random number with sigma = 1.0
 *
 *  \param index The index of the array that will get a random number
 */
float get_gaussian_random_number(int index)
{
#if (RANDOM_NUMBER_TABLE > 0)
	return RndTableGaussian[(index % RANDOM_NUMBER_TABLE)];
#else
	return gsl_ran_gaussian(rng_gaussian,1.0);
#endif
}


/*! \brief This function returns a uniform random number between 0 and 1
 *
 *  \param index The index of the array that will get a random number
 */
float get_uniform_random_number(int index)
{
#if (RANDOM_NUMBER_TABLE > 0)
	return RndTableUniform[(index % RANDOM_NUMBER_TABLE)];
#else
	return gsl_rng_uniform(rng_uniform);
#endif
}


/*! \brief This function returns the time in seconds
 */
double second(void)
{
	return MPI_Wtime();
}
