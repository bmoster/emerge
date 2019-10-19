///////////////////////////////////////////////////////////////////////////////////////////////////
// Emerge code - File clustering.c                                                               //
///////////////////////////////////////////////////////////////////////////////////////////////////
///
/// \file clustering.c
/// \brief Contains functions that compute the two-point correlation function and its projection
///
/// This file contains all functions that are needed to compute the two-point galaxy correlation
/// function xi(r) and the projected galaxy correlation function wp(rp). For this a kd-tree is
/// built which is then used to find pairs in parallel for each univere.
///
///////////////////////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_spline.h>
#include "allvars.h"
#include "proto.h"


/*! \brief This function computes the projected correlation functions
 *
 *  The function first allocates memory for a an array that stores the properties of all galaxies on this
 *  task that are used for the correlation functions. Going through all haloes on this task the galaxy
 *  array is filled with positions, masses and SFRs. If the memory is not enough more memory is allocated
 *  dynamically. Using the galaxy arrays on each task the projected correlation functions in each stellar
 *  mass bin is computed.
 */
int compute_wp(void)
{
	int i, j;
	int igalaxy = 0;
	int ngaltot = NALLOCATION;
	struct galaxy *localgal = emalloc("LocalGalaxies", ngaltot * sizeof(struct galaxy));

	//Go through all haloes on this task
	for (i = 0; i < Nhalos; i++)
	{
		//Check if the galaxy is above the minimum mass and below the maximum mass for clustering
		if (H[i].iscale == All.Wpscale && H[i].mstar >= All.wpmmin && H[i].mstar < All.wpmmax)
		{
			//Write all galaxy properties to a local array (localgal)
			localgal[igalaxy].x[0] = H[i].pos[0];
			localgal[igalaxy].x[1] = H[i].pos[1];
			localgal[igalaxy].x[2] = H[i].pos[2];
			localgal[igalaxy].mass = H[i].mstar;
			localgal[igalaxy].sfr  = H[i].sfr;
			igalaxy++;
			//If the allocated array becomes too small we reallocate
			if (igalaxy >= ngaltot)
			{
				ngaltot += NALLOCATION;
				localgal = (struct galaxy *) erealloc(localgal, ngaltot * sizeof(struct galaxy));
			}
		}
	}

	//Record the total number of galaxies for clustering on this task
	ngaltot = igalaxy;

	//For each mass bin determine the projected correlation function using the local galaxy array
	for (j = 0; j < All.Nwpset; j++) compute_wp_bin(j,ngaltot,localgal);

	//Free the galaxy array
	efree(localgal);

	return 1;

}


/*! \brief This function computes the projected correlation functions
 *
 *  The function first collects all galaxies in this stellar mass bin from the local galaxy array. Then
 *  the kd-tree is built. The radii for the 3d-correlation function are set logarithmically from
 *  \f$ r_\mathrm{min} \f$ to \f$ \sqrt{r_\mathrm{max}^2 + \pi_\mathrm{max}^2} \f$, where
 *  \f$ r_\mathrm{min} \f$ is the smallest radius of the corresonding observed projected correlation
 *  function, and \f$ r_\mathrm{max} \f$ is defined by the maximum of the largest radius of the corresonding observed
 *  projected correlation function and #WP_RMAX times the box length: \f$ r_\mathrm{max}  = \mathrm{min}
 *  (r_\mathrm{p, max}, \mathrm{WP\_RMAX} \cdot L_\mathrm{box})\f$.
 *
 *  The 3d-auto-correlation function is then computed as \f[ \xi(r) = \frac{\mathrm{DD}(r)}{\mathrm{RR}(r)} -1 ,\f]
 *  where \f$ \mathrm{DD}(r)\f$ is calculated using the parallel pair counter funcion #tree_pairs, and
 *  \f$ \mathrm{RR}(r) \f$ is given by the natural estimator (removing the double counting): \f[ \mathrm{RR}(r) =
 *  \frac{N(N-1)}{2}\frac{V}{L_\mathrm{box}^3} = N(N-1) \frac{2\pi (r^2 \mathrm{d}r +
 *  \mathrm{d}r/12)}{L_\mathrm{box}^3} , \f] where V is the volume of the radial bin with width
 *  \f$ \mathrm{d}r \f$, and \f$ N \f$ is the number of galaxies in this mass bin.
 *
 *  Next, the computed correlation function is interpolated using #WP_RBINS_INT bins. Finally the projected
 *  correlation function at each observed bin is calculated: \f[ w_\mathrm{p} (r_\mathrm{p}) =
 *  \int_0^{\pi_\mathrm{max}} 2 \xi(r_\mathrm{p},\pi) \mathrm{d}\pi = 2 \int_{r_\mathrm{p}}^{\sqrt{r_\mathrm{p}^2
 *  + \pi_\mathrm{max}^2}} \frac{r \xi(r)}{\sqrt{r^2-r_\mathrm{p}^2}} \; \mathrm{d}r , \f] where \f$
 *  \pi_\mathrm{max} \f$ is the maximum distance along the line of sight.
 *
 *  \param iwpbin The index of the stellar mass bin of the projected correlation functions
 *  \param ngaltot The number of galaxies on this task in the array localgal
 *  \param *localgal Local galaxy array that holds the properties of all galaxies that are used for correlation functions
 */
int compute_wp_bin(int iwpbin, int ngaltot, struct galaxy *localgal)
{
	int i, j, ngal, nspline;
	float rmin, rmax, dlogr, logr, r1, r2, norm, rr, dd, splinemin, splinemax;
	float *r,*xi,*xierr;
	double *x,*ya,*yb;
	struct kd_tree tree;
	gsl_spline *splineX, *splineE;
	gsl_interp_accel *xacc;

	//Get all galaxies in the corresponding stellar mass bin on every task
	ngal = collect_galaxies(iwpbin,ngaltot,localgal);
	//There are no pairs for 1 or 0 galaxies (and negative ngal states an error)
	if (ngal <= 1)
	{ //Set all model values to zero and all errors to a large number
		for (i = WpSet[iwpbin].offset; i < WpSet[iwpbin].offset + WpSet[iwpbin].ndata; i++)
		{
			Wp[i].mod_y     = 0.0;
			Wp[i].mod_sigma = 1.e10;
		}
		if (ngal >= 0) efree(Galpos);
		return 0;
	}

	//Construct the kd-tree for the galaxies in this mass bin (return 0 if not possible)
	if (!(construct_tree(&tree, Galpos, ngal, TREEDIM)))
	{ //Set all model values to zero and all errors to a large number
		for (i = WpSet[iwpbin].offset; i < WpSet[iwpbin].offset + WpSet[iwpbin].ndata; i++)
		{
			Wp[i].mod_y     = 0.0;
			Wp[i].mod_sigma = 1.e10;
		}
		efree(tree.nodes);
		efree(Galpos);
		return 0;
	}

	//Check if we have enough memory left - otherwise return 0
	if ((long long)((3*WP_RBINS*sizeof(float)) * sizeof(struct kd_pos)) > (long long)(FreeBytes))
	{ //Set all model values to zero and all errors to a large number
		for (i = WpSet[iwpbin].offset; i < WpSet[iwpbin].offset + WpSet[iwpbin].ndata; i++)
		{
			Wp[i].mod_y     = 0.0;
			Wp[i].mod_sigma = 1.e10;
		}
		efree(tree.nodes);
		efree(Galpos);
		return 0;
	}

	//Allocate the radius, correlation function, and correlation function error arrays
	r     = emalloc("RADIUS", WP_RBINS * sizeof(float));
	xi    = emalloc("XI    ", WP_RBINS * sizeof(float));
	xierr = emalloc("XI_ERR", WP_RBINS * sizeof(float));

	//Set the minimum radius, maximum radius, log bin size, and the normalisation for the correlation function
	rmin  = Wp[WpSet[iwpbin].offset].obs_x;
	rmax  = sqrt(min(sqr(Wp[WpSet[iwpbin].offset+WpSet[iwpbin].ndata-1].obs_x),sqr(WP_RMAX*All.Lbox)) + sqr(WpSet[iwpbin].cut));
	dlogr = (log10(rmax)-log10(rmin))/(WP_RBINS-1);
	norm  =  2.0 * M_PI / 3.0 / cub(All.Lbox);

	//Go through all radial bins (double counting accounted for in both dd and rr)
	for (i = 0; i < WP_RBINS; i++)
	{
		//Set the current log radius, radius, and bin edges
		logr  = log10(rmin) + (float)(i) * dlogr;
		r[i]  = pow(10.0, logr);
		r1    = pow(10.0, log10(rmin)+((float)(i)-0.5)*dlogr);
		r2    = pow(10.0, log10(rmin)+((float)(i)+0.5)*dlogr);
		//Compute the random pair count from the natural estimator
		rr    = norm * (float)(ngal * (ngal-1)) * (cub(r2) - cub(r1));
		//Compute the pair count using the parallel version of the dualtree algorithm
		dd    = (float)(tree_pairs(tree,Galpos,r1,r2));
		//Get the 3d auto-correlation function
		xi[i] = dd/rr - 1.0;
		//Set the poisson error for the correlation function
		if (dd > 0) xierr[i] = xi[i]/sqrt(dd);
		else
		{ //If the correlation function is zero or negative (error) set it to zero and its error to 1e10
			xi[i]    = 0.0;
			xierr[i] = 1.0e10;
		}
	}

	//Check if we have enough memory left - otherwise return 0
	if ((long long)((3*WP_RBINS*sizeof(float)) * sizeof(struct kd_pos)) > (long long)(FreeBytes))
	{ //Set all model values to zero and all errors to a large number
		for (i = WpSet[iwpbin].offset; i < WpSet[iwpbin].offset + WpSet[iwpbin].ndata; i++)
		{
			Wp[i].mod_y     = 0.0;
			Wp[i].mod_sigma = 1.e10;
		}
		efree(xierr);
		efree(xi);
		efree(r);
		efree(tree.nodes);
		efree(Galpos);
		return 0;
	}

	//Allocate arrays for logarithmic spline interpolation
	x  = emalloc("X",  WP_RBINS * sizeof(double));
	ya = emalloc("YA", WP_RBINS * sizeof(double));
	yb = emalloc("YB", WP_RBINS * sizeof(double));

	//Initialise the number of points for the spline interpolation
	nspline   =  0;
	splinemin =  1.e30;
	splinemax = -1.e30;
	//Set the interpolation arrays to the logarithmic radius, xi and xi error
	for (i = 0; i < WP_RBINS; i++)
	{ //Check if the correlation function is positive
		if (xi[i] > 0)
		{ //Set the interplotation to logarithmic
			x[nspline]  = log10(r[i]);
			ya[nspline] = (double)(log10(xi[i]));
			yb[nspline] = (double)(log10(xierr[i]));
			//Determine the minimum and maximum radius using all valid points
			if (x[nspline] < splinemin) splinemin = x[nspline];
			if (x[nspline] > splinemax) splinemax = x[nspline];
			nspline++;
		}
	}

	//If there are too few points to create the splines
	if (nspline < 3)
	{ //Set all model values to zero and all errors to a large number
		for (i = WpSet[iwpbin].offset; i < WpSet[iwpbin].offset + WpSet[iwpbin].ndata; i++)
		{
			Wp[i].mod_y     = 0.0;
			Wp[i].mod_sigma = 1.e10;
		}
		//Free all arrays
		efree(yb);
		efree(ya);
		efree(x);
		efree(xierr);
		efree(xi);
		efree(r);
		efree(tree.nodes);
		efree(Galpos);
		//And return
		return 0;
	}

	//Set the logarithmic bin size and allocate the splines and acceleration array
	dlogr   = (log10(rmax)-log10(rmin))/(WP_RBINS_INT-1);
	splineX = gsl_spline_alloc(gsl_interp_cspline,nspline);
	splineE = gsl_spline_alloc(gsl_interp_cspline,nspline);
	xacc    = gsl_interp_accel_alloc();

	//Initialise the splines for the correlation function and its errors
	gsl_spline_init(splineX, x, ya, nspline);
	gsl_spline_init(splineE, x, yb, nspline);

	//Go through all radii, evaluate the splines and go back to non-logarithmic values
	for (i = WP_RBINS_INT * iwpbin; i < WP_RBINS_INT * (iwpbin+1); i++)
	{ //Set the logarithmic radius
		Radius[i]     = log10(rmin) + (float)(i - WP_RBINS_INT * iwpbin) * dlogr;
		//Check if the radius is in the interpolation range
		if (Radius[i] > splinemin && Radius[i] < splinemax)
		{ //If it is evaluate the splines and go back to non-logarithmic values
			Modelxi[i]    = gsl_spline_eval(splineX, Radius[i], xacc);
			Modelxierr[i] = gsl_spline_eval(splineE, Radius[i], xacc);
			Modelxi[i]    = pow(10.0, Modelxi[i]);
			Modelxierr[i] = pow(10.0, Modelxierr[i]);
		}
		else
		{ //Otherwise we can't evaluate this radius and set the correlation function to zero and its error to 1e10
			Modelxi[i]    = 0.0;
			Modelxierr[i] = 1.e10;
		}
		//Go back to non-logarithmic radii
		Radius[i]     = pow(10.0, Radius[i]);
	}

	// Integrate xi to get projected correlation function wp
	for (i = WpSet[iwpbin].offset; i < WpSet[iwpbin].offset + WpSet[iwpbin].ndata; i++)
	{ //Set the initial values to zero
		Wp[i].mod_y     = 0.0;
		Wp[i].mod_sigma = 0.0;
		for (j = WP_RBINS_INT * iwpbin; j < WP_RBINS_INT * (iwpbin+1); j++)
		{ //Integrate from rp to sqrt(rp^2+pimax^2)
			if (Radius[j] > Wp[i].obs_x && sqr(Radius[j]) <= sqr(WpSet[iwpbin].cut) + sqr(Wp[i].obs_x))
			{ //Logarithmic integration, hence the factor log(10)
				Wp[i].mod_y     += 2.0 * log(10.) * sqr(Radius[j]) * dlogr * Modelxi[j] / sqrt(sqr(Radius[j])-sqr(Wp[i].obs_x));
				Wp[i].mod_sigma += 2.0 * log(10.) * sqr(Radius[j]) * dlogr * Modelxierr[j] / sqrt(sqr(Radius[j])-sqr(Wp[i].obs_x));
			}
		}
	}

	//Free gsl structures
	gsl_interp_accel_free(xacc);
	gsl_spline_free(splineE);
	gsl_spline_free(splineX);

	//Free all arrays
	efree(yb);
	efree(ya);
	efree(x);

	efree(xierr);
	efree(xi);
	efree(r);

	//Free the tree nodes and the position array
	efree(tree.nodes);
	efree(Galpos);

	return 1;

}


/*! \brief This function collects all galaxies in a stellar mass bin and copies them to all tasks
 *
 *  The function first counts the number of galaxies on each task that are in this mass bin. This number is
 *  then sent by each worker to the master task, where the numbers are summed up and sent back to the workers.
 *  Each task allocates the global positions array and a helper array to send data. The positions of all galaxies
 *  in this mass bin are then sent to the master task using the helper array, where they are stored in the
 *  global positions array. Finally the global positions array is sent back to each worker, and the total
 *  number of galaxies in this mass bin is returned.
 *
 *  \param iwpbin The index of the stellar mass bin of the projected correlation functions
 *  \param ngaltot The number of galaxies on this task in the array localgal
 *  \param *localgal Local galaxy array that holds the properties of all galaxies that are used for correlation functions
 */
int collect_galaxies(int iwpbin, int ngaltot, struct galaxy *localgal)
{

	int i, j, task, ngalaxies, nthistask, nrecv, offset, Nmax;
	struct kd_pos *possend;
	MPI_Status status;

	//Initialise the number of galaxies in this bin to zero
	nthistask = 0;
	//Go through all galaxies in the local array
	for (i = 0; i < ngaltot; i++)
	{
		//If they are in this mass bin increment the counter
		if (localgal[i].mass >= WpSet[iwpbin].min && localgal[i].mass < WpSet[iwpbin].max)
			nthistask++;
	}

	//Initialise the total number of galaxies in this mass bin with the number on this task
	ngalaxies = nthistask;
	Nmax      = nthistask;

	//Send the number of galaxies on each task to the master task
	if (ThisTask != MasterTask)
		//Each worker task sends nthistask
		MPI_Ssend(&nthistask, 1, MPI_INT, MasterTask, TAG_NHALOS, MPI_COMM_WORLD);
	else
	{
		//The Master task goes through all workers
		for (task = MasterTask + 1; task < MasterTask + All.NTaskPerUniverse; task++)
		{
			//Receive the number of galaxies on this worker and increment the total number
			MPI_Recv(&nrecv, 1, MPI_INT, task, TAG_NHALOS, MPI_COMM_WORLD, &status);
			ngalaxies += nrecv;
			if (nrecv > Nmax) Nmax = nrecv;
		}
	}

	//Send the total number of galaxies back to the workers
	if (ThisTask != MasterTask)
		//Each worker task receives the total number of galaxies
		MPI_Recv(&ngalaxies, 1, MPI_INT, MasterTask, TAG_NHALOS, MPI_COMM_WORLD, &status);
	else
	{
		//The Master task go through all workers and sends them the total number of galaxies
		for (task = MasterTask + 1; task < MasterTask + All.NTaskPerUniverse; task++)
			MPI_Ssend(&ngalaxies, 1, MPI_INT, task, TAG_NHALOS, MPI_COMM_WORLD);
	}

	//Check if we have enough memory left - otherwise return -1
	if ((long long)((ngalaxies+Nmax) * sizeof(struct kd_pos)) > (long long)(FreeBytes)) return -1;

	//Allocate the array holding the positions of all galaxies in this mass bin and a helper array to send data
	Galpos  = emalloc("KDPositions", ngalaxies * sizeof(struct kd_pos));
	if (Nmax > 0) possend = emalloc("PosSend", Nmax * sizeof(struct kd_pos));

	//The master task first writes all local galaxies that are in this bin to the position array
	if (ThisTask == MasterTask)
	{ //Go through the local galaxy array
		for (i = 0, j = 0; i < ngaltot; i++)
		{ //If the galaxy has the right stellar mass
			if (localgal[i].mass >= WpSet[iwpbin].min && localgal[i].mass < WpSet[iwpbin].max)
			{ //Write it to the position array and increment the counter
				Galpos[j].x[0] = localgal[i].x[0];
				Galpos[j].x[1] = localgal[i].x[1];
				Galpos[j].x[2] = localgal[i].x[2];
				j++;
			}
		}
		// Initialise the offset in the position array to the number of galaxies that have just been copied
		offset = nthistask;
	}

	//Each worker task writes all local galaxies that are in this bin to the helper array and sends it to the master
	if (ThisTask != MasterTask)
	{ //Go through the local galaxy array
		for (i = 0, j = 0; i < ngaltot; i++)
		{ //If the galaxy has the right stellar mass
			if (localgal[i].mass >= WpSet[iwpbin].min && localgal[i].mass < WpSet[iwpbin].max && Nmax > 0)
			{ //Write it to the helper array and increment the counter
				possend[j].x[0] = localgal[i].x[0];
				possend[j].x[1] = localgal[i].x[1];
				possend[j].x[2] = localgal[i].x[2];
				j++;
			}
		}
		//Send the number of galaxies in the helper array and then the helper array to the master task
		MPI_Ssend(&nthistask, 1, MPI_INT, MasterTask, TAG_NHALOS, MPI_COMM_WORLD);
		MPI_Ssend(possend, nthistask * sizeof(struct kd_pos), MPI_BYTE, MasterTask, TAG_HDATA, MPI_COMM_WORLD);
	}
	else
	{ //The master task goes through all worker tasks and receives the number of galaxies and the helper array
		for (task = MasterTask + 1; task < MasterTask + All.NTaskPerUniverse; task++)
		{ //Receive the number of galaxies and the helper array
			MPI_Recv(&nrecv, 1, MPI_INT, task, TAG_NHALOS, MPI_COMM_WORLD, &status);
			MPI_Recv(possend, nrecv * sizeof(struct kd_pos), MPI_BYTE, task, TAG_HDATA, MPI_COMM_WORLD, &status);
			//Write all positions in the helper array to the global positions array (using the right offset)
			for (i = offset, j = 0; j < nrecv; i++, j++) Galpos[i] = possend[j];
			//Increment the offset by the number of galaxies just received and copied
			offset += nrecv;
		}
	}

	//Free the helper array
	if (Nmax > 0) efree(possend);

	//Now all worker arrays receive the total global positions array from the master task
	if (ThisTask != MasterTask)
		MPI_Recv(Galpos, ngalaxies * sizeof(struct kd_pos), MPI_BYTE, MasterTask, TAG_HDATA, MPI_COMM_WORLD, &status);
	else
	{ //The master task goes through all workers and sends them the total global positions array
		for (task = MasterTask + 1; task < MasterTask + All.NTaskPerUniverse; task++)
			MPI_Ssend(Galpos, ngalaxies * sizeof(struct kd_pos), MPI_BYTE, task, TAG_HDATA, MPI_COMM_WORLD);
	}

	//Return the total number of galaxies in this bin
	return ngalaxies;
}




/*! \brief This function constructs a kd-tree from an array of positions
 *
 *  The function first allocates memory for the worst-case scenario, which is a maximally imbalanced tree.
 *  Each node has two integers, i0 and i1, which point either to the child nodes, or if the node is a leaf
 *  to the index of the first and last object in the position array. The initial count for each node is set
 *  to zero. Using the helper structure kd_build_task that holds the information about the objects on each node
 *  the nodes are split. For each node the bounding box is determined, and from this the widest dimension is
 *  found. The node is then split in this dimension at the centre. The splitting is stopped with a leaf if
 *  either the node has left fewer than #WP_NLEAF_MIN objects, or the widest dimension is smaller than
 *  #WP_NODE_WIDTH_MIN times the box length. Once the tree has been built, the minimum starting level (where
 *  there is at least one node per task), the maximum starting level (up to where nodes exist), and the mean
 *  level are calculated. Finally it is recursively determined how many nodes there are up to the starting level.
 *
 *
 *  \param *tree The kd-tree that is being built
 *  \param *pos The array of positions that is used to build the kd-tree
 *  \param count The number of objects that the tree contains
 *  \param dim The dimension of the tree
 */
int construct_tree(struct kd_tree *tree, struct kd_pos *pos, int count, int dim)
{

	int i, j, i0, i1, inode, n, widestdim, result;
	float maxwidth;
	struct kd_build_task task, task0, task1;

	//If there are no objects return zero
	if (count == 0) return 0;

	//Set the tree's dimensionality and number of objects
	tree->dim   = dim;
	tree->count = count;

	//Set the maximum number of nodes that will be needed if the tree is maximally imbalanced
	const int node_count_max = max(2 * (count - WP_NLEAF_MIN) + 1, 1);

	//Allocate all nodes that may be needed
	tree->nodes = emalloc("KD-TREE", (node_count_max) * sizeof(struct kd_node));

	//Set all initial node values
	for (i = 0; i < node_count_max; i++)
	{
		tree->nodes[i].i0    = -1;
		tree->nodes[i].i1    = -1;
		tree->nodes[i].leaf  = 0;
		tree->nodes[i].cnt   = 0;
	}

	//Set the starting node index (beyond the root)
	inode = 1;

	//Define the stack size used to create the tree (typically log2(count) is needed). Allow for tolerance.
	const int stack_size = (int)(max(log2((double)(count)) / (1.0 - 2.0 * ALLOC_TOLERANCE),1.0));

	//Declare and allocate the helper task array tasks with stack_size elements
	struct kd_build_task tasks[stack_size];
	//Initialise the helper tasks (use root values for all).
	for (i = 0; i < stack_size; i++)
	{
		tasks[i].first = 0;
		tasks[i].last  = count - 1;
		tasks[i].node  = 0;
	}

	//Initialise the current task to zero
	int current_task = 0;

	//Loop over all nodes until all nodes are leaves
	do
	{
		//To simplify the code set task to the current task and n to the current node index
		task = tasks[current_task];
		n    = tasks[current_task].node;

		//Set the initial bounding box limits
		for (j = 0; j < tree->dim; j++)
		{
			tree->nodes[n].min[j] = pos[task.first].x[j];
			tree->nodes[n].max[j] = pos[task.first].x[j];
		}

		//Determine the number of objects in this node
		tree->nodes[n].cnt = task.last - task.first + 1;

		//Go through all particles in the node and determine the bounding box
		for (i = task.first; i <= task.last; i++)
			for (j = 0; j < tree->dim; j++)
			{
				if (pos[i].x[j] <= tree->nodes[n].min[j]) tree->nodes[n].min[j] = pos[i].x[j];
				if (pos[i].x[j] >= tree->nodes[n].max[j]) tree->nodes[n].max[j] = pos[i].x[j];
			}

		//Initialise the default widest dimension as zero, and the maximum width as zero
		widestdim = 0;
		maxwidth  = 0.0;
		//Go through all dimensions
		for (j = 0; j < tree->dim; j++)
		{
			//If the width of this dimension is wider than the currently stored value
			if (tree->nodes[n].max[j] - tree->nodes[n].min[j] > maxwidth)
			{
				//Set the maximum width to the width in this dimension, and set the widest dimension to this one
				maxwidth = tree->nodes[n].max[j] - tree->nodes[n].min[j];
				widestdim = j;
			}
		}

		//Stop the splitting if we have a leaf, i.e. the number of particles is lower than WP_NLEAF_MIN
		//or the maximum width is smaller than WP_NODE_WIDTH_MIN times the box size
		if (task.last - task.first < WP_NLEAF_MIN || maxwidth < WP_NODE_WIDTH_MIN * All.Lbox)
		{
			//Set the indices i0 and i1 to the first and last object in the node, and set leaf to 1
			tree->nodes[n].i0   = task.first;
			tree->nodes[n].i1   = task.last;
			tree->nodes[n].leaf = 1;
			//As this task resulted in a leaf quit it and continue with the previous node/task
			current_task--;
			continue;
		}

		//Set the split value to half the bounding box size (in the widest dimension)
		const float splitvalue = (tree->nodes[n].min[widestdim] + tree->nodes[n].max[widestdim]) / 2.0;

		//Partition the object array such that there are only positions in the widest dimension that are
		//smaller than the split value in the first part and only positions that are larger in the second
		//part. The index of the first particle in the second part is given by k.
		int k = partition(pos, task.first, task.last, splitvalue, widestdim);

		// Set the index of the left child node to i0 and the index of the right child node to i1
		i0 = inode;
		i1 = i0 + 1;
		tree->nodes[n].i0  = i0;
		tree->nodes[n].i1  = i1;
		//Set the values for the left child node (i0) from first to k-1
		task0.first = task.first;
		task0.last  = k - 1;
		task0.node  = i0;
		//Set the values for the right child node (i1) from k to last
		task1.first = k;
		task1.last  = task.last;
		task1.node  = i1;

		//Let the next task deal with the left child node (task0) and increment the node index and current task
		tasks[current_task] = task0;
		inode++;
		current_task++;
		//If the task index is larger than the stack size the tree was more imbalanced than assumed, so return zero
		if (current_task >= stack_size)
		{
			if (ThisTask == MasterTask)
				printf("KD-tree for universe %d cannot be built. Problem with the stack size!\n",MasterTask/All.NUniverses);
			return 0;
		}
		//Let the next task deal with the left child node (task0) and increment the node index
		tasks[current_task] = task1;
		inode++;

		//Loop until the current task number is negative, i.e. all tasks have been resolved
	} while (current_task != -1);

	//Set the number of nodes to the final node index and reallocate the nodes to the needed size
	const int node_count = inode;
	tree->nodes = erealloc(tree->nodes, node_count * sizeof(struct kd_node));
	tree->node_count = node_count;

	//Define the intitial minimum level as the log2 of the number of tasks per universe
	tree->lmin = (int)(log2(All.NTaskPerUniverse));
	//Define the intitial maximum level as the log2 of the number of objects in the tree
	tree->lmax = (int)(log2(tree->node_count)) - 2;
	//Define the starting level for the parallel pair counting as the intermediate level between min and max
	tree->midlevel = (int)(0.5 * (float)(tree->lmin) + 0.5 * (float)(tree->lmax));

	//Loop until the number of starting nodes is at least the number of tasks per universe
	do
	{
		//Initialise the counter
		i = 0;
		//Compute the number of nodes up the the starting level
		result = find_nodes(*(tree), &i, 0, tree->midlevel, 0, &i, 0);
		//If the starting nodes could not be determined return zero
		if (result == 0) return 0;
		//If this number is lower than the number of tasks per universe increase the staring level by and repeat
		if (i < All.NTaskPerUniverse) (tree->midlevel)++;
		//Loop until the number of nodes at the starting level is at least the number of tasks per universe
	} while (i < All.NTaskPerUniverse && tree->midlevel <= tree->lmax);

	//Store the number nodes at the starting level
	tree->nstartnodes = i;

	//If all went fine return 1
	return 1;
}


/*! \brief This function identifies all starting nodes for the parallel pair counting
 *
 *  This function recursively determines the number of nodes up to a specific level l, and if the flag
 *  save is set to 1, writes all node indices to the array nodeindex. It checks if the current node is
 *  either at the required level or the node is a leaf. In this case it increments the counter and if
 *  needed writes the index to nodeindex. Otherwise is continues one level further down.
 *
 *  \param tree kd-tree that is being used to count the pairs
 *  \param *nodeindex Array containing all starting nodes
 *  \param inode Index of the current node
 *  \param l Level at which the starting nodes are located
 *  \param level Level of the current node
 *  \param *counter Counter that stores the total number of nodes at the level l and the leaves above
 *  \param save Flag that indicates if the array containing the starting nodes will be filled (1) or not (0)
 */
int find_nodes(struct kd_tree tree, int *nodeindex, int inode, int l, int level, int *counter, int save)
{
	//The resulting values for the left and right nodes
	int result1, result2;

	//If we are at the desired level or hit a leaf before add this node to the array (if save=1) and increment the counter - then return 1 (success)
	if (level == l || tree.nodes[inode].leaf == 1)
	{
		//if (ThisTask==0) printf("%d %d %d\n",(*counter),inode,save);
		if (save) nodeindex[(*counter)] = inode;
		(*counter)++;
		return 1;
	}
	//If we are at a level before that we go down the left and right child nodes with level+1
	else
	{
		result1 = find_nodes(tree, nodeindex, tree.nodes[inode].i0, l, level+1, counter, save);
		result2 = find_nodes(tree, nodeindex, tree.nodes[inode].i1, l, level+1, counter, save);
		return result1 * result2;
	}

	//We should never be here - if yes return zero (failure)
	return 0;
}


/*! \brief This function partitions the arrays of positions around the pivot value
 *
 *  The position array pos is partitioned such that all positions along the dimension dim that are smaller
 *  than the specified pivot value are moved to the first half of the array, and all positions in that dimension
 *  that are larger are moved to the second half. The function returns the index of the first element of the
 *  second half, i.e. the element at this index and all following elements are larger than the pivot index.
 *  Both the first and the second halves are otherwise unordered.
 *
 *  \param *pos The array of positions that is used to build the kd-tree
 *  \param left First element of the array to do the partitioning
 *  \param right Last element of the array to do the partitioning
 *  \param pivotValue Pivot value used to partition the array (smaller elements are to the left, larger elements are to the right)
 *  \param dim The dimension of the position that is used for the partitioning
 */
int partition(struct kd_pos *pos, int left, int right, double pivotValue, int dim)
{
//This swap function exchanges two elements in the position array
#define SWAP(a,b) {struct kd_pos swap=(a); (a)=(b); (b)=swap;}
	int i;
	//We start at the left index
	int storeIndex = left;
	//Now we move through all elements of the array from left to right
	for (i = left; i <= right; i++)
	{
		//If the position along the chosen dimension is smaller than the pivot value
		if (pos[i].x[dim] < pivotValue)
		{
			//Swap the position at the store index with the current position and increment the store index
			SWAP(pos[storeIndex],pos[i]);
			storeIndex++;
		}
	}
	//Return the store index where all elements to the left are smaller than the pivot value and the element at the index and all elements on the right are larger
	return storeIndex;
//Undefine the swap function
#undef SWAP
}


/*! \brief This function computes the squared minimum distance between two nodes
 *
 *  First the mean position of both nodes along each dimension is determined. If the first node has a larger
 *  mean position, the distance is determined as the difference between the maximum value of the bounding box
 *  of node 1 and the minimum of the bounding box of node 2, or by the difference between the minimum value of
 *  the bounding box of node 1 and the maximum of the bounding box of node 2, whatever is smaller (periodic
 *  boundary conditions are taken into account). If the first node has a smaller mean position, the distance is
 *  determined as the difference between the maximum value of the bounding box of node 2 and the minimum of the
 *  bounding box of node 1, or by the difference between the minimum value of the bounding box of node 2 and the
 *  maximum of the bounding box of node 1, whatever is smaller. The squared distances are added for each dimension
 *  and returned.
 *
 *  \param node1 The first node in the distance computation
 *  \param node2 The second node in the distance computation
 *  \param dim The dimensionality of the kd-tree
 */
float mindist2_node_node(struct kd_node node1, struct kd_node node2, int dim)
{
	int i;
	float x1,x2,m1,m2,d;
	//Initialise the distance to zero
	d = 0.0;
	//Go through all dimensions
	for (i = 0; i < dim; i++)
	{
		//Compute the positions of the box centres
		m1 = 0.5 * node1.min[i] + 0.5 * node1.max[i];
		m2 = 0.5 * node2.min[i] + 0.5 * node2.max[i];
		//If the central position of the first box is larger than the second one
		if (m1 >= m2)
		{
			//Compute the distance between the max of node 1 and the min of node 2 (periodic boundaries)
			x1 = (node1.max[i] - node2.min[i]);
			if (x1 < 0.0) x1 = 0.0;
			if (x1 > 0.5 * All.Lbox) x1 = All.Lbox - x1;
			//Compute the distance between the min of node 1 and the max of node 2 (periodic boundaries)
			x2 = (node1.min[i] - node2.max[i]);
			if (x2 < 0.0) x2 = 0.0;
			if (x2 > 0.5 * All.Lbox) x2 = All.Lbox - x2;
		}
		//If the central position of the second box is larger than the first one
		else
		{
			//Compute the distance between the max of node 2 and the min of node 1 (periodic boundaries)
			x1 = (node2.max[i] - node1.min[i]);
			if (x1 < 0.0) x1 = 0.0;
			if (x1 > 0.5 * All.Lbox) x1 = All.Lbox - x1;
			//Compute the distance between the min of node 2 and the max of node 1 (periodic boundaries)
			x2 = (node2.min[i] - node1.max[i]);
			if (x2 < 0.0) x2 = 0.0;
			if (x2 > 0.5 * All.Lbox) x2 = All.Lbox - x2;
		}
		if (x1*x1 < x2*x2) d += sqr(x1);
		else d += sqr(x2);
	}
	return d;
}


/*! \brief This function computes the squared maximum distance between two nodes
 *
 *  First the mean position of both nodes along each dimension is determined. If the first node has a larger
 *  mean position, the distance is determined as the difference between the maximum value of the bounding box
 *  of node 1 and the minimum of the bounding box of node 2 (periodic boundary conditions are taken into account).
 *  If the first node has a smaller mean position, the distance is determined as the difference between the
 *  maximum value of the bounding box of node 2 and the minimum of the bounding box of node 1. The squared
 *  distances are added for each dimension and returned.
 *
 *  \param node1 The first node in the distance computation
 *  \param node2 The second node in the distance computation
 *  \param dim The dimensionality of the kd-tree
 */
float maxdist2_node_node(struct kd_node node1, struct kd_node node2, int dim)
{
	int i;
	float x,m1,m2,d;
	//Initialise the distance to zero
	d = 0.0;
	//Go through all dimensions
	for (i = 0; i < dim; i++)
	{
		//Compute the positions of the box centres
		m1 = 0.5 * node1.min[i] + 0.5 * node1.max[i];
		m2 = 0.5 * node2.min[i] + 0.5 * node2.max[i];
		//If the central position of the first box is larger than the second one
		if (m1 >= m2)
		{
			//If the distance between the centres is smaller than half the box size the maximum distance is the distance between the max of node 1 and the min of node 2
			if (m1 - m2 <= 0.5 * All.Lbox) x = node1.max[i] - node2.min[i];
			//Otherwise the maximum distance is box size minus the distance between the min of node 1 and the max of node 2
			else x = All.Lbox + node2.max[i] - node1.min[i];
		}
		//If the central position of the second box is larger than the first one
		else
		{
			//If the distance between the centres is smaller than half the box size the maximum distance is the distance between the max of node 2 and the min of node 1
			if (m2 - m1 <= 0.5 * All.Lbox) x = node2.max[i] - node1.min[i];
			//Otherwise the maximum distance is box size minus the distance between the min of node 2 and the max of node 1
			else x = All.Lbox + node1.max[i] - node2.min[i];
		}
		d += x*x;
	}
	return d;
}


/*! \brief This function computes the number of pairs between two nodes
 *
 *  The algorithm follows the method presented in <a href="https://arxiv.org/abs/astro-ph/0012333">
 *  Moore et al. (2001)</a>.
 *  The number of pairs between two nodes are counted efficiently by pruning nodes and subnodes if they are
 *  either completely inside or outside of the distance bin. If the minimum distance of the nodes is larger than
 *  the inner radius rmin and the maximum distance is smaller than the outer radius rmax the product of the
 *  counts of the nodes is returned. Otherwise, if both nodes are leaves the number of pairs is computed in a
 *  brute force way, by computing the distance of each pair and checking if it is in the distance range, in
 *  which case the pair counter is incremented. If one or more nodes are parent nodes and cannot be pruned
 *  the function is called recursively for the child nodes and the results are added up and returned.
 *
 *  \param tree The kd-tree that is used to compute the number of pairs
 *  \param node1 The first node in the pair counting
 *  \param node2 The second node in the pair counting
 *  \param *pos Position array of the objects in the tree
 *  \param rmin2 The squared minimum distance for the pair counting
 *  \param rmax2 The squared maximum distance for the pair counting
 */
int dualtree(struct kd_tree tree, struct kd_node node1, struct kd_node node2, struct kd_pos *pos, float rmin2, float rmax2)
{
	int i,j,k;
	float mindist2, maxdist2, dist2;
	int cleft, cright, cnt;

	//Compute the squared minimum distance between the two nodes
	mindist2 = mindist2_node_node(node1,node2,tree.dim);
	//Return zero if it is larger than the squared maximum distance (this and all subnodes are too far away)
	if (mindist2 > rmax2) return 0;

	//Compute the squared maximum distance between the two nodes
	maxdist2 = maxdist2_node_node(node1,node2,tree.dim);
	//Return zero if it is smaller than the squared minimum distance (this and all subnodes are too close)
	if (maxdist2 <= rmin2) return 0;

	//If the nodes are further away than the minimum distance but closer than the maximum distance the number of pairs is given by the product of the particle numbers in each node
	if (rmin2 <= mindist2 && mindist2 <= maxdist2 && maxdist2 < rmax2)
		return (node1.cnt * node2.cnt);

	//Both nodes are leaves
	if (node1.leaf == 1 && node2.leaf == 1)
	{
		//Initialise the pair counter to zero
		cnt = 0;
		//Loop through all objects in the first leaf
		for (i = node1.i0; i <= node1.i1; i++)
		{
			//Loop through all objects in the second leaf
			for (j = node2.i0; j <= node2.i1; j++)
			{
				//Compute the distance between the pair of objects (i in node 1 and j in node 2)
				dist2 = 0.0;
				for (k = 0; k < tree.dim; k++) dist2 += sqr(NEAREST(pos[i].x[k] - pos[j].x[k]));
				//If the distance is larger than the minimum and smaller than the maximum increment the pair counter
				if (dist2 > rmin2 && dist2 < rmax2) cnt++;
			}
		}
		//Return all computed pairs between the two leaves
		return cnt;
	}
	//If not both nodes are leaves
	else
	{
		//Initialise the number of pairs for the left and right child node of the larger leaf
		cleft = cright = 0;
		//If the first node is a leaf (i.e. the second node is a parent node)
		if (node1.leaf == 1)
		{
			//Compute the pair counts for the left and right child node of node 2 recursively
			cleft  = dualtree(tree,tree.nodes[node2.i0],node1,pos,rmin2,rmax2);
			cright = dualtree(tree,tree.nodes[node2.i1],node1,pos,rmin2,rmax2);
		}
		//If the second node is a leaf (i.e. the first node is a parent node)
		else if (node2.leaf == 1)
		{
			//Compute the pair counts for the left and right child node of node 1 recursively
			cleft  = dualtree(tree,tree.nodes[node1.i0],node2,pos,rmin2,rmax2);
			cright = dualtree(tree,tree.nodes[node1.i1],node2,pos,rmin2,rmax2);
		}
		//Both nodes are parent nodes
		else
		{
			//If node 1 has more objects go further down node 1, otherwise go further down node 2
			if (node1.cnt > node2.cnt)
			{
				//Compute the pair counts for the left and right child node of node 1 recursively
				cleft  = dualtree(tree,tree.nodes[node1.i0],node2,pos,rmin2,rmax2);
				cright = dualtree(tree,tree.nodes[node1.i1],node2,pos,rmin2,rmax2);
			}
			//Node 2 has more (or the same number of) objects than node 1
			else
			{
				//Compute the pair counts for the left and right child node of node 2 recursively
				cleft  = dualtree(tree,tree.nodes[node2.i0],node1,pos,rmin2,rmax2);
				cright = dualtree(tree,tree.nodes[node2.i1],node1,pos,rmin2,rmax2);
			}
		}
		//If an error has occured in the subnodes and a negative value has been returned return -1 as well
		if (cleft  < 0) return -1;
		if (cright < 0) return -1;
		//If everything went well return the sum of the left and right child nodes
		return cleft+cright;
	}
	//We should never arrive here, if we do return -1
	return -1;
}


/*! \brief This function computes the number of pairs in a given distance bin
 *
 *  This algorithm is a parallelisation of the kd-tree pair counter following the method presented by
 *  <a href="http://www.linuxclustersinstitute.org/conferences/archive/2008/PDF/Dolence_98279.pdf">
 *  Dolence & Brunner (2008)</a>. If the minimum level is larger than the maximum level, the starting level
 *  is smaller than the minimum level, the starting level is larger than the maximum level, or the number of
 *  tasks per universe is one, the pairs are not searched in parallel and dualtree is called from the root
 *  nodes. Otherwise the pairs are searched in parallel. First the starting nodes are determined by calling
 *  find_nodes with the flag save = 1. Then the master task processes request using processRequests, and
 *  gives starting nodes to worker tasks. If there are currently no requests for work and there are still
 *  starting nodes left, the master task computes pairs for the next node. If there are no more nodes left
 *  the master task answers all work requests with a terminate message. The worker tasks request work
 *  from the master task and once they receive a new node index they compute the pairs for this node. This
 *  is repeated until a terminate message has been sent by the master task. In the end all counts are
 *  added up on all tasks using a message ring. The total number of pairs is then returned.
 *
 *  \param tree The kd-tree that is used to compute the number of pairs
 *  \param *pos Position array of the objects in the tree
 *  \param rmin The minimum distance for the pair counting
 *  \param rmax The maximum distance for the pair counting
 */
int tree_pairs(struct kd_tree tree, struct kd_pos *pos, float rmin, float rmax)
{
	int i, listpos, thispos, nextpos, nterm, totcnt, cnt, cntrecv, left, right;
	float rmax2 = sqr(rmax);
	float rmin2 = sqr(rmin);
	MPI_Status status;
	MPI_Status statuses[2];
	MPI_Request request[2];

	//If the starting level (midlevel) is smaller than the minimum level or larger than the maximum level
	//Or if the minimum level is larger than the maximum level or if we just have one task
	//Return the pair count (not done in parallel) starting from the root node
	//Divide by two to avoid double counting
	if (tree.lmin > tree.lmax || tree.midlevel < tree.lmin || tree.midlevel > tree.lmax || All.NTaskPerUniverse == 1 || tree.nstartnodes <= 1)
		return dualtree(tree,tree.nodes[0],tree.nodes[0],pos,rmin2,rmax2)/2;

	//Set W as the number of nodes at this level
	const int W = tree.nstartnodes;

	//Declare and allocate the list of nodes that will be used as starting points for the parallel computation
	int nodelist[W];
	for (i = 0; i < W; i++) nodelist[i] = 0;

	//Initialise counter to zero and identify the starting nodes by going down the tree recursively
	i = 0;
	find_nodes(tree, nodelist, 0, tree.midlevel, 0, &i, 1);

	//Initialise the pair counter to zero
	cnt = 0;

	//Each master task now distributes the work and if it has time does some computation as well
	if (ThisTask == MasterTask)
	{
		listpos = All.NTaskPerUniverse - 1;
		thispos = 0;
		nterm   = 0;
		//Loop until all worker tasks have terminated
		while (nterm < All.NTaskPerUniverse - 1)
		{
			//Process all requests (if a worker has finished computing pairs for one node give it a new node)
			processRequests(&listpos,&nterm,W);
			//If there are still uncomputed nodes do some pair counting yourself
			if (listpos < W)
			{
				//If this is not a leaf compute the pairs for both its child nodes and process requests between
				if (tree.nodes[nodelist[thispos]].leaf == 0)
				{
					cnt += dualtree(tree,tree.nodes[tree.nodes[nodelist[thispos]].i0],tree.nodes[0],pos,rmin2,rmax2);
					processRequests(&listpos,&nterm,W);
					cnt += dualtree(tree,tree.nodes[tree.nodes[nodelist[thispos]].i1],tree.nodes[0],pos,rmin2,rmax2);
				}
				//If this is a leaf just compute the counts for this node
				else
				{
					cnt += dualtree(tree,tree.nodes[nodelist[thispos]],tree.nodes[0],pos,rmin2,rmax2);
				}
				//As we have computed pairs for a node increment the index in the node list
				listpos++;
				//Set the next node to be computed to the just incremented node index
				thispos = listpos;
			}//Done with computing the current node on the master task
		}//Done with looping (all worker tasks have terminated)
	}//Done with master task
	//If this is a worker task
	else
	{
		//Set the initial node to be computed
		nextpos = ThisTask - MasterTask;
		//Loop while the next node is positive (if terminated nextpos is -1)
		while (nextpos >= 0)
		{
			//Ask the master task for the next node (send this task's index to the master task)
			MPI_Isend(&ThisTask,1,MPI_INT,MasterTask,TAG_REQUEST,MPI_COMM_WORLD,&request[0]);
			//Set the node to be computed to the next one in the queue
			thispos = nextpos;
			//Receive new node that this worker is responsible for
			MPI_Irecv(&nextpos,1,MPI_INT,MasterTask,TAG_REQUEST,MPI_COMM_WORLD,&request[1]);
			//Compute the pairs for this node
			cnt += dualtree(tree,tree.nodes[nodelist[thispos]],tree.nodes[0],pos,rmin2,rmax2);
			//Wait for all communication to finish (so we can use the memory again)
			MPI_Waitall(2,request,statuses);
		}//Done looping (no more nodes to be computed)
	}//Done with the worker tasks

	//Initialise the total number of pairs over all tasks with the one on this task
	totcnt = cnt;

	//Repeat for all other tasks on this universe
	for (i = 1; i < All.NTaskPerUniverse; i++)
	{
		//Set the left and right tasks (all within their own universe)
		left  = ThisTask - 1 < MasterTask ? MasterTask + All.NTaskPerUniverse - 1 : ThisTask - 1;
		right = ThisTask + 1 > MasterTask + All.NTaskPerUniverse - 1 ? MasterTask : ThisTask + 1;
		//Send the count on this task to the task on the left and receive the count from the right
		MPI_Sendrecv(&cnt, 1, MPI_INT, left, TAG_COUNT, &cntrecv, 1, MPI_INT, right, TAG_COUNT, MPI_COMM_WORLD, &status);
		//Set the count on this task to what has been received from the right
		cnt = cntrecv;
		//Add this count to the total number of pairs
		totcnt += cnt;
	}

	//Return half the total number of pairs (to avoid double counting)
	return totcnt/2;
}


/*! \brief This function processes all requests for work from worker tasks for the pair counting
 *
 *  This function is called by the master task to process all messages and work requests. A worker task sends
 *  a message containing its task index. If a starting node is left the master task responds with the next
 *  worknode and increments the position in the node list. If there are no more nodes left the master task sends
 *  a terminate message. In the end the function checks if there are new work request and if so answers them
 *  as well until no more requests are present.
 *
 *  \param *listpos The array of starting nodes for the parallel pair counting
 *  \param *nterm Number of tasks that have already terminated
 *  \param W Number of starting nodes for the parallel pair counting
 */
void processRequests(int *listpos, int *nterm, int W)
{
	int isMessage,source,Terminate;
	MPI_Status status;
	//Set the value for the terminate message
	Terminate = -1;
	//Check if there is a message from a worker task waiting (i.e. a request for work)
	MPI_Iprobe(MPI_ANY_SOURCE,TAG_REQUEST,MPI_COMM_WORLD,&isMessage,&status);
	//If there is a message
	while(isMessage != 0)
	{
		//Process messages until queue is empty - receive the task index
		MPI_Recv(&source,1,MPI_INT,status.MPI_SOURCE,TAG_REQUEST,MPI_COMM_WORLD,&status);
		//If there is still a node available
		if((*listpos)+1 < W)
		{ //Increment the node index
			(*listpos)++;
			//Send the next node in the list to the worker
			MPI_Send(listpos,1,MPI_INT,status.MPI_SOURCE,TAG_REQUEST,MPI_COMM_WORLD);
		}
		else
		{ //Out of work, send a terminate message and increment the number of terminated workers
			MPI_Send(&Terminate,1,MPI_INT,status.MPI_SOURCE,TAG_REQUEST,MPI_COMM_WORLD);
			(*nterm)++;
		}
		//Check if there is a new message from a worker
		MPI_Iprobe(MPI_ANY_SOURCE,TAG_REQUEST,MPI_COMM_WORLD,&isMessage,&status);
	}
	return;
}
