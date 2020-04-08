///////////////////////////////////////////////////////////////////////////////////////////////////
// Emerge code - File read_trees.c                                                               //
// Parts of these functions have been adapted from the GADGET code developed by Volker Springel  //
///////////////////////////////////////////////////////////////////////////////////////////////////
///
/// \file read_trees.c
/// \brief Contains functions that read halo merger trees and distribute them on the processors
///
/// This file contains all functions that are needed to read the halo merger trees and distribute
/// them to the different processors. This includes a function that identifies the forests and
/// joins all trees in a given forest on one task, a function that distributes trees or forests
/// on the tasks such that the load is balanced, a function that identifies the time steps, a
/// function that creates orphans, a function that computes halo properties based on their history,
/// a function that identifies each halo's host halo, and a function that copies the trees to the
/// tasks containing all other universes.
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


#define SKIP  fread(&blksize,sizeof(int),1,ifp);  ///< Reads in a block with the size of an integer


/*! \brief This function reads all merger trees and distributes them on the tasks
 *
 *  This function reads all merger trees from several files and distributes them smoothly
 *  on NTaskPerUniverse tasks. Then the routines for setting up the trees are called.
 *
 *  \param fname Filename to be read
 */
int read_trees(char *fname)
{
	int i, rep, num_files, rest_files, ngroups, master, filenr, masterTask, lastTask, max_load;
	int *nhalosInTask, *ntreesInTask;
	char buf[NSTRING];
	unsigned long long NtreesLocal, NhalosLocal;

	All.TotNhalos = 0;
	All.TotNtrees = 0;

	//Get number of files
	num_files  = find_files(fname);

	//Repeat reading the headers of the files twice. In the first iteration, only the tree and halo numbers ending
	//up on each task are assembled, followed by memory allocation. In the second iteration, the data is actually read in.
	for (rep = 0; rep < 2; rep++)
	{

		Nhalos = 0;
		Ntrees = 0;

		//Set the files left to read equal to the number of files
		rest_files = num_files;

		//First do all files one per task
		while(rest_files > All.NTaskPerUniverse)
		{
			//Set the file name (here there are always multiple files)
			sprintf(buf, "%s.%d", fname, ThisTask + (rest_files - All.NTaskPerUniverse));

			//Number of groups and master task
			ngroups = All.NTaskPerUniverse / All.NumFilesInParallel;
			if((All.NTaskPerUniverse % All.NumFilesInParallel)) ngroups++;
			master = (ThisTask / ngroups) * ngroups;

			//Loop through all tasks in the first universe
			for(i = 0; i < ngroups; i++)
			{
				if ((ThisTask == (master + i)) && (ThisTask < All.NTaskPerUniverse))
				{
					//This processor's turn
					if (rep == 0) share_halo_number_in_file(buf, ThisTask, ThisTask);
					else read_tree_file(buf, ThisTask, ThisTask);
				}
				MPI_Barrier(MPI_COMM_WORLD);
			}

			//Now the number of files is reduced by the number files just read
			rest_files -= All.NTaskPerUniverse;

		}

		//Then do remaining files and distribute to tasks
		if(rest_files > 0)
		{
			//Since there are fewer files than read tasks we need to ditribute the tasks
			distribute_file(All.NTaskPerUniverse, rest_files, 0, 0, All.NTaskPerUniverse-1, &filenr, &masterTask, &lastTask);

			//Get file name either for one or multiple files
			if(num_files > 1) sprintf(buf, "%s.%d", fname, filenr);
			else sprintf(buf, "%s", fname);

			//Number of groups
			ngroups = rest_files / All.NumFilesInParallel;
			if((rest_files % All.NumFilesInParallel)) ngroups++;

			//Loop through all files left
			for(i = 0; i < ngroups; i++)
			{
				if ((filenr / All.NumFilesInParallel == i) && (ThisTask < All.NTaskPerUniverse))
				{
					//This processor's turn
					if (rep == 0) share_halo_number_in_file(buf, masterTask, lastTask);
					else read_tree_file(buf, masterTask, lastTask);
				}
				MPI_Barrier(MPI_COMM_WORLD);
			}
		}

		//Determine total number of trees and haloes and allocate memory
		if (rep == 0)
		{
			//Calculate total number of trees and haloes
			NtreesLocal = (unsigned long long)(Ntrees);
			NhalosLocal = (unsigned long long)(Nhalos);
			MPI_Allreduce(&NtreesLocal, &All.TotNtrees, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
			MPI_Allreduce(&NhalosLocal, &All.TotNhalos, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

			//Determine the maximum number of haloes per task including some tolerance
			MPI_Allreduce(&Nhalos, &max_load, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
			All.MaxHalos = max_load / (1.0 - 2 * ALLOC_TOLERANCE);

			//Determine the maximum number of trees per task including some tolerance
			MPI_Allreduce(&Ntrees, &max_load, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
			All.MaxTrees = max_load / (1.0 - 2 * ALLOC_TOLERANCE) * 2;

			//Allocate the number of haloes in each tree and the haloes
			NhalosInTree   = (int *) emalloc_movable(&NhalosInTree, "NhalosInTree", All.MaxTrees * sizeof(int));
			OffsetHalos    = (int *) emalloc_movable(&OffsetHalos,  "OffsetHalos",  All.MaxTrees * sizeof(int));
#ifdef COMPUTE_ICM
			NtreesInForest = (int *) emalloc_movable(&NtreesInForest, "NtreesInForest", All.MaxTrees * sizeof(int));
			NhalosInForest = (int *) emalloc_movable(&NhalosInForest, "NhalosInForest", All.MaxTrees * sizeof(int));
#endif
			H              = (struct halo *) emalloc_movable(&H, "Halos", All.MaxHalos * sizeof(struct halo));
			TreeIDInTree   = (IDType *) emalloc_movable(&TreeIDInTree, "TreeIDInTree", All.MaxTrees * sizeof(IDType));

			//Set all to zero
			memset(NhalosInTree, 0, All.MaxTrees * sizeof(int));
			memset(OffsetHalos,  0, All.MaxTrees * sizeof(int));
#ifdef COMPUTE_ICM
			memset(NtreesInForest, 0, All.MaxTrees * sizeof(int));
			memset(NhalosInForest, 0, All.MaxTrees * sizeof(int));
#endif
			memset(H, 0, All.MaxHalos * sizeof(struct halo));
			memset(TreeIDInTree, 0, All.MaxTrees * sizeof(IDType));

			//Print what is done...
			if (ThisTask == 0) printf("%s\n%s Reading %llu trees with %llu haloes on %d tasks...\n", All.fullline,All.startline,All.TotNtrees,All.TotNhalos,All.NumFilesInParallel);

			//Allocate the CommBuffer
			CommBuffer = emalloc_movable(&CommBuffer,"CommBuffer", All.BufferSize * 1024 * 1024);

			MPI_Barrier(MPI_COMM_WORLD);
		}

		if (rep == 1) efree_movable(CommBuffer);

	}

	//If we want to compute the ICM we need to collect all trees in forests
#ifdef COMPUTE_ICM
	find_forests();
#endif

	//Now we can free TreeIDInTree;
	efree_movable(TreeIDInTree);

	//Sort the trees or forests by their size, so that smaller trees/forests can be moved to other tasks
	sort_by_size();

	//Since the number of haloes is not yet the same on all tasks they are distributed more evenly
	distribute_trees(All.NTaskPerUniverse, (double)(Nhalos));

	//Report memory usage
	report_memory_usage(&HighMark, "Trees_Read");

	//Identify the scale factors and the number of timesteps
	get_timesteps();

	//Add orphans to the trees - each leaf continues to the last timestep
	add_orphans();

	//Determine total number of trees and haloes
	NtreesLocal = (unsigned long long)(Ntrees);
	NhalosLocal = (unsigned long long)(Nhalos);
	MPI_Allreduce(&NtreesLocal, &All.TotNtrees, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&NhalosLocal, &All.TotNhalos, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

	//Sort the trees or forests by their size, so that smaller trees/forests can be moved to other tasks
	sort_by_size();

	//Since we now have added orphans we need to re-distribute the trees more evenly again
	distribute_trees(All.NTaskPerUniverse, (double)(Nhalos));

#ifndef COMPUTE_ICM
	//Setup all haloes (progenitor/descendant connection) tree-by-tree
	setup_haloes_by_tree();
#else
	//Setup all haloes (progenitor/descendant connection) forest-by-forest
	setup_haloes_by_forest();
	//Find the parent of each halo in each forest
	find_parents();
#endif

	//Compute peak mass and accretion rate
	compute_halo_history();

	//For communication purposed
	ntreesInTask = emalloc("NtreesInTask", NTask * sizeof(int));
	nhalosInTask = emalloc("NhalosInTask", NTask * sizeof(int));

	//Gather the number of trees and haloes locally
	MPI_Gather(&Ntrees, 1, MPI_INT, ntreesInTask, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Gather(&Nhalos, 1, MPI_INT, nhalosInTask, 1, MPI_INT, 0, MPI_COMM_WORLD);

	//And print
	if (ThisTask == 0)
	{
		printf("%s\n",All.fullline);
		for (i=0; i<All.NTaskPerUniverse; i++) printf("%s Task %6d contains %8d trees and %12d haloes (%f%%).\n",All.startline,i,ntreesInTask[i],nhalosInTask[i],(double)(nhalosInTask[i])/(double)(All.TotNhalos));
	}

	//Free communication arrays
	efree(nhalosInTask);
	efree(ntreesInTask);

	//If we have more than one universe we copy all data to the other universes
	if (All.NUniverses > 1) copy_trees_to_other_universes(All.NTaskPerUniverse,All.NUniverses);

	return 1;

}


/*! \brief This function determines on how many files the halo merger trees are distributed.
 *
 *  The root task first checks if there is a file that has a file name equal to fname. In this
 *  case a number of files equal to 1 is returned. Otherwise it checks if there are files that
 *  have a file name equal to fname.number and counts how many files there are. This number is
 *  then returned.
 *
 *  \param fname Filename to be read
 */
int find_files(char *fname)
{
	FILE *fd;
	char buf[NSTRING];
	int num_files, loop;

	num_files = 0;
	loop      = 1;
	sprintf(buf, "%s", fname);

	if(ThisTask == 0)
	{
		if((fd = fopen(buf, "r"))) {
			fclose(fd);
			num_files = 1;
		}
		else
		{
			while (loop)
			{
				sprintf(buf, "%s.%d", fname,num_files);
				if((fd = fopen(buf, "r"))) num_files++;
				else loop = 0;
				if (loop > 0) fclose(fd);
			}
		}
	}

	MPI_Bcast(&num_files, 1, MPI_INT, 0, MPI_COMM_WORLD);

	return num_files;

}


/*! \brief This function assigns a certain number of tasks to each file.
 *
 *  These tasks are containing the content of that file after the trees have been read
 *  The number of tasks per file is as homogeneous as possible.
 *  The number of files may at most be equal to the number of tasks.
 *
 *  \param ntask Number of tasks on which the files will be distributed
 *  \param nfiles Number of files on which the trees are located
 *  \param firstfile Index of the first file
 *  \param firsttask Index of the first task
 *  \param lasttask Index of the last task
 *  \param filenr Contains the file number to which this task belongs
 *  \param master The number of the task responsible to read the file
 *  \param last Number of the last task belonging to the same file as this task
 */
void distribute_file(int ntask, int nfiles, int firstfile, int firsttask, int lasttask, int *filenr, int *master, int *last)
{
	int i, group;
	int tasks_per_file = ntask / nfiles;
	int tasks_left = ntask % nfiles;

	if(tasks_left == 0)
	{
		group = ThisTask / tasks_per_file;
		*master = group * tasks_per_file;
		*last = (group + 1) * tasks_per_file - 1;
		*filenr = group;
		return;
	}

	for(i = 0, *last = -1; i < nfiles; i++)
	{
		*master = *last + 1;
		*last = (i + 1) * ((double) ntask) / nfiles;
		if(*last >= ntask) *last = *last - 1;
		if(*last < *master) endrun("last < master");
		*filenr = i;
		if(i == nfiles - 1) *last = ntask - 1;
		if(ThisTask >= *master && ThisTask <= *last) return;
	}
}


/*! \brief This function reads the tree and halo numbers in the file fname and distributes them
 *  to tasks 'readTask' to 'lastTask'.
 *
 *  First the number of trees in each file is read and communicated to all tasks. Then the array
 *  containing the number of haloes for each tree is read and communicated. The trees are divided
 *  among tasks such that the number of haloes for each task is as close as possible. This information
 *  is stored in Ntrees and Nhalos for each task.
 *
 *  \param fname Filename to be read
 *  \param readTask Task responsible for reading the file fname
 *  \param lastTask Last Task which gets data contained in the file
 */
void share_halo_number_in_file(const char *fname, int readTask, int lastTask)
{
	int i, itree, ntask, task, n_in_task;
	int ntrees_in_file, totNhalos_in_file;
	int *nhalos_in_file;
	int blksize;
	MPI_Status status;
	FILE *ifp;
	char buf[500];

	//First let the readTask read the number of trees in this file and communicate
	if(ThisTask == readTask)
	{
		//Open the file
		if(!(ifp = fopen(fname, "r")))
		{
			sprintf(buf, "Can't open file `%s' for reading merger trees.", fname);
			endrun(buf);
		}

		//Read the number of trees in this file
		SKIP;
		fread(&ntrees_in_file, sizeof(int), 1, ifp);
		SKIP;

		//Check if this was an integer
		if(blksize != 4)
			endrun("incorrect file format\n");

		//Send this to all other tasks in the group
		for(task = readTask + 1; task <= lastTask; task++)
		{
			MPI_Ssend(&ntrees_in_file, 1, MPI_INT, task, TAG_NTREES, MPI_COMM_WORLD);
		}
	}
	//All other tasks receive the number of trees in this file here
	else
	{
		MPI_Recv(&ntrees_in_file, 1, MPI_INT, readTask, TAG_NTREES, MPI_COMM_WORLD, &status);
	}

	//Allocate nhalos_in_file
	if(!(nhalos_in_file = (int *) emalloc("Nhalos_in_file",ntrees_in_file*sizeof(int))))
		endrun("Failed to allocate memory for nhalos_in_file.");

	//Now read the nhalos_in_file array
	if(ThisTask == readTask)
	{
		//Read the array containing the number of haloes for each tree
		SKIP;
		fread(nhalos_in_file, sizeof(int), ntrees_in_file, ifp);
		SKIP;

		//Send this to all other tasks in the group
		for(task = readTask + 1; task <= lastTask; task++)
		{
			MPI_Ssend(nhalos_in_file, ntrees_in_file, MPI_INT, task, TAG_NHALOS, MPI_COMM_WORLD);
		}
	}
	//All other tasks receive nhalos_in_file here
	else
	{
		MPI_Recv(nhalos_in_file, ntrees_in_file, MPI_INT, readTask, TAG_NHALOS, MPI_COMM_WORLD, &status);
	}

	//Get the total number of tasks in this group
	ntask  = lastTask - readTask + 1;

	//Count total number of haloes in this file
	totNhalos_in_file = 0;
	for (i = 0; i < ntrees_in_file; i++) totNhalos_in_file += nhalos_in_file[i];

	//Initialise counters
	itree = n_in_task = 0;
	task = readTask;
	//Loop through all halos in file
	for (i = 0; i < ntrees_in_file; i++)
	{
		//Add number of haloes in this tree to current number
		n_in_task += nhalos_in_file[i];
		itree++;
		//If current number of haloes is larger than the equal fraction
		if (n_in_task > totNhalos_in_file / ntask)
		{
			//Add number of trees and haloes to this task
			if (ThisTask == task)
			{
				Ntrees += itree;
				Nhalos += n_in_task;
			}
			//Now do next task and reset counters
			task++;
			itree = 0;
			n_in_task = 0;
		}
		//Go to next tree
	}
	//Since the last task has fewer haloes than the equal fraction, do it seperately
	if (ThisTask == lastTask)
	{
		Ntrees += itree;
		Nhalos += n_in_task;
	}

	efree(nhalos_in_file);

	//Close the file again
	if(ThisTask == readTask) fclose(ifp);

}


/*! \brief This function reads the halos in the file fname and distributes them to tasks 'readTask' to 'lastTask'.
 *
 *  First the number of trees in each file is read and communicated to all tasks. Then the array
 *  containing the number of haloes for each tree is read and communicated. The trees are divided
 *  among tasks such that the number of haloes for each task is as close as possible. This information
 *  is stored in Ntrees and Nhalos for each task.
 *
 *  \param fname Filename to be read
 *  \param readTask Task responsible for reading the file fname
 *  \param lastTask Last Task which gets data contained in the file
 */
void read_tree_file(char *fname, int readTask, int lastTask)
{
	int i, itree, hc, ntask, task, n_in_task, offset;
	int ntrees_in_file, totNhalos_in_file;
	int *nhalos_in_file, *nhalo_for_this_task, *ntree_for_this_task;
	int blksize, maxlen;
	IDType *treeID_in_file;
	MPI_Status status;
	FILE *ifp;
	char buf[500];

	//First let the readTask read the number of trees in this file and communicate
	if(ThisTask == readTask)
	{
		//Open the file
		if(!(ifp = fopen(fname, "r")))
		{
			sprintf(buf, "Can't open file `%s' for reading merger trees.", fname);
			endrun(buf);
		}

		//Read the number of trees in this file
		SKIP;
		fread(&ntrees_in_file, sizeof(int), 1, ifp);
		SKIP;

		//Check if this was an integer
		if(blksize != 4)
			endrun("incorrect file format\n");

		//Send this to all other tasks in the group
		for(task = readTask + 1; task <= lastTask; task++)
		{
			MPI_Ssend(&ntrees_in_file, 1, MPI_INT, task, TAG_NTREES, MPI_COMM_WORLD);
		}
	}
	//All other tasks receive the number of trees in this file here
	else
	{
		MPI_Recv(&ntrees_in_file, 1, MPI_INT, readTask, TAG_NTREES, MPI_COMM_WORLD, &status);
	}

	//Allocate nhalos_in_file and treeID_in_file
	nhalos_in_file = emalloc("Nhalos_in_file",ntrees_in_file*sizeof(int));
	treeID_in_file = emalloc("TreeID_in_file",ntrees_in_file*sizeof(IDType));

	//Now read the nhalos_in_file array
	if(ThisTask == readTask)
	{
		//Read the array containing the number of haloes for each tree
		SKIP;
		fread(nhalos_in_file, sizeof(int), ntrees_in_file, ifp);
		SKIP;

		SKIP;
		fread(treeID_in_file, sizeof(IDType), ntrees_in_file, ifp);
		SKIP;

		//Send this to all other tasks in the group
		for(task = readTask + 1; task <= lastTask; task++)
		{
			MPI_Ssend(nhalos_in_file, ntrees_in_file, MPI_INT, task, TAG_NHALOS, MPI_COMM_WORLD);
			MPI_Ssend(treeID_in_file, ntrees_in_file * sizeof(IDType), MPI_BYTE, task, TAG_TREEID, MPI_COMM_WORLD);
		}
	}
	//All other tasks receive nhalos_in_file here
	else
	{
		MPI_Recv(nhalos_in_file, ntrees_in_file, MPI_INT, readTask, TAG_NHALOS, MPI_COMM_WORLD, &status);
		MPI_Recv(treeID_in_file, ntrees_in_file * sizeof(IDType), MPI_BYTE, readTask, TAG_TREEID, MPI_COMM_WORLD, &status);
	}

	//Get the total number of tasks in this group
	ntask  = lastTask - readTask + 1;

	//Allocate ntree_for_this_task
	ntree_for_this_task = emalloc("ntree_for_this_task",ntask*sizeof(int));

	//Allocate nhalo_for_this_task
	nhalo_for_this_task = (int *) emalloc("nhalo_for_this_task",ntask*sizeof(int));

	//Count total number of haloes in this file
	totNhalos_in_file = 0;
	for (i = 0; i < ntrees_in_file; i++) totNhalos_in_file += nhalos_in_file[i];

	//Initialise counters
	task = itree = n_in_task = 0;
	//Loop through all halos in file
	for (i = 0; i < ntrees_in_file; i++)
	{
		//Add number of haloes in this tree to current number
		n_in_task += nhalos_in_file[i];
		itree++;
		//If current number of haloes is larger than the equal fraction
		if (n_in_task > totNhalos_in_file / ntask)
		{
			//Add number of trees and haloes to task
			nhalo_for_this_task[task] = n_in_task;
			ntree_for_this_task[task] = itree;
			//Now do next task and reset counters
			task++;
			itree = 0;
			n_in_task = 0;
		}
		//Go to next tree
	}
	//Since the last task has fewer haloes than the equal fraction, do it seperately
	nhalo_for_this_task[ntask-1] = n_in_task;
	ntree_for_this_task[ntask-1] = itree;

	//Print the file details
	if(ThisTask == readTask)
	{
		printf("%s Reading file `%s' on task=%d (contains %d trees and %d haloes).\n"
		       "%s Distributing this file to tasks %d-%d\n",
		       All.startline, fname, ThisTask, ntrees_in_file, totNhalos_in_file, All.startline, readTask, lastTask);
		fflush(stdout);
	}

	//Move NhalosInTree and H by ntree_for_this_task so there is exactly the right amount of space to hold the new data
	memmove(&NhalosInTree[ntree_for_this_task[ThisTask-readTask]], &NhalosInTree[0], Ntrees * sizeof(int));
	memmove(&TreeIDInTree[ntree_for_this_task[ThisTask-readTask]], &TreeIDInTree[0], Ntrees * sizeof(IDType));
	memmove(&H[nhalo_for_this_task[ThisTask-readTask]], &H[0], Nhalos * sizeof(struct halo));

	//Get the maximum number of haloes the CommBuffer can hold
	maxlen = ((int) (All.BufferSize * 1024 * 1024)) / (sizeof(struct haloread));
	itree  = 0;
	offset = 0;

	//Start reading the haloes
	if(ThisTask == readTask) SKIP;

	//Loop over tasks in this group
	for(task = readTask; task <= lastTask; task++)
	{
		//Check if this task has enough memory to store the new haloes
		if(task == ThisTask)
			if(Nhalos + nhalo_for_this_task[task-readTask] > All.MaxHalos)
			{
				sprintf(buf, "too many haloes. %d %d %d\n", Nhalos, nhalo_for_this_task[task-readTask], All.MaxHalos);
				endrun(buf);
			}

		//Copy the number of haloes for each tree for each task
		for (i=0; i<ntree_for_this_task[task-readTask]; i++,itree++)
			if (ThisTask == task)
			{
				NhalosInTree[i] = nhalos_in_file[itree];
				TreeIDInTree[i] = treeID_in_file[itree];
			}

		//Loop until all haloes have been sent
		do
		{
			//Number of haloes to send is the minimum of all haloes and the maximum number that fits onto the Buffer
			hc = nhalo_for_this_task[task-readTask];
			if(hc > maxlen) hc = maxlen;

			//If it is the reading task read hc haloes from the file and write it to the buffer
			if(ThisTask == readTask) fread(CommBuffer, sizeof(struct haloread), hc, ifp);

			//If it is the reading task and there are haloes to read send the CommBuffer to other task
			if(ThisTask == readTask && task != readTask && hc > 0)
				MPI_Ssend(CommBuffer, (sizeof(struct haloread)) * hc, MPI_BYTE, task, TAG_HDATA, MPI_COMM_WORLD);

			//If it is not the reading task and there are haloes to read receive the CommBuffer from the reading task
			if(ThisTask != readTask && task == ThisTask && hc > 0)
				MPI_Recv(CommBuffer, (sizeof(struct haloread)) * hc, MPI_BYTE, readTask, TAG_HDATA, MPI_COMM_WORLD, &status);

			//If this is the task where the haloes should end
			if(ThisTask == task)
			{
				//Read out the haloes from the CommBuffer and increment the offset by the number of haloes just read
				empty_read_buffer(offset, hc);
				offset += hc;
			}

			//Remove the number of haloes just read from the number of haloes left to read
			nhalo_for_this_task[task-readTask] -= hc;
		}
		//Loop until all haloes are read
		while(nhalo_for_this_task[task-readTask] > 0);
	}

	//When every halo has been read finish reading the file
	if(ThisTask == readTask) SKIP;

	//Initialise counters
	itree = n_in_task = 0;
	task = readTask;
	//Loop through all halos in file
	for (i = 0; i < ntrees_in_file; i++)
	{
		//Add number of haloes in this tree to current number
		n_in_task += nhalos_in_file[i];
		itree++;
		//If current number of haloes is larger than the equal fraction
		if (n_in_task > totNhalos_in_file / ntask)
		{
			//Add number of trees and haloes to this task
			if (ThisTask == task)
			{
				Ntrees += itree;
				Nhalos += n_in_task;
			}
			//Now do next task and reset counters
			task++;
			itree = 0;
			n_in_task = 0;
		}
		//Go to next tree
	}
	//Since the last task has fewer haloes than the equal fraction, do it seperately
	if (ThisTask == lastTask)
	{
		Ntrees += itree;
		Nhalos += n_in_task;
	}

	efree(nhalo_for_this_task);
	efree(ntree_for_this_task);
	efree(treeID_in_file);
	efree(nhalos_in_file);

	//Close the file again
	if(ThisTask == readTask) fclose(ifp);

}


/*! \brief This function reads out the CommBuffer that was filled with merger tree data.
 *
 *  The CommBuffer is cast into the halo structure and hc haloes are read from it.
 *  Masses are converted from Msol/h to code units (typically 1e9Msol).
 *  Positions are converted from Mpc/h to code units (typically Mpc).
 *  Radii are converted from Mpc/h to 1/1000 code units (typically kpc).
 *
 *  \param offset Halo corresponding to the first element in the CommBuffer
 *  \param hc Number of haloes in the CommBuffer
 */
void empty_read_buffer(int offset, int hc)
{
	int n;
	struct haloread *h;

	//Cast CommBuffer into halo structure
	h = (struct haloread *) CommBuffer;

	//Copy values for each halo to H
#if defined(_OPENMP) && (OPENMPTHREADS > 1)
#pragma omp parallel for schedule(static,10)
#endif
	for(n = 0; n < hc; n++)
	{ //In case the positions are out of bounds put all haloes back into the box
		if (h[n].pos[0] >= (float)(All.Lbox)) h[n].pos[0] = fmod(h[n].pos[0], (float)(All.Lbox));
		if (h[n].pos[1] >= (float)(All.Lbox)) h[n].pos[1] = fmod(h[n].pos[1], (float)(All.Lbox));
		if (h[n].pos[2] >= (float)(All.Lbox)) h[n].pos[2] = fmod(h[n].pos[2], (float)(All.Lbox));
		if (h[n].pos[0] < 0) h[n].pos[0] = fmod(h[n].pos[0], (float)(All.Lbox)) + (float)(All.Lbox);
		if (h[n].pos[1] < 0) h[n].pos[1] = fmod(h[n].pos[1], (float)(All.Lbox)) + (float)(All.Lbox);
		if (h[n].pos[2] < 0) h[n].pos[2] = fmod(h[n].pos[2], (float)(All.Lbox)) + (float)(All.Lbox);

		H[offset + n].haloid = h[n].haloid;
		H[offset + n].descid = h[n].descid;
		H[offset + n].upid   = h[n].upid;

		H[offset + n].np     = h[n].np;
		H[offset + n].mmp    = h[n].mmp;

		H[offset + n].a      = h[n].a;
		H[offset + n].mvir   = h[n].mvir/All.h_100/All.m_unit;;
		H[offset + n].rvir   = h[n].rvir/All.h_100/All.x_unit/1000.0*h[n].a;
		H[offset + n].c      = h[n].c;
		H[offset + n].lambda = h[n].lambda;

		H[offset + n].pos[0] = h[n].pos[0]/All.h_100/All.x_unit;;
		H[offset + n].pos[1] = h[n].pos[1]/All.h_100/All.x_unit;;
		H[offset + n].pos[2] = h[n].pos[2]/All.h_100/All.x_unit;;
		H[offset + n].vel[0] = h[n].vel[0];
		H[offset + n].vel[1] = h[n].vel[1];
		H[offset + n].vel[2] = h[n].vel[2];

		H[offset + n].idesc    = -1;
		H[offset + n].iprog    = -1;
		H[offset + n].icoprog  = -1;
	}

}


#ifdef COMPUTE_ICM
/*! \brief This function identifies all forests and collects them on one task.
 *
 *  First the forest file is read and the forest IDs belonging to each tree and halo are stored.
 *  Then the function checks which halo needs to go to which task to collect all trees/haloes in
 *  a forest and stores this info (temporarily in imain). The halo array is then sorted according
 *  to where the haloes need to go. Going through all task pairs the trees and haloes are exchanged.
 *  Finally the halo array is sorted according to the forest ID (and then according to the tree ID),
 *  and the number of trees and haloes in each forest is computed.
 */
void find_forests(void)
{
	int i, j, read, iforest, itree, ihalo, isend, irecv, maxlen, send_this_turn, left_to_send;
	int nids, ntreesInTask, nhalosInTask, ntrees_send, nhalos_send;
	int *ntreesToTask, *nhalosToTask, *ntreesOffset, *nhalosOffset, *ip;
	unsigned long long NtreesLocal, NhalosLocal, trees_to_read;
	IDType *forestID, *treeID, *forestIDTree;

	struct search_list_IDs *id_list;
	struct halo *hp;
	MPI_Status status;

	FILE *ifp;
	const int linesize = 1000;
	char line[linesize], buf[500], forestfilename[NSTRING];

	//Print what is done...
	if (ThisTask == 0) printf("%s\n%s Identifying the forests...\n", All.fullline,All.startline);

	//Allocate arrays
	nids = ((int) (All.BufferSize * 1024 * 1024 / 10)) / (sizeof(IDType));
	id_list      = emalloc_movable(&id_list,      "id_list",      All.MaxTrees * sizeof(struct search_list_IDs));
	forestIDTree = emalloc_movable(&forestIDTree, "ForestIDTree", All.MaxTrees * sizeof(IDType));
	forestID     = emalloc("ForestID",     nids * sizeof(IDType));
	treeID       = emalloc("TreeID",       nids * sizeof(IDType));

	//Write all tree IDs in this task in id_list
#if defined(_OPENMP) && (OPENMPTHREADS > 1)
#pragma omp parallel for
#endif
	for (i = 0; i < Ntrees; i++)
	{
		id_list[i].ID       = TreeIDInTree[i];
		id_list[i].position = i;
	}

	//Sort all tree IDs on this task
	qsort(id_list, Ntrees, sizeof(struct search_list_IDs), compare_id);

	//Get the forest file name
	sprintf(forestfilename, "%s.forests", All.treefile_name);

	if (ThisTask==0)
	{
		//Open the file on task 0
		if(!(ifp = fopen(forestfilename, "r")))
		{
			sprintf(buf, "Can't open file `%s' for reading the forests.", forestfilename);
			endrun(buf);
		}
	}

	//Read the forest file
	i             = 0;
	read          = 1;
	trees_to_read = All.TotNtrees;

	if (ThisTask==0) fgets(line,linesize,ifp); //1st line is a header

	//Go through each line
	do
	{
		//Only task 0 reads
		if (ThisTask==0)
		{
			//Read line
			if (fgets(line,linesize,ifp))
			{
				//Write to treeID and forestID
#ifndef LONGIDS
				sscanf(line,"%u %u\n",&treeID[i],&forestID[i]);
#else
				sscanf(line,"%llu %llu\n",&treeID[i],&forestID[i]);
#endif
				read = 1;
			}
			else read = 0;
		}

		//Increment the index
		if (ThisTask == 0)
		{
			//If this is the reading task (0) and we still read increment i by 1
			if (read == 1) i++;
		}
		else
		{
			//If this is not the reading task and we have more trees to read than one set increment by nids
			if (trees_to_read >= (unsigned long long)(nids)) i += nids;
			//Otherwise increment by the number of trees left to read and set read to 0
			else
			{
				i += (int)(trees_to_read);
				read = 0;
			}
		}

		//If one set has been read link the ids
		if (i == nids)
		{
			//If we still need to read then tell this to the other tasks
			MPI_Bcast(&read, 1, MPI_INT, 0, MPI_COMM_WORLD);

			//Send this set to all other tasks
#ifndef LONGIDS
			MPI_Bcast(treeID,   nids, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
			MPI_Bcast(forestID, nids, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
#else
			MPI_Bcast(treeID,   nids, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
			MPI_Bcast(forestID, nids, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
#endif
			//Each task searches for each tree ID in the set and if found sets its forest ID
			for (j = 0; j < nids; j++)
			{
				itree = binary_search_id(id_list, Ntrees, treeID[j]);
				if (itree >= 0) forestIDTree[itree] = forestID[j];
			}
			//Reset the index so we can fill the buffer again
			i = 0;
			//Decrease the number of trees left to read by nids
			trees_to_read -= (unsigned long long)(nids);
		}
		//Continute reading
	} while (read);

	//set the remaining IDs if there are some left
	if (i > 0)
	{
		//Send this set to all other tasks
#ifndef LONGIDS
		MPI_Bcast(treeID,   i, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
		MPI_Bcast(forestID, i, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
#else
		MPI_Bcast(treeID,   i, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
		MPI_Bcast(forestID, i, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
#endif

		//Each task searches for each tree ID in the set and if found sets its forestID
		for (j = 0; j < i; j++)
		{
			itree = binary_search_id(id_list, Ntrees, treeID[j]);
			if (itree >= 0) forestIDTree[itree] = forestID[j];
		}
	}

	//Free the read arrays
	efree(treeID);
	efree(forestID);

	//Set the forestID of each halo
	for (itree = 0, ihalo = 0; itree < Ntrees; itree++)
		for (i = 0; i < NhalosInTree[itree]; i++, ihalo++)
			H[ihalo].forestid = forestIDTree[itree];

	//Now join all forests by moving them from task to task

	//In order to sort the haloes easily we misuse iscale as original index and idesc as number of haloes in tree
	for (itree = 0, ihalo = 0; itree < Ntrees; itree++)
		for (i = 0; i < NhalosInTree[itree]; i++, ihalo++)
		{
			H[ihalo].iscale = ihalo;
			H[ihalo].idesc  = NhalosInTree[itree];
		}

	//Allocate the arrays that store how many trees/haloes to move to other tasks and where they are
	ntreesToTask = emalloc_movable(&ntreesToTask,"NtreesToTask", NTask * sizeof(int));
	nhalosToTask = emalloc_movable(&nhalosToTask,"NhalosToTask", NTask * sizeof(int));
	ntreesOffset = emalloc_movable(&ntreesOffset,"NtreesOffset", NTask * sizeof(int));
	nhalosOffset = emalloc_movable(&nhalosOffset,"NhalosOffset", NTask * sizeof(int));
	//Initialise as zero
	for (i = 0; i < NTask; i++)
	{
		ntreesToTask[i] = 0;
		nhalosToTask[i] = 0;
	}

	//Allocate forestID array and set it to zero
	forestID = emalloc_movable(&forestID,"ForestID", All.MaxTrees * sizeof(IDType));
	memset(forestID,0,All.MaxTrees * sizeof(IDType));
	//treeID   = emalloc("TreeID",   All.MaxTrees * sizeof(IDType));

	//Loop over all tasks that have to send trees
	for (isend = 0; isend < All.NTaskPerUniverse; isend++)
	{
		if (ThisTask == isend)
		{
			//Copy the forest IDs to forestID and set the number of haloes in this task
			for (i = 0; i < Ntrees; i++) forestID[i] = forestIDTree[i];
			ntreesInTask = Ntrees;
		}
		//Send the number of haloes in the send task to all tasks
		MPI_Bcast(&ntreesInTask, 1, MPI_INT, isend, MPI_COMM_WORLD);
		//Send the forestID array to all other tasks
#ifndef LONGIDS
		MPI_Bcast(forestID, ntreesInTask, MPI_UNSIGNED, isend, MPI_COMM_WORLD);
#else
		MPI_Bcast(forestID, ntreesInTask, MPI_UNSIGNED_LONG_LONG, isend, MPI_COMM_WORLD);
#endif

		//Go through all other tasks that can receive trees
		for (irecv = 0; irecv < All.NTaskPerUniverse; irecv++)
		{
			//Go through all haloes/trees on this task
			for (itree = 0, ihalo = 0; itree < ntreesInTask; itree++, ihalo+=H[ihalo].idesc)
			{
				//If this is the receiving task search the forest ID in the tree ID list
				if (ThisTask == irecv) iforest = binary_search_id(id_list, Ntrees, forestID[itree]);
				//Send this index to all other tasks
				MPI_Bcast(&iforest, 1, MPI_INT, irecv, MPI_COMM_WORLD);
				//If it has been found and this is the sending task...
				if (iforest >= 0 && ThisTask == isend)
				{
					//... increment the number of tree going to the receiving task
					ntreesToTask[irecv]++;
					//imain now stores the task where the tree needs to go
					for (j = ihalo; j < ihalo + H[ihalo].idesc; j++) H[j].imain = irecv;
				}
			}
		}
	}

	//Print what is done...
	if (ThisTask == 0) printf("%s Sorting haloes by the task to which they will be moved...\n",All.startline);

	//Sort all haloes according to which task they need to be send
	qsort(H, Nhalos, sizeof(struct halo), compare_totask);

	//As the sorting destroyed the old one - get the new NhalosInTrees array
	for (itree = 0, ihalo = 0; itree < Ntrees; itree++)
	{
		NhalosInTree[itree]  = H[ihalo].idesc;
		ihalo               += NhalosInTree[itree];
	}

	//Count how many haloes go to each task
	for (i = 0; i < Nhalos; i++) nhalosToTask[H[i].imain]++;

	//Compute the offsets that is the starting indices of each set of trees/haloes
	ntreesOffset[0] = 0;
	nhalosOffset[0] = 0;
	for (i = 1; i < NTask; i++)
	{
		ntreesOffset[i] = ntreesOffset[i-1] + ntreesToTask[i-1];
		nhalosOffset[i] = nhalosOffset[i-1] + nhalosToTask[i-1];
	}

	//Compute how many haloes are sent to other tasks (for printing)
	ntrees_send = 0;
	nhalos_send = 0;
	for (i = 0; i < All.NTaskPerUniverse; i++)
	{
		if (i != ThisTask)
		{
			ntrees_send += ntreesToTask[i];
			nhalos_send += nhalosToTask[i];
		}
	}

	//Print what is done...
	if (MasterTask == 0)
		printf("%s Task %6d is moving %6d trees with %8d haloes to other tasks.\n",All.startline,ThisTask,ntrees_send,nhalos_send);

	//Allocate the CommBuffer
	CommBuffer = emalloc_movable(&CommBuffer,"CommBuffer", All.BufferSize * 1024 * 1024);

	//Go through all pairs and exchange the trees and haloes
	for (isend = 1; isend < All.NTaskPerUniverse; isend++)
	{
		//irecv is always at least 1 smaller than isend
		for (irecv = 0; irecv < isend; irecv++)
		{
			//First send from isend to irecv

			//Get number of trees that will be sent and broadcast
			ntrees_send = ntreesToTask[irecv];
			MPI_Bcast(&ntrees_send, 1, MPI_INT, isend, MPI_COMM_WORLD);
			//Get number of trees that are already at the irecv task and broadcast
			ntreesInTask = Ntrees;
			MPI_Bcast(&ntreesInTask, 1, MPI_INT, irecv, MPI_COMM_WORLD);

			//Check if the receiving task has enough space for the new trees
			if(ntreesInTask + ntrees_send > All.MaxTrees)
			{
				//If we need more trees than currently allocated then reallocate
				All.MaxTrees = ntreesInTask + ntrees_send;
				NhalosInTree = (int *) erealloc_movable(NhalosInTree, All.MaxTrees * sizeof(int));
			}

			//The isend task sends all NhalosInTree to the irecv task
			if (ThisTask == isend && ntrees_send > 0)
			{
				ip = (int *) CommBuffer;
				//Write CommBuffer from NhalosInTree
				for (i = ntreesOffset[irecv]; i < ntreesOffset[irecv]+ntrees_send; i++)
					*ip++ = NhalosInTree[i];
				//If there are trees to send then send the CommBuffer to the receiving task
				MPI_Ssend(CommBuffer, (sizeof(int)) * ntrees_send, MPI_INT, irecv, TAG_NHALOS, MPI_COMM_WORLD);
				//Move the remainder of the tree array forward to fill the gap that has been sent
				memmove(&NhalosInTree[ntreesOffset[irecv]], &NhalosInTree[ntreesOffset[irecv]+ntrees_send], (Ntrees - ntreesOffset[irecv] - ntrees_send) * sizeof(int));
				//Remove the number of trees that have been sent from Ntrees and set ntreesToTask to zero
				Ntrees -= ntrees_send;
				ntreesToTask[irecv] = 0;
				//Remove the number of sent trees from the offsets for the remaining tasks
				for (i = irecv + 1; i < All.NTaskPerUniverse; i++) ntreesOffset[i] -= ntrees_send;
			}

			//The irecv task now receives all NhalosInTree
			if (ThisTask == irecv && ntrees_send > 0) {
				MPI_Recv(CommBuffer, (sizeof(int)) * ntrees_send, MPI_INT, isend, TAG_NHALOS, MPI_COMM_WORLD, &status);
				ip = (int *) CommBuffer;
				//Read CommBuffer and write to NhalosInTree
				for (j = 0, i = Ntrees; j < ntrees_send; i++,j++) NhalosInTree[i] = ip[j];
				//Add the number of trees that have been sent to Ntrees
				Ntrees += ntrees_send;
			}

			//Now the haloes from isend to irecv
			//Get number of haloes that will be sent and broadcast
			nhalos_send = nhalosToTask[irecv];
			MPI_Bcast(&nhalos_send, 1, MPI_INT, isend, MPI_COMM_WORLD);
			//Get number of haloes that are already at the irecv task and broadcast
			nhalosInTask = Nhalos;
			MPI_Bcast(&nhalosInTask, 1, MPI_INT, irecv, MPI_COMM_WORLD);

			//Check if the receiving task has enough space for the new haloes
			if (nhalos_send + nhalosInTask > All.MaxHalos)
			{
				//If we can make space for the new trees then increase MaxHaloes and reallocate
				if ( (long long)((nhalos_send + nhalosInTask - All.MaxHalos) * sizeof(struct halo) / (1.0 - 2 * ALLOC_TOLERANCE)) < (long long)(FreeBytes))
				{
					All.MaxHalos += (int)(FreeBytes * (1.0 - 2 * ALLOC_TOLERANCE) / sizeof(struct halo));
					H  = (struct halo *) erealloc_movable(H, All.MaxHalos * sizeof(struct halo));
				}
				//Otherwise crash
				else endrun("Cannot allocate enough memory for the halos to move the forests.");
			}

			//Set the maximum number of haloes that can be fitted onto the CommBuffer
			maxlen       = ((int) (All.BufferSize * 1024 * 1024)) / (sizeof(struct halo));
			left_to_send = nhalos_send;

			//Loop until all haloes have been sent
			do
			{
				//Send this turn the minimum of all haloes left to send or the maximum number that fits onto the Buffer
				send_this_turn = left_to_send;
				if (send_this_turn > maxlen) send_this_turn = maxlen;

				//Sending task
				if (ThisTask == isend && send_this_turn > 0)
				{
					// Write send_this_turn haloes to CommBuffer
					hp = (struct halo *) CommBuffer;
					for (i = nhalosOffset[irecv] + nhalosToTask[irecv] - left_to_send; i < nhalosOffset[irecv] + nhalosToTask[irecv] - left_to_send + send_this_turn; i++) *hp++ = H[i];
					MPI_Ssend(CommBuffer, (sizeof(struct halo)) * send_this_turn, MPI_BYTE, irecv, TAG_HDATA, MPI_COMM_WORLD);
				}

				//Receiving task
				if(ThisTask == irecv && send_this_turn > 0)
				{
					//Receive CommBuffer from isend
					MPI_Recv(CommBuffer, (sizeof(struct halo)) * send_this_turn, MPI_BYTE, isend, TAG_HDATA, MPI_COMM_WORLD, &status);
					//Write send_this_turn haloes in CommBuffer to end of halo array
					hp = (struct halo *) CommBuffer;
					for (j = 0, i = Nhalos + nhalos_send - left_to_send; j < send_this_turn; i++,j++) H[i] = hp[j];
				}

				//Haloes left to send are reduced by what has been sent this turn
				left_to_send -= send_this_turn;
			}
			//Loop until all haloes are sent
			while (left_to_send > 0);

			//The isend task now cleans up
			if (ThisTask == isend && nhalos_send > 0)
			{
				//Move the remainder of the halo array forward to fill the gap that has been sent
				memmove(&H[nhalosOffset[irecv]], &H[nhalosOffset[irecv]+nhalos_send], (Nhalos - nhalosOffset[irecv] - nhalos_send) * sizeof(struct halo));
				//Remove the number of haloes that have been sent from Ntrees and set nhalosToTask to zero
				Nhalos -= nhalos_send;
				nhalosToTask[irecv] = 0;
				for (i = irecv + 1; i < All.NTaskPerUniverse; i++) nhalosOffset[i] -= nhalos_send;
				//Set the rest of the halo array to zero
				memset(&H[Nhalos], 0, (All.MaxHalos - Nhalos) * sizeof(struct halo));
			}

			//The irecv task adds the number of haloes that have been sent to Nhalos
			if (ThisTask == irecv && nhalos_send > 0) Nhalos += nhalos_send;

			//Now send back from irecv to isend
			//Get number of trees that will be sent and broadcast
			ntrees_send = ntreesToTask[isend];
			MPI_Bcast(&ntrees_send, 1, MPI_INT, irecv, MPI_COMM_WORLD);
			//Get number of trees that are already at the isend task and broadcast
			ntreesInTask = Ntrees;
			MPI_Bcast(&ntreesInTask, 1, MPI_INT, isend, MPI_COMM_WORLD);

			//Check if the receiving task has enough space for the new trees
			if(ntreesInTask + ntrees_send > All.MaxTrees)
			{
				//If we need more trees than currently allocated then reallocate
				All.MaxTrees = ntreesInTask + ntrees_send;
				NhalosInTree = (int *) erealloc_movable(NhalosInTree, All.MaxTrees * sizeof(int));
			}

			//The irecv task sends all NhalosInTree to the isend task
			if (ThisTask == irecv && ntrees_send > 0)
			{
				ip = (int *) CommBuffer;
				//Write CommBuffer from NhalosInTree
				for (i = ntreesOffset[isend]; i < ntreesOffset[isend]+ntrees_send; i++) *ip++ = NhalosInTree[i];
				//If there are trees to send then send the CommBuffer to the receiving task
				MPI_Ssend(CommBuffer, (sizeof(int)) * ntrees_send, MPI_INT, isend, TAG_NHALOS, MPI_COMM_WORLD);
				//Move the remainder of the tree array forward to fill the gap that has been sent
				memmove(&NhalosInTree[ntreesOffset[isend]], &NhalosInTree[ntreesOffset[isend]+ntrees_send], (Ntrees - ntreesOffset[isend] - ntrees_send) * sizeof(int));
				//Remove the number of trees that have been sent from Ntrees and set ntreesToTask to zero
				Ntrees -= ntrees_send;
				ntreesToTask[isend] = 0;
				//Remove the number of sent trees from the offsets for the remaining tasks
				for (i = isend + 1; i < All.NTaskPerUniverse; i++) ntreesOffset[i] -= ntrees_send;
			}

			//The isend task now receives all NhalosInTree
			if (ThisTask == isend && ntrees_send > 0) {
				MPI_Recv(CommBuffer, (sizeof(int)) * ntrees_send, MPI_INT, irecv, TAG_NHALOS, MPI_COMM_WORLD, &status);
				ip = (int *) CommBuffer;
				//Read CommBuffer and write to NhalosInTree
				for (j = 0, i = Ntrees; j < ntrees_send; i++,j++) NhalosInTree[i] = ip[j];
				//Add the number of trees that have been sent to Ntrees
				Ntrees += ntrees_send;
			}

			//Now the haloes from irecv to isend
			//Get number of haloes that will be sent and broadcast
			nhalos_send = nhalosToTask[isend];
			MPI_Bcast(&nhalos_send, 1, MPI_INT, irecv, MPI_COMM_WORLD);
			//Get number of haloes that are already at the irecv task and broadcast
			nhalosInTask = Nhalos;
			MPI_Bcast(&nhalosInTask, 1, MPI_INT, isend, MPI_COMM_WORLD);

			//Check if the receiving task has enough space for the new haloes
			if (nhalos_send + nhalosInTask > All.MaxHalos)
			{
				//If we can make space for the new trees then increase MaxHaloes and reallocate
				if ( (long long)((nhalos_send + nhalosInTask - All.MaxHalos) * sizeof(struct halo) / (1.0 - 2 * ALLOC_TOLERANCE)) < (long long)(FreeBytes))
				{
					All.MaxHalos += (int)(FreeBytes * (1.0 - 2 * ALLOC_TOLERANCE) / sizeof(struct halo));
					H  = (struct halo *) erealloc_movable(H, All.MaxHalos * sizeof(struct halo));
				}
				//Otherwise free all arrays and skip the re-distribution
				else endrun("Cannot allocate enough memory for the halos to move the forests.");
			}

			//Set the maximum number of haloes that can be fitted onto the CommBuffer
			maxlen       = ((int) (All.BufferSize * 1024 * 1024)) / (sizeof(struct halo));
			left_to_send = nhalos_send;

			//Loop until all haloes have been sent
			do
			{
				//Send this turn the minimum of all haloes left to send or the maximum number that fits onto the Buffer
				send_this_turn = left_to_send;
				if (send_this_turn > maxlen) send_this_turn = maxlen;

				//Sending task
				if (ThisTask == irecv && send_this_turn > 0)
				{
					// Write send_this_turn haloes to CommBuffer
					hp = (struct halo *) CommBuffer;
					for (i = nhalosOffset[isend] + nhalosToTask[isend] - left_to_send; i < nhalosOffset[isend] + nhalosToTask[isend] - left_to_send + send_this_turn; i++) *hp++ = H[i];
					MPI_Ssend(CommBuffer, (sizeof(struct halo)) * send_this_turn, MPI_BYTE, isend, TAG_HDATA, MPI_COMM_WORLD);
				}

				//Receiving task
				if(ThisTask == isend && send_this_turn > 0)
				{
					//Receive CommBuffer from taskSorted[isend]
					MPI_Recv(CommBuffer, (sizeof(struct halo)) * send_this_turn, MPI_BYTE, irecv, TAG_HDATA, MPI_COMM_WORLD, &status);
					//Write send_this_turn haloes in CommBuffer to end of halo array
					hp = (struct halo *) CommBuffer;
					for (j = 0, i = Nhalos + nhalos_send - left_to_send; j < send_this_turn; i++,j++) H[i] = hp[j];
				}

				//Haloes left to send are reduced by what has been sent this turn
				left_to_send -= send_this_turn;
			}
			//Loop until all haloes are sent
			while (left_to_send > 0);

			//The irecv task now cleans up
			if (ThisTask == irecv && nhalos_send > 0)
			{
				//Move the remainder of the halo array forward to fill the gap that has been sent
				memmove(&H[nhalosOffset[isend]], &H[nhalosOffset[isend]+nhalos_send], (Nhalos - nhalosOffset[isend] - nhalos_send) * sizeof(struct halo));
				//Remove the number of haloes that have been sent from Ntrees and set nhalosToTask to zero
				Nhalos -= nhalos_send;
				nhalosToTask[isend] = 0;
				for (i = isend + 1; i < All.NTaskPerUniverse; i++) nhalosOffset[i] -= nhalos_send;
				//Set the rest of the halo array to zero
				memset(&H[Nhalos], 0, (All.MaxHalos - Nhalos) * sizeof(struct halo));
			}

			//The isend task adds the number of haloes that have been sent to Nhalos
			if (ThisTask == isend && nhalos_send > 0) Nhalos += nhalos_send;

		}//End irecv
	}//End isend

	//Free all arrays that are no longer needed
	efree_movable(CommBuffer);
	efree_movable(forestID);
	efree_movable(nhalosOffset);
	efree_movable(ntreesOffset);
	efree_movable(nhalosToTask);
	efree_movable(ntreesToTask);
	efree_movable(forestIDTree);
	efree_movable(id_list);

	//Set the values of all haloes and trees beyond the limit to zero
	memset(&NhalosInTree[Ntrees], 0, (All.MaxTrees - Ntrees) * sizeof(int));
	memset(&H[Nhalos], 0, (All.MaxHalos - Nhalos) * sizeof(struct halo));

	//Compute the total number of trees and haloes
	NtreesLocal = (unsigned long long)(Ntrees);
	NhalosLocal = (unsigned long long)(Nhalos);
	MPI_Allreduce(&NtreesLocal, &All.TotNtrees, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&NhalosLocal, &All.TotNhalos, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

	//Copy the maximum number of trees to maxlen
	MPI_Allreduce(&Ntrees, &maxlen, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
	//If there are more trees than what we have now reallocate the NhalosInTree array
	if (maxlen > All.MaxTrees)
	{
		//Leave some tolerance
		All.MaxTrees = maxlen / (1.0 - 2 * ALLOC_TOLERANCE) * 2;
		NhalosInTree = (int *) erealloc_movable(NhalosInTree, All.MaxTrees * sizeof(int));
	}

	//For sorting we again misuse iscale as the original index
#if defined(_OPENMP) && (OPENMPTHREADS > 1)
#pragma omp parallel for schedule(static,100)
#endif
	for (i = 0; i < Nhalos; i++) H[i].iscale = i;

	//Print what is done...
	if (ThisTask == 0) printf("%s Sorting haloes by their forest ID...\n",All.startline);

	//Sort all haloes according to their forest ID
	qsort(H, Nhalos, sizeof(struct halo), compare_forest);

	//Get new NhalosInTrees array as the old one was destroyed again in the sorting
	for (itree = 0, ihalo = 0; itree < Ntrees; itree++)
	{
		NhalosInTree[itree]  = H[ihalo].idesc;
		ihalo               += NhalosInTree[itree];
	}

	//Initialise the number of forests on each task
	Nforests = 0;
	nhalosInTask = 0; //Used as haloes in this forest
	ntreesInTask = 0; //Used as trees in this forest
	//Go through all haloes and count the number of trees and haloes in each forest
	for (i = 1; i < Nhalos; i++)
	{
		//Increment the number of haloes in this forest
		nhalosInTask++;
		//If we are at the root then increment the number of trees in this forest
		if (H[i].descid == 0) ntreesInTask++;
		//If we arrive at a new forest
		if (H[i].forestid != H[i-1].forestid)
		{
			//copy the number of trees and haloes to this forest
			NhalosInForest[Nforests] = nhalosInTask;
			NtreesInForest[Nforests] = ntreesInTask;
			//Set the number of trees and haloes per forest back to zero
			nhalosInTask = 0;
			ntreesInTask = 0;
			//Go to new forest
			Nforests++;
		}
	}
	//For the last forest we need to copy the number of trees and haloes
	NhalosInForest[Nforests] = nhalosInTask + 1;
	NtreesInForest[Nforests] = ntreesInTask;
	//And increment the number of forests one last time
	Nforests++;

	//Reset iscale imain and idesc to the default values
#if defined(_OPENMP) && (OPENMPTHREADS > 1)
#pragma omp parallel for schedule(static,100)
#endif
	for (i = 0; i < Nhalos; i++)
	{
		H[i].iscale =  0;
		H[i].imain  =  i;
		H[i].idesc  = -1;
	}

	//Report Memory Usage
	HighMark = 0;
	report_memory_usage(&HighMark, "find_forests");

}
#endif


/*! \brief This function sorts the trees or forests by their size (largest go first)
 *
 *  Since trees and forests can get moved to other tasks to distribute the load more evenly it is best to order
 *  them according to their size. In this way the smallest trees or forests are in the end of the arrays and
 *  the number of haloes that are moved to other tasks can be set closer to the target number.
 *  First some extra information is added to each halo so they can be sorted. The original index is stored
 *  in icoprog, the size of the enitity that is being sorted (trees or forests) is stored in idescgal, and
 *  if the forests are being sorted the tree size is stored in imass. Then all haloes are sorted with
 *  compare_size according to their idescgal value and if tied according to their icoprog value.
 *  The arrays #NhalosInTree, #NtreesInForest, and #NhalosInForest are restored from the information given
 *  in the halo array. Finally the variables used for the sorting are restored (icoprog, idescgal and
 *  possibly imass).
 */
void sort_by_size(void)
{
	int i, itree, ihalo;
#ifdef COMPUTE_ICM
	int j, imass, iforest, nHInForest, nTInForest;
#endif

	//If we do not need the ICM we can go tree-by-tree
#ifndef COMPUTE_ICM
	//Loop through all trees
	for (itree = 0, ihalo = 0; itree < Ntrees; itree++)
	{
		//Go through all haloes in this tree
		for (i = 0; i < NhalosInTree[itree]; i++, ihalo++)
		{
			H[ihalo].icoprog  = ihalo;                    //Use icoprog as original index
			H[ihalo].idescgal = NhalosInTree[itree];      //Use idescgal as number of haloes per tree
		}
	}

	//Print what is done...
	if (ThisTask == 0) printf("%s\n%s Sorting trees by their size...\n",All.fullline,All.startline);

	//If we do need the ICM we must go forest-by-forest
#else
	//Loop through all forests
	for (iforest = 0, itree = 0, ihalo = 0; iforest < Nforests; iforest++)
	{
		//Loop through all trees
		for (i = 0; i < NtreesInForest[iforest]; i++, itree++)
		{
			//Go through all haloes in this tree
			for (j = 0; j < NhalosInTree[itree]; j++, ihalo++)
			{
				H[ihalo].icoprog  = ihalo;                    //Use icoprog as original index
				H[ihalo].idescgal = NhalosInForest[iforest];  //Use idescgal as number of haloes per forest
				H[ihalo].imass    = NhalosInTree[itree];      //Use imass as number of haloes per tree
			}
		}
	}

	//Print what is done...
	if (ThisTask == 0) printf("%s\n%s Sorting forests by their size...\n",All.fullline,All.startline);
#endif

	//Sort all haloes according to their tree or forest length
	qsort(H, Nhalos, sizeof(struct halo), compare_size);

	//Get new NhalosInTrees array as the old one was destroyed in the sorting
	for (itree = 0, ihalo = 0; itree < Ntrees; itree++)
	{
#ifndef COMPUTE_ICM
		//If we have sorted by trees NhalosInTree is stored in idescgal
		NhalosInTree[itree]  = H[ihalo].idescgal;
#else
		//If we have sorted by forests NhalosInTree is stored in imass
		NhalosInTree[itree]  = H[ihalo].imass;
#endif
		ihalo               += NhalosInTree[itree];
	}
#ifdef COMPUTE_ICM
	//Initialise the number of forests on each task
	iforest = 0;
	nHInForest = 0;
	nTInForest = 0;
	//Go through all haloes and count the number of trees and haloes in each forest
	for (i = 1; i < Nhalos; i++)
	{
		//Increment the number of haloes in this forest
		nHInForest++;
		//If we are at the root then increment the number of trees in this forest
		if (H[i].descid == 0 && H[i].type != 2) nTInForest++;
		//If we arrive at a new forest
		if (H[i].forestid != H[i-1].forestid)
		{
			//copy the number of trees and haloes to this forest
			NhalosInForest[iforest] = nHInForest;
			NtreesInForest[iforest] = nTInForest;
			//Set the number of trees and haloes per forest back to zero
			nHInForest = 0;
			nTInForest = 0;
			//Go to new forest
			iforest++;
		}
	}
	//For the last forest we need to copy the number of trees and haloes
	NhalosInForest[iforest] = nHInForest + 1;
	NtreesInForest[iforest] = nTInForest;
	//And increment the number of forests one last time
	iforest++;

	if (iforest != Nforests) endrun("The number of forests after sorting is not the same as the stored number");
#endif

	//Finally restore all the values that have been used for the sorting (icoprog, idescgal and possibly imass)
#ifndef COMPUTE_ICM
	//Loop through all trees
	for (itree = 0, ihalo = 0; itree < Ntrees; itree++)
	{
		//Go through all haloes in this tree
		for (i = 0; i < NhalosInTree[itree]; i++, ihalo++)
		{
#else
	//Loop through all forests
	for (iforest = 0, ihalo = 0; iforest < Nforests; iforest++)
	{
		//Go through all haloes in this tree
		for (i = 0, imass = 0; i < NhalosInForest[iforest]; i++, ihalo++)
		{
#endif
			//Reset idescgal and icoprog
			H[ihalo].icoprog  = -1;
			H[ihalo].idescgal = H[ihalo].idesc;
			//If we have also used imass restore that as well
#ifdef COMPUTE_ICM
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
#endif
		}
	}
	//For formating purposes add the two braces
#ifdef NEVERDEF
}
}
#endif

}


/*! \brief This function distributes the trees to the tasks as evenly as possible
 *
 *  The tasks are ordered by their load. Then the task containing the heaviest load sends haloes
 *  to tasks with the lightest load. If the sending task has sent enough haloes the next high-loaded
 *  task in line sends haloes. If the receiving task has received enough haloes the next low-loaded
 *  task in line receives haloes. Haloes are sent in packeges of full trees.
 *
 *  \param ntask Number of tasks that the trees will be distributed over
 *  \param load Load of each of tasks that the trees will be distributed over
 */
void distribute_trees(int ntask, double load)
{
	int i,j,maxlen;
	int nhalos_send_target,nhalos_send_actual,ntreesInTask,ntrees_send;
	int isend, irecv;
	int send_this_turn, left_to_send;
	int *nhalosInTask, *taskSorted;
	int *ip;
	unsigned long long tothalos, NtreesLocal, NhalosLocal;
	double totload, load_average;
	double *loadInTask;
	char buf[500];
	struct halo *hp;
	struct task_list *taskList;
	MPI_Status status;

#ifdef COMPUTE_ICM
	int Nforests_send, NforestsInTask;
#endif

	//Allocate the CommBuffer
	CommBuffer = emalloc_movable(&CommBuffer,"CommBuffer", All.BufferSize * 1024 * 1024);

	//Allocate arrays to hold the information on the task number and the number of haloes it contains
	nhalosInTask = emalloc_movable(&nhalosInTask,"NhalosInTask", NTask * sizeof(int));
	loadInTask   = emalloc_movable(&loadInTask,  "LoadInTask",   NTask * sizeof(double));
	taskSorted   = emalloc_movable(&taskSorted,  "TaskSorted",   NTask * sizeof(int));
	taskList     = emalloc("TaskList", NTask * sizeof(struct task_list));

	//Gather the number of haloes from each task and write to local task in nhalosInTask
	MPI_Allgather(&Nhalos, 1, MPI_INT,    nhalosInTask, 1, MPI_INT,    MPI_COMM_WORLD);
	MPI_Allgather(&load,   1, MPI_DOUBLE, loadInTask,   1, MPI_DOUBLE, MPI_COMM_WORLD);

	//Get the total number of haloes and the total load
	NhalosLocal = (unsigned long long)(Nhalos);
	MPI_Allreduce(&NhalosLocal, &tothalos, 1, MPI_UNSIGNED_LONG_LONG,    MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&load,        &totload,  1, MPI_DOUBLE,                MPI_SUM, MPI_COMM_WORLD);

	//Compute the average load each task should ideally hold
	load_average = totload / (double)(ntask);

	//Assign the tasks to the task_list
	for (i=0; i<ntask; i++)
	{
		taskList[i].NhalosInTask = nhalosInTask[i];
		taskList[i].LoadInTask   = loadInTask[i];
		taskList[i].TaskSorted   = i;
		//Check if each task can make enough space for the new trees and haloes
		if ((long long)(((int)((max(load_average-loadInTask[i],0.0)) * (double)(tothalos)/totload) + nhalosInTask[i] - All.MaxHalos) * sizeof(struct halo)) > (long long)(FreeBytes))
		{
			//If not then free all arrays and skip the re-distribution
			if (ThisTask==0) printf("%s\n%s Cannot distribute the haloes equally given this load. Skipping the re-distribution.\n",All.fullline,All.startline);
			efree(taskList);
			efree_movable(taskSorted);
			efree_movable(loadInTask);
			efree_movable(nhalosInTask);
			efree_movable(CommBuffer);
			return;
		}
	}
	//Sort the task_list according to their load
	qsort(taskList, ntask, sizeof(struct task_list), compare_load);
	//Assign the values in the task_list back to the arrays
	for (i=0; i<ntask; i++)
	{
		nhalosInTask[i] = taskList[i].NhalosInTask;
		loadInTask[i]   = taskList[i].LoadInTask;
		taskSorted[i]   = taskList[i].TaskSorted;
	}
	efree(taskList);

	//Write the total number of trees and haloes
	if (ThisTask == 0)
	{
		printf("%s\n",All.fullline);
		printf("%s Distributing all %llu trees with %llu haloes over %d tasks.\n",All.startline,All.TotNtrees,All.TotNhalos,ntask);
	}
	MPI_Barrier(MPI_COMM_WORLD);

	//Set the number of trees in this task to its initial value
	ntreesInTask = Ntrees;

#ifdef COMPUTE_ICM
	//Set the number of forests in this task to its initial value
	NforestsInTask = Nforests;
#endif

	//Start with sending from the fullest task and receiving at the emptiest task
	isend = 0;
	irecv = ntask - 1;

	//Loop while there are still haloes to send
	while (isend < irecv)
	{
		//Compute the target number of haloes to send ignoring that haloes come in trees
		//Set this number to whatever is smaller: number of haloes the receiving task wants to take
		//or the number of haloes the sending task wants to send
		//All numbers are determined from the load of each task
		nhalos_send_target = min((int)((loadInTask[isend]-load_average) * (double)(tothalos)/totload),(int)((load_average-loadInTask[irecv]) * (double)(tothalos)/totload));

		//Sending task
		//Send closest to nhalos_send_target haloes (in tree increments)
		if (ThisTask == taskSorted[isend])
		{
			//Start counting haloes from last tree
			i = ntreesInTask - 1;
#ifdef COMPUTE_ICM
			//Start counting haloes from last forest
			j = NforestsInTask - 1;
#endif
			nhalos_send_actual = 0;
			//Loop until the actual number of haloes to send exceeds the target number
			while (nhalos_send_actual < nhalos_send_target)
			{
#ifndef COMPUTE_ICM
				//Increment the number of haloes we will send by the number of haloes in this tree
				nhalos_send_actual += NhalosInTree[i];
				i--;
#else
				//Increment the number of haloes we will send by the number of haloes in this forest
				nhalos_send_actual += NhalosInForest[j];
				//Increment the number of trees by the number of trees in this forest
				i -= NtreesInForest[j];
				j--;
#endif
			}
			//Record the number of trees that will be sent
			ntrees_send = ntreesInTask - i - 1;
#ifdef COMPUTE_ICM
			//Record the number of forests that will be sent
			Nforests_send = NforestsInTask - j - 1;
#endif
		}

		//Send nhalos_send_actual and ntrees_send to all tasks
		MPI_Bcast(&nhalos_send_actual, 1, MPI_INT, taskSorted[isend], MPI_COMM_WORLD);
		MPI_Bcast(&ntrees_send, 1, MPI_INT, taskSorted[isend], MPI_COMM_WORLD);
		//Decrease number of trees for sending task and increase it for receiving task
		if (ThisTask == taskSorted[isend]) ntreesInTask -= ntrees_send;
		if (ThisTask == taskSorted[irecv]) ntreesInTask += ntrees_send;
#ifdef COMPUTE_ICM
		//Decrease number of forests for sending task and increase it for receiving task
		MPI_Bcast(&Nforests_send, 1, MPI_INT, taskSorted[isend], MPI_COMM_WORLD);
		if (ThisTask == taskSorted[isend]) NforestsInTask -= Nforests_send;
		if (ThisTask == taskSorted[irecv]) NforestsInTask += Nforests_send;
#endif

		//If zero haloes would be sent, break the loop
		if (nhalos_send_actual <= 0) break;

		//Print what will be sent
		if (ThisTask == 0)
			printf("%s Moving %8d trees with %10d haloes from task %d to task %d.\n",All.startline,ntrees_send,nhalos_send_actual,taskSorted[isend],taskSorted[irecv]);

		//Now send all tree information and then all haloes

		//Check if the CommBuffer is large enough to hold ntrees_send
		if (ThisTask == taskSorted[isend])
		{
			if(ntrees_send * sizeof(int) > All.BufferSize * 1024 * 1024)
			{
				sprintf(buf, "CommBuffer not large enough. %d %d\n", (int) (ntrees_send * sizeof(int)), (int) (All.BufferSize * 1024 * 1024));
				endrun(buf);
			}
		}

		//Check if the receiving task has enough space for the new trees
		if(ntreesInTask > All.MaxTrees)
		{
			//If we need more trees than currently allocated then reallocate
			All.MaxTrees = ntreesInTask;
			NhalosInTree = (int *) erealloc_movable(NhalosInTree, All.MaxTrees * sizeof(int));
#ifdef COMPUTE_ICM
			NtreesInForest = (int *) erealloc_movable(NtreesInForest, All.MaxTrees * sizeof(int));
			NhalosInForest = (int *) erealloc_movable(NhalosInForest, All.MaxTrees * sizeof(int));
#endif
		}

		//Check if the receiving task has enough space for the new haloes
		if (nhalos_send_actual + nhalosInTask[irecv] > All.MaxHalos)
		{
			//If we can make space for the new trees then increase MaxHaloes and reallocate
			if ((long long)((nhalos_send_actual + nhalosInTask[irecv] - All.MaxHalos) * sizeof(struct halo)) < (long long)(FreeBytes))
			{
				All.MaxHalos = nhalos_send_actual + nhalosInTask[irecv];
				H  = (struct halo *) erealloc_movable(H, All.MaxHalos * sizeof(struct halo));
				if (ThisTask==0) printf("%s Reallocating the haloes with a size of %d\n",All.startline,All.MaxHalos);
			}
			//Otherwise free all arrays and skip the re-distribution
			else
			{
				if (ThisTask==0) printf("%s\n%s Cannot distribute the haloes equally given this load. Skipping the re-distribution.\n",All.fullline,All.startline);
				efree_movable(taskSorted);
				efree_movable(loadInTask);
				efree_movable(nhalosInTask);
				efree_movable(CommBuffer);
				return;
			}
		}

		if (ThisTask == taskSorted[isend] && ntrees_send > 0)
		{
			//Load NhalosInTree for the trees that will be sent into the CommBuffer
			ip = (int *) CommBuffer;
			for (i = ntreesInTask; i < ntreesInTask+ntrees_send; i++)
			{
				//Write CommBuffer from NhalosInTree
				*ip++ = NhalosInTree[i];
				NhalosInTree[i] = 0;
			}
			//If there are trees to send then send the CommBuffer to the receiving task
			MPI_Ssend(CommBuffer, (sizeof(int)) * ntrees_send, MPI_INT, taskSorted[irecv], TAG_NHALOSINTREE, MPI_COMM_WORLD);

#ifdef COMPUTE_ICM
			//Load NtreesInForest for the trees that will be sent into the CommBuffer
			ip = (int *) CommBuffer;
			for (i = NforestsInTask; i < NforestsInTask+Nforests_send; i++)
			{
				//Write CommBuffer from NtreesInForest
				*ip++ = NtreesInForest[i];
				NtreesInForest[i] = 0;
			}
			//If there are forests to send then send the CommBuffer to the receiving task
			MPI_Ssend(CommBuffer, (sizeof(int)) * Nforests_send, MPI_INT, taskSorted[irecv], TAG_NTREESINFOREST, MPI_COMM_WORLD);

			//Load NtreesInForest for the trees that will be sent into the CommBuffer
			ip = (int *) CommBuffer;
			for (i = NforestsInTask; i < NforestsInTask+Nforests_send; i++)
			{
				//Write CommBuffer from NhalosInForest
				*ip++ = NhalosInForest[i];
				NhalosInForest[i] = 0;
			}
			//If there are forests to send then send the CommBuffer to the receiving task
			MPI_Ssend(CommBuffer, (sizeof(int)) * Nforests_send, MPI_INT, taskSorted[irecv], TAG_NHALOSINFOREST, MPI_COMM_WORLD);
#endif
		}

		//If there are trees to send then receive the CommBuffer from the sending task
		if (ThisTask == taskSorted[irecv] && ntrees_send > 0) {
			MPI_Recv(CommBuffer, (sizeof(int)) * ntrees_send, MPI_INT, taskSorted[isend], TAG_NHALOSINTREE, MPI_COMM_WORLD, &status);
			ip = (int *) CommBuffer;
			//Read CommBuffer and write to NhalosInTree
			for (j = 0, i = ntreesInTask - ntrees_send; j < ntrees_send; i++,j++) NhalosInTree[i] = ip[j];

#ifdef COMPUTE_ICM
			MPI_Recv(CommBuffer, (sizeof(int)) * Nforests_send, MPI_INT, taskSorted[isend], TAG_NTREESINFOREST, MPI_COMM_WORLD, &status);
			ip = (int *) CommBuffer;
			//Read CommBuffer and write to NtreesInForest
			for (j = 0, i = NforestsInTask - Nforests_send; j < Nforests_send; i++,j++) NtreesInForest[i] = ip[j];

			MPI_Recv(CommBuffer, (sizeof(int)) * Nforests_send, MPI_INT, taskSorted[isend], TAG_NHALOSINFOREST, MPI_COMM_WORLD, &status);
			ip = (int *) CommBuffer;
			//Read CommBuffer and write to NhalosInForest
			for (j = 0, i = NforestsInTask - Nforests_send; j < Nforests_send; i++,j++) NhalosInForest[i] = ip[j];
#endif
		}

		//Set the maximum number of haloes that can be fitted onto the CommBuffer
		maxlen       = ((int) (All.BufferSize * 1024 * 1024)) / (sizeof(struct halo));
		left_to_send = nhalos_send_actual;

		//Loop until all haloes have been sent
		do
		{
			//Send this turn the minimum of all haloes left to send or the maximum number that fits onto the Buffer
			send_this_turn = left_to_send;
			if (send_this_turn > maxlen) send_this_turn = maxlen;

			//if (ThisTask == 0) printf("%d %d %d\n",nhalos_send_actual,left_to_send,send_this_turn);

			//Sending task
			if (ThisTask == taskSorted[isend] && send_this_turn > 0)
			{
				// Write send_this_turn haloes to CommBuffer
				hp = (struct halo *) CommBuffer;
				for (i = nhalosInTask[isend]-left_to_send; i < nhalosInTask[isend]-left_to_send+send_this_turn; i++)
					*hp++ = H[i];
				memset(&H[nhalosInTask[isend]-left_to_send], 0, send_this_turn * sizeof(struct halo));
				MPI_Ssend(CommBuffer, (sizeof(struct halo)) * send_this_turn, MPI_BYTE, taskSorted[irecv], TAG_HDATA, MPI_COMM_WORLD);
			}

			//Receiving task
			if(ThisTask == taskSorted[irecv] && send_this_turn > 0)
			{
				//Receive CommBuffer from taskSorted[isend]
				MPI_Recv(CommBuffer, (sizeof(struct halo)) * send_this_turn, MPI_BYTE, taskSorted[isend], TAG_HDATA, MPI_COMM_WORLD, &status);
				//Write send_this_turn haloes in CommBuffer to end of halo array
				hp = (struct halo *) CommBuffer;
				for (j = 0, i = nhalosInTask[irecv] + nhalos_send_actual - left_to_send; j < send_this_turn; i++,j++) H[i] = hp[j];
			}

			//Haloes left to send are reduced by what has been sent this turn
			left_to_send -= send_this_turn;
		}
		//Loop until all haloes are sent
		while (left_to_send > 0);

		//Update the number of haloes and the load in the sending and receiving task locally
		nhalosInTask[isend] -= nhalos_send_actual;
		nhalosInTask[irecv] += nhalos_send_actual;
		loadInTask[isend] -= (double)(nhalos_send_actual)/(double)(tothalos) * totload;
		loadInTask[irecv] += (double)(nhalos_send_actual)/(double)(tothalos) * totload;
		//Update the number of forests/trees/haloes in the sending and receiving task
		if (ThisTask == taskSorted[isend]) Ntrees = ntreesInTask;
		if (ThisTask == taskSorted[irecv]) Ntrees = ntreesInTask;
		if (ThisTask == taskSorted[isend]) Nhalos = nhalosInTask[isend];
		if (ThisTask == taskSorted[irecv]) Nhalos = nhalosInTask[irecv];
#ifdef COMPUTE_ICM
		if (ThisTask == taskSorted[isend]) Nforests = NforestsInTask;
		if (ThisTask == taskSorted[irecv]) Nforests = NforestsInTask;
#endif
		//If sending task has a lighter load than average go to next task
		//If receiving task has a heavier load than average go to next task
		if (loadInTask[isend] <= load_average * 1.000001) isend++;
		if (loadInTask[irecv] >= load_average / 1.000001) irecv--;
	}

	efree_movable(taskSorted);
	efree_movable(loadInTask);
	efree_movable(nhalosInTask);
	efree_movable(CommBuffer);

	// Find maximum number of trees and haloes for all tasks
	MPI_Allreduce(&Ntrees, &All.MaxTrees, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
	MPI_Allreduce(&Nhalos, &All.MaxHalos, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

	// Find total number of trees and haloes for all tasks
	NtreesLocal = (unsigned long long)(Ntrees);
	NhalosLocal = (unsigned long long)(Nhalos);
	MPI_Allreduce(&NtreesLocal, &All.TotNtrees, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&NhalosLocal, &All.TotNhalos, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

	//Reallocate to decrease the memory
	NhalosInTree   = (int *) erealloc_movable(NhalosInTree, All.MaxTrees * sizeof(int));
#ifdef COMPUTE_ICM
	NtreesInForest = (int *) erealloc_movable(NtreesInForest, All.MaxTrees * sizeof(int));
	NhalosInForest = (int *) erealloc_movable(NhalosInForest, All.MaxTrees * sizeof(int));
#endif
	H  = (struct halo *) erealloc_movable(H, (int)(All.MaxHalos) * sizeof(struct halo));

	//Set the offsets of the first halo in the trees/forests
	OffsetHalos[0] = 0;
#ifndef COMPUTE_ICM
	for (i = 1; i < Ntrees; i++) OffsetHalos[i] = OffsetHalos[i-1] + NhalosInTree[i-1];
#else
	for (i = 1; i < Nforests; i++) OffsetHalos[i] = OffsetHalos[i-1] + NhalosInForest[i-1];
#endif

	//Report memory usage
	HighMark = 0;
	report_memory_usage(&HighMark, "distribute_trees");

	MPI_Barrier(MPI_COMM_WORLD);

}


/*! \brief This function calculates the number of timesteps for the trees and creates a list
 *
 *  The scale factors of all haloes in each task are sorted. Then the number of timesteps is
 *  determined by counting how often it increases. Further a list containing all unique timesteps
 *  on this task is stored. The maximum number of timesteps over all tasks is then computed and
 *  all timestep lists are gathered on each task. This combined list is sorted again and the
 *  global number of timesteps is counted as above including a global list containing all unique
 *  timesteps. Also the table containing the fraction of mass left between scale factors is computed.
 */
void get_timesteps(void)
{

	int i, j, nTimesteps, maxnTimesteps;
	float *a, *alist, *alistall;
	char buf[500];

	//Print what is done...
	if (ThisTask == 0) printf("%s\n%s Identifying the number of timesteps and the list of scale factors from the trees...\n", All.fullline,All.startline);

	MPI_Barrier(MPI_COMM_WORLD);

	//First the arrays are allocated (MAXTIMESTEPS is used for ScaleFactor initially)
	ScaleFactor = emalloc_movable(&ScaleFactor,"ScaleFactor", MAXTIMESTEPS * sizeof(float));
	a           = emalloc("ScaleFactorHalo", All.MaxHalos * sizeof(float));
	alist       = emalloc("ScaleFactorList", sizeof(float));

	//Write scale factors of all trees in array
#if defined(_OPENMP) && (OPENMPTHREADS > 1)
#pragma omp parallel for
#endif
	for (i=0; i<Nhalos; i++) a[i] = H[i].a;

	//Sort scale factors
	qsort(a, Nhalos, sizeof(float), compare_float);

	//Initialise number of timesteps and first scale factor
	nTimesteps = 1;
	alist[0]     = a[0];

	//Check for each entry if a increases
	for (i=1; i<Nhalos; i++)
	{
		//If a increases increment the number of timesteps
		if (a[i]>a[i-1]) {
			nTimesteps++;
			//Increase the size of alist
			alist = (float *) erealloc(alist, nTimesteps * sizeof(float));
			//Write the new unique scale factor to alist
			alist[nTimesteps-1] = a[i];
		}
	}

	//Get the maximum number of timesteps across all tasks
	MPI_Allreduce(&nTimesteps, &maxnTimesteps, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

	//Increase the size of alist to maxnTimesteps and fill with zeros
	alist = (float *) erealloc(alist, maxnTimesteps * sizeof(float));
	for (i=nTimesteps; i<maxnTimesteps; i++) alist[i] = 0.0;

	//Allocate the array containing all scale factors from all tasks
	alistall = emalloc("ScaleFactorListAll", maxnTimesteps * NTask * sizeof(float));

	//Gather all scale factors from the other tasks in alistall
	MPI_Allgather(alist, maxnTimesteps, MPI_FLOAT, alistall, maxnTimesteps, MPI_FLOAT, MPI_COMM_WORLD);

	//Sort this list
	qsort(alistall, maxnTimesteps * All.NTaskPerUniverse, sizeof(float), compare_float);

	//Initialise number of timesteps
	All.NTimesteps = 1;

	//Skip over all zeros at the beginning and initialise the first element of ScaleFactor
	i = 0;
	while (alistall[i]<=0.0) i++;
	ScaleFactor[0] = alistall[i];

	//Continue at this position and check for each entry if a increases
	for (i = i+1; i < maxnTimesteps * NTask; i++)
	{
		//If a increases increment the global number of timesteps
		if (alistall[i]>alistall[i-1])
		{
			All.NTimesteps++;
			//Write the new unique scale factor to the global list
			ScaleFactor[All.NTimesteps-1] = alistall[i];
		}
	}

	//Print the number of timesteps and the starting and final scale factor
	if (ThisTask == 0) printf("%s There are %d timesteps ranging from a=%f to a=%f\n",All.startline,All.NTimesteps,ScaleFactor[0],ScaleFactor[All.NTimesteps-1]);

	//Check that the number of Timesteps does not exceed the maximum
	if (All.NTimesteps > MAXTIMESTEPS)
	{
		sprintf(buf, "The number of timesteps is %d which is larger than the maximum allowed number of %d.", nTimesteps, MAXTIMESTEPS);
		endrun(buf);
	}

	//Check that the number of Timesteps is not too low (at least 4 for the splines - though in practice much higher)
	if (All.NTimesteps < 4)
	{
		sprintf(buf, "The number of timesteps is %d which is lower than the minimim allowed number of 4.", nTimesteps);
		endrun(buf);
	}

	efree(alistall);
	efree(alist);
	efree(a);

	//Reallocate with a size of the list of the actual global number of timesteps
	ScaleFactor = (float *) erealloc_movable(ScaleFactor, All.NTimesteps * sizeof(float));

	CosmicTime = emalloc_movable(&CosmicTime,"CosmicTime",    All.NTimesteps * sizeof(float));
	Timestep   = emalloc_movable(&Timestep,  "Timestep",      All.NTimesteps * sizeof(float));
	DynTime    = emalloc_movable(&DynTime,   "DynamicalTime", All.NTimesteps * sizeof(float));

	for (i = 0; i < All.NTimesteps; i++)
	{
		CosmicTime[i] = cosmictime(ScaleFactor[i]);
		if (i == 0) Timestep[i] = CosmicTime[i];
		else Timestep[i] = CosmicTime[i] - CosmicTime[i-1];
		DynTime[i] = tdyn(ScaleFactor[i]);
	}

	//Loop through all haloes and set its iscale to the corresponding index in ScaleFactor
#if defined(_OPENMP) && (OPENMPTHREADS > 1)
#pragma omp parallel for
#endif
	for (i=0; i<Nhalos; i++)
	{
		//Search through the list and return the index
		H[i].iscale = binary_search_float(ScaleFactor, All.NTimesteps, H[i].a);
		//If the index cannot be found or is larger than allowed end programme
		if (H[i].iscale < 0 || H[i].iscale >= All.NTimesteps)
			endrun("Halo scale factor cannot be found in the list.");
	}

	//Go through all output redshifts and find the corresponding snapshot number
	for (i = 0; i < All.Noutputredshifts; i++)
	{ //Find the index for this output scale factor
		for (j = 0; j < All.NTimesteps; j++)
		{
			if (ScaleFactor[j] >= 1./(OutputRedshifts[i]+1.0)) break;
		}
		if (j == 0) Output_iscale[i] = 0;
		else
		{
			if (ScaleFactor[j]-1./(OutputRedshifts[i]+1.0) >= 1./(OutputRedshifts[i]+1.0)-ScaleFactor[j-1])
				Output_iscale[i] = j-1;
			else Output_iscale[i] = j;
		}
	}

#ifdef WRITE_MAINBRANCH
	for (j = 0; j < All.NTimesteps; j++)
	{
		if (ScaleFactor[j] >= 1./(All.mainBranchRedshift+1.0)) break;
	}
	if (j == 0) All.MainBranch_iscale = 0;
	else
	{
		if (ScaleFactor[j]-1./(All.mainBranchRedshift+1.0) >= 1./(All.mainBranchRedshift+1.0)-ScaleFactor[j-1])
			All.MainBranch_iscale = j-1;
		else All.MainBranch_iscale = j;
	}
#endif

#ifdef READ_WP
	//Find the index for the correlation function scale factor
	for (i = 0; i < All.NTimesteps; i++)
	{
		if (ScaleFactor[i] >= 1./(All.wpredshift+1.0)) break;
	}
	if (i == 0) All.Wpscale = 0;
	else
	{
		if (ScaleFactor[i]-1./(All.wpredshift+1.0) >= 1./(All.wpredshift+1.0)-ScaleFactor[i-1]) All.Wpscale = i-1;
		else All.Wpscale = i;
	}
#endif

	//Compute the MassLeft table which stores the fraction of mass left after any time interval
	MassLeft = emalloc_movable(&MassLeft,"MassLeft", All.NTimesteps * All.NTimesteps * sizeof(float));
	for (i = 0; i < All.NTimesteps; i++)
	{
		for (j = 0; j < All.NTimesteps; j++)
		{
			if (j < i) MassLeft[i*All.NTimesteps+j] = 1.0 - 0.05 * log( ((CosmicTime[i]-CosmicTime[j])*All.t_unit/1.4e6) + 1.0);
			else MassLeft[i*All.NTimesteps+j] = 1.0;
		}
	}

	//Broadcast so that every universe has the joy
	MPI_Bcast(ScaleFactor,     All.NTimesteps, MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Bcast(CosmicTime,      All.NTimesteps, MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Bcast(Timestep,        All.NTimesteps, MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Bcast(DynTime,         All.NTimesteps, MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Bcast(MassLeft,        All.NTimesteps * All.NTimesteps, MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Bcast(Output_iscale,   All.Noutputredshifts, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(OutputRedshifts, All.Noutputredshifts, MPI_FLOAT, 0, MPI_COMM_WORLD);

}


/*! \brief This function adds orphans to each tree
 *
 *  First the total number of orphans on each task is computed by summing the number of
 *  timesteps left for each leaf. Then the halo arrays is reallocated to make space for
 *  the orphans. The maximum ID is calculated throughout the simulation and is used as
 *  the first orphan ID. After the number of orphans in each tree is computed we shift
 *  the haloes in the memory for each tree to make space for the orphans. We then go
 *  through all leaves and detach all merging haloes from their branch and add orphans
 *  as their descendants. Finally the orphans are sorted into each tree by scale factor.
 */
void add_orphans(void)
{
	int i,j,ihalo,iorphan,itree,idesc,ithis;
	int nhalos_cumulative, norphans_cumulative;
	int *nOrphansInTree, *nOrphansInTask;
	IDType maxid, maxhaloIDalltasks;

#ifdef COMPUTE_ICM
	int iforest = 0;
#endif

	//Calculate the number of orphans on this task and the maximum over all tasks
	Norphans = 0;
	for (i = 0; i < Nhalos; i++) if (H[i].mmp==0) Norphans += All.NTimesteps - H[i].iscale - 1;
	MPI_Allreduce(&Norphans, &All.MaxOrphans, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

	//Print what is done...
	if (ThisTask == 0) printf("%s\n%s Inserting %12d possible orphans in each tree...\n", All.fullline,All.startline,All.MaxOrphans);

	//Calculate the maximum number of haloes including the orphans and reallocate the halo array
	All.MaxHalos = All.MaxHalos + All.MaxOrphans;
	H  = (struct halo *) erealloc_movable(H, All.MaxHalos * sizeof(struct halo));
	memset(&H[Nhalos], 0, (All.MaxHalos-Nhalos) * sizeof(struct halo));

	//Check maximum halo ID in tree
	maxid = H[0].haloid;
	for (i = 0; i < Nhalos; i++) if (H[i].haloid > maxid) maxid = H[i].haloid;

#ifndef LONGIDS
	//If the new IDs cannot be represented by IDType end the programme
	if ((unsigned long long)(maxid) + (unsigned long long)(All.MaxOrphans) > (1ULL << (sizeof(IDType)*8)) - 1ULL && ThisTask==0) endrun("The bit size of the IDType is not large enough to hold all IDs. Use LONGIDS!");
#endif

	// Find maximum halo ID across all tasks
	maxhaloIDalltasks = 0; //probably not needed
#ifndef LONGIDS
	MPI_Allreduce(&maxid, &maxhaloIDalltasks, 1, MPI_UNSIGNED, MPI_MAX, MPI_COMM_WORLD);
#else
	MPI_Allreduce(&maxid, &maxhaloIDalltasks, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
#endif
	maxid = maxhaloIDalltasks;

	//Allocate and compute array to store the number of orphans in each task
	nOrphansInTask = emalloc("NOrphansInTask", NTask * sizeof(int));
	MPI_Allgather(&Norphans, 1, MPI_INT, nOrphansInTask, 1, MPI_INT, MPI_COMM_WORLD);

	//For each task determine the starting orphan id
	for (i=0; i<All.NTaskPerUniverse; i++)
	{
		if (ThisTask>i) maxid += nOrphansInTask[i];
	}

	//Allocate and compute array to store the number of orphans in each tree
	nOrphansInTree = emalloc("NOrphansInTree", Ntrees * sizeof(int));
	for (itree = 0, ihalo = 0; itree < Ntrees; itree++)
	{
		nOrphansInTree[itree] = 0;
		//Now all the haloes in the tree
		for (i = 0; i < NhalosInTree[itree]; i++, ihalo++)
		{
			if (H[ihalo].mmp==0) nOrphansInTree[itree] += All.NTimesteps - H[ihalo].iscale - 1;
		}
	}

	//Initialise the cumulative numbers of haloes and orphans
	nhalos_cumulative = norphans_cumulative = 0;

	//Go through the trees backwards
	for (itree = Ntrees-1; itree >= 0; itree--)
	{
		//Get the cumulative numbers of haloes and orphans up to this tree
		nhalos_cumulative += NhalosInTree[itree];
		norphans_cumulative += nOrphansInTree[itree];
		//Move the haloes in this tree towards the end of the array and leave space for the orphans
		memmove(&H[Nhalos + Norphans - nhalos_cumulative - norphans_cumulative], &H[Nhalos-nhalos_cumulative], NhalosInTree[itree] * sizeof(struct halo));
		//If there are orphans
		if (nOrphansInTree[itree]>0)
		{
			//Set all orphans to zero
			memset(&H[Nhalos + Norphans - nhalos_cumulative - norphans_cumulative + NhalosInTree[itree]], 0, nOrphansInTree[itree] * sizeof(struct halo));
		}
	}

	//Set the orphan index to zero
	iorphan = 0;

	//Now go through all trees separately
	for (itree = 0, ihalo = 0; itree < Ntrees; itree++)
	{
		//increment the orphan index to be the first after all regular haloes
		iorphan += NhalosInTree[itree];

		//Go through all haloes in this tree
		for (i = 0; i < NhalosInTree[itree]; i++, ihalo++)
		{
			//Set the halo type to 0 for main haloes (upid==0) and to 1 for subhaloes (otherwise)
			if (H[ihalo].upid == 0) H[ihalo].type = 0;
			else H[ihalo].type = 1;

			//Create a new orphan if the halo is not the most massive progenitor i.e. merging
			if (H[ihalo].mmp==0)
			{
				//Go from this halo down to its descendant for which the halo ID is this halo's descendant ID
				idesc = ihalo;
				while (H[idesc].haloid != H[ihalo].descid) idesc++;

				//Go from the descendant up to its progenitor for which the descendant ID is this halo's ID
				ithis = idesc;
				while (H[ithis].descid != H[idesc].haloid) ithis--;

				//Detach the halos from the larger branch
				//Descendant ID is now the orphan which is set t maxid+1
				H[ihalo].descid     = maxid+1;
				//The halo is now the MMP (of the orphan)
				H[ihalo].mmp        = 1;
				//The previous descendant has now only one progenitor
				H[idesc].np         = 1;

				//Now at the next timestep we create the orphan
				H[iorphan].haloid   = maxid+1;
				//Unless its true descendant is 0 create the next orphan descendant
				if (H[idesc].descid==0) H[iorphan].descid = 0;
				else H[iorphan].descid = maxid+2;
				//We now use the upid of the orphan to record the true descendants halo ID
				H[iorphan].upid     = H[idesc].haloid;

				//Copy all relevant information from the true descendant
				H[iorphan].iscale   = H[idesc].iscale;
				H[iorphan].type     = 2; //All orphans are type 2
				H[iorphan].np       = 1; //The orphans can never be leaves and nothing merges into them
				H[iorphan].mmp      = 1; //Always the main progenitor
				H[iorphan].a        = H[idesc].a;
#ifdef ORPHAN_MASSLOSS
				H[iorphan].mvir     = 0;
#else
				H[iorphan].mvir     = H[ihalo].mvir;
#endif
				H[iorphan].rvir     = H[ihalo].rvir;
				H[iorphan].c        = H[ihalo].c;
				H[iorphan].lambda   = H[ihalo].lambda;
				H[iorphan].vel[0]   = H[idesc].vel[0];
				H[iorphan].vel[1]   = H[idesc].vel[1];
				H[iorphan].vel[2]   = H[idesc].vel[2];
#ifdef COMPUTE_ICM
				H[iorphan].forestid = H[idesc].forestid;
#endif

				//The first orphan is set - increment the iorphan index and ID
				iorphan++;
				maxid++;

				//Now do create all remaining orphans until the final timestep
				for (j=H[ihalo].iscale+2; j<All.NTimesteps; j++)
				{
					//Record the current descendant ID
					ithis = idesc;
					//Go to the next descendant which has an ID equal to the this one's descendant ID
					while (H[idesc].haloid != H[ithis].descid) idesc++;
					//Create next orphan descendant until the descid is 0
					H[iorphan].haloid   = maxid+1;
					if (H[idesc].descid==0) H[iorphan].descid = 0;
					else H[iorphan].descid = maxid+2;
					//We still use the upid of the orphan to record the true descendants halo ID
					H[iorphan].upid     = H[idesc].haloid;
					//Copy all relevant information from the true descendant
					H[iorphan].iscale   = H[idesc].iscale;
					H[iorphan].type     = 2;
					H[iorphan].np       = 1;
					H[iorphan].mmp      = 1;
					H[iorphan].a        = H[idesc].a;
#ifdef ORPHAN_MASSLOSS
					H[iorphan].mvir     = 0;
#else
					H[iorphan].mvir     = H[ihalo].mvir;
#endif
					H[iorphan].rvir     = H[ihalo].rvir;
					H[iorphan].c        = H[ihalo].c;
					H[iorphan].lambda   = H[ihalo].lambda;
					H[iorphan].vel[0]   = H[idesc].vel[0];
					H[iorphan].vel[1]   = H[idesc].vel[1];
					H[iorphan].vel[2]   = H[idesc].vel[2];
#ifdef COMPUTE_ICM
					H[iorphan].forestid = H[idesc].forestid;
#endif

					//This orphan is set - increment the iorphan index and ID
					iorphan++;
					maxid++;
				}//Done looping through this orphan's descendant

			}//Done with this leaf
		}//Done with all haloes in this tree

		//increment the halo index by the number of orphans we just created
		ihalo += nOrphansInTree[itree];
		//Add the number of all orphans in this tree to the haloes in this tree
		NhalosInTree[itree] += nOrphansInTree[itree];

#ifdef COMPUTE_ICM
		//Add the number of all orphans in this tree to the haloes in this forest
		NhalosInForest[iforest] += nOrphansInTree[itree];
		if (H[ihalo].forestid != H[ihalo-1].forestid) iforest++;
#endif

		//Sort the halos in this tree by the scale factor and Halo ID
		//Sorts in the orphans at a given scale factor after all real haloes
		//Might be able to improve the sorting here...
		qsort(&H[ihalo - NhalosInTree[itree]], NhalosInTree[itree], sizeof(struct halo), compare_scale);

	}//Done with this tree

	//Print the number of orphans created on each task
	if (ThisTask < All.NTaskPerUniverse)
		printf("%s Inserted %7d orphans on task %6d.\n",All.startline,Norphans,ThisTask);

	//The new number of haloes is now the old number plus the number of orphans
	Nhalos = Nhalos + Norphans;

	//Free the arrays
	efree(nOrphansInTree);
	efree(nOrphansInTask);

	//And print the memory usage
	report_memory_usage(&HighMark, "Add_Orphans");

}


/*! \brief This function links progenitors and descendants tree-by-tree
 *
 *  All merger trees are processed one at a time. All IDs of haloes in a given merger tree are
 *  sorted in a list. If a halo has no progenitors its iprog is set to -1. If it has no descendant
 *  its idesc is set to -1. Otherwise the index of its descendant is determined by searching the
 *  list. If the halo is the most massive progenitor of its descendant the descendant's iprog is
 *  set to the halo's index.
 */
void setup_haloes_by_tree(void)
{
	int i,ihalo,itree;
	int nHalosInTreeMax;
	struct search_list_IDs *id_list;

	//Print what is done...
	if (ThisTask == 0) printf("%s\n%s Connecting progenitors to descendants...\n", All.fullline,All.startline);
	MPI_Barrier(MPI_COMM_WORLD);

	//Initialise and find the maximum number of haloes per tree on this task
	nHalosInTreeMax = 0;
	for (itree = 0; itree < Ntrees; itree++) nHalosInTreeMax = max(nHalosInTreeMax,NhalosInTree[itree]);

	//Allocate memory for the id list for this tree which can be searched for progenitors and descendants
	id_list = emalloc("ID_LIST", nHalosInTreeMax * sizeof(struct search_list_IDs));

	//Loop through all trees
	for (itree = 0, ihalo = 0; itree < Ntrees; itree++)
	{
		//Write all IDs and indices of the haloes in this tree to a list
		for (i = 0; i < NhalosInTree[itree]; i++)
		{
			id_list[i].position = ihalo + i;
			id_list[i].ID       = H[ihalo+i].haloid;
		}

		//Sort the ID list using the function compare_id
		qsort(id_list, NhalosInTree[itree], sizeof(struct search_list_IDs), compare_id);

		//Go through all haloes
		for (i = 0; i < NhalosInTree[itree]; i++, ihalo++)
		{
			//If this halo does not have any progenitors set its iprog to -1 and its impeak to this halo's index
			if (H[ihalo].np < 1) H[ihalo].iprog  = -1;

			//If this halo does not have a descendant set its idesc and icoprog to -1
			if (H[ihalo].descid == 0)
			{
				H[ihalo].idesc   = -1;
			}
			//If this halo has a descendant
			else
			{
				//Search the list for the descid and write the index found to idesc
				H[ihalo].idesc = binary_search_id(id_list, NhalosInTree[itree], H[ihalo].descid);
				if (H[ihalo].idesc < 0) endrun("Cannot find descendant in this tree.\n");
				//If this is the most massive progenitor then set its descendants's iprog to the halo index
				H[H[ihalo].idesc].iprog = ihalo;
			}//End all haloes with descendants

			//As we used the upid for orphans to point to the halo its progenitor would have merged with
			//We use this ID to identify imain as the index of this halo
			if (H[ihalo].type == 2) H[ihalo].imain = binary_search_id(id_list, NhalosInTree[itree], H[ihalo].upid);
			//For non-orphnas we set imain to its own index
			else H[ihalo].imain = ihalo;

			//Now for orphans set its upid to the id of its main halo or to this halo's upid
			if (H[ihalo].type == 2)
			{
				if (H[H[ihalo].imain].upid == 0) H[ihalo].upid = H[H[ihalo].imain].haloid;
				else H[ihalo].upid = H[H[ihalo].imain].upid;
			}

		}//End loop through all haloes - increment ihalo

	}//End loop through all trees

	//Set the offsets of the first halo in the trees
	OffsetHalos[0] = 0;
	for (itree = 1; itree < Ntrees; itree++) OffsetHalos[itree] = OffsetHalos[itree-1] + NhalosInTree[itree-1];

	//Free the ID list
	efree(id_list);

	//And print the memory usage
	report_memory_usage(&HighMark, "Setup_Halos_by_Tree");

}


#ifdef COMPUTE_ICM
/*! \brief This function links progenitors and descendants forest-by-forest
 *
 *  All merger trees are processed one at a time. All IDs of haloes in a given merger tree are
 *  sorted in a list. If a halo has no progenitors its iprog is set to -1. If it has no descendant
 *  its idesc is set to -1. Otherwise the index of its descendant is determined by searching the
 *  list. If the halo is the most massive progenitor of its descendant the descendant's iprog is
 *  set to the halo's index.
 */
void setup_haloes_by_forest(void)
{
	int i,ihalo,iforest;
	int nHalosInForestMax;
	struct search_list_IDs *id_list;

	//Print what is done...
	if (ThisTask == 0) printf("%s\n%s Connecting progenitors to descendants...\n", All.fullline,All.startline);
	MPI_Barrier(MPI_COMM_WORLD);

	//Initialise and find the maximum number of haloes per forest on this task
	nHalosInForestMax = 0;
	for (iforest = 0; iforest < Nforests; iforest++) nHalosInForestMax = max(nHalosInForestMax,NhalosInForest[iforest]);

	//Allocate memory for the id list for this forest which can be searched for progenitors and descendants
	id_list = emalloc("ID_LIST", nHalosInForestMax * sizeof(struct search_list_IDs));

	//Loop through all forests
	for (iforest = 0, ihalo = 0; iforest < Nforests; iforest++)
	{
		//Sort the halos in this forest by the scale factor and Halo ID
		qsort(&H[ihalo], NhalosInForest[iforest], sizeof(struct halo), compare_scale);

		//Write all IDs and indices of the haloes in this forest to a list
		for (i = 0; i < NhalosInForest[iforest]; i++)
		{
			id_list[i].position = ihalo + i;
			id_list[i].ID       = H[ihalo+i].haloid;
		}

		//Sort the ID list using the function compare_id
		qsort(id_list, NhalosInForest[iforest], sizeof(struct search_list_IDs), compare_id);

		//Go through all haloes
		for (i = 0; i < NhalosInForest[iforest]; i++, ihalo++)
		{
			//If this halo does not have any progenitors set its iprog to -1 and its impeak to this halo's index
			if (H[ihalo].np < 1) H[ihalo].iprog  = -1;

			//If this halo does not have a descendant set its idesc and icoprog to -1
			if (H[ihalo].descid == 0)
			{
				H[ihalo].idesc   = -1;
			}
			//If this halo has a descendant
			else
			{
				//Search the list for the descid and write the index found to idesc
				H[ihalo].idesc = binary_search_id(id_list, NhalosInForest[iforest], H[ihalo].descid);
				if (H[ihalo].idesc < 0) endrun("Cannot find descendant in this forest.\n");
				//If this is the most massive progenitor then set its descendants's iprog to the halo index
				H[H[ihalo].idesc].iprog = ihalo;
			}//End all haloes with descendants

			//As we used the upid for orphans to point to the halo its progenitor would have merged with
			//We use this ID to identify imain as the index of this halo
			if (H[ihalo].type == 2) H[ihalo].imain = binary_search_id(id_list, NhalosInForest[iforest], H[ihalo].upid);
			//For non-orphnas we set imain to its own index
			else H[ihalo].imain = ihalo;

			//Now for orphans set its upid to the id of its main halo or to this halo's upid
			if (H[ihalo].type == 2)
			{
				if (H[H[ihalo].imain].upid == 0) H[ihalo].upid = H[H[ihalo].imain].haloid;
				else H[ihalo].upid = H[H[ihalo].imain].upid;
			}

		}//End loop through all haloes - increment ihalo

	}//End loop through all forests

	//Set the offsets of the first halo in the forests
	OffsetHalos[0] = 0;
	for (iforest = 1; iforest < Nforests; iforest++) OffsetHalos[iforest] = OffsetHalos[iforest-1] + NhalosInForest[iforest-1];

	//Free the ID list
	efree(id_list);

	//And print the memory usage
	report_memory_usage(&HighMark, "Setup_Halos_by_Forest");

}
#endif


/*! \brief This function computes the peak mass and redshift and the baryonic accretion rate
 *
 *  The index where the halo reaches peak mass is stored throughout the merger tree.
 *  The baryonic accretion rate averaged over one dynamical time (forwards) is computed as well.
 */
void compute_halo_history(void)
{
	int ihalo,iprog,ithis;

	//Print what is done...
	if (ThisTask == 0) printf("%s\n%s Computing maximum mass and accretion rate...\n", All.fullline,All.startline);
	MPI_Barrier(MPI_COMM_WORLD);

	//Go through all haloes
	for (ihalo = 0; ihalo < Nhalos; ihalo++)
	{
		//If this halo does not have any progenitors set its iprog to -1 and its impeak to this halo's index
		if (H[ihalo].np < 1) H[ihalo].impeak = ihalo;
		//Copy the index of the peak halo mass from the progenitor
		else H[ihalo].impeak = H[H[ihalo].iprog].impeak;
		//If the current halo mass is larger than the peak mass set the peak index to this halo's index
		if (H[ihalo].mvir > H[H[ihalo].impeak].mvir) H[ihalo].impeak = ihalo;

		//Determine the baryonic accretion rate
		if (H[ihalo].iprog == -1)
		{
			//If the halo does not have a progenitor the rate is mvir/t_cosmic*fb
			H[ihalo].mdotbary = H[ihalo].mvir / CosmicTime[H[ihalo].iscale] * All.f_baryon;
		}
		else
		//If there is a progenitor go back to the progenitor that has a cosmic time + dynamical time equal to this halos cosmic time
		{
			//Identify the progenitor (cannot be -1)
			iprog = H[ihalo].iprog;
			//Check if we are at least one dynamical time back
			while (CosmicTime[H[iprog].iscale] + DynTime[H[iprog].iscale] >= CosmicTime[H[ihalo].iscale])
			{
				//Save current progenitor index and compute next progenitor index
				ithis = iprog;
				iprog = H[iprog].iprog;
				//If we are at a leave set the progenitor index to the leaf index and stop
				if (iprog < 0)
				{
					iprog = ithis;
					break;
				}
			}
			//Set the baryonic accretion rate to (mvir_j-mvir_i)/(t_j-t_i)*fbary where i_j>t_i+tdyn_i
			H[ihalo].mdotbary = (H[ihalo].mvir - H[iprog].mvir) / (CosmicTime[H[ihalo].iscale] - CosmicTime[H[iprog].iscale]) * All.f_baryon;
			//If the accretion rate is negative set it to 0
			H[ihalo].mdotbary = max(H[ihalo].mdotbary,0.0);
		}//End computation of accretion rate

	}//End loop through all haloes

	//And print the memory usage
	report_memory_usage(&HighMark, "Compute_Halo_History");

}


#ifdef COMPUTE_ICM
/*! \brief This function links all haloes to their parents
 *
 *  All halo IDs and indices are stored in an id_list and sorted by their ID. For each halo it
 *  is checked if it is a main halo (then the host index is its own index) or if it is a subhalo.
 *  In the latter case the id_list on this task is searched for the UPID. If it is found then the
 *  host id is set equal to the index of the parent.
 */
void find_parents(void)
{
	int i,ihalo,iforest,ihost,nlist;
	struct search_list_IDs *id_list;

	//Print what is done...
	if (ThisTask == 0) printf("%s\n%s Linking the haloes to their parents...\n", All.fullline,All.startline);
	MPI_Barrier(MPI_COMM_WORLD);

	//Allocate the memory for some arrays. The list id_on_other_task is allocated in increments of allocsize
	id_list = emalloc("ID_LIST", All.MaxHalos * sizeof(struct search_list_IDs));

	//Go through each forest
	for (iforest = 0, ihalo = 0; iforest < Nforests; iforest++)
	{

		//Write all IDs and indices of the haloes in this forest to a list
		for (i = 0, nlist = 0; i < NhalosInForest[iforest]; i++)
		{
			if (H[ihalo+i].type==0 || H[ihalo+i].type==1)
			{
				id_list[nlist].position = ihalo + i;
				id_list[nlist].ID       = H[ihalo+i].haloid;
				nlist++;
			}
		}

		//Sorting the ID list on each task using the function compare_id
		qsort(id_list, nlist, sizeof(struct search_list_IDs), compare_id);

		//Go through all haloes in this forest
		for (i = 0; i < NhalosInForest[iforest]; i++, ihalo++)
		{
			//If the upid is 0 it is a host halo - set the host index equal to the current index
			if (H[ihalo].upid == 0) H[ihalo].ihost = ihalo;
			//This is a subhalo so search for its parent
			else
			{
				//Search the local ID list for the upid and return the index of the parent
				ihost = binary_search_id(id_list, nlist, H[ihalo].upid);
				//If we find a parent set the host index equal to the parent's index
				if (ihost >= 0) H[ihalo].ihost    = ihost;
				//If we cannot find the parent on this task crash
				else endrun("Cannot find parent on this task.");
			}

		}//End haloes
	}//End forests

	efree(id_list);

	report_memory_usage(&HighMark, "Find_Parents");

}
#endif


/*! \brief This function copies the merger trees to the tasks for all universes
 *
 *  As the merger trees are only read by and distributed over the tasks for the first universe,
 *  we need to copy the trees to all other tasks. This needs to be done such that the haloes are
 *  on the same tasks modulo the universe. It is first checked that the number of tasks * universes
 *  is equal to or less than the total number of tasks. All tasks on the first universe then copy
 *  their data (Ntrees, Nhalos, NhalosInTree, H) to the tasks that hold universe nuniverses/2.
 *  The maximum universe that needs to be addressed in following calls is copied as well. Each
 *  universe that received data then copies this data to other universes as well, until all universes
 *  have received all data. The Haloes are copied into the CommBuffer which is then sent to the
 *  other task.
 *
 *  \param ntask Number of tasks needed for one universe
 *  \param nuniverses Number of universes that shall be processed
 */
void copy_trees_to_other_universes(int ntask, int nuniverses)
{
	int i, j, maxlen, send_this_turn, left_to_send, has_data, to_universe, to_task, from_task, max_universe;
	struct halo *hp;
	MPI_Status status;

	//Check if there are enough tasks to hold all universes
	if (ntask * nuniverses > NTask) endrun("Not enough tasks to hold all universes.");

	//Print what is done...
	if (ThisTask == 0) printf("%s\n%s Copying the merger trees to other universes...\n", All.fullline,All.startline);
	MPI_Barrier(MPI_COMM_WORLD);

	//Communicate all global parameters to the other processes
	MPI_Bcast(&All, sizeof(struct global_data), MPI_BYTE, 0, MPI_COMM_WORLD);

	//Allocate the CommBuffer
	CommBuffer = emalloc("CommBuffer", All.BufferSize * 1024 * 1024);

	//Set the maximum number of haloes that can be fitted onto the CommBuffer
	maxlen = ((int) (All.BufferSize * 1024 * 1024)) / (sizeof(struct halo));

	//Initialise values for all universe
	has_data = max_universe = to_universe = 0;

	//Set values for the first universe that holds all data initially
	if (MasterTask == 0)
	{
		has_data     = 1;
		max_universe = nuniverses;
		to_universe  = max_universe/2;
	}

	//Loop until the tasks would send their data to themselves
	while (to_universe * ntask != MasterTask)
	{ //Sending universes (tasks)
		if (has_data)
		{ //Set index of task that this task is sending the data to
			to_task = to_universe * ntask + ThisTask - MasterTask;
			//Print which universe is dealt with now
			if (ThisTask == MasterTask) printf("%s Copying trees from universe %4d to universe %4d on tasks %5d to %5d.\n", All.startline,MasterTask/ntask,to_universe,to_universe*ntask,(to_universe+1)*ntask-1);
			MPI_Ssend(&max_universe, 1, MPI_INT, to_task, TAG_REQUEST, MPI_COMM_WORLD);
			//Send Ntrees, Nhalos, NhalosInTree, Nforests, NtreesInForest, and NhalosInForest
			MPI_Ssend(&Ntrees, 1, MPI_INT, to_task, TAG_NTREES, MPI_COMM_WORLD);
			MPI_Ssend(&Nhalos, 1, MPI_INT, to_task, TAG_NHALOS, MPI_COMM_WORLD);
			MPI_Ssend(NhalosInTree, Ntrees, MPI_INT, to_task, TAG_NHALOSINTREE, MPI_COMM_WORLD);
			printf("%s Task %d to   task %d: %d trees, %d haloes\n",All.startline,ThisTask,to_task,Ntrees,Nhalos);
#ifdef COMPUTE_ICM
			MPI_Ssend(&Nforests, 1, MPI_INT, to_task, TAG_NFORESTS, MPI_COMM_WORLD);
			MPI_Ssend(NtreesInForest, Nforests, MPI_INT, to_task, TAG_NTREESINFOREST, MPI_COMM_WORLD);
			MPI_Ssend(NhalosInForest, Nforests, MPI_INT, to_task, TAG_NHALOSINFOREST, MPI_COMM_WORLD);
#endif
			MPI_Ssend(&Norphans, 1, MPI_INT, to_task, TAG_NORPHANS, MPI_COMM_WORLD);
#ifdef COMPUTE_ICM
			MPI_Ssend(OffsetHalos, Nforests, MPI_INT, to_task, TAG_OFFSETHALOS, MPI_COMM_WORLD);
#else
			MPI_Ssend(OffsetHalos, Ntrees, MPI_INT, to_task, TAG_OFFSETHALOS, MPI_COMM_WORLD);
#endif
			//Number of haloes that we want to send (all of them)
			left_to_send = Nhalos;
			//Loop until all haloes have been sent
			do
			{ //Send this turn the minimum of all haloes left to send or the maximum number that fits onto the Buffer
				send_this_turn = left_to_send;
				if (send_this_turn > maxlen) send_this_turn = maxlen;
				//Sending task
				if (send_this_turn > 0)
				{ // Write send_this_turn haloes to CommBuffer and send to universe * ntask + task
					hp = (struct halo *) CommBuffer;
					for (i = Nhalos - left_to_send; i < Nhalos - left_to_send + send_this_turn; i++)
						*hp++ = H[i];
					MPI_Send(CommBuffer, (sizeof(struct halo)) * send_this_turn, MPI_BYTE, to_universe*ntask+ThisTask-MasterTask, TAG_HDATA, MPI_COMM_WORLD);
				}
				//Haloes left to send are reduced by what has been sent this turn
				left_to_send -= send_this_turn;
			}
			//Loop until all haloes are sent
			while (left_to_send > 0);
			//Set new maximum universe and compute the new universe the data is sent to
			max_universe = to_universe;
			to_universe  = (MasterTask/ntask+max_universe)/2;
		}
		//This universe has no data yet and therefore receives data
		else
		{ //Wait until we get a message from the sending universe / task
			MPI_Probe(MPI_ANY_SOURCE, TAG_REQUEST, MPI_COMM_WORLD, &status);
			//Determine task that is sending the data
			from_task = status.MPI_SOURCE;
			//Receive the maximum universe this universe / task has to deal with
			MPI_Recv(&max_universe, 1, MPI_INT, from_task, TAG_REQUEST, MPI_COMM_WORLD, &status);
			//Receive Ntrees, Nhalos, NhalosInTree, Nforests, NtreesInForest, and NhalosInForest
			MPI_Recv(&Ntrees, 1, MPI_INT, from_task, TAG_NTREES, MPI_COMM_WORLD, &status);
			MPI_Recv(&Nhalos, 1, MPI_INT, from_task, TAG_NHALOS, MPI_COMM_WORLD, &status);
			MPI_Recv(NhalosInTree, Ntrees, MPI_INT, from_task, TAG_NHALOSINTREE, MPI_COMM_WORLD, &status);
#ifdef COMPUTE_ICM
			MPI_Recv(&Nforests, 1, MPI_INT, from_task, TAG_NFORESTS, MPI_COMM_WORLD, &status);
			MPI_Recv(NtreesInForest, Nforests, MPI_INT, from_task, TAG_NTREESINFOREST, MPI_COMM_WORLD, &status);
			MPI_Recv(NhalosInForest, Nforests, MPI_INT, from_task, TAG_NHALOSINFOREST, MPI_COMM_WORLD, &status);
#endif
			MPI_Recv(&Norphans, 1, MPI_INT, from_task, TAG_NORPHANS, MPI_COMM_WORLD, &status);
#ifdef COMPUTE_ICM
			MPI_Recv(OffsetHalos, Nforests, MPI_INT, from_task, TAG_OFFSETHALOS, MPI_COMM_WORLD, &status);
#else
			MPI_Recv(OffsetHalos, Ntrees, MPI_INT, from_task, TAG_OFFSETHALOS, MPI_COMM_WORLD, &status);
#endif
			//Number of haloes that we want to send (all of them)
			left_to_send = Nhalos;
			//Loop until all haloes have been sent
			do
			{ //Send this turn the minimum of all haloes left to send or the maximum number that fits onto the Buffer
				send_this_turn = left_to_send;
				if (send_this_turn > maxlen) send_this_turn = maxlen;
				//Receiving task
				if(send_this_turn > 0)
				{ //Receive CommBuffer from task
					MPI_Recv(CommBuffer, (sizeof(struct halo)) * send_this_turn, MPI_BYTE, from_task, TAG_HDATA, MPI_COMM_WORLD, &status);
					//Write send_this_turn haloes in CommBuffer to end of halo array
					hp = (struct halo *) CommBuffer;
					for (i = Nhalos - left_to_send, j = 0; i < Nhalos - left_to_send + send_this_turn; i++, j++)
						H[i] = hp[j];
				}
				//Haloes left to send are reduced by what has been sent this turn
				left_to_send -= send_this_turn;
			}
			//Loop until all haloes are sent
			while (left_to_send > 0);
			//Compute the new universe the data is sent to and remember that this universe / task now has data
			to_universe  = (MasterTask/ntask+max_universe)/2;
			has_data = 1;
		} //End sending and receiving options
	} //End loop through universes

	//Free the CommBuffer
	efree(CommBuffer);

}
