///////////////////////////////////////////////////////////////////////////////////////////////////
// Emerge code - File output.c                                                                 //
///////////////////////////////////////////////////////////////////////////////////////////////////
///
/// \file galaxies.c
/// \brief Contains functions that store computed data in files
///
/// This file contains all functions that print computed galaxy and halo data to output files,
/// except for that function that prints the statistics (see statistics.c).
/// A function that prints all galaxies at a specific scale factor to output files is also included.
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


#ifdef WRITE_GALAXY_CATALOG
/*! \brief This function collects all galaxy data and prints it
 *
 *  This function collects the galaxy data from all tasks and writes it to All.NumOutputFiles files. There
 *  are All.NumFilesInParallel files written in parallel. Each master task collects the galaxy/halo data from
 *  its slave tasks and writes it. If there is one output file it will be called galaxies.out otherwise there
 *  will be n files called galaxies.i where i < n.
 */
void output_galaxies(void)
{
	int i, j, ihalo, iprog, filenr, masterTask, lastTask, ngroups, maxlen, task;
	int left_to_send, p, send_this_turn, offset, nhaloout;
	float sigmaobs, mstarobs, sfrobs, scatter, amhalf;
	char buf[NSTRING];
	FILE *fd;
	MPI_Status status;

	//Structure that stores all halo/galaxy properties that will be saved
	struct haloout {
		float mh,mdot,mhp,mpd,mhh,amp,amh,r,c,l,ms,sfr,icm,mso,sfro,x,y,z,u,v,w;
		unsigned short t;
		long long hid, did, uid;
	} *hp, hsend;

#ifdef HDF5_SUPPORT
	int k;
	hid_t gal_file, hdf5_dataspace_memory, hdf5_dataspace_in_file, hdf5_dtype, hdf5_ftype, hdf5_filespace;
	hid_t hdf5_properties, hdf5_dataset;
	hid_t hdf5_att_space, hdf5_att_redshift, hdf5_att_lbox, hdf5_att_hubble;
	hid_t hdf5_att_omega0, hdf5_att_omegal, hdf5_att_omegab, hdf5_att_mmin;
	herr_t hdf5_status;
	hsize_t dims[1], maxdims[1], count[1], start[1];
	int sendsum = 0;
	double tmp;
	char path[NSTRING];
#endif

	//Print what is done...
	if (ThisTask == 0) printf("%s\n%s Writing galaxies in 1st universe to %d files...\n", All.fullline,All.startline,All.NumOutputFiles * All.Noutputredshifts);

	//Check if there are enough tasks to write all output files
	if (All.NTaskPerUniverse < All.NumOutputFiles && ThisTask == 0)
		endrun("Fatal error. Number of processors must be larger or equal than All.NumOutputFiles.");

	//Allocate CommBuffer
	CommBuffer = emalloc("CommBuffer", All.BufferSize * 1024 * 1024);

	//Get maximum number of haloes that fit in the CommBuffer
	maxlen = ((int) (All.BufferSize * 1024 * 1024)) / sizeof(struct halo);

	//Assign processors to output files
	distribute_file(All.NTaskPerUniverse, All.NumOutputFiles, 0, 0, All.NTaskPerUniverse - 1, &filenr, &masterTask, &lastTask);

	//Go through all output redshifts
	for (i = 0; i < All.Noutputredshifts; i++)
	{
		//If the output format has been set to HDF5
		if (All.OutputFormat == 2)
		{
			//If HDF5 libraries are set
#ifdef HDF5_SUPPORT
			if (All.NumOutputFiles > 1)
				sprintf(buf, "%s/galaxies/galaxies.S%d.%d.h5", All.OutputDir, Output_iscale[i], filenr);
			else
				sprintf(buf, "%s/galaxies/galaxies.S%d.h5", All.OutputDir, Output_iscale[i]);
			//If not return
#else
			if (ThisTask == 0) printf("%s Output format has been set to 2 but HDF5 support was not enabled.\n",All.startline);
			return;
#endif
		}
		//Otherwise use standard ascii files
		else
		{
			if (All.NumOutputFiles > 1)
				sprintf(buf, "%s/galaxies/galaxies.S%d.%d", All.OutputDir, Output_iscale[i], filenr);
			else
				sprintf(buf, "%s/galaxies/galaxies.S%d.out", All.OutputDir, Output_iscale[i]);
		}

		//Get number of groups
		ngroups = All.NumOutputFiles / All.NumFilesInParallel;
		if ((All.NumOutputFiles % All.NumFilesInParallel)) ngroups++;

		//For each group do...
		for (j = 0; j < ngroups; j++)
		{
			//This task will be processed now
			if ((filenr / All.NumFilesInParallel) == j && MasterTask == 0)
			{ //The masterTask opens the file
				if (ThisTask == masterTask)
				{ //If the hdf5 output format has been selected
					if (All.OutputFormat == 2)
					{ //Check if the libraries have been included
#ifdef HDF5_SUPPORT
						//Create hdf5 file
						gal_file = H5Fcreate(buf, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
						//Print to screen
						printf("%s Writing output file at z = %f: '%s' (file %d of %d)\n",All.startline, OutputRedshifts[i], buf, filenr+1, All.NumOutputFiles);
						//Specify data type
						hdf5_dtype  = H5Tcreate(H5T_COMPOUND, sizeof(struct haloout));
						k = 0;
						hdf5_status = H5Tinsert(hdf5_dtype, "Halo_mass",          HOFFSET(struct haloout, mh),   H5T_NATIVE_FLOAT); k+=4;
						hdf5_status = H5Tinsert(hdf5_dtype, "Halo_growth_rate",   HOFFSET(struct haloout, mdot), H5T_NATIVE_FLOAT); k+=4;
						hdf5_status = H5Tinsert(hdf5_dtype, "Halo_mass_peak",     HOFFSET(struct haloout, mhp),  H5T_NATIVE_FLOAT); k+=4;
						hdf5_status = H5Tinsert(hdf5_dtype, "Halo_growth_peak",   HOFFSET(struct haloout, mpd),  H5T_NATIVE_FLOAT); k+=4;
						hdf5_status = H5Tinsert(hdf5_dtype, "Halo_mass_host",     HOFFSET(struct haloout, mhh),  H5T_NATIVE_FLOAT); k+=4;
						hdf5_status = H5Tinsert(hdf5_dtype, "Scale_peak_mass",    HOFFSET(struct haloout, amp),  H5T_NATIVE_FLOAT); k+=4;
						hdf5_status = H5Tinsert(hdf5_dtype, "Scale_half_mass",    HOFFSET(struct haloout, amh),  H5T_NATIVE_FLOAT); k+=4;
						hdf5_status = H5Tinsert(hdf5_dtype, "Halo_radius",        HOFFSET(struct haloout, r),    H5T_NATIVE_FLOAT); k+=4;
						hdf5_status = H5Tinsert(hdf5_dtype, "Concentration",      HOFFSET(struct haloout, c),    H5T_NATIVE_FLOAT); k+=4;
						hdf5_status = H5Tinsert(hdf5_dtype, "Halo_spin",          HOFFSET(struct haloout, l),    H5T_NATIVE_FLOAT); k+=4;
						hdf5_status = H5Tinsert(hdf5_dtype, "Stellar_mass",       HOFFSET(struct haloout, ms),   H5T_NATIVE_FLOAT); k+=4;
						hdf5_status = H5Tinsert(hdf5_dtype, "SFR",                HOFFSET(struct haloout, sfr),  H5T_NATIVE_FLOAT); k+=4;
						hdf5_status = H5Tinsert(hdf5_dtype, "Intra_cluster_mass", HOFFSET(struct haloout, icm),  H5T_NATIVE_FLOAT); k+=4;
						hdf5_status = H5Tinsert(hdf5_dtype, "Stellar_mass_obs",   HOFFSET(struct haloout, mso),  H5T_NATIVE_FLOAT); k+=4;
						hdf5_status = H5Tinsert(hdf5_dtype, "SFR_obs",            HOFFSET(struct haloout, sfro), H5T_NATIVE_FLOAT); k+=4;
						hdf5_status = H5Tinsert(hdf5_dtype, "X_pos",              HOFFSET(struct haloout, x),    H5T_NATIVE_FLOAT); k+=4;
						hdf5_status = H5Tinsert(hdf5_dtype, "Y_pos",              HOFFSET(struct haloout, y),    H5T_NATIVE_FLOAT); k+=4;
						hdf5_status = H5Tinsert(hdf5_dtype, "Z_pos",              HOFFSET(struct haloout, z),    H5T_NATIVE_FLOAT); k+=4;
						hdf5_status = H5Tinsert(hdf5_dtype, "X_vel",              HOFFSET(struct haloout, u),    H5T_NATIVE_FLOAT); k+=4;
						hdf5_status = H5Tinsert(hdf5_dtype, "Y_vel",              HOFFSET(struct haloout, v),    H5T_NATIVE_FLOAT); k+=4;
						hdf5_status = H5Tinsert(hdf5_dtype, "Z_vel",              HOFFSET(struct haloout, w),    H5T_NATIVE_FLOAT); k+=4;
						hdf5_status = H5Tinsert(hdf5_dtype, "Type",               HOFFSET(struct haloout, t),    H5T_NATIVE_USHORT); k+=2;
						hdf5_status = H5Tinsert(hdf5_dtype, "Halo_ID",            HOFFSET(struct haloout, hid),  H5T_NATIVE_LLONG); k+=8;
						hdf5_status = H5Tinsert(hdf5_dtype, "Desc_ID",            HOFFSET(struct haloout, did),  H5T_NATIVE_LLONG); k+=8;
						hdf5_status = H5Tinsert(hdf5_dtype, "Up_ID",              HOFFSET(struct haloout, uid),  H5T_NATIVE_LLONG); k+=8;

						//Specify file type
						hdf5_ftype  = H5Tcreate(H5T_COMPOUND, k);
						k = 0;
						hdf5_status = H5Tinsert(hdf5_ftype, "Halo_mass",          k,  H5T_IEEE_F32LE); k+=4;
						hdf5_status = H5Tinsert(hdf5_ftype, "Halo_growth_rate",   k,  H5T_IEEE_F32LE); k+=4;
						hdf5_status = H5Tinsert(hdf5_ftype, "Halo_mass_peak",     k,  H5T_IEEE_F32LE); k+=4;
						hdf5_status = H5Tinsert(hdf5_ftype, "Halo_growth_peak",   k,  H5T_IEEE_F32LE); k+=4;
						hdf5_status = H5Tinsert(hdf5_ftype, "Halo_mass_host",     k,  H5T_IEEE_F32LE); k+=4;
						hdf5_status = H5Tinsert(hdf5_ftype, "Scale_peak_mass",    k,  H5T_IEEE_F32LE); k+=4;
						hdf5_status = H5Tinsert(hdf5_ftype, "Scale_half_mass",    k,  H5T_IEEE_F32LE); k+=4;
						hdf5_status = H5Tinsert(hdf5_ftype, "Halo_radius",        k,  H5T_IEEE_F32LE); k+=4;
						hdf5_status = H5Tinsert(hdf5_ftype, "Concentration",      k,  H5T_IEEE_F32LE); k+=4;
						hdf5_status = H5Tinsert(hdf5_ftype, "Halo_spin",          k,  H5T_IEEE_F32LE); k+=4;
						hdf5_status = H5Tinsert(hdf5_ftype, "Stellar_mass",       k,  H5T_IEEE_F32LE); k+=4;
						hdf5_status = H5Tinsert(hdf5_ftype, "SFR",                k,  H5T_IEEE_F32LE); k+=4;
						hdf5_status = H5Tinsert(hdf5_ftype, "Intra_cluster_mass", k,  H5T_IEEE_F32LE); k+=4;
						hdf5_status = H5Tinsert(hdf5_ftype, "Stellar_mass_obs",   k,  H5T_IEEE_F32LE); k+=4;
						hdf5_status = H5Tinsert(hdf5_ftype, "SFR_obs",            k,  H5T_IEEE_F32LE); k+=4;
						hdf5_status = H5Tinsert(hdf5_ftype, "X_pos",              k,  H5T_IEEE_F32LE); k+=4;
						hdf5_status = H5Tinsert(hdf5_ftype, "Y_pos",              k,  H5T_IEEE_F32LE); k+=4;
						hdf5_status = H5Tinsert(hdf5_ftype, "Z_pos",              k,  H5T_IEEE_F32LE); k+=4;
						hdf5_status = H5Tinsert(hdf5_ftype, "X_vel",              k,  H5T_IEEE_F32LE); k+=4;
						hdf5_status = H5Tinsert(hdf5_ftype, "Y_vel",              k,  H5T_IEEE_F32LE); k+=4;
						hdf5_status = H5Tinsert(hdf5_ftype, "Z_vel",              k,  H5T_IEEE_F32LE); k+=4;
						hdf5_status = H5Tinsert(hdf5_ftype, "Type",               k,  H5T_STD_U16LE); k+=2;
						hdf5_status = H5Tinsert(hdf5_ftype, "Halo_ID",            k,  H5T_STD_I64LE); k+=8;
						hdf5_status = H5Tinsert(hdf5_ftype, "Desc_ID",            k,  H5T_STD_I64LE); k+=8;
						hdf5_status = H5Tinsert(hdf5_ftype, "Up_ID",              k,  H5T_STD_I64LE); k+=8;
#endif
					}
					//Otherwise open standard ascii file
					else
					{ //Open file
						if (!(fd = fopen(buf, "w")))
						{ //If not possible print to screen and abort
							printf("%s Can't open file `%s' for writing output.\n",All.startline,buf);
							endrun("file open error");
						}
						//If file can be opened print to screen
						printf("%s Writing output file at z = %f: '%s' (file %d of %d)\n",All.startline, OutputRedshifts[i], buf, filenr+1, All.NumOutputFiles);
						//Write ascii header
						task = 0;
						fprintf(fd,"#Halo_mass(%d)",task); task++;
						fprintf(fd," Halo_growth_rate(%d)",task); task++;
						fprintf(fd," Halo_mass_peak(%d)",task); task++;
						fprintf(fd," Halo_growth_peak(%d)",task); task++;
						fprintf(fd," Halo_mass_host(%d)",task); task++;
						fprintf(fd," Scale_peak_mass(%d)",task); task++;
						fprintf(fd," Scale_half_mass(%d)",task); task++;
						fprintf(fd," Halo_radius(%d)",task); task++;
						fprintf(fd," Concentration(%d)",task); task++;
						fprintf(fd," Halo_spin(%d)",task); task++;
						fprintf(fd," Stellar_mass(%d)",task); task++;
						fprintf(fd," SFR(%d)",task); task++;
						fprintf(fd," Intra_cluster_mass(%d)",task); task++;
						fprintf(fd," Stellar_mass_obs(%d)",task); task++;
						fprintf(fd," SFR_obs(%d)",task); task++;
						fprintf(fd," X_pos(%d)",task); task++;
						fprintf(fd," Y_pos(%d)",task); task++;
						fprintf(fd," Z_pos(%d)",task); task++;
						fprintf(fd," X_vel(%d)",task); task++;
						fprintf(fd," Y_vel(%d)",task); task++;
						fprintf(fd," Z_vel(%d)",task); task++;
						fprintf(fd," Type(%d)",task); task++;
						fprintf(fd," Halo_ID(%d)",task); task++;
						fprintf(fd," Desc_ID(%d)",task); task++;
						fprintf(fd," Up_ID(%d)",task); task++;
						fprintf(fd,"\n");
						fprintf(fd,"#Emerge Branch: %s\n",BRANCH);
						fprintf(fd,"#Emerge Version: %s\n",VERSION);
						fprintf(fd,"#Model Name: %s\n",All.model_name);
						fprintf(fd,"#Scale Factor: %f\n",ScaleFactor[Output_iscale[i]]);
						fprintf(fd,"#Box Size: %f (Mpc)\n",All.Lbox*All.x_unit);
						fprintf(fd,"#Hubble Parameter: %f\n",All.h_100);
						fprintf(fd,"#Omega_0: %f\n",All.Omega_0);
						fprintf(fd,"#Omega_Lambda: %f\n",All.Omega_Lambda_0);
						fprintf(fd,"#Omega_Baryon: %f\n",All.Omega_Baryon_0);
						fprintf(fd,"#Minimum Stellar Mass: %f (log Msun)\n",log10(All.minmass*All.m_unit));
						fprintf(fd,"#Halo_mass: Current virial mass of the halo (log Msun).\n");
						fprintf(fd,"#Halo_growth_rate: Growth rate of the halo (Msun/yr)\n");
						fprintf(fd,"#Halo_mass_peak: Peak virial mass of the halo through its history (log Msun).\n");
						fprintf(fd,"#Halo_growth_peak: Growth rate of the halo when peak mass is reached (Msun/yr).\n");
						fprintf(fd,"#Halo_mass_host: Current virial mass of the host halo (log Msun).\n");
						fprintf(fd,"#Scale_peak_mass: Scale Factor when halo had its peak mass.\n");
						fprintf(fd,"#Scale_half_mass: Scale Factor when halo first had half of its peak mass.\n");
						fprintf(fd,"#Halo_radius: Current virial radius of the halo (kpc).\n");
						fprintf(fd,"#Concentration: Virial radius divided by halo scale length.\n");
						fprintf(fd,"#Halo_spin: Halo spin parameter.\n");
						fprintf(fd,"#Stellar_mass: Intrinsic stellar mass of the galaxy (log Msun).\n");
						fprintf(fd,"#SFR: Star formation rate of the galaxy (Msun/yr).\n");
						fprintf(fd,"#Intra_cluster_mass: Stellar mass in the halo (log Msun).\n");
						fprintf(fd,"#Stellar_mass_obs: Observed stellar mass including observation error (log Msun).\n");
						fprintf(fd,"#SFR_obs: Observed star formation rate including observation error (Msun/yr).\n");
						fprintf(fd,"#X_pos: Comoving galaxy position along the X axis (Mpc).\n");
						fprintf(fd,"#Y_pos: Comoving galaxy position along the Y axis (Mpc).\n");
						fprintf(fd,"#Z_pos: Comoving galaxy position along the Z axis (Mpc).\n");
						fprintf(fd,"#X_vel: Peculiar/physical velocity of the galaxy along the X axis (km/s).\n");
						fprintf(fd,"#Y_vel: Peculiar/physical velocity of the galaxy along the Y axis (km/s).\n");
						fprintf(fd,"#Z_vel: Peculiar/physical velocity of the galaxy along the Z axis (km/s).\n");
						fprintf(fd,"#Type: Whether the galaxy is a central (0), satellite (1), or orphan (2).\n");
						fprintf(fd,"#Halo_ID: Halo ID that corresponds to the ID in the original halo catalogue.\n");
						fprintf(fd,"#Desc_ID: ID of the descendant halo. If there is none then Desc_ID = -1.\n");
						fprintf(fd,"#Up_ID: ID of the most massive halo. If this is a main halo then Up_ID = -1.\n");
					}
				}

				//Determine number of systems that will be printed on each task (depending on redshift and mass)
				nhaloout = 0;
				for (ihalo = 0; ihalo < Nhalos; ihalo++)
					if (H[ihalo].iscale == Output_iscale[i] && H[ihalo].gone == 0 && H[ihalo].mstar >= All.minmass) nhaloout++;

#ifdef HDF5_SUPPORT
				if (ThisTask == masterTask && All.OutputFormat == 2)
				{
					//Set the path name for the galaxy data set
					sprintf(path,"/Galaxies");
					//Initial size is the number of haloes on the master task
					dims[0] = nhaloout;
					//Allow for variable total size
					maxdims[0] = H5S_UNLIMITED;
					//Create the data space
					hdf5_dataspace_in_file = H5Screate_simple(1, dims, maxdims);
					//Create the data properties
					hdf5_properties = H5Pcreate(H5P_DATASET_CREATE);
					hdf5_status = H5Pset_chunk(hdf5_properties, 1, dims);   // set chunk size
					hdf5_status = H5Pset_shuffle(hdf5_properties);          // reshuffle bytes to get better compression ratio
					hdf5_status = H5Pset_deflate(hdf5_properties, 9);       // gzip compression level 9
					hdf5_status = H5Pset_fletcher32(hdf5_properties);       // Fletcher32 checksum on dataset
					//If filters could be set use them to create the data set
					if (H5Pall_filters_avail(hdf5_properties))
						hdf5_dataset = H5Dcreate(gal_file, path, hdf5_ftype, hdf5_dataspace_in_file, H5P_DEFAULT, hdf5_properties, H5P_DEFAULT);
					//Otherwise create the default data set
					else
						hdf5_dataset = H5Dcreate(gal_file, path, hdf5_ftype, hdf5_dataspace_in_file, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
					//Initialise the number of haloes in the set to zero
					sendsum = 0;
					//Define Scalar for Attributes
					hdf5_att_space = H5Screate(H5S_SCALAR);
					//Write Redshift as attribute to dataset
					hdf5_att_redshift = H5Acreate(hdf5_dataset, "Scale Factor", H5T_NATIVE_FLOAT, hdf5_att_space, H5P_DEFAULT, H5P_DEFAULT);
					H5Awrite(hdf5_att_redshift, H5T_IEEE_F32LE, &ScaleFactor[Output_iscale[i]]);
					H5Aclose(hdf5_att_redshift);
					//Write Box Size as attribute to dataset
					tmp               = All.Lbox*All.x_unit;
					hdf5_att_lbox     = H5Acreate(hdf5_dataset, "Box Size", H5T_NATIVE_DOUBLE, hdf5_att_space, H5P_DEFAULT, H5P_DEFAULT);
					H5Awrite(hdf5_att_lbox, H5T_IEEE_F64LE, &tmp);
					H5Aclose(hdf5_att_lbox);
					//Write Hubble parameter as attribute to dataset
					hdf5_att_hubble     = H5Acreate(hdf5_dataset, "Hubble Parameter", H5T_NATIVE_DOUBLE, hdf5_att_space, H5P_DEFAULT, H5P_DEFAULT);
					H5Awrite(hdf5_att_hubble, H5T_IEEE_F64LE, &All.h_100);
					H5Aclose(hdf5_att_hubble);
					//Write Omega_0 as attribute to dataset
					hdf5_att_omega0     = H5Acreate(hdf5_dataset, "Omega_0", H5T_NATIVE_DOUBLE, hdf5_att_space, H5P_DEFAULT, H5P_DEFAULT);
					H5Awrite(hdf5_att_omega0, H5T_IEEE_F64LE, &All.Omega_0);
					H5Aclose(hdf5_att_omega0);
					//Write Omega_Lambda as attribute to dataset
					hdf5_att_omegal     = H5Acreate(hdf5_dataset, "Omega_Lambda", H5T_NATIVE_DOUBLE, hdf5_att_space, H5P_DEFAULT, H5P_DEFAULT);
					H5Awrite(hdf5_att_omegal, H5T_IEEE_F64LE, &All.Omega_Lambda_0);
					H5Aclose(hdf5_att_omegal);
					//Write Omega_Baryon as attribute to dataset
					hdf5_att_omegab     = H5Acreate(hdf5_dataset, "Omega_Baryon", H5T_NATIVE_DOUBLE, hdf5_att_space, H5P_DEFAULT, H5P_DEFAULT);
					H5Awrite(hdf5_att_omegab, H5T_IEEE_F64LE, &All.Omega_Baryon_0);
					H5Aclose(hdf5_att_omegab);
					//Write minimum stellar mass as attribute to dataset
					tmp               = log10(All.minmass*All.m_unit);
					hdf5_att_mmin     = H5Acreate(hdf5_dataset, "Minimum Stellar Mass", H5T_NATIVE_DOUBLE, hdf5_att_space, H5P_DEFAULT, H5P_DEFAULT);
					H5Awrite(hdf5_att_mmin, H5T_IEEE_F64LE, &tmp);
					H5Aclose(hdf5_att_mmin);
				}
#endif

				//Now go through all tasks that are member of this group
				for (task = masterTask, offset = 0; task <= lastTask; task++)
				{
					//if this task is processed
					if (task == ThisTask)
					{
						//We need to get those haloes
						left_to_send = nhaloout;

						//Tell this to the other tasks
						for (p = masterTask; p <= lastTask; p++)
							if (p != ThisTask)
								MPI_Send(&left_to_send, 1, MPI_INT, p, TAG_NHALOS, MPI_COMM_WORLD);
					}
					//Reveive the number of haloes from active task
					else MPI_Recv(&left_to_send, 1, MPI_INT, task, TAG_NHALOS, MPI_COMM_WORLD, &status);

					//Now send the haloes in increments
					while (left_to_send > 0)
					{
						//Compute how many haloes we can send this turn
						send_this_turn = left_to_send;
						//If they do not all fit into the CommBuffer set to the maximum
						if (send_this_turn > maxlen) send_this_turn = maxlen;

						//Sending task
						if (ThisTask == task)
						{
							//Use the haloout struct for the CommBuffer
							hp = (struct haloout *) CommBuffer;
							//Initialise the number of systems sent to zero
							ihalo = 0;
							//Go through all systems until the maximum that can be sent
							while (ihalo < send_this_turn)
							{   //If we are at the right scale factor and the halo is not gone send it to the master
								if (H[offset].iscale == Output_iscale[i] && H[offset].gone == 0 && H[offset].mstar >= All.minmass)
								{
									//Set the scatter in stellar mass
									scatter = get_gaussian_random_number(offset);
									//Compute the observational error and the observed stellar mass
									sigmaobs = All.obssigma0+(1./H[offset].a-1.)*All.obssigmaz;
									if (1./H[offset].a-1. > ZMAX_SMFERR) sigmaobs = All.obssigma0+(ZMAX_SMFERR)*All.obssigmaz;
									mstarobs = H[offset].mstar * pow(10.,sigmaobs*scatter);
									//Set the scatter in the star formation rate
									scatter = get_gaussian_random_number(offset + RANDOM_NUMBER_TABLE/2);
									//Compute the observed star formation rate
									sfrobs = H[offset].sfr * pow(10.,sigmaobs*scatter);

									//Find scale factor when virial mass was half of maximum value
									amhalf = ScaleFactor[H[H[offset].impeak].iscale];
									iprog = H[H[offset].impeak].iprog;
									if (iprog >=0 )
									{
										while (H[iprog].mvir > 0.5 * H[H[offset].impeak].mvir)
										{
											if (H[iprog].iprog < 0) break;
											iprog = H[iprog].iprog;
										}
										amhalf = ScaleFactor[H[iprog].iscale];
									}

									//Set all halo/galaxy properties
									hsend.mh   = log10(H[offset].mvir*All.m_unit);
									hsend.mdot = H[offset].mdotbary*All.m_unit/All.t_unit/All.f_baryon;
									hsend.mhp  = log10(H[H[offset].impeak].mvir*All.m_unit);
									hsend.mpd  = H[H[offset].impeak].mdotbary*All.m_unit/All.t_unit/All.f_baryon;
#ifdef COMPUTE_ICM
									hsend.mhh  = log10(H[H[offset].ihost].mvir*All.m_unit);
#else
									hsend.mhh  = 0.0;
#endif
									hsend.amp  = ScaleFactor[H[H[offset].impeak].iscale];
									hsend.amh  = amhalf;
									hsend.r    = H[offset].rvir*All.x_unit*1.e3;
									hsend.c    = H[offset].c;
									hsend.l    = H[offset].lambda;
									hsend.ms   = log10(H[offset].mstar*All.m_unit);
									hsend.sfr  = H[offset].sfr*All.m_unit/All.t_unit;
									hsend.icm  = log10(H[offset].icm*All.m_unit);
									hsend.mso  = log10(mstarobs*All.m_unit);
									hsend.sfro = sfrobs*All.m_unit/All.t_unit;
									hsend.x    = H[offset].pos[0];
									hsend.y    = H[offset].pos[1];
									hsend.z    = H[offset].pos[2];
									hsend.u    = H[offset].vel[0];
									hsend.v    = H[offset].vel[1];
									hsend.w    = H[offset].vel[2];
									hsend.t    = H[offset].type;
									hsend.hid  = (long long)(H[offset].haloid)-1;
									hsend.did  = (long long)(H[offset].descid)-1;
									hsend.uid  = (long long)(H[offset].upid)-1;

									//Write this system to the CommBuffer
									*hp++ = hsend;
									//Increment the number of systems that are being sent
									ihalo++;
								}
								//Increment the halo index
								offset++;
							}
						}

						//Receive haloes
						if (ThisTask == masterTask && task != masterTask) MPI_Recv(CommBuffer, (sizeof(struct halo)) * send_this_turn, MPI_BYTE, task, TAG_HDATA, MPI_COMM_WORLD, &status);

						//Send haloes
						if (ThisTask != masterTask && task == ThisTask) MPI_Ssend(CommBuffer, (sizeof(struct halo)) * send_this_turn, MPI_BYTE, masterTask, TAG_HDATA, MPI_COMM_WORLD);

						//Collect haloes and write to file
						if (ThisTask == masterTask)
						{

							//If the output format is hdf5
							if (All.OutputFormat == 2)
							{
#ifdef HDF5_SUPPORT
								//Set the starting index and the count
								start[0] = sendsum;
								count[0] = send_this_turn;
								//Add the number of haloes that are sent this turn to the total number and set the size to the total
								sendsum += send_this_turn;
								dims[0]  = sendsum;
								//If the total size is larger than the initial set size for the data set extend it
								if (sendsum > nhaloout) hdf5_status = H5Dset_extent(hdf5_dataset, dims);
								//Get the data space of the data set
								hdf5_filespace = H5Dget_space(hdf5_dataset);
								//Select the hyperslab that corresponds to the current offset
								hdf5_status = H5Sselect_hyperslab(hdf5_filespace, H5S_SELECT_SET, start, NULL, count, NULL);
								//Set the size to the number of systems sent this turn
								dims[0] = send_this_turn;
								//Create a data space for the memory
								hdf5_dataspace_memory = H5Screate_simple(1, dims, NULL);
								//Write the data set
								hdf5_status = H5Dwrite(hdf5_dataset, hdf5_dtype, hdf5_dataspace_memory, hdf5_filespace, H5P_DEFAULT, CommBuffer);
								//Close the data space of the data set
								H5Sclose(hdf5_filespace);
								//Close the data space for the memory
								H5Sclose(hdf5_dataspace_memory);
#endif
							}  //End HDF5 output
							//Otherwise simply write the properties to an ascii file
							else
							{   //Use the haloout struct for the CommBuffer
								hp = (struct haloout *) CommBuffer;
								//Go through all systems that were sent this turn
								for (ihalo = 0; ihalo < send_this_turn; ihalo++)
								{   //Write a new line in the ascii file
									fprintf(fd,"%f %e %f %e %f %f %f %f %f %f %f %e %f %f %e %f %f %f %f %f %f %hu %lld %lld %lld",
									        hp[ihalo].mh,
									        hp[ihalo].mdot,
									        hp[ihalo].mhp,
									        hp[ihalo].mpd,
									        hp[ihalo].mhh,
									        hp[ihalo].amp,
									        hp[ihalo].amh,
									        hp[ihalo].r,
									        hp[ihalo].c,
									        hp[ihalo].l,
									        hp[ihalo].ms,
									        hp[ihalo].sfr,
									        hp[ihalo].icm,
									        hp[ihalo].mso,
									        hp[ihalo].sfro,
									        hp[ihalo].x,
									        hp[ihalo].y,
									        hp[ihalo].z,
									        hp[ihalo].u,
									        hp[ihalo].v,
									        hp[ihalo].w,
									        hp[ihalo].t,
									        hp[ihalo].hid,
									        hp[ihalo].did,
									        hp[ihalo].uid);
									fprintf(fd,"\n");
								}  //End loop through all haloes
							}  //End ascii output
						}  //End master task (writer task)
						//Decrease the number of haloes left to send by the number sent this turn
						left_to_send -= send_this_turn;
					}  //End halo increments
				}  //End loop through all tasks


				//If this is the master task (writer task) close the file (and other objects)
				if (ThisTask == masterTask)
				{
					//For HDF5 output
					if (All.OutputFormat == 2)
					{ //Check if libraries have been included
#ifdef HDF5_SUPPORT
						if (nhaloout > 0)
						{ // Close the dataset
							hdf5_status = H5Dclose(hdf5_dataset);
							// Close the dataset
							hdf5_status = H5Pclose(hdf5_properties);
							// Close the data space
							hdf5_status = H5Sclose(hdf5_dataspace_in_file);
						}
						// Close the file type
						hdf5_status = H5Tclose(hdf5_ftype);
						// Close the data type
						hdf5_status = H5Tclose(hdf5_dtype);
						//Close the file and print to screen if writing was not  successful
						hdf5_status = H5Fclose(gal_file);
						if (hdf5_status != 0) printf("%s Task %5d failed to write to file %s\n",All.startline,ThisTask,buf);
#endif
					}//End HDF5
					//Otherwise close the ascii file
					else
					{ //Close the file
						fclose(fd);
					}//End ascii
				}//End master/writer task

			}//End tasks in the group

		}//End output redshifts

		//Set MPI barrier so that only one task writes to a given file
		MPI_Barrier(MPI_COMM_WORLD);
	}//End loop through all groups

	//Free the CommBuffer
	efree(CommBuffer);

}
#endif


#ifdef WRITE_HALO_CATALOG
/*! \brief This function collects all halo data and prints it
 *
 *  This function collects the halo data from all tasks and writes it to All.NumOutputFiles files. There
 *  are All.NumFilesInParallel files written in parallel. Each master task collects the halo data from
 *  its slave tasks and writes it. If there is one output file it will be called haloes.out otherwise there
 *  will be n files called haloes.i where i < n.
 */
void output_halos(void)
{
	int i, j, k, ihalo, iprog, filenr, masterTask, lastTask, ngroups, maxlen, task;
	int left_to_send, p, send_this_turn, offset, nhaloout;
	double amhalf;
	char buf[NSTRING];
	FILE *fd;
	MPI_Status status;

	//Structure that stores all halo/galaxy properties that will be saved
	struct haloout {
		float mvir,mdot,mpeak,mpeakdot,ampeak,amhalf,rvir,c,l,x,y,z,u,v,w;
		unsigned short t;
		long long hid, did, uid;
	} *hp, hsend;

#ifdef HDF5_SUPPORT
	hid_t halo_file, hdf5_dataspace_memory, hdf5_dataspace_in_file, hdf5_dtype, hdf5_ftype, hdf5_filespace;
	hid_t hdf5_properties, hdf5_dataset;
	hid_t hdf5_att_space, hdf5_att_redshift, hdf5_att_lbox, hdf5_att_hubble;
	hid_t hdf5_att_omega0, hdf5_att_omegal, hdf5_att_omegab;
	herr_t hdf5_status;
	hsize_t dims[1], maxdims[1], count[1], start[1];
	int sendsum = 0;
	double tmp;
	char path[NSTRING];
#endif

	//Print what is done...
	if (ThisTask == 0) printf("%s\n%s Writing haloes in 1st universe to %d files...\n", All.fullline,All.startline,All.NumOutputFiles * All.Noutputredshifts);

	//Check if there are enough tasks to write all output files
	if (All.NTaskPerUniverse < All.NumOutputFiles && ThisTask == 0)
		endrun("Fatal error. Number of processors must be larger or equal than All.NumOutputFiles.");

	//Allocate CommBuffer
	CommBuffer = emalloc("CommBuffer", All.BufferSize * 1024 * 1024);

	//Get maximum number of haloes that fit in the CommBuffer
	maxlen = ((int) (All.BufferSize * 1024 * 1024)) / sizeof(struct halo);

	//Assign processors to output files
	distribute_file(All.NTaskPerUniverse, All.NumOutputFiles, 0, 0, All.NTaskPerUniverse - 1, &filenr, &masterTask, &lastTask);

	//Go through all output redshifts
	for (i = 0; i < All.Noutputredshifts; i++)
	{
		//If the output format has been set to HDF5
		if (All.OutputFormat == 2)
		{
			//If HDF5 libraries are set
#ifdef HDF5_SUPPORT
			if (All.NumOutputFiles > 1)
				sprintf(buf, "%s/haloes/haloes.S%d.%d.h5", All.OutputDir, Output_iscale[i], filenr);
			else
				sprintf(buf, "%s/haloes/haloes.S%d.h5", All.OutputDir, Output_iscale[i]);
			//If not return
#else
			if (ThisTask == 0) printf("%s Output format has been set to 2 but HDF5 support was not enabled.\n",All.startline);
			return;
#endif
		}
		//Otherwise use standard ascii files
		else
		{
			if (All.NumOutputFiles > 1)
				sprintf(buf, "%s/haloes/haloes.S%d.%d", All.OutputDir, Output_iscale[i], filenr);
			else
				sprintf(buf, "%s/haloes/haloes.S%d.out", All.OutputDir, Output_iscale[i]);
		}

		//Get number of groups
		ngroups = All.NumOutputFiles / All.NumFilesInParallel;
		if ((All.NumOutputFiles % All.NumFilesInParallel)) ngroups++;

		//For each group do...
		for (j = 0; j < ngroups; j++)
		{
			//This task will be processed now
			if ((filenr / All.NumFilesInParallel) == j && MasterTask == 0)
			{ //The masterTask opens the file
				if (ThisTask == masterTask)
				{ //If the hdf5 output format has been selected
					if (All.OutputFormat == 2)
					{ //Check if the libraries have been included
#ifdef HDF5_SUPPORT
						//Create hdf5 file
						halo_file = H5Fcreate(buf, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
						//Print to screen
						printf("%s Writing output file at z = %f: '%s' (file %d of %d)\n",All.startline, OutputRedshifts[i], buf, filenr+1, All.NumOutputFiles);
						//Specify data type
						hdf5_dtype  = H5Tcreate(H5T_COMPOUND, sizeof(struct haloout));
						k = 0;
						hdf5_status = H5Tinsert(hdf5_dtype, "Halo_mass",          HOFFSET(struct haloout, mvir),     H5T_NATIVE_FLOAT); k+=4;
						hdf5_status = H5Tinsert(hdf5_dtype, "Halo_growth_rate",   HOFFSET(struct haloout, mdot),     H5T_NATIVE_FLOAT); k+=4;
						hdf5_status = H5Tinsert(hdf5_dtype, "Halo_mass_peak",     HOFFSET(struct haloout, mpeak),    H5T_NATIVE_FLOAT); k+=4;
						hdf5_status = H5Tinsert(hdf5_dtype, "Halo_growth_peak",   HOFFSET(struct haloout, mpeakdot), H5T_NATIVE_FLOAT); k+=4;
						hdf5_status = H5Tinsert(hdf5_dtype, "Scale_peak_mass",    HOFFSET(struct haloout, ampeak),   H5T_NATIVE_FLOAT); k+=4;
						hdf5_status = H5Tinsert(hdf5_dtype, "Scale_half_mass",    HOFFSET(struct haloout, amhalf),   H5T_NATIVE_FLOAT); k+=4;
						hdf5_status = H5Tinsert(hdf5_dtype, "Halo_radius",        HOFFSET(struct haloout, rvir),     H5T_NATIVE_FLOAT); k+=4;
						hdf5_status = H5Tinsert(hdf5_dtype, "Concentration",      HOFFSET(struct haloout, c),        H5T_NATIVE_FLOAT); k+=4;
						hdf5_status = H5Tinsert(hdf5_dtype, "Halo_spin",          HOFFSET(struct haloout, l),        H5T_NATIVE_FLOAT); k+=4;
						hdf5_status = H5Tinsert(hdf5_dtype, "X_pos",              HOFFSET(struct haloout, x),        H5T_NATIVE_FLOAT); k+=4;
						hdf5_status = H5Tinsert(hdf5_dtype, "Y_pos",              HOFFSET(struct haloout, y),        H5T_NATIVE_FLOAT); k+=4;
						hdf5_status = H5Tinsert(hdf5_dtype, "Z_pos",              HOFFSET(struct haloout, z),        H5T_NATIVE_FLOAT); k+=4;
						hdf5_status = H5Tinsert(hdf5_dtype, "X_vel",              HOFFSET(struct haloout, u),        H5T_NATIVE_FLOAT); k+=4;
						hdf5_status = H5Tinsert(hdf5_dtype, "Y_vel",              HOFFSET(struct haloout, v),        H5T_NATIVE_FLOAT); k+=4;
						hdf5_status = H5Tinsert(hdf5_dtype, "Z_vel",              HOFFSET(struct haloout, w),        H5T_NATIVE_FLOAT); k+=4;
						hdf5_status = H5Tinsert(hdf5_dtype, "Type",               HOFFSET(struct haloout, t),        H5T_NATIVE_USHORT); k+=2;
						hdf5_status = H5Tinsert(hdf5_dtype, "Halo_ID",            HOFFSET(struct haloout, hid),      H5T_NATIVE_LLONG); k+=8;
						hdf5_status = H5Tinsert(hdf5_dtype, "Desc_ID",            HOFFSET(struct haloout, did),      H5T_NATIVE_LLONG); k+=8;
						hdf5_status = H5Tinsert(hdf5_dtype, "Up_ID",              HOFFSET(struct haloout, uid),      H5T_NATIVE_LLONG); k+=8;
						//Specify file type
						hdf5_ftype  = H5Tcreate(H5T_COMPOUND, k);
						k = 0;
						hdf5_status = H5Tinsert(hdf5_ftype, "Halo_mass",          k,  H5T_IEEE_F32LE); k+=4;
						hdf5_status = H5Tinsert(hdf5_ftype, "Halo_growth_rate",   k,  H5T_IEEE_F32LE); k+=4;
						hdf5_status = H5Tinsert(hdf5_ftype, "Halo_mass_peak",     k,  H5T_IEEE_F32LE); k+=4;
						hdf5_status = H5Tinsert(hdf5_ftype, "Halo_growth_peak",   k,  H5T_IEEE_F32LE); k+=4;
						hdf5_status = H5Tinsert(hdf5_ftype, "Scale_peak_mass",    k,  H5T_IEEE_F32LE); k+=4;
						hdf5_status = H5Tinsert(hdf5_ftype, "Scale_half_mass",    k,  H5T_IEEE_F32LE); k+=4;
						hdf5_status = H5Tinsert(hdf5_ftype, "Halo_radius",        k,  H5T_IEEE_F32LE); k+=4;
						hdf5_status = H5Tinsert(hdf5_ftype, "Concentration",      k,  H5T_IEEE_F32LE); k+=4;
						hdf5_status = H5Tinsert(hdf5_ftype, "Halo_spin",          k,  H5T_IEEE_F32LE); k+=4;
						hdf5_status = H5Tinsert(hdf5_ftype, "X_pos",              k,  H5T_IEEE_F32LE); k+=4;
						hdf5_status = H5Tinsert(hdf5_ftype, "Y_pos",              k,  H5T_IEEE_F32LE); k+=4;
						hdf5_status = H5Tinsert(hdf5_ftype, "Z_pos",              k,  H5T_IEEE_F32LE); k+=4;
						hdf5_status = H5Tinsert(hdf5_ftype, "X_vel",              k,  H5T_IEEE_F32LE); k+=4;
						hdf5_status = H5Tinsert(hdf5_ftype, "Y_vel",              k,  H5T_IEEE_F32LE); k+=4;
						hdf5_status = H5Tinsert(hdf5_ftype, "Z_vel",              k,  H5T_IEEE_F32LE); k+=4;
						hdf5_status = H5Tinsert(hdf5_ftype, "Type",               k,  H5T_STD_U16LE); k+=2;
						hdf5_status = H5Tinsert(hdf5_ftype, "Halo_ID",            k,  H5T_STD_I64LE); k+=8;
						hdf5_status = H5Tinsert(hdf5_ftype, "Desc_ID",            k,  H5T_STD_I64LE); k+=8;
						hdf5_status = H5Tinsert(hdf5_ftype, "Up_ID",              k,  H5T_STD_I64LE); k+=8;
#endif
					}
					//Otherwise open standard ascii file
					else
					{ //Open file
						if (!(fd = fopen(buf, "w")))
						{ //If not possible print to screen and abort
							printf("%s Can't open file `%s' for writing output.\n",All.startline,buf);
							endrun("file open error");
						}
						//If file can be opened print to screen
						printf("%s Writing output file at z = %f: '%s' (file %d of %d)\n",All.startline, OutputRedshifts[i], buf, filenr+1, All.NumOutputFiles);
						//Write ascii header
						task = 0;
						fprintf(fd,"#Halo_mass(%d)",task); task++;
						fprintf(fd," Halo_growth_rate(%d)",task); task++;
						fprintf(fd," Halo_mass_peak(%d)",task); task++;
						fprintf(fd," Halo_growth_peak(%d)",task); task++;
						fprintf(fd," Scale_peak_mass(%d)",task); task++;
						fprintf(fd," Scale_half_mass(%d)",task); task++;
						fprintf(fd," Halo_radius(%d)",task); task++;
						fprintf(fd," Concentration(%d)",task); task++;
						fprintf(fd," Halo_spin(%d)",task); task++;
						fprintf(fd," X_pos(%d)",task); task++;
						fprintf(fd," Y_pos(%d)",task); task++;
						fprintf(fd," Z_pos(%d)",task); task++;
						fprintf(fd," X_vel(%d)",task); task++;
						fprintf(fd," Y_vel(%d)",task); task++;
						fprintf(fd," Z_vel(%d)",task); task++;
						fprintf(fd," Type(%d)",task); task++;
						fprintf(fd," Halo_ID(%d)",task); task++;
						fprintf(fd," Desc_ID(%d)",task); task++;
						fprintf(fd," Up_ID(%d)",task); task++;
						fprintf(fd,"\n");
						fprintf(fd,"#Emerge Branch: %s\n",BRANCH);
						fprintf(fd,"#Emerge Version: %s\n",VERSION);
						fprintf(fd,"#Model Name: %s\n",All.model_name);
						fprintf(fd,"#Scale Factor: %f\n",ScaleFactor[Output_iscale[i]]);
						fprintf(fd,"#Box Size: %f (Mpc)\n",All.Lbox*All.x_unit);
						fprintf(fd,"#Hubble Parameter: %f\n",All.h_100);
						fprintf(fd,"#Omega_0: %f\n",All.Omega_0);
						fprintf(fd,"#Omega_Lambda: %f\n",All.Omega_Lambda_0);
						fprintf(fd,"#Omega_Baryon: %f\n",All.Omega_Baryon_0);
						fprintf(fd,"#Halo_mass: Current virial mass of the halo (log Msun).\n");
						fprintf(fd,"#Halo_growth_rate: Growth rate of the halo (Msun/yr)\n");
						fprintf(fd,"#Halo_mass_peak: Peak virial mass of the halo through its history (log Msun).\n");
						fprintf(fd,"#Halo_growth_peak: Growth rate of the halo when peak mass is reached (Msun/yr).\n");
						fprintf(fd,"#Scale_peak_mass: Scale Factor when halo had its peak mass.\n");
						fprintf(fd,"#Scale_half_mass: Scale Factor when halo first had half of its peak mass.\n");
						fprintf(fd,"#Halo_radius: Current virial radius of the halo (kpc).\n");
						fprintf(fd,"#Concentration: Virial radius divided by halo scale length.\n");
						fprintf(fd,"#Halo_spin: Halo spin parameter.\n");
						fprintf(fd,"#X_pos: Comoving halo position along the X axis (Mpc).\n");
						fprintf(fd,"#Y_pos: Comoving halo position along the Y axis (Mpc).\n");
						fprintf(fd,"#Z_pos: Comoving halo position along the Z axis (Mpc).\n");
						fprintf(fd,"#X_vel: Peculiar/physical velocity of the halo along the X axis (km/s).\n");
						fprintf(fd,"#Y_vel: Peculiar/physical velocity of the halo along the Y axis (km/s).\n");
						fprintf(fd,"#Z_vel: Peculiar/physical velocity of the halo along the Z axis (km/s).\n");
						fprintf(fd,"#Type: Whether the halo is a main halo (0), subhalo (1), or lost halo with orphan galaxy (2).\n");
						fprintf(fd,"#Halo_ID: Halo ID that corresponds to the ID in the original halo catalogue.\n");
						fprintf(fd,"#Desc_ID: ID of the descendant halo. If there is none then Desc_ID = -1.\n");
						fprintf(fd,"#Up_ID: ID of the most massive halo. If this is a main halo then Up_ID = -1.\n");
					}
				}

				//Determine number of systems that will be printed on each task (depending on redshift and mass)
				nhaloout = 0;
				for (ihalo = 0; ihalo < Nhalos; ihalo++)
					if (H[ihalo].iscale == Output_iscale[i] && H[ihalo].gone == 0 && H[H[ihalo].impeak].mvir >= All.minmass) nhaloout++;

				if (nhaloout > 0)
				{

#ifdef HDF5_SUPPORT
					if (ThisTask == masterTask && All.OutputFormat == 2)
					{
						//Set the path name for the halo data set
						sprintf(path,"/Haloes");
						//Initial size is the number of haloes on the master task
						dims[0] = nhaloout;
						//Allow for variable total size
						maxdims[0] = H5S_UNLIMITED;
						//Create the data space
						hdf5_dataspace_in_file = H5Screate_simple(1, dims, maxdims);
						//Create the data properties
						hdf5_properties = H5Pcreate(H5P_DATASET_CREATE);
						hdf5_status = H5Pset_chunk(hdf5_properties, 1, dims); // set chunk size
						hdf5_status = H5Pset_shuffle(hdf5_properties);        // reshuffle bytes to get better compression ratio
						hdf5_status = H5Pset_deflate(hdf5_properties, 9);     // gzip compression level 9
						hdf5_status = H5Pset_fletcher32(hdf5_properties);     // Fletcher32 checksum on dataset
						//If filters could be set use them to create the data set
						if (H5Pall_filters_avail(hdf5_properties))
							hdf5_dataset = H5Dcreate(halo_file, path, hdf5_ftype, hdf5_dataspace_in_file, H5P_DEFAULT, hdf5_properties, H5P_DEFAULT);
						//Otherwise create the default data set
						else
							hdf5_dataset = H5Dcreate(halo_file, path, hdf5_ftype, hdf5_dataspace_in_file, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
						//Initialise the number of haloes in the set to zero
						sendsum = 0;
						//Define Scalar for Attributes
						hdf5_att_space = H5Screate(H5S_SCALAR);
						//Write Redshift as attribute to dataset
						hdf5_att_redshift = H5Acreate(hdf5_dataset, "Scale Factor", H5T_NATIVE_FLOAT, hdf5_att_space, H5P_DEFAULT, H5P_DEFAULT);
						H5Awrite(hdf5_att_redshift, H5T_IEEE_F32LE, &ScaleFactor[Output_iscale[i]]);
						H5Aclose(hdf5_att_redshift);
						//Write Box Size as attribute to dataset
						tmp               = All.Lbox*All.x_unit;
						hdf5_att_lbox     = H5Acreate(hdf5_dataset, "Box Size", H5T_NATIVE_DOUBLE, hdf5_att_space, H5P_DEFAULT, H5P_DEFAULT);
						H5Awrite(hdf5_att_lbox, H5T_IEEE_F64LE, &tmp);
						H5Aclose(hdf5_att_lbox);
						//Write Hubble parameter as attribute to dataset
						hdf5_att_hubble     = H5Acreate(hdf5_dataset, "Hubble Parameter", H5T_NATIVE_DOUBLE, hdf5_att_space, H5P_DEFAULT, H5P_DEFAULT);
						H5Awrite(hdf5_att_hubble, H5T_IEEE_F64LE, &All.h_100);
						H5Aclose(hdf5_att_hubble);
						//Write Omega_0 as attribute to dataset
						hdf5_att_omega0     = H5Acreate(hdf5_dataset, "Omega_0", H5T_NATIVE_DOUBLE, hdf5_att_space, H5P_DEFAULT, H5P_DEFAULT);
						H5Awrite(hdf5_att_omega0, H5T_IEEE_F64LE, &All.Omega_0);
						H5Aclose(hdf5_att_omega0);
						//Write Omega_Lambda as attribute to dataset
						hdf5_att_omegal     = H5Acreate(hdf5_dataset, "Omega_Lambda", H5T_NATIVE_DOUBLE, hdf5_att_space, H5P_DEFAULT, H5P_DEFAULT);
						H5Awrite(hdf5_att_omegal, H5T_IEEE_F64LE, &All.Omega_Lambda_0);
						H5Aclose(hdf5_att_omegal);
						//Write Omega_Baryon as attribute to dataset
						hdf5_att_omegab     = H5Acreate(hdf5_dataset, "Omega_Baryon", H5T_NATIVE_DOUBLE, hdf5_att_space, H5P_DEFAULT, H5P_DEFAULT);
						H5Awrite(hdf5_att_omegab, H5T_IEEE_F64LE, &All.Omega_Baryon_0);
						H5Aclose(hdf5_att_omegab);
					}
#endif

					//Now go through all tasks that are member of this group
					for (task = masterTask, offset = 0; task <= lastTask; task++)
					{
						//if this task is processed
						if (task == ThisTask)
						{
							//We need to get those haloes
							left_to_send = nhaloout;

							//Tell this to the other tasks
							for (p = masterTask; p <= lastTask; p++)
								if (p != ThisTask)
									MPI_Send(&left_to_send, 1, MPI_INT, p, TAG_NHALOS, MPI_COMM_WORLD);
						}
						//Reveive the number of haloes from active task
						else MPI_Recv(&left_to_send, 1, MPI_INT, task, TAG_NHALOS, MPI_COMM_WORLD, &status);

						//Now send the haloes in increments
						while (left_to_send > 0)
						{
							//Compute how many haloes we can send this turn
							send_this_turn = left_to_send;
							//If they do not all fit into the CommBuffer set to the maximum
							if (send_this_turn > maxlen) send_this_turn = maxlen;

							//Sending task
							if (ThisTask == task)
							{
								//Use the haloout struct for the CommBuffer
								hp = (struct haloout *) CommBuffer;
								//Initialise the number of systems sent to zero
								ihalo = 0;
								//Go through all systems until the maximum that can be sent
								while (ihalo < send_this_turn)
								{ //If we are at the right scale factor and the halo is not gone send it to the master
									if (H[offset].iscale == Output_iscale[i] && H[offset].gone == 0 && H[H[ihalo].impeak].mvir >= All.minmass)
									{
										//Find scale factor when virial mass was half of current value
										amhalf = ScaleFactor[H[H[offset].impeak].iscale];
										iprog = H[H[offset].impeak].iprog;
										if (iprog >=0 )
										{
											while (H[iprog].mvir > 0.5 * H[H[offset].impeak].mvir)
											{
												if (H[iprog].iprog < 0) break;
												iprog = H[iprog].iprog;
											}
											amhalf = ScaleFactor[H[iprog].iscale];
										}

										//Set all halo properties
										hsend.mvir     = log10(H[offset].mvir*All.m_unit);
										hsend.mdot     = H[offset].mdotbary*All.m_unit/All.t_unit/All.f_baryon;
										hsend.mpeak    = log10(H[H[offset].impeak].mvir*All.m_unit);
										hsend.mpeakdot = H[H[offset].impeak].mdotbary*All.m_unit/All.t_unit/All.f_baryon;
										hsend.ampeak   = ScaleFactor[H[H[offset].impeak].iscale];
										hsend.amhalf   = amhalf;
										hsend.rvir     = H[offset].rvir*All.x_unit*1.e3;
										hsend.c        = H[offset].c;
										hsend.l        = H[offset].lambda;
										hsend.x        = H[offset].pos[0];
										hsend.y        = H[offset].pos[1];
										hsend.z        = H[offset].pos[2];
										hsend.u        = H[offset].vel[0];
										hsend.v        = H[offset].vel[1];
										hsend.w        = H[offset].vel[2];
										hsend.t        = H[offset].type;
										hsend.hid      = (long long)(H[offset].haloid)-1;
										hsend.did      = (long long)(H[offset].descid)-1;
										hsend.uid      = (long long)(H[offset].upid)-1;

										//Write this system to the CommBuffer
										*hp++ = hsend;
										//Increment the number of systems that are being sent
										ihalo++;
									}
									//Increment the halo index
									offset++;
								}
							}

							//Receive haloes
							if (ThisTask == masterTask && task != masterTask) MPI_Recv(CommBuffer, (sizeof(struct halo)) * send_this_turn, MPI_BYTE, task, TAG_HDATA, MPI_COMM_WORLD, &status);

							//Send haloes
							if (ThisTask != masterTask && task == ThisTask) MPI_Ssend(CommBuffer, (sizeof(struct halo)) * send_this_turn, MPI_BYTE, masterTask, TAG_HDATA, MPI_COMM_WORLD);

							//Collect haloes and write to file
							if (ThisTask == masterTask)
							{
								//If the output format is hdf5
								if (All.OutputFormat == 2)
								{
#ifdef HDF5_SUPPORT
									//Set the starting index and the count
									start[0] = sendsum;
									count[0] = send_this_turn;
									//Add the number of haloes that are sent this turn to the total number and set the size to the total
									sendsum += send_this_turn;
									dims[0]  = sendsum;
									//If the total size is larger than the initial set size for the data set extend it
									if (sendsum > nhaloout) hdf5_status = H5Dset_extent(hdf5_dataset, dims);
									//Get the data space of the data set
									hdf5_filespace = H5Dget_space(hdf5_dataset);
									//Select the hyperslab that corresponds to the current offset
									hdf5_status = H5Sselect_hyperslab(hdf5_filespace, H5S_SELECT_SET, start, NULL, count, NULL);
									//Set the size to the number of systems sent this turn
									dims[0] = send_this_turn;
									//Create a data space for the memory
									hdf5_dataspace_memory = H5Screate_simple(1, dims, NULL);
									//Write the data set
									hdf5_status = H5Dwrite(hdf5_dataset, hdf5_dtype, hdf5_dataspace_memory, hdf5_filespace, H5P_DEFAULT, CommBuffer);
									//Close the data space of the data set
									H5Sclose(hdf5_filespace);
									//Close the data space for the memory
									H5Sclose(hdf5_dataspace_memory);
#endif
								}//End HDF5 output
								//Otherwise simply write the properties to an ascii file
								else
								{ //Use the haloout struct for the CommBuffer
									hp = (struct haloout *) CommBuffer;
									//Go through all systems that were sent this turn
									for (ihalo = 0; ihalo < send_this_turn; ihalo++)
									{ //Write a new line in the ascii file
										fprintf(fd,
										        "%f %e %f %f %f %f %f %f %f %f %f %f %f %f %f %hu %lld %lld %lld",
										        hp[ihalo].mvir,
										        hp[ihalo].mdot,
										        hp[ihalo].mpeak,
										        hp[ihalo].mpeakdot,
										        hp[ihalo].ampeak,
										        hp[ihalo].amhalf,
										        hp[ihalo].rvir,
										        hp[ihalo].c,
										        hp[ihalo].l,
										        hp[ihalo].x,
										        hp[ihalo].y,
										        hp[ihalo].z,
										        hp[ihalo].u,
										        hp[ihalo].v,
										        hp[ihalo].w,
										        hp[ihalo].t,
										        hp[ihalo].hid,
										        hp[ihalo].did,
										        hp[ihalo].uid);
										fprintf(fd,"\n");
									}//End loop through all haloes
								}//End ascii output
							}//End master task (writer task)
							//Decrease the number of haloes left to send by the number sent this turn
							left_to_send -= send_this_turn;
						}//End halo increments
					}//End loop through all tasks

				}// nhaloout > 0

				//If this is the master task (writer task) close the file (and other objects)
				if (ThisTask == masterTask)
				{
					//For HDF5 output
					if (All.OutputFormat == 2)
					{ //Check if libraries have been included
#ifdef HDF5_SUPPORT
						if (nhaloout > 0)
						{ // Close the dataset
							hdf5_status = H5Dclose(hdf5_dataset);
							// Close the dataset
							hdf5_status = H5Pclose(hdf5_properties);
							// Close the data space
							hdf5_status = H5Sclose(hdf5_dataspace_in_file);
						}
						// Close the file type
						hdf5_status = H5Tclose(hdf5_ftype);
						// Close the data type
						hdf5_status = H5Tclose(hdf5_dtype);
						//Close the file and print to screen if writing was not  successful
						hdf5_status = H5Fclose(halo_file);
						if (hdf5_status != 0) printf("%s Task %5d failed to write to file %s\n",All.startline,ThisTask,buf);
#endif
					}//End HDF5
					//Otherwise close the ascii file
					else
					{ //Close the file
						fclose(fd);
					}//End ascii
				}//End master/writer task

			}//End tasks in the group

		}//End output redshifts

		//Set MPI barrier so that only one task writes to a given file
		MPI_Barrier(MPI_COMM_WORLD);
	}//End loop through all groups

	//Free the CommBuffer
	efree(CommBuffer);

}
#endif


#ifdef WRITE_MAINBRANCH
/*! \brief This function prints the main branch of galaxies with a specified halo or stellar mass at a specified redshift
 *
 *  This function checks each tree is the main progenitor has a specific stellar or halo mass at a specific
 *  redshift. If the tree fulfils this criterion it is sent to the master task and then written to a file.
 *  Output can be generated in either ascii format, where main branch data is written in blocks, or in HDF5
 *  format, where main branch data sets are stored individually in groups according to their mass bin.
 *  There are a total of All.NumOutputFiles files and All.NumFilesInParallel files are written in parallel.
 *  Each master task collects the halo data from its slave tasks and writes it. If there is one output file
 *  it will be called mainbranches.out otherwise there will be n files called mainbranches.i where i < n.
 */
void output_mainbranch(void)
{

	int i, j, k, filenr, masterTask, lastTask, ngroups, task;
	int iroot, ihalo, iprog, present_at_z, Nmainbranch, Nrootloop, mbflag, imb, Nroot;
	int *roothalos;
	double mass;
	char buf[NSTRING];
	FILE *fd;
	MPI_Status status;

	//Structure that stores all halo/galaxy properties that will be saved
	struct mbout {
		float a,mh,mdot,mhp,amp,r,c,l,ms,sfr,icm,x,y,z,u,v,w;
		unsigned short t;
		long long hid, did, uid;
	} *mainbranch, mb;

#ifdef HDF5_SUPPORT
	hid_t mb_file, hdf5_dtype, hdf5_ftype, set_group;
	hid_t hdf5_properties, hdf5_dataset, hdf5_dataspace_in_file, atype;
	hid_t hdf5_att_space, hdf5_att_redshift, hdf5_att_massbin, hdf5_att_binsize, hdf5_att_masstag;
	herr_t hdf5_status;
	hsize_t dims[1], maxdims[1];
	int asize;
	char path[NSTRING], masstag[NSTRING];
#endif

	//Find number of root haloes (at z=0) for each task
	Nroot = 0;
	for (ihalo = 0; ihalo < Nhalos; ihalo++)
		if (H[ihalo].iscale == All.NTimesteps - 1)
			Nroot++;

	//Allocate memory for the array containing indices of root haloes
	roothalos = (int *) emalloc_movable(&roothalos,  "RootHalos", Nroot *sizeof(int));

	//Fill array containing indices of root haloes
	iroot = 0;
	for (ihalo = 0; ihalo < Nhalos; ihalo++)
	{      //Each halo at z=0 is a root halo (also orphans)
		if (H[ihalo].iscale == All.NTimesteps - 1)
		{            //Fill array
			roothalos[iroot] = ihalo;
			iroot++;
		}
	}

	//Print what is done...
	if (ThisTask == 0) printf("%s\n%s Writing main branches in 1st universe to %d files...\n", All.fullline,All.startline,All.NumOutputFiles);

	//Check if there are enough tasks to write all output files
	if (All.NTaskPerUniverse < All.NumOutputFiles && ThisTask == 0)
		endrun("Fatal error. Number of processors must be larger or equal than All.NumOutputFiles.");

	//Assign processors to output files
	distribute_file(All.NTaskPerUniverse, All.NumOutputFiles, 0, 0, All.NTaskPerUniverse - 1, &filenr, &masterTask, &lastTask);

	//If the output format has been set to HDF5
	if (All.OutputFormat == 2)
	{
		//Check if HDF5 libraries are set
#ifdef HDF5_SUPPORT
		if (All.NumOutputFiles > 1)
			sprintf(buf, "%s/mainbranches/mainbranches.S%d.%d.h5", All.OutputDir, All.MainBranch_iscale, filenr);
		else
			sprintf(buf, "%s/mainbranches/mainbranches.S%d.h5", All.OutputDir, All.MainBranch_iscale);
		//If not return
#else
		if (ThisTask == 0) printf("%s Output format has been set to 2 but HDF5 support was not enabled.\n",All.startline);
		return;
#endif
	}
	//Otherwise use standard ascii files
	else
	{
		if (All.NumOutputFiles > 1)
			sprintf(buf, "%s/mainbranches/mainbranches.S%d.%d", All.OutputDir, All.MainBranch_iscale, filenr);
		else
			sprintf(buf, "%s/mainbranches/mainbranches.S%d.out", All.OutputDir, All.MainBranch_iscale);
	}

	//Get number of groups
	ngroups = All.NumOutputFiles / All.NumFilesInParallel;
	if ((All.NumOutputFiles % All.NumFilesInParallel)) ngroups++;

	//For each group do...
	for (j = 0; j < ngroups; j++)
	{
		//This task will be processed now
		if ((filenr / All.NumFilesInParallel) == j && MasterTask == 0)
		{ //The masterTask opens the file
			if (ThisTask == masterTask)
			{ //If the hdf5 output format has been selected
				if (All.OutputFormat == 2)
				{ //Check if the libraries have been included
#ifdef HDF5_SUPPORT
					//Create hdf5 file
					mb_file = H5Fcreate(buf, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
					//Print to screen
					printf("%s Writing main branch for galaxies at z = %3.2f: '%s' (file %d of %d)\n",All.startline, All.mainBranchRedshift, buf, filenr+1, All.NumOutputFiles);
					//Specify data type
					hdf5_dtype  = H5Tcreate(H5T_COMPOUND, sizeof(struct mbout));
					k = 0;
					hdf5_status = H5Tinsert(hdf5_dtype, "Scale_factor",       HOFFSET(struct mbout, a),    H5T_NATIVE_FLOAT); k+=4;
					hdf5_status = H5Tinsert(hdf5_dtype, "Halo_mass",          HOFFSET(struct mbout, mh),   H5T_NATIVE_FLOAT); k+=4;
					hdf5_status = H5Tinsert(hdf5_dtype, "Halo_growth_rate",   HOFFSET(struct mbout, mdot), H5T_NATIVE_FLOAT); k+=4;
					hdf5_status = H5Tinsert(hdf5_dtype, "Halo_mass_peak",     HOFFSET(struct mbout, mhp),  H5T_NATIVE_FLOAT); k+=4;
					hdf5_status = H5Tinsert(hdf5_dtype, "Scale_peak_mass",    HOFFSET(struct mbout, amp),  H5T_NATIVE_FLOAT); k+=4;
					hdf5_status = H5Tinsert(hdf5_dtype, "Halo_radius",        HOFFSET(struct mbout, r),    H5T_NATIVE_FLOAT); k+=4;
					hdf5_status = H5Tinsert(hdf5_dtype, "Concentration",      HOFFSET(struct mbout, c),    H5T_NATIVE_FLOAT); k+=4;
					hdf5_status = H5Tinsert(hdf5_dtype, "Halo_spin",          HOFFSET(struct mbout, l),    H5T_NATIVE_FLOAT); k+=4;
					hdf5_status = H5Tinsert(hdf5_dtype, "Stellar_mass",       HOFFSET(struct mbout, ms),   H5T_NATIVE_FLOAT); k+=4;
					hdf5_status = H5Tinsert(hdf5_dtype, "SFR",                HOFFSET(struct mbout, sfr),  H5T_NATIVE_FLOAT); k+=4;
					hdf5_status = H5Tinsert(hdf5_dtype, "Intra_cluster_mass", HOFFSET(struct mbout, icm),  H5T_NATIVE_FLOAT); k+=4;
					hdf5_status = H5Tinsert(hdf5_dtype, "X_pos",              HOFFSET(struct mbout, x),    H5T_NATIVE_FLOAT); k+=4;
					hdf5_status = H5Tinsert(hdf5_dtype, "Y_pos",              HOFFSET(struct mbout, y),    H5T_NATIVE_FLOAT); k+=4;
					hdf5_status = H5Tinsert(hdf5_dtype, "Z_pos",              HOFFSET(struct mbout, z),    H5T_NATIVE_FLOAT); k+=4;
					hdf5_status = H5Tinsert(hdf5_dtype, "X_vel",              HOFFSET(struct mbout, u),    H5T_NATIVE_FLOAT); k+=4;
					hdf5_status = H5Tinsert(hdf5_dtype, "Y_vel",              HOFFSET(struct mbout, v),    H5T_NATIVE_FLOAT); k+=4;
					hdf5_status = H5Tinsert(hdf5_dtype, "Z_vel",              HOFFSET(struct mbout, w),    H5T_NATIVE_FLOAT); k+=4;
					hdf5_status = H5Tinsert(hdf5_dtype, "Type",               HOFFSET(struct mbout, t),    H5T_NATIVE_USHORT); k+=2;
					hdf5_status = H5Tinsert(hdf5_dtype, "Halo_ID",            HOFFSET(struct mbout, hid),  H5T_NATIVE_LLONG); k+=8;
					hdf5_status = H5Tinsert(hdf5_dtype, "Desc_ID",            HOFFSET(struct mbout, did),  H5T_NATIVE_LLONG); k+=8;
					hdf5_status = H5Tinsert(hdf5_dtype, "Up_ID",              HOFFSET(struct mbout, uid),  H5T_NATIVE_LLONG); k+=8;
					//Specify file type
					hdf5_ftype  = H5Tcreate(H5T_COMPOUND, k);
					k = 0;
					hdf5_status = H5Tinsert(hdf5_ftype, "Scale_factor",       k,  H5T_IEEE_F32LE); k+=4;
					hdf5_status = H5Tinsert(hdf5_ftype, "Halo_mass",          k,  H5T_IEEE_F32LE); k+=4;
					hdf5_status = H5Tinsert(hdf5_ftype, "Halo_growth_rate",   k,  H5T_IEEE_F32LE); k+=4;
					hdf5_status = H5Tinsert(hdf5_ftype, "Halo_mass_peak",     k,  H5T_IEEE_F32LE); k+=4;
					hdf5_status = H5Tinsert(hdf5_ftype, "Scale_peak_mass",    k,  H5T_IEEE_F32LE); k+=4;
					hdf5_status = H5Tinsert(hdf5_ftype, "Halo_radius",        k,  H5T_IEEE_F32LE); k+=4;
					hdf5_status = H5Tinsert(hdf5_ftype, "Concentration",      k,  H5T_IEEE_F32LE); k+=4;
					hdf5_status = H5Tinsert(hdf5_ftype, "Halo_spin",          k,  H5T_IEEE_F32LE); k+=4;
					hdf5_status = H5Tinsert(hdf5_ftype, "Stellar_mass",       k,  H5T_IEEE_F32LE); k+=4;
					hdf5_status = H5Tinsert(hdf5_ftype, "SFR",                k,  H5T_IEEE_F32LE); k+=4;
					hdf5_status = H5Tinsert(hdf5_ftype, "Intra_cluster_mass", k,  H5T_IEEE_F32LE); k+=4;
					hdf5_status = H5Tinsert(hdf5_ftype, "X_pos",              k,  H5T_IEEE_F32LE); k+=4;
					hdf5_status = H5Tinsert(hdf5_ftype, "Y_pos",              k,  H5T_IEEE_F32LE); k+=4;
					hdf5_status = H5Tinsert(hdf5_ftype, "Z_pos",              k,  H5T_IEEE_F32LE); k+=4;
					hdf5_status = H5Tinsert(hdf5_ftype, "X_vel",              k,  H5T_IEEE_F32LE); k+=4;
					hdf5_status = H5Tinsert(hdf5_ftype, "Y_vel",              k,  H5T_IEEE_F32LE); k+=4;
					hdf5_status = H5Tinsert(hdf5_ftype, "Z_vel",              k,  H5T_IEEE_F32LE); k+=4;
					hdf5_status = H5Tinsert(hdf5_ftype, "Type",               k,  H5T_STD_U16LE); k+=2;
					hdf5_status = H5Tinsert(hdf5_ftype, "Halo_ID",            k,  H5T_STD_I64LE); k+=8;
					hdf5_status = H5Tinsert(hdf5_ftype, "Desc_ID",            k,  H5T_STD_I64LE); k+=8;
					hdf5_status = H5Tinsert(hdf5_ftype, "Up_ID",              k,  H5T_STD_I64LE); k+=8;
					//Add attribute space to file
					hdf5_att_space = H5Screate(H5S_SCALAR);
					//Write selection redshift as attribute to file
					hdf5_att_redshift = H5Acreate(mb_file, "Selection_Redshift", H5T_NATIVE_DOUBLE, hdf5_att_space, H5P_DEFAULT, H5P_DEFAULT);
					H5Awrite(hdf5_att_redshift, H5T_IEEE_F64LE, &All.mainBranchRedshift);
					H5Aclose(hdf5_att_redshift);
					//Close attribute space
					H5Sclose(hdf5_att_space);
#endif
				}
				//Otherwise open standard ascii file
				else
				{ //Open file
					if (!(fd = fopen(buf, "w")))
					{ //If not possible print to screen and abort
						printf("%s Can't open file `%s' for writing output.\n",All.startline,buf);
						endrun("file open error");
					}
					//If file can be opened print to screen
					printf("%s Writing main branch for galaxies at z = %3.2f: '%s' (file %d of %d)\n",All.startline, All.mainBranchRedshift, buf, filenr+1, All.NumOutputFiles);
					//Print column names
					task = 0;
					fprintf(fd,"#Scale_factor(%d)",task); task++;
					fprintf(fd," Halo_mass(%d)",task); task++;
					fprintf(fd," Halo_growth_rate(%d)",task); task++;
					fprintf(fd," Halo_mass_peak(%d)",task); task++;
					fprintf(fd," Scale_peak_mass(%d)",task); task++;
					fprintf(fd," Halo_radius(%d)",task); task++;
					fprintf(fd," Concentration(%d)",task); task++;
					fprintf(fd," Halo_spin(%d)",task); task++;
					fprintf(fd," Stellar_mass(%d)",task); task++;
					fprintf(fd," SFR(%d)",task); task++;
					fprintf(fd," Intra_cluster_mass(%d)",task); task++;
					fprintf(fd," X_pos(%d)",task); task++;
					fprintf(fd," Y_pos(%d)",task); task++;
					fprintf(fd," Z_pos(%d)",task); task++;
					fprintf(fd," X_vel(%d)",task); task++;
					fprintf(fd," Y_vel(%d)",task); task++;
					fprintf(fd," Z_vel(%d)",task); task++;
					fprintf(fd," Type(%d)",task); task++;
					fprintf(fd," Halo_ID(%d)",task); task++;
					fprintf(fd," Desc_ID(%d)",task); task++;
					fprintf(fd," Up_ID(%d)",task); task++;
					fprintf(fd,"\n");
					fprintf(fd,"#Emerge Branch: %s\n",BRANCH);
					fprintf(fd,"#Emerge Version: %s\n",VERSION);
					fprintf(fd,"#Model Name: %s\n",All.model_name);
					fprintf(fd,"#Box Size: %f (Mpc)\n",All.Lbox*All.x_unit);
					fprintf(fd,"#Hubble Parameter: %f\n",All.h_100);
					fprintf(fd,"#Omega_0: %f\n",All.Omega_0);
					fprintf(fd,"#Omega_Lambda: %f\n",All.Omega_Lambda_0);
					fprintf(fd,"#Omega_Baryon: %f\n",All.Omega_Baryon_0);
					//Write a header including the selection redshift to the first line
					fprintf(fd,"#This file contains all main branches with a mass of %s at z = %4.2f\n",All.output_mass_mb,All.mainBranchRedshift);
					//Add which mass has been used for selecting the trees
					if (All.MainBranchMassType == 1) fprintf(fd,"#Stellar mass has been used for this selection.\n");
					else fprintf(fd,"#Halo mass has been used for this selection.\n");
					fprintf(fd,"#Halo_mass: Current virial mass of the halo (log Msun).\n");
					fprintf(fd,"#Halo_growth_rate: Growth rate of the halo (Msun/yr)\n");
					fprintf(fd,"#Halo_mass_peak: Peak virial mass of the halo through its history (log Msun).\n");
					fprintf(fd,"#Scale_peak_mass: Scale Factor when halo had its peak mass.\n");
					fprintf(fd,"#Halo_radius: Current virial radius of the halo (kpc).\n");
					fprintf(fd,"#Concentration: Virial radius divided by halo scale length.\n");
					fprintf(fd,"#Halo_spin: Halo spin parameter.\n");
					fprintf(fd,"#Stellar_mass: Intrinsic stellar mass of the galaxy (log Msun).\n");
					fprintf(fd,"#SFR: Star formation rate of the galaxy (Msun/yr).\n");
					fprintf(fd,"#Intra_cluster_mass: Stellar mass in the halo (log Msun).\n");
					fprintf(fd,"#X_pos: Comoving galaxy position along the X axis (Mpc).\n");
					fprintf(fd,"#Y_pos: Comoving galaxy position along the Y axis (Mpc).\n");
					fprintf(fd,"#Z_pos: Comoving galaxy position along the Z axis (Mpc).\n");
					fprintf(fd,"#X_vel: Peculiar/physical velocity of the galaxy along the X axis (km/s).\n");
					fprintf(fd,"#Y_vel: Peculiar/physical velocity of the galaxy along the Y axis (km/s).\n");
					fprintf(fd,"#Z_vel: Peculiar/physical velocity of the galaxy along the Z axis (km/s).\n");
					fprintf(fd,"#Type: Whether the galaxy is a central (0), satellite (1), or orphan (2).\n");
					fprintf(fd,"#Halo_ID: Halo ID that corresponds to the ID in the original halo catalogue.\n");
					fprintf(fd,"#Desc_ID: ID of the descendant halo. If there is none then Desc_ID = -1.\n");
					fprintf(fd,"#Up_ID: ID of the most massive halo. If this is a main halo then Up_ID = -1.\n");
					//Print an empty line
					fprintf(fd,"##################################################\n\n");
				}//End ascii
			}//End master task

			// Loop over all selected mass bins for the main branch output
			for (i = 0; i < All.Noutputbranch; i++)
			{ //If this is the master task start with a new group
				if (ThisTask == masterTask)
				{ //If we write an HDF5 file, create a group and write the bin's mass and size, and a tag which mass has been used as an attribute
					if (All.OutputFormat == 2)
					{
#ifdef HDF5_SUPPORT
						//Set the path name for the mass bin
						sprintf(path,"/MainBranch_M%03d",i);
						//Create new group for the SMFs
						set_group = H5Gcreate(mb_file, path, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
						//Add attribute space to group
						hdf5_att_space = H5Screate(H5S_SCALAR);
						//Write mass bin as attribute to group
						hdf5_att_massbin = H5Acreate(set_group, "Mass_bin", H5T_NATIVE_DOUBLE, hdf5_att_space, H5P_DEFAULT, H5P_DEFAULT);
						mass = (double)(MainBranchMasses[i]);
						H5Awrite(hdf5_att_massbin, H5T_IEEE_F64LE, &mass);
						H5Aclose(hdf5_att_massbin);
						//Write mass bin size as attribute to group
						hdf5_att_binsize = H5Acreate(set_group, "Mass_bin_size", H5T_NATIVE_DOUBLE, hdf5_att_space, H5P_DEFAULT, H5P_DEFAULT);
						mass = (double)(MainBranchBinSize[i]);
						H5Awrite(hdf5_att_binsize, H5T_IEEE_F64LE, &mass);
						H5Aclose(hdf5_att_binsize);
						//Write mass tag as attribute to group
						atype = H5Tcopy(H5T_C_S1);
						//Select which mass has been used
						if (All.MainBranchMassType == 1)
						{ //For stellar mass we need 12 characters
							asize = 12;
							sprintf(masstag, "Stellar_Mass");
						}
						else
						{ //For halo mass we need 9 characters
							asize = 9;
							sprintf(masstag, "Halo_Mass");
						}
						//Set size
						H5Tset_size(atype, asize);
						H5Tset_strpad(atype,H5T_STR_NULLTERM);
						hdf5_att_masstag = H5Acreate(set_group, "Mass_Type", atype, hdf5_att_space, H5P_DEFAULT, H5P_DEFAULT);
						H5Awrite(hdf5_att_masstag, atype, masstag);
						H5Aclose(hdf5_att_masstag);
						//Close String Type
						H5Tclose(atype);
						//Close attribute space
						H5Sclose(hdf5_att_space);
#endif
					}
					else //If we use an ascii file
					{ //Print the mass and width of this bin
						fprintf(fd,"# The following trees are for masses of %7.5f with a bin size of %7.5f\n",MainBranchMasses[i],MainBranchBinSize[i]);
						//Print line
						fprintf(fd,"##################################################\n");
					}
				} //Done for the master task

				//Allocate memory for one main branch over the full possible number of snapshots
				mainbranch = emalloc_movable(&mainbranch,"Mainbranch", All.NTimesteps * sizeof(struct mbout));
				//Set main branch index to 0
				imb = 0;

				//Now go through all tasks that are member of this group
				for (task = masterTask; task <= lastTask; task++)
				{ //If this task is processed
					if (task == ThisTask)
					{ //If we're beyond the master task, send number of root haloes to master task
						if (ThisTask != masterTask) MPI_Ssend(&Nroot, 1, MPI_INT, masterTask, TAG_NTREES, MPI_COMM_WORLD);
						//Go through all root haloes on this task
						for (iroot = 0; iroot < Nroot; iroot++)
						{ //Initialise main branch flag to 0
							mbflag = 0;
							//Set flag for presence at selected redshift to 1
							present_at_z = 1;
							//Set iprog to root halo
							iprog = roothalos[iroot];
							//Loop until main branch progenitor is at selected redshift
							while (H[iprog].iscale > All.MainBranch_iscale)
							{ //If progenitor does not exist
								if (H[iprog].iprog < 0)
								{ //Set presence flag to 0 and break
									present_at_z = 0;
									break;
								}
								//Otherwise go to main progenitor
								else iprog = H[iprog].iprog;
							}//Done with looping over progenitors
							//Select which mass is used for selection
							if (All.MainBranchMassType == 1) mass = log10(H[iprog].mstar*All.m_unit);
							else mass = log10(H[iprog].mvir*All.m_unit);
							//If mass is in the bin and halo is still present
							if (mass >= MainBranchMasses[i] - 0.5 * MainBranchBinSize[i] &&
							    mass <  MainBranchMasses[i] + 0.5 * MainBranchBinSize[i] &&
							    present_at_z == 1)
							{ //System is relevant - so deal with it and set flag to 1
								mbflag = 1;
								//Tell master task that this is a relevant tree
								if (ThisTask != masterTask) MPI_Ssend(&mbflag, 1, MPI_INT, masterTask, TAG_COUNT, MPI_COMM_WORLD);
								//Start again at z = 0
								iprog = roothalos[iroot];
								//Set counter for the number of entries in this tree to 0
								Nmainbranch = 0;
								//Loop over progenitors until there is no more progenitor
								while (H[iprog].iprog >= 0)
								{ //Set all properties for the mb structure
									mb.a    = H[iprog].a;
									mb.mh   = log10(H[iprog].mvir*All.m_unit);
									mb.mdot = H[iprog].mdotbary*All.m_unit/All.t_unit/All.f_baryon;
									mb.mhp  = log10(H[H[iprog].impeak].mvir*All.m_unit);
									mb.amp  = ScaleFactor[H[H[iprog].impeak].iscale];
									mb.r    = H[iprog].rvir*All.x_unit*1.e3;
									mb.c    = H[iprog].c;
									mb.l    = H[iprog].lambda;
									mb.ms   = log10(H[iprog].mstar*All.m_unit);
									mb.sfr  = H[iprog].sfr*All.m_unit/All.t_unit;
									mb.icm  = log10(H[iprog].icm*All.m_unit);
									mb.x    = H[iprog].pos[0];
									mb.y    = H[iprog].pos[1];
									mb.z    = H[iprog].pos[2];
									mb.u    = H[iprog].vel[0];
									mb.v    = H[iprog].vel[1];
									mb.w    = H[iprog].vel[2];
									mb.t    = H[iprog].type;
									mb.hid  = (long long)(H[iprog].haloid)-1;
									mb.did  = (long long)(H[iprog].descid)-1;
									mb.uid  = (long long)(H[iprog].upid)-1;
									//Write it all to next main branch entry
									mainbranch[Nmainbranch] = mb;
									//Go to next main progenitor
									iprog = H[iprog].iprog;
									//And increment the counter
									Nmainbranch++;
								} //End loop through main branch
								//If we're beyond the master task, send number of haloes in main branch and then send full main branch
								if (ThisTask != masterTask)
								{
									MPI_Ssend(&Nmainbranch, 1, MPI_INT, masterTask, TAG_NHALOS, MPI_COMM_WORLD);
									MPI_Ssend(mainbranch, (sizeof(struct mbout)) * Nmainbranch, MPI_BYTE, masterTask, TAG_HDATA, MPI_COMM_WORLD);
								}
								else //If this is the master task, write the main branch in a file
								{ //Check if we write an HDF5 file
									if (All.OutputFormat == 2)
									{
#ifdef HDF5_SUPPORT
										//Set the path name for the halo data set
										sprintf(path,"/MainBranch_M%03d/Tree_%07d",i,imb);
										//Initial size is the number of haloes on the master task
										dims[0] = Nmainbranch;
										//Allow for variable total size
										maxdims[0] = H5S_UNLIMITED;
										//Create the data space
										hdf5_dataspace_in_file = H5Screate_simple(1, dims, maxdims);
										//Create the data properties
										hdf5_properties = H5Pcreate(H5P_DATASET_CREATE);
										hdf5_status = H5Pset_chunk(hdf5_properties, 1, dims); // set chunk size
										hdf5_status = H5Pset_shuffle(hdf5_properties);        // reshuffle bytes to get better compression ratio
										hdf5_status = H5Pset_deflate(hdf5_properties, 9);     // gzip compression level 9
										hdf5_status = H5Pset_fletcher32(hdf5_properties);     // Fletcher32 checksum on dataset
										//If filters could be set use them to create the data set
										if (H5Pall_filters_avail(hdf5_properties))
											hdf5_dataset = H5Dcreate(mb_file, path, hdf5_ftype, hdf5_dataspace_in_file, H5P_DEFAULT, hdf5_properties, H5P_DEFAULT);
										//Otherwise create the default data set
										else
											hdf5_dataset = H5Dcreate(mb_file, path, hdf5_ftype, hdf5_dataspace_in_file, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
										//Write main branch to data set
										hdf5_status = H5Dwrite(hdf5_dataset, hdf5_dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, mainbranch);
										//Close the dataset
										hdf5_status = H5Dclose(hdf5_dataset);
										//Close the properties
										hdf5_status = H5Pclose(hdf5_properties);
										//Close the data space
										hdf5_status = H5Sclose(hdf5_dataspace_in_file);
#endif
									}
									else //Write main branch to the ascii file
									{ //First we write the index number of this main branch
										fprintf(fd,"#%03d_%07d\n",i,imb);
										//Then write all systems in the main branch
										for (ihalo = 0; ihalo < Nmainbranch; ihalo++)
										{
											fprintf(fd,
											        "%f %f %e %f %f %f %f %f %f %e %f %f %f %f %f %f %f %hu %lld %lld %lld",
											        mainbranch[ihalo].a,
											        mainbranch[ihalo].mh,
											        mainbranch[ihalo].mdot,
											        mainbranch[ihalo].mhp,
											        mainbranch[ihalo].amp,
											        mainbranch[ihalo].r,
											        mainbranch[ihalo].c,
											        mainbranch[ihalo].l,
											        mainbranch[ihalo].ms,
											        mainbranch[ihalo].sfr,
											        mainbranch[ihalo].icm,
											        mainbranch[ihalo].x,
											        mainbranch[ihalo].y,
											        mainbranch[ihalo].z,
											        mainbranch[ihalo].u,
											        mainbranch[ihalo].v,
											        mainbranch[ihalo].w,
											        mainbranch[ihalo].t,
											        mainbranch[ihalo].hid,
											        mainbranch[ihalo].did,
											        mainbranch[ihalo].uid);
											fprintf(fd,"\n");
										}
										//And separate with an empty line
										fprintf(fd,"\n");
									} //Done writing the main branch for master task
								} //Done with the master task
								//Increment index for the main branch in this mass bin
								imb++;
							} //Done with this selected system
							else
							//Tree does not have the right mass at the selected redshift
							{ //System is not relevant - set flag to 0
								mbflag = 0;
								//Tell master task that this tree is not relevant
								if (ThisTask != masterTask) MPI_Ssend(&mbflag, 1, MPI_INT, masterTask, TAG_COUNT, MPI_COMM_WORLD);
							}
						}//Done with the loop over all trees
					}//Done with active task task
					//Now deal with master task while other task is active
					else if (ThisTask == masterTask && task != masterTask)
					{ //Receive number of root haloes on the active task
						MPI_Recv(&Nrootloop, 1, MPI_INT, task, TAG_NTREES, MPI_COMM_WORLD, &status);
						//Go through all root haloes in the active task
						for (iroot = 0; iroot < Nrootloop; iroot++)
						{ //Receive the flag whether the halo is selected or not
							MPI_Recv(&mbflag, 1, MPI_INT, task, TAG_COUNT, MPI_COMM_WORLD, &status);
							//If this tree is relevant we need to print it
							if (mbflag == 1)
							{ //Receive number of haloes in main branch and then full main branch
								MPI_Recv(&Nmainbranch, 1, MPI_INT, task, TAG_NHALOS, MPI_COMM_WORLD, &status);
								MPI_Recv(mainbranch, (sizeof(struct mbout)) * Nmainbranch, MPI_BYTE, task, TAG_HDATA, MPI_COMM_WORLD, &status);
								//Now we need to print it, either to HDF5 or ascii - check if it's HDF5 first
								if (All.OutputFormat == 2)
								{ //If format is set to HDF5 then write it to the current data group
#ifdef HDF5_SUPPORT
									//Set the path name for the halo data set
									sprintf(path,"/MainBranch_M%03d/Tree_%07d",i,imb);
									//Initial size is the number of haloes on the master task
									dims[0] = Nmainbranch;
									//Allow for variable total size
									maxdims[0] = H5S_UNLIMITED;
									//Create the data space
									hdf5_dataspace_in_file = H5Screate_simple(1, dims, maxdims);
									//Create the data properties
									hdf5_properties = H5Pcreate(H5P_DATASET_CREATE);
									hdf5_status = H5Pset_chunk(hdf5_properties, 1, dims); // set chunk size
									hdf5_status = H5Pset_shuffle(hdf5_properties);        // reshuffle bytes to get better compression ratio
									hdf5_status = H5Pset_deflate(hdf5_properties, 9);     // gzip compression level 9
									hdf5_status = H5Pset_fletcher32(hdf5_properties);     // Fletcher32 checksum on dataset
									//If filters could be set use them to create the data set
									if (H5Pall_filters_avail(hdf5_properties))
										hdf5_dataset = H5Dcreate(mb_file, path, hdf5_ftype, hdf5_dataspace_in_file, H5P_DEFAULT, hdf5_properties, H5P_DEFAULT);
									//Otherwise create the default data set
									else
										hdf5_dataset = H5Dcreate(mb_file, path, hdf5_ftype, hdf5_dataspace_in_file, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
									//Write main branch to data set
									hdf5_status = H5Dwrite(hdf5_dataset, hdf5_dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, mainbranch);
									//Close the dataset
									hdf5_status = H5Dclose(hdf5_dataset);
									//Close the dataset
									hdf5_status = H5Pclose(hdf5_properties);
									//Close the data space
									hdf5_status = H5Sclose(hdf5_dataspace_in_file);
#endif
								}
								else //Write main branch to the ascii file
								{ //First we write the index number of this main branch
									fprintf(fd,"#%03d_%07d\n",i,imb);
									//Then write all systems in the main branch
									for (ihalo = 0; ihalo < Nmainbranch; ihalo++)
									{
										fprintf(fd,
										        "%f %f %e %f %f %f %f %f %f %e %f %f %f %f %f %f %f %hu %lld %lld %lld",
										        mainbranch[ihalo].a,
										        mainbranch[ihalo].mh,
										        mainbranch[ihalo].mdot,
										        mainbranch[ihalo].mhp,
										        mainbranch[ihalo].amp,
										        mainbranch[ihalo].r,
										        mainbranch[ihalo].c,
										        mainbranch[ihalo].l,
										        mainbranch[ihalo].ms,
										        mainbranch[ihalo].sfr,
										        mainbranch[ihalo].icm,
										        mainbranch[ihalo].x,
										        mainbranch[ihalo].y,
										        mainbranch[ihalo].z,
										        mainbranch[ihalo].u,
										        mainbranch[ihalo].v,
										        mainbranch[ihalo].w,
										        mainbranch[ihalo].t,
										        mainbranch[ihalo].hid,
										        mainbranch[ihalo].did,
										        mainbranch[ihalo].uid);
										fprintf(fd,"\n");
									}
									//And separate with an empty line
									fprintf(fd,"\n");
								} //Done writing the main branch for master task
								//Increment index for the main branch in this mass bin
								imb++;
							}//Done with this relevant tree - now go to next tree
						}//Loop over all trees on active task
					}//Done with the master task while another task is active
				}//End loop through all tasks

				//Free the memomry for the main branch structure
				efree_movable(mainbranch);

				//If this is the master task (writer task) close the group
#ifdef HDF5_SUPPORT
				if (ThisTask == masterTask && All.OutputFormat == 2)
				{ // Close set goup
					hdf5_status = H5Gclose(set_group);
				}
#endif

			}//Done with this mass bins - go to next

			//If this is the master task (writer task) close the file (and other objects)
			if (ThisTask == masterTask)
			{ //For HDF5 output
				if (All.OutputFormat == 2)
				{ //Check if libraries have been included
#ifdef HDF5_SUPPORT
					// Close the file type
					hdf5_status = H5Tclose(hdf5_ftype);
					// Close the data type
					hdf5_status = H5Tclose(hdf5_dtype);
					//Close the file and print to screen if writing was not  successful
					hdf5_status = H5Fclose(mb_file);
					if (hdf5_status != 0) printf("%s Task %5d failed to write to file %s\n",All.startline,ThisTask,buf);
#endif
				}//End HDF5
				//Otherwise close the ascii file
				else
				{ //Close the file
					fclose(fd);
				}//End ascii
			}//End master/writer task

		}//End task in the group

	}//End groups

	//Free array with root haloes
	efree_movable(roothalos);

}
#endif

