///////////////////////////////////////////////////////////////////////////////////////////////////
// Emerge code - File read_data.c                                                                //
///////////////////////////////////////////////////////////////////////////////////////////////////
///
/// \file read_data.c
/// \brief Contains functions that read statistical galaxy data
///
/// This file contains all functions that read the statistical galaxy data from various files.
/// This includes the stellar mass functions, the quenched fractions, the cosmic star formation
/// rate density, and the specific star formation rates.
///
///////////////////////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include "allvars.h"
#include "proto.h"


/*! \brief This function calls the read functions for all required data
 *
 *  This function is called at the start of the programme and will call all read functions
 *  for which the READ_XYZ options have been selected. All data is read into galaxy_data structures
 */
void read_data(void)
{
	//Print what is done...
	if (ThisTask == 0) printf("%s\n%s Reading observed data...\n", All.fullline,All.startline);

#ifdef READ_SMF
	read_smf();
	All.Nobs++;
#endif

#ifdef READ_FQ
	read_fq();
	All.Nobs++;
#endif

#ifdef READ_CSFRD
	read_csfrd();
	All.Nobs++;
#endif

#ifdef READ_SSFR
	read_ssfr();
	All.Nobs++;
#endif

#ifdef READ_WP
	read_wp();
	All.Nobs++;
#endif

	MPI_Barrier(MPI_COMM_WORLD);

}


#ifdef READ_SMF
/*! \brief This function reads all observed stellar mass functions from a file
 *
 *  First the galaxy_data structure #smf is allocated. The function then reads all lines and if it is
 *  a header line (#) then the redshift and the corrections for IMF and Hubble_h are read. Then all
 *  data points in this set are read and corrected. Redshifts are stored in smf.bin.
 */
void read_smf(void)
{
	FILE *ifp;
	const int linesize = 500;
	char line[linesize], buf[linesize], *buffer;
	int i,ismf,nset,iset;
	float imf,hubble,sigmaplus,sigmaminus;

	//Task 0 reads all data
	if (ThisTask == 0)
	{
		//Print what is done...
		printf("%s Reading observed stellar mass functions.\n",All.startline);

		//Open file
		if(!(ifp=fopen(All.smffile_name,"r")))
		{
			sprintf(buf, "Cannot open file `%s' for reading the stellar mass functions.", All.smffile_name);
			endrun(buf);
		}
		//Skip the header
		fgets(line,linesize,ifp);
		//Initialise counter
		All.Nsmfset = 0;
		All.Nsmf    = 0;
		while (fgets(line,linesize,ifp))
		{
			if (line[0] == '#')
			{
				sscanf(line,"# %d\n",&nset);
				All.Nsmfset++;
				All.Nsmf += nset;
			}
		}
		//Close file
		fclose(ifp);
	}

	MPI_Bcast(&All.Nsmf, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&All.Nsmfset, 1, MPI_INT, 0, MPI_COMM_WORLD);
	SmfSet = emalloc("SMFSET", All.Nsmfset * sizeof(struct data_set));
	Smf    = emalloc("SMF",    All.Nsmf    * sizeof(struct galaxy_data));

	//Task 0 reads all data
	if (ThisTask == 0)
	{
		buffer = emalloc("BUFFER", SSTRING * sizeof(char));

		//Open file
		if(!(ifp=fopen(All.smffile_name,"r")))
		{
			sprintf(buf, "Cannot open file `%s' for reading the stellar mass functions.", All.smffile_name);
			endrun(buf);
		}
		//Skip the header
		fgets(line,linesize,ifp);
		//Initialise counter
		ismf = iset = 0;
		//Loop over all data sets
		while (fgets(line,linesize,ifp))
		{
			//Check if this is a new header line
			if (line[0] == '#')
			{
				//Read the header of each set
				sscanf(line,"# %d %f %f %f %f %[^\n]\n",&SmfSet[iset].ndata,&SmfSet[iset].min,&SmfSet[iset].max,&imf,&hubble,buffer);
				strcpy(SmfSet[iset].tag,buffer);
				//Read all data points in this set and compute sigma and correct for imf / hubble
				for (i = 0; i < SmfSet[iset].ndata; i++)
				{
					fgets(line,linesize,ifp);
					sscanf(line,"%f %f %f %f %d\n",&Smf[ismf].obs_x,&Smf[ismf].obs_y,&sigmaplus,&sigmaminus,&Smf[ismf].nfit);
					Smf[ismf].obs_sigma  =  0.5 * fabsf(sigmaplus - sigmaminus);
					Smf[ismf].bin        =  0.5 * (1.0/(SmfSet[iset].min+1.0) + 1.0/(SmfSet[iset].max+1.0));
					Smf[ismf].obs_x     += -imf - 2.0 * log10(All.h_100/hubble);
					Smf[ismf].obs_y     +=  3.0 * log10(All.h_100/hubble);
					ismf++;
				}//End set
				iset++;
			}//End header
		}//End reading lines
		//Close file
		fclose(ifp);
		efree(buffer);

		//Set the offset for each data set
		if (All.Nsmfset > 0) SmfSet[0].offset = 0;
		for (iset = 1; iset < All.Nsmfset; iset++) SmfSet[iset].offset = SmfSet[iset-1].offset + SmfSet[iset-1].ndata;

#ifdef GLOBAL_SIGMA_SMF
		//Add global error to SMF
		for (i = 0; i < All.Nsmf; i++)
		{
			if (Smf[i].bin >= 1./(SIGMASMFZTHRESH+1.0)) Smf[i].obs_sigma = sqrt(Smf[i].obs_sigma*Smf[i].obs_sigma+(float)(All.GlobalSigmaSmfLz*All.GlobalSigmaSmfLz));
			else Smf[i].obs_sigma = sqrt(Smf[i].obs_sigma*Smf[i].obs_sigma+(float)(All.GlobalSigmaSmfHz*All.GlobalSigmaSmfHz));
		}
#endif
	}

	//Broadcast to other tasks
	MPI_Bcast(Smf, All.Nsmf * sizeof(struct galaxy_data), MPI_BYTE, 0, MPI_COMM_WORLD);
	MPI_Bcast(SmfSet, All.Nsmfset * sizeof(struct data_set), MPI_BYTE, 0, MPI_COMM_WORLD);
}
#endif


#ifdef READ_FQ
/*! \brief This function reads all observed quenched fractions from a file
 *
 *  First the galaxy_data structure #fq is allocated. The function then reads all lines and if it is
 *  a header line (#) then the redshift and the corrections for IMF and Hubble_h are read. Then all
 *  data points in this set are read and corrected.
 */
void read_fq(void)
{
	FILE *ifp;
	const int linesize = 500;
	char line[linesize], buf[linesize], *buffer;
	int i,ifq,iset,nset;
	float imf,hubble;

	//Task 0 reads all data
	if (ThisTask == 0)
	{
		//Print what is done...
		printf("%s Reading observed quenched fractions.\n",All.startline);

		//Open file
		if(!(ifp=fopen(All.fqfile_name,"r")))
		{
			sprintf(buf, "Cannot open file `%s' for reading the quenched fractions.", All.fqfile_name);
			endrun(buf);
		}
		//Skip the header
		fgets(line,linesize,ifp);
		//Initialise counter
		All.Nfqset = 0;
		All.Nfq    = 0;
		while (fgets(line,linesize,ifp))
		{
			if (line[0] == '#')
			{
				sscanf(line,"# %d\n",&nset);
				All.Nfqset++;
				All.Nfq += nset;
			}
		}
		//Close file
		fclose(ifp);
	}

	MPI_Bcast(&All.Nfq, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&All.Nfqset, 1, MPI_INT, 0, MPI_COMM_WORLD);
	FqSet = emalloc("FQSET", All.Nfqset * sizeof(struct data_set));
	Fq    = emalloc("FQ",    All.Nfq    * sizeof(struct galaxy_data));

	//Task 0 reads all data
	if (ThisTask == 0)
	{
		buffer = emalloc("BUFFER", SSTRING * sizeof(char));

		//Open file
		if(!(ifp=fopen(All.fqfile_name,"r")))
		{
			sprintf(buf, "Cannot open file `%s' for reading the quenched fractions.", All.fqfile_name);
			endrun(buf);
		}
		//Skip the header
		fgets(line,linesize,ifp);
		//Initialise counter
		ifq = iset = 0;
		//Loop over all data sets
		while (fgets(line,linesize,ifp))
		{
			//Check if this is a new header line
			if (line[0] == '#')
			{
				//Read the header of each set
				sscanf(line,"# %d %f %f %f %f %[^\n]\n",&FqSet[iset].ndata,&FqSet[iset].min,&FqSet[iset].max,&imf,&hubble,buffer);
				strcpy(FqSet[iset].tag,buffer);
				//Read all data points in this set and compute sigma and correct for imf / hubble
				for (i = 0; i < FqSet[iset].ndata; i++)
				{
					fgets(line,linesize,ifp);
					sscanf(line,"%f %f %f %d\n",&Fq[ifq].obs_x,&Fq[ifq].obs_y,&Fq[ifq].obs_sigma,&Fq[ifq].nfit);
					Fq[ifq].bin        =  0.5 * (1.0/(FqSet[iset].min+1.0) + 1.0/(FqSet[iset].max+1.0));
					Fq[ifq].obs_x     += -imf - 2.0 * log10(All.h_100/hubble);
					ifq++;
				}//End set
				iset++;
			}//End header
		}//End reading
		//Close file
		fclose(ifp);
		efree(buffer);

		//Set the offset for each data set
		if (All.Nfqset > 0) FqSet[0].offset = 0;
		for (iset = 1; iset < All.Nfqset; iset++) FqSet[iset].offset = FqSet[iset-1].offset + FqSet[iset-1].ndata;

#ifdef GLOBAL_SIGMA_FQ
		for (i = 0; i < All.Nfq; i++)
			Fq[i].obs_sigma = sqrt(Fq[i].obs_sigma*Fq[i].obs_sigma+(float)(All.GlobalSigmaFq*All.GlobalSigmaFq));
#endif
	}

	//Broadcast to other tasks
	MPI_Bcast(Fq, All.Nfq * sizeof(struct galaxy_data), MPI_BYTE, 0, MPI_COMM_WORLD);
	MPI_Bcast(FqSet, All.Nfqset * sizeof(struct data_set), MPI_BYTE, 0, MPI_COMM_WORLD);
}
#endif


#ifdef READ_CSFRD
/*! \brief This function reads all observed cosmic star formation rate densities from a file
 *
 *  First the galaxy_data structure #csfrd is allocated. The function then reads all lines and if it is
 *  a header line (#) then the corrections for IMF and Hubble_h are read. Then all data points in this
 *  set are read and corrected.
 */
void read_csfrd(void)
{
	FILE *ifp;
	const int linesize = 500;
	char line[linesize], buf[linesize], *buffer;
	int i,icsfrd,iset,nset;
	float imf,hubble,sigmaplus,sigmaminus;

	//Task 0 reads all data
	if (ThisTask == 0)
	{
		//Print what is done...
		printf("%s Reading observed cosmic star formation rate densities.\n",All.startline);

		//Open file
		if(!(ifp=fopen(All.csfrdfile_name,"r")))
		{
			sprintf(buf, "Cannot open file `%s' for reading the cosmic star formation rate density.", All.csfrdfile_name);
			endrun(buf);
		}
		//Skip the header
		fgets(line,linesize,ifp);
		//Initialise counter
		All.Ncsfrdset = 0;
		All.Ncsfrd    = 0;
		while (fgets(line,linesize,ifp))
		{
			if (line[0] == '#')
			{
				sscanf(line,"# %d\n",&nset);
				All.Ncsfrdset++;
				All.Ncsfrd += nset;
			}
		}
		//Close file
		fclose(ifp);
	}

	MPI_Bcast(&All.Ncsfrd, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&All.Ncsfrdset, 1, MPI_INT, 0, MPI_COMM_WORLD);
	CsfrdSet = emalloc("CSFRDSET", All.Ncsfrdset * sizeof(struct data_set));
	Csfrd    = emalloc("CSFRD",    All.Ncsfrd    * sizeof(struct galaxy_data));

	//Task 0 reads all data
	if (ThisTask == 0)
	{
		buffer = emalloc("BUFFER", SSTRING * sizeof(char));

		//Open file
		if(!(ifp=fopen(All.csfrdfile_name,"r")))
		{
			sprintf(buf, "Cannot open file `%s' for reading the cosmic star formation rate density.", All.csfrdfile_name);
			endrun(buf);
		}
		//Skip the header
		fgets(line,linesize,ifp);
		//Initialise counter
		icsfrd = iset = 0;
		//Loop over all data sets
		while (fgets(line,linesize,ifp))
		{
			//Check if this is a new header line
			if (line[0] == '#')
			{
				//Read the header of each set
				sscanf(line,"# %d %f %f %[^\n]\n",&CsfrdSet[iset].ndata,&imf,&hubble,buffer);
				strcpy(CsfrdSet[iset].tag,buffer);
				//Read all data points in this set and compute sigma and correct for imf / hubble
				for (i = 0; i < CsfrdSet[iset].ndata; i++)
				{
					fgets(line,linesize,ifp);
					sscanf(line,"%f %f %f %f %d\n",&Csfrd[icsfrd].obs_x,&Csfrd[icsfrd].obs_y,&sigmaplus,&sigmaminus,&Csfrd[icsfrd].nfit);
					Csfrd[icsfrd].obs_sigma  =  0.5 * fabsf(sigmaplus + sigmaminus);
					Csfrd[icsfrd].obs_y     += -imf + log10(All.h_100/hubble/All.m_unit*All.t_unit*All.x_unit*All.x_unit*All.x_unit);
					icsfrd++;
				}//End set
				iset++;
			}//End header
		}//End reading
		//Close file
		fclose(ifp);
		efree(buffer);

		//Set the offset for each data set
		if (All.Ncsfrdset > 0) CsfrdSet[0].offset = 0;
		for (iset = 1; iset < All.Ncsfrdset; iset++) CsfrdSet[iset].offset = CsfrdSet[iset-1].offset + CsfrdSet[iset-1].ndata;

#ifdef GLOBAL_SIGMA_CSFRD
		for (i = 0; i < All.Ncsfrd; i++)
			Csfrd[i].obs_sigma = sqrt(Csfrd[i].obs_sigma*Csfrd[i].obs_sigma+(float)(All.GlobalSigmaCsfrd*All.GlobalSigmaCsfrd));
#endif
	}

	//Broadcast to other tasks
	MPI_Bcast(Csfrd, All.Ncsfrd * sizeof(struct galaxy_data), MPI_BYTE, 0, MPI_COMM_WORLD);
	MPI_Bcast(CsfrdSet, All.Ncsfrdset * sizeof(struct data_set), MPI_BYTE, 0, MPI_COMM_WORLD);
}
#endif


#ifdef READ_SSFR
/*! \brief This function reads all observed specific star formation rates from a file
 *
 *  First the galaxy_data structure #ssfr is allocated. The function then reads all lines and if it is
 *  a header line (#) then the corrections for IMF and Hubble_h are read. Then all data points in this
 *  set are read and corrected. Stellar masses are stored in ssfr.bin.
 */
void read_ssfr(void)
{
	FILE *ifp;
	const int linesize = 500;
	char line[linesize], buf[linesize], *buffer;
	int i,issfr,iset,nset;
	float imf,hubble,sigmaplus,sigmaminus;

	//Task 0 reads all data
	if (ThisTask == 0)
	{
		//Print what is done...
		printf("%s Reading observed specific star formation rates.\n",All.startline);

		//Open file
		if(!(ifp=fopen(All.ssfrfile_name,"r")))
		{
			sprintf(buf, "Cannot open file `%s' for reading the cosmic star formation rate density.", All.ssfrfile_name);
			endrun(buf);
		}
		//Skip the header
		fgets(line,linesize,ifp);
		//Initialise counter
		All.Nssfrset = 0;
		All.Nssfr    = 0;
		while (fgets(line,linesize,ifp))
		{
			if (line[0] == '#')
			{
				sscanf(line,"# %d\n",&nset);
				All.Nssfrset++;
				All.Nssfr += nset;
			}
		}
		//Close file
		fclose(ifp);
	}

	MPI_Bcast(&All.Nssfr, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&All.Nssfrset, 1, MPI_INT, 0, MPI_COMM_WORLD);
	SsfrSet = emalloc("SSFRSET", All.Nssfrset * sizeof(struct data_set));
	Ssfr    = emalloc("SSFR",    All.Nssfr    * sizeof(struct galaxy_data));

	//Task 0 reads all data
	if (ThisTask == 0)
	{
		buffer = emalloc("BUFFER", SSTRING * sizeof(char));

		//Open file
		if(!(ifp=fopen(All.ssfrfile_name,"r")))
		{
			sprintf(buf, "Cannot open file `%s' for reading the cosmic star formation rate density.", All.ssfrfile_name);
			endrun(buf);
		}
		//Skip the header
		fgets(line,linesize,ifp);
		//Initialise counter
		issfr = iset = 0;
		//Loop over all data sets
		while (fgets(line,linesize,ifp))
		{
			//Check if this is a new header line
			if (line[0] == '#')
			{
				//Read the header of each set
				sscanf(line,"# %d %f %f %[^\n]\n",&SsfrSet[iset].ndata,&imf,&hubble,buffer);
				strcpy(SsfrSet[iset].tag,buffer);
				//Read all data points in this set and compute sigma and correct for imf / hubble
				for (i = 0; i < SsfrSet[iset].ndata; i++)
				{
					fgets(line,linesize,ifp);
					sscanf(line,"%f %f %f %f %f %d\n",&Ssfr[issfr].obs_x,&Ssfr[issfr].bin,&Ssfr[issfr].obs_y,&sigmaplus,&sigmaminus,&Ssfr[issfr].nfit);
					Ssfr[issfr].bin       += -imf - 2.0 * log10(All.h_100/hubble);
					Ssfr[issfr].obs_sigma  =  0.5 * fabsf(sigmaplus + sigmaminus);
					issfr++;
				}//End set
				iset++;
			}//End header
		}//End reading
		//Close file
		fclose(ifp);
		efree(buffer);

		//Set the offset for each data set
		if (All.Nssfrset > 0) SsfrSet[0].offset = 0;
		for (iset = 1; iset < All.Nssfrset; iset++) SsfrSet[iset].offset = SsfrSet[iset-1].offset + SsfrSet[iset-1].ndata;

#ifdef GLOBAL_SIGMA_CSFRD
		for (i = 0; i < All.Nssfr; i++)
			Ssfr[i].obs_sigma = sqrt(Ssfr[i].obs_sigma*Ssfr[i].obs_sigma+(float)(All.GlobalSigmaSsfr*All.GlobalSigmaSsfr));
#endif
	}

	//Broadcast to other tasks
	MPI_Bcast(Ssfr, All.Nssfr * sizeof(struct galaxy_data), MPI_BYTE, 0, MPI_COMM_WORLD);
	MPI_Bcast(SsfrSet, All.Nssfrset * sizeof(struct data_set), MPI_BYTE, 0, MPI_COMM_WORLD);
}
#endif


#ifdef READ_WP
/*! \brief This function reads all observed projected correlation functions from a file
 *
 *  First the function reads all header lines (#) including the redshift for all correlation functions.
 *  Then the galaxy_data structure #Wp and the data_set structure #WpSet are allocated. The function then
 *  reads all lines and if it is a header line (#) then the minimum and maximum bin masses, and the
 *  corrections for IMF and Hubble_h are read. Then all data points in this set are read and corrected.
 */
void read_wp(void)
{
	FILE *ifp;
	const int linesize = 500;
	char line[linesize], buf[linesize], *buffer;
	int i,iwp,iwpbin,nset;
	float imf,hubble;

	//Task 0 first checks how many correlation functions there are
	if (ThisTask == 0)
	{
		//Print what is done...
		printf("%s Reading observed projected correlation functions.\n",All.startline);

		//Open file
		if(!(ifp=fopen(All.wpfile_name,"r")))
		{
			sprintf(buf, "Cannot open file `%s' for reading the projected correlation functions.", All.wpfile_name);
			endrun(buf);
		}

		//Skip the header
		fgets(line,linesize,ifp);
		//Read the global redshift
		fgets(line,linesize,ifp);
		sscanf(line,"%lf\n",&All.wpredshift);
		//Skip the header
		fgets(line,linesize,ifp);
		//Count how many different corerlation functions there are (lines starting with #)
		All.Nwpset = 0;
		All.Nwp    = 0;
		while (fgets(line,linesize,ifp))
		{
			if (line[0] == '#')
			{
				sscanf(line,"# %d\n",&nset);
				All.Nwpset++;
				All.Nwp += nset;
			}
		}
		//Close file
		fclose(ifp);
	}

	MPI_Bcast(&All.Nwp, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&All.Nwpset, 1, MPI_INT, 0, MPI_COMM_WORLD);
	WpSet = emalloc("WPBIN", All.Nwpset * sizeof(struct data_set));
	Wp    = emalloc("WP", All.Nwp    * sizeof(struct galaxy_data));

	//Task 0 reads all data
	if (ThisTask == 0)
	{
		buffer = emalloc("BUFFER", SSTRING * sizeof(char));

		//Open file
		if(!(ifp=fopen(All.wpfile_name,"r")))
		{
			sprintf(buf, "Cannot open file `%s' for reading the projected correlation functions.", All.wpfile_name);
			endrun(buf);
		}
		//Skip the header
		fgets(line,linesize,ifp);
		fgets(line,linesize,ifp);
		fgets(line,linesize,ifp);
		//Initialise counters
		iwp    = 0;
		iwpbin = 0;
		//Loop over all data sets
		while (fgets(line,linesize,ifp))
		{
			//Check if this is a new header line
			if (line[0] == '#')
			{
				//Read the header of each set
				sscanf(line,"# %d %f %f %f %f %f %[^\n]\n",&WpSet[iwpbin].ndata,&WpSet[iwpbin].min,&WpSet[iwpbin].max,&WpSet[iwpbin].cut,&imf,&hubble,buffer);
				strcpy(WpSet[iwpbin].tag,buffer);
				//Set the offset
				if (iwpbin == 0) WpSet[iwpbin].offset = 0;
				else WpSet[iwpbin].offset = WpSet[iwpbin-1].ndata + WpSet[iwpbin-1].offset;
				//Correct the mass bins for imf/hubble
				WpSet[iwpbin].min += -imf - 2.0 * log10(All.h_100/hubble);
				WpSet[iwpbin].max += -imf - 2.0 * log10(All.h_100/hubble);
				WpSet[iwpbin].min  =  pow(10.,WpSet[iwpbin].min)/All.m_unit;
				WpSet[iwpbin].max  =  pow(10.,WpSet[iwpbin].max)/All.m_unit;
				WpSet[iwpbin].cut *= hubble/All.h_100/All.x_unit;
				//Read all data points in this set and compute sigma and correct for imf / hubble
				for (i = 0; i < WpSet[iwpbin].ndata; i++)
				{
					fgets(line,linesize,ifp);
					sscanf(line,"%f %f %f %d\n",&Wp[iwp].obs_x,&Wp[iwp].obs_y,&Wp[iwp].obs_sigma,&Wp[iwp].nfit);
					Wp[iwp].obs_x     *= hubble/All.h_100/All.x_unit;
					Wp[iwp].obs_y     *= hubble/All.h_100/All.x_unit;
					Wp[iwp].obs_sigma *= hubble/All.h_100/All.x_unit;
					iwp++;
				}//End set
				iwpbin++;
			}//End header
		}//End reading
		//Close file
		fclose(ifp);
		efree(buffer);

#ifdef GLOBAL_SIGMA_WP
		for (i = 0; i < All.Nwp; i++)
			Wp[i].obs_sigma = sqrt(Wp[i].obs_sigma*Wp[i].obs_sigma+(float)(All.GlobalSigmaWp*All.GlobalSigmaWp)*Wp[i].obs_y*Wp[i].obs_y);
#endif
		//Identify the minimum and maximum mass of the correlation functions
		All.wpmmin = WpSet[0].min;
		All.wpmmax = WpSet[0].max;
		for (i = 0; i < All.Nwpset; i++)
		{
			if (WpSet[i].min < All.wpmmin) All.wpmmin = WpSet[i].min;
			if (WpSet[i].max > All.wpmmax) All.wpmmax = WpSet[i].max;
		}
	}

	//Broadcast to other tasks
	MPI_Bcast(&All.wpredshift, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&All.wpmmin, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&All.wpmmax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(Wp, All.Nwp * sizeof(struct galaxy_data), MPI_BYTE, 0, MPI_COMM_WORLD);
	MPI_Bcast(WpSet, All.Nwpset * sizeof(struct data_set), MPI_BYTE, 0, MPI_COMM_WORLD);
}
#endif
