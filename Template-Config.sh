####################################################################################################
#!/bin/bash            # this line only there to enable syntax highlighting in this file
####################################################################################################
#  EMERGE Config file - Enable/Disable compile-time options as needed                              #
####################################################################################################
#
#
# --------------------------------------------------------------------------------------------------
#                            CODE OPTIONS
# --------------------------------------------------------------------------------------------------
#OPENMPTHREADS=2              # Enables OpenMP support and sets the number of threads per task
RANDOM_NUMBER_TABLE=1000000  # Enables the random number table and sets its size
#LONGIDS                      # Use 64 bit integers for IDs (default is 32 bit) for large runs
#DISABLE_MEMORY_MANAGER       # Disables the default memory manager
#
# --------------------------------------------------------------------------------------------------
#                            MODEL OPTIONS
# --------------------------------------------------------------------------------------------------
SFE_MPEAK_ZEVOLV             # Use redshift evolution for M1
SFE_NORM_ZEVOLV              # Use redshift evolution for Epsilon
SFE_BETA_ZEVOLV              # Use redshift evolution for Beta
#SFE_GAMMA_ZEVOLV             # Use redshift evolution for Gamma
#SFE_SAME_SLOPE               # Use the same slope for Beta and Gamma
SAT_QUENCH_MASS_DEPENDENT    # Use mass dependent quenching timescale for sats
#SAT_STRIP_USE_MSTAR          # Use fraction of mstar in stripping
#SAT_SFR_EXP_DECAY            # Exponential decline in str after peak mass
#COMPUTE_ICM                  # Compute the ICM for each halo
DF_USE_BK08                  # Use the dynamical friction fitting function by Boylan-Kolchin (2008)
ORPHAN_MASSLOSS              # Let the halo mass of orphans decay on merging timescale
#ORPHAN_NONRADIAL_INFALL      # Recompute the position angles of orphans at each time step
#
# --------------------------------------------------------------------------------------------------
#                            DATA OPTIONS
# --------------------------------------------------------------------------------------------------
READ_SMF                     # Read the stellar mass functions and use in fit
READ_FQ                      # Read the quenched fractions and use in fit
READ_CSFRD                   # Read the cosmic SFR density and use in fit
READ_SSFR                    # Read the specific SFRs and use in fit
READ_WP                      # Read the projected correlation functions and use in fit
GLOBAL_SIGMA_SMF             # Specify a global error for the SMF and add it
GLOBAL_SIGMA_FQ              # Specify a global error for the FQ and add it
GLOBAL_SIGMA_CSFRD           # Specify a global error for the CSFRD and add it
GLOBAL_SIGMA_SSFR            # Specify a global error for the SSFR and add it
GLOBAL_SIGMA_WP              # Specify a global error for the WP and add it
#
# --------------------------------------------------------------------------------------------------
#                            CORRELATION FUNCTION OPTIONS
# --------------------------------------------------------------------------------------------------
#WP_RBINS=20                  # Number of bins for which the 3d correlation function is computed
#WP_RBINS_INT=1000            # Number of interpolation bins for the 3d correlation function
#WP_RMAX=0.1                  # Maximum radius (fraction of box size) of the correlation function
#WP_NLEAF_MIN=4               # Minimum number of objects per kd-node (if more the node is split)
#WP_NODE_WIDTH_MIN=0.01       # Minimum kd-node size (fraction of box size - no more splitting)
#
# --------------------------------------------------------------------------------------------------
#                            FIT OPTIONS
# --------------------------------------------------------------------------------------------------
#HYBRID_ALPHA=0.4             # Hybrid optimization: slope for function g
#HYBRID_BETA=1.0              # Hybrid optimization: slope for function f (larger chi^2)
#HYBRID_GAMMA=1.0             # Hybrid optimization: slope for function f (smaller chi^2)
#PT_COLD_FRACTION=0.25        # Parallel Tempering: Fraction of cold walkers (T=1)
#PT_TARGET_ACCEPTANCE=0.3     # Parallel Tempering: Target acceptance rate
#PT_STEPS_SCALE_ADJUST=10     # Parallel Tempering: Chain steps for adjusting scale
#
# --------------------------------------------------------------------------------------------------
#                           OUTPUT OPTIONS
# --------------------------------------------------------------------------------------------------
HDF5_SUPPORT                 # Enables HDF5 support for the output files (Use OutputFormat = 2)
WRITE_GALAXY_CATALOG         # Writes galaxy catalogues at the specified redshifts
#WRITE_HALO_CATALOG           # Writes halo catalogues at the specified redshifts
#WRITE_MAINBRANCH             # Writes main branch history of specified systems
#
####################################################################################################
