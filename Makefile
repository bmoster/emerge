#######################################################################
#                                                                     #
#   Makefile for the Empirical ModEl for the foRmation of GalaxiEs    #
#                                                                     #
#######################################################################
#                                                                     #
# To build the code, do the following:                                #
# (1) Copy the file "Template-Config.sh" to "Config.sh"               #
# (2) Edit "Config.sh" as needed for your application                 #
# (3) Copy the file "Template-TypeOfSystem.sh" to "TypeOfSystem.sh"   #
# (4) Uncomment your system in  "TypeOfSystem.sh"                     #
# (3) Run "make"                                                      #
#                                                                     #
#######################################################################
#                                                                     #
# New compile-time options should be added to the file                #
# "Template-Config.sh" only. Usually, they should be added there in   #
# the disabled/default version. New types of systems can be added     #
# below and should also be added to "Template-TypeOfSystem.sh" in     #
# the disabled version.                                               #
#                                                                     #
# Neither "Config.sh" nor "TypeOfSystem.sh" should be checked in to   #
# the repository.                                                     #
#                                                                     #
#######################################################################
#                                                                     #
# Note: It is possible to override the default name of the Config.sh  #
# file, if desired, as well as the name of the executable.            #
# For example:                                                        #
#                                                                     #
# make CONFIG=NewConfig.sh EXEC=emerge_new                            #
#                                                                     #
#######################################################################
#
#
#
#
# ----------------------------------------------------------------
EXEC         = emerge
CONFIG       = Config.sh
BUILD_DIR    = build
SRC_DIR      = src
TOOLS_DIR    = tools
# ----------------------------------------------------------------
include TypeOfSystem.sh
# ----------------------------------------------------------------
MAKEFILES    = Makefile config-makefile TypeOfSystem.sh
PERL	     = /usr/bin/perl
RESULT     := $(shell CONFIG=$(CONFIG) PERL=$(PERL) BUILD_DIR=$(BUILD_DIR) make -f config-makefile)
CONFIGVARS := $(shell cat $(BUILD_DIR)/codeoptions.h)
# ----------------------------------------------------------------
#
#
#
#######################################################################
# Default options                                                     #
#######################################################################
CC       = mpicc
MPICHLIB = -lmpich
GSLLIB   = -lgsl -lgslcblas
MATHLIB  = -lm
#######################################################################
#
#
#
#######################################################################
# Define the types of systems here
#######################################################################
#
# MacBook with OS X Mojave, OpenMPI 3.1.1, and GSL 2.5
#
ifeq ($(TYPEOFSYSTEM),"MacBook")
CC        = /usr/local/bin/mpicc
OPTIMIZE  = -O3 -g -Wall -Wno-uninitialized -Wformat-overflow=0
ifeq (OPENMPTHREADS,$(findstring OPENMPTHREADS,$(CONFIGVARS)))
OPTIMIZE += -fopenmp
OMP_LIBS  = -lgomp
else
OPTIMIZE += -Wno-unknown-pragmas
endif
GSL_INCL  =
GSL_LIBS  =
MPICHLIB  =
HDF5_INCL =
HDF5_LIB  = -lhdf5 -lz
endif
#
# ---------------------------------------------------------------------
#
# ODIN @ MPA
#
# module load gsl/2.1
# module load impi
# module load hdf5-serial
#
ifeq ($(TYPEOFSYSTEM),"Odin")
CC        = mpiicc
OPTIMIZE  = -O3 -g -m64
ifeq (OPENMPTHREADS,$(findstring OPENMPTHREADS,$(CONFIGVARS)))
OPTIMIZE += -openmp
OMP_LIBS  = -openmp
else
OPTIMIZE += -Wno-unknown-pragmas
endif
GSL_INCL  = -I$(GSL_HOME)/include/
GSL_LIBS  = -L$(GSL_HOME)/lib/ -Xlinker -R -Xlinker $(GSL_HOME)/lib
MPICHLIB  =
HDF5_INCL = -I$(HDF5_HOME)/include
HDF5_LIB  = -L$(HDF5_HOME)/lib -Xlinker -R -Xlinker $(HDF5_HOME)/lib -lhdf5 -lz
OPT      +=  -DNOCALLSOFSYSTEM
endif
#
# ---------------------------------------------------------------------
#
# FREYA @ MPA
#
# module load gsl/2.2
# module load impi
# module load hdf5-serial
#
ifeq ($(TYPEOFSYSTEM),"Freya")
CC        = mpiicc
OPTIMIZE  = -O3 -g -m64
ifeq (OPENMPTHREADS,$(findstring OPENMPTHREADS,$(CONFIGVARS)))
OPTIMIZE += -qopenmp
OMP_LIBS  = -qopenmp
else
OPTIMIZE += -Wno-unknown-pragmas
endif
GSL_INCL  = -I$(GSL_HOME)/include/
GSL_LIBS  = -L$(GSL_HOME)/lib/ -Xlinker -R -Xlinker $(GSL_HOME)/lib
MPICHLIB  =
HDF5_INCL = -I$(HDF5_HOME)/include
HDF5_LIB  = -L$(HDF5_HOME)/lib -Xlinker -R -Xlinker $(HDF5_HOME)/lib -lhdf5 -lz
OPT      +=  -DNOCALLSOFSYSTEM
endif
#
# ---------------------------------------------------------------------
#
# LRZ: SuperMuc / LinuxCluster
#
# module load gsl/2.4
# module load hdf5/serial
#
ifeq ($(TYPEOFSYSTEM),"LRZ")
CC        = mpiicc
OPTIMIZE  = -O3 -g -m64
ifeq (OPENMPTHREADS,$(findstring OPENMPTHREADS,$(CONFIGVARS)))
OPTIMIZE += -qopenmp
OMP_LIBS  = -qopenmp
else
OPTIMIZE += -Wno-unknown-pragmas
endif
GSL_INCL  = $(GSL_INC)
GSL_LIBS  = -L$(GSL_LIBDIR)
MPICHLIB  =
HDF5_INCL = $(HDF5_INC)
HDF5_LIB  = $(HDF5_SHLIB) -lhdf5 -lz
OPT      +=  -DNOCALLSOFSYSTEM
endif
#
# ---------------------------------------------------------------------
#
# Dorc @ USM
#
# module load intel
# module load mpi.intel
#
ifeq ($(TYPEOFSYSTEM),"Dorc")
CC        = mpiicc
OPTIMIZE  = -O3 -g -Wall -Wno-uninitialized
ifeq (OPENMPTHREADS,$(findstring OPENMPTHREADS,$(CONFIGVARS)))
OPTIMIZE += -fopenmp
OMP_LIBS  = -lgomp
else
OPTIMIZE += -Wno-unknown-pragmas
endif
GSL_INCL  = -I/usr/local/include/
GSL_LIBS  = -L/usr/local/lib64/
MPICHLIB  =
HDF5_INCL =
HDF5_LIB  = -lhdf5 -lz
endif
#
# ---------------------------------------------------------------------
#
# C2PAP
#
# module load gsl/2.4
# module load hdf5/serial
#
ifeq ($(TYPEOFSYSTEM),"C2PAP")
CC        = mpicc
OPTIMIZE  = -O3 -g -m64
ifeq (OPENMPTHREADS,$(findstring OPENMPTHREADS,$(CONFIGVARS)))
OPTIMIZE += -qopenmp
OMP_LIBS  = -qopenmp
else
OPTIMIZE += -Wno-unknown-pragmas
endif
GSL_INCL  = $(GSL_INC)
GSL_LIBS  = -L$(GSL_LIBDIR)
MPICHLIB  =
HDF5_INCL = $(HDF5_INC)
HDF5_LIB  = $(HDF5_SHLIB) -lhdf5 -lz
OPT      +=  -DNOCALLSOFSYSTEM
endif
#
#######################################################################
# Determine libraries that are not needed
#
ifneq (HDF5_SUPPORT,$(findstring HDF5_SUPPORT,$(CONFIGVARS)))
HDF5_LIB  =
endif
#
#######################################################################
#
#
#
#
#######################################################################
# Define LINKER
#######################################################################
ifndef LINKER
LINKER = $(CC)
endif
#######################################################################
#
#
#
#
#######################################################################
# Object and header files
#######################################################################
OBJS = allocate.o allvars.o clustering.o fit.o functions.o galaxies.o \
  main.o output.o read_data.o read_trees.o setup.o statistics.o

INCL += allvars.h proto.h



#######################################################################
# Compile options
#######################################################################

CCFLAGS = $(OPTIMIZE) $(OPT) $(GSL_INCL) $(HDF5_INCL) -I$(BUILD_DIR)
LDFLAGS = $(MATHLIB) $(OMP_LIBS) $(MPICHLIB) $(GSL_LIBS) $(GSLLIB) $(HDF5_LIB)

OBJS := $(addprefix $(BUILD_DIR)/,$(OBJS)) $(BUILD_DIR)/compile_info.o
INCL := $(addprefix $(SRC_DIR)/,$(INCL)) $(BUILD_DIR)/codeoptions.h



#######################################################################
# Build rules
#######################################################################

all: $(EXEC)

$(EXEC): $(OBJS)
	$(LINKER) $(OBJS) $(LDFLAGS) -o $(EXEC)
	$(CC) $(TOOLS_DIR)/convert_CT_to_emerge.c -o $(TOOLS_DIR)/convert_CT_to_emerge -lm

clean:
	rm -f $(OBJS) $(EXEC)
	rm -f $(BUILD_DIR)/compile_info.c $(BUILD_DIR)/codeoptions.h
	rm -f $(TOOLS_DIR)/convert_CT_to_emerge

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c $(INCL) $(MAKEFILES)
	$(CC) $(CCFLAGS) -c $< -o $@

$(BUILD_DIR)/compile_info.o: $(BUILD_DIR)/compile_info.c $(MAKEFILES)
	$(CC) $(CCFLAGS) -c $< -o $@
