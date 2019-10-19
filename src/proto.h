///////////////////////////////////////////////////////////////////////////////////////////////////
// Emerge code - File proto.h                                                                    //
///////////////////////////////////////////////////////////////////////////////////////////////////
///
/// \file proto.h
/// \brief Contains the declarations of all functions
///
/// This file declares all functions in the code. All new functions should be added here as well.
/// The function headers should be declared in the block that corresponds to the file the function
/// is in.
///
///////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef PROTO
#define PROTO
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "allvars.h"

#if defined(_OPENMP) && (OPENMPTHREADS > 1)
#include <omp.h>
#endif


//allocate.c
void malloc_init(void);
void report_detailed_memory_usage_of_largest_task(size_t * OldHighMarkBytes, const char *label, const char *func, const char *file, int line);
void dump_memory_table(void);
int dump_memory_table_buffer(char *p);
void *malloc_fullinfo(const char *varname, size_t n, const char *func, const char *file, int line);
void *malloc_movable_fullinfo(void *ptr, const char *varname, size_t n, const char *func, const char *file, int line);
void free_fullinfo(void *p, const char *func, const char *file, int line);
void free_movable_fullinfo(void *p, const char *func, const char *file, int line);
void *realloc_fullinfo(void *p, size_t n, const char *func, const char *file, int line);
void *realloc_movable_fullinfo(void *p, size_t n, const char *func, const char *file, int line);

//clustering.c
int compute_wp(void);
int compute_wp_bin(int iwpbin, int Ngaltot, struct galaxy *localgal);
int collect_galaxies(int iwpbin, int Ngaltot, struct galaxy *localgal);
int construct_tree(struct kd_tree *tree, struct kd_pos *p, int count, int dim);
int find_nodes(struct kd_tree tree, int *nodeindex, int inode, int l, int level, int *counter, int save);
int partition(struct kd_pos *p, int left, int right, double pivotValue, int dim);
float mindist2_node_node(struct kd_node node1, struct kd_node node2, int dim);
float maxdist2_node_node(struct kd_node node1, struct kd_node node2, int dim);
int dualtree(struct kd_tree tree, struct kd_node node1, struct kd_node node2, struct kd_pos *p, float rmin2, float rmax2);
int tree_pairs(struct kd_tree tree, struct kd_pos *p, float rmin, float rmax);
void processRequests(int *listpos, int *nterm, int W);

//compile_info.c
void output_code_options(void);

//fit.c
void aies_mcmc(int restart);
void optimize_hybrid(int restart);
void parallel_tempering(int restart);
void print_mcmc_header(void);
void init_param(double *p);
void init_sigma(double *s);
float lnprob(double *p);
float print_chi2(void);
void get_chi2(float *chi2);
float chi2_smf(void);
float chi2_fq(void);
float chi2_csfrd(void);
float chi2_ssfr(void);
float chi2_wp(void);

//functions.c
int compare_float(const void *a, const void *b);
int compare_id(const void *a, const void *b);
int compare_load(const void *a, const void *b);
int compare_scale(const void *a, const void *b);
int compare_totask(const void *a, const void *b);
int compare_forest(const void *a, const void *b);
int compare_size(const void *a, const void *b);
int binary_search_id(struct search_list_IDs *a, int n, IDType key);
int binary_search_float(float *a, int n, float key);
float cosmictime(float a);
float Omega(float z);
float Epeebles(float z);
float tdyn(float a);
float get_gaussian_random_number(int index);
float get_uniform_random_number(int index);
double second(void);

//galaxies.c
void init_galaxies(void);
void finish_galaxies(void);
void make_galaxies(void);
void reset_galaxies(void);
void assemble_galaxies(void);
void assemble_tree(int itree, int thisthread);
float sfe(float hmass, float a);
float stellar_mass(int imasstab, int scalemax, int thisthread);
float intra_cluster_mass(int imass, int scalemax, int thisthread);
float get_tdf(int ihalo, int imain);

//output.c
void output_galaxies(void);
void output_halos(void);
void output_mainbranch(void);

//read_data.c
void read_data(void);
void read_smf(void);
void read_fq(void);
void read_csfrd(void);
void read_ssfr(void);
void read_wp(void);

//read_trees.c
int read_trees(char *fname);
int find_files(char *fname);
void distribute_file(int ntask, int nfiles, int firstfile, int firsttask, int lasttask, int *filenr, int *master, int *last);
void share_halo_number_in_file(const char *fname, int readTask, int lastTask);
void read_tree_file(char *fname, int readTask, int lastTask);
void empty_read_buffer(int offset, int hc);
void find_forests(void);
void sort_by_size(void);
void distribute_trees(int ntask, double load);
void get_timesteps(void);
void add_orphans(void);
void setup_haloes_by_tree(void);
void setup_haloes_by_forest(void);
void compute_halo_history(void);
void find_parents(void);
void copy_trees_to_other_universes(int ntask, int nuniverses);

//setup.c
void read_parameterfile(char *fname);
void setup(void);
void finalize (void);
void open_logfiles(void);
void close_logfiles(void);
void init_random_numbers(void);
void print_banner(void);

//statistics.c
void add_galaxy_to_statistics(int ihalo, int thisthread);
void get_statistics(void);
void write_statistics(void);
void write_statistics_init_ascii(void);
void write_statistics_init_hdf5(void);
void write_statistics_ascii(char *outdir, int iuniverse);
void write_statistics_hdf5(char *outfname, int iuniverse);
void write_statistics_chi2(char *outdir, int step);
void write_statistics_chi2_ascii(char *outdir, int step);
void write_walker_statistics(void);
void write_walker_statistics_init_ascii(char *outfname);
void write_walker_statistics_init_hdf5(char *outfname);
void write_walker_statistics_ascii(char *outfname, int iuniverse);
void write_walker_statistics_hdf5(char *outfname, int iuniverse);

#endif
