#ifndef CLEAN_H_
#define CLEAN_H_
#include <stdint.h>
#include <inttypes.h>
#include <glib.h>
#include "bwase.h"

#define MAX_N16			4294967295
#define LOW_KMER		2
#define UNEVEN_THRE		1
#define LOW_RANGE		1.2
#define MAX_TIME	    128	
#define SOLID_PERCERN	0.15

typedef struct {
	int kmer;
	int mode;
	double stop_thre;
	char *lib_name;
	int n_threads;
} clean_opt;

typedef struct {
	int key;
	uint16_t count;
} kmer;

typedef struct {
	uint32_t read_id;
	double k_freq;
	double k_sd;
	double k_mean;
	uint8_t checked;
} counter;

typedef struct {
	uint32_t start;
	uint32_t end;
	int rm_low_kmer;
	bwa_seq_t *seqs;
	int n_seqs;
	counter *counter_list;
	uint16_t *kmer_list;
	clean_opt *opt;
} clean_thread_aux;

typedef struct {
	uint32_t start;
	uint32_t end;
	bwa_seq_t *seqs;
	int n_seqs;
	counter	*sorted_counters;
	uint16_t *kmer_list;
	int *n_solid;
	int n_needed;
	int n_iterate;
	GPtrArray *solid_reads;
	clean_opt *opt;
} iterate_thread_aux;

int clean_reads(int argc, char *argv[]);
clean_opt *init_clean_opt();
void set_kmer_index(const bwa_seq_t *read, int k, uint16_t *kmer_list);
int cmp_kmer(const void *a, const void *b);
void set_k_freq(bwa_seq_t *read, counter *k_count, uint16_t *kmer_list,
		const int k);
GPtrArray *calc_solid_reads(bwa_seq_t *seqs, const int n_seqs, clean_opt *opt,
		const int n_needed, const int by_coverage, const int rm_low_kmer);
void pe_clean_core(char *fa_fn, clean_opt *opt);
#endif
