/*
 * kmer.h
 *
 *  Created on: 13-Mar-2013
 *      Author: carl
 */

#ifndef KMER_H_
#define KMER_H_

#include <stdint.h>

typedef struct {
	uint64_t s;
	uint32_t count;
	uint8_t status;
} mer;

typedef struct {
	uint64_t n_seqs;
	uint64_t n_kmers;
	int k;
} mer_meta;

void ext_by_kmers(char *lib_file, const char *solid_file,
		const char *kmer_file, const int insert_size, const int sd_insert_size,
		const int n_threads);
void build_kmers(const char *fa_fn, const char *out_fn, const int k);
void build_kmers_gvdb(const char *fa_fn, const char *out_fn, const int k);
GHashTable *load_kmers(const char *kmer_file, GPtrArray *kmer_list,
		mer_meta *meta);

#endif /* KMER_H_ */
