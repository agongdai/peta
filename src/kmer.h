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
} mer;

void ext_by_kmers(char *lib_file, const char *solid_file,
		const char *kmer_file, const int insert_size, const int sd_insert_size,
		const int n_threads);
void build_kmers(const char *fa_fn, const char *out_fn, const int k);

#endif /* KMER_H_ */
