/*
 * kmer.h
 *
 *  Created on: 13-Mar-2013
 *      Author: carl
 */

#ifndef KMER_H_
#define KMER_H_

#include "bwase.h"

typedef struct {
	bwa_seq_t *s;
	int fre;
} mer;

void ext_by_kmers(char *lib_file, const char *solid_file,
		const char *kmer_file, const int insert_size, const int sd_insert_size,
		const int n_threads);

#endif /* KMER_H_ */
