/*
 * kmers.h
 *
 *  Created on: 02-Apr-2013
 *      Author: carl
 */

#ifndef KMERS_H_
#define KMERS_H_

#include <unordered_map>
#include <stdint.h>

typedef struct mer {
	uint64_t s;
	uint32_t count;
	uint8_t status;
} mer;

typedef struct {
	uint64_t n_seqs;
	uint64_t n_kmers;
	int k;
} mer_meta;

#ifdef __cplusplus
extern "C" {
#endif

	void test_sdarray();
	void load_sda_kmers(const char *kmer_file, GPtrArray *kmer_list, mer_meta *meta);
	void build_kmers(const char *fa_fn, const char *out_fn, const int k);
	mer *new_mer();
	int pe_kmer(int argc, char *argv[]);

#ifdef __cplusplus
}
#endif

#endif /* KMERS_H_ */
