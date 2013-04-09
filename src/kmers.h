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
#include "pehash.h"
#include "bwtaln.h"

using namespace std;
typedef unordered_map<uint64_t, uint64_t*> mer_hash;

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

typedef struct {
	hash_opt *o;
	uint64_t n_seqs;
	bwa_seq_t *seqs;
	mer_hash *hash;
} hash_map;

typedef unordered_map<uint64_t, mer*> mer_map;
typedef unordered_map<uint64_t, uint32_t> mer_counter;

#ifdef __cplusplus
extern "C" {
#endif

	void test_sdarray();
	void build_kmers(const char *fa_fn, const char *out_fn, const int k);
	mer *new_mer();
	int pe_kmer(int argc, char *argv[]);

#ifdef __cplusplus
}
#endif

#endif /* KMERS_H_ */
