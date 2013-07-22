/*
 * hash.h
 *
 *  Created on: 09-May-2013
 *      Author: carl
 */

#ifndef HASH_H_
#define HASH_H_

#include <unordered_map>
#include <stdint.h>
#include <inttypes.h>
#include <glib.h>
#include "bwtaln.h"
#include "utils.h"
#include "k_hash.h"

using namespace std;
typedef unordered_map<uint64_t, uint64_t*> mer_hash;
typedef unordered_map<uint64_t, uint64_t> mer_counter;

typedef struct {
	uint16_t k;
	uint64_t n_reads;
	uint64_t n_k_mers;
	uint64_t n_valid_k_mers;
	uint64_t n_pos;
	uint32_t read_len;
} map_opt;

typedef struct {
	map_opt *o;
	uint64_t n_reads;
	uint64_t *kmers_ordered;
	bwa_seq_t *seqs;
	mer_hash *hash;
} hash_map;

typedef struct {
	uint64_t kmer;
	uint64_t count;
} kmer_counter;

#ifdef __cplusplus
extern "C" {
#endif

	void destroy_hm(hash_map *hm);
	void build_kmers_hash(const char *fa_fn, const int k, const int with_reads);
	void build_tpl_hash(mer_hash& hash, GPtrArray *tpls, const int k);
	hash_map *load_hash_map(const char *fa_fn, const int with_reads, mer_hash& kmers);
	void test_kmer_hash(const char *fa_fn);
	int build_kmer_hash(int argc, char *argv[]);
	int export_frequency(int argc, char *argv[]);


#ifdef __cplusplus
}
#endif

#endif /* HASH_H_ */
