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
#include "bwtaln.h"
#include "utils.h"

#define N_POS_BITS			16
#define HASH_VALUE_HIGHER	281474976710655 // 48 1's
#define HASH_VALUE_LOWER	65535

typedef uint64_t hash_key;
typedef uint64_t hash_value;
typedef int16_t read_pos;

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

	hash_key get_hash_key(const ubyte_t *seq, const int start,
			const int interleaving, const int len);
	uint64_t get_kmer_int(const ubyte_t *seq, const int start,
				const int interleaving, const int len);
	hash_value get_hash_value(const index64 seq_id, const int pos_start);
	void read_hash_value(index64 *seq_id, int *pos_start, hash_value value);
	void destroy_hm(hash_map *hm);
	void build_kmers_hash(const char *fa_fn, const int k, const int with_reads);
	hash_map *load_hash_map(const char *fa_fn, const int with_reads, mer_hash& kmers);
	void test_kmer_hash(const char *fa_fn);
	int build_kmer_hash(int argc, char *argv[]);


#ifdef __cplusplus
}
#endif

#endif /* HASH_H_ */
