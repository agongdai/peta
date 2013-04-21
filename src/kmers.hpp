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
#include <glib.h>
#include "pehash.h"
#include "bwtaln.h"

using namespace std;
typedef unordered_map<uint64_t, uint64_t*> mer_hash;
typedef unordered_map<uint64_t, uint64_t> mer_counter;

typedef struct {
	uint16_t k;
	uint64_t n_reads;
	uint64_t n_k_mers;
	uint64_t n_valid_k_mers;
	uint64_t n_pos;
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
	uint32_t count;
} kmer_counter;

#ifdef __cplusplus
extern "C" {
#endif

uint64_t rev_comp_kmer(uint64_t kmer, const int n);
void destroy_hm(hash_map *hm);
uint64_t get_kmer_int(const ubyte_t *seq, const int start,
		const int interleaving, const int len);
void build_kmers_hash(const char *fa_fn, const int k, const int with_reads);
hash_map *load_hash_map(const char *fa_fn, const int with_reads, mer_hash& kmers);
bwa_seq_t *get_kmer_seq(uint64_t kmer, const int k);
int next_char_by_kmers(mer_hash *kmers, const int k, bwa_seq_t *query,
		const int ori);
void test_kmer_hash(const char *fa_fn);
uint64_t get_kmer_count(const uint64_t kmer_int, hash_map *hm);
void mark_kmer_used(const uint64_t kmer_int, const hash_map *hm);
void kmer_aln_query(const bwa_seq_t *query, const hash_map *hm, GPtrArray *hits);

#ifdef __cplusplus
}
#endif

#endif /* KMERS_H_ */
