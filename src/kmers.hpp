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
#include "bwtaln.h"
#include "hash.hpp"

#define		LOWER_ONES_32		4294967295				// 32 0's, followed by 32 1's
#define		LOWER_ONES_24		16777215				// 40 0's, followed by 24 1's
#define		LOWER_ONES_16		65535					// 48 0's, 16 11's

using namespace std;

#ifdef __cplusplus
extern "C" {
#endif

	int *count_next_kmers(hash_map *hm, uint64_t query, const int fresh_only, const int ori);
	uint64_t rev_comp_kmer(uint64_t kmer, const int n);
	bwa_seq_t *get_kmer_seq(uint64_t kmer, const int k);
	GPtrArray *kmer_find_reads(const bwa_seq_t *query, const hash_map *hm,
			const int mismatch, const uint32_t n_part_only);
	int next_char_by_kmers(hash_map *hm, uint64_t kmer_int, const int fresh_only, const int ori);
	int next_char_max_freq(hash_map *hm, uint64_t kmer_int, const int fresh_only,
			const int ori);
	uint64_t get_kmer_count(const uint64_t kmer_int, hash_map *hm, const int fresh);
	uint64_t get_kmer_rf_count(const uint64_t query, hash_map *hm, const int fresh_only);
	void mark_kmer_used(const uint64_t kmer_int, const hash_map *hm, const int tpl_id, const int locus, const int tpl_len);
	void mark_kmer_not_used(const uint64_t kmer_int, const hash_map *hm);
	int kmer_is_used(const uint64_t kmer_int, hash_map *hm);
	void read_tpl_using_kmer(const uint64_t kmer_int, mer_hash *hash,
			int *tpl_id, int *locus, uint64_t *value);
	void kmer_aln_query(const bwa_seq_t *query, const hash_map *hm,
			const int query_is_part, GPtrArray *hits);
	GPtrArray *align_full_seq(const bwa_seq_t *query, const hash_map *hm, const int mismatch);

#ifdef __cplusplus
}
#endif

#endif /* KMERS_H_ */
