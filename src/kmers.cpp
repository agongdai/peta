#include <vector>
#include <unordered_map>
#include <cassert>
#include <unistd.h>
#include <fstream>
#include <string>
#include <ctime>
#include <algorithm>
#include <stdint.h>
#include <glib.h>
#include <time.h>
#define __STDC_FORMAT_MACROS
#include <inttypes.h>
#include <string.h>
#include <stdio.h>
#include "kmers.hpp"
#include "utils.h"
#include "peseq.h"
#include "rnaseq.h"
#include "bwtaln.h"
#include "pechar.h"
#include "hash.hpp"

/**
 * Return a reverse complement value of a kmer int
 */
uint64_t rev_comp_kmer(uint64_t kmer, const int n) {
	int i = 0, c = 0;
	uint64_t copy = kmer;
	uint64_t rev_comp = 0;
	for (i = 1; i <= n; i++) {
		c = 3 & copy;
		c = 3 - c; // complement of next char
		rev_comp <<= 2;
		rev_comp |= c;
		copy >>= 2;
	}
	return rev_comp;
}

bwa_seq_t *get_kmer_seq(uint64_t kmer, const int k) {
	ubyte_t *seq = NULL;
	bwa_seq_t *read = NULL;
	int i = 0;
	uint64_t copy = kmer, index = 0;
	read = blank_seq(k);
	seq = read->seq;
	for (i = 0; i < k; i++) {
		index = copy;
		index &= 3; // Keep only the last two bits
		seq[k - 1 - i] = index;
		copy >>= 2;
	}
	read->len = k;
	read->status = FRESH;
	read->name = (char*) malloc(64);
	sprintf(read->name, "%" ID64, kmer);
	set_rev_com(read);
	return read;
}

/**
 * Set the upper 40~64 bits to be template id, 24~40 bits to be locus
 */
void mark_kmer_used(const uint64_t kmer_int, const hash_map *hm,
		const int tpl_id, const int locus) {
	uint64_t *freq = NULL, rev_kmer_int = 0, count = 0;
	mer_hash *hash = hm->hash;
	mer_hash::iterator it = hash->find(kmer_int);
	if (it != hash->end()) {
		freq = it->second;
		count = tpl_id;
		count <<= 16;
		count += locus;
		count <<= 24;
		// Reset the used template to be none
		freq[0] <<= 40;
		freq[0] >>= 40;
		freq[0] += count;
	}
	rev_kmer_int = rev_comp_kmer(kmer_int, hm->o->k);
	it = hash->find(rev_kmer_int);
	if (it != hash->end()) {
		freq = it->second;
		count = tpl_id;
		count <<= 16;
		count += locus;
		count <<= 24;
		// Reset the used template to be none
		freq[0] <<= 40;
		freq[0] >>= 40;
		freq[0] += count;
	}
}

/**
 * Mark a kmer as not used.
 */
void mark_kmer_not_used(const uint64_t kmer_int, const hash_map *hm) {
	uint64_t *freq = NULL, rev_kmer_int = 0, count = 0;
	mer_hash *hash = hm->hash;
	mer_hash::iterator it = hash->find(kmer_int);
	if (it != hash->end()) {
		freq = it->second;
		freq[0] <<= 40;
		freq[0] >>= 40;
	}
	rev_kmer_int = rev_comp_kmer(kmer_int, hm->o->k);
	it = hash->find(rev_kmer_int);
	if (it != hash->end()) {
		freq = it->second;
		freq[0] <<= 40;
		freq[0] >>= 40;
	}
}

/**
 * Check whether a kmer is used.
 */
int kmer_is_used(const uint64_t kmer_int, hash_map *hm) {
	uint64_t *freq = NULL, count = 0;
	mer_hash *hash = hm->hash;
	mer_hash::iterator it = hash->find(kmer_int);
	if (it != hash->end()) {
		freq = it->second;
		count = freq[0];
		count >>= 24;
		if (count > 0)
			return 1;
		return 0;
	}
	return -1;
}

/**
 * Read the teamplate id and locus using some kmer
 */
void read_tpl_using_kmer(const uint64_t kmer_int, const hash_map *hm,
		int *tpl_id, int *locus, uint64_t *value) {
	uint64_t *freq = NULL, count = 0, count_copy = 0;
	mer_hash *hash = hm->hash;
	mer_hash::iterator it = hash->find(kmer_int);
	if (it != hash->end()) {
		freq = it->second;
		count = freq[0];
		count_copy = count;
		count_copy &= LOWER_ONES_24;
		*value = count_copy;
		count >>= 24;
		count_copy = count;
		count >>= 16;
		*tpl_id = count;
		count_copy &= LOWER_ONES_16;
		*locus = count_copy;
	}
}

/**
 * Get how many kmers in the hash map
 * if fresh_only is 1, return the frequency only if the kmer is not used.
 */
uint64_t get_kmer_count(const uint64_t kmer_int, hash_map *hm,
		const int fresh_only) {
	uint64_t *freq = NULL, count = 0, count_copy = 0;
	mer_hash *hash = hm->hash;
	mer_hash::iterator it = hash->find(kmer_int);
	if (it != hash->end()) {
		freq = it->second;
		count = freq[0];
		count_copy = freq[0];
		// The upper 32 bits stores the template using this kmer,
		//	reset the upper 32 bits to be 0, to get the frequency
		if (fresh_only) {
			count_copy >>= 24;
			if (count_copy > 0)
				return 0;
		}
		count <<= 40;
		count >>= 40;
		//show_debug_msg(__func__, "Kmer: %" ID64 "=>%" ID64 "\n", kmer_int, count);
		return count;
	}
	return 0;
}

int *count_next_kmers(hash_map *hm, uint64_t query, const int fresh_only,
		const int ori) {
	int *counters = (int*) calloc(4, sizeof(int)), i = 0;
	uint64_t next_probable_kmer = 0;
	for (i = 0; i < 4; i++) {
		counters[i] = 0;
		next_probable_kmer = shift_bit(query, i, hm->o->k, ori);
		counters[i] = get_kmer_count(next_probable_kmer, hm, fresh_only);
		// Get its reverse complement as well
		next_probable_kmer = rev_comp_kmer(next_probable_kmer, hm->o->k);
		counters[i] += get_kmer_count(next_probable_kmer, hm, fresh_only);
	}
	return counters;
}

/**
 *  query_is_part=1: the query length is shorter than read length
 *  query_is_part=0: the returned reads should be on the query
 */
void parse_hit_ints(const uint64_t *occs, const int query_i,
		const int query_is_part, const int query_len, GPtrArray *hits,
		bwa_seq_t *seqs, const uint8_t ori) {
	int pos = 0;
	uint64_t j = 0, n_hits = 0, hit_id = 0, read_id = 0;
	bwa_seq_t *r = NULL;

	if (occs) {
		n_hits = occs[0];
		n_hits &= LOWER_ONES_24;
		for (j = 0; j < n_hits; j++) {
			hit_id = occs[j + 1];
			read_hash_value(&read_id, &pos, hit_id);
			r = &seqs[read_id];
			// Make sure the read can be put into the query
			// Query: --------|-----------|------
			//            i        kmer
			// Read:     -----|-----------|---
			//            pos
			//p_query(__func__, r);
			//show_debug_msg(__func__, "query_i: %d; pos: %d; r->len: %d \n",
			//		query_i, pos, r->len);
			if (r->pos == -1) { // To avoid adding the same read repeatly.
				if (query_is_part || (query_i - pos >= 0 && (query_i - pos
						+ r->len) <= query_len)) {
					r->pos = query_i - pos;
					r->rev_com = ori;
					g_ptr_array_add(hits, r);
				}
			}
		}
	}
}

/**
 * Find common reads of two lists of reads
 */
GPtrArray *interset_reads(GPtrArray *list_1, GPtrArray *list_2, GPtrArray *set) {
	uint32_t index_1 = 0, index_2 = 0;
	bwa_seq_t *read_1 = NULL, *read_2 = NULL;
	g_ptr_array_sort(list_1, (GCompareFunc) cmp_reads_by_name);
	g_ptr_array_sort(list_2, (GCompareFunc) cmp_reads_by_name);
	if (!set)
		set = g_ptr_array_sized_new(list_1->len);
	while (index_1 < list_1->len && index_2 < list_2->len) {
		read_1 = (bwa_seq_t*) g_ptr_array_index(list_1, index_1);
		//p_query("READ 1", read_1);
		while (index_2 < list_2->len) {
			read_2 = (bwa_seq_t*) g_ptr_array_index(list_2, index_2);
			//p_query("READ 2", read_2);
			if (read_1 == read_2) {
				g_ptr_array_add(set, read_1);
				//show_debug_msg(__func__, "INTERSET\n");
				index_2++;
				index_1++;
				break;
			}
			if (atoi(read_1->name) < atoi(read_2->name)) {
				read_1->pos = -1;
				index_1++;
				break;
			} else {
				read_1->pos = -1;
				index_2++;
			}
		}
	}
	return set;
}

/**
 * Get reads containing some kmer, including forward and reverse
 */
void kmer_aln_query(const bwa_seq_t *query, const hash_map *hm,
		const int query_is_part, GPtrArray *hits) {
	uint64_t kmer_int = 0, *occs = NULL;
	int i = 0;
	map_opt *opt = hm->o;
	mer_hash *hash = hm->hash;
	mer_hash::iterator it;

	for (i = 0; i <= query->len - opt->k; i++) {
		kmer_int = get_kmer_int(query->seq, i, 1, opt->k);
		//show_debug_msg(__func__, "Kmer int: %" ID64 "\n", kmer_int);
		it = hash->find(kmer_int);
		if (it != hash->end()) {
			occs = it->second;
			//show_debug_msg(__func__, "Kmer count: %" ID64 "\n", occs[0]);
			parse_hit_ints(occs, i, query_is_part, query->len, hits, hm->seqs,
					0);
		}
		kmer_int = get_kmer_int(query->rseq, i, 1, opt->k);
		//show_debug_msg(__func__, "Kmer int: %" ID64 "\n", kmer_int);
		it = hash->find(kmer_int);
		if (it != hash->end()) {
			occs = it->second;
			//show_debug_msg(__func__, "Kmer count: %" ID64 "\n", occs[0]);
			parse_hit_ints(occs, i, query_is_part, query->len, hits, hm->seqs,
					1);
		}
	}
}

/**
 * Get reads with both head and tail kmers, query length is the same as read length
 */
GPtrArray *head_tail_kmer_reads(const bwa_seq_t *query, const hash_map *hm,
		GPtrArray *hits) {
	bwa_seq_t *kmer_seq = NULL, *read = NULL;
	uint32_t i = 0;
	if (query->len <= hm->o->k)
		return NULL;
	GPtrArray *head_hits = g_ptr_array_sized_new(0);
	GPtrArray *tail_hits = g_ptr_array_sized_new(0);
	if (!hits)
		hits = g_ptr_array_sized_new(32);

	kmer_seq = new_seq(query, hm->o->k, 0);
	kmer_aln_query(kmer_seq, hm, 1, head_hits);
	bwa_free_read_seq(1, kmer_seq);
	for (i = 0; i < head_hits->len; i++) {
		read = (bwa_seq_t*) g_ptr_array_index(head_hits, i);
		read->pos = -1;
	}

	kmer_seq = new_seq(query, hm->o->k, query->len - hm->o->k);
	kmer_aln_query(kmer_seq, hm, 1, tail_hits);
	for (i = 0; i < tail_hits->len; i++) {
		read = (bwa_seq_t*) g_ptr_array_index(tail_hits, i);
		read->pos = -1;
	}
	bwa_free_read_seq(1, kmer_seq);

	interset_reads(head_hits, tail_hits, hits);
	g_ptr_array_free(head_hits, TRUE);
	g_ptr_array_free(tail_hits, TRUE);
	return hits;
}

/**
 * Align all reads to the 'seq' with 'mismatch'
 * The length of 'seq' should be longer than read length
 */
GPtrArray *align_full_seq(const bwa_seq_t *query, const hash_map *hm,
		const int mismatch) {
	GPtrArray *hits = g_ptr_array_sized_new(32);
	uint32_t i = 0, j = 0;
	bwa_seq_t *read = NULL, *part = NULL;
	if (query->len < hm->o->read_len)
		return g_ptr_array_sized_new(0);
	for (i = 0; i <= query->len - hm->o->read_len; i++) {
		part = new_seq(query, hm->o->read_len, i);
		head_tail_kmer_reads(query, hm, hits);
		//p_readarray(hits, 1);
		for (j = 0; j < hits->len; j++) {
			read = (bwa_seq_t*) g_ptr_array_index(hits, j);
			if (read->rev_com)
				switch_fr(part);
			if (!seq_ol(read, part, read->len, mismatch)) {
				read->pos = -1;
				g_ptr_array_remove_index_fast(hits, j--);
			}
			if (read->rev_com)
				switch_fr(part);
		}
		bwa_free_read_seq(1, part);
	}
	return hits;
}

/**
 * Get reads that can be aligned to parameter 'query'.
 * If n_part_only > 0, return max n_part_only hits.
 */
GPtrArray *kmer_find_reads(const bwa_seq_t *query, const hash_map *hm,
		const int mismatch, const uint32_t n_part_only) {
	GPtrArray *hits = NULL, *part_hits = NULL;
	bwa_seq_t *read = NULL, *part = NULL, *r = NULL;
	uint32_t i = 0;
	hits = g_ptr_array_sized_new(64);
	if (n_part_only)
		part_hits = g_ptr_array_sized_new(n_part_only);
	kmer_aln_query(query, hm, 0, hits);
	for (i = 0; i < hits->len; i++) {
		read = (bwa_seq_t*) g_ptr_array_index(hits, i);
		part = new_seq(query, read->len, read->pos);
		//show_debug_msg(__func__, "pos: %d \n", read->pos);
		if (read->rev_com)
			switch_fr(part);
		//p_query(__func__, read);
		//p_query(__func__, part);
		//show_debug_msg(__func__, "---\n");
		if (!seq_ol(read, part, read->len, mismatch)) {
			read->pos = -1;
			g_ptr_array_remove_index_fast(hits, i--);
		} else {
			if (n_part_only) {
				g_ptr_array_add(part_hits, read);
				if (part_hits->len >= n_part_only) {
					bwa_free_read_seq(1, part);
					break;
				}
			}
		}
		bwa_free_read_seq(1, part);
	}
	// Reset the hit positions
	for (i = 0; i < hits->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(hits, i);
		r->pos = -1;
	}
	if (n_part_only) {
		g_ptr_array_free(hits, TRUE);
		return part_hits;
	}
	return hits;
}

int next_char_by_kmers(hash_map *hm, uint64_t kmer_int, const int fresh_only,
		const int ori) {
	int *counters = count_next_kmers(hm, kmer_int, fresh_only, ori);
	//show_debug_msg(__func__, "Counters: [%d, %d, %d, %d]\n", counters[0],
	//		counters[1], counters[2], counters[3]);
	int max = get_max_index(counters);
	free(counters);
	return max;
}

int next_char_max_freq(hash_map *hm, uint64_t kmer_int, const int fresh_only,
		const int ori) {
	int *counters = count_next_kmers(hm, kmer_int, fresh_only, ori);
	//show_debug_msg(__func__, "Counters: [%d, %d, %d, %d]\n", counters[0],
	//		counters[1], counters[2], counters[3]);
	int max_index = get_max_index(counters);
	int max = counters[max_index];
	free(counters);
	return max;
}

int find_kmer(int argc, char *argv[]) {
	int c;
	while ((c = getopt(argc, argv, "")) >= 0) {
	}
	if (optind + 1 > argc) {
		return 1;
	}
	return 0;
}
