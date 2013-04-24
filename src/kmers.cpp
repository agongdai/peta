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
#include "pehash.h"
#include "peseq.h"
#include "rnaseq.h"
#include "bwtaln.h"
#include "pechar.h"

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

map_opt *new_map_opt() {
	map_opt *o = (map_opt*) malloc(sizeof(map_opt));
	o->k = 0;
	o->n_k_mers = 0;
	o->n_pos = 0;
	o->n_reads = 0;
	o->n_valid_k_mers = 0;
	o->read_len = 0;
	return o;
}

void destroy_hm(hash_map *hm) {
	if (hm) {
		bwa_free_read_seq(hm->n_reads, hm->seqs);
		free(hm->kmers_ordered);
		free(hm->o);
		hm->hash->clear();
		free(hm);
	}
}

uint64_t get_kmer_int(const ubyte_t *seq, const int start,
		const int interleaving, const int len) {
	uint64_t key = 0;
	int i = 0;
	for (i = 0; i < len * interleaving; i += interleaving) {
		if (seq[i + start] < 0 || seq[i + start] > 4) {
			return 0;
		}
		key <<= 2;
		// printf("%d \n", seq[i + start]);
		// Keep the last two digits and perform 'or' oepration
		key = key | (3 & seq[i + start]);
		// printf("%d: key = %" ID64 "\n", i, key);
	}
	return key;
}

/**
 * Structure of hash_map:
 * 		mer_opt
 * 		kmer occurrences: kmer int; # of kmers; [read int]
 * 			a 'read int' is a 'hash_value': upper 48 bits store read id; lower 16 bits store the kmer position on read
 *      list of kmer ints decreasing by '# of kmers'.
 */
void build_kmers_hash(const char *fa_fn, const int k, const int with_reads) {
	clock_t t = clock();
	uint32_t n_reads = 0, i = 0;
	int j = 0;
	uint64_t mer_v = 0, *kmer_freq = NULL;
	mer_counter counter;
	mer_hash hash;
	kmer_counter *mer = NULL;
	bwa_seq_t *reads = NULL, *r = NULL;
	FILE *hash_fp;
	char *hash_fn = (char*) malloc(BUFSIZE);
	map_opt *opt = (map_opt*) malloc(sizeof(map_opt));
	opt->k = k;

	show_msg(__func__, "Hashing library %s ...\n", fa_fn);
	reads = load_reads(fa_fn, &n_reads);
	if (n_reads <= 0)
		err_fatal(__func__, "No reads are loaded. Terminated.\n");
	opt->n_reads = n_reads;
	opt->read_len = reads->len;

	// Count frequencies of every kmer
	show_msg(__func__, "Counting kmer frequencies ...\n");
	for (i = 0; i < n_reads; i++) {
		r = &reads[i];
		for (j = 0; j <= r->len - k; j++) {
			mer_v = get_kmer_int(r->seq, j, 1, k);
			if (counter[mer_v]) {
				counter[mer_v]++;
			} else {
				counter[mer_v] = 1;
				opt->n_k_mers++;
			}
			opt->n_pos++;
		}
	}

	// Set the value of hash map: first value stores how many reads having this kmer

	show_msg(__func__, "Allocating space for %" ID64 " kmers, %" ID64 " distinct kmers ...\n", opt->n_pos, opt->n_k_mers);
	if (with_reads) {
		for (mer_counter::iterator it = counter.begin(); it != counter.end(); ++it) {
			kmer_freq = (uint64_t*) calloc(it->second + 1, sizeof(uint64_t));
			kmer_freq[0] = 0;
			hash[it->first] = kmer_freq;
			mer = (kmer_counter*) malloc(sizeof(kmer_counter));
			mer->kmer = it->first;
			mer->count = it->second;
			if (mer->count > 1) {
				opt->n_valid_k_mers++;
			}
		}

		show_msg(__func__, "Setting hash map ...\n");
		for (i = 0; i < n_reads; i++) {
			r = &reads[i];
			for (j = 0; j <= r->len - opt->k; j++) {
				mer_v = get_kmer_int(r->seq, j, 1, k);
				kmer_freq = hash[mer_v];
				kmer_freq[++kmer_freq[0]] = get_hash_value(atoi(r->name), j);
			}
		}
		sprintf(hash_fn, "%s.map.reads", fa_fn);
	} else {
		for (mer_counter::iterator it = counter.begin(); it != counter.end(); ++it) {
			kmer_freq = (uint64_t*) calloc(1, sizeof(uint64_t));
			kmer_freq[0] = it->second;
			hash[it->first] = kmer_freq;
		}
		sprintf(hash_fn, "%s.map", fa_fn);
	}

	show_msg(__func__, "Saving hash map of %d reads: %.2f sec ... \n", n_reads,
			(float) (clock() - t) / CLOCKS_PER_SEC);
	hash_fp = xopen(hash_fn, "w");
	fwrite(opt, sizeof(map_opt), 1, hash_fp);
	for (mer_hash::iterator m = hash.begin(); m != hash.end(); ++m) {
		mer_v = m->first;
		fwrite(&mer_v, sizeof(uint64_t), 1, hash_fp);
		kmer_freq = m->second;
		if (with_reads)
			fwrite(kmer_freq, sizeof(uint64_t), kmer_freq[0] + 1, hash_fp);
		else
			fwrite(kmer_freq, sizeof(uint64_t), 1, hash_fp);
	}
	show_msg(__func__, "%d reads hashed: %.2f sec\n", n_reads, (float) (clock()
			- t) / CLOCKS_PER_SEC);
}

hash_map *load_hash_map(const char *fa_fn, const int with_reads,
		mer_hash& kmers) {
	clock_t t = clock();
	FILE *map_fp = NULL;
	map_opt *opt = NULL;
	uint64_t i = 0, *count = NULL, *freq = NULL, *kmer_int = NULL;
	uint32_t n_reads = 0;
	int rs = 0;
	bwa_seq_t *seqs = NULL;
	char *map_fn = (char*) malloc(BUFSIZE);
	hash_map *hm = (hash_map*) malloc(sizeof(hash_map));

	seqs = load_reads(fa_fn, &n_reads);
	hm->seqs = seqs;
	hm->n_reads = n_reads;

	sprintf(map_fn, "%s.map.reads", fa_fn);
	show_msg(__func__, "Loading hash map from %s ...\n", map_fn);
	map_fp = xopen(map_fn, "rb");
	opt = new_map_opt();
	if (!fread(opt, sizeof(map_opt), 1, map_fp)) {
		err_fatal(__func__, "Unable to read from the map file %s! \n", map_fn);
	}
	for (i = 0; i < opt->n_k_mers; i++) {
		kmer_int = (uint64_t*) malloc(sizeof(uint64_t));
		rs = fread(kmer_int, sizeof(uint64_t), 1, map_fp);
		count = (uint64_t*) malloc(sizeof(uint64_t));
		rs = fread(count, sizeof(uint64_t), 1, map_fp);
		freq = (uint64_t*) calloc(*count + 1, sizeof(uint64_t));
		rs = fread(freq + 1, sizeof(uint64_t), *count, map_fp);
		freq[0] = *count;
		free(count);
		kmers[*kmer_int] = freq;
	}

	fclose(map_fp);
	free(map_fn);
	opt->n_reads = n_reads;
	hm->o = opt;
	hm->hash = &kmers;
	show_msg(__func__, "%" ID64 " distinct kmers loaded\n", opt->n_k_mers);
	show_msg(__func__, "%" ID64 " kmer occurrences loaded\n", opt->n_pos);
	show_msg(__func__, "Hash map loaded: %.2f sec ... \n",
			(float) (clock() - t) / CLOCKS_PER_SEC);
	return hm;
}

void test_kmer_hash(const char *fa_fn) {
	mer_hash map, *map_copy = NULL;
	bwa_seq_t *query = NULL;
	hash_map *hm = load_hash_map(fa_fn, 0, map);
	map_copy = hm->hash;

	show_msg(__func__, "Count of 3922922232675962: %d \n", get_kmer_count(
			3922922232675962, hm, 1));
	query = get_kmer_seq(3922922232675962, 25);
	p_query(__func__, query);
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

void parse_hit_ints(const uint64_t *occs, GPtrArray *hits, bwa_seq_t *seqs,
		const uint8_t ori) {
	int pos = 0;
	uint64_t j = 0, n_hits = 0, hit_id = 0, read_id = 0;
	bwa_seq_t *r = NULL;

	if (occs) {
		n_hits = occs[0];
		n_hits &= LOWER_ONES_32;
		for (j = 0; j < n_hits; j++) {
			hit_id = occs[j + 1];
			read_hash_value(&read_id, &pos, hit_id);
			r = &seqs[read_id];
			if (r->pos == -1) {
				r->pos = pos;
				r->rev_com = ori;
				g_ptr_array_add(hits, r);
			}
		}
	}
}

/**
 * Set the upper 48~64 bits to be template id, 32~48 bits to be locus
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
		count <<= 32;
		freq[0] += count;
	}
	rev_kmer_int = rev_comp_kmer(kmer_int, hm->o->k);
	it = hash->find(rev_kmer_int);
	if (it != hash->end()) {
		freq = it->second;
		count = tpl_id;
		count <<= 16;
		count += locus;
		count <<= 32;
		freq[0] += count;
	}
}

int kmer_is_used(const uint64_t kmer_int, hash_map *hm) {
	uint64_t *freq = NULL, count = 0;
	mer_hash *hash = hm->hash;
	mer_hash::iterator it = hash->find(kmer_int);
	if (it != hash->end()) {
		freq = it->second;
		count = freq[0];
		count >>= 32;
		if (count > 0)
			return 1;
		return 0;
	}
	return -1;
}

void read_tpl_using_kmer(const uint64_t kmer_int, const hash_map *hm,
		int *tpl_id, int *locus, uint64_t *value) {
	uint64_t *freq = NULL, count = 0, count_copy = 0;
	mer_hash *hash = hm->hash;
	mer_hash::iterator it = hash->find(kmer_int);
	if (it != hash->end()) {
		freq = it->second;
		count = freq[0];
		count_copy = count;
		count_copy &= LOWER_ONES_32;
		*value = count_copy;
		count >>= 32;
		count_copy = count;
		count >>= 16;
		*tpl_id = count;
		count_copy &= MIDDLE_ONES_16;
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
			count_copy >>= 32;
			if (count_copy > 0)
				return 0;
		}
		count &= LOWER_ONES_32;
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
 * Get reads containing some kmer, including forward and reverse
 */
void kmer_aln_query(const bwa_seq_t *query, const hash_map *hm, GPtrArray *hits) {
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
			parse_hit_ints(occs, hits, hm->seqs, 0);
		}
		kmer_int = get_kmer_int(query->rseq, i, 1, opt->k);
		//show_debug_msg(__func__, "Kmer int: %" ID64 "\n", kmer_int);
		it = hash->find(kmer_int);
		if (it != hash->end()) {
			occs = it->second;
			//show_debug_msg(__func__, "Kmer count: %" ID64 "\n", occs[0]);
			parse_hit_ints(occs, hits, hm->seqs, 1);
		}
	}
}

GPtrArray *kmer_find_reads(const bwa_seq_t *query, const hash_map *hm,
		const int mismatch) {
	GPtrArray *hits = NULL;
	bwa_seq_t *read = NULL, *part = NULL;
	uint32_t i = 0;
	hits = g_ptr_array_sized_new(64);
	kmer_aln_query(query, hm, hits);
	for (i = 0; i < hits->len; i++) {
		read = (bwa_seq_t*) g_ptr_array_index(hits, i);
		part = new_seq(query, read->len, read->shift);
		if (read->rev_com)
			switch_fr(part);
		if (!seq_ol(read, part, read->len, mismatch)) {
			g_ptr_array_remove_index_fast(hits, i--);
		}
		bwa_free_read_seq(1, part);
	}
	return hits;
}

int next_char_by_kmers(hash_map *hm, uint64_t kmer_int, const int fresh_only,
		const int ori) {
	int *counters = count_next_kmers(hm, kmer_int, fresh_only, ori);
	int max = get_max_index(counters);
	free(counters);
	return max;
}
