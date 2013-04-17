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
void build_kmers_hash(const char *fa_fn, const int k) {
	clock_t t = clock();
	uint32_t n_reads = 0, i = 0, j = 0;
	uint64_t mer_v = 0, *kmer_freq = NULL, m = 0;
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
	opt->n_reads = n_reads;

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
	show_msg(__func__, "Saving hash map of %d reads: %.2f sec ... \n", n_reads,
			(float) (clock() - t) / CLOCKS_PER_SEC);
	sprintf(hash_fn, "%s.map", fa_fn);
	hash_fp = xopen(hash_fn, "w");
	fwrite(opt, sizeof(map_opt), 1, hash_fp);
	for (mer_hash::iterator m = hash.begin(); m != hash.end(); ++m) {
		mer_v = m->first;
		fwrite(&mer_v, sizeof(uint64_t), 1, hash_fp);
		kmer_freq = m->second;
		fwrite(kmer_freq, sizeof(uint64_t), kmer_freq[0] + 1, hash_fp);
	}
	show_msg(__func__, "%d reads hashed: %.2f sec\n", n_reads, (float) (clock()
			- t) / CLOCKS_PER_SEC);
}

hash_map *load_hash_map(const char *fa_fn, mer_hash& kmers) {
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

	sprintf(map_fn, "%s.map", fa_fn);
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
	uint64_t i = 0, j = 0;
	GPtrArray *hits = g_ptr_array_sized_new(BUFSIZ);
	bwa_seq_t *query = NULL;
	hash_map *hm = load_hash_map(fa_fn, map);
	map_copy = hm->hash;

	show_msg(__func__, "Count of 3922922232675962: %d \n", get_kmer_count(
			3922922232675962, hm));
	query = get_kmer_seq(3922922232675962, 25);
	p_query(__func__, query);
}

bwa_seq_t *get_kmer_seq(uint64_t kmer, const int k) {
	ubyte_t *seq = NULL;
	bwa_seq_t *read = NULL;
	int i = 0;
	uint64_t copy = kmer, index = 0;
	read = blank_seq();
	free(read->seq);
	read->seq = (ubyte_t*) calloc(k + 1, sizeof(ubyte_t));
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
 * Free the memory allocated to a kmer;
 * Remove the kmer from the hashmap
 */
void mark_kmer_used(const uint64_t kmer_int, const hash_map *hm) {
	uint64_t *freq = NULL;
	mer_hash *hash = hm->hash;
	mer_hash::iterator it = hash->find(kmer_int);
	if (it != hash->end()) {
		freq = it->second;
		freq[0] = 0;
	}
}

/**
 * Get how many kmers in the hash map
 */
uint64_t get_kmer_count(const uint64_t kmer_int, hash_map *hm) {
	uint64_t *freq = NULL;
	mer_hash *hash = hm->hash;
	mer_hash::iterator it = hash->find(kmer_int);
	if (it != hash->end()) {
		freq = it->second;
		return freq[0];
	}
	return 0;
}

/**
 * Get reads containing some kmer, including forward and reverse
 */
void kmer_aln_query(const bwa_seq_t *query, const hash_map *hm, GPtrArray *hits) {
	GPtrArray *ret_hits = NULL;
	uint64_t kmer_int = 0, i = 0, *occs = NULL;
	map_opt *opt = hm->o;
	mer_hash *hash = hm->hash;
	bwa_seq_t *r = NULL, *r_pre = NULL;
	mer_hash::iterator it;

	ret_hits = g_ptr_array_sized_new(0);
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

/**
 * Check kmer frequencies of all possible next chars, get the highest one
 * Return: ACGTN -> 01234;
 *         None -> -1;
 */
int next_char_by_kmers(mer_hash *kmers, const int k, bwa_seq_t *query,
		const int ori) {
	bwa_seq_t *next_q = NULL;
	int counters[5], i = 0, max = 0;
	uint64_t kmer_int = 0, *freq = NULL;
	for (i = 0; i < 5; i++) {
		next_q = new_seq(query, query->len, 0);
		ext_que(next_q, i, ori);
		counters[i] = 0;
		// Check the forward kmer
		//p_query("SEQ", next_q);
		kmer_int = get_kmer_int(next_q->seq, 0, 1, k);
		freq = (*kmers)[kmer_int];
		if (freq) {
			counters[i] = freq[0];
		}
		// Check the reverse kmer
		kmer_int = get_kmer_int(next_q->rseq, 0, 1, k);
		freq = (*kmers)[kmer_int];
		if (freq) {
			counters[i] += freq[0];
		}
		bwa_free_read_seq(1, next_q);
	}
	// Get max occurrence of the chars
	for (i = 0; i < 5; i++) {
		max = (counters[i] > max) ? counters[i] : max;
	}
	//show_debug_msg(__func__, "Max: %d \n", max);
	if (max == 0)
		return -1;
	for (i = 0; i < 5; i++) {
		if (max == counters[i])
			return i;
	}
	return -1;
}
