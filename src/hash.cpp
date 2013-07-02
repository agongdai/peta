/*
 * hash.c
 *
 *  Created on: 09-May-2013
 *      Author: carl
 */
#include <vector>
#include <unordered_map>
#include <stdint.h>
#include <glib.h>
#include <ctime>
#include <stdlib.h>
#include <stdio.h>
#define __STDC_FORMAT_MACROS
#include <inttypes.h>
#include "hash.hpp"
#include "utils.h"
#include "rnaseq.h"
#include "tpl.h"
#include "kmers.hpp"
#include "peseq.h"

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
 * The elements in pos array are 64bits, where upper 48 bits representing the rna-seq id,
 * 	and lower 16 bits representing the starting position of hashing sub-sequence.
 */
hash_value get_hash_value(const index64 seq_id, const int pos_start) {
	hash_value value = seq_id;
	value <<= N_POS_BITS;
	value = value | pos_start;
	// printf("[get_hash_value] (%" ID64 ", %d) => %" ID64 "\n", seq_id, pos_start, value);
	return value;
}

void read_hash_value(index64 *seq_id, int *pos_start, hash_value value) {
	*seq_id = value >> N_POS_BITS;
	*seq_id = *seq_id & HASH_VALUE_HIGHER;
	*pos_start = value & HASH_VALUE_LOWER;
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

void build_tpl_hash(mer_hash& hash, GPtrArray *tpls, const int k) {
	uint32_t i = 0, j = 0, start = 0, end = 0;
	uint64_t kmer = 0, *kmer_freq = NULL;
	tpl *eg = NULL;
	mer_counter counter;
	for (i = 0; i < tpls->len; i++) {
		eg = (tpl*)g_ptr_array_index(tpls, i);
		if (eg->len < k)
			continue;
		end = eg->len > 2 * k ? k : (eg->len - k);
		for (j = 0; j < end; j++) {
			kmer = get_kmer_int(eg->ctg->seq, j, 1, k);
			counter[kmer]++;
		}
		start = eg->len > 2 * k ? (eg->len - 2 * k) : 0;
		for (j = start; j < eg->len - k; j++) {
			kmer = get_kmer_int(eg->ctg->seq, j, 1, k);
			counter[kmer]++;
		}
	}
	for (mer_counter::iterator it = counter.begin(); it != counter.end(); ++it) {
		kmer_freq = (uint64_t*) calloc(it->second + 1, sizeof(uint64_t));
		hash[it->first] = kmer_freq;
	}
	for (i = 0; i < tpls->len; i++) {
		eg = (tpl*)g_ptr_array_index(tpls, i);
		if (eg->len < k)
			continue;
		end = eg->len > 2 * k ? k : (eg->len - k);
		for (j = 0; j < end; j++) {
			kmer = get_kmer_int(eg->ctg->seq, j, 1, k);
			kmer_freq = hash[kmer];
			kmer_freq[++kmer_freq[0]] = get_hash_value(eg->id, j);
		}
		start = eg->len > 2 * k ? (eg->len - 2 * k) : 0;
		for (j = start; j < eg->len - k; j++) {
			kmer = get_kmer_int(eg->ctg->seq, j, 1, k);
			kmer_freq = hash[kmer];
			kmer_freq[++kmer_freq[0]] = get_hash_value(eg->id, j);
		}
	}
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

void export_frequency(const char *kmer_fa, const char *contigs_fn, const int k) {
	uint32_t n_ctgs = 0;
	uint64_t kmer = 0;
	uint64_t count = 0, i = 0, j = 0, x = 0;
	mer_hash map;
	bwa_seq_t *contigs = load_reads(contigs_fn, &n_ctgs);
	bwa_seq_t *r = NULL, *h = NULL;
	hash_map *hm = load_hash_map(kmer_fa, 1, map);
	GPtrArray *hits = NULL;
	char entry[BUFSIZE];
	sprintf(entry, "%s.bedgraph", kmer_fa);
	FILE *freq = xopen(entry, "w");
	for (i = 0; i < n_ctgs; i++) {
		r = &contigs[i];
		show_msg(__func__, "Processing %s...\n", r->name);
		for (j = 0; j <= r->len - k; j++) {
			kmer = get_kmer_int(r->seq, j, 1, k);
			count = get_kmer_rf_count(kmer, hm, 0);

			bwa_seq_t *debug = get_kmer_seq(kmer, k);
			//p_query(__func__, debug);
			hits = kmer_find_reads(debug, hm, 0, 0);
			printf("%s\t%d\t%d\t%" ID64 "\t", r->name, j, j + k, count);
			for (x = 0; x < hits->len; x++) {
				h = (bwa_seq_t*) g_ptr_array_index(hits, x);
				printf("%s,", h->name);
			}
			printf("\n");
			bwa_free_read_seq(1, debug);
			g_ptr_array_free(hits, TRUE);

			sprintf(entry, "%s\t%d\t%d\t%" ID64 "\n", r->name, j, j + k, count);
			fputs(entry, freq);
		}
	}
	fclose(freq);
}

int build_kmer_hash(int argc, char *argv[]) {
	build_kmers_hash(argv[optind], atoi(argv[optind + 1]), 1);
	return 0;
}

int export_frequency(int argc, char *argv[]) {
	export_frequency(argv[optind], argv[optind + 1], atoi(argv[optind + 2]));
	return 0;
}
