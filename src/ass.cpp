/*
 * ass.cpp
 *
 *  Created on: 09-Apr-2013
 *      Author: carl
 */

#include <unordered_map>
#include <stdint.h>
#include <glib.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#define __STDC_FORMAT_MACROS
#include <inttypes.h>
#include "bwtaln.h"
#include "edge.h"
#include "edgelist.h"
#include "mate.h"
#include "pool.h"
#include "merge.h"
#include "readrm.h"
#include "kmers.hpp"
#include "ass.hpp"
#include "rnaseq.h"
#include "pechar.h"
#include "pelib.h"
#include "utils.h"
#include "pehash.h"
#include "peseq.h"
#include "clean.h"

using namespace std;

int kmer_ctg_id = 1;
int ins_size = 0;
int sd_ins_size = 0;
int kmer_n_threads = 0;
int kmer_len = 0;
uint32_t n_used_reads = 0;
char *kmer_out = NULL;

GMutex *kmer_id_mutex;
struct timespec kmer_start_time, kmer_finish_time;

edge *blank_edge(uint64_t query_int, int init_len, int ori) {
	bwa_seq_t *kmer = NULL;
	edge *eg = new_eg();
	g_mutex_lock(kmer_id_mutex);
	eg->id = kmer_ctg_id++;
	g_mutex_unlock(kmer_id_mutex);
	// Get a copy of the kmer
	kmer = get_kmer_seq(query_int, init_len);
	if (ori)
		eg->contig = new_seq(kmer, init_len, 0);
	else
		eg->contig = new_seq(kmer, init_len, kmer->len - init_len);
	eg->len = eg->contig->len;
	eg->tid = atoi(kmer->name);
	eg->name = strdup(kmer->name);
	bwa_free_read_seq(1, kmer);
	return eg;
}

gint cmp_kmers_by_count(gpointer a, gpointer b) {
	kmer_counter *c_a = *((kmer_counter**) a);
	kmer_counter *c_b = *((kmer_counter**) b);
	return ((c_b->count) - c_a->count);
}

/**
 * Extend the short template as int64
 * max_len: should be smaller than 64
 */
uint64_t try_short_tpl_ext(hash_map *hm, uint64_t query, const int ori,
		const int max_len) {
	uint64_t branch_template = 0;
	int i = 0, c = 0;
	for (i = 0; i < max_len; i++) {
		c = next_char_by_kmers(hm, query, ori);
		if (c == -1)
			break;
		else {
			branch_template = shift_bit(branch_template, c, hm->o->k, ori);
			query = shift_bit(query, c, hm->o->k, ori);
		}
	}
	return branch_template;
}

/**
 * Extend the short template as ubyte_t's, max_len may be larger than 64
 */
ubyte_t *try_short_tpl_byte_ext(hash_map *hm, uint64_t query, const int ori,
		const int max_len) {
	ubyte_t *branch_seq = NULL, c = 0;
	int i = 0;
	branch_seq = (ubyte_t*) calloc(max_len, sizeof(ubyte_t));
	for (i = 0; i < max_len; i++) {
		c = next_char_by_kmers(hm, query, ori);
		if (c == -1)
			break;
		else
			branch_seq[i] = c;
	}
	return branch_seq;
}

int find_junction_reads(hash_map *hm, uint64_t query_int, int c, const int ori,
		const int max_len) {
	uint64_t query = 0, i = 0, is_valid = 1;
	ubyte_t *branch_seq = NULL;
	query = shift_bit(query_int, c, hm->o->k, ori);
	branch_seq = try_short_tpl_byte_ext(hm, query, ori, max_len);
	free(branch_seq);
	return is_valid;
}

void val_short_tpl(hash_map *hm, uint64_t query_int, int max_c, int second_c,
		const int ori) {
	uint64_t main_query = 0, main_branch = 0;
	uint64_t second_query = 0, second_branch = 0;
	main_query = shift_bit(query_int, max_c, hm->o->k, ori);
	main_branch = try_short_tpl_ext(hm, main_query, ori, SHORT_BRANCH_LEN);
	second_query = shift_bit(query_int, second_c, hm->o->k, ori);
	second_branch = try_short_tpl_ext(hm, main_query, ori, SHORT_BRANCH_LEN);
	if (main_branch == second_branch)
		return 1;
	return 0;
}

void one_step_forward(edge *eg, uint64_t query, int next_c, hash_map *hm,
		const int ori) {
	mark_kmer_used(query, hm);
	ext_con(eg->contig, next_c, 0);
	eg->len = eg->contig->len;
	return shift_bit(query, next_c, hm->o->k, ori);
}

void kmer_ext_edge(edge *eg, uint64_t query_int, hash_map *hm, const int ori) {
	int max_c = 0, second_c = 0, *counters = NULL;
	int main_is_valid = 0, branch_is_valid = 0;
	uint64_t rev_kmer_int = 0, second_query = 0;
	edge *branch_eg = NULL;
	//	bwa_seq_t *debug = NULL;

	if (ori)
		seq_reverse(eg->len, eg->contig->seq, 0);
	while (1) {
		counters = count_next_kmers(hm, query_int, ori);
		max_c = get_max_index(counters);
		second_c = get_second_freq_char(counters, max_c);
		free(counters);

		/**
		 if (ori)
		 seq_reverse(eg->len, eg->contig->seq, 0);
		 debug = get_kmer_seq(query_int, 25);
		 p_query(__func__, debug);
		 bwa_free_read_seq(1, debug);
		 show_debug_msg(__func__,
		 "Ori %d, Edge %d, length %d, Next char: %d \n", ori, eg->id,
		 eg->len, c);
		 p_ctg_seq("Contig", eg->contig);
		 if (ori)
		 seq_reverse(eg->len, eg->contig->seq, 0);
		 **/
		if (max_c == -1) {
			show_debug_msg(__func__, "[%d, %d] No hits, stop here. \n", eg->id,
					eg->len);
			break;
		}
		if (second_c != -1) {
			if (val_short_tpl(hm, query_int, max_c, second_c, ori)) {
				main_is_valid = find_junction_reads(hm, query_int, max_c, ori,
						hm->o->read_len - SHORT_BRANCH_LEN);
				branch_is_valid = find_junction_reads(hm, query_int, second_c,
						ori, hm->o->read_len - SHORT_BRANCH_LEN);
				if (main_is_valid && branch_is_valid) {
					query_int = one_step_forward(eg, query_int, max_c, hm, ori);
					kmer_ext_edge(eg, query_int, hm, ori);
					branch_eg = blank_edge(query_int, 1, ori);
				}
			} else {
				// Mark the branch kmer as used
				second_query = shift_bit(query_int, second_c, hm->o->k, ori);
				mark_kmer_used(second_query, hm);
			}
		}
		query_int = one_step_forward(eg, query_int, max_c, hm, ori);
		if (eg->len % 100 == 0)
			show_debug_msg(__func__, "Ori %d, Edge %d, length %d \n", ori,
					eg->id, eg->len);
	}
	if (ori)
		seq_reverse(eg->len, eg->contig->seq, 0);
}

/**
 * Extend a kmer
 */
void *kmer_ext_thread(gpointer data, gpointer thread_params) {
	edge *eg = NULL;
	bwa_seq_t *kmer = NULL;
	int round_1_len = 0, round_2_len = 0;
	uint64_t kmer_int = 0;
	map_opt *opt = NULL;
	kmer_counter *counter = NULL;
	kmer_t_meta *params = (kmer_t_meta*) thread_params;
	opt = params->hm->o;

	counter = (kmer_counter*) data;
	kmer_int = counter->kmer;
	if (get_kmer_count(kmer_int, params->hm) <= 1) {
		return NULL;
	}
	eg = blank_edge(kmer_int, opt->k, 0);

	show_debug_msg(__func__, "============= %s: %d ============ \n",
			kmer->name, counter->count);
	kmer_ext_edge(eg, kmer_int, params->hm, 0);
	round_1_len = eg->len;
	show_debug_msg(__func__, "Edge %d with length: %d \n", eg->id, eg->len);

	kmer_int = get_kmer_int(kmer->seq, 0, 1, kmer->len);
	bwa_free_read_seq(1, kmer);
	kmer_ext_edge(eg, kmer_int, params->hm, 1);
	show_debug_msg(__func__, "Edge %d with length: %d \n", eg->id, eg->len);

	if (eg->len - round_1_len > 2) {
		kmer = new_seq(eg->contig, opt->k, eg->len - opt->k);
		kmer_int = get_kmer_int(kmer->seq, 0, 1, kmer->len);
		bwa_free_read_seq(1, kmer);
		round_2_len = eg->len;
		kmer_ext_edge(eg, kmer_int, params->hm, 0);
		show_debug_msg(__func__, "Edge %d with length: %d \n", eg->id, eg->len);
		if (eg->len - round_2_len > 2) {
			kmer = new_seq(eg->contig, opt->k, 0);
			kmer_int = get_kmer_int(kmer->seq, 0, 1, kmer->len);
			bwa_free_read_seq(1, kmer);
			kmer_ext_edge(eg, kmer_int, params->hm, 1);
			show_debug_msg(__func__, "Edge %d with length: %d \n", eg->id,
					eg->len);
		}
	}

	eg->start_kmer_int = *((uint64_t*) data);
	if (eg->len < 50) {
		destroy_eg(eg);
	} else {
		g_mutex_lock(kmer_id_mutex);
		g_ptr_array_add(params->all_edges, eg);
		g_mutex_unlock(kmer_id_mutex);
	}
	free(counter);
	return NULL;
}

void kmer_threads(kmer_t_meta *params, GPtrArray *solid_reads) {
	GThreadPool *thread_pool = NULL;
	uint64_t i = 0, *count = 0;
	hash_map *hm = params->hm;
	mer_hash *hash = hm->hash;
	kmer_counter *counter = NULL;
	GPtrArray *start_kmers = g_ptr_array_sized_new(BUFSIZ);

	thread_pool = g_thread_pool_new((GFunc) kmer_ext_thread, params, 1, TRUE,
			NULL);

	show_msg(__func__, "Extending by kmers...\n");
	for (mer_hash::iterator it = hash->begin(); it != hash->end(); ++it) {
		count = it->second;
		if (count[0] > 1) {
			counter = (kmer_counter*) malloc(sizeof(kmer_counter));
			counter->kmer = it->first;
			counter->count = count[0];
			g_ptr_array_add(start_kmers, counter);
		}
	}
	g_ptr_array_sort(start_kmers, (GCompareFunc) cmp_kmers_by_count);
	for (i = 0; i < start_kmers->len; i++) {
		counter = (kmer_counter*) g_ptr_array_index(start_kmers, i);
		g_thread_pool_push(thread_pool, (gpointer) counter, NULL);
		//if (i > 100)
		//	break;
	}
	g_thread_pool_free(thread_pool, 0, 1);
}

void pick_unused_kmers(kmer_t_meta *params) {
	GThreadPool *thread_pool = NULL;
	hash_map *hm = params->hm;
	mer_hash *hash = hm->hash;
	map_opt *opt = hm->o;
	uint64_t mer_v = 0, n = 0, i = 0;
	bwa_seq_t *seq = NULL;
	GPtrArray *remaining = g_ptr_array_sized_new(BUFSIZ);
	show_msg(__func__, "Going to assemble remaining reads: %d/%" ID64 " ...\n", n_used_reads, opt->n_reads);
	for (i = 0; i < hm->n_reads; i++) {
		seq = &hm->seqs[i];
		if (seq->status == FRESH)
			g_ptr_array_add(remaining, seq);
	}
	show_msg(__func__, "%d reads remaining \n", remaining->len);
	thread_pool = g_thread_pool_new((GFunc) kmer_ext_thread, params, 1, TRUE,
			NULL);
	for (i = 0; i < remaining->len; i++) {
		seq = (bwa_seq_t*) g_ptr_array_index(remaining, i);
		g_thread_pool_push(thread_pool, (gpointer) seq, NULL);
	}
	g_thread_pool_free(thread_pool, 0, 1);
	g_ptr_array_free(remaining, TRUE);
}

void test_kmer_ext(kmer_t_meta *params) {
	uint64_t kmer_int = 289104232493298;
	edge *eg = new_eg();
	FILE *contigs = NULL;
	GPtrArray *all_edges = g_ptr_array_sized_new(12);
	bwa_seq_t *kmer = get_kmer_seq(kmer_int, params->hm->o->k);

	eg->contig = new_seq(kmer, kmer->len, 0);
	eg->len = eg->contig->len;
	eg->tid = atoi(kmer->name);
	eg->name = strdup(kmer->name);
	kmer_ext_edge(eg, kmer_int, params->hm, 0);
	kmer_ext_edge(eg, kmer_int, params->hm, 1);
	g_ptr_array_add(all_edges, eg);

	contigs = xopen(get_output_file("single.fa", kmer_out), "w");
	save_edges(all_edges, contigs, 0, 0, 0);
}

void ext_by_kmers_core(char *lib_file, const char *solid_file) {
	hash_map *hm = NULL;
	mer_hash map;
	GPtrArray *all_edges = NULL;
	GPtrArray *solid_reads = NULL;
	FILE *contigs = NULL;
	kmer_t_meta *params = (kmer_t_meta*) calloc(1, sizeof(kmer_t_meta));

	show_msg(__func__, "Library: %s \n", lib_file);
	show_msg(__func__, "Solid Reads: %s \n", solid_file);

	hm = load_hash_map(lib_file, 0, map);
	hm->hash = &map;

	clock_gettime(CLOCK_MONOTONIC, &kmer_finish_time);
	show_msg(__func__, "Done loading kmers hash map: %.2f sec\n",
			(float) (kmer_finish_time.tv_sec - kmer_start_time.tv_sec));

	all_edges = g_ptr_array_sized_new(BUFSIZ);
	params->hm = hm;
	params->all_edges = all_edges;
	solid_reads = load_solid_reads(solid_file, hm->seqs, hm->n_reads);

	test_kmer_ext(params);
	exit(1);
	kmer_threads(params, solid_reads);
	//pick_unused_kmers(params);

	contigs = xopen(get_output_file("paired.fa", kmer_out), "w");
	save_edges(all_edges, contigs, 0, 0, 0);
	fflush(contigs);
	fclose(contigs);

	reset_edge_ids(all_edges);
	merge_ol_edges(all_edges, ins_size, sd_ins_size, hm->seqs, kmer_n_threads);
	free(params);
	destroy_hm(hm);

	contigs = xopen(get_output_file("merged.fa", kmer_out), "w");
	save_edges(all_edges, contigs, 0, 0, 100);
	fflush(contigs);
	fclose(contigs);
	contigs = xopen(get_output_file("merged_all.fa", kmer_out), "w");
	save_edges(all_edges, contigs, 0, 0, 0);
	fflush(contigs);
	fclose(contigs);
}

int pe_kmer(int argc, char *argv[]) {
	int c = 0;
	clock_gettime(CLOCK_MONOTONIC, &kmer_start_time);
	while ((c = getopt(argc, argv, "k:m:s:o:t:")) >= 0) {
		switch (c) {
		case 'k':
			kmer_len = atoi(optarg);
			break;
		case 'o':
			kmer_out = optarg;
			break;
		case 'm':
			ins_size = atoi(optarg);
			break;
		case 's':
			sd_ins_size = atoi(optarg);
			break;
		case 't':
			kmer_n_threads = atoi(optarg);
			break;
		}
	}
	if (optind + 2 > argc) {
		show_msg(__func__, "Parameters error! \n");
		return 1;
	}
	if (!g_thread_supported())
		g_thread_init(NULL);
	kmer_id_mutex = g_mutex_new();

	ext_by_kmers_core(argv[optind], argv[optind + 1]);

	clock_gettime(CLOCK_MONOTONIC, &kmer_finish_time);
	show_msg(__func__, "Done: %.2f sec\n", (float) (kmer_finish_time.tv_sec
			- kmer_start_time.tv_sec));
	return 0;
}
