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

int kmer_ctg_id = 0;
int ins_size = 0;
int sd_ins_size = 0;
int kmer_n_threads = 0;
int kmer_len = 0;
uint32_t n_used_reads = 0;
char *kmer_out = NULL;

GMutex *kmer_id_mutex;
struct timespec kmer_start_time, kmer_finish_time;

gint cmp_kmers_by_count(gpointer a, gpointer b) {
	kmer_counter *c_a = *((kmer_counter**) a);
	kmer_counter *c_b = *((kmer_counter**) b);
	return ((c_b->count) - c_a->count);
}

int next_char_by_kmers(hash_map *hm, uint64_t kmer_int, const int ori) {
	int counters[4], i = 0, max = 0;
	uint64_t next_kmer_int = 0;
	for (i = 0; i < 4; i++) {
		counters[i] = 0;
		// Check the forward kmer
		next_kmer_int = shift_bit(kmer_int, i, hm->o->k, ori);
		counters[i] = get_kmer_count(next_kmer_int, hm);
		// Check the reverse kmer
		next_kmer_int = rev_comp_kmer(next_kmer_int, hm->o->k);
		counters[i] += get_kmer_count(next_kmer_int, hm);
	}
	// Get max occurrence of the chars
	for (i = 0; i < 4; i++) {
		max = (counters[i] > max) ? counters[i] : max;
	}
	show_debug_msg(__func__, "Next [%d:%d:%d:%d]\n", counters[0], counters[1],
			counters[2], counters[3]);
	show_debug_msg(__func__, "Max: %d \n", max);
	if (max == 0)
		return -1;
	for (i = 0; i < 4; i++) {
		if (max == counters[i])
			return i;
	}
	return -1;
}

void kmer_pool(GPtrArray *hits, const hash_map *hm, pool *cur_pool, edge *eg,
		bwa_seq_t *query, int *next, const int ori) {
	uint32_t i = 0;
	bwa_seq_t *s = NULL, *mate = NULL, *seqs = hm->seqs;
	// Add aligned reads to current pool
	for (i = 0; i < hits->len; i++) {
		s = (bwa_seq_t*) g_ptr_array_index(hits, i);
		mate = get_mate(s, seqs);
		//p_query(__func__, s);
		//p_query(__func__, mate);
		//pre_cursor = s->cursor;

		if (s->status != FRESH || s->tid != -1)
			continue;

		mate->rev_com = s->rev_com;
		if (s->rev_com)
			s->cursor = ori ? (s->len - query->len - 1 - s->pos) : (s->len
					- s->pos);
		else
			s->cursor = ori ? (s->pos - 1) : (s->pos + query->len);

		if (s->cursor >= s->len || s->cursor < 0) {
			s->cursor = 0;
			continue;
		}
		//if ((!ori && pre_cursor > s->cursor) || (ori && pre_cursor < s->cursor))
		//	continue;
		pool_add(cur_pool, s, eg->tid);
	}
	//show_debug_msg(__func__, "Removing partial...\n");
	rm_partial(eg, cur_pool, ori, seqs, query, MISMATCHES);
	//p_pool("POOL after removing partial", cur_pool, next);
	check_next_char(cur_pool, eg, next, ori);
}

void kmer_ext_edge(edge *eg, uint64_t query_int, hash_map *hm, const int ori) {
	int c = 0;
	uint64_t rev_kmer_int = 0;
	bwa_seq_t *debug = NULL;

	if (ori)
		seq_reverse(eg->len, eg->contig->seq, 0);
	while (1) {
		c = next_char_by_kmers(hm, query_int, ori);

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

		if (c == -1) {
			show_debug_msg(__func__, "[%d, %d] No hits, stop here. \n", eg->id,
					eg->len);
			break;
		}
		mark_kmer_used(query_int, hm);
		rev_kmer_int = rev_comp_kmer(query_int, hm->o->k);
		mark_kmer_used(rev_kmer_int, hm);
		ext_con(eg->contig, c, 0);
		eg->len = eg->contig->len;
		query_int = shift_bit(query_int, c, hm->o->k, ori);
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
	eg = new_eg();
	g_mutex_lock(kmer_id_mutex);
	eg->id = kmer_ctg_id++;
	g_mutex_unlock(kmer_id_mutex);
	// Get a copy of the kmer
	kmer = get_kmer_seq(kmer_int, opt->k);
	eg->contig = new_seq(kmer, kmer->len, 0);
	eg->len = eg->contig->len;
	eg->tid = atoi(kmer->name);
	eg->name = strdup(kmer->name);

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
