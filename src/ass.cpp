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
#include "pool.h"
#include "edge.h"
#include "edgelist.h"
#include "pool.h"
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
	rm_partial(eg, cur_pool, ori, seqs, query, MISMATCHES);
	//p_pool("POOL after removing partial", cur_pool, next);
	check_next_char(cur_pool, eg, next, ori);
}

void kmer_ext_edge(edge *eg, bwa_seq_t *query, hash_map *hm, const int ori) {
	pool *cur_pool = NULL;
	int c = 0, *next = NULL;
	uint64_t kmer_int = 0;
	GPtrArray *hits = NULL;

	if (ori)
		seq_reverse(eg->len, eg->contig->seq, 0);
	next = (int*) calloc(5, sizeof(int));
	cur_pool = new_pool();
	p_query(__func__, query);
	while (1) {
		//p_query(__func__, query);
		reset_c(next, NULL); // Reset the counter
		hits = kmer_aln_query(query, hm);
		//p_readarray(hits, 1);
		kmer_pool(hits, hm, cur_pool, eg, query, next, ori);
		if (cur_pool->reads->len <= 0 && eg->len > ins_size - SD_TIMES
				* sd_ins_size) {
			//show_debug_msg(__func__, "Trying mate pool... \n");
			add_mates_by_ol(hm->seqs, eg, cur_pool, RELAX_MATE_OL_THRE,
					SHORT_MISMATCH, query, ori, ins_size, sd_ins_size);
			reset_c(next, NULL); // Reset the counter
			check_next_char(cur_pool, eg, next, ori);
			//p_pool("After adding mates", cur_pool, next);
		}
		c = get_pure_most(next);
		//show_debug_msg(__func__, "Ori %d, Edge %d, length %d \n", ori, eg->id,
		//		eg->len);
		//p_ctg_seq("Contig", eg->contig);
		//p_pool("Current Pool", cur_pool, next);
		//show_debug_msg(__func__, "Ori: %d, Next char: %d \n", ori, c);
		if (cur_pool->n <= 0) {
			show_debug_msg(__func__, "[%d, %d] No hits, stop here. \n", eg->id,
					eg->len);
			break;
		}
		kmer_int = get_kmer_int(query->seq, 0, 1, query->len);
		mark_kmer_used(kmer_int, hm);
		forward(cur_pool, c, eg, ori);
		ext_con(eg->contig, c, 0);
		eg->len = eg->contig->len;
		ext_que(query, c, ori);
		g_ptr_array_free(hits, TRUE);
		if (eg->len % 100 == 0)
			show_debug_msg(__func__, "Ori %d, Edge %d, length %d \n", ori,
					eg->id, eg->len);
	}
	g_ptr_array_free(hits, TRUE);
	if (ori)
		seq_reverse(eg->len, eg->contig->seq, 0);
	bwa_free_read_seq(1, query);
	free_pool(cur_pool);
	free(next);
}

/**
 * Extend a kmer
 */
void *kmer_ext_thread(gpointer data, gpointer thread_params) {
	edge *eg = NULL;
	bwa_seq_t *query = NULL, *kmer = NULL;
	int round_1_len = 0, round_2_len = 0;
	map_opt *opt = NULL;
	kmer_t_meta *params = (kmer_t_meta*) thread_params;

	query = (bwa_seq_t*) data;
	if (query->status != FRESH || has_n(query, 4) || is_biased_q(query)) {
		return NULL;
	}
	opt = params->hm->o;
	eg = new_eg();
	g_mutex_lock(kmer_id_mutex);
	eg->id = kmer_ctg_id++;
	g_mutex_unlock(kmer_id_mutex);
	// Get a copy of the kmer
	kmer = new_seq(query, opt->k, query->len - opt->k);
	eg->contig = new_seq(query, query->len, 0);
	eg->len = eg->contig->len;
	eg->tid = atoi(query->name);
	eg->name = strdup(query->name);

	query->shift = 0;
	readarray_add(eg, query);
	show_debug_msg(__func__, "============= %s ============ \n", query->name);
	kmer_ext_edge(eg, kmer, params->hm, 0);
	round_1_len = eg->len;
	show_debug_msg(__func__, "Edge %d length: %d \n", eg->id, eg->len);

	kmer = new_seq(eg->contig, opt->k, 0);
	rev_reads_pos(eg);
	kmer_ext_edge(eg, kmer, params->hm, 1);
	rev_reads_pos(eg);
	show_debug_msg(__func__, "Edge %d length: %d \n", eg->id, eg->len);

	if (eg->len - round_1_len > 2) {
		kmer = new_seq(eg->contig, opt->k, eg->len - opt->k);
		round_2_len = eg->len;
		kmer_ext_edge(eg, kmer, params->hm, 0);
		show_debug_msg(__func__, "Edge %d length: %d \n", eg->id, eg->len);
		if (eg->len - round_2_len > 2) {
			kmer = new_seq(eg->contig, opt->k, 0);
			rev_reads_pos(eg);
			kmer_ext_edge(eg, kmer, params->hm, 1);
			rev_reads_pos(eg);
			show_debug_msg(__func__, "Edge %d length: %d \n", eg->id, eg->len);
		}
	}

	eg->start_kmer_int = get_kmer_int(kmer->seq, 0, 1, kmer->len);
	if (eg->len < 50)
		destroy_eg(eg);
	else {
		g_mutex_lock(kmer_id_mutex);
		g_ptr_array_add(params->all_edges, eg);
		n_used_reads += eg->reads->len;
		g_mutex_unlock(kmer_id_mutex);
	}
	show_debug_msg(__func__, "Consumed reads: %d/%" ID64 " \n", n_used_reads, params->hm->n_reads);
	return NULL;
}

void kmer_threads(kmer_t_meta *params, GPtrArray *solid_reads) {
	GThreadPool *thread_pool = NULL;
	uint64_t i = 0;
	bwa_seq_t *query = NULL;
	thread_pool = g_thread_pool_new((GFunc) kmer_ext_thread, params, 1, TRUE,
			NULL);

	show_msg(__func__, "Extending by kmers...\n");
	for (i = 0; i < solid_reads->len; i++) {
		query = (bwa_seq_t*) g_ptr_array_index(solid_reads, i);
		g_thread_pool_push(thread_pool, (gpointer) query, NULL);
		//break;
	}
	g_thread_pool_free(thread_pool, 0, 1);
}

void pick_unused_kmers(hash_map *hm) {
	mer_hash *hash = hm->hash;
	uint64_t mer_v = 0, *freq = NULL, n_kmers = 0, n = 0;
	bwa_seq_t *seq = NULL;
	for (mer_hash::iterator m = hash->begin(); m != hash->end(); ++m) {
		mer_v = m->first;
		freq = m->second;
		n = get_kmer_count(mer_v, hm);
		if (n > 0) {
			n_kmers++;
			seq = get_kmer_seq(mer_v, hm->o->k);
			p_query(__func__, seq);
			show_debug_msg(__func__, "Count: %" ID64 "\n", n);
		}
	}
	show_msg(__func__, "Remaining kmers: %" ID64 "\n", n_kmers);
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

	hm = load_hash_map(lib_file, map);
	hm->hash = &map;

	clock_gettime(CLOCK_MONOTONIC, &kmer_finish_time);
	show_msg(__func__, "Done loading kmers hash map: %.2f sec\n",
			(float) (kmer_finish_time.tv_sec - kmer_start_time.tv_sec));

	all_edges = g_ptr_array_sized_new(BUFSIZ);
	params->hm = hm;
	params->all_edges = all_edges;
	solid_reads = load_solid_reads(solid_file, hm->seqs, hm->n_reads);

	//test_kmer_ext(params, kmer_list);
	//exit(1);
	kmer_threads(params, solid_reads);
	g_ptr_array_free(solid_reads, TRUE);

	contigs = xopen("../SRR097897_out/paired.fa", "w");
	save_edges(all_edges, contigs, 0, 0, 100);
	fflush(contigs);
	fclose(contigs);

	pick_unused_kmers(hm);

	reset_edge_ids(all_edges);
	merge_ol_edges(all_edges, ins_size, sd_ins_size, hm->seqs, kmer_n_threads);
	free(params);
	destroy_hm(hm);

	contigs = xopen("../SRR097897_out/merged.fa", "w");
	save_edges(all_edges, contigs, 0, 0, 100);
	fflush(contigs);
	fclose(contigs);
	contigs = xopen("../SRR097897_out/merged_all.fa", "w");
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
