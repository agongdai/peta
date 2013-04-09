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
#include "bwtaln.h"
#include "pool.h"
#include "edge.h"
#include "edgelist.h"
#include "pool.h"
#include "kmers.hpp"
#include "mate.h"
#include "pool.h"

using namespace std;

int kmer_ctg_id = 0;
int ins_size = 0;
int sd_ins_size = 0;
int kmer_n_threads = 0;
int kmer_len = 0;
char *kmer_out = NULL;

GMutex *kmer_id_mutex;
struct timespec kmer_start_time, kmer_finish_time;

typedef struct {
	hash_map *hm;
	GPtrArray *edge_list;
} kmer_thread_meta;

/**
 * Extend a kmer
 */
void *kmer_ext_thread(gpointer data, gpointer thread_params) {
	edge *eg = NULL;
	mer *m = NULL;
	bwa_seq_t *query = NULL, *kmer = NULL;
	int round_1_len = 0, round_2_len = 0;
	kmer_thread_meta *params = (kmer_thread_meta*) thread_params;
	mer_map *kmers = (params->kmers);
	m = (mer*) data;
	kmer = get_kmer_seq(m, params->meta->k);

	if (has_n(kmer, 4) || is_biased_q(kmer)) {
		bwa_free_read_seq(1, kmer);
		return NULL;
	}show_debug_msg(__func__,
			"========= Kmer: %" ID64 " ========= \n", get_kmer_int(kmer->seq, 0, 1, kmer->len));
	eg = new_eg();
	g_mutex_lock(kmer_id_mutex);
	eg->id = kmer_ctg_id++;
	g_mutex_unlock(kmer_id_mutex);
	// Get a copy of the kmer
	query = new_seq(kmer, kmer->len, 0);
	eg->contig = new_seq(query, query->len, 0);
	eg->len = eg->contig->len;
	//	p_query(__func__, kmer);
	//	p_query(__func__, query);
	kmer_ext_edge(eg, query, kmers, params->meta, params->ht, 0);
	round_1_len = eg->len;
	//p_ctg_seq("Contig", eg->contig);
	show_debug_msg(__func__, "Edge %d length: %d \n", eg->id, eg->len);

	query = new_seq(eg->contig, kmer->len, 0);
	kmer_ext_edge(eg, query, kmers, params->meta, params->ht, 1);
	//p_ctg_seq("Contig", eg->contig);
	show_debug_msg(__func__, "Edge %d length: %d \n", eg->id, eg->len);

	if (eg->len - round_1_len > 2) {
		query = new_seq(eg->contig, kmer->len, 0);
		round_2_len = eg->len;
		kmer_ext_edge(eg, query, kmers, params->meta, params->ht, 0);
		show_debug_msg(__func__, "Edge %d length: %d \n", eg->id, eg->len);
		if (eg->len - round_2_len > 2) {
			query = new_seq(eg->contig, kmer->len, 0);
			kmer_ext_edge(eg, query, kmers, params->meta, params->ht, 1);
			show_debug_msg(__func__, "Edge %d length: %d \n", eg->id, eg->len);
		}
	}

	eg->start_kmer_int = get_kmer_int(kmer->seq, 0, 1, kmer->len);
	if (eg->len < 100000)
		destroy_eg(eg);
	else {
		g_mutex_lock(kmer_id_mutex);
		g_ptr_array_add(params->edge_list, eg);
		g_mutex_unlock(kmer_id_mutex);
	}
	bwa_free_read_seq(1, kmer);

	return NULL;
}

void kmer_threads(kmer_thread_meta *params, GPtrArray *kmer_list) {
	GThreadPool *thread_pool = NULL;
	mer *m = NULL;
	uint64_t i = 0;
	thread_pool = g_thread_pool_new((GFunc) kmer_ext_thread, params, 1, TRUE,
			NULL);

	show_msg(__func__, "Extending by kmers...\n");
	for (i = 0; i < kmer_list->len; i++) {
		m = (mer*) g_ptr_array_index(kmer_list, i);
		if (m->status != USED && m->count > 1) {
			g_thread_pool_push(thread_pool, (gpointer) m, NULL);
		}
	}
	g_thread_pool_free(thread_pool, 0, 1);
}

void test_kmer_ext(kmer_thread_meta *params, GPtrArray *kmer_list) {
	FILE *kmer_contigs = NULL;
	bwa_seq_t *kmer_seq = NULL;
	mer_map *map = params->kmers;

	mer *m = (mer*) (*map)[53073660363574];
	kmer_seq = get_kmer_seq(m, params->meta->k);
	kmer_ext_thread((gpointer) kmer_seq, (gpointer) params);

	m = (mer*) (*map)[857693345222861];
	kmer_seq = get_kmer_seq(m, params->meta->k);
	kmer_ext_thread((gpointer) kmer_seq, (gpointer) params);

	kmer_contigs = xopen("../SRR097897_out/single.fa", "w");
	save_edges(params->edge_list, kmer_contigs, 0, 0, 0);
	fflush(kmer_contigs);
	fclose(kmer_contigs);
}

void ext_by_kmers_core(char *lib_file, const char *solid_file,
		const char *kmer_file, const int insert_size, const int sd_insert_size,
		const int n_threads) {
	hash_table *ht = NULL;
	GPtrArray *all_edges = NULL, *kmer_list = NULL;
	mer_meta *meta = NULL;
	FILE *kmer_contigs = NULL;
	mer_map kmers;
	kmer_thread_meta *params = (kmer_thread_meta*) calloc(1,
			sizeof(kmer_thread_meta));

	clock_gettime(CLOCK_MONOTONIC, &kmer_start_time);

	show_msg(__func__, "Library: %s \n", lib_file);
	show_msg(__func__, "Solid Reads: %s \n", solid_file);
	ins_size = insert_size;
	sd_ins_size = sd_insert_size;

	kmer_list = g_ptr_array_sized_new(BUFSIZ);
	meta = (mer_meta*) malloc(sizeof(mer_meta));
	load_kmers(kmer_file, kmers, kmer_list, meta);
	clock_gettime(CLOCK_MONOTONIC, &kmer_finish_time);
	show_msg(__func__, "Done loading kmers: %.2f sec\n",
			(float) (kmer_finish_time.tv_sec - kmer_start_time.tv_sec));

	ht = pe_load_hash(lib_file);
	all_edges = g_ptr_array_sized_new(BUFSIZ);
	params->ht = ht;
	params->kmers = &kmers;
	params->meta = meta;
	params->edge_list = all_edges;

	//test_kmer_ext(params, kmer_list);
	//exit(1);
	kmer_threads(params, kmer_list);

	kmer_contigs = xopen("../SRR097897_out/kmer_contigs.fa", "w");
	save_edges(all_edges, kmer_contigs, 0, 0, 100);
	fflush(kmer_contigs);
	fclose(kmer_contigs);
	g_ptr_array_free(kmer_list, TRUE);

	reset_edge_ids(all_edges);
	merge_ol_edges(all_edges, insert_size, sd_insert_size, ht, n_threads);
	destroy_ht(ht);
	free(params);

	kmer_contigs = xopen("../SRR097897_out/merged.fa", "w");
	save_edges(all_edges, kmer_contigs, 0, 0, 100);
	fflush(kmer_contigs);
	fclose(kmer_contigs);
	kmer_contigs = xopen("../SRR097897_out/merged_all.fa", "w");
	save_edges(all_edges, kmer_contigs, 0, 0, 0);
	fflush(kmer_contigs);
	fclose(kmer_contigs);
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

	hash_opt *opt = (hash_opt*) malloc(sizeof(hash_opt));
	opt->k = kmer_len;
	test_kmer_hash(
			"/home/carl/Projects/peta/rnaseq/Spombe/SRR097897/SRR097897_corrected.fa");
	//build_kmers_hash(argv[optind], opt);
	//ext_by_kmers_core(argv[optind], argv[optind + 1], argv[optind + 2],
	//		ins_size, sd_ins_size, kmer_n_threads);

	clock_gettime(CLOCK_MONOTONIC, &kmer_finish_time);
	show_msg(__func__, "Done: %.2f sec\n", (float) (kmer_finish_time.tv_sec
			- kmer_start_time.tv_sec));
	return 0;
}
