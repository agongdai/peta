/*
 * kmer.c
 *
 *  Created on: 13-Mar-2013
 *      Author: carl
 */

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <glib.h>
#include <time.h>
#include "utils.h"
#include "pehash.h"
#include "bwtaln.h"
#include "kmer.h"
#include "pehash.h"
#include "rnaseq.h"
#include "list.h"
#include "peseq.h"
#include "roadmap.h"
#include "readrm.h"
#include "mate.h"
#include "merge.h"

int kmer_ctg_id = 0;
int ins_size = 0;
int sd_ins_size = 0;
struct timespec kmer_start_time, kmer_finish_time;

slist *prepare_kmers(bwa_seq_t *kmers, const int n_kmers) {
	int i = 0;
	bwa_seq_t *m = NULL;
	slist *kmer_list = NULL;
	kmer_list = new_slist();
	for (i = 0; i < n_kmers; i++) {
		m = &kmers[i];
		g_ptr_array_add(kmer_list->values, m);
		if (atoi(m->name) > 1)
			g_ptr_array_add(kmer_list->kmers, m);
	}
	show_msg(__func__, "Sorting kmers by sequences...\n");
	g_ptr_array_sort(kmer_list->values, (GCompareFunc) cmp_kmer_by_seq);
	show_msg(__func__, "Sorting kmers by frequency...\n");
	g_ptr_array_sort(kmer_list->kmers, (GCompareFunc) cmp_kmer_by_freq);
	kmer_list->order = 0;
	return kmer_list;
}

/**
 * Check kmer frequencies of all possible next chars, get the highest one
 * Return: ACGTN -> 01234;
 *         None -> -1;
 */
int next_char_by_kmers(slist *kmer_list, bwa_seq_t *query, const int ori) {
	bwa_seq_t *next_q = NULL, *mer = NULL;
	int counters[5], i = 0, max = 0, exists = 0;
	for (i = 0; i < 5; i++) {
		next_q = new_seq(query, query->len, 0);
		ext_que(next_q, i, ori);
		counters[i] = 0;
		// Check the forward kmer
		//p_query("SEQ", next_q);
		exists = slist_binary(kmer_list, next_q);
		if (exists != -1) {
			mer = g_ptr_array_index(kmer_list->values, exists);
			//p_query("FORWARD", mer);
			if (mer->status != USED)
				counters[i] = atoi(mer->name);
		}
		// Check the reverse kmer
		switch_ubyte(next_q);
		exists = slist_binary(kmer_list, next_q);
		//p_query("REV", next_q);
		if (exists != -1) {
			mer = g_ptr_array_index(kmer_list->values, exists);
			//p_query("BACKWARD", mer);
			if (mer->status != USED)
				counters[i] += atoi(mer->name);
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

/**
 * Mark a kmer as used.
 * The parameter 'kmer' is a 'copy' only, need to mark the kmer on the list
 */
void mark_kmer_used(slist *kmer_list, bwa_seq_t *kmer) {
	bwa_seq_t *ori_kmer = NULL;
	int pos = 0;
	pos = slist_binary(kmer_list, kmer);
	if (pos > -1) {
		ori_kmer = g_ptr_array_index(kmer_list->values, pos);
		ori_kmer->status = USED;
	}
	switch_ubyte(kmer);
	pos = slist_binary(kmer_list, kmer);
	if (pos > -1) {
		ori_kmer = g_ptr_array_index(kmer_list->values, pos);
		ori_kmer->status = USED;
	}
	switch_ubyte(kmer); // Switch back for the query
}

int ext_eg_query_n(slist *kmer_list, bwa_seq_t *query, edge *eg,
		bwa_seq_t *read, const int shift, const int ori) {
	int i = 0, pos = 0, ext_len = 0;
	ubyte_t *s = NULL, next_c = 0;
	bwa_seq_t *mer = NULL;
	if (shift > query->len)
		return 0;
	s = read->seq;
	for (i = shift; i < read->len; i++) {
		next_c = s[i];
		ext_que(query, next_c, ori);
//		pos = slist_binary(kmer_list, query);
//		if (pos >= 0) {
//			mer = g_ptr_array_index(kmer_list->values, pos);
//			if (mer->status = USED)
//				break;
//		}
		ext_con(eg->contig, next_c, 0);
		//p_ctg_seq(__func__, eg->contig);
		mark_kmer_used(kmer_list, query);
	}
	ext_len = eg->contig->len - eg->len;
	eg->len = eg->contig->len;
	return ext_len;
}

int ext_by_mates(slist *kmer_list, edgearray *mates, bwa_seq_t *query,
		edge *eg, const int ori) {
	int ol = 0, i = 0, ext_len = 0;
	bwa_seq_t *m = NULL, *tmp = NULL, *template = NULL;
	template = new_seq(eg->contig, query->len, eg->len - query->len);
	if (ori) {
		seq_reverse(template->len, template->seq, 0);
	}
	for (i = 0; i < mates->len; i++) {
		m = g_ptr_array_index(mates, i);
		tmp = m;
		if (m->rev_com)
			tmp = new_mem_rev_seq(m, m->len, 0);
		ol = find_ol_within_k(tmp, template, 1, 0, template->len, ori);
		//p_ctg_seq(__func__, tmp);
		//p_ctg_seq(__func__, template);
		//show_debug_msg(__func__, "OL: %d \n", ol);
		if (ol >= 8) {
			//p_query("Mate", m);
			ext_len = ext_eg_query_n(kmer_list, query, eg, tmp, ol, ori);
		}

		if (m->rev_com)
			bwa_free_read_seq(1, tmp);
		if (ext_len > 0)
			break;
	}
	bwa_free_read_seq(1, template);
	return ext_len;
}

void kmer_ext_edge(edge *eg, bwa_seq_t *query, slist *kmer_list,
		hash_table *ht, const int ori) {
	int next_c = 0;
	pool *mate_pool = NULL;
	if (ori)
		seq_reverse(eg->len, eg->contig->seq, 0);
	p_query(__func__, query);
	while (1) {
		if (eg->len % 100 == 0)
			show_debug_msg(__func__, "Edge [%d, %d] \n", eg->id, eg->len);
		//p_query(__func__, query);
		//p_ctg_seq("Contig", eg->contig);
		next_c = next_char_by_kmers(kmer_list, query, ori);
		//show_debug_msg(__func__, "Next char: %d \n", next_c);
		if (next_c == -1) {
			show_debug_msg(__func__, "Edge [%d, %d] \n", eg->id, eg->len);
			show_debug_msg(__func__, "Realigning \n");
			realign_reads_by_ht(ht, eg, MISMATCHES, ori);
			//p_readarray(eg->reads, 1);
			show_debug_msg(__func__, "Getting mate pool \n");
			mate_pool = get_mate_pool_from_edge(eg, ht, ori, ins_size,
					sd_ins_size);
			//p_readarray(mate_pool->reads, 1);
			show_debug_msg(__func__, "Extending by mates \n");
			if (ext_by_mates(kmer_list, mate_pool->reads, query, eg, ori))
				continue;
			else
				break;
		}
		ext_con(eg->contig, next_c, 0);
		eg->len = eg->contig->len;
		ext_que(query, next_c, ori);
		mark_kmer_used(kmer_list, query);
	}
	if (ori)
		seq_reverse(eg->len, eg->contig->seq, 0);
	bwa_free_read_seq(1, query);
}

/**
 * Extend a kmer
 */
edge *kmer_ext(slist *kmer_list, bwa_seq_t *kmer, hash_table *ht) {
	edge *eg = NULL;
	bwa_seq_t *query = NULL;
	int round_1_len = 0, round_2_len = 0;

	eg = new_eg();
	eg->id = kmer_ctg_id++;
	// Get a copy of the kmer
	query = new_seq(kmer, kmer->len, 0);
	eg->contig = new_seq(query, query->len, 0);
	kmer_ext_edge(eg, query, kmer_list, ht, 0);
	round_1_len = eg->len;
	//p_ctg_seq("Contig", eg->contig);
	show_debug_msg(__func__, "Edge %d length: %d \n", eg->id, eg->len);

	query = new_seq(eg->contig, kmer->len, 0);
	kmer_ext_edge(eg, query, kmer_list, ht, 1);
	//p_ctg_seq("Contig", eg->contig);
	show_debug_msg(__func__, "Edge %d length: %d \n", eg->id, eg->len);

	if (eg->len - round_1_len > 2) {
		query = new_seq(eg->contig, kmer->len, 0);
		round_2_len = eg->len;
		kmer_ext_edge(eg, query, kmer_list, ht, 0);
		show_debug_msg(__func__, "Edge %d length: %d \n", eg->id, eg->len);
		if (eg->len - round_2_len > 2) {
			query = new_seq(eg->contig, kmer->len, 0);
			kmer_ext_edge(eg, query, kmer_list, ht, 1);
			show_debug_msg(__func__, "Edge %d length: %d \n", eg->id, eg->len);
		}
	}

	return eg;
}

void ext_by_kmers(char *lib_file, const char *solid_file,
		const char *kmer_file, const int insert_size, const int sd_insert_size,
		const int n_threads) {
	hash_table *ht = NULL;
	bwa_seq_t *kmers = NULL, *mer = NULL;
	uint32_t n_kmers = 0;
	GPtrArray *all_edges = NULL;
	slist *kmer_list = NULL;
	int i = 0;
	edge *eg = NULL;
	FILE *kmer_contigs = NULL;
	clock_gettime(CLOCK_MONOTONIC, &kmer_start_time);

	show_msg(__func__, "Library: %s \n", lib_file);
	show_msg(__func__, "Solid Reads: %s \n", solid_file);
	ins_size = insert_size;
	sd_ins_size = sd_insert_size;

	kmers = load_reads(kmer_file, &n_kmers);
	all_edges = g_ptr_array_sized_new(BUFSIZ); // kmer_list->kmers->len

	show_msg(__func__, "Preparing kmer list...\n");
	kmer_list = prepare_kmers(kmers, n_kmers);

	clock_gettime(CLOCK_MONOTONIC, &kmer_finish_time);
	show_msg(__func__, "Done preparation: %.2f sec\n",
			(float) (kmer_finish_time.tv_sec - kmer_start_time.tv_sec));

	ht = pe_load_hash(lib_file);
	show_msg(__func__, "Extending by kmers...\n");
	for (i = 0; i < kmer_list->kmers->len; i++) {
		mer = g_ptr_array_index(kmer_list->kmers, i);
		if (mer->status != USED && !same_bytes(mer)) {
			show_debug_msg(__func__,
					"========================== %d ===================== \n", i);
			eg = kmer_ext(kmer_list, mer, ht);
			if (eg->len > 100)
				g_ptr_array_add(all_edges, eg);
		}
	}
	kmer_contigs = xopen("../SRR097897_out/kmer_contigs.fa", "w");
	save_edges(all_edges, kmer_contigs, 0, 0, 100);
	fflush(kmer_contigs);
	fclose(kmer_contigs);
	free_slist(kmer_list);
	bwa_free_read_seq(n_kmers, kmers);

	reset_edge_ids(all_edges);
	merge_ol_edges(all_edges, insert_size, sd_insert_size, ht, n_threads);
	destroy_ht(ht);
	kmer_contigs = xopen("../SRR097897_out/merged.fa", "w");
	save_edges(all_edges, kmer_contigs, 0, 0, 100);

	clock_gettime(CLOCK_MONOTONIC, &kmer_finish_time);
	show_msg(__func__, "Done: %.2f sec\n", (float) (kmer_finish_time.tv_sec
			- kmer_start_time.tv_sec));
}
