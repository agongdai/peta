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
#include "pehash.h"
#include "bwase.h"
#include "kmer.h"
#include "pehash.h"
#include "rnaseq.h"
#include "list.h"
#include "peseq.h"
#include "roadmap.h"

slist *prepare_kmers(bwa_seq_t *kmers, const int n_kmers) {
	int i = 0;
	bwa_seq_t *m = NULL;
	slist *kmer_list = NULL;
	kmer_list = new_slist((GCompareFunc) cmp_kmer_by_seq);
	for (i = 0; i < n_kmers; i++) {
		m = &kmers[i];
		m = new_seq(m, m->len, 0);
		g_ptr_array_add(kmer_list->values, m);
	}
	g_ptr_array_sort(kmer_list->values, kmer_list->cmp_f);
	kmer_list->order = 0;
	return kmer_list;
}

/**
 * Check kmer frequencies of all possible next chars, get the highest one
 * Return: ACGTN -> 01234;
 *         None -> -1;
 */
ubyte_t next_char_by_kmers(slist *kmer_list, bwa_seq_t *query) {
	bwa_seq_t *next_q = NULL, *mer = NULL;
	int counters[5], i = 0, max = -1, exists = 0;
	for (i = 0; i < 5; i++) {
		next_q = new_seq(query, query->len, 0);
		ext_que(next_q, i, 0);
		p_query(__func__, next_q);
		exists = slist_binary(kmer_list, next_q);
		if (exists == -1) {
			counters[i] = -1;
		} else {
			mer = g_ptr_array_index(kmer_list->values, exists);
			p_query(__func__, mer);
			// The times this kmer appears
			if (mer->status != USED)
				counters[i] = atoi(mer->name);
		}
		bwa_free_read_seq(1, next_q);
	}
	// Get max occurrence of the chars
	for (i = 0; i < 5; i++) {
		max = (counters[i] > max) ? counters[i] : max;
	}
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
}

/**
 * Extend a kmer
 */
edge *kmer_ext(slist *kmer_list, bwa_seq_t *kmer) {
	edge *eg = NULL;
	int next_c = 0;
	bwa_seq_t *query = NULL;

	eg = new_eg();
	// Get a copy of the kmer
	query = new_seq(kmer, kmer->len, 0);
	eg->contig = new_seq(query, query->len, 0);
	while (1) {
		p_query(__func__, query);
		next_c = next_char_by_kmers(kmer_list, query);
		if (next_c == -1)
			break;
		ext_con(eg->contig, next_c, 0);
		eg->len = eg->contig->len;
		ext_que(query, next_c, 0);
		mark_kmer_used(kmer_list, query);
	}
	p_ctg_seq("Contig", eg->contig);
	return eg;
}

void ext_by_kmers(const char *lib_file, const char *solid_file,
		const char *kmer_file) {
	hash_table *ht = NULL;
	bwa_seq_t *kmers = NULL, *mer = NULL;
	uint32_t n_kmers = 0;
	GPtrArray *all_edges = NULL;
	slist *kmer_list = NULL;
	int i = 0;

	show_msg(__func__, "Library: %s \n", lib_file);
	show_msg(__func__, "Solid Reads: %s \n", solid_file);

	kmers = load_reads(kmer_file, &n_kmers);
	//	all_edges = g_ptr_array_sized_new(BUFSIZ);
	//	ht = pe_load_hash(lib_file);

	show_msg(__func__, "Preparing kmer list...\n");
	kmer_list = prepare_kmers(kmers, n_kmers);

	show_msg(__func__, "Extending by kmers...\n");
	mer = g_ptr_array_index(kmer_list->values, 0);
	kmer_ext(kmer_list, mer);
}
