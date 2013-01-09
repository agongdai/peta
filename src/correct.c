/*
 * correct.c
 *
 *  Created on: 06-Jan-2013
 *      Author: carl
 */

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <glib.h>
#include <math.h>
#include "bwase.h"
#include "utils.h"
#include "peseq.h"
#include "pechar.h"
#include "pealn.h"
#include "edge.h"
#include "edgelist.h"
#include "pehash.h"
#include "pool.h"
#include "clean.h"

int rm_repetitive_reads(bwa_seq_t *seqs, const int n_seqs) {
	bwa_seq_t *s;
	int i = 0, n_most = 0, k = 0, n_rep = 0;
	ubyte_t c = 0;
	int *counter = NULL;
	counter = (int*) calloc(5, sizeof(int));
	for (k = 0; k < n_seqs; k++) {
		s = &seqs[k];
		n_most = 0;
		for (i = 0; i < s->len; i++) {
			c = s->seq[i];
			counter[c]++;
		}
		for (i = 0; i < 5; i++) {
			if (counter[i] > n_most) {
				n_most = counter[i];
			}
			counter[i] = 0;
		}
		if (n_most >= s->len - 2) {
			s->status = DEAD;
			n_rep++;
		}
	}
	free(counter);
	return n_rep;
}

void correct_bases(bwa_seq_t *seqs, bwa_seq_t *ori_read, alignarray *aligns,
		const int tid) {
	pool *p = NULL;
	int j = 0, cursor = 0, i = 0, index = 0, has_hit = 0;
	alg *a = NULL;
	bwa_seq_t *s = NULL;
	int *counter = NULL;

	//p_query("ORI", ori_read);
	p = new_pool();
	ori_read->cursor = 0;
	for (i = 0; i < aligns->len; i++) {
		a = g_ptr_array_index(aligns, i);
		index = a->r_id;
		s = &seqs[index];
		if (s->is_in_c_pool > 0 && s->is_in_c_pool != tid)
			continue;
		if (s->status == USED || strcmp(s->name, ori_read->name) == 0)
			continue;
		s->rev_com = a->rev_comp;
		if (s->rev_com)
			s->cursor = s->len - ori_read->len - a->pos;
		else
			s->cursor = a->pos;
		if (s->is_in_c_pool == 0)
			pool_add(p, s, tid);
	}
	if (p->n >= 4) {
		//p_pool(__func__, p, NULL);
		counter = (int*) calloc(5, sizeof(int));
		for (j = 0; j < ori_read->len; j++) {
			reset_c(counter, NULL);
			has_hit = 0;
			for (i = 0; i < p->n; i++) {
				s = g_ptr_array_index(p->reads, i);
				cursor = s->cursor + j;
				if (cursor >= 0 || cursor < s->len) {
					has_hit = 1;
					if (s->rev_com) {
						counter[s->rseq[cursor]]++;
					} else {
						counter[s->seq[cursor]]++;
					}
				}
			}
			//show_debug_msg(__func__, "Correcting %d: %d:%d:%d:%d\n", j, counter[0],
			//		counter[1], counter[2], counter[3]);
			if (has_hit) {
				ori_read->seq[j] = get_pure_most(counter);
				ori_read->rseq[ori_read->len - 1 - j] = 3 - ori_read->seq[j];
			}
		}
		free(counter);
	}
	free_pool(p);
	//p_query("AFT", ori_read);
}

GPtrArray *get_low_kmer_reads(hash_table *ht, const int k, const double thre) {
	int i = 0, j = 0, h = 0, lowest = 0, highest = 0, key = 0;
	bwa_seq_t *seqs = NULL, *s = NULL;
	uint16_t *kmer_list = NULL;
	counter *max_kmer_diff = NULL, *c = NULL;
	uint32_t n_kmers = 0;
	GPtrArray *low_kmer_reads = NULL;

	seqs = ht->seqs;
	n_kmers = (1 << (k * 2)) + 1;
	kmer_list = (uint16_t*) calloc(n_kmers, sizeof(uint16_t));
	show_msg(__func__, "Counting %d-mers of %d reads... \n", k, ht->n_seqs);
	for (i = 0; i < ht->n_seqs; i++) {
		s = &seqs[i];
		if (s->status == USED || s->status == DEAD)
			continue;
		set_kmer_index(s, k, kmer_list);
	}
	low_kmer_reads = g_ptr_array_sized_new(ht->n_seqs / 2);
	max_kmer_diff = (counter*) calloc(ht->n_seqs, sizeof(counter));
	show_msg(__func__,
			"Getting maximum difference between highest/lowest %d-mers... \n",
			k);
	for (i = 0; i < ht->n_seqs; i++) {
		s = &seqs[i];
		lowest = 0;
		highest = 0;
		if (s->status == USED || s->status == DEAD)
			continue;
		for (j = 0; j <= s->len - k; j++) {
			key = 0;
			for (h = 0; h < k; h++) {
				key *= 4;
				key = key | s->seq[h + j];
			}
			if (kmer_list[key] > highest)
				highest = kmer_list[key];
			if (j == 0)
				lowest = kmer_list[key];
			if (kmer_list[key] < lowest)
				lowest = kmer_list[key];
		}
		c = &max_kmer_diff[i];
		c->read_id = atoi(s->name);
		c->k_freq = highest - lowest;
	}
	show_msg(__func__, "Sorting... \n");
	qsort(max_kmer_diff, ht->n_seqs, sizeof(counter), cmp_kmer);
	for (i = 0; i < ht->n_seqs; i++) {
		if (low_kmer_reads->len >= ht->n_seqs * thre)
			break;
		c = &max_kmer_diff[i];
		s = &seqs[c->read_id];
		if (has_n(s))
			continue;
		g_ptr_array_add(low_kmer_reads, s);
		// show_debug_msg(__func__, "%s: %f \n", s->name, c->k_freq);
	}
	show_msg(__func__, "# of low kmer reads: %d \n", low_kmer_reads->len);
	free(kmer_list);
	free(max_kmer_diff);
	return low_kmer_reads;
}

typedef struct {
	readarray *low_kmer_reads;
	hash_table *ht;
	int start;
	int end;
	int tid;
	double stop_thre;
} correct_aux_t;

static void *correct_thread(void *data) {
	correct_aux_t *d = (correct_aux_t*) data;
	int i = 0;
	bwa_seq_t *s = NULL, *query = NULL, *seqs = d->ht->seqs;
	readarray *low_kmer_reads = d->low_kmer_reads;
	alignarray *aligns = NULL;

	aligns = g_ptr_array_sized_new(N_DEFAULT_ALIGNS);
	for (i = d->start; i < d->end; i++) {
		if (i % 10000 == 0)
			show_msg(__func__,
					"Thread %d correction progress: [%d,%d,%d]... \n", d->tid,
					d->start, i, d->end);
		s = g_ptr_array_index(low_kmer_reads, i);
		if (is_repetitive_q(s)) {
			s->status = USED;
			continue;
		}
		// Only the fresh reads, or the reads tried once would be corrected.
		if (s->status != FRESH)
			continue;
		query = new_seq(s, s->len - 8, 0);
		pe_aln_query(s, s->seq, d->ht, MISMATCHES, s->len, 0, aligns);
		pe_aln_query(s, s->rseq, d->ht, MISMATCHES, s->len, 1, aligns);
		if (aligns->len >= 4)
			correct_bases(seqs, s, aligns, d->tid);
		s->status = TRIED;
		reset_alg(aligns);
		bwa_free_read_seq(1, query);
		//if (i > 10000)
		//	break;
	}
	free_alg(aligns);
	show_msg(__func__, "Thread %d finished. \n", d->tid);
}

int correct_reads(hash_table *ht, const int n_threads) {
	hash_opt *opt = NULL;
	bwa_seq_t *seqs = NULL, *s = NULL, *query = NULL;
	hash_key *k_mers_occ_acc = NULL;
	correct_aux_t *data;
	int i = 0, n_per_threads = 0;
	GThread *threads[n_threads];
	GPtrArray *low_kmer_reads = NULL;

	data = (correct_aux_t*) calloc(n_threads, sizeof(correct_aux_t));
	low_kmer_reads = get_low_kmer_reads(ht, 15, 0.1);

	n_per_threads = low_kmer_reads->len / n_threads;
	for (i = 0; i < n_threads; ++i) {
		data[i].low_kmer_reads = low_kmer_reads;
		data[i].ht = ht;
		data[i].start = n_per_threads * i;
		data[i].end = n_per_threads * (i + 1);
		if (data[i].end >= low_kmer_reads->len)
			data[i].end = low_kmer_reads->len - 1;
		data[i].tid = i + 1;
		data[i].stop_thre = 0;
		threads[i]
				= g_thread_create((GThreadFunc) correct_thread, data + i, TRUE, NULL);
	}

	/* wait for threads to finish */
	for (i = 0; i < n_threads; ++i) {
		g_thread_join(threads[i]);
	}
	free(data);
	show_msg(__func__, "Resetting the status... \n");

	// Reset the status to FRESH
	seqs = ht->seqs;
	for (i = 0; i < ht->n_seqs; i++) {
		s = &seqs[i];
		if (s->status != DEAD) {
			s->status = FRESH;
			s->is_in_c_pool = 0;
			s->tid = 0;
		}
	}
	g_ptr_array_free(low_kmer_reads, TRUE);
	return 0;
}
