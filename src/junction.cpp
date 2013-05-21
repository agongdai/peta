/*
 * junction.c
 *
 *  Created on: 09-May-2013
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
#include "kmers.hpp"
#include "ass.hpp"
#include "rnaseq.h"
#include "pechar.h"
#include "utils.h"
#include "peseq.h"
#include "junction.hpp"

junction *new_junction(edge *main_tpl, edge *branch_tpl, uint64_t kmer, int locus,
		int ori, int weight) {
	junction *j = (junction*) malloc(sizeof(junction));
	j->main_tpl = main_tpl;
	j->branch_tpl = branch_tpl;
	j->locus = main_tpl->len;
	j->ori = ori;
	j->locus = locus;
	j->weight = weight;
	j->kmer = kmer;
	return j;
}

/**
 * Check whether there are reads in the junction area
 */
int find_junc_reads(hash_map *hm, bwa_seq_t *left, bwa_seq_t *right,
		const int max_len, int *weight) {
	int left_len = 0, right_len = 0, n_reads = 0;
	GPtrArray *reads = NULL;
	bwa_seq_t *junction_seq = blank_seq(max_len);

	left_len = (left->len > max_len / 2) ? (max_len / 2) : left->len;
	memcpy(junction_seq->seq, left->seq + (left->len - left_len),
			sizeof(ubyte_t) * left_len);
	right_len = (right->len) > (max_len / 2) ? (max_len / 2) : (right->len);
	memcpy(junction_seq->seq + left_len, right->seq, sizeof(ubyte_t)
			* right_len);
	junction_seq->len = left_len + right_len;
	set_rev_com(junction_seq);
	// p_query("Junction seq", junction_seq);
	reads = kmer_find_reads(junction_seq, hm, 0, 0);
	n_reads = reads->len;
	// show_debug_msg(__func__, "# of junction reads: %d \n", n_reads);
	// p_readarray(reads, 1);
	*weight = n_reads;
	bwa_free_read_seq(1, junction_seq);
	g_ptr_array_free(reads, TRUE);
	if (n_reads > 0)
		return 1;
	return 0;
}

int find_junc_reads_w_tails(hash_map *hm, edge *left, edge *right,
		const int r_shift, const int max_len, int *weight) {
	bwa_seq_t *left_seq = NULL, *right_seq = NULL;
	int is_valid = 0;
	left_seq = cut_edge_tail(left, left->len, max_len / 2, 0);
	right_seq = cut_edge_tail(right, r_shift, max_len / 2, 1);
	is_valid = find_junc_reads(hm, left_seq, right_seq, max_len, weight);
	bwa_free_read_seq(1, left_seq);
	bwa_free_read_seq(1, right_seq);
	return is_valid;
}

/**
 * Update the junction locus for those edges right connected to itself.
 * Because the locus is not correct when the junction is recorded.
 */
void upd_tpl_jun_locus(edge *eg, GPtrArray *branching_events,
		const int kmer_len) {
	int i = 0;
	uint64_t query_int = 0, k = 0;
	junction *jun = NULL;
	if (eg->len < kmer_len)
		return;
	for (k = 0; k < branching_events->len; k++) {
		jun = (junction*) g_ptr_array_index(branching_events, k);
		if (jun->main_tpl == eg && jun->branch_tpl == eg) {
			for (i = 0; i < eg->len - kmer_len; i++) {
				query_int = get_kmer_int(eg->ctg->seq, i, 1, kmer_len);
				if (query_int == jun->kmer) {
					jun->locus = i;
					break;
				}
			}
		}
	}
}

gint cmp_junctions_by_id(gpointer a, gpointer b) {
	junction *c_a = *((junction**) a);
	junction *c_b = *((junction**) b);
	return ((c_a->main_tpl->id) - c_b->main_tpl->id);
}

void store_junctions(char *name, GPtrArray *branching_events) {
	junction *jun = NULL;
	uint64_t i = 0;
	char entry[BUFSIZE];
	FILE *f = xopen(name, "w");
	sprintf(entry, "Main\tBranch\tLocus\tWeight\tDirection\n");
	fputs(entry, f);
	for (i = 0; i < branching_events->len; i++) {
		jun = (junction*) g_ptr_array_index(branching_events, i);
		sprintf(entry, "[%d, %d]\t[%d, %d]\t%d\t%d\t%d\n", jun->main_tpl->id,
				jun->main_tpl->len, jun->branch_tpl->id, jun->branch_tpl->len,
				jun->locus, jun->weight, jun->ori);
		fputs(entry, f);
	}
	fclose(f);
}
