/*
 * edge.c
 *
 *  Created on: 06-Mar-2012
 *      Author: carl
 */

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <glib.h>
#include "utils.h"
#include "edge.h"
#include "peseq.h"
#include "bwtaln.h"
#include "pehash.h"

eg_gap *init_gap(int s_index, int size, int ori) {
	eg_gap *gap = (eg_gap*) malloc (sizeof(eg_gap));
	gap->s_index = s_index;
	gap->size = size;
	gap->ori = ori;
	return gap;
}

void free_eg_gap(eg_gap *gap) {
	if (gap)
		free(gap);
}

edge *new_eg() {
	edge *eg = (edge*) malloc(sizeof(edge));
	eg->contig = 0;
	eg->in_egs = g_ptr_array_sized_new(0);
	eg->out_egs = g_ptr_array_sized_new(0);
	eg->reads = g_ptr_array_sized_new(0);
	eg->pairs = g_ptr_array_sized_new(0);
	eg->name = NULL;
	eg->right_ctg = NULL;
	eg->left_ctg = NULL;
	eg->len = 0;
	eg->r_shift = 0;
	eg->l_shift = 0;
	eg->id = 0;
	eg->visited = 0;
	eg->alive = 1;
	eg->is_root = 0;
	eg->ori = 0;
	eg->gaps = g_ptr_array_sized_new(0);
	eg->level = -1;
	eg->comp_id = -1;
	eg->start_kmer_int = 0;
	return eg;
}

/**
 * Reset thread id of all edges
 */
void reset_tid(edge *eg) {
	int i = 0;
	bwa_seq_t *r = NULL;
	for (i = 0; i < eg->reads->len; i++) {
		r = g_ptr_array_index(eg->reads, i);
		if (r->status != USED && r->status != DEAD)
			r->tid = -1;
		else
			r->tid = eg->tid;
	}
}

void free_readarray(readarray *ra) {
	if (!ra)
		return;
	g_ptr_array_free(ra, TRUE);
}

void destroy_eg(edge *eg) {
	eg_gap *gap = NULL;
	bwa_seq_t *read = NULL;
	int i = 0;
	if (eg) {
		bwa_free_read_seq(1, eg->contig); // bug if free it
		g_ptr_array_free(eg->in_egs, TRUE);
		if (!eg->right_ctg) {
			// If eg's right contig is not null, its out_egs is set to be right contig's out_egs
			g_ptr_array_free(eg->out_egs, TRUE);
		}
		for (i = 0; i < eg->reads->len; i++) {
			read = g_ptr_array_index(eg->reads, i);
			read->status = TRIED;
			read->contig_id = UNUSED_CONTIG_ID;
		}
		for (i = 0; i < eg->pairs->len; i++) {
			read = g_ptr_array_index(eg->pairs, i);
			read->status = TRIED;
			read->contig_id = UNUSED_CONTIG_ID;
		}
		while (eg->reads->len > 0) {
			g_ptr_array_remove_index_fast(eg->reads, 0);
		}
		while (eg->pairs->len > 0) {
			g_ptr_array_remove_index_fast(eg->pairs, 0);
		}
		free_readarray(eg->reads);
		free_readarray(eg->pairs);
		for (i = 0; i < eg->gaps->len; i++) {
			gap = g_ptr_array_index(eg->gaps, i);
			free_eg_gap(gap);
		}
		free_readarray(eg->gaps);
		eg->alive = 0;
		free(eg->name);
		free(eg);
	}
}
