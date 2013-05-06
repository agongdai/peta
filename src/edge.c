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
	eg_gap *gap = (eg_gap*) malloc(sizeof(eg_gap));
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
	eg->contig = NULL;
	eg->l_tail = NULL;
	eg->r_tail = NULL;
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
	eg->tid = 0;
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

/**
 * Get virtual tail of an edge.
 * This is used when a branch edge is connectted to the locus of '^'.
 * When another third-layer branch edge is connectted to this branch edge,
 * 	its own length may be not long enough, then it cuts some length from the main template.
 * If the length of the branch edge is long enough, the virtual tail is not used.
 *
 * Edge: 	==============================
 * shift: 	                    ^
 * ori: 	0 (to the right)
 * tail_len:	                --------------
 * Return:                      ==========++++
 * '++++' is the partial virtual tail of current edge
 */
bwa_seq_t *cut_edge_tail(edge *eg, const int shift, const int tail_len,
		const int ori) {
	bwa_seq_t *tail = NULL, *partial = NULL, *main_tail = NULL;
	int v_tail_len = 0;
	// Get partial edge at the locus 'shift'.
	show_debug_msg(__func__, "[%d, %d] Cutting edge tail at shift %d for tail length %d from ori %d...\n", eg->id, eg->len, shift, tail_len, ori);
	p_ctg_seq(__func__, eg->contig);
	if (ori) {
		partial = new_seq(eg->contig, eg->len - shift, shift);
		main_tail = eg->r_tail;
	} else {
		partial = new_seq(eg->contig, shift, 0);
		main_tail = eg->l_tail;
	}
	p_ctg_seq("MAIN_TA", main_tail);
	// If the edge is long, cut the tail directly
	if (partial->len >= tail_len) {
		if (ori)
			tail = new_seq(partial, tail_len, 0);
		else
			tail = new_seq(partial, tail_len, partial->len - tail_len);
	} else {
		// If the edge has a virtual tail, try to get it
		if (main_tail && main_tail->len > 0) {
			v_tail_len = main_tail->len + partial->len;
			v_tail_len = (v_tail_len > tail_len) ? tail_len : v_tail_len;
			tail = blank_seq(v_tail_len);
			if (ori) {
				memcpy(tail->seq, partial->seq, sizeof(ubyte_t) * partial->len);
				memcpy(tail->seq + partial->len, tail->seq, sizeof(ubyte_t)
						* (v_tail_len - partial->len));
			} else {
				show_debug_msg(__func__, "Virtual tail length %d\n", v_tail_len);
				p_ctg_seq("PARTIAL", partial);
				p_ctg_seq("MAIN_TA", main_tail);
				memcpy(tail->seq, main_tail->seq + (main_tail->len
						+ partial->len - v_tail_len), sizeof(ubyte_t)
						* (v_tail_len - partial->len));
				memcpy(tail->seq + (v_tail_len - partial->len), partial->seq,
						sizeof(ubyte_t) * partial->len);
			}
			set_rev_com(tail);
		} else
			tail = new_seq(partial, partial->len, 0);
	}
	bwa_free_read_seq(1, partial);
	p_ctg_seq("TAIL", tail);
	return tail;
}

void set_tail(edge *eg, const int shift, const int tail_len, const int ori) {
	if (ori) {
		bwa_free_read_seq(1, eg->r_tail);
		eg->r_tail = cut_edge_tail(eg, shift, tail_len, ori);
	} else {
		bwa_free_read_seq(1, eg->l_tail);
		eg->l_tail = cut_edge_tail(eg, shift, tail_len, ori);
	}
}

void destroy_eg(edge *eg) {
	eg_gap *gap = NULL;
	bwa_seq_t *read = NULL;
	int i = 0;
	if (eg) {
		bwa_free_read_seq(1, eg->contig);
		bwa_free_read_seq(1, eg->r_tail);
		bwa_free_read_seq(1, eg->l_tail);
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
