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
	eg->ctg = NULL;
	eg->l_tail = NULL;
	eg->r_tail = NULL;
	eg->in_egs = NULL;
	eg->out_egs = NULL;
	eg->len = 0;
	eg->id = 0;
	eg->alive = 1;
	eg->is_root = 0;
	eg->ori = 0;
	eg->comp_id = -1;
	eg->start_kmer = 0;
	eg->tid = 0;
	eg->coverage = 0.0;
	eg->kmer_freq = 0;
	return eg;
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
 * ori: 	1 (to the left), using the right virtual tail
 * tail_len:	                --------------
 * Return:                      ==========++++
 * '++++' is the partial virtual tail of current edge
 */
bwa_seq_t *cut_edge_tail(edge *eg, const int shift, const int tail_len,
		const int ori) {
	bwa_seq_t *tail = NULL, *partial = NULL, *main_tail = NULL;
	int v_tail_len = 0;
	if (eg->len < shift)
		return new_seq(eg->ctg, eg->len, 0);
	if (ori) {
		partial = new_seq(eg->ctg, eg->len - shift, shift);
		main_tail = eg->r_tail;
	} else {
		partial = new_seq(eg->ctg, shift, 0);
		main_tail = eg->l_tail;
	}
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
				memcpy(tail->seq, main_tail->seq + (main_tail->len
						+ partial->len - v_tail_len), sizeof(ubyte_t)
						* (v_tail_len - partial->len));
				memcpy(tail->seq + (v_tail_len - partial->len), partial->seq,
						sizeof(ubyte_t) * partial->len);
			}
			tail->len = v_tail_len;
			set_rev_com(tail);
		} else
			tail = new_seq(partial, partial->len, 0);
	}
	bwa_free_read_seq(1, partial);
	return tail;
}

void set_tail(edge *branch, edge *parent_eg, const int shift,
		const int tail_len, const int ori) {
	bwa_seq_t *tmp = NULL;
	// The right/left tail would be used in cut_edge_tail function,
	//	but will be replaced later. So save to tmp first, free the old one later.
	if (ori) {
		tmp = branch->r_tail;
		branch->r_tail = cut_edge_tail(parent_eg, shift, tail_len, ori);
	} else {
		tmp = branch->l_tail;
		branch->l_tail = cut_edge_tail(parent_eg, shift, tail_len, ori);
	}
	bwa_free_read_seq(1, tmp);
}

void save_edges(edgearray *pfd_ctg_ids, FILE *ass_fa, const int ori,
		const int p_all, const int min_len) {
	int i = 0;
	edge *eg;
	char *h;
	bwa_seq_t *contig;
	if (!pfd_ctg_ids || pfd_ctg_ids->len <= 0)
		return;
	h = malloc(BUFSIZE);
	for (i = 0; i < pfd_ctg_ids->len; i++) {
		eg = (edge*) g_ptr_array_index(pfd_ctg_ids, i);
		//show_debug_msg(__func__, "Saving edge %d length %d, alive %d \n", eg->id, eg->len, eg->alive);
		if (p_all || (eg && eg->alive && eg->ctg && eg->len > min_len)) {
			contig = eg->ctg;
			if (ori)
				seq_reverse(contig->len, contig->seq, 0);
			sprintf(h, ">%"ID64" length: %d start: %" ID64 "\n", eg->id,
					contig->len, eg->start_kmer);
			save_con(h, contig, ass_fa);
		}
	}
	free(h);
}

void destroy_eg(edge *eg) {
	if (eg) {
		show_debug_msg(__func__, "Freeing edge [%d, %d] \n", eg->id, eg->len);
		bwa_free_read_seq(1, eg->ctg);
		bwa_free_read_seq(1, eg->r_tail);
		bwa_free_read_seq(1, eg->l_tail);
		if (eg->in_egs)
			g_ptr_array_free(eg->in_egs, TRUE);
		if (eg->out_egs)
			g_ptr_array_free(eg->out_egs, TRUE);
		free(eg);
	}
}

/**
 * Free the edge contigs first, will be destoryed later.
 */
void free_eg_seq(edge *eg) {
	if (eg) {
		free_read_seq(eg->ctg);
		eg->ctg = NULL;
		free_read_seq(eg->r_tail);
		eg->r_tail = NULL;
		free_read_seq(eg->l_tail);
		eg->l_tail = NULL;
		if (eg->in_egs)
			g_ptr_array_free(eg->in_egs, TRUE);
		if (eg->out_egs)
			g_ptr_array_free(eg->out_egs, TRUE);
		eg->in_egs = NULL;
		eg->out_egs = NULL;
	}
}
