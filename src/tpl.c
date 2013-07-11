/*
 * tpl.c
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
#include "tpl.h"
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

tpl *new_eg() {
	tpl *t = (tpl*) malloc(sizeof(tpl));
	t->ctg = NULL;
	t->l_tail = NULL;
	t->r_tail = NULL;
	t->len = 0;
	t->id = 0;
	t->alive = 1;
	t->is_root = 0;
	t->ori = 0;
	t->comp_id = -1;
	t->start_kmer = 0;
	t->tid = 0;
	t->coverage = 0.0;
	t->kmer_freq = 0;
	t->in_connect = 0;
	t->vertexes = g_ptr_array_sized_new(0);
	return t;
}

void free_readarray(readarray *ra) {
	if (!ra)
		return;
	g_ptr_array_free(ra, TRUE);
}

/**
 * Get virtual tail of an tpl.
 * This is used when a branch tpl is connected to the locus of '^'.
 * When another third-layer branch tpl is connected to this branch tpl,
 * 	its own length may be not long enough, then it cuts some length from the main template.
 * If the length of the branch tpl is long enough, the virtual tail is not used.
 *
 * tpl: 	==============================
 * shift: 	                    ^
 * ori: 	1 (to the left), using the right virtual tail
 * tail_len:	                --------------
 * Return:                      ==========++++
 * '++++' is the partial virtual tail of current tpl
 */
bwa_seq_t *cut_tpl_tail(tpl *t, const int pos, const int tail_len,
		const int ori) {
	bwa_seq_t *tail = NULL, *partial = NULL, *main_tail = NULL;
	int v_tail_len = 0;
	if (t->len < pos)
		return new_seq(t->ctg, t->len, 0);
	if (ori) {
		partial = new_seq(t->ctg, t->len - pos, pos);
		main_tail = t->r_tail;
	} else {
		partial = new_seq(t->ctg, pos, 0);
		main_tail = t->l_tail;
	}
	// If the tpl is long, cut the tail directly
	if (partial->len >= tail_len) {
		if (ori)
			tail = new_seq(partial, tail_len, 0);
		else
			tail = new_seq(partial, tail_len, partial->len - tail_len);
	} else {
		// If the tpl has a virtual tail, try to get it
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

void set_tail(tpl *branch, tpl *parent_eg, const int shift,
		const int tail_len, const int ori) {
	bwa_seq_t *tmp = NULL;
	// The right/left tail would be used in cut_tpl_tail function,
	//	but will be replaced later. So save to tmp first, free the old one later.
	if (ori) {
		tmp = branch->r_tail;
		branch->r_tail = cut_tpl_tail(parent_eg, shift, tail_len, ori);
	} else {
		tmp = branch->l_tail;
		branch->l_tail = cut_tpl_tail(parent_eg, shift, tail_len, ori);
	}
	bwa_free_read_seq(1, tmp);
}

void save_tpls(tplarray *pfd_ctg_ids, FILE *ass_fa, const int ori,
		const int p_all, const int min_len) {
	int i = 0;
	tpl *t;
	char *h;
	bwa_seq_t *contig;
	if (!pfd_ctg_ids || pfd_ctg_ids->len <= 0)
		return;
	h = malloc(BUFSIZE);
	for (i = 0; i < pfd_ctg_ids->len; i++) {
		t = (tpl*) g_ptr_array_index(pfd_ctg_ids, i);
		//show_debug_msg(__func__, "Saving tpl %d length %d, alive %d \n", t->id, t->len, t->alive);
		if (p_all || (t && t->alive && t->ctg && t->len > min_len)) {
			contig = t->ctg;
			if (ori)
				seq_reverse(contig->len, contig->seq, 0);
			sprintf(h, ">%"ID64" length: %d start: %" ID64 "\n", t->id,
					contig->len, t->start_kmer);
			save_con(h, contig, ass_fa);
		}
	}
	free(h);
}

void destroy_eg(tpl *t) {
	if (t) {
		show_debug_msg(__func__, "Freeing tpl [%d, %d] \n", t->id, t->len);
		bwa_free_read_seq(1, t->ctg);
		bwa_free_read_seq(1, t->r_tail);
		bwa_free_read_seq(1, t->l_tail);
		if (t->vertexes)
			g_ptr_array_free(t->vertexes, TRUE);
		free(t);
	}
}

/**
 * Free the tpl contigs first, will be destoryed later.
 */
void free_eg_seq(tpl *t) {
	if (t) {
		free_read_seq(t->ctg);
		t->ctg = NULL;
		free_read_seq(t->r_tail);
		t->r_tail = NULL;
		free_read_seq(t->l_tail);
		t->l_tail = NULL;
	}
}