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
#include "tpl.hpp"
#include "peseq.h"
#include "bwtaln.h"
#include "k_hash.h"

using namespace std;

eg_gap *init_gap(int s_index, int size, int ori) {
	eg_gap *gap = (eg_gap*) malloc(sizeof(eg_gap));
	gap->s_index = s_index;
	gap->size = size;
	gap->ori = ori;
	return gap;
}

void p_tpl(tpl *t) {
	if (!t) {
		show_debug_msg(__func__, "---- Template is NULL ----\n");
	}
	show_debug_msg(__func__, "---- Template %d ----\n", t->id);
	show_debug_msg(__func__, "\t Length: %d \n", t->len);
	if (t->reads)
		show_debug_msg(__func__, "\t Reads: %d\n", t->reads->len);
	if (t->tried)
		show_debug_msg(__func__, "\t Tried: %d\n", t->tried->len);
	if (t->m_juncs)
		show_debug_msg(__func__, "\t Junctions as main: %d \n", t->m_juncs->len);
	if (t->b_juncs)
		show_debug_msg(__func__, "\t Junctions as branch: %d \n",
				t->b_juncs->len);
	if (t->start_read)
		p_query("start", t->start_read);
	if (t->l_tail)
		p_query("ltail", t->l_tail);
	if (t->r_tail)
		p_query("rtail", t->r_tail);
	if (t->ctg)
		p_ctg_seq("contg", t->ctg);
	show_debug_msg(__func__, "---- Template %d ---- \n", t->id);
}

void free_eg_gap(eg_gap *gap) {
	if (gap)
		free(gap);
}

tpl *new_tpl() {
	tpl *t = (tpl*) malloc(sizeof(tpl));
	t->ctg = NULL;
	t->l_tail = NULL;
	t->r_tail = NULL;
	t->m_juncs = NULL;
	t->b_juncs = NULL;
	t->len = 0;
	t->id = 0;
	t->alive = 1;
	t->is_root = 0;
	t->ori = 0;
	t->start_read = NULL;
	t->tid = 0;
	t->kmer_freq = 0;
	t->in_connect = 0;
	t->vertexes = g_ptr_array_sized_new(0);
	t->reads = g_ptr_array_sized_new(0);
	t->tried = g_ptr_array_sized_new(0);
	return t;
}

/**
 * Get the tail for extension.
 * The only chance that the template is shorter than len is after timmed by junction.
 */
bwa_seq_t *get_tail(tpl *t, int len, const int ori) {
	bwa_seq_t *tail = NULL;
	// Must be something wrong!
	if (!t || t->len < 0 || len < 0)
		return NULL;
	if (t->len >= len) {
		return ori ? new_seq(t->ctg, len, 0) : new_seq(t->ctg, len, t->len
				- len);
	} else {
		return ori ? new_seq(t->start_read, len, 0) : new_seq(t->start_read,
				len, t->start_read->len - len);
	}
}

/**
 * When connecting, need to check whether its reverse-complement or not
 * So need to get the subseq on branch and main to compare base-by-base
 * Its length is read_len - 1.
 */
bwa_seq_t *get_ol_with_connector(tpl *branch, const int read_len, const int ori) {
	bwa_seq_t *ol = NULL;
	if (!branch || !branch->start_read)
		return blank_seq(0);
	if (branch->len < read_len) {
		return ori ? new_seq(branch->start_read, read_len - 1, 0) : new_seq(
				branch->start_read, read_len - 1, 1);
	} else {
		return ori ? new_seq(branch->ctg, read_len - 1, 0) : new_seq(
				branch->ctg, read_len - 1, branch->len - (read_len - 1));
	}
}

gint cmp_tpl_by_id(gpointer a, gpointer b) {
	tpl *c_a = *((tpl**) a);
	tpl *c_b = *((tpl**) b);
	return (c_a->id - c_b->id);
}

/**
 * Remove duplicates in an array
 */
GPtrArray *rm_dup_reads_on_tpl(GPtrArray *reads) {
	int i = 0;
	bwa_seq_t *pre = NULL, *post = NULL;
	GPtrArray *uni_reads = NULL;
	if (!reads || reads->len < 2) {
		return reads;
	}
	uni_reads = g_ptr_array_sized_new(reads->len);
	g_ptr_array_sort(reads, (GCompareFunc) cmp_reads_by_contig_id);
	pre = (bwa_seq_t*) g_ptr_array_index(reads, 0);
	g_ptr_array_add(uni_reads, pre);
	for (i = 1; i < reads->len; i++) {
		post = (bwa_seq_t*) g_ptr_array_index(reads, i);
		if (post->contig_id != pre->contig_id) {
			g_ptr_array_add(uni_reads, post);
			pre = post;
		}
	}
	g_ptr_array_free(reads, TRUE);
	return uni_reads;
}

/**
 * Get the nearest index of read with name as the query
 */
int binary_search_read(GPtrArray *reads, bwa_seq_t *q) {
	int start = 0, end = reads->len, middle = (start + end) / 2;
	bwa_seq_t *r = NULL, *m = NULL;
	//p_readarray(reads, 1);
	//p_query(__func__, q);
	while (start < end - 1) {
		middle = (start + end) / 2;
		//show_debug_msg(__func__, "[%d, %d, %d]\n", start, middle, end);
		m = (bwa_seq_t*) g_ptr_array_index(reads, middle);
		//p_query("MIDDLE", m);
		if (atoi(m->name) - atoi(q->name) == 0) {
			return middle;
		} else {
			if (atoi(m->name) - atoi(q->name) > 0) {
				end = middle;
			} else {
				start = middle + 1;
			}
		}
		//show_debug_msg(__func__, "[%d, %d, %d]\n", start, middle, end);
	}
	return middle;
}

/**
 * Find pairs between two lists
 * If t1_id is not 0, reads on reads_1 must has contig_id == t1_id
 */
int find_pairs(GPtrArray *reads_1, GPtrArray *reads_2, int t1_id, int t2_id,
		int start_2, int end_2, const int min_n_pairs) {
	bwa_seq_t *r1 = NULL, *r2 = NULL;
	uint32_t i = 0, j = 0, index_1 = 0, index_2 = 0;
	int n_pairs = 0;
	if (!reads_1 || !reads_2 || reads_1->len == 0 || reads_2->len == 0)
		return 0;
	r1 = (bwa_seq_t*) g_ptr_array_index(reads_2, 0);
	//show_debug_msg(__func__, "Searching %s in reads_1 \n", r1->name);
	index_1 = binary_search_read(reads_1, r1);
	//show_debug_msg(__func__, "Index_1: %d\n", index_1);

	r2 = (bwa_seq_t*) g_ptr_array_index(reads_1, 0);
	//show_debug_msg(__func__, "Searching %s in reads_2 \n", r2->name);
	index_2 = binary_search_read(reads_2, r2);
	//show_debug_msg(__func__, "Index_1: %d\n", index_2);

	//show_debug_msg(__func__, "Reads 1 \n");
	//p_readarray(reads_1, 1);
	//show_debug_msg(__func__, "Reads 2 \n");
	//p_readarray(reads_2, 1);

	for (i = index_1; i < reads_1->len; i++) {
		r1 = (bwa_seq_t*) g_ptr_array_index(reads_1, i);
		// Only check template id when the id is not 0
		if (t1_id && r1->contig_id != t1_id)
			continue;
		if (index_2 >= reads_2->len - 1)
			break;
		//p_query("READ_1", r1);
		for (j = index_2; j < reads_2->len; j++) {
			r2 = (bwa_seq_t*) g_ptr_array_index(reads_2, j);
			if (t2_id && (r2->contig_id != t2_id))
				continue;
			if (r2->contig_locus < start_2 || r2->contig_locus > end_2)
				continue;
			//p_query("READ_2", r2);
			if (is_mates(r1->name, r2->name)) {
				//p_query(__func__, r1);
				//p_query(__func__, r2);
				// Ignore if either one of them is like: "AAATAAAAAA"
				if (is_biased_q(r1) || is_biased_q(r2))
					continue;
				n_pairs++;
			}
			if (atoi(r2->name) - atoi(r1->name) >= 1) {
				index_2 = j;
				break;
			}
		}
		if (n_pairs >= min_n_pairs)
			return 1;
	}
	return 0;
}

/**
 * Find whether the two templates contain at least min_n_pairs pairs of mates.
 * Assume the reads on them are already sorted by name
 */
int vld_tpl_mates(tpl *t1, tpl *t2, int start_2, int end_2,
		const int min_n_pairs) {
	int is_valid = 0;
	is_valid = find_pairs(t1->reads, t2->reads, t1->id, t2->id, start_2, end_2,
			min_n_pairs);
	return is_valid;
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

void set_tail(tpl *branch, tpl *parent_eg, const int shift, const int tail_len,
		const int ori) {
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

/**
 * Update the contig_locus value after left extension
 */
void upd_locus_on_tpl(tpl *t, int pre_t_len, int pre_n_reads) {
	index64 i = 0;
	bwa_seq_t *r = NULL;

	if (!t || !t->reads || pre_n_reads < 0 || pre_t_len < 0)
		return;
	for (i = 0; i < pre_n_reads; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(t->reads, i);
		r->contig_locus += t->len - pre_t_len;
	}
	for (i = pre_n_reads; i < t->reads->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(t->reads, i);
		r->contig_locus = t->len - r->contig_locus;
	}
}

/**
 * Add a read to the template
 */
void add2tpl(tpl *t, bwa_seq_t *r, const int locus) {
	r->contig_id = t->id;
	r->contig_locus = locus;
	r->pos = -1;
	r->status = USED;
	g_ptr_array_add(t->reads, r);
}

/**
 * Reset the reads on the template
 */
void refresh_tpl_reads(hash_table *ht, tpl *t, int mismatches) {
	bwa_seq_t *r = NULL, *seq = NULL, *window = NULL;
	int left_len = 0, counted_len = 0, right_len = 0;
	int i = 0, j = 0;
	GPtrArray *refresh = NULL, *hits = NULL;
	if (!t || !t->reads || t->reads->len <= 0 || t->len < 0)
		return;
	unfrozen_tried(t);
	refresh = g_ptr_array_sized_new(t->reads->len);
	unhold_reads_array(t->reads);
	t->reads = refresh;

	seq = get_tpl_ctg_wt(t, &left_len, &right_len, &counted_len);

	// If it happens, means something wrong
	if (seq->len < ht->o->read_len) {
		show_debug_msg(
				"[WARNING]",
				"The sequence with tails shorter than read length: [%d, %d] \n",
				t->id, t->len);
		p_tpl(t);
		bwa_free_read_seq(1, seq);
		return;
	}

	p_tpl(t);
	for (i = 0; i <= seq->len - ht->o->read_len; i++) {
		window = new_seq(seq, ht->o->read_len, i);
		hits = g_ptr_array_sized_new(4);
		hits = find_both_fr_full_reads(ht, window, hits, mismatches);
		for (j = 0; j < hits->len; j++) {
			r = (bwa_seq_t*) g_ptr_array_index(hits, j);
			// For reads partially on left tail, the locus is negative
			if (r->status == FRESH)
				add2tpl(t, r, i - left_len);
		}
		g_ptr_array_free(hits, TRUE);
		bwa_free_read_seq(1, window);
	}
	bwa_free_read_seq(1, seq);
	g_ptr_array_sort(t->reads, (GCompareFunc) cmp_reads_by_name);
}

/**
 * After connecting to existing, need to set the junction reads to USED,
 * 	to avoid multiple junctions
 */
void refresh_reads_on_tail(hash_table *ht, tpl *t, int mismatches) {
	bwa_seq_t *tail = NULL, *s = NULL, *window = NULL;
	bwa_seq_t *r = NULL;
	int i = 0, j = 0, len = 0, borrow_len = 0;
	GPtrArray *hits = NULL;
	if (t->l_tail) {
		borrow_len = (t->len >= (ht->o->read_len - 1)) ? (ht->o->read_len)
				: t->len;
		len = t->l_tail->len + borrow_len;
		s = blank_seq(len);
		memcpy(s->seq, t->l_tail->seq, t->l_tail->len * sizeof(ubyte_t));
		s->len += t->l_tail->len;
		memcpy(s->seq + s->len, t->ctg->seq, borrow_len * sizeof(ubyte_t));
		s->len += borrow_len;

		for (i = 0; i <= s->len - ht->o->read_len; i++) {
			window = new_seq(s, ht->o->read_len, i);
			hits = g_ptr_array_sized_new(4);
			hits = find_both_fr_full_reads(ht, window, hits, mismatches);
			for (j = 0; j < hits->len; j++) {
				r = (bwa_seq_t*) g_ptr_array_index(hits, j);
				// For reads partially on left tail, the locus is negative
				if (r->status == FRESH)
					add2tpl(t, r, i - t->l_tail->len);
			}
			g_ptr_array_free(hits, TRUE);
			bwa_free_read_seq(1, window);
		}
		bwa_free_read_seq(1, s);
	}

	if (t->r_tail) {
		borrow_len = (t->len >= (ht->o->read_len - 1)) ? (ht->o->read_len)
				: t->len;
		len = t->r_tail->len + borrow_len;
		s = blank_seq(len);
		memcpy(s->seq, t->ctg->seq + (t->len - borrow_len), borrow_len
				* sizeof(ubyte_t));
		s->len = borrow_len;
		memcpy(s->seq + s->len, t->r_tail->seq, t->r_tail->len
				* sizeof(ubyte_t));
		s->len += t->r_tail->len;

		for (i = 0; i <= s->len - ht->o->read_len; i++) {
			window = new_seq(s, ht->o->read_len, i);
			hits = g_ptr_array_sized_new(4);
			hits = find_both_fr_full_reads(ht, window, hits, mismatches);
			for (j = 0; j < hits->len; j++) {
				r = (bwa_seq_t*) g_ptr_array_index(hits, j);
				// For reads partially on left tail, the locus is negative
				add2tpl(t, r, t->len - borrow_len + i);
			}
			g_ptr_array_free(hits, TRUE);
			bwa_free_read_seq(1, window);
		}
		bwa_free_read_seq(1, s);
	}
}

/**
 * Remove a read from the pool and reset the the read status
 */
void rm_from_tpl(tpl *t, int index) {
	if (index < 0 || !t || !t->reads || index >= t->reads->len)
		return;
	bwa_seq_t * r = (bwa_seq_t*) g_ptr_array_index(t->reads, index);
	reset_to_fresh(r);
	g_ptr_array_remove_index_fast(t->reads, index);
}

/**
 * Add to TRIED, pretend from using for this template
 */
void add2tried(tpl *t, bwa_seq_t *r) {
	r->contig_id = t->id;
	r->contig_locus = -1;
	r->status = TRIED;
	r->pos = -1;
	g_ptr_array_add(t->tried, r);
}

void reset_reads_to_fresh(GPtrArray *reads) {
	index64 i = 0;
	bwa_seq_t *r = NULL;
	for (i = 0; i < reads->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(reads, i);
		reset_to_fresh(r);
	}
}

void unfrozen_tried(tpl *t) {
	reset_reads_to_fresh(t->tried);
	while (t->tried->len > 0)
		g_ptr_array_remove_index_fast(t->tried, 0);
}

/**
 * If some read is already used some template, do not reset
 */
void unhold_reads_array(GPtrArray *reads) {
	int i = 0;
	bwa_seq_t *r = NULL;
	for (i = 0; i < reads->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(reads, i);
		if (r->status == FRESH)
			reset_to_fresh(r);
	}
	g_ptr_array_free(reads, TRUE);
}

/**
 * Read extension starts from some read, mark this read and its close reads as used first
 */
void mark_init_reads_used(hash_table *ht, tpl *t, bwa_seq_t *read,
		int mismatches) {
	GPtrArray *hits = NULL;
	bwa_seq_t *r = NULL;
	index64 i = 0;

	hits = find_both_fr_full_reads(ht, read, hits, mismatches);
	for (i = 0; i < hits->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(hits, i);
		//p_query(__func__, r);
		if (r->status == FRESH) {
			add2tpl(t, r, 0);
		}
	}
	//p_readarray(t->reads, 1);
	g_ptr_array_free(hits, TRUE);
}

/**
 * Align the tail of template to find fresh reads
 */
GPtrArray *align_tpl_tail(hash_table *ht, tpl *t, bwa_seq_t *tail,
		int mismatches, int8_t status, int ori) {
	int64_t i = 0, cursor = 0, ol = 0;
	int n_mis = 0, added = 0;
	bwa_seq_t *r = NULL, *tpl_seq = NULL;
	GPtrArray *hits = align_query(ht, tail, status, mismatches);
	GPtrArray *fresh_reads = g_ptr_array_sized_new(hits->len);

	// Smaller than read length happens when it is trimmed to connect at the right/left side
	tpl_seq = (t->len < ht->o->read_len) ? t->start_read : t->ctg;
	if (t->len < ht->o->read_len) {
		tpl_seq = blank_seq(ht->o->read_len);
		memcpy(tpl_seq->seq, t->ctg->seq, t->len * sizeof(ubyte_t));
		tpl_seq->len += t->len;
		memcpy(tpl_seq->seq + tpl_seq->len, t->start_read->seq,
				(ht->o->read_len - tpl_seq->len) * sizeof(ubyte_t));
		tpl_seq->len = ht->o->read_len;
		set_rev_com(tpl_seq);
	} else {
		tpl_seq = t->ctg;
	}

	// These reads are not duplicated
	for (i = 0; i < hits->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(hits, i);
		added = 0;
		// pos is the kmer position on the read
		cursor = ori ? (r->pos - 1) : (r->pos + tail->len);

		//p_query(__func__, r);
		//show_debug_msg(__func__, "CURSOR: %d\n", cursor);

		if (cursor >= 0 && cursor <= r->len - 1) {
			ol = ori ? (r->len - cursor - 1) : cursor;
			if (ol >= tail->len) {
				if (ori) {
					n_mis = seq_ol(r, tpl_seq, ol, mismatches);
				} else {
					n_mis = seq_ol(tpl_seq, r, ol, mismatches);
				}
				// if (!ori) {
				//p_query(__func__, r);
				//p_ctg_seq(__func__, t->ctg);
				//show_debug_msg(__func__, "Should have overlap: %d\n", ol);
				//show_debug_msg(__func__, "N_MISMATCHES: %d\n", n_mis);
				//	}
				// n_mis >= 0 means similar with n_mis mismatches; -1 means not similar
				if (n_mis >= 0) {
					//show_debug_msg(__func__, "Cursor: %d\n", cursor);
					r->cursor = cursor;
					// In the pool, 'pos' stores how many mismatches between read and template
					r->pos = n_mis;
					//p_query("ADDED", r);
					g_ptr_array_add(fresh_reads, r);
					added = 1;
				}
			}
		}
		// If not added, reset the pos to be -1
		if (!added) {
			r->pos = -1;
			r->rev_com = 0;
		}
	}
	//show_debug_msg(__func__, "Reads with the tail: %d\n", fresh_reads->len);
	g_ptr_array_free(hits, TRUE);
	if (t->len < ht->o->read_len) {
		bwa_free_read_seq(1, tpl_seq);
	}
	return fresh_reads;
}

void save_tpls(tplarray *pfd_ctg_ids, FILE *ass_fa, const int ori,
		const int p_all, const int min_len) {
	int i = 0;
	tpl *t;
	char *h;
	bwa_seq_t *contig;
	if (!pfd_ctg_ids || pfd_ctg_ids->len <= 0)
		return;
	h = (char*) malloc(BUFSIZE);
	for (i = 0; i < pfd_ctg_ids->len; i++) {
		t = (tpl*) g_ptr_array_index(pfd_ctg_ids, i);
		//show_debug_msg(__func__, "Saving tpl %d length %d, alive %d \n", t->id, t->len, t->alive);
		if (p_all || (t && t->alive && t->ctg && t->len > min_len)) {
			contig = t->ctg;
			if (ori)
				seq_reverse(contig->len, contig->seq, 0);
			sprintf(h, ">%d length: %d start: %s\n", t->id, contig->len,
					t->start_read->name);
			save_con(h, contig, ass_fa);
		}
	}
	free(h);
}

void destroy_tpl(tpl *t) {
	index64 i = 0;
	bwa_seq_t *r = NULL;
	if (t) {
		show_debug_msg(__func__, "Freeing tpl [%d, %d] \n", t->id, t->len);
		bwa_free_read_seq(1, t->ctg);
		bwa_free_read_seq(1, t->r_tail);
		bwa_free_read_seq(1, t->l_tail);
		if (t->vertexes)
			g_ptr_array_free(t->vertexes, TRUE);
		if (t->reads) {
			for (i = 0; i < t->reads->len; i++) {
				r = (bwa_seq_t*) g_ptr_array_index(t->reads, i);
				r->contig_id = -1;
				r->contig_locus = -1;
				r->pos = -1;
				r->status = TRIED;
			}
			g_ptr_array_free(t->reads, TRUE);
		}

		if (t->tried) {
			for (i = 0; i < t->tried->len; i++) {
				r = (bwa_seq_t*) g_ptr_array_index(t->tried, i);
				r->contig_id = -1;
				r->contig_locus = -1;
				r->pos = -1;
				r->status = TRIED;
			}
			g_ptr_array_free(t->tried, TRUE);
		}

		if (t->m_juncs)
			g_ptr_array_free(t->m_juncs, TRUE);
		if (t->b_juncs)
			g_ptr_array_free(t->b_juncs, TRUE);
		free(t);
	}
}

/**
 * Get the sequence with left/right tails.
 */
bwa_seq_t *get_tpl_ctg_wt(tpl *t, int *l_len, int *r_len, int *t_len) {
	bwa_seq_t *s = NULL;
	int counted_len = 0;
	*l_len = (t->l_tail) ? (t->l_tail->len) : 0;
	*r_len = (t->r_tail) ? (t->r_tail->len) : 0;
	*t_len = *l_len + t->len + *r_len;
	s = blank_seq(*t_len);
	if (t->l_tail) {
		memcpy(s->seq, t->l_tail->seq, t->l_tail->len * sizeof(ubyte_t));
		s->len = *l_len;
	}
	memcpy(s->seq + s->len, t->ctg->seq, t->len);
	s->len += t->len;
	if (t->r_tail) {
		memcpy(s->seq + s->len, t->r_tail->seq, t->r_tail->len
				* sizeof(ubyte_t));
		s->len += t->r_tail->len;
	}
	set_rev_com(s);
	return s;
}
