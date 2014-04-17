/*
 * merge.cpp
 *
 *  Created on: Jul 29, 2013
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
#include "junction.hpp"
#include "merge.hpp"
#include "k_hash.h"

void mv_reads_bt_tpls(tpl *from, tpl *to, int ol, int ori) {
	bwa_seq_t *r = NULL;
	int i = 0, new_locus = 0;
	if (ori) {
		// 'from' is at left, 'to' is at right
		for (i = 0; i < to->reads->len; i++) {
			r = (bwa_seq_t*) g_ptr_array_index(to->reads, i);
			if (r->contig_locus < ol) {
				p_query("DEAD", r);
				reset_to_dead(r);
				g_ptr_array_remove_index_fast(to->reads, i--);
			} else
				r->contig_locus += from->len - ol;
		}
		for (i = 0; i < from->reads->len; i++) {
			r = (bwa_seq_t*) g_ptr_array_index(from->reads, i);
			// Reset to fresh if out of range; otherwise, template correction problematic
			add2tpl(to, r, r->contig_locus);
		}
	} else {
		// 'from' is at right, 'to' is at left
		for (i = 0; i < to->reads->len; i++) {
			r = (bwa_seq_t*) g_ptr_array_index(to->reads, i);
			if (r->contig_locus + r->len >= to->len - ol) {
				p_query("DEAD", r);
				reset_to_dead(r);
				g_ptr_array_remove_index_fast(to->reads, i--);
			}
		}
		for (i = 0; i < from->reads->len; i++) {
			r = (bwa_seq_t*) g_ptr_array_index(from->reads, i);
			if (r->status == USED && r->contig_id == from->id) {
				new_locus = r->contig_locus + to->len - ol;
				add2tpl(to, r, new_locus);
			}
		}
	}
	while (from->reads->len > 0)
		g_ptr_array_remove_index_fast(from->reads, 0);
}

/**
 * Merge jumped template to the main template;
 * Main template is at the left
 */
void merge_tpl_to_right(tpl *jumped, tpl *t, int ol, int rev_com) {
	int i = 0, l_len = 0, new_locus = 0;
	bwa_seq_t *r = NULL;
	ubyte_t c = 0;
	junction *jun = 0;
	show_debug_msg(__func__, "Merging Right: [%d, %d] Vs. [%d, %d]; ol: %d\n",
			jumped->id, jumped->len, t->id, t->len, ol);
	if (!t || !jumped || ol >= t->len || ol >= jumped->len)
		return;
	if (t->l_tail) {
		show_debug_msg(
				__func__,
				"[WARNING]: jumped template [%d, %d] cannot be merged to [%d, %d] \n",
				jumped->id, jumped->len, t->id, t->len);
		return;
	}

	l_len = t->len;
	// Copy the template sequence from right to left
	for (i = jumped->len - ol - 1; i >= 0; i--) {
		c = rev_com ? jumped->ctg->rseq[i] : jumped->ctg->seq[i];
		ext_con(t->ctg, c, 1);
	}
	t->len = t->ctg->len;
	set_rev_com(t->ctg);

	// Update the read locus on the main template
	for (i = 0; i < t->reads->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(t->reads, i);
		r->contig_locus += jumped->len - ol;
	}

	if (t->m_juncs) {
		for (i = 0; i < t->m_juncs->len; i++) {
			jun = (junction*) g_ptr_array_index(t->m_juncs, i);
			jun->locus += jumped->len - ol;
		}
	}

	// Add the reads on jumped template to the main template
	for (i = 0; i < jumped->reads->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(jumped->reads, i);
		if (r->status == USED && r->contig_id == jumped->id) {
			new_locus = r->contig_locus;
			if (rev_com) {
				r->rev_com = r->rev_com ? 0 : 1;
				new_locus = jumped->len - r->contig_locus - r->len;
			}
			add2tpl(t, r, new_locus);
		}
	}
	while (jumped->reads->len > 0)
		g_ptr_array_remove_index_fast(jumped->reads, 0);
	jumped->alive = 0;
}

/**
 * Merge jumped template to the main template;
 * Main template is at the left
 */
void merge_tpl_to_left(tpl *t, tpl *jumped, int ol, int rev_com) {
	int i = 0, l_len = 0, new_locus = 0;
	bwa_seq_t *r = NULL;
	ubyte_t c = 0;
	if (!t || !jumped || ol >= t->len || ol >= jumped->len)
		return;
	if (t->r_tail) {
		show_debug_msg(
				__func__,
				"[WARNING]: jumped template [%d, %d] cannot be merged to [%d, %d] \n",
				jumped->id, jumped->len, t->id, t->len);
		return;
	}
	l_len = t->len;
	show_debug_msg(__func__, "Merging Left: [%d, %d] Vs. [%d, %d]; ol: %d\n",
			t->id, t->len, jumped->id, jumped->len, ol);
	// Copy the template sequence from right to left
	for (i = ol; i < jumped->len; i++) {
		c = rev_com ? jumped->ctg->rseq[i] : jumped->ctg->seq[i];
		ext_con(t->ctg, c, 0);
		//p_ctg_seq(__func__, t->ctg);
	}
	t->len = t->ctg->len;
	set_rev_com(t->ctg);

	// Add the reads on right template to the left template
	for (i = 0; i < jumped->reads->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(jumped->reads, i);
		if (r->status == USED && r->contig_id == jumped->id) {
			// Update the locus
			new_locus = rev_com ? (jumped->len - r->contig_locus - r->len
					+ l_len - ol) : (r->contig_locus + l_len - ol);
			// If need to reverse complement the right template,
			//	need to reverse the reverse complement flag on the reads
			if (rev_com) {
				r->rev_com = r->rev_com ? 0 : 1;
			}
			add2tpl(t, r, new_locus);
		}
	}
	while (jumped->reads->len > 0)
		g_ptr_array_remove_index_fast(jumped->reads, 0);
	jumped->alive = 0;
}

int similar_to_merge(bwa_seq_t *from, bwa_seq_t *jumped, int unique_len) {
	int from_s = 0, from_e = 0, jumped_s = 0, jumped_e = 0;
	int score = 0, similar = 0;
	int min_score = unique_len;

	score = smith_waterman_simple(from, jumped, &from_s, &from_e, &jumped_s,
			&jumped_e, min_score);
	if (score > min_score && score + 4 > from_e - from_e)
		return 1;
}

/**
 * Before merging, check whether there are pairs at the head/tail area.
 */
int paired_at_head_tail(bwa_seq_t *seqs, tpl *left, tpl *right) {
	int range = INS_SIZE - 2 * SD_INS_SIZE;
	GPtrArray *tail_reads = get_supporting_reads(left, left->len - range, left->len);
	int n = 0, i = 0, total = 0;
	bwa_seq_t *m = NULL, *r = NULL;
	for (i = 0; i < tail_reads->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(tail_reads, i);
		m = get_mate(r, seqs);
		if (m->status == USED && m->contig_id == right->id && m->contig_locus < range)
			n++;
	}
	g_ptr_array_free(tail_reads, TRUE);
	if (n >= 1)
		return 1;
}

/**
 * Try to merge two templates.
 * The 'jumped' template does not have any junctions on it
 */
int merged_jumped(hash_table *ht, tpl *from, tpl *jumped,
		bwa_seq_t *jumping_read, int mis) {
	int rev_com = 0, n_mis = 0, in_paired = 0, spanning = 0;
	int from_s = 0, from_e = 0, jumped_s = 0, jumped_e = 0;
	int score = 0, similar = 0, ori_len = 0, side = 0;
	int max_ol = ht->o->read_len * 1.5;
	float from_cov = 0.0, jumped_cov = 0.0;
	bwa_seq_t *from_seq = NULL, *jumped_seq = NULL, *r = NULL;
	//	p_tpl_reads(t);
	//	p_tpl_reads(jumped);
	if (!paired_by_reads(ht->seqs, from, jumped, 1) && jumped->len
			> ht->o->read_len + 2)
		return 0;

	side = should_at_which_side(ht->seqs, from, jumping_read);
	show_debug_msg(__func__, "Read %s at side %d; %s \n", jumping_read->name,
			side, side == 0 ? "LEFT" : "RIGHT");

	if (!from->r_tail && (side == RIGHT_SIDE || side == UNKNOWN_SIDE)) {
		from_seq = new_seq(from->ctg, min(max_ol, from->len), from->len
				- min(max_ol, from->len));
		jumped_seq = new_seq(jumped->ctg, min(jumped->len, max_ol), 0);
		score = smith_waterman_simple(from_seq, jumped_seq, &from_s, &from_e,
				&jumped_s, &jumped_e, ht->o->k);
		in_paired = paired_at_head_tail(ht->seqs, from, jumped);
		/**
		p_ctg_seq("FROM", from->ctg);
		p_ctg_seq("FROM_PART END", from_seq);
		p_ctg_seq("JUMPED", jumped->ctg);
		p_ctg_seq("JUMPED_PART HEAD", jumped_seq);
		show_debug_msg(__func__, "SCORE: %d \n", score);
		show_debug_msg(__func__, "FROM: [%d, %d] \n", from_s, from_e);
		show_debug_msg(__func__, "JUMPED: [%d, %d] \n", jumped_s, jumped_e);
		//p_tpl_reads(jumped);
		//**/

		if (score >= ht->o->k - 1 && in_paired) {
			spanning = pairs_spanning_locus(ht->seqs, from, from_e);
			if (spanning < 2) {
				// For the overlapped region, pick the one with higher coverage
				from_cov = calc_tpl_cov(from, from->len - from_seq->len + from_s,
						from->len - from_seq->len + from_e, ht->o->read_len);
				jumped_cov = calc_tpl_cov(jumped, jumped_s, jumped_e,
						ht->o->read_len);
				if (from_cov > jumped_cov) {
					truncate_tpl(from, from_seq->len - from_e, 0);
					truncate_tpl(jumped, jumped_e, 1);
				} else {
					truncate_tpl(from, from_seq->len - from_s, 0);
					truncate_tpl(jumped, jumped_s, 1);
				}
				bwa_free_read_seq(1, from_seq);
				bwa_free_read_seq(1, jumped_seq);
				ori_len = from->len;
				merge_tpl_to_left(from, jumped, 0, rev_com);
				refresh_tpl_reads(ht, from, N_MISMATCHES);
				correct_tpl_base(ht->seqs, from, ht->o->read_len, ori_len
						- ht->o->read_len, ori_len + ht->o->read_len);
				return 1;
			}
		}
		bwa_free_read_seq(1, from_seq);
		bwa_free_read_seq(1, jumped_seq);
	}

	if (!from->l_tail && (side == UNKNOWN_SIDE || side == LEFT_SIDE)) {
		// From right template jump to left
		from_seq = new_seq(from->ctg, min(max_ol, from->len), 0);
		jumped_seq = new_seq(jumped->ctg, min(jumped->len, max_ol), jumped->len
				- min(jumped->len, max_ol));
		score = smith_waterman_simple(from_seq, jumped_seq, &from_s, &from_e,
				&jumped_s, &jumped_e, ht->o->k);
		in_paired = paired_at_head_tail(ht->seqs, jumped, from);
		/**
		p_ctg_seq("FROM", from->ctg);
		p_ctg_seq("FROM_PART END", from_seq);
		p_ctg_seq("JUMPED", jumped->ctg);
		p_ctg_seq("JUMPED_PART HEAD", jumped_seq);
		show_debug_msg(__func__, "Score: %d \n", score);
		show_debug_msg(__func__, "FROM: [%d, %d] \n", from_s, from_e);
		show_debug_msg(__func__, "JUMPED: [%d, %d] \n", jumped_s, jumped_e);
		// **/
		if (score >= ht->o->k - 1 && in_paired) {
			spanning = pairs_spanning_locus(ht->seqs, from, from_s);
			if (spanning < 2) {
				from_cov = calc_tpl_cov(from, from_s, from_e, ht->o->read_len);
				jumped_cov = calc_tpl_cov(jumped, jumped->len - jumped_seq->len
						+ jumped_s, jumped->len - jumped_seq->len + jumped_e,
						ht->o->read_len);
				if (from_cov > jumped_cov) {
					truncate_tpl(from, from_s, 1);
					truncate_tpl(jumped, jumped_seq->len - jumped_s, 0);
				} else {
					truncate_tpl(from, from_e, 1);
					truncate_tpl(jumped, jumped_seq->len - jumped_e, 0);
				}
				bwa_free_read_seq(1, from_seq);
				bwa_free_read_seq(1, jumped_seq);
				ori_len = jumped->len;
				merge_tpl_to_right(jumped, from, 0, rev_com);
				refresh_tpl_reads(ht, from, N_MISMATCHES);
				correct_tpl_base(ht->seqs, from, ht->o->read_len, ori_len
						- ht->o->read_len, ori_len + ht->o->read_len);
				//p_ctg_seq("MERGED", t->ctg);
				return 1;
			}

		}
		bwa_free_read_seq(1, from_seq);
		bwa_free_read_seq(1, jumped_seq);
	}
	return 0;
}

/**
 * Merge the right template to the left template
 * The right template is not destroyed (its attributes are not altered), it is simply marked as 'DEAD'.
 * The right template is destroyed at hash_to_array in ass.cpp.
 * If param rev_com is 1, will merge the reverse complement of right template to left.
 */
int merge_tpls(tpl *left, tpl *right, int ol, int rev_com) {
	bwa_seq_t *r = NULL;
	int i = 0, connected = 0;
	int new_locus = 0, l_len = 0;
	junction *jun = NULL;
	ubyte_t c = 0;
	if (!left->alive || !right->alive || left->r_tail || right->l_tail) {
		show_debug_msg(
				__func__,
				"[WARNING] Merging templates not alive or with wrong ends: [%d, %d] and [%d, %d] \n",
				left->id, left->len, right->id, right->len);
		p_tpl(left);
		p_tpl(right);
		return 0;
	}

	// If they have junctions, just ignore
	if (tpls_have_junction(left, right)) {
		show_debug_msg(__func__,
				"[WARNING] These two templates have some junction, cannot be merged. \n");
		show_debug_msg(__func__, "Left  template: [%d, %d] \n", left->id,
				left->len);
		show_debug_msg(__func__, "RIght template: [%d, %d] \n", right->id,
				right->len);
		return 0;
	}

	// If trying to reverse complement it, there cannot be any connection on right
	if (rev_com) {
		if ((right->m_juncs && right->m_juncs->len > 0) || (right->b_juncs
				&& right->b_juncs->len > 0)) {
			show_debug_msg(
					__func__,
					"[WARNING] The right template is connected, cannot be reverse complement merged. \n");
			show_debug_msg(__func__, "Left  template: [%d, %d] \n", left->id,
					left->len);
			show_debug_msg(__func__, "RIght template: [%d, %d] \n", right->id,
					right->len);
			return 0;
		}
	}

	show_debug_msg(
			__func__,
			"Merging templates [%d, %d] and [%d, %d] with reverse complement: %d \n",
			left->id, left->len, right->id, right->len, rev_com);

	//p_tpl(left);
	//p_tpl(right);

	l_len = left->len;
	// Copy the template sequence from right to left
	for (i = ol; i < right->len; i++) {
		c = rev_com ? right->ctg->rseq[i] : right->ctg->seq[i];
		ext_con(left->ctg, c, 0);
	}
	left->len = left->ctg->len;
	set_rev_com(left->ctg);

	// Copy the tail of right to left
	if (right->r_tail) {
		left->r_tail = new_seq(right->r_tail, right->r_tail->len, 0);
	}

	// Add the reads on right template to the left template
	for (i = 0; i < right->reads->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(right->reads, i);
		if (r->status == USED && r->contig_id == right->id) {
			// Update the locus
			new_locus = rev_com ? (right->len - r->contig_locus - r->len
					+ l_len - ol) : (r->contig_locus + l_len - ol);
			// If need to reverse complement the right template,
			//	need to reverse the reverse complement flag on the reads
			if (rev_com) {
				r->rev_com = r->rev_com ? 0 : 1;
			}
			add2tpl(left, r, new_locus);
		}
	}
	while (right->reads->len > 0)
		g_ptr_array_remove_index_fast(right->reads, 0);

	// The junctions on the right now go to the left template
	if (right->b_juncs) {
		for (i = 0; i < right->b_juncs->len; i++) {
			jun = (junction*) g_ptr_array_index(right->b_juncs, i);
			if (jun->branch_tpl == right) {
				jun->branch_tpl = left;
				// If right template connects to somewhere on the right,
				//  now the left templates connects to it.
				if (!left->b_juncs)
					left->b_juncs = g_ptr_array_sized_new(2);
				g_ptr_array_add(left->b_juncs, jun);
			}
		}
	}
	if (right->m_juncs) {
		for (i = 0; i < right->m_juncs->len; i++) {
			jun = (junction*) g_ptr_array_index(right->m_juncs, i);
			if (jun->main_tpl == right) {
				jun->main_tpl = left;
				// Update the junction locus
				jun->locus += l_len;
				if (!left->m_juncs)
					left->m_juncs = g_ptr_array_sized_new(2);
				g_ptr_array_add(left->m_juncs, jun);
			}
		}
	}

	//p_tpl(left);
	right->alive = 0;
	return 1;
}

