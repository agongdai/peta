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
	//p_tpl(t);

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
	//p_tpl(t);
}

/**
 * Try to merge two templates.
 * The 'jumped' template does not have any junctions on it
 */
int merged_jumped(hash_table *ht, tpl *from, tpl *jumped, int mis) {
	int from_s = 0, from_e = 0, jumped_s = 0, jumped_e = 0;
	int score = 0, similar = 0, len = 0, i = 0;
	int max_ol = ht->o->k;
	float from_cov = 0.0, jumped_cov = 0.0;
	bwa_seq_t *from_seq = NULL, *jumped_seq = NULL, *r = NULL;
	junction *jun = NULL;

	p_tpl(from);
	p_tpl(jumped);
	from_seq = new_seq(from->ctg, min(max_ol, from->len), from->len
			- min(max_ol, from->len));
	jumped_seq = new_seq(jumped->ctg, min(jumped->len, max_ol), 0);
	score = smith_waterman_simple(from_seq, jumped_seq, &from_s, &from_e,
			&jumped_s, &jumped_e);
	len = from_seq->len;
	///**
	printf("\n");
	p_ctg_seq("FROM", from->ctg);
	p_ctg_seq("FROM_PART END", from_seq);
	p_ctg_seq("JUMPED", jumped->ctg);
	p_ctg_seq("JUMPED_PART HEAD", jumped_seq);
	show_debug_msg(__func__, "SCORE: %d \n", score);
	show_debug_msg(__func__, "FROM: [%d, %d] \n", from_s, from_e);
	show_debug_msg(__func__, "JUMPED: [%d, %d] \n", jumped_s, jumped_e);

	bwa_free_read_seq(1, from_seq);
	bwa_free_read_seq(1, jumped_seq);

	//p_tpl_reads(jumped);
	//**/
	if (score >= 4) {
		// For the overlapped region, pick the one with higher coverage
		from_cov = calc_tpl_cov(from, from->len - len + from_s,
				from->len - len + from_e, ht->o->read_len);
		jumped_cov = calc_tpl_cov(jumped, jumped_s, jumped_e,
				ht->o->read_len);
		if (from_cov > jumped_cov) {
			truncate_tpl(from, len - from_e, 1, 0);
			truncate_tpl(jumped, jumped_e, 1, 1);
		} else {
			truncate_tpl(from, len - from_s, 1, 0);
			truncate_tpl(jumped, jumped_s, 1, 1);
		}
		len = from->len;
		merge_tpl_to_left(from, jumped, 0, 0);
		refresh_tpl_reads(ht, from, 0, from->len, N_MISMATCHES);
		correct_tpl_base(ht->seqs, from, ht->o->read_len, len
				- ht->o->read_len, len + ht->o->read_len);
		return 1;
	}
	printf("\n");
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
		g_ptr_array_free(right->b_juncs, TRUE);
		right->b_juncs = NULL;
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
		g_ptr_array_free(right->m_juncs, TRUE);
		right->m_juncs = NULL;
	}

	//p_tpl(left);
	right->alive = 0;
	return 1;
}

/**
 * Check whether the unpaired reads on the 'right' template fall on template 'left'
 */
int left_tpl_is_paired(bwa_seq_t *seqs, tpl *left, tpl *right, float pair_pc) {
	float pc = 0.0, n = 0.0, n_right = 0.0;
	bwa_seq_t *r = NULL, *m = NULL;
	int i = 0, dist = 0;
	int max_range = INS_SIZE + GRACE_TIMES * SD_INS_SIZE;
	if (!left->alive || !right->alive) return 0;
	for (i = 0; i < right->reads->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(right->reads, i);
		m = get_mate(r, seqs);
		if (r->rev_com != m->rev_com) continue;
		if (read_on_tpl(left, m)) {
			dist = left->len - m->contig_locus + r->contig_locus;
			if (good_insert_size(dist)) n_right++;
		} else n++;
	}
	pc = n_right / n;
	show_debug_msg(__func__, "Pair percentage: %.2f/%.2f \n", pc, pair_pc);
	return pc >= pair_pc ? 1 : 0;
}

/**
 * At the right end of template 'left', if there are single reads, try
 * to get its right template id to merge.
 * If the mates of 'pair_pc' of the single reads are on right template 'r',
 * return r->id
 */
int right_tpl_to_merge(bwa_seq_t *seqs, tpl *left, float pair_pc) {
	bwa_seq_t *r = NULL, *m = NULL, *pre = NULL;
	int i = 0, dist = 0, max_tpl_id = 0;
	float max_n = 0.0, n_left = 0.0, pc = 0.0;
	int max_range = INS_SIZE + GRACE_TIMES * SD_INS_SIZE;
	GPtrArray *reads = g_ptr_array_sized_new(4);
	for (i = 0; i < left->reads->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(left->reads, i);
		m = get_mate(r, seqs);
		dist = left->len - r->contig_locus;
		if (dist > max_range) continue;
		if (!read_on_tpl(left, m)) {
			//p_query("USED", r);
			//p_query("MATE", m);
			g_ptr_array_add(reads, m);
		}
	}
	g_ptr_array_sort(reads, (GCompareFunc) cmp_reads_by_contig_id);
	for (i = 0; i < reads->len; i++) {
		m = (bwa_seq_t*) g_ptr_array_index(reads, i);
		//p_query("PAIR", m);
		if (m->contig_id <= 0) continue;
		if (!pre) {
			pre = m; n_left++; max_tpl_id = m->contig_id; continue;
		}
		if (m->contig_id == pre->contig_id) {
			r = get_mate(m, seqs);
			dist = left->len - r->contig_locus + m->contig_locus;
			//show_debug_msg(__func__, "Dist: %d \n", dist);
			if (r->rev_com == m->rev_com && good_insert_size(dist)) n_left++;
		} else {
			if (n_left > max_n) {
				max_tpl_id = pre->contig_id;
				max_n = n_left;
			}
			n_left = 1;
		}
		pre = m;
	}
	if (n_left > max_n) {
		//p_query("MAX", r);
		max_tpl_id = m->contig_id;
		max_n = n_left;
	}
	if (reads->len > 0) pc = max_n / ((float) reads->len);
	show_debug_msg(__func__, "Pair percentage: %.2f/%.2f \n", pc, pair_pc);
	g_ptr_array_free(reads, TRUE);
	return pc >= pair_pc ? max_tpl_id : -1;
}

/**
 * Try to connect the template 'b' to template 't' at locus 't_locus'.
 */
int connect_at_locus_right(hash_table *ht, tpl *t, tpl *b, int t_locus, int b_locus) {
	if (!b->alive || b->l_tail) return 0;
	int i = 0, con_pos = t_locus, cut_pos = b_locus;
	float n_not_paired = 0.0, n_spanning = 0.0;
	int t_part_len = 0, dist = 0, ori_len = 0;
	//show_debug_msg(__func__, "Left template [%d, %d] @ %d \n", t->id, t->len, t_locus);
	//show_debug_msg(__func__, "Right template [%d, %d] @ %d \n", b->id, b->len, b_locus);
	if (has_junction_at_locus(t, con_pos, 0)) return 0;
	for (i = 0; i < min(t->len - t_locus, b->len - b_locus); i++) {
		if (t->ctg->seq[i + t_locus] != b->ctg->seq[i + b_locus]) break;
		con_pos++; cut_pos++;
	}
	//show_debug_msg(__func__, "con_pos: %d; cut_pos: %d \n", con_pos, cut_pos);
	bwa_seq_t *r = NULL, *m = NULL;
	// Count pairs spanning the two templates
	for (i = 0; i < t->reads->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(t->reads, i);
		m = get_mate(r, ht->seqs);
		if (read_on_tpl(t, m)) continue;
		t_part_len = con_pos - r->contig_locus;
		if (t_part_len <= 0 || t_part_len > INS_SIZE + GRACE_TIMES * SD_INS_SIZE) continue;
		if (m->status != USED || m->contig_id != t->id) n_not_paired++;
		if (read_on_tpl(b, m)) {
			dist = (con_pos - r->contig_locus) + (m->contig_locus - cut_pos);
			if (r->rev_com == m->rev_com && good_insert_size(dist)) n_spanning++;
		}
	}
	if (n_spanning > 0 && n_spanning >= n_not_paired * PAIR_PERCENTAGE) {
		if (b_locus < ht->o->k * 4 && (t->len - con_pos) < 5 * ht->o->k) {
			p_tpl(t);
			show_debug_msg(__func__, "Left template [%d, %d] @ %d \n", t->id, t->len, t_locus);
			show_debug_msg(__func__, "Right template [%d, %d] @ %d \n", b->id, b->len, b_locus);
			show_debug_msg(__func__, "Spanning reads: %.2f/%.2f.\n", n_spanning, n_not_paired);
			show_debug_msg(__func__, "Merging: left [%d, %d] @ %d; right [%d, %d] %d. \n",
				t->id, t->len, con_pos, b->id, b->len, cut_pos);

			truncate_tpl(t, t->len - con_pos, 0, 0);
			truncate_tpl(b, cut_pos, 0, 1);
			ori_len = t->len;
			merge_tpl_to_left(t, b, 0, 0);
			refresh_tpl_reads(ht, t, ori_len - ht->o->read_len, ori_len + ht->o->read_len, N_MISMATCHES);
			correct_tpl_base(ht->seqs, t, ht->o->read_len, ori_len
					- ht->o->read_len, ori_len + ht->o->read_len);

			p_tpl(t);
			p_tpl_reads(t);
			return 1;
		}
	}
	return 0;
}

/**
 * Connect left end of the branch template to the main template;
 * Assumption: a->locus <= ht->o->k
 */
int connect_left_end(hash_table *ht, anchor *a, tpl *t) {
	int has_juc_read = 0, has_spanning_read = 0, connected = 0, similar = 0;
	int j = 0, len = 0, dist = 0;
	float n_spanning = 0.0, n = 0.0;
	tpl *b = a->b;
	bwa_seq_t *r = NULL, *m = NULL, *t_seq = NULL, *b_seq = NULL;
	len = min(t->len - a->from - a->size, b->len - a->from - a->size);
	t_seq = new_seq(t->ctg, len, a->from + a->size);
	b_seq = new_seq(b->ctg, len, a->locus + a->size);
	similar = similar_seqs(t_seq, b_seq, len / 10, MAX_GAPS,
					MATCH_SCORE, MISMATCH_SCORE, INDEL_SCORE);
	bwa_free_read_seq(1, t_seq);
	bwa_free_read_seq(1, b_seq);
	// If the branch is at least 90% similar with the main template,
	// It is likely to be constructed from 'squeezed' reads,
	// then merge it to the main template
	if (similar) {
		n_spanning = pairs_connect_left_locus(ht->seqs, t, b, a->from, a->locus, &n);
		if (n_spanning >= n * PAIR_PERCENTAGE) {
			mv_reads_to_main_tpl(t, b, a->from);
			b->alive = 0;
			return 1;
		}
		return 0;
	}

	// Make sure there are at least one pair at junction, and one pair spanning them
	for (j = 0; j < b->reads->len; j++) {
		r = (bwa_seq_t*) g_ptr_array_index(b->reads, j);
		m = get_mate(r, ht->seqs);
		dist = -1; has_juc_read = 0; has_spanning_read = 0;
		if (r->contig_locus < a->locus + a->size && r->contig_locus > a->locus + a->size) {
			if (read_on_tpl(t, m)) dist = (a->from - m->contig_locus) + (r->contig_locus - a->locus);
			if (read_on_tpl(b, m)) dist = (m->contig_locus - a->locus);
			if (good_insert_size(dist)) has_juc_read = 1;
			n++;
		}
		if (r->contig_locus >= a->locus + a->size && read_on_tpl(t, m)) {
			dist = (a->from - m->contig_locus) + (r->contig_locus - a->locus);
			if (good_insert_size(dist)) has_spanning_read = 1;
		}
	}
	if (has_spanning_read && has_juc_read) {
		truncate_tpl(b, a->locus, 0, 1);
		add_a_junction(t, b, NULL, a->from + a->size, 0, n);
		return 1;
	}
	return 0;
}

/**
 * Connect one end of a template to some location of the main template.
 * Create a junction.
 */
int connect_one_end(hash_table *ht, GPtrArray *anchors, tpl *t) {
	if (!anchors || anchors->len <= 0 || !t || !t->alive) return 0;
	tpl *b = NULL;
	int connected = 0;
	int i = 0, j = 0;
	anchor *a = NULL;
	bwa_seq_t *b_seq = NULL, *t_seq = NULL, *r = NULL, *m = NULL;
	show_debug_msg(__func__, "Anchors with template [%d, %d] \n", t->id, t->len);
	for (i = 0; i < anchors->len; i++) {
		a = (anchor*) g_ptr_array_index(anchors, i);
		p_anchor("BRANCHING", a);
		b = a->b;
		if (t == b || !b->alive || a->size > ht->o->read_len) continue;
		if (a->locus <= ht->o->k) { // Try to connect left end of the branch to the main
			connected = connect_left_end(ht, a, t);
		}
	}
	return connected;
}

/**
 * Description: http://caishaojiang.com/2014/05/01/peta-connect-both-ends-of-branch-template-to-main-template/
 * t: ----------------------------------------------------
 * b:               a1 ---------         ----------- a2
 * Assumption: the anchors are already sorted by 'hash' attribute;
 * 	That means, first sorted by t->id, then sorted by 'locus'
 */
int connect_both_ends(hash_table *ht, GPtrArray *anchors, tpl *t) {
	if (!anchors || anchors->len <= 0 || !t || !t->alive) return 0;
	tpl *b = NULL;
	int i = 0, j = 0, x = 0, k = ht->o->k, tmp = 0, adjusted = 0;
	int connected = 0, dist = 0, n_left = 0, n_right = 0, con_left = 0, con_right = 0;
	int t_mid = 0, b_mid = 0, similar = 0;
	float n_pairs = 0.0, n_single_reads = 0.0;
	anchor *a = NULL, *a1 = NULL, *a2 = NULL;
	bwa_seq_t *r = NULL, *m = NULL, *t_seq = NULL, *b_seq = NULL;
	GPtrArray *junc_reads = NULL;
	show_debug_msg(__func__, "Anchors shared with [%d, %d] \n", t->id, t->len);
//	for (i = 0; i < anchors->len; i++) {
//		a1 = (anchor*) g_ptr_array_index(anchors, i);
//		p_anchor("ANCHOR", a1);
//	}
	for (i = 0; i < anchors->len; i++) {
		a1 = (anchor*) g_ptr_array_index(anchors, i);
		if (a1->size <= 0 || a1->b == t) continue;
		b = a1->b;
		if (!b->alive || has_any_junction(b)) continue;
		for (j = i; j < anchors->len; j++) {
			a2 = (anchor*) g_ptr_array_index(anchors, j);
			if (a1->b != a2->b) break;
			if (a1->size <= 0 || a2->size <= 0) continue;
			if (a1 == a2 || !a2->b->alive || a1->from >= a2->from) continue;
			if (a1->size + a2->size < ht->o->read_len || has_any_junction(b)) continue;
			if (abs(a1->from + a1->size - a2->from) <= 2 * k && abs(a2->locus - a1->locus - a1->size) <= 2 * k) continue;
			if (a1->locus >= ht->o->k * 2 || b->len - a2->locus - a2->size >= ht->o->k * 2) continue;
			p_tpl(t); p_tpl(b);
			p_anchor("A1", a1); p_anchor("A2", a2);
			show_debug_msg(__func__, "Connecting branch template [%d, %d] to main template [%d, %d] ... \n",
					b->id, b->len, t->id, t->len);
			con_left = a1->from + a1->size;
			con_right = a2->from;
			adjusted = 0;
			if (con_right - con_left < 0 || a2->locus < a1->locus + a1->size) {
				show_debug_msg(__func__, "Left position %d larger than right position %d. \n", con_left, con_right);
				// Either the main branch anchors overlap or branch anchors overlap
				dist = max(a1->from + a1->size - a2->from, a1->locus + a1->size - a2->locus);
				tmp = a2->locus + dist;
				for (x = 1; x <= dist; x++) {
					show_debug_msg(__func__, "LEFT: %d; RIGHT: %d \n", tmp-x, a1->locus + a1->size - x);
					show_debug_msg(__func__, "LEFT: %c; RIGHT: %c \n", "ACGTN"[b->ctg->seq[tmp - x]],
							"ACGTN"[b->ctg->seq[a1->locus + a1->size - x]]);
					if (b->ctg->seq[tmp - x] == b->ctg->seq[a1->locus + a1->size - x]) {
						adjusted++;
					} else break;
				}
				show_debug_msg(__func__, "Adjusted: left %d; right %d. \n", con_left, con_right + adjusted);
			}
			con_right += adjusted;
			// Not allow circle in the graph
			if (con_left > con_right) continue;
			if ((b->len - a2->locus - adjusted) + (a1->locus + a1->size) > b->len) {
				show_debug_msg(__func__, "[WARNING] Branch template [%d, %d] is too short to be truncated: Right %d; Left %d\n",
						b->id, b->len, (b->len - a2->locus - adjusted), (a1->locus + a1->size));
				continue;
			}

			t_mid = a2->from + adjusted - a1->from - a1->size;
			b_mid = a2->locus + adjusted - a1->locus - a1->size;
			if (abs(t_mid - b_mid) <= k) {
				t_seq = new_seq(t->ctg, t_mid, a1->from + a1->size);
				b_seq = new_seq(b->ctg, b_mid, a1->locus + a1->size);
				p_query("MAIN_MIDDLE", t_seq);
				p_query("BRAN MIDDLE", b_seq);
				similar = similar_seqs(b_seq, b_seq, (1 - SM_SIMILARY) * min(t_mid, b_mid),
						MAX_GAPS, MATCH_SCORE, MISMATCH_SCORE, INDEL_SCORE);
				bwa_free_read_seq(1, t_seq);
				bwa_free_read_seq(1, b_seq);
				if (similar) {
					show_debug_msg(__func__, "The middle sequences are similar, not connected. \n");
					continue;
				}
			}

			junc_reads = g_ptr_array_sized_new(4);
			n_left = 0; n_right = 0;
			for (x = 0; x < b->reads->len; x++) {
				r = (bwa_seq_t*) g_ptr_array_index(b->reads, x);
				if (r->contig_locus < a1->locus + a1->size && r->contig_locus + r->len > a1->locus + a1->size) {
					n_left++; g_ptr_array_add(junc_reads, r);
				}
				if (r->contig_locus < a2->locus + adjusted && r->contig_locus + r->len > a2->locus + adjusted) {
					n_right++; g_ptr_array_add(junc_reads, r);
				}
			}
			show_debug_msg(__func__, "Junction reads: %d \n", junc_reads->len);
			if (n_left >= MIN_PAIRS && n_right >= MIN_PAIRS) {
				n_pairs = 0.0; n_single_reads = 0.0;
				for (x = 0; x < junc_reads->len; x++) {
					r = (bwa_seq_t*) g_ptr_array_index(junc_reads, x);
					m = get_mate(r, ht->seqs);
					if (r->rev_com != m->rev_com) continue;
					if (r->contig_locus < a1->locus || r->contig_locus + r->len > a2->locus + a2->size) continue;
					//p_query(__func__, r);
					//p_query(__func__, m);
					if (read_on_tpl(b, m)) continue;
					n_single_reads++; dist = -1;
					if (read_on_tpl(t, m)) {
						if (m->contig_locus < a1->from + a1->size)
							dist = (a1->from - m->contig_locus) + (r->contig_locus - a1->size);
						if (m->contig_locus >= a2->from + adjusted)
							dist = (a2->locus + adjusted - r->contig_locus) + (m->contig_locus - (a2->from + adjusted));
						if (good_insert_size(dist)) n_pairs++;
						else g_ptr_array_remove_index_fast(junc_reads, x--);
					}
					//show_debug_msg(__func__, "Pairs: %.2f/%.2f; Distance: %d \n", n_pairs, n_single_reads, dist);
				}
			}
			g_ptr_array_free(junc_reads, TRUE);
			if (n_pairs >= MIN_PAIRS && n_pairs >= n_single_reads * PAIR_PERCENTAGE) {
				shrink_anchor_right(a2, adjusted);
				add_a_junction(t, b, NULL, a2->from, 1, n_right);
				add_a_junction(t, b, NULL, a1->from + a1->size, 0, n_left);

				show_debug_msg(__func__, "N_LEFT: %d \n", n_left);
				show_debug_msg(__func__, "N_RIGHT: %d \n", n_right);
				p_tpl(t); p_tpl(b);
				show_debug_msg(__func__, "Truncating template [%d, %d] at %d; right end \n",
						b->id, b->len, b->len - a2->locus);

				truncate_tpl(b, b->len - a2->locus, 0, 0);
				set_tail(b, t, a1->from + a1->size, ht->o->read_len - 1, 0);

				show_debug_msg(__func__, "Truncating template [%d, %d] at %d; left end \n",
						b->id, b->len, a1->locus + a1->size);

				truncate_tpl(b, a1->locus + a1->size, 0, 1);
				set_tail(b, t, a2->from, ht->o->read_len - 1, 1);
				disable_anchor(a1); disable_anchor(a2);
				connected = 1;
				//p_tpl(b); p_tpl_reads(b);
			}
		}
	}
	return connected;
}
