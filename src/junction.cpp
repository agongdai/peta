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
#include "tpl.hpp"
#include "ass.hpp"
#include "rnaseq.h"
#include "pechar.h"
#include "utils.h"
#include "peseq.h"
#include "junction.hpp"
#include "graph.hpp"
#include "k_hash.h"

using namespace std;

junction *new_junction(tpl *main_tpl, tpl *branch_tpl, bwa_seq_t * connector,
		int locus, int ori, int weight) {
	junction *j = (junction*) malloc(sizeof(junction));
	j->main_tpl = main_tpl;
	j->branch_tpl = branch_tpl;
	j->locus = main_tpl->len;
	j->ori = ori;
	j->locus = locus;
	j->weight = weight;
	j->connector = connector;
	j->reads = g_ptr_array_sized_new(0);
	j->status = 0;
	return j;
}

/**
 * Add a new junction;
 * The new junction exists on both main and branch templates
 */
junction *add_a_junction(tpl *main_tpl, tpl *branch_tpl, bwa_seq_t *connector,
		int locus, int ori, int weight) {
	// Indicating the templates are in-connect, cannot be reverse-complement
	junction *new_j = new_junction(main_tpl, branch_tpl, connector, locus, ori,
			weight);
	// Add the junction to the templates. For the 'existing connect' later.
	if (!main_tpl->m_juncs)
		main_tpl->m_juncs = g_ptr_array_sized_new(2);
	g_ptr_array_add(main_tpl->m_juncs, new_j);
	if (!branch_tpl->b_juncs)
		branch_tpl->b_juncs = g_ptr_array_sized_new(2);
	g_ptr_array_add(branch_tpl->b_juncs, new_j);
	return new_j;
}

void destroy_junction(junction *j) {
	tpl *main_tpl = NULL, *branch_tpl = NULL;
	int i = 0;
	junction *jun = NULL;
	if (j) {
		main_tpl = j->main_tpl;
		branch_tpl = j->branch_tpl;

		/**
		 show_debug_msg(__func__, "Destroying junction ... \n");
		 p_junction(j);
		 p_junctions(main_tpl->m_juncs);
		 p_junctions(main_tpl->b_juncs);
		 p_junctions(branch_tpl->m_juncs);
		 p_junctions(branch_tpl->b_juncs);
		 **/

		if (j->reads)
			g_ptr_array_free(j->reads, TRUE);
		if (j->locus == 1 && j->branch_tpl->r_tail) {
			bwa_free_read_seq(1, j->branch_tpl->r_tail);
			j->branch_tpl->r_tail = NULL;
		}
		if (j->locus == 0 && j->branch_tpl->l_tail) {
			bwa_free_read_seq(1, j->branch_tpl->l_tail);
			j->branch_tpl->l_tail = NULL;
		}
		if (main_tpl->m_juncs) {
			for (i = main_tpl->m_juncs->len - 1; i >= 0; i--) {
				jun = (junction*) g_ptr_array_index(main_tpl->m_juncs, i);
				if (jun == j) {
					g_ptr_array_remove_index_fast(main_tpl->m_juncs, i--);
					break;
				}
			}
		}
		if (branch_tpl->b_juncs) {
			for (i = branch_tpl->b_juncs->len - 1; i >= 0; i--) {
				jun = (junction*) g_ptr_array_index(branch_tpl->b_juncs, i);
				if (jun == j) {
					g_ptr_array_remove_index_fast(branch_tpl->b_juncs, i--);
					break;
				}
			}
		}
		free(j);
	}
}

gint cmp_junc_by_locus(gpointer a, gpointer b) {
	junction *c_a = *((junction**) a);
	junction *c_b = *((junction**) b);
	return ((c_a->locus) - c_b->locus);
}

gint cmp_junc_by_branch_id(gpointer a, gpointer b) {
	junction *c_a = *((junction**) a);
	junction *c_b = *((junction**) b);
	return ((c_a->branch_tpl->id) - c_b->branch_tpl->id);
}

int same_main_juncs(junction *jun_1, junction *jun_2) {
	if (!jun_1 || !jun_2)
		return 0;
	if (jun_1->main_tpl == jun_2->main_tpl && jun_1->locus == jun_2->locus
			&& jun_1->ori == jun_2->ori)
		return 1;
}

/**
 * Get junctions whose branch tpl is the given tpl branch
 */
GPtrArray *find_branch_junctions(GPtrArray *all, tpl *branch) {
	uint32_t i = 0;
	junction *j = NULL;
	GPtrArray *hits = g_ptr_array_sized_new(2);
	for (i = 0; i < all->len; i++) {
		j = (junction*) g_ptr_array_index(all, i);
		if (j->branch_tpl == branch)
			g_ptr_array_add(hits, j);
	}
	return hits;
}

/**
 * Concatenate all reads on nearby templates;
 * For paired reads finding
 * The attribute is_root is used as an indicator temporarily.
 */
GPtrArray *get_nearby_tpls(tpl *t, GPtrArray *tpls) {
	int j = 0, i = 0;
	junction *jun = NULL;
	if (!tpls)
		tpls = g_ptr_array_sized_new(4);
	//p_tpl(t);
	if (t->is_root)
		return tpls;
	// Set to 1, means that the reads on it have been visited.
	t->is_root = 1;
	g_ptr_array_add(tpls, t);
	//show_debug_msg(__func__, "Nearby template: [%d, %d] \n", t->id, t->len);
	if (t->b_juncs && t->b_juncs->len > 0) {
		for (j = 0; j < t->b_juncs->len; j++) {
			jun = (junction*) g_ptr_array_index(t->b_juncs, j);
			//p_junction(jun);
			//p_tpl(jun->main_tpl);
			if (jun->main_tpl->is_root == 0)
				get_nearby_tpls(jun->main_tpl, tpls);
		}
	}
	return tpls;
}

void reset_is_root(GPtrArray *tpls) {
	int i = 0;
	tpl *t = NULL;
	for (i = 0; i < tpls->len; i++) {
		t = (tpl*) g_ptr_array_index(tpls, i);
		t->is_root = 0;
	}
}

GPtrArray *nearby_tpls(tpl *t) {
	GPtrArray *near = g_ptr_array_sized_new(4);
	near = get_nearby_tpls(t, near);
	reset_is_root(near);
	return near;
}

GPtrArray *find_junc_reads(hash_table *ht, bwa_seq_t *left, bwa_seq_t *right,
		const int max_len, int *weight) {
	int left_len = 0, right_len = 0, n_reads = 0;
	GPtrArray *reads = g_ptr_array_sized_new(0), *hits = NULL;
	bwa_seq_t *window = NULL, *junc_seq = blank_seq(max_len), *r = NULL;
	int i = 0, j = 0;

	left_len = (left->len > max_len / 2) ? (max_len / 2) : left->len;
	memcpy(junc_seq->seq, left->seq + (left->len - left_len), sizeof(ubyte_t)
			* left_len);
	right_len = (right->len) > (max_len / 2) ? (max_len / 2) : (right->len);
	memcpy(junc_seq->seq + left_len, right->seq, sizeof(ubyte_t) * right_len);
	junc_seq->len = left_len + right_len;
	set_rev_com(junc_seq);
	//p_query("Left  seq", left);
	//p_query("Right seq", right);
	//p_query("Junction seq", junc_seq);

	for (i = 0; i <= junc_seq->len - ht->o->read_len; i++) {
		window = new_seq(junc_seq, ht->o->read_len, i);
		hits = g_ptr_array_sized_new(4);
		hits = find_both_fr_full_reads(ht, window, hits, N_MISMATCHES);
		for (j = 0; j < hits->len; j++) {
			r = (bwa_seq_t*) g_ptr_array_index(hits, j);
			g_ptr_array_add(reads, r);
		}
		g_ptr_array_free(hits, TRUE);
		bwa_free_read_seq(1, window);
	}

	n_reads = reads->len;
	//show_debug_msg(__func__, "# of junction reads: %d \n", n_reads);
	*weight = n_reads;
	bwa_free_read_seq(1, junc_seq);
	return reads;
}

/**
 * Align the reads to a junction sequence and return the number
 */
int count_jun_reads(hash_table *ht, junction *jun) {
	tpl *main_tpl = NULL, *branch = NULL;
	bwa_seq_t *left = NULL, *right = NULL, *branch_seq = NULL, *main_seq = NULL;
	int n_reads = 0, b_l_len = 0, b_r_len = 0, b_t_len = 0;
	int m_l_len = 0, m_r_len = 0, m_t_len = 0;
	GPtrArray *j_reads = NULL;
	if (!jun || jun->status != 0)
		return 0;
	//show_debug_msg(__func__, "Setting junction reads...\n");
	main_tpl = jun->main_tpl;
	branch = jun->branch_tpl;
	branch_seq = get_tpl_ctg_wt(branch, &b_l_len, &b_r_len, &b_t_len);
	main_seq = get_tpl_ctg_wt(main_tpl, &m_l_len, &m_r_len, &m_t_len);
	//p_junction(jun);
	//p_tpl(main_tpl);
	//p_tpl(branch);
	left = jun->ori ? new_seq(branch_seq, branch_seq->len - b_r_len, 0) : new_seq(
			main_seq, jun->locus + m_l_len, 0);
	right = jun->ori ? new_seq(main_seq, main_seq->len - jun->locus - m_l_len,
			jun->locus + m_l_len) : new_seq(branch_seq, branch->len + b_r_len, b_l_len);
	//p_ctg_seq("Main", main_tpl->ctg);
	//p_ctg_seq("Bran", branch->ctg);
	j_reads = find_junc_reads(ht, left, right, (ht->o->read_len
			- JUNCTION_BOUNDARY_BASE) * 2, &n_reads);
	//p_readarray(j_reads, 1);
	g_ptr_array_free(j_reads, TRUE);
	bwa_free_read_seq(1, left);
	bwa_free_read_seq(1, right);
	bwa_free_read_seq(1, branch_seq);
	bwa_free_read_seq(1, main_seq);
	return n_reads;
}

bwa_seq_t *get_junc_seq(tpl *left, int l_pos, int *left_len, tpl *right,
		int r_pos, int *right_len, int max_len) {
	bwa_seq_t *junc_seq = blank_seq(max_len);
	bwa_seq_t *left_seq = NULL, *right_seq = NULL;
	left_seq = cut_tpl_tail(left, l_pos, max_len / 2, 0);
	right_seq = cut_tpl_tail(right, r_pos, max_len / 2, 1);
	*left_len = left_seq->len;
	*right_len = right_seq->len;
	memcpy(junc_seq->seq, left_seq->seq, left_seq->len);
	memcpy(junc_seq->seq + left_seq->len, right_seq->seq, right_seq->len);
	junc_seq->len = left_seq->len + right_seq->len;
	bwa_free_read_seq(1, left_seq);
	bwa_free_read_seq(1, right_seq);
	return junc_seq;
}

GPtrArray *find_junc_reads_w_tails(hash_table *ht, tpl *main_tpl,
		tpl *branch_tpl, const int pos, const int max_len, const int ori,
		int *weight) {
	bwa_seq_t *left_seq = NULL, *right_seq = NULL;
	tpl *left_eg = branch_tpl, *right_eg = main_tpl;
	int l_pos = branch_tpl->len, r_pos = pos;
	GPtrArray *reads = g_ptr_array_sized_new(4);
	if (ori) {
		left_eg = main_tpl;
		l_pos = pos;
		right_eg = branch_tpl;
		r_pos = 0;
	}
	//show_debug_msg(__func__, "left pos: %d; right pos: %d\n", l_pos, r_pos);
	left_seq = cut_tpl_tail(left_eg, l_pos, max_len / 2, 0);
	//p_query("Left  seq", left_eg->ctg);
	//p_query("Left tail", left_seq);
	right_seq = cut_tpl_tail(right_eg, r_pos, max_len / 2, 1);
	//p_query("Right  seq", right_eg->ctg);
	//p_query("Right tail", right_seq);
	reads = find_junc_reads(ht, left_seq, right_seq, max_len, weight);
	bwa_free_read_seq(1, left_seq);
	bwa_free_read_seq(1, right_seq);
	return reads;
}

/**
 * Validate the junction by checking mate pairs.
 * Depending on ori and con_pos, only partial reads on the main template are counted
 */
int vld_junc_by_mates(tpl *main_tpl, tpl *branch_tpl, GPtrArray *junc_reads,
		hash_table *ht, const int con_pos, const int ins_size, const int ori) {
	int start = con_pos, end = main_tpl->len;
	int is_valid = 0, seg_len = 0;
	GPtrArray *tmp = NULL;
	//return 1;
	seg_len += branch_tpl->len;
	if (ori) {
		seg_len += con_pos;
		if (main_tpl->l_tail)
			seg_len += main_tpl->l_tail->len;
		if (branch_tpl->r_tail)
			seg_len += branch_tpl->r_tail->len;
	} else {
		seg_len += main_tpl->len - con_pos;
		if (main_tpl->r_tail)
			seg_len += main_tpl->r_tail->len;
		if (branch_tpl->l_tail)
			seg_len += branch_tpl->l_tail->len;
	}
	// If the length is too short, do not do the validation
	if (seg_len < ins_size + 100)
		return 1;
	if (ori) {
		start = 0;
		end = con_pos - ht->o->read_len;
	}

	/**
	 show_debug_msg(__func__, "Main: [%d, %d]; Branch: [%d, %d]; Start~End: [%d, %d]\n",
	 main_tpl->id, main_tpl->len, branch_tpl->id, branch_tpl->len, start, end);
	 show_debug_msg(__func__, "Range [%d, %d]; ORI: %d\n", start, end, ori);
	 show_debug_msg(__func__, "Main template %d reads: %d\n", main_tpl->id,
	 main_tpl->reads->len);
	 p_readarray(main_tpl->reads, 1);
	 show_debug_msg(__func__, "Branch template %d reads: %d\n", branch_tpl->id,
	 branch_tpl->reads->len);
	 p_ctg_seq("BRANCH", branch_tpl->ctg);
	 p_readarray(branch_tpl->reads, 1);
	 **/

	g_ptr_array_sort(branch_tpl->reads, (GCompareFunc) cmp_reads_by_name);
	is_valid = vld_tpl_mates(branch_tpl, main_tpl, start, end, MIN_PAIRS);
	if (!is_valid) {
		/**
		 show_debug_msg(__func__, "Tag 1 \n");
		 p_readarray(main_tpl->reads, 1);
		 show_debug_msg(__func__, "Tag 2\n");
		 p_readarray(junc_reads, 1);
		 **/
		// Maybe exon shorter than read length, the mates located at the junction
		g_ptr_array_sort(junc_reads, (GCompareFunc) cmp_reads_by_name);
		is_valid = find_pairs(junc_reads, main_tpl->reads, 0, main_tpl->id,
				start, end, MIN_PAIRS);
		if (!is_valid) {
			is_valid = find_pairs(junc_reads, branch_tpl->reads, 0,
					branch_tpl->id, 0, branch_tpl->len, MIN_PAIRS);
		}
	}
	return is_valid;
}

/**
 * Update the junction locus for those tpls connected to itself.
 * Because the locus is not correct when the junction is recorded.
 * Example:
 * Firstly extending to the right and connect to itself:
 *           <<<<<
 *          /     \
 *         /       \
 *         |        ^
 * ------------------
 * Locus:  ^ (value 8)
 * But, later it extends to the left:
 *               <<<<<
 *              /     \
 *             /       \
 *             |        ^
 * ====------------------
 * Locus should be 8 + 4
 */
void upd_tpl_jun_locus(tpl *t, GPtrArray *branching_events, const int extra_len) {
	int i = 0;
	uint64_t query_int = 0, k = 0;
	junction *jun = NULL;
	if (!t->b_juncs || t->b_juncs->len <= 0 || !branching_events
			|| branching_events->len <= 0)
		return;
	jun = (junction*) g_ptr_array_index(t->b_juncs, 0);
	if (jun->branch_tpl == jun->main_tpl)
		jun->locus += extra_len;
}

gint cmp_junc_by_id(gpointer a, gpointer b) {
	junction *c_a = *((junction**) a);
	junction *c_b = *((junction**) b);
	return ((c_a->main_tpl->id) - c_b->main_tpl->id);
}

void store_features(char *name, GPtrArray *branching_events,
		GPtrArray *all_tpls) {
	junction *jun = NULL;
	uint64_t i = 0;
	tpl *t = NULL;
	char entry[BUFSIZE];
	FILE *f = xopen(name, "w");
	sprintf(entry, "Main\tBranch\tLocus\tWeight\tDirection\n");
	fputs(entry, f);
	//	for (i = 0; i < all_tpls->len; i++) {
	//		t = (tpl*) g_ptr_array_index(all_tpls, i);
	//		sprintf(entry, "[%d, %d]\t[%d, %d]\t0\t%d\t-1\n", t->id, t->len, t->id,
	//				t->len, t->kmer_freq);
	//		fputs(entry, f);
	//	}
	for (i = 0; i < branching_events->len; i++) {
		jun = (junction*) g_ptr_array_index(branching_events, i);
		p_junction(jun);
		sprintf(entry, "[%d, %d]\t[%d, %d]\t%d\t%d\t%d\n", jun->main_tpl->id,
				jun->main_tpl->len, jun->branch_tpl->id, jun->branch_tpl->len,
				jun->locus, jun->weight, jun->ori);
		fputs(entry, f);
	}
	fclose(f);
}

void p_junction(junction *jun) {
	if (!jun)
		show_debug_msg(__func__, "Junction is NULL.\n");
	else
		show_debug_msg(__func__,
				"[%d%c, %d]\t[%d%c, %d]\t%d\t%d\t%d\tStatus:%d\n",
				jun->main_tpl->id, jun->main_tpl->alive ? '@' : '!',
				jun->main_tpl->len, jun->branch_tpl->id,
				jun->branch_tpl->alive ? '@' : '!', jun->branch_tpl->len,
				jun->locus, jun->weight, jun->ori, jun->status);
}

void p_junctions(GPtrArray *juns) {
	junction *jun = NULL;
	int i = 0;
	if (juns) {
		show_debug_msg(__func__, "------------- %d junctions ------------- \n",
				juns->len);
		for (i = 0; i < juns->len; i++) {
			jun = (junction*) g_ptr_array_index(juns, i);
			printf("%d\t", i);
			p_junction(jun);
		}
		show_debug_msg(__func__, "------------- %d junctions ------------- \n",
				juns->len);
	} else {
		show_debug_msg(__func__, "0 junctions \n");
	}
}

/**
 * Remove those junctions with status not equal to 0 from the list
 */
void remove_dead_junctions(GPtrArray *junctions) {
	junction *j = NULL;
	uint32_t i = 0;
	for (i = 0; i < junctions->len; i++) {
		j = (junction*) g_ptr_array_index(junctions, i);
		p_junction(j);
		//p_tpl(j->main_tpl);
		//p_tpl(j->branch_tpl);
		if (j->status != 0) {
			destroy_junction(j);
			g_ptr_array_remove_index_fast(junctions, i--);
		}
	}
}

/**
 * Remove those junctions with not-alive templates
 */
void rm_junc_w_dead_tpls(GPtrArray *junctions, tpl *t) {
	int i = 0;
	junction *jun = NULL;
	for (i = 0; i < junctions->len; i++) {
		jun = (junction*) g_ptr_array_index(junctions, i);
		if (jun->main_tpl == t || jun->branch_tpl == t || jun->status != 0) {
			jun->status = 1;
			//destroy_junction(jun);
			g_ptr_array_remove_index_fast(junctions, i--);
		}
	}
}

void destory_tpl_junctions(tpl *t) {
	int i = 0, j = 0;
	junction *jun = NULL, *jun2 = NULL;
	if (!t)
		return;
	if (t->m_juncs) {
		while (t->m_juncs->len > 0) {
			jun = (junction*) g_ptr_array_index(t->m_juncs, 0);
			destroy_junction(jun);
		}
	}
	if (t->b_juncs) {
		while (t->b_juncs->len > 0) {
			jun = (junction*) g_ptr_array_index(t->b_juncs, 0);
			destroy_junction(jun);
		}
	}
}

/**
 * Remove duplicate junctions and junctions with dead templates
 */
void get_junction_arr(GPtrArray *read_tpls, GPtrArray *junctions) {
	junction *junc = NULL, *pre = NULL;
	uint32_t i = 0, j = 0;
	int not_alive = 0;
	tpl *t = NULL;
	for (i = 0; i < read_tpls->len; i++) {
		t = (tpl*) g_ptr_array_index(read_tpls, i);
		if (t->m_juncs) {
			for (j = 0; j < t->m_juncs->len; j++) {
				junc = (junction*) g_ptr_array_index(t->m_juncs, j);
				if (junc->status == 0) {
					g_ptr_array_add(junctions, junc);
					junc->status = 1;
				}
			}
		}
		if (t->b_juncs) {
			for (j = 0; j < t->b_juncs->len; j++) {
				junc = (junction*) g_ptr_array_index(t->b_juncs, j);
				if (junc->status == 0) {
					g_ptr_array_add(junctions, junc);
					junc->status = 1;
				}
			}
		}
	}
	for (i = 0; i < junctions->len; i++) {
		junc = (junction*) g_ptr_array_index(junctions, i);
		junc->status = 0;
	}
}

/**
 * Mark a template not alive;
 * Mark all junctions regarding this template as not alive
 */
void mark_tpl_dead(tpl *t) {
	destory_tpl_junctions(t);
	t->alive = 0;
}

/**
 * Main:  ================================
 * Shift:                 ^
 * Ori: 0 (to the right)
 * Branch:                 ---------
 * Check:                  =========
 *                         |||||||||
 *                         ---------
 */
int branch_on_main(tpl *main_tpl, tpl *branch, const int pos,
		const int mismatches, const int exist_ori) {
	bwa_seq_t *sub = NULL, *full = NULL;
	int similar = 0, len_not_valid = 0;
	int l_len = 0, r_len = 0, t_len = 0;
	full = get_tpl_ctg_wt(main_tpl, &l_len, &r_len, &t_len);
	if (exist_ori) {
		if (pos + l_len < branch->len)
			len_not_valid = 1;
		else
			sub = new_seq(full, branch->len, pos - branch->len + l_len);
	} else {
		if (full->len - pos - l_len < branch->len)
			len_not_valid = 1;
		else
			sub = new_seq(full, branch->len, pos + l_len);
	}

	//p_ctg_seq("FULL", full);
	//p_ctg_seq(__func__, sub);
	//p_ctg_seq(__func__, branch->ctg);
	if (!len_not_valid) {
		//similar = (seq_ol(sub, branch, branch->len, mismatches) == -1) ? 0 : 1;
		similar = similar_seqs(sub, branch->ctg, mismatches, MAX_GAPS,
				MATCH_SCORE, MISMATCH_SCORE, INDEL_SCORE);
		//show_debug_msg(__func__, "Mismatches: %d; similar: %d\n", mismatches,
		//		similar);
	} else {
		similar = 0;
	}
	bwa_free_read_seq(1, sub);
	bwa_free_read_seq(1, full);
	return similar;
}

/**
 * Get junctions with the template t as main/branch
 * all_juncs are supposed to be ordered in main/branch template id increasingly
 */
GPtrArray *tpl_junctions(tpl *t, GPtrArray *all_juncs, int start_index,
		int to_get_main) {
	uint32_t i = 0;
	junction *j = NULL;
	GPtrArray *juncs = g_ptr_array_sized_new(32);
	for (i = start_index; i < all_juncs->len; i++) {
		j = (junction*) g_ptr_array_index(all_juncs, i);
		if (to_get_main) {
			if (j->main_tpl == t)
				g_ptr_array_add(juncs, j);
			else
				break;
		} else {
			if (j->branch_tpl == t)
				g_ptr_array_add(juncs, j);
			else
				break;
		}
	}
	return juncs;
}

/**
 * Mark junctions as 'dead' if:
 * 1. Branch can be mereged into the main template
 * 2. Two junctions are too close (<= 2bp)
 */
void filter_branches(GPtrArray *junctions, const int read_len) {
	uint32_t i = 0;
	int start = 0, end = 0, similar = 0;
	junction *cur = NULL, *pre = NULL;
	bwa_seq_t *main_seq = NULL, *branch_seq = NULL;
	tpl *main_tpl = NULL, *branch_tpl = NULL;
	for (i = 0; i < junctions->len; i++) {
		cur = (junction*) g_ptr_array_index(junctions, i);
		if (cur->status != 0)
			continue;
		if (!pre) {
			pre = (junction*) g_ptr_array_index(junctions, 0);
			continue;
		}
		if (cur->branch_tpl != pre->branch_tpl) {
			pre = cur;
			continue;
		}
		if (abs(cur->locus - pre->locus) <= 2 && cur->branch_tpl->len <= 2) {
			pre->status = 1;
			cur->status = 1;
		} else {
			if (abs(cur->locus - pre->locus) < read_len && pre->ori != cur->ori) {
				branch_seq = cur->branch_tpl->ctg;
				start = min(cur->locus, pre->locus);
				end = max(cur->locus, pre->locus);
				main_seq = new_seq(cur->main_tpl->ctg, end - start, start);
				similar = seq_ol(main_seq, branch_seq, main_seq->len,
						N_MISMATCHES * 2);
				if (similar >= 0) {
					pre->status = 1;
					cur->status = 1;
				}
				bwa_free_read_seq(1, main_seq);
				/**
				 p_junction(pre);
				 p_junction(cur);
				 p_ctg_seq(__func__, main_seq);
				 p_ctg_seq(__func__, branch_seq);
				 show_debug_msg(__func__, "SIMILAR: %d\n", similar);
				 **/
			}
		}
		pre = cur;
	}
}

/**
 * Mark dead the junctions which are short and cannot be connected to other templates
 */
void prune_short_branches(GPtrArray *junctions, GPtrArray *tpls,
		const int max_len) {
	junction *jun = NULL;
	uint32_t i = 0, m = 0, start = 0;
	tpl *t = NULL;
	int id = 0;
	junc_count main_count, branch_count;
	g_ptr_array_sort(tpls, (GCompareFunc) cmp_tpl_by_id);
	g_ptr_array_sort(junctions, (GCompareFunc) cmp_junc_by_id);

	for (i = 0; i < tpls->len; i++) {
		t = (tpl*) g_ptr_array_index(tpls, i);
		for (m = start; m < junctions->len; m++) {
			jun = (junction*) g_ptr_array_index(junctions, m);
			//p_junction(jun);
			if (jun->main_tpl != t)
				break;
		}
		main_count[t->id] = m - start;
		//show_debug_msg(__func__, "Main of %d: %d \n", t->id, m-start);
		start = m;
	}

	start = 0;
	g_ptr_array_sort(junctions, (GCompareFunc) cmp_junc_by_branch_id);
	for (i = 0; i < tpls->len; i++) {
		t = (tpl*) g_ptr_array_index(tpls, i);
		for (m = start; m < junctions->len; m++) {
			jun = (junction*) g_ptr_array_index(junctions, m);
			if (jun->branch_tpl != t)
				break;
		}
		branch_count[t->id] = m - start;
		start = m;
	}

	for (i = 0; i < junctions->len; i++) {
		jun = (junction*) g_ptr_array_index(junctions, i);
		if (jun->branch_tpl->len > max_len)
			continue;
		id = jun->branch_tpl->id;
		if (main_count[id] == 0 && branch_count[id] == 1) {
			jun->status = 1;
		}
	}
	main_count.clear();
	branch_count.clear();
}

/**
 * Check whether there is any junction between the two templates
 * Return: 1 means 'yes'
 */
int tpls_have_junction(tpl *left, tpl *right) {
	int i = 0, have = 0;
	junction *jun = NULL;
	if (right->b_juncs) {
		for (i = 0; i < right->b_juncs->len; i++) {
			jun = (junction*) g_ptr_array_index(right->b_juncs, i);
			if ((jun->branch_tpl == right && jun->main_tpl == left)
					|| (jun->branch_tpl == left && jun->main_tpl == right)) {
				have = 1;
			}
		}
	}
	if (!have && right->m_juncs) {
		for (i = 0; i < right->m_juncs->len; i++) {
			jun = (junction*) g_ptr_array_index(right->m_juncs, i);
			if ((jun->branch_tpl == right && jun->main_tpl == left)
					|| (jun->branch_tpl == left && jun->main_tpl == right)) {
				have = 1;
			}
		}
	}
	if (!have && left->m_juncs) {
		for (i = 0; i < left->m_juncs->len; i++) {
			jun = (junction*) g_ptr_array_index(left->m_juncs, i);
			if ((jun->branch_tpl == right && jun->main_tpl == left)
					|| (jun->branch_tpl == left && jun->main_tpl == right)) {
				have = 1;
			}
		}
	}
	if (!have && left->b_juncs) {
		for (i = 0; i < left->b_juncs->len; i++) {
			jun = (junction*) g_ptr_array_index(left->b_juncs, i);
			if ((jun->branch_tpl == right && jun->main_tpl == left)
					|| (jun->branch_tpl == left && jun->main_tpl == right)) {
				have = 1;
			}
		}
	}
	return have;
}

void clean_junctions(GPtrArray *junctions) {
	int i = 0;
	junction *jun = NULL;
	for (i = 0; i < junctions->len; i++) {
		jun = (junction*) g_ptr_array_index(junctions, i);
		if (jun->status != 0) {
			g_ptr_array_remove_index_fast(junctions, i--);
			destroy_junction(jun);
		}
	}
}

/**
 * Filter out some junctions if:
 * 	- template length 0, connect to the same template at same locus
 *  - hanging branch template can be merged to the main template
 */
void filter_junctions(GPtrArray *junctions, GPtrArray *tpls, hash_table *ht) {
	uint32_t i = 0, j = 0, start_index = 0;
	junction *junc = NULL;
	tpl *t = NULL;
	GPtrArray *main_junctions = NULL;
	clean_junctions(junctions);
	while (start_index < junctions->len - 1) {
		//show_debug_msg(__func__, "======== Started at %d =======\n",
		//		start_index);
		junc = (junction*) g_ptr_array_index(junctions, start_index);
		t = junc->main_tpl;
		main_junctions = tpl_junctions(t, junctions, start_index, 1);
		filter_branches(main_junctions, ht->o->read_len);
		start_index += main_junctions->len;
		g_ptr_array_free(main_junctions, TRUE);
	}
	remove_dead_junctions(junctions);
	prune_short_branches(junctions, tpls, 10);
	remove_dead_junctions(junctions);
	g_ptr_array_sort(junctions, (GCompareFunc) cmp_junc_by_id);
}
