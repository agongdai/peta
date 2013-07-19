/*
 * ass.cpp
 *
 *  Created on: 09-Apr-2013
 *      Author: carl
 */

#include <unordered_map>
#include <stdint.h>
#include <glib.h>
#include <cstdio>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#define __STDC_FORMAT_MACROS
#include <inttypes.h>
#include "bwtaln.h"
#include "tpl.hpp"
#include "rnaseq.h"
#include "pechar.h"
#include "utils.h"
#include "peseq.h"
#include "kmers.hpp"
#include "ass.hpp"
#include "junction.hpp"
#include "graph.hpp"

using namespace std;

int kmer_ctg_id = 1;
int ins_size = 0;
int sd_ins_size = 0;
int kmer_n_threads = 0;
int kmer_len = 0;
uint32_t n_used_reads = 0;
char *kmer_out = NULL;
GPtrArray *branching_events = NULL;

GMutex *kmer_id_mutex;
struct timespec kmer_start_time, kmer_finish_time;

void add_a_junction(tpl *main_tpl, tpl *branch_tpl, uint64_t kmer, int locus,
		int ori, int weight) {
	// Indicating the templates are in-connect, cannot be reverse-complement
	main_tpl->in_connect = 1;
	branch_tpl->in_connect = 1;
	junction *new_j = new_junction(main_tpl, branch_tpl, kmer, locus, ori,
			weight);
	g_mutex_lock(kmer_id_mutex);
	g_ptr_array_add(branching_events, new_j);
	g_mutex_unlock(kmer_id_mutex);
	// Add the junction to the templates. For the 'existing connect' later.
	if (!main_tpl->m_juncs)
		main_tpl->m_juncs = g_ptr_array_sized_new(2);
	g_ptr_array_add(main_tpl->m_juncs, new_j);
	if (!branch_tpl->b_juncs)
		branch_tpl->b_juncs = g_ptr_array_sized_new(2);
	g_ptr_array_add(branch_tpl->b_juncs, new_j);
}

tpl *blank_tpl(uint64_t query_int, int kmer_len, int init_len, int ori) {
	bwa_seq_t *kmer = NULL;
	tpl *t = new_eg();
	g_mutex_lock(kmer_id_mutex);
	t->id = kmer_ctg_id++;
	g_mutex_unlock(kmer_id_mutex);
	// Get a copy of the kmer
	kmer = get_kmer_seq(query_int, kmer_len);
	if (ori)
		t->ctg = new_seq(kmer, init_len, 0);
	else
		t->ctg = new_seq(kmer, init_len, kmer->len - init_len);
	t->len = t->ctg->len;
	t->tid = atoi(kmer->name);
	t->start_kmer = query_int;
	bwa_free_read_seq(1, kmer);
	return t;
}

GPtrArray *hash_to_array(tpl_hash *all_tpls) {
	uint64_t id = 0;
	tpl *t = NULL;
	GPtrArray *tpls = g_ptr_array_sized_new(all_tpls->size());
	for (tpl_hash::iterator m = all_tpls->begin(); m != all_tpls->end(); ++m) {
		id = m->first;
		t = (tpl*) m->second;
		g_ptr_array_add(tpls, t);
		/**
		 show_debug_msg(__func__, "Tails of template %d\n", t->id);
		 p_ctg_seq(__func__, t->l_tail);
		 p_ctg_seq(__func__, t->r_tail);
		 **/
	}
	g_ptr_array_sort(tpls, (GCompareFunc) cmp_tpl_by_id);
	return tpls;
}

gint cmp_kmers_by_count(gpointer a, gpointer b) {
	kmer_counter *c_a = *((kmer_counter**) a);
	kmer_counter *c_b = *((kmer_counter**) b);
	return ((c_b->count) - c_a->count);
}

void mark_tpl_kmers_used(tpl *t, hash_map *hm, const int kmer_len,
		const int shift) {
	int i = 0;
	uint64_t query_int = 0;
	if (t->len < kmer_len)
		return;
	for (i = 0; i <= t->len - kmer_len; i++) {
		query_int = get_kmer_int(t->ctg->seq, i, 1, kmer_len);
		mark_kmer_used(query_int, hm, t->id, i + shift, t->len);
	}
}

void mark_tpl_kmers_fresh(tpl *t, hash_map *hm, const int kmer_len) {
	int i = 0;
	uint64_t query_int = 0;
	for (i = 0; i < t->len - kmer_len; i++) {
		query_int = get_kmer_int(t->ctg->seq, i, 1, kmer_len);
		mark_kmer_not_used(query_int, hm);
	}
}

/**
 * Mark reads on template as used.
 * If there are virtual tails, concat them as well.
 * To save time, just get reads on partial template
 */
void mark_reads_on_tpl(tpl *t, hash_map *hm, int start, int end) {
	int i = 0, j = 0, k = 0, ctg_len = 0, l_len = 0;
	uint64_t kmer_int = 0;
	bwa_seq_t *seq = NULL, *read = NULL, *ctg = NULL;
	GPtrArray *hits = NULL;
	int kmer_len = hm->o->k, read_len = hm->o->read_len;
	if (!t->reads)
		t->reads = g_ptr_array_sized_new(4);
	// Clear the whole reads array
	while (t->reads->len)
		g_ptr_array_remove_index_fast(t->reads, 0);
	ctg_len = t->len;
	if (t->l_tail) {
		ctg_len += t->l_tail->len;
		l_len = t->l_tail->len;
	}
	if (t->r_tail) {
		ctg_len += t->r_tail->len;
	}
	// show_debug_msg(__func__, "ctg_len: %d \n", ctg_len);
	if (ctg_len < read_len)
		return;
	// Concat the tails if there is any
	ctg = blank_seq(ctg_len);
	if (t->l_tail) {
		memcpy(ctg->seq, t->l_tail->seq, t->l_tail->len);
		ctg->len = t->l_tail->len;
	}
	memcpy(ctg->seq + ctg->len, t->ctg->seq, t->ctg->len);
	ctg->len += t->len;
	if (t->r_tail) {
		memcpy(ctg->seq + ctg->len, t->r_tail->seq, t->r_tail->len);
		ctg->len += t->r_tail->len;
	}
	set_rev_com(ctg);

	start -= read_len;
	start = start < 0 ? 0 : start;
	end = (end > ctg->len - read_len) ? (ctg->len - read_len) : end;
    /**
	 show_debug_msg(__func__, "Start~End: [%d, %d] \n", start, end);
	 p_ctg_seq(__func__, t->ctg);
	 p_ctg_seq(__func__, ctg);
    **/
	for (i = start; i <= end; i++) {
		seq = new_seq(ctg, read_len, i);
		hits = align_full_seq(seq, hm, N_MISMATCHES);
		for (j = 0; j < hits->len; j++) {
			read = (bwa_seq_t*) g_ptr_array_index(hits, j);
			// The locus could be <0
			add_read_to_tpl(t, read, i - l_len);
			//p_query(__func__, read);
			// Mark less frequent kmers on the read as used
			for (k = 0; k <= read->len - kmer_len; k++) {
				kmer_int = get_kmer_int(read->seq, k, 1, kmer_len);
				if (get_kmer_count(kmer_int, hm, 0) <= MIN_WEIGHT)
					mark_kmer_used(kmer_int, hm, t->id, i + k, t->len);
			}
		}
		g_ptr_array_free(hits, TRUE);
		bwa_free_read_seq(1, seq);
	}
	bwa_free_read_seq(1, ctg);
	g_ptr_array_sort(t->reads, (GCompareFunc) cmp_reads_by_name);
}

/**
 * Called after try_short_tpl_byte_ext.
 * If the tpl after short extension is too short but right connected,
 * 	borrow some sequence from the right connected tpl.
 */
void borrow_existing_seq(bwa_seq_t *short_seq, tpl *connected_tpl,
		const int shift, const int kmer_len, const int max_len, const int ori) {
	bwa_seq_t *tail = NULL;
	int borrow_len = 0, existing_ori = 0, existing_shift = 0;
	if (short_seq->len >= max_len)
		return;
	existing_ori = ori ? 0 : 1;
	existing_shift = existing_ori ? (shift + kmer_len - 1) : shift;
	tail = cut_tpl_tail(connected_tpl, existing_shift, max_len, existing_ori);
	borrow_len = max_len - short_seq->len;
	borrow_len = borrow_len > tail->len ? tail->len : borrow_len;
	// The space has be allocated to max_len in try_short_tpl_byte_ext already.
	if (existing_ori)
		memcpy(short_seq->seq + short_seq->len, tail->seq, borrow_len
				* sizeof(ubyte_t));
	else {
		memmove(short_seq->seq + borrow_len, short_seq->seq, short_seq->len);
		memcpy(short_seq->seq, tail->seq + (tail->len - borrow_len), borrow_len
				* sizeof(ubyte_t));
	}

	short_seq->len += borrow_len;
	free_read_seq(tail);
}

/**
 * Try to get the right connected
 */
tpl *get_connected_tpl(hash_map *hm, tpl_hash *all_tpls, uint64_t query_int,
		int *locus, const int ori) {
	int next_c = 0, eg_id = 0;
	uint64_t value = 0;
	tpl *right_tpl = NULL;

	next_c = next_char_by_kmers(hm, query_int, 0, ori);
	query_int = shift_bit(query_int, next_c, hm->o->k, ori);
	read_tpl_using_kmer(query_int, hm->hash, &eg_id, locus, &value);

	tpl_hash::iterator it = all_tpls->find(eg_id);
	if (it != all_tpls->end()) {
		right_tpl = (tpl*) it->second;
	}
	return right_tpl;
}

/**
 * Extend the short template as bwa_seq_t, max_len may be larger than 64
 * The query is updated. In case the extended seq is shorter than k.
 */
bwa_seq_t *try_short_tpl_byte_ext(hash_map *hm, uint64_t *query, int first_c,
		const int ori, const int max_len) {
	bwa_seq_t *branch_seq = NULL;
	int i = 0, c = 0;
	branch_seq = blank_seq(max_len);
	branch_seq->seq[0] = first_c;
	branch_seq->len = 1;
	for (i = 1; i < max_len; i++) {
		c = next_char_by_kmers(hm, *query, 1, ori);
		if (c == -1)
			break;
		branch_seq->seq[i] = c;
		branch_seq->len = i + 1;
		*query = shift_bit(*query, c, hm->o->k, ori);
		//p_ctg_seq(__func__, branch_seq);
	}
	if (ori)
		seq_reverse(branch_seq->len, branch_seq->seq, 0);
	set_rev_com(branch_seq);
	return branch_seq;
}

/**
 * Validate a branching is valid or not.
 * Need to tailor and extend the two templates if they are not long enough
 **/
int val_branching(hash_map *hm, tpl *main_tpl, tpl_hash *all_tpls,
		const int shift, uint64_t query_int, int c, const int ori, int *weight,
		const int max_len) {
	uint64_t query = 0;
	int n_reads = 0, existing_shift = 0;
	bwa_seq_t *branch_seq = NULL, *main_seq = NULL;
	tpl *connected_tpl = NULL;
	GPtrArray *reads = NULL;

	query = shift_bit(query_int, c, hm->o->k, ori);
	branch_seq = try_short_tpl_byte_ext(hm, &query, c, ori, max_len / 2);
	connected_tpl
			= get_connected_tpl(hm, all_tpls, query, &existing_shift, ori);
	// If the branch is right connected.
	if (connected_tpl && branch_seq->len < max_len / 2) {
		borrow_existing_seq(branch_seq, connected_tpl, existing_shift,
				hm->o->k, max_len / 2, ori);
	}

	if (ori) {
		main_seq = cut_tpl_tail(main_tpl, shift, max_len / 2, ori);
		reads = find_junc_reads(hm, branch_seq, main_seq, max_len, weight);
	} else {
		main_seq = cut_tpl_tail(main_tpl, shift + hm->o->k, max_len / 2, ori);
		reads = find_junc_reads(hm, main_seq, branch_seq, max_len, weight);
	}
	n_reads = reads->len;
	g_ptr_array_free(reads, TRUE);
	//p_ctg_seq("MAIN", main_tpl->ctg);
	//p_ctg_seq("SUB", main_seq);
	//p_ctg_seq("BRANCH", branch_seq);
	//show_debug_msg(__func__, "N_READS: %d \n", n_reads);
	bwa_free_read_seq(1, main_seq);
	bwa_free_read_seq(1, branch_seq);
	return n_reads;
}

/**
 * Validate the junction by checking mate pairs.
 * Depending on ori and con_pos, only partial reads on the main template are counted
 */
int vld_junc_by_mates(tpl *main_tpl, tpl *branch_tpl, GPtrArray *junc_reads,
		hash_map *hm, const int con_pos, const int ori) {
	int start = con_pos, end = main_tpl->len;
	int is_valid = 0, seg_len = 0;
	GPtrArray *tmp = NULL;
	return 1;
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
		end = con_pos - hm->o->read_len;
	}

	mark_reads_on_tpl(main_tpl, hm, start, end);
	mark_reads_on_tpl(branch_tpl, hm, 0, branch_tpl->len);

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

	is_valid = vld_tpl_mates(branch_tpl, main_tpl, start, end, MIN_PAIRS);
	if (!is_valid) {
		g_ptr_array_sort(junc_reads, (GCompareFunc) cmp_reads_by_name);
		/**
		 show_debug_msg(__func__, "Tag 1 \n");
		 p_readarray(main_tpl->reads, 1);
		 show_debug_msg(__func__, "Tag 2\n");
		 p_readarray(junc_reads, 1);
		 **/
		p_readarray(junc_reads, 1);
		// Maybe exon shorter than read length, the mates located at the junction
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
 * To connect to templates with length smaller than k (25 by default)
 * The input template short_tpl is supposed to be shorter than k
 */
tpl *connect_to_small_tpl(hash_map *hm, uint64_t query_int, tpl *branch_tpl,
		tpl *short_tpl, int *parent_locus, int *borrow_bases, int ori) {
	uint32_t i = 0;
	uint64_t kmer_int = 0;
	int max_len = hm->o->k - 1;
	junction *j = NULL, *left_j = NULL, *right_j = NULL;
	tpl *left = NULL, *right = NULL, *parent = NULL;
	bwa_seq_t *left_seq = NULL, *right_seq = NULL, *junc_seq = NULL;
	GPtrArray *branch_juncs = short_tpl->b_juncs;
	if (!branch_juncs || branch_juncs->len != 2) {
		return NULL;
	}
	show_debug_msg(__func__, "Branch: [%d, %d] \n", branch_tpl->id,
			branch_tpl->len);
	show_debug_msg(__func__, "Short: [%d, %d] \n", short_tpl->id,
			short_tpl->len);
	// Assign the left template  and right template
	for (i = 0; i < branch_juncs->len; i++) {
		j = (junction*) g_ptr_array_index(branch_juncs, i);
		if (j->ori == 0) {
			left_j = j;
			p_junction(left_j);
		} else {
			right_j = j;
			p_junction(right_j);
		}
	}
	if (!left_j || !right_j)
		return NULL;

	// Concatenate the three short sequences, to find where the branch_tpl should be connected to
	//p_ctg_seq(__func__, left_j->main_tpl->ctg);
	//show_debug_msg(__func__, "Locus: %d\n", left_j->locus);
	left_seq = cut_tpl_tail(left_j->main_tpl, left_j->locus, max_len, 0);
	right_seq = cut_tpl_tail(right_j->main_tpl, right_j->locus, max_len, 1);
	junc_seq = blank_seq(left_seq->len + short_tpl->len + right_seq->len);
	memcpy(junc_seq->seq, left_seq->seq, left_seq->len);
	memcpy(junc_seq->seq + left_seq->len, short_tpl->ctg->seq, short_tpl->len);
	memcpy(junc_seq->seq + left_seq->len + short_tpl->len, right_seq->seq,
			right_seq->len);
	junc_seq->len = left_seq->len + short_tpl->len + right_seq->len;
	set_rev_com(junc_seq);
	//p_ctg_seq("JUNCTION", junc_seq);

	for (i = 0; i <= junc_seq->len - hm->o->k; i++) {
		kmer_int = get_kmer_int(junc_seq->seq, i, 1, hm->o->k);

		/**
		show_debug_msg(__func__, "Debug %d: left_len: %d; right_len: %d \n", i,
				left_seq->len, right_seq->len);
		bwa_seq_t *debug = get_kmer_seq(kmer_int, 25);
		p_query(__func__, debug);
		bwa_free_read_seq(1, debug);
		**/

		if (kmer_int == query_int) {
			bwa_free_read_seq(1, junc_seq);
			if (ori) {
				if (i >= left_seq->len) {
					*parent_locus = i - left_seq->len;
					parent = short_tpl;
				} else {
					*borrow_bases = hm->o->k - 1 - i;
					*parent_locus = left_j->locus - *borrow_bases;
					parent = left_j->main_tpl;
				}
			} else {
				if (i > left_seq->len) {
					*parent_locus = i - left_seq->len;
					parent = short_tpl;
				} else {
					*borrow_bases = hm->o->k - short_tpl->len - (left_seq->len
							- i) - 1;
					*parent_locus = right_j->locus + *borrow_bases;
					parent = right_j->main_tpl;
				}
			}
			break;
		}
	}
	bwa_free_read_seq(1, left_seq);
	bwa_free_read_seq(1, right_seq);
	return parent;
}

/**
 * If the template reaches some kmer which is used, stop the extension and add a branching event
 */
int connect(tpl *branch, hash_map *hm, tpl_hash *all_tpls, uint64_t query_int,
		const int ori) {
	int *counters = NULL, locus = 0, i = 0, connected = 0, eg_id = 0,
			valid = 0, weight = 0, con_pos = 0, exist_ori = 0,
			parent_locus = 0, borrow_bases = 0, on_main = 0, to_destroy = 0;
	int loop_len = 0, j = 0;
	uint64_t value = 0, query_copy = query_int;
	tpl *existing = NULL, *parent_existing = NULL;
	GPtrArray *junc_reads = NULL;
	junction *right_junc = NULL;

	counters = count_next_kmers(hm, query_int, 0, ori);

	/**
	 p_ctg_seq(__func__, branch->ctg);
	 bwa_seq_t *debug = get_kmer_seq(query_int, 25);
	 p_query(__func__, debug);
	 bwa_free_read_seq(1, debug);
	 show_debug_msg(__func__, "Counters: %d,%d,%d,%d\n", counters[0],
	 counters[1], counters[2], counters[3]);
	 **/

	// In case connecting to the template itself
	mark_tpl_kmers_used(branch, hm, hm->o->k, 0);
	for (i = 0; i < 4; i++) {
		if (counters[i] < MIN_WEIGHT)
			continue;
		show_debug_msg(__func__,
				"---------- Connecting to existing, ori %d ----------\n", ori);

		query_int = shift_bit(query_copy, i, hm->o->k, ori);
		read_tpl_using_kmer(query_int, hm->hash, &eg_id, &locus, &value);

		tpl_hash::iterator it = all_tpls->find(eg_id);
		if (it != all_tpls->end()) {
			existing = (tpl*) it->second;
			// It happens when 'existing' and 'branch' are the same.
			// And the same template has been trimmed before.

			/**
			 show_debug_msg(__func__, "Existing tpl %d: %d\n", existing->id,
			 existing->len);
			 show_debug_msg(__func__,
			 "Locus: %d; existing->len: %d; hm->o->k: %d\n", locus,
			 existing->len, hm->o->k);
			 p_ctg_seq("MAIN", existing->ctg);
			 p_ctg_seq("BRANCH", branch->ctg);
			 **/

			//if (locus > existing->len - hm->o->k)
			//	continue;
			con_pos = ori ? (locus + 1) : (locus + hm->o->k - 1);
			exist_ori = ori ? 0 : 1;

			/**
			bwa_seq_t *debug = get_kmer_seq(query_int, 25);
			p_query(__func__, debug);
			bwa_free_read_seq(1, debug);
			show_debug_msg(__func__, "connect pos: %d; locus: %d \n", con_pos,
					locus);
			if (branch->b_juncs && branch->b_juncs->len > 0) {
				right_junc = (junction*) g_ptr_array_index(branch->b_juncs, 0);
				p_junction(right_junc);
			}
			**/

			// If the existing is too short, maybe need to connect to its connector instead.
			if (existing->len < hm->o->k) {
				parent_existing = connect_to_small_tpl(hm, query_int, branch,
						existing, &parent_locus, &borrow_bases, ori);
				if (parent_existing) {
					/**
					show_debug_msg(
							__func__,
							"TAG: Existing changed from [%d, %d] to [%d, %d]\n",
							existing->id, existing->len, parent_existing->id,
							parent_existing->len);
					show_debug_msg(__func__,
							"TAG: Connect position changed from %d to %d\n",
							con_pos, parent_locus);
					**/
					existing = parent_existing;
					con_pos = parent_locus;
				}
			}

			// If the branch is can be merged to main template, skip
			// Only check when it is extending to the left and
			//		its right end is not connected to anywhere.
			if (branch->len < 100 && branch->r_tail) {
				/**
				 show_debug_msg(__func__, "Main template: %d; pos: %d \n",
				 existing->id, con_pos);
				 show_debug_msg(__func__, "Branch template: %d; ori: %d \n",
				 branch->id, exist_ori);
				 p_ctg_seq("Right Tail", branch->r_tail);
				 **/

				right_junc = (junction*) g_ptr_array_index(branch->b_juncs, 0);
				if (get_abs(right_junc->locus - con_pos) <= N_MISMATCHES
						&& branch->len <= N_MISMATCHES) {
					branch->alive = 0;
					g_ptr_array_remove_index_fast(branching_events,
							branching_events->len - 1);
					break;
				}
				loop_len = get_abs(branch->len - (right_junc->locus - con_pos));
				if (loop_len <= N_MISMATCHES) {
					on_main = branch_on_main(existing->ctg, branch->ctg,
							con_pos, (branch->len / hm->o->k + 2) * 3,
							exist_ori);
					valid = on_main ? 0 : 1;
					if (!valid) {
						// Mark it as 'dead', will be destroyed in kmer_ext_thread.
						branch->alive = 0;
						// Remove its right junction from the global junctions
						// @TODO: here it works only if it is single thread
						g_ptr_array_remove_index_fast(branching_events,
								branching_events->len - 1);
						break;
					}
				}
			}

			// If no enough junction reads, skip
			junc_reads = find_junc_reads_w_tails(hm, existing, branch, con_pos,
					(hm->o->read_len - SHORT_BRANCH_SHIFT) * 2, ori, &weight);
			valid = (junc_reads->len >= MIN_JUNCTION_READS) ? 1 : 0;
			//p_readarray(junc_reads, 0);
			if (!valid) {
				show_debug_msg(__func__, "No enough junction reads\n");
				continue;
			}

			// If no pairs, skip
			valid = vld_junc_by_mates(existing, branch, junc_reads, hm,
					con_pos, ori);
			g_ptr_array_free(junc_reads, TRUE);
			if (!valid) {
				show_debug_msg(__func__, "Not passed the pair validation\n");
				continue;
			}

			// Trim the branch
			if (branch->len < hm->o->k) {
				if (exist_ori)
					con_pos -= branch->len;
				else
					con_pos = locus + branch->len + 1;
				branch->len = 0;
				branch->ctg->len = 0;
				set_rev_com(branch->ctg);
			} else {
				/**
				p_ctg_seq("EXISTING", existing->ctg);
				p_ctg_seq("BRANCH", branch->ctg);
				show_debug_msg(__func__, "Borrow bases: %d\n", borrow_bases);
				show_debug_msg(__func__, "Branch length: %d \n", branch->len);
				show_debug_msg(__func__, "ORI: %d \n", ori);
				**/
				// Make the branch not sharing a 24-mer with the main
				if (borrow_bases) {
					if (exist_ori) {
						con_pos -= borrow_bases;
						branch->len -= borrow_bases;
					} else {
						con_pos += borrow_bases;
						memmove(branch->ctg->seq, branch->ctg->seq
								+ borrow_bases, sizeof(ubyte_t) * (branch->len
								- borrow_bases));
						branch->len -= borrow_bases;
					}
					//show_debug_msg(__func__, "Connecting position: %d\n",
					//		con_pos);
				} else {
					if (exist_ori) {
						con_pos -= (hm->o->k - 1);
					} else {
						con_pos += (hm->o->k - 1);
						memmove(branch->ctg->seq, branch->ctg->seq + (hm->o->k
								- 1), sizeof(ubyte_t) * (branch->len
								- (hm->o->k - 1)));
					}
					branch->len -= (hm->o->k - 1);
				}

				branch->ctg->len = branch->len;
				//p_ctg_seq("TRUNCATED", branch->ctg);
				set_rev_com(branch->ctg);
				//if (borrow_bases && branch->id > 1472)
				//	exit(1);
			}
			// If there is a small loop, erase it.
			if (ori && branch->b_juncs && branch->b_juncs->len == 1) {
				right_junc = (junction*) g_ptr_array_index(branch->b_juncs, 0);
				if (right_junc->main_tpl == existing && right_junc->ori == 1) {
					loop_len = con_pos - right_junc->locus;
					show_debug_msg(__func__, "Erasing small loop... \n");
					p_junction(right_junc);
					if (loop_len > 0 && loop_len < hm->o->k) {
						// Copy the small loop to the head of the branch, adjust the connecting position
						for (j = right_junc->locus + loop_len - 1; j
								>= right_junc->locus; j--) {
							ext_con(branch->ctg, existing->ctg->seq[j], 1);
						}
						branch->len = branch->ctg->len;
						con_pos = right_junc->locus;
					}
				}
			}

			show_debug_msg(__func__,
					"Connect existing [%d, %d] to [%d, %d] at %d. \n",
					branch->id, branch->len, existing->id, existing->len,
					con_pos);
			set_tail(branch, existing, con_pos, hm->o->read_len
					- SHORT_BRANCH_SHIFT, exist_ori);
			//p_ctg_seq("Right tail", branch->r_tail);
			//p_ctg_seq("Left  tail", branch->l_tail);
			add_a_junction(existing, branch, query_int, con_pos, exist_ori,
					weight);
			connected = 1;
		}
	}
	free(counters);
	return connected;
}

/**
 * If next kmer is used by another template, try to connect to it.
 * Return: 0: cannot be connected
 *         1: forward connected
 *         2: reverse complement connected
 **/
int existing_connect(tpl *branch, hash_map *hm, tpl_hash *all_tpls,
		uint64_t query_int, int ori) {
	int connected = 0, rev_ori = 0;
	uint64_t rev_query = 0;
	// During extension, the sequence is actually reversed, here reverse back temp
	if (ori)
		seq_reverse(branch->len, branch->ctg->seq, 0);
	set_rev_com(branch->ctg);
	connected = connect(branch, hm, all_tpls, query_int, ori);
	// Try the reverse complement of the branch and connect
	//	, if there is no other template connecting to it currently
	if (!connected && !branch->in_connect) {
		//show_debug_msg(__func__, "ATTENTION: going to connect reverse complement\n");
		switch_fr(branch->ctg);
		rev_ori = ori ? 0 : 1;
		rev_query = rev_comp_kmer(query_int, hm->o->k);
		connected = connect(branch, hm, all_tpls, rev_query, rev_ori);
		// If connected, no need to reverse back, because the extending will be always stopped
		if (!connected)
			switch_fr(branch->ctg);
		else {
			connected = 2; // Indicates that its rev-comp connectes to existing template
			//show_debug_msg(__func__, "ATTENTION: connected reverse complement\n");
		}
	}
	if (ori)
		seq_reverse(branch->len, branch->ctg->seq, 0);
	return connected;
}

/**
 * Extend the short template as int64
 * max_len: should be smaller than 64
 */
uint64_t try_short_tpl_ext(hash_map *hm, uint64_t query, const int ori,
		const int max_len) {
	uint64_t branch_template = 0;
	int i = 0, c = 0;
	for (i = 0; i < max_len; i++) {
		c = next_char_by_kmers(hm, query, 1, ori);
		if (c == -1)
			return LOWER_ONES_32;
		else {
			branch_template = shift_bit(branch_template, c, hm->o->k, ori);
			query = shift_bit(query, c, hm->o->k, ori);
			//show_debug_msg(__func__, "Branch int: %"ID64"\n", branch_template);
		}
	}
	return branch_template;
}

/**
 * Trim the kmers whose support is only 1.
 */
void trim_weak_tails(tpl *t, hash_map *hm, const int ori) {
	uint64_t query_int = 0;
	int i = 0, support = 0, kmer_len = hm->o->k, shift = 0, sum_sup = 0;
	if (ori) {
		for (i = 0; i <= t->len - kmer_len; i++) {
			query_int = get_kmer_int(t->ctg->seq, i, 1, kmer_len);
			support = get_kmer_count(query_int, hm, 0);
			if (support <= 1) {
				shift++;
				sum_sup += support;
			} else
				break;
		}
		if (shift > 0) {
			//p_ctg_seq(__func__, t->ctg);
			memcpy(t->ctg->seq, t->ctg->seq + shift, sizeof(ubyte_t) * (t->len
					- shift));
			t->ctg->len = t->len - shift;
			t->len = t->ctg->len;
			t->kmer_freq -= sum_sup;
			set_rev_com(t->ctg);
			//p_ctg_seq(__func__, t->ctg);
		}
	} else {
		//p_ctg_seq(__func__, t->ctg);
		for (i = t->len - kmer_len; i >= 0; i--) {
			query_int = get_kmer_int(t->ctg->seq, i, 1, kmer_len);
			support = get_kmer_count(query_int, hm, 0);
			if (support <= 1) {
				sum_sup += support;
				t->ctg->len--;
			} else
				break;
		}
		t->kmer_freq -= sum_sup;
		t->len = t->ctg->len;
		set_rev_com(t->ctg);
		//p_ctg_seq(__func__, t->ctg);
	}
}

/**
 * For an existing tpl, try to assemble the branches
 */
void kmer_ext_branch(tpl *t, hash_map *hm, tpl_hash *all_tpls, const int ori) {
	int *counters = NULL, weight = 0, c = 0, i = 0, j = 0, con_pos = 0;
	int kmer_len = hm->o->k, branch_is_valid = 0, max_freq = 0, con_existing =
			0, valid = 0;
	uint64_t query_int = 0, branch_query = 0;
	tpl *branch = NULL;
	if (t->len <= kmer_len)
		return;
	//show_debug_msg(__func__, "^^^^^ Branching [%d, %d] to ori %d ^^^^^^\n",
	//		t->id, t->len, ori);
	for (i = 0; i < t->len - kmer_len; i++) {
		query_int = get_kmer_int(t->ctg->seq, i, 1, kmer_len);
		counters = count_next_kmers(hm, query_int, 1, ori);

		for (j = 0; j < 4; j++) {
			c = counters[j];
			max_freq = next_char_max_freq(hm, query_int, 0, ori);
			// At the branching point, the frequency should be at least BRANCH_THRE of the max frequency.
			if (c < MIN_WEIGHT)
				continue;

			/**
			 bwa_seq_t *debug = get_kmer_seq(query_int, 25);
			 p_query(__func__, debug);
			 bwa_free_read_seq(1, debug);
			 show_debug_msg(__func__,
			 "[%d, %d] %d pos counters ori %d: [%d, %d, %d, %d]\n",
			 t->id, t->len, i, ori, counters[0], counters[1],
			 counters[2], counters[3]);
			 show_debug_msg(__func__, "Counters[]: %d, Max frequency: %d \n",
			 counters[j], max_freq);
			 **/

			weight = 0;
			branch_query = shift_bit(query_int, j, kmer_len, ori);
			branch_is_valid = val_branching(hm, t, all_tpls, i, query_int, j,
					ori, &weight, (hm->o->read_len - SHORT_BRANCH_SHIFT) * 2);
			if (!branch_is_valid || weight < MIN_WEIGHT) {
				continue;
			}
			show_debug_msg(__func__, "Weight: %d \n", weight);
			branch = blank_tpl(branch_query, kmer_len, 1, ori);
			branch->kmer_freq = get_kmer_count(branch_query, hm, 1);

			con_pos = ori ? i : i + kmer_len;
			set_tail(branch, t, con_pos, hm->o->read_len - SHORT_BRANCH_SHIFT,
					ori);
			// Insert first, in case it connects to itself during extension
			g_mutex_lock(kmer_id_mutex);
			all_tpls->insert(make_pair<int, tpl*> ((int) branch->id,
					(tpl*) branch));
			g_mutex_unlock(kmer_id_mutex);
			con_existing
					= kmer_ext_tpl(branch, branch_query, hm, all_tpls, ori);
			// If the branch can be merged into main template, erase the branch
			if (branch_on_main(t->ctg, branch->ctg, con_pos, (branch->len
					/ hm->o->k + 2) * 3, ori)) {
				show_debug_msg(__func__, "Branch is not valid: [%d, %d] \n",
						branch->id, branch->len);
				// Remove the branch from the global hash table
				g_mutex_lock(kmer_id_mutex);
				all_tpls->erase(branch->id);
				g_mutex_unlock(kmer_id_mutex);
				// Although the branch will be removed, its kmers are marked as used by the main template.
				branch->id = t->id;
				// The kmer frequencies are added to the main tpl now.
				t->kmer_freq += branch->kmer_freq;
				if (ori) {
					mark_tpl_kmers_used(branch, hm, kmer_len, t->len
							- branch->len);
				} else
					mark_tpl_kmers_used(branch, hm, kmer_len, t->len);
				// Do not destroy first, because it may be right connected to some tpl. Mark it as not alive
				// But free most of the memory it occupies, since there are many invalid branches.
				// Destory not-alive junctions and tpls in function clean_junctions.
				free_eg_seq(branch);
				branch->alive = 0;
			} else {
				//				valid = vld_junc_by_mates(t, branch, hm, con_pos, ori);
				//				if (valid) {
				add_a_junction(t, branch, query_int, con_pos, ori, weight);
				mark_tpl_kmers_used(branch, hm, kmer_len, 0);
				upd_tpl_jun_locus(branch, branching_events, kmer_len);
				// Try to extend branches of current branch
				kmer_ext_branch(branch, hm, all_tpls, 0);
				kmer_ext_branch(branch, hm, all_tpls, 1);
				//				}
			}
		}
		free(counters);
		counters = NULL;
	}
	free(counters);
}

int kmer_ext_tpl(tpl *t, uint64_t query_int, hash_map *hm, tpl_hash *all_tpls,
		const int ori) {
	int max_c = 0, *counters = NULL, weight = 0, con_existing = 0;
	int max_c_all = 0, *counters_all = NULL;

	show_debug_msg(__func__, "------ Started extending tpl %d to ori %d...\n",
			t->id, ori);
	bwa_seq_t *debug = get_kmer_seq(query_int, 25);
	p_query(__func__, debug);
	bwa_free_read_seq(1, debug);
	if (ori)
		seq_reverse(t->len, t->ctg->seq, 0);
	while (1) {
		weight = 0;
		counters = count_next_kmers(hm, query_int, 1, ori);
		max_c = get_max_index(counters);
		counters_all = count_next_kmers(hm, query_int, 0, ori);
		max_c_all = get_max_index(counters_all);

		/**
		 if (t->id == 3) {
		 if (ori)
		 seq_reverse(t->len, t->ctg->seq, 0);
		 bwa_seq_t *debug = get_kmer_seq(query_int, 25);
		 p_query(__func__, debug);
		 bwa_free_read_seq(1, debug);
		 show_debug_msg(__func__, "Next chars: [%d, %d, %d, %d] \n",
		 counters[0], counters[1], counters[2], counters[3]);
		 show_debug_msg(__func__,
		 "Ori %d, tpl %d, length %d, Next char: %d \n", ori, t->id,
		 t->len, max_c);
		 p_ctg_seq("Contig", t->ctg);
		 if (ori)
		 seq_reverse(t->len, t->ctg->seq, 0);
		 }
		 **/

		// If the max direction is taken already, connect to it.
		if (max_c_all != max_c) {
			if (existing_connect(t, hm, all_tpls, query_int, ori)) {
				con_existing = 1;
				free(counters_all);
				break;
			}
		}
		// Maybe marked as not alive in existing_connect
		if (!t->alive) {
			free(counters);
			break;
		}
		if (max_c == -1) {
			if (!existing_connect(t, hm, all_tpls, query_int, ori))
				show_debug_msg(__func__, "[%d, %d] No hits, stop here. \n",
						t->id, t->len);
			else {
				con_existing = 1;
			}
			free(counters);
			break;
		}
		t->kmer_freq += counters[max_c];
		free(counters);
		free(counters_all);
		ext_con(t->ctg, max_c, 0);
		t->len = t->ctg->len;
		// In case the template has repeats, so mark kmers used TEMPERARIALY. The locus is incorrect.
		mark_kmer_used(query_int, hm, t->id, t->len, t->len);
		query_int = shift_bit(query_int, max_c, hm->o->k, ori);
		if (t->len % 100 == 0)
			show_debug_msg(__func__, "Ori %d, tpl %d, length %d \n", ori,
					t->id, t->len);
	}
	if (ori)
		seq_reverse(t->len, t->ctg->seq, 0);
	return con_existing;
}

/**
 * Extend a kmer
 */
void *kmer_ext_thread(gpointer data, gpointer thread_params) {
	tpl *t = NULL;
	int round_1_len = 0, round_2_len = 0, connected = 0, ori = 1, flag = 0;
	uint64_t kmer_int = 0, rev_kmer_int = 0, c = 0;
	map_opt *opt = NULL;
	kmer_counter *counter = NULL;
	bwa_seq_t *query = NULL;
	kmer_t_meta *params = (kmer_t_meta*) thread_params;
	tpl_hash *all_tpls = params->all_tpls;
	opt = params->hm->o;

	counter = (kmer_counter*) data;
	kmer_int = counter->kmer;
	rev_kmer_int = rev_comp_kmer(kmer_int, opt->k);
	c = get_kmer_rf_count(kmer_int, params->hm, 0);
	query = get_kmer_seq(kmer_int, opt->k);

	/**
	 show_debug_msg("TAG",
	 "============= %" ID64 ": %" ID64 " ============ \n", kmer_int, c);
	 p_query(__func__, query);
	 **/

	if (kmer_int < 64 || rev_kmer_int < 64 || c <= 1 || is_repetitive_q(query)
			|| is_biased_q(query) || kmer_is_used(kmer_int, params->hm)) {
		bwa_free_read_seq(1, query);
		return NULL;
	}
	bwa_free_read_seq(1, query);
	show_debug_msg(__func__,
			"============= %" ID64 ": %" ID64 " ============ \n", kmer_int, c);
	t = blank_tpl(kmer_int, opt->k, opt->k, 0);
	t->kmer_freq = get_kmer_count(kmer_int, params->hm, 1);
	// Insert first, in case it connects to itself during extension
	g_mutex_lock(kmer_id_mutex);
	all_tpls->insert(make_pair<int, tpl*> ((int) t->id, (tpl*) t));
	g_mutex_unlock(kmer_id_mutex);

	connected = kmer_ext_tpl(t, kmer_int, params->hm, all_tpls, 0);
	round_1_len = t->len;
	show_debug_msg(__func__, "tpl %d with length: %d \n", t->id, t->len);

	// Its reverse complement is connected to an existing template
	if (connected == 2) {
		ori = 0;
		kmer_int = rev_kmer_int;
	}
	connected |= kmer_ext_tpl(t, kmer_int, params->hm, all_tpls, ori);
	show_debug_msg(__func__, "tpl %d with length: %d \n", t->id, t->len);

	/**
	 if (!connected && t->len - round_1_len > 2) {
	 kmer_int = get_kmer_int(t->ctg->seq, t->len - opt->k, 1, opt->k);
	 round_2_len = t->len;
	 flag = kmer_ext_tpl(t, kmer_int, params->hm, all_tpls, 0);
	 show_debug_msg(__func__, "tpl %d with length: %d \n", t->id, t->len);
	 if (t->len - round_2_len > 2) {
	 kmer_int = get_kmer_int(t->ctg->seq, 0, 1, opt->k);
	 connected |= kmer_ext_tpl(t, kmer_int, params->hm, all_tpls, 1);
	 show_debug_msg(__func__, "tpl %d with length: %d \n", t->id,
	 t->len);
	 }
	 }
	 **/

	if (!t->alive || (!connected && t->len <= opt->k)) {
		g_mutex_lock(kmer_id_mutex);
		all_tpls->erase(t->id);
		g_mutex_unlock(kmer_id_mutex);
		destroy_eg(t);
	} else {
		mark_tpl_kmers_used(t, params->hm, opt->k, 0);
		mark_reads_on_tpl(t, params->hm, 0, t->len);
		upd_tpl_jun_locus(t, branching_events, opt->k);
		t->start_kmer = *((uint64_t*) data);
	}
	return NULL;
}

void kmer_threads(kmer_t_meta *params) {
	GThreadPool *thread_pool = NULL;
	uint64_t i = 0, *freq = 0, count = 0;
	hash_map *hm = params->hm;
	mer_hash *hash = hm->hash;
	kmer_counter *counter = NULL;
	GPtrArray *start_kmers = g_ptr_array_sized_new(BUFSIZ);

	for (mer_hash::iterator it = hash->begin(); it != hash->end(); ++it) {
		freq = it->second;
		count = freq[0];
		count <<= 40;
		count >>= 40;
		if (count > 1 && it->first > 64 && rev_comp_kmer(it->first, hm->o->k)
				> 64) {
			counter = (kmer_counter*) malloc(sizeof(kmer_counter));
			counter->kmer = it->first;
			counter->count = get_kmer_rf_count(it->first, hm, 0);
			g_ptr_array_add(start_kmers, counter);
		}
	}
	show_msg(__func__, "Sorting %d initial kmers... \n", start_kmers->len);
	g_ptr_array_sort(start_kmers, (GCompareFunc) cmp_kmers_by_count);
	show_msg(__func__, "Extending by kmers...\n");
	thread_pool = g_thread_pool_new((GFunc) kmer_ext_thread, params, 1, TRUE,
			NULL);
	for (i = 0; i < start_kmers->len; i++) {
		if (i % 100000 == 0)
		show_debug_msg(__func__, "Extending %" ID64 "-th kmer... \n ", i);
		counter = (kmer_counter*) g_ptr_array_index(start_kmers, i);
		kmer_ext_thread(counter, params);
		free(counter);
		//g_thread_pool_push(thread_pool, (gpointer) counter, NULL);
		//if (i >= 100100)
		//	break;
	}
	g_ptr_array_free(start_kmers, TRUE);
	g_thread_pool_free(thread_pool, 0, 1);
}

void test_kmer_ext(kmer_t_meta *params) {
	uint64_t kmer_int = 957139786396488;
	int round_1_len = 0, round_2_len = 0;
	tpl *t = new_eg();
	FILE *contigs = NULL;
	GPtrArray *all_tpls = NULL;
	map_opt *opt = params->hm->o;
	tpl_hash *tpl_h = params->all_tpls;

	t = blank_tpl(kmer_int, 25, 25, 0);
	tpl_h->insert(make_pair<int, tpl*> ((int) t->id, (tpl*) t));
	kmer_ext_tpl(t, kmer_int, params->hm, params->all_tpls, 0);
	round_1_len = t->len;
	show_debug_msg(__func__, "Extending to the right: [%d, %d]... \n", t->id,
			t->len);
	kmer_ext_tpl(t, kmer_int, params->hm, params->all_tpls, 1);
	show_debug_msg(__func__, "Done Extending to the right: [%d, %d]... \n",
			t->id, t->len);

	mark_tpl_kmers_used(t, params->hm, opt->k, 0);
	upd_tpl_jun_locus(t, branching_events, opt->k);
	//kmer_ext_branch(t, params->hm, all_tpls, 0);
	//kmer_ext_branch(t, params->hm, all_tpls, 1);

	kmer_int = 502458011902035;
	t = blank_tpl(kmer_int, 25, 25, 0);
	tpl_h->insert(make_pair<int, tpl*> ((int) t->id, (tpl*) t));
	kmer_ext_tpl(t, kmer_int, params->hm, params->all_tpls, 0);
	round_1_len = t->len;
	show_debug_msg(__func__, "Extending to the right: [%d, %d]... \n", t->id,
			t->len);
	kmer_ext_tpl(t, kmer_int, params->hm, params->all_tpls, 1);
	show_debug_msg(__func__, "Done Extending to the right: [%d, %d]... \n",
			t->id, t->len);

	mark_tpl_kmers_used(t, params->hm, opt->k, 0);
	upd_tpl_jun_locus(t, branching_events, opt->k);
	//	kmer_ext_branch(t, params->hm, all_tpls, 0);
	//	kmer_ext_branch(t, params->hm, all_tpls, 1);

	all_tpls = hash_to_array(tpl_h);
	contigs = xopen(get_output_file("single.fa", kmer_out), "w");
	save_tpls(all_tpls, contigs, 0, 0, 0);
	fflush(contigs);
	fclose(contigs);

	clean_junctions(branching_events);
	g_ptr_array_sort(branching_events, (GCompareFunc) cmp_junc_by_id);
	store_features(get_output_file("single.junctions", kmer_out),
			branching_events, all_tpls);
}

/**
 * Start branching after the frequent kmers are consumed already.
 */
void start_branching(tpl_hash *all_tpls, kmer_t_meta *params) {
	uint64_t id = 0;
	tpl *t = NULL;
	for (tpl_hash::iterator m = all_tpls->begin(); m != all_tpls->end(); ++m) {
		id = m->first;
		t = (tpl*) m->second;
		kmer_ext_branch(t, params->hm, all_tpls, 0);
		kmer_ext_branch(t, params->hm, all_tpls, 1);
	}
}

void ext_by_kmers_core(char *lib_file, const char *solid_file) {
	hash_map *hm = NULL;
	mer_hash map;
	FILE *contigs = NULL;
	kmer_t_meta *params = (kmer_t_meta*) calloc(1, sizeof(kmer_t_meta));
	tpl_hash all_tpls;
	GPtrArray *kmer_tpls = NULL;

	show_msg(__func__, "Library: %s \n", lib_file);
	show_msg(__func__, "Solid Reads: %s \n", solid_file);

	hm = load_hash_map(lib_file, 0, map);
	hm->hash = &map;

	clock_gettime(CLOCK_MONOTONIC, &kmer_finish_time);
	show_msg(__func__, "Done loading kmers hash map: %.2f sec\n",
			(float) (kmer_finish_time.tv_sec - kmer_start_time.tv_sec));

	branching_events = g_ptr_array_sized_new(BUFSIZ);
	params->hm = hm;
	params->all_tpls = &all_tpls;

	//test_kmer_ext(params);
	//exit(1);
	kmer_threads(params);
	// Start branching after the frequent kmers are consumed already.
	//start_branching(&all_tpls, params);

	show_msg(__func__, "Saving contigs: %.2f sec\n",
			(float) (kmer_finish_time.tv_sec - kmer_start_time.tv_sec));
	contigs = xopen(get_output_file("paired.fa", kmer_out), "w");
	kmer_tpls = hash_to_array(&all_tpls);
	save_tpls(kmer_tpls, contigs, 0, 0, 0);
	fflush(contigs);
	fclose(contigs);

	g_ptr_array_sort(branching_events, (GCompareFunc) cmp_junc_by_id);
	clean_junctions(branching_events);
	store_features(get_output_file("paired.junctions", kmer_out),
			branching_events, kmer_tpls);
	process_graph(kmer_tpls, branching_events, hm);
	//	reset_tpl_ids(all_kmer_tpls);
	//	merge_ol_tpls(all_kmer_tpls, ins_size, sd_ins_size, hm->seqs,
	//			kmer_n_threads);
	//	free(params);
	//	destroy_hm(hm);
	//
	//	contigs = xopen(get_output_file("merged.fa", kmer_out), "w");
	//	save_tpls(all_kmer_tpls, contigs, 0, 0, 100);
	//	fflush(contigs);
	//	fclose(contigs);
	//	contigs = xopen(get_output_file("merged_all.fa", kmer_out), "w");
	//	save_tpls(all_kmer_tpls, contigs, 0, 0, 0);
	//	fflush(contigs);
	//	fclose(contigs);
}

void read_juncs_from_file(char *junc_fn, char *pair_fa, GPtrArray *all_tpls,
		GPtrArray *all_junctions) {
	FILE *junc_fp = xopen(junc_fn, "r");
	bwa_seq_t *seqs = NULL, *ctg = NULL;
	uint32_t n_ctgs = 0, i = 0, id = 0;
	seqs = load_reads(pair_fa, &n_ctgs);
	tpl *t = NULL, *main_tpl = NULL, *branch = NULL;
	tpl_hash tpls;
	char buf[BUFSIZ];
	char *attr[18];
	junction *jun = NULL;
	for (i = 0; i < n_ctgs; i++) {
		ctg = &seqs[i];
		t = new_eg();
		if (atoi(ctg->name) - id == 1) {
			ctg = &seqs[i];
			t->id = atoi(ctg->name);
			t->ctg = new_seq(ctg, ctg->len, 0);
			id = t->id;
		} else {
			t->ctg = blank_seq(0);
			t->id = ++id;
			i--;
		}
		tpls[t->id] = t;
		t->len = t->ctg->len;
		t->alive = 1;
		g_ptr_array_add(all_tpls, t);
	}

	while (fgets(buf, sizeof(buf), junc_fp)) {
		i = 0;
		attr[0] = strtok(buf, "\t");
		while (attr[i] != NULL) { //ensure a pointer was found
			//			printf("fields[%d] = %s\n", i, fields[i]);
			attr[++i] = strtok(NULL, "\t"); //continue to tokenize the string
		}
		if (atoi(attr[3]) < 2)
			continue;
		main_tpl = tpls[atoi(attr[0])];
		branch = tpls[atoi(attr[1])];
		jun = new_junction(main_tpl, branch, 0, atoi(attr[2]), atoi(attr[4]),
				atoi(attr[3]));
		g_ptr_array_add(all_junctions, jun);
	}
	bwa_free_read_seq(n_ctgs, seqs);
}

void process_only(char *junc_fn, char *pair_fa, char *hash_fn) {
	GPtrArray *all_tpls = g_ptr_array_sized_new(32);
	GPtrArray *all_junctions = g_ptr_array_sized_new(32);
	junction *j = NULL;
	read_juncs_from_file(junc_fn, pair_fa, all_tpls, all_junctions);
	uint32_t i = 0;
	mer_hash map;
	hash_map *hm = load_hash_map(hash_fn, 1, map);
	filter_junctions(all_junctions, all_tpls, hm);
	process_graph(all_tpls, all_junctions, hm);
}

int pe_kmer(int argc, char *argv[]) {
	//	process_only("../SRX011545_out/paired.junctions.nolen",
	//			"../SRX011545_out/paired.fa",
	//			"/home/ariyaratnep/shaojiang/peta/rnaseq/hg19/SRX011545/both.fa");
	//	return 0;

	int c = 0;
	clock_gettime(CLOCK_MONOTONIC, &kmer_start_time);
	while ((c = getopt(argc, argv, "k:m:s:o:t:")) >= 0) {
		switch (c) {
		case 'k':
			kmer_len = atoi(optarg);
			break;
		case 'o':
			kmer_out = optarg;
			break;
		case 'm':
			ins_size = atoi(optarg);
			break;
		case 's':
			sd_ins_size = atoi(optarg);
			break;
		case 't':
			kmer_n_threads = atoi(optarg);
			break;
		}
	}
	if (optind + 2 > argc) {
		show_msg(__func__, "Parameters error! \n");
		return 1;
	}
	if (!g_thread_supported())
		g_thread_init( NULL);
	kmer_id_mutex = g_mutex_new();

	ext_by_kmers_core(argv[optind], argv[optind + 1]);

	clock_gettime(CLOCK_MONOTONIC, &kmer_finish_time);
	show_msg(__func__, "Done: %.2f sec\n", (float) (kmer_finish_time.tv_sec
			- kmer_start_time.tv_sec));
	return 0;
}
