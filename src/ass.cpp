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
#include "ass.hpp"
#include "junction.hpp"
#include "graph.hpp"
#include "pool.hpp"
#include "merge.hpp"
#include "hash.hpp"

using namespace std;

int kmer_ctg_id = 1;
int ins_size = 0;
int sd_ins_size = 0;
int kmer_n_threads = 0;
int kmer_len = 0;
int stage = 1;
int fresh_trial = 0;
uint32_t n_used_reads = 0;
char *kmer_out = NULL;
GPtrArray *hang_reads = NULL;

bwa_seq_t *TEST = NULL;

GMutex *kmer_id_mutex;
struct timespec kmer_start_time, kmer_finish_time;

/**
 * Initialize a template
 */
tpl *blank_tpl(bwa_seq_t *start_read, int len, int ori) {
	tpl *t = new_tpl();
	g_mutex_lock(kmer_id_mutex);
	t->id = kmer_ctg_id++;
	g_mutex_unlock(kmer_id_mutex);
	if (start_read->rev_com)
		switch_fr(start_read);
	if (len == 0) {
		t->ctg = blank_seq(start_read->len);
	} else {
		if (len < 0 || len > start_read->len) {
			t->ctg = new_seq(start_read, start_read->len, 0);
		} else {
			t->ctg = ori ? new_seq(start_read, len, 0) : new_seq(start_read,
					len, start_read->len - len);
		}
	}
	if (start_read->rev_com)
		switch_fr(start_read);
	t->start_read = start_read;
	t->len = t->ctg->len;
	return t;
}

GPtrArray *hash_to_array(tpl_hash *all_tpls) {
	tpl *t = NULL;
	GPtrArray *tpls = g_ptr_array_sized_new(all_tpls->size());
	for (tpl_hash::iterator m = all_tpls->begin(); m != all_tpls->end(); ++m) {
		t = (tpl*) m->second;
		if (t->alive) {
			show_debug_msg(__func__, "LAST READ [%d, %d] \n", t->id, t->len);
			p_query("LAST READ", t->last_read);
			//p_tpl_reads(t);
			//p_junctions(t->m_juncs);
			//p_junctions(t->b_juncs);
			g_ptr_array_add(tpls, t);
		} else {
			//destory_tpl_junctions(t);
			destroy_tpl(t, TRIED);
		}
	}
	g_ptr_array_sort(tpls, (GCompareFunc) cmp_tpl_by_id);
	return tpls;
}

/**
 * Reset the rev_com flag on used reads
 */
void reset_rev_com(tpl_hash *all_tpls, bwa_seq_t *r) {
	int rev_com = 0, n_mis = 0, i = 0;
	tpl *t = NULL;
	tpl_hash::iterator it = all_tpls->find(r->contig_id);
	if (it == all_tpls->end()) {
		return;
	}
	t = (tpl*) it->second;
	for (i = 0; i < r->len; i++) {
		if (r->contig_locus + i < 0 || r->contig_locus + i >= t->len)
			continue;
		if (r->seq[i] != t->ctg->seq[r->contig_locus + i])
			n_mis++;
		if (n_mis > N_MISMATCHES + 2)
			rev_com = 1;
	}
	r->rev_com = rev_com ? 1 : 0;
	r->cursor = -1;
	r->pos = IMPOSSIBLE_NEGATIVE;
}

/**
 * Use read-length tail to search,
 * 	find those templates could be connected to current branch
 */
GPtrArray *find_connected_reads(hash_table *ht, tpl_hash *all_tpls,
		tpl *branch, const int ori) {
	bwa_seq_t *tail = NULL, *r = NULL, *copy_r = NULL, *branch_seq = NULL;
	int ol_len = 0, n_mis = 0;
	int i = 0, j = 0, branch_main_same_ori = 0;
	GPtrArray *mains = g_ptr_array_sized_new(0);
	GPtrArray *hits = NULL;
	tpl *m_tpl = NULL;

	branch->ctg->rev_com = 0;
	tail = get_tail(branch, kmer_len, ori);

	// If the tail is like 'AAAAAAATAAAA', ignore
	if (is_biased_q(tail) || is_repetitive_q(tail) || tail->len < kmer_len) {
		bwa_free_read_seq(1, tail);
		return mains;
	}

	//	p_test();
	//	p_query(__func__, tail);
	hits = align_query(ht, tail, USED, N_MISMATCHES);

	//	if (branch->id == 4) {
	//		p_tpl(branch);
	//		show_debug_msg(__func__, "To the %s \n", ori ? "left" : "right");
	//		p_query(__func__, tail);
	//		p_readarray(hits, 1);
	//	}
	// Check whether the branch and main are reverse complement or not
	for (i = 0; i < hits->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(hits, i);
		n_mis = 0;
		if (r->contig_id != branch->id) {
			tpl_hash::iterator it = all_tpls->find(r->contig_id);
			if (it != all_tpls->end()) {
				m_tpl = (tpl*) it->second;
				branch_main_same_ori
						= rev_com_on_tpl(m_tpl, r->contig_locus, r) ? 0 : 1;
				break;
			}
		}
	}
	//show_debug_msg(__func__, "branch_main_same_ori: %d \n",
	//		branch_main_same_ori);
	branch_seq = blank_seq(ht->o->read_len);
	set_rev_com(branch->ctg);
	for (i = 0; i < hits->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(hits, i);
		if (r->contig_id == branch->id) {
			reset_rev_com(all_tpls, r);
			continue;
		}
		n_mis = -1;
		if (ori) {
			ol_len = r->len - r->pos;
			if (branch->len < ol_len) {
				reset_rev_com(all_tpls, r);
				continue;
			}
			copy_partial(branch->ctg, branch_seq, 0, ol_len);
			n_mis = seq_ol(r, branch_seq, ol_len, N_MISMATCHES);
			if (n_mis >= 0)
				r->cursor = r->pos - 1;
		} else {
			ol_len = r->pos + tail->len;
			if (branch->len < ol_len) {
				reset_rev_com(all_tpls, r);
				continue;
			}
			copy_partial(branch->ctg, branch_seq, branch->len - ol_len, ol_len);

			n_mis = seq_ol(branch_seq, r, ol_len, N_MISMATCHES);
			if (n_mis >= 0)
				r->cursor = ol_len;
			//p_query(__func__, branch_seq);
			//p_query(__func__, r);
			//show_debug_msg(__func__, "Overlap length: %d\n", ol_len);
			//show_debug_msg(__func__, "Mismatches: %d\n---\n", n_mis);
		}
		if (n_mis >= 0) {
			g_ptr_array_add(mains, r);
		} else {
			reset_rev_com(all_tpls, r);
		}
	}
	bwa_free_read_seq(1, branch_seq);

	for (i = 0; i < mains->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(mains, i);
		//p_query(__func__, r);
	}
	if (ori)
		g_ptr_array_sort(mains, (GCompareFunc) cmp_reads_by_cursor);
	else
		g_ptr_array_sort(mains, (GCompareFunc) cmp_reads_by_rev_cursor);

	bwa_free_read_seq(1, tail);
	g_ptr_array_free(hits, TRUE);
	return mains;
}

/**
 * Check whether the branch can be merged to the main template
 */
int merge_branch_to_main(hash_table *ht, tpl *branch) {
	junction *exist_junc = NULL, *right = NULL, *left = NULL;
	int i = 0;
	int l_len = 0, r_len = 0, t_len = 0;
	int on_main = 0, merged = 0, similar = 0;
	int exist_ori = 0;
	tpl *main_tpl = NULL;
	bwa_seq_t *main_seq = NULL, *sub = NULL, *r = NULL;
	// If right and left connections are too close, just ignore.
	//p_tpl(branch);
	if (!branch->alive || !branch->b_juncs || branch->b_juncs->len <= 0)
		return 0;
	if (branch->b_juncs->len == 1) {
		exist_junc = (junction*) g_ptr_array_index(branch->b_juncs, 0);
		exist_ori = exist_junc->ori;
		main_tpl = exist_junc->main_tpl;
		on_main = branch_on_main(main_tpl, branch, exist_junc->locus,
				branch->len * (1 - BRANCH_SIMILARITY), exist_ori);
		if (on_main) {
			// Mark it as 'dead', will be destroyed in kmer_ext_thread.
			show_debug_msg(__func__,
					"Template [%d, %d] merged to template [%d, %d] \n",
					branch->id, branch->len, main_tpl->id, main_tpl->len);
			branch->alive = 0;
		}
	}
	if (branch->b_juncs->len == 2) {
		for (i = 0; i < branch->b_juncs->len; i++) {
			exist_junc = (junction*) g_ptr_array_index(branch->b_juncs, i);
			if (exist_junc->ori)
				right = exist_junc;
			else
				left = exist_junc;
		}
		if (!right || !left || right == left || left->locus > right->locus) {
			show_debug_msg(__func__, "WARNING: wrong with junctions. \n");
			p_junction(left);
			p_junction(right);
			branch->alive = 0;
			return 0;
		}
		if (left->main_tpl != right->main_tpl) {
			return 0;
		}
		if ((get_abs(left->locus - right->locus) <= IGNORE_DIFF && branch->len
				<= IGNORE_DIFF + ht->o->k * 3)) {
			branch->alive = 0;
			show_debug_msg(__func__,
					"Ignored the template [%d, %d] because too short\n",
					branch->id, branch->len);
			return 0;
		}
		main_tpl = left->main_tpl;
		main_seq = get_tpl_ctg_wt(main_tpl, &l_len, &r_len, &t_len);
		if (left->locus + l_len < 0) {
			show_debug_msg(__func__, "WARNING: locus not correct\n");
			p_junction(left);
		} else {
			sub = new_seq(main_seq, right->locus - left->locus, left->locus
					+ l_len);
			similar = similar_seqs(sub, branch->ctg, branch->len * (1
					- BRANCH_SIMILARITY), MAX_GAPS, MATCH_SCORE,
					MISMATCH_SCORE, INDEL_SCORE);
			//p_ctg_seq("MAIN", sub);
			//p_ctg_seq("BRANCH", branch->ctg);
			if (similar) {
				merged = 1;
				for (i = 0; i < branch->reads->len; i++) {
					r = (bwa_seq_t*) g_ptr_array_index(branch->reads, i);
					add2tpl(main_tpl, r, r->contig_locus + left->locus);
				}
				while (branch->reads->len > 0)
					g_ptr_array_remove_index_fast(branch->reads, 0);
				branch->alive = 0;
				show_debug_msg(__func__,
						"Template [%d, %d] merged to template [%d, %d] \n",
						branch->id, branch->len, main_tpl->id, main_tpl->len);
			}
			bwa_free_read_seq(1, sub);
		}
		bwa_free_read_seq(1, main_seq);
	}
	return merged;
}

/**
 * For branches on a template, count the reads at the junctions
 */
void set_jun_reads(hash_table *ht, tpl *t) {
	GPtrArray *juncs = t->m_juncs;
	junction *jun = NULL;
	int i = 0, changed = 0;
	tpl *branch = NULL;
	ubyte_t *new_seq = NULL;
	bwa_seq_t *jun_seq = NULL;
	if (!juncs || !t->alive)
		return;
	g_ptr_array_sort(juncs, (GCompareFunc) cmp_junc_by_locus);
	for (i = 0; i < juncs->len; i++) {
		jun = (junction*) g_ptr_array_index(juncs, i);
		if (jun->status != 0 || jun->weight > 0)
			continue;
		//p_junction(jun);
		//p_tpl_reads(jun->main_tpl);
		//p_tpl_reads(jun->branch_tpl);
		//p_tpl(jun->main_tpl);
		//p_tpl(jun->branch_tpl);
		jun->weight = count_jun_reads(ht, jun);
		//		if (jun->weight < MIN_JUNCTION_READS) {
		//			destroy_junction(jun);
		//			i--;
		//		}
		//		p_junction(jun);
	}
}

void p_tpl_junctions(tpl_hash *all_tpls, int id) {
	tpl_hash::iterator it = all_tpls->find(id);
	if (it == all_tpls->end()) {
		show_debug_msg(__func__, "No template id %d \n", id);
		return;
	}
	tpl *t = (tpl*) it->second;
	printf("\n");
	p_tpl(t);
	show_debug_msg(__func__, "Template [%d, %d]: as main branches: \n", t->id,
			t->len);
	p_junctions(t->m_juncs);
	show_debug_msg(__func__, "Template [%d, %d]: as branch branches: \n",
			t->id, t->len);
	p_junctions(t->b_juncs);
	show_debug_msg(__func__, "--- End of printing junctions. --- \n\n");
}

/**
 * Remove a dead template from the global template hash
 */
void rm_global_tpl(tpl_hash *all_tpls, tpl *t, int status) {
	if (!t->alive) {
		show_debug_msg(__func__, "Template [%d, %d] is destroyed.\n", t->id,
				t->len);
		g_mutex_lock(kmer_id_mutex);
		all_tpls->erase(t->id);
		g_mutex_unlock(kmer_id_mutex);
		destory_tpl_junctions(t);
		destroy_tpl(t, status);
	}
}

tpl *add_global_tpl(tpl_hash *all_tpls, bwa_seq_t *branch_read, int len,
		int ori) {
	tpl *branch = blank_tpl(branch_read, len, ori);
	g_mutex_lock(kmer_id_mutex);
	all_tpls->insert(make_pair<int, tpl*> ((int) branch->id, (tpl*) branch));
	g_mutex_unlock(kmer_id_mutex);
	return branch;
}

void p_test_read() {
//	p_query("TEST", TEST);
//	if (TEST->status == FRESH && TEST->pos != IMPOSSIBLE_NEGATIVE) {
//		show_debug_msg(__func__, "ID: %d \n", kmer_ctg_id);
//		exit(1);
//	}
}

/**
 * Mark the reads on a template as TRIED temporarily
 */
void mark_as_hang_tmp(tpl *t) {
	bwa_seq_t *r = 0;
	int i = 0;
	if (t->reads) {
		for (i = 0; i < t->reads->len; i++) {
			r = (bwa_seq_t*) g_ptr_array_index(t->reads, i);
			r->status = HANG;
			r->pos = IMPOSSIBLE_NEGATIVE;
			r->cursor = -1;
			r->contig_id = -1;
			r->contig_locus = -1;
			r->rev_com = 0;
			g_ptr_array_add(hang_reads, r);
		}
	}
	if (t->tried) {
		for (i = 0; i < t->tried->len; i++) {
			r = (bwa_seq_t*) g_ptr_array_index(t->tried, i);
			r->status = HANG;
			r->pos = IMPOSSIBLE_NEGATIVE;
			r->cursor = -1;
			r->contig_id = -1;
			r->contig_locus = -1;
			r->rev_com = 0;
			g_ptr_array_add(hang_reads, r);
		}
	}
}

/**
 * Unfrozen all TRIED reads
 */
void unfrozen_hang_reads() {
	bwa_seq_t *r = NULL;
	while (hang_reads->len > 0) {
		r = (bwa_seq_t*) g_ptr_array_index(hang_reads, 0);
		//if (r->status == HANG)
		reset_to_fresh(r);
		g_ptr_array_remove_index_fast(hang_reads, 0);
	}
}

/**
 * Check whether a branch has paired reads with current component.
 */
int val_branch_by_pairs(hash_table *ht, tpl *main_tpl, tpl *branch_tpl) {
	int has_pairs = 0;
	if (branch_tpl->alive == 0)
		return 0;
	GPtrArray *near_tpls = NULL;
	near_tpls = nearby_tpls(branch_tpl, 0);
	has_pairs = has_nearby_pairs(ht, near_tpls, branch_tpl, MIN_PAIRS);
	g_ptr_array_free(near_tpls, TRUE);
	if (!has_pairs) {
		branch_tpl->alive = 0;
		return 0;
	}
	return 1;
}

/**
 * During extension, if it reaches some read which is used already, try to connect to it.
 */
int connect_by_full_reads(hash_table *ht, tpl_hash *all_tpls, tpl *branch,
		const int ori) {
	GPtrArray *con_reads = NULL;
	int i = 0, j = 0, tmp = 0;
	bwa_seq_t *r = NULL, *tail = NULL;
	tpl *main_tpl = NULL;
	ubyte_t branch_c = 0, main_c = 0;
	// Positions, lengths, etc.
	int con_pos = 0, exist_ori = 0, loop_len = 0;
	int ori_len = 0, truncated_len = 0, ol_len = 0;
	// Indicators, etc
	int valid = 0, weight = 0;
	int connected = 0, is_rev = 0;
	int rev_com = 0, n_mis = 0;
	// If the branch is reverse complement connected, the direction needs to be switch
	int adj_ori = 0;
	int max_trial = 0;
	junction *exist_junc = NULL;

	if (!branch || !branch->alive || !branch->reads || branch->reads->len == 0)
		return 0;

	// If extending to the left, and it's not connected to any template, mark it 'dead'
	if (ori && (!branch->b_juncs || branch->b_juncs->len == 0) && branch->len
			<= ht->o->read_len) {
		mark_tpl_dead(branch);
		return 0;
	}
	printf("\n");
	show_debug_msg(__func__, "Trying to connect template [%d, %d] to %s ...\n",
			branch->id, branch->len, ori ? "left" : "right");
	con_reads = find_connected_reads(ht, all_tpls, branch, ori);
	//show_debug_msg(__func__, "Connecting reads: \n");
	//p_readarray(con_reads, 1);
	max_trial = con_reads->len > 4 ? 4 : con_reads->len;
	for (i = 0; i < max_trial; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(con_reads, i);

		tpl_hash::iterator it = all_tpls->find(r->contig_id);
		if (it == all_tpls->end()) {
			continue;
		}

		adj_ori = ori;
		// The candidate template to connect
		main_tpl = (tpl*) it->second;

		// If the main template is too short, just ignore
		if (main_tpl->len <= ht->o->k || main_tpl == branch
				|| !val_branch_by_pairs(ht, main_tpl, branch) || main_tpl->cov
				* MIN_BRANCH_MAIN_COV > (branch->reads->len * ht->o->read_len
				/ max(branch->len, 1))) {
			continue;
		}

		// If the branch is reversed complement in last round,
		//	but it is not successfully connected, need to reverse it back for this round.
		if (is_rev) {
			switch_tpl_fr(branch);
			is_rev = 0;
		}

		//p_query("CONNECTOR", r);
		//p_tpl(main_tpl);
		//p_tpl(branch);
		/**
		 * When connecting, need to check whether it's reverse-complement or not
		 * So need to get the subseq on branch and main to compare base-by-base
		 */

		rev_com = rev_com_on_tpl(main_tpl, r->contig_locus, r);
		if (rev_com != r->rev_com) {
			if ((branch->b_juncs && branch->b_juncs->len > 0)
					|| (branch->m_juncs && branch->m_juncs->len > 0)) {
				show_debug_msg(
						__func__,
						"[WARNING] Template [%d, %d] is already connected, should not reverse. \n",
						branch->id, branch->len);
				p_tpl_junctions(all_tpls, branch->id);
				continue;
			}
			//p_tpl(branch);
			//p_tpl(main_tpl);
			show_debug_msg(__func__,
					"Branch [%d, %d] is reverse complemented \n", branch->id,
					branch->len);
			// Reverse the branch and the direction
			switch_tpl_fr(branch);
			adj_ori = ori ? 0 : 1;
			is_rev = 1;
			r->pos = r->len - kmer_len - r->pos;
			r->cursor = adj_ori ? r->pos - 1 : r->pos + kmer_len;
		}
		//p_ctg_seq("TO CONNECT", branch->ctg);

		p_query("CONNECTOR", r);
		con_pos = ori ? r->contig_locus + 1 + r->cursor : r->contig_locus;
		show_debug_msg(__func__, "ORI: %d; CONPOS: %d\n", ori, con_pos);

		exist_ori = adj_ori ? 0 : 1;
		//p_ctg_seq("MAIN", main_tpl->ctg);
		//p_ctg_seq("BRANCH", branch->ctg);
		ol_len = adj_ori ? (r->len - r->cursor) : r->cursor;
		find_reads_ahead(branch, ht->o->read_len, ol_len, &weight, ori);
		valid = weight >= MIN_JUNCTION_READS ? 1 : 0;
		//show_debug_msg(__func__, "BRANCH READS: %d \n", branch->reads->len);
		//p_readarray(branch->reads, 1);
		//show_debug_msg(__func__, "Junction reads: %d \n", weight);

		if (!valid) {
			show_debug_msg(__func__, "No enough junction reads. Ignore. \n");
			break;
		} // End of checking junction reads and pairs

		// Make sure the truncated subsequence is exactly the same on branch and main
		truncated_len = 0;
		if (exist_ori) {
			tmp = 0;
			for (j = branch->len - 1; j >= branch->len - r->cursor; j--) {
				branch_c = branch->ctg->seq[j];
				main_c = main_tpl->ctg->seq[r->contig_locus + r->cursor
						- (++tmp)];
				if (branch_c == main_c)
					truncated_len++;
				else
					break;
			}
			con_pos = r->contig_locus + (r->cursor - truncated_len);
		} else {
			// Determine how many bases to truncate
			for (j = 0; j < r->len - r->cursor - 1; j++) {
				branch_c = branch->ctg->seq[j];
				main_c
						= main_tpl->ctg->seq[r->contig_locus + r->cursor + j
								+ 1];
				//show_debug_msg(__func__, "Main_c vs Branch_c: %c vs %c \n", "ACGTN"[main_c], "ACGTN"[branch_c]);
				if (branch_c == main_c)
					truncated_len++;
				else
					break;
			}
			con_pos = r->contig_locus + r->cursor + truncated_len + 1;
		}

		show_debug_msg(__func__, "Truncated length: %d \n", truncated_len);

		// Here we are sure they are going to connect,
		show_debug_msg(__func__,
				"Trying to connect [%d, %d] and [%d, %d] at %d ori %d...\n",
				main_tpl->id, main_tpl->len, branch->id, branch->len, con_pos,
				adj_ori);

		if (con_pos < 0 || con_pos >= main_tpl->len + 1) {
			show_debug_msg(__func__,
					"WARNING: the connecting position %d is not valid \n",
					con_pos);
			continue;
		}

		// If the there is a small loop, determine the loop len
		loop_len = 0;
		if (branch->b_juncs && branch->b_juncs->len == 1) {
			exist_junc = (junction*) g_ptr_array_index(branch->b_juncs, 0);
			if (exist_junc->status == 0 && exist_junc->main_tpl == main_tpl
					&& exist_junc->ori != exist_ori) {
				if (exist_ori == 1 && con_pos < exist_junc->locus) {
					loop_len = exist_junc->locus - con_pos;
					loop_len = (branch->len - r->cursor + loop_len < 0) ? 0
							: loop_len;
				}
				if (exist_ori == 0 && con_pos > exist_junc->locus) {
					loop_len = con_pos - exist_junc->locus;
					loop_len = (branch->len - (r->len - r->cursor - 1)
							- loop_len) < 0 ? 0 : loop_len;
				}
				// In case too long loop length
				loop_len = loop_len > 8 ? loop_len : 0;
				con_pos = loop_len > 0 ? exist_junc->locus : con_pos;
			}
		}
		show_debug_msg(__func__,
				"Loop length: %d; con_pos: %d; exist_ori: %d \n", loop_len,
				con_pos, exist_ori);

		// Trim the branch
		if (branch->len - truncated_len + loop_len > branch->len)
			loop_len = 0;
		if (exist_ori) {
			branch->len = branch->len - truncated_len + loop_len;
		} else {
			show_debug_msg(__func__, "Original length: %d\n", branch->len);
			show_debug_msg(__func__, "Truncated length: %d\n", truncated_len);
			ori_len = branch->len;
			branch->len = branch->len - truncated_len + loop_len;
			memmove(branch->ctg->seq, branch->ctg->seq + ori_len - branch->len,
					sizeof(ubyte_t) * branch->len);
			upd_reads_after_truncate(branch, truncated_len - loop_len);
		}
		branch->ctg->len = branch->len;
		//p_ctg_seq("TRUNCATED", branch->ctg);
		set_rev_com(branch->ctg);

		// Finally! Go to add the junction!
		show_debug_msg(__func__,
				"Connect existing [%d, %d] to [%d, %d] at %d with ori %d. \n",
				branch->id, branch->len, main_tpl->id, main_tpl->len, con_pos,
				exist_ori);
		set_tail(branch, main_tpl, con_pos, ht->o->read_len - 1, exist_ori);
		//p_ctg_seq("Right tail", branch->r_tail);
		//p_ctg_seq("Left  tail", branch->l_tail);
		add_a_junction(main_tpl, branch, 0, con_pos, exist_ori, weight);
		if (is_rev) {
			// Means reverse complement sequence connected
			connected = 2;
		} else {
			// Means forward sequence connected
			connected = 1;
		}
		set_jun_reads(ht, main_tpl);
		break;
	} // End of connecting all probable templates
	// If it is not connected and is reverse complemented, just get it back
	if (is_rev && !connected) {
		switch_tpl_fr(branch);
	}
	// Reset the attributes of the branching reads.
	for (i = 0; i < con_reads->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(con_reads, i);
		r->pos = IMPOSSIBLE_NEGATIVE;
		r->cursor = -1;
		tpl_hash::iterator it = all_tpls->find(r->contig_id);
		if (it == all_tpls->end()) {
			continue;
		}
		main_tpl = (tpl*) it->second;
		rev_com = rev_com_on_tpl(main_tpl, r->contig_locus, r);
		if (rev_com != r->rev_com) {
			r->rev_com = r->rev_com ? 0 : 1;
		}
	}
	g_ptr_array_free(con_reads, TRUE);
	return connected;
}

/**
 * Look harder for mate reads at head/tail;
 * If there is another direction supported by a paired read, switch to that direction
 */
bwa_seq_t *look_harder_at_tail(hash_table *ht, pool *p,
		GPtrArray *near_tpls, tpl *t, int ori) {
	bwa_seq_t *r = NULL, *tail = NULL;
	int i = 0, j = 0, start = 1, end = ht->o->read_len;
	if (!t || t->len < ht->o->read_len)
		return NULL;
	if (ori == 0) {
		start = t->len - ht->o->read_len;
		end = t->len;
	}
	tail = new_seq(t->ctg, kmer_len, start);
	//show_debug_msg(__func__, "Looking harder at [%d, %d] ...\n", start, end);
	//p_ctg_seq("ORI", t->ctg);
	for (i = start; i < end - kmer_len; i++) {
		//p_query(__func__, tail);
		find_match_mates(ht, p, near_tpls, t, tail, N_MISMATCHES, ori);
		if (p->reads->len > 0) {
			//p_pool("HARDER", p, NULL);
			if (ori) {
				t->len -= i;
				t->ctg->len = t->len;
				memmove(t->ctg->seq, t->ctg->seq + i, sizeof(ubyte_t) * t->len);
				set_rev_com(t->ctg);
				for (j = t->reads->len - 1; j >= 0; j--) {
					r = (bwa_seq_t*) g_ptr_array_index(t->reads, j);
					if (r->contig_locus >= t->len) {
						reset_to_dead(r);
						g_ptr_array_remove_index_fast(t->reads, j);
					} else
						break;
				}
			} else {
				t->len = i;
				t->ctg->len = t->len;
				// Remove reads without range
				for (j = t->reads->len - 1; j >= 0; j--) {
					r = (bwa_seq_t*) g_ptr_array_index(t->reads, j);
					if (r->contig_locus + r->len >= i) {
						reset_to_dead(r);
						g_ptr_array_remove_index_fast(t->reads, j);
					} else
						break;
				}
			}
			//p_ctg_seq("TRUNCATED", t->ctg);
			return tail;
		}
		ext_que(tail, t->ctg->seq[i + kmer_len], 0);
	}
	bwa_free_read_seq(1, tail);
	return NULL;
}

/**
 * Extend a template until no next kmer
 */
int kmer_ext_tpl(hash_table *ht, tpl_hash *all_tpls, pool *p, tpl *t,
		bwa_seq_t *query, const int ori) {
	int max_c = -1;
	int to_connect = 0;
	int ext_len = 0, no_read_len = 0;
	GPtrArray *near_tpls = g_ptr_array_sized_new(4);
	bwa_seq_t *tail = new_seq(query, query->len, 0), *adj_tail = NULL;
	bwa_seq_t *last_read = NULL;

	show_debug_msg(__func__,
			"------ Started extending tpl %d to ori %d ... ------\n", t->id,
			ori);
	p_query(__func__, tail);
	near_tpls = nearby_tpls(t, 1);
	p_test_read();
	while (1) {
		// If the query is bad, stop
		if (is_bad_query(tail)) {
			show_debug_msg(__func__, "Repetitive tail, stop.\n");
			p_query(__func__, tail);
			mark_pool_reads_tried(p, t);
			to_connect = 0;
			break;
		}

		// If any long-enough subsequence is not covered by any read, try to connect
		if (no_read_len >= ht->o->read_len * 1.5) {
			show_debug_msg(__func__,
					"[%d, %d] %d bases not covered by a read. \n", t->id,
					t->len, no_read_len);
			show_debug_msg(__func__, "No hits, stop ori %d: [%d, %d] \n", ori,
					t->id, t->len);
			to_connect = 1;
			break;
		}

		if (p->reads->len > 0)
			last_read
					= (bwa_seq_t*) g_ptr_array_index(p->reads, p->reads->len - 1);

		if (t->id == -1) {
		p_query(__func__, tail);
		p_ctg_seq("TEMPLATE", t->ctg);
		p_pool("CURRENT POOL", p, NULL);
		}

		max_c = get_next_char(ht, p, near_tpls, t, ori);

		if (t->id == -1)
		show_debug_msg(__func__,
				"Ori: %d, Template [%d, %d], Next char: %c \n", ori, t->id,
				t->len, "ZACGTN"[max_c + 1]);

		// If cannot extend, try to add mates into the pool
		if (max_c == -1) {
			//p_ctg_seq("TEMPLATE", t->ctg);
			//show_debug_msg(__func__, "Looking for mates on [%d, %d] ...\n",
			//		t->id, t->len);
			//p_tpl_reads(t);
			find_hashed_mates(ht, p, near_tpls, t, tail, N_MISMATCHES,
					ori);
			max_c = get_next_char(ht, p, near_tpls, t, ori);
			if (max_c == -1) {
				//adj_tail = look_harder_at_tail(ht, p, near_tpls, t, ori);
				if (adj_tail) {
					bwa_free_read_seq(1, tail);
					tail = adj_tail;
					max_c = get_next_char(ht, p, near_tpls, t, ori);
				}
				if (max_c == -1) {
					show_debug_msg(__func__, "No hits, stop ori %d: [%d, %d] \n",
							ori, t->id, t->len);
					p_query("LAST READ", last_read);
					if (last_read)
						t->last_read = last_read;
					to_connect = 1;
					break;
				}
			} else {
				p_query("TAIL", tail);
				show_debug_msg(__func__, "Added mates: ori %d \n", ori);
				p_pool("MATE_POOL", p, NULL);
			}
		}

		ext_con(t->ctg, max_c, ori);
		t->len = t->ctg->len;
		ext_que(tail, max_c, ori);
		ext_len++;
		// If the extended length is save long enough, refresh the frozen reads.
		// @Desperate: risky, may produce infinite extension
		//if (ext_len % ht->o->read_len == 0) {
		//	unfrozen_tried(t);
		//}
		// If the overlapped region between t and r has too many mismatches, remove from pool
		rm_half_clip_reads(p, t, max_c, N_MISMATCHES, ori);
		if (forward(p, t, ori)) {
			no_read_len = 0;
		} else {
			no_read_len++;
		}
		// Try to align the tail only if:
		//	1. once for every 4bp
		//	2. the reads in pool is less than 4
		//	3. not consuming the pool
		if (t->len % 1 == 0 || p->reads->len <= 4)
			next_pool(ht, p, t, tail, LESS_MISMATCH, ori);

		if (t->len % 100 == 0)
			show_debug_msg(__func__, "Ori %d, tpl %d, length %d \n", ori,
					t->id, t->len);
		if (t->len > 100000) {
			p_tpl(t);
			err_fatal(
					__func__,
					"[ERROR] Too long contig [%d, %d]. Please report to caishaojiang@gmail.com.\n",
					t->id, t->len);
		}
	}
	set_rev_com(t->ctg);
	g_ptr_array_free(near_tpls, TRUE);
	bwa_free_read_seq(1, tail);
	p_test_read();
	return to_connect;
}

/**
 * Check whether need to jump, starting from the 'read'
 * At least MIN_GAPPED_PAIRS pairs should be found
 */
tpl *please_jump(hash_table *ht, tpl_hash *all_tpls, tpl *from, bwa_seq_t *read) {
	pool *p = NULL;
	tpl *to = NULL;
	int i = 0;
	bwa_seq_t *r = NULL, *m = NULL;
	to = add_global_tpl(all_tpls, read, read->len, 0);
	return to;
	p = new_pool();
	init_pool(ht, p, to, kmer_len, N_MISMATCHES, 0);
	init_pool(ht, p, to, kmer_len, N_MISMATCHES, 1);
	g_ptr_array_add(p->reads, read);
	//p_readarray(p->reads, 1);
	for (i = 0; i < p->reads->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(p->reads, i);
		m = get_mate(r, ht->seqs);
		if (m->status == USED && m->contig_id == from->id) {
			//p_query("USED ", m);
			//p_query("FRESH", r);
		} else {
			reset_to_fresh(r);
			g_ptr_array_remove_index_fast(p->reads, i--);
		}
	}
	if (p->reads->len < MIN_GAPPED_PAIRS) {
		to->alive = 0;
		rm_global_tpl(all_tpls, to, FRESH);
		destroy_pool(p);
		return NULL;
	} else {
		destroy_pool(p);
		return to;
	}
}

/**
 * Reaching this step, the initial pool cannot be empty
 */
void do_jumping(hash_table *ht, tpl_hash *all_tpls, tpl *t, bwa_seq_t *r) {
	pool *p = new_pool();
	int to_con_left = 0, to_con_right = 0;
	bwa_seq_t *query = NULL;
	int pre_len = t->len;
	int pre_n_reads = t->reads->len;

	show_debug_msg(__func__, "Jumping to read %s as template %d\n", r->name,
			t->id);

	init_pool(ht, p, t, kmer_len, N_MISMATCHES, 1);
	correct_init_tpl_base(p, t, 1);
	query = get_tail(t, kmer_len, 1);
	to_con_left = kmer_ext_tpl(ht, all_tpls, p, t, query, 1);
	bwa_free_read_seq(1, query);
	destroy_pool(p);
	set_rev_com(t->ctg);

	// Maybe marked as not alive in last extension
	if (!t->alive)
		return;

	upd_locus_on_tpl(t, pre_len, pre_n_reads);

	// Then extend to the left
	query = get_tail(t, kmer_len, 0);
	p = new_pool();
	init_pool(ht, p, t, kmer_len, N_MISMATCHES, 0);
	//p_query(__func__, query);
	//p_pool("INITIAL_POOL", p, NULL);
	to_con_right = kmer_ext_tpl(ht, all_tpls, p, t, query, 0);
	set_rev_com(t->ctg);
	destroy_pool(p);
	bwa_free_read_seq(1, query);

	show_debug_msg(__func__, "Jumping tpl %d with length: %d \n", t->id, t->len);
	//p_tpl_reads(t);
}

void tpl_jumping(hash_table *ht, tpl_hash *all_tpls, tpl *from) {
	int i = 0, merged = 0;
	bwa_seq_t *r = NULL, *m = NULL;
	tpl *to = NULL;
	pool *p = NULL;
	if (!from || !from->alive)
		return;
	//	g_ptr_array_sort(from->reads, (GCompareFunc) cmp_reads_by_contig_locus);
	for (i = 0; i < from->reads->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(from->reads, i);
		m = get_mate(r, ht->seqs);
		//p_test();
		//		show_debug_msg(__func__, "i = %d\n", i);
		//		p_query(__func__, r);
		//		p_query(__func__, m);
		// If the read is in the middle of the template, must do not jump from its mate
		if (m->status != FRESH || (r->contig_locus > ins_size + 100
				&& r->contig_locus < from->len - ins_size - 100))
			continue;
		show_debug_msg(__func__, "Trying to jump to %s ...\n", m->name);
		to = please_jump(ht, all_tpls, from, m);
		if (!to)
			continue;
		mark_init_reads_used(ht, to, m, N_MISMATCHES);
		do_jumping(ht, all_tpls, to, m);
		unfrozen_tried(to);
		if (!to->alive) {
			mark_as_hang_tmp(to);
			to->alive = 0;
			rm_global_tpl(all_tpls, to, HANG);
			continue;
		}
		merged = merged_jumped(ht, from, to, MORE_MISMATCH);
		if (merged) {
			show_debug_msg(__func__,
					"Jumped to read %s [%d, %d] as [%d, %d]...\n", m->name,
					to->id, to->len, from->id, from->len);
			//			g_ptr_array_sort(from->reads, (GCompareFunc) cmp_reads_by_contig_locus);
			//p_tpl_reads(from);
			i = 0;
			unfrozen_hang_reads();
		}
		// For TRIED reads; to->reads is empty after merging
		mark_as_hang_tmp(to);
		to->alive = 0;
		rm_global_tpl(all_tpls, to, HANG);
		show_debug_msg(__func__, "--- End of jumping %s\n\n", m->name);
	}
	unfrozen_hang_reads();
}

/**
 * Mark a template not alive if:
 * 1. Shorter than 100bp
 * 2. Coverage too low
 * 3. No pairs on the template
 */
int try_destroy_tpl(hash_table *ht, tpl_hash *all_tpls, tpl *t, int read_len) {
	int i = 0;
	int is_valid = 1, with_pairs = 1;
	float branch_cov = 0.0, main_cov = 0.0;
	junction *jun = NULL;
	tpl *main_tpl = NULL;
	GPtrArray *near_tpls = NULL;
	if (!t)
		return 1;
	unfrozen_tried(t);
	if (!t->alive)
		return 1;
	// If it is a short leaf branch, check the coverage
	if (t->len <= 100) {
		is_valid = 0;
		if (t->b_juncs) {
			if (t->b_juncs->len == 2)
				is_valid = 1;
			else if (t->b_juncs->len == 1) {
				jun = (junction*) g_ptr_array_index(t->b_juncs, 0);
				main_tpl = jun->main_tpl;
				if (main_tpl->len > 0 && t->len > 0) {
					branch_cov = ((float) t->reads->len) / ((float) t->len);
					main_cov = ((float) main_tpl->len)
							/ ((float) main_tpl->len);
					if (branch_cov >= main_cov * 1000)
						is_valid = 1;
				}
			}
		}
		if (!is_valid)
			show_debug_msg(__func__,
					"Template [%d, %d] is shorter than 100bp \n", t->id, t->len);
	}
	//p_tpl_reads(t);
	if (is_valid && t->cov <= MIN_BRANCH_COV)
		is_valid = 0;
	with_pairs = has_pairs_on_tpl(ht, t, MIN_PAIRS);
	//show_debug_msg(__func__, "Has pairs on template [%d, %d]: %d\n", t->id, t->len, with_pairs);
	// At any stage, if the template is longer than insert size, require some pairs
	if (is_valid && t->len >= ins_size + 100 && !with_pairs) {
		show_debug_msg(__func__, "No pairs on template [%d, %d]\n", t->id,
				t->len);
		is_valid = 0;
	}
	if (!is_valid) {
		t->alive = 0;
		return 1;
	}
	return 0;
}

/**
 * For those branches at the head/tail, choose the longer branch
 */
int prune_tpl_tails(hash_table *ht, tpl_hash *all_tpls, tpl *t) {
	GPtrArray *branches = t->m_juncs;
	junction *jun = NULL, *jun2 = NULL;
	int i = 0, j = 0, changed = 0, pruned_right = 0;
	int added_len_to_left = 0, freed = 1;
	tpl *branch = NULL;
	ubyte_t *new_seq = NULL;
	float main_cov = 0.0, branch_cov = 0.0;
	if (!branches || branches->len <= 0 || !t->alive)
		return 0;

	show_debug_msg(__func__, "Pruning head/tails of template [%d, %d] ... \n",
			t->id, t->len);
	g_ptr_array_sort(branches, (GCompareFunc) cmp_junc_by_locus);

	//p_tpl(t);
	//p_junctions(branches);
	//p_tpl_reads(t);
	while (branches->len > 0 && freed && !(added_len_to_left && pruned_right)) {
		freed = 0;
		branches = t->m_juncs;
		for (i = 0; i < branches->len; i++) {
			jun = (junction*) g_ptr_array_index(branches, i);
			branch = jun->branch_tpl;
			if (jun->status != 0 || branch == t)
				continue;
			if (!branch->alive) {
				jun->status = 1;
				continue;
			}
			if (jun->ori == 0 && t->r_tail)
				continue;
			if (jun->ori == 1 && t->l_tail)
				continue;
			if (jun->locus < 0 || jun->locus >= t->len) {
				show_debug_msg(__func__,
						"[WARNING] Junction locus not correct. Removed. \n");
				p_junction(jun);
				destroy_junction(jun);
				freed = 1;
				break;
			}
			//p_tpl(branch);
			//p_tpl_reads(branch);
			if (added_len_to_left == 0 && jun->locus < ht->o->read_len
					&& jun->ori == 1) {
				main_cov = calc_tpl_cov(t, 0, jun->locus, ht->o->read_len);
				show_debug_msg(__func__,
						"Branch coverage: %.2f; main coverage: %.2f\n",
						branch->cov, main_cov);
				if (branch->len > jun->locus && main_cov < LOW_PART_COV
						&& branch->cov > main_cov) {
					p_junction(jun);
					mv_reads_bt_tpls(branch, t, jun->locus, 1);
					new_seq = (ubyte_t*) calloc(branch->len + t->len,
							sizeof(ubyte_t));
					memcpy(new_seq, branch->ctg->seq, sizeof(ubyte_t)
							* branch->len);
					memcpy(new_seq + branch->len, t->ctg->seq + jun->locus,
							sizeof(ubyte_t) * (t->ctg->len - jun->locus));
					added_len_to_left = branch->len - jun->locus;
					free(t->ctg->seq);
					t->ctg->seq = new_seq;
					t->len = t->len - jun->locus + branch->len;
					t->ctg->len = t->len;
					t->ctg->full_len = t->len;
					set_rev_com(t->ctg);
					changed = 1;
					mark_tpl_dead(branch);
					show_debug_msg(__func__,
							"Branch [%d, %d] merged to Main [%d, %d]\n",
							branch->id, branch->len, t->id, t->len);
					// Update the locus of all other junctions
					for (j = 0; j < branches->len; j++) {
						jun2 = (junction*) g_ptr_array_index(branches, j);
						if (jun2->status == 0) {
							// Remove all junctions at the trimmed head
							if (jun2->locus < jun->locus || jun2->branch_tpl
									== branch)
								jun2->status = 1;
							else
								jun2->locus += added_len_to_left;
						}
					}
				}
			} else if (!pruned_right && jun->locus > (t->len - ht->o->read_len)
					&& jun->ori == 0) {
				p_junction(jun);
				main_cov = calc_tpl_cov(t, jun->locus, t->len, ht->o->read_len);
				show_debug_msg(__func__,
						"Branch coverage: %.2f; main coverage: %.2f\n",
						branch->cov, main_cov);
				if (branch->len > (t->len - jun->locus) && main_cov
						< LOW_PART_COV && branch->cov > main_cov) {
					mv_reads_bt_tpls(branch, t, t->len - jun->locus, 0);
					new_seq = (ubyte_t*) calloc(branch->len + t->len,
							sizeof(ubyte_t));
					memcpy(new_seq, t->ctg->seq, sizeof(ubyte_t) * jun->locus);
					memcpy(new_seq + jun->locus, branch->ctg->seq, branch->len);
					free(t->ctg->seq);
					t->ctg->seq = new_seq;
					t->len = jun->locus + branch->len;
					t->ctg->len = t->len;
					t->ctg->full_len = t->len;
					set_rev_com(t->ctg);
					changed = 1;
					pruned_right = 1;
					mark_tpl_dead(branch);
					show_debug_msg(__func__,
							"Branch [%d, %d] merged to Main [%d, %d]\n",
							branch->id, branch->len, t->id, t->len);
					// Update the locus of all other junctions
					for (j = 0; j < branches->len; j++) {
						jun2 = (junction*) g_ptr_array_index(branches, j);
						if (jun2->status == 0) {
							// Remove all junctions at the trimmed head
							if (jun2->locus > jun->locus || jun2->branch_tpl
									== branch)
								jun2->status = 1;
						}
					}
				}
			}
			if (branch->cov < t->cov * MIN_BRANCH_MAIN_COV)
				branch->alive = 0;
			merge_branch_to_main(ht, branch);
			val_branch_by_pairs(ht, t, branch);
			try_destroy_tpl(ht, all_tpls, branch, ht->o->read_len);

			//p_tpl(t);
			if (!branch->alive) {
				rm_global_tpl(all_tpls, branch, FRESH);
				freed = 1;
				break;
			}
			//p_tpl(t);
		}
	}
	return changed;
}

/**
 * Perform branching without restriction
 */
void branching(hash_table *ht, tpl_hash *all_tpls, tpl *t, int mismatches,
		int ori) {
	bwa_seq_t *tail = NULL, *branch_read = NULL, *jun_read = NULL, *query =
			NULL;
	int i = 0, j = 0, x = 0, shift = 0, cursor = 0, pos = 0, read_status = DEAD;
	int con_pos = 0, n_junc_reads = 0;
	int exist_ori = ori, dead = 0, to_connect = 0, connected = 0;
	tpl *branch = NULL;
	pool *p = NULL;
	int least_ol_len = kmer_len;
	GPtrArray *b_reads = NULL;
	junction *jun = NULL;

	if (!t || !t->alive || !t->ctg || t->len <= least_ol_len)
		return;
	unfrozen_tried(t);
	printf("\n");
	show_debug_msg(__func__,
			"===== Branching template [%d, %d] to %s ===== \n", t->id, t->len,
			ori ? "left" : "right");
	tail = ori ? new_seq(t->ctg, least_ol_len, t->len - least_ol_len)
			: new_seq(t->ctg, least_ol_len, 0);
	p_test_read();
	for (i = 0; i <= t->len - least_ol_len; i++) {
		if (i > 0) {
			if (ori)
				ext_que(tail, t->ctg->seq[t->len - i - least_ol_len], 1);
			else
				ext_que(tail, t->ctg->seq[i - 1 + least_ol_len], 0);
		}

		shift = ori ? t->len - least_ol_len - i : i;

		//show_debug_msg(__func__, "Template [%d, %d] at %d (i: %d) \n", t->id,
		//					t->len, shift, i);
		//p_query(__func__, tail);

		b_reads = check_branch_tail(ht, t, tail, shift, mismatches, FRESH, ori);
		p_test_read();
		if (b_reads->len > 0) {
			printf("\n ---- \n");
			show_debug_msg(__func__, "Template [%d, %d] at %d \n", t->id,
					t->len, shift);
			p_query(__func__, tail);
			//p_readarray(b_reads, 1);
		}
		if (ori)
			g_ptr_array_sort(b_reads, (GCompareFunc) cmp_reads_by_cursor);
		else
			g_ptr_array_sort(b_reads, (GCompareFunc) cmp_reads_by_rev_cursor);

		for (j = 0; j < b_reads->len; j++) {
			// For later truncate the branch template
			branch_read = (bwa_seq_t*) g_ptr_array_index(b_reads, j);

			cursor = branch_read->cursor;
			pos = branch_read->pos;

			// Add the branching junction first;
			// Later may add connection junction
			con_pos = ori ? shift - (pos - cursor - 1) : cursor - pos + shift;
			// In case the junctions create some small loop
			if (has_nearby_junc(t, con_pos)) {
				break;
			}
			// Create a new template
			branch = add_global_tpl(all_tpls, branch_read, branch_read->len,
					ori);
			branch->len = 1;
			branch->ctg->seq[0] = branch->ctg->seq[cursor];
			branch->ctg->len = branch->len;
			branch->ctg->rseq[0] = 3 - branch->ctg->seq[0];
			set_tail(branch, t, con_pos, ht->o->read_len - 1, exist_ori);
			jun = add_a_junction(t, branch, NULL, con_pos, exist_ori,
					n_junc_reads);

			p = new_pool();
			for (x = 0; x < b_reads->len; x++) {
				jun_read = (bwa_seq_t*) g_ptr_array_index(b_reads, x);
				jun_read->cursor = ori ? jun_read->cursor - 1
						: jun_read->cursor + 1;
				// !Important: pos in the pool means how many mismatches
				jun_read->pos = 0;
				add2pool(p, jun_read);
			}

			jun_read = get_tail(t, ht->o->read_len, ori);
			mark_init_reads_used(ht, branch, jun_read, mismatches);
			bwa_free_read_seq(1, jun_read);
			for (x = 0; x < branch->reads->len; x++) {
				jun_read = (bwa_seq_t*) g_ptr_array_index(branch->reads, x);
				jun_read->contig_locus = 0 - (jun_read->len - branch->len);
			}

			p_query("BRANCH_QUERY", branch_read);
			//show_debug_msg(__func__, "shift: %d; POS: %d; CURSOR: %d\n", shift,
			//		pos, cursor);
			//show_debug_msg(__func__, "Branching at %d \n", con_pos);
			//p_pool(__func__, p, NULL);
			p_ctg_seq("AFTER TRUNCATE", branch->ctg);

			// Initialization
			if (branch_read->rev_com)
				switch_fr(branch_read);
			query = ori ? new_seq(branch_read, kmer_len, cursor) : new_seq(
					branch_read, kmer_len, cursor + 1 - kmer_len);
			query->rev_com = 0;
			if (branch_read->rev_com)
				switch_fr(branch_read);

			// Perform extension
			printf("\n");
			show_debug_msg(
					__func__,
					"===== Branching template [%d, %d] to %s at %d; Locus %d, Started %d ===== \n",
					t->id, t->len, ori ? "left" : "right", shift, con_pos,
					branch->id);

			to_connect = kmer_ext_tpl(ht, all_tpls, p, branch, query, ori);
			destroy_pool(p);
			bwa_free_read_seq(1, query);
			if (ori)
				upd_locus_on_tpl(branch, 0, 0);

			if (to_connect && stage == 1) {
				connected = connect_by_full_reads(ht, all_tpls, branch, ori);
			}
			set_rev_com(branch->ctg);

			//unfrozen_tried(branch);
			refresh_tpl_reads(ht, branch, mismatches);
			dead = 0;
			if ((branch->len <= branch_read->len && !connected && !is_tail_junction(jun))
					|| (!branch->alive) || (t->cov * MIN_BRANCH_MAIN_COV)
					> branch->cov)
				dead = 1;
			// If because of on pairs, mark as FRESH, not DEAD.
			if (!dead) {
				if (!val_branch_by_pairs(ht, t, branch)) {
					dead = 1;
				}
			}
			if (dead) mark_as_hang_tmp(branch);

			if (!dead) {
				unfrozen_tried(branch);
				correct_tpl_base(ht->seqs, branch, ht->o->read_len);
				//p_tpl(branch);
				//p_junctions(branch->b_juncs);
			} else {
				branch->alive = 0;
				rm_global_tpl(all_tpls, branch, read_status);
			}
			p_test_read();
			break;
		}
		g_ptr_array_free(b_reads, TRUE);
	}
	unfrozen_hang_reads();
	p_test_read();
	bwa_free_read_seq(1, tail);
}

void try_connect(hash_table *ht, tpl_hash *all_tpls, int to_con_left,
		int to_con_right, tpl *t) {
	int connected = 0;
	int ori = 0;
	if (!t->alive)
		return;
	//if (t->id == 3)
	//	p_tpl_reads(t);
	if (to_con_left) {
		connected = connect_by_full_reads(ht, all_tpls, t, 1);
		ori = (connected == 2) ? 1 : 0;
	}
	if (to_con_right) {
		connect_by_full_reads(ht, all_tpls, t, ori);
	}
	//if (t->id == 3)
	//	p_tpl_reads(t);
}

/**
 * Validate a template;
 * Branch the template
 */
void finalize_tpl(hash_table *ht, tpl_hash *all_tpls, tpl *t, int to_branching,
		int to_con_left, int to_con_right) {
	int changed = 1, pre_n_m_juncs = 0;
	unfrozen_tried(t);
	if (t->alive) {
		set_rev_com(t->ctg);
		refresh_tpl_reads(ht, t, N_MISMATCHES);
	}
	try_connect(ht, all_tpls, to_con_left, to_con_right, t);
	// Reactive the TRIED reads to FRESH, for other starting reads
	//show_debug_msg(__func__, "Finalizing template [%d, %d] \n", t->id, t->len);
	//p_tpl(t);
	try_destroy_tpl(ht, all_tpls, t, ht->o->read_len);
	//p_tpl(t);
	if (t->alive) {
		//show_debug_msg(__func__, "Reads before refresh: %d \n", t->reads->len);
		merge_branch_to_main(ht, t);
		try_destroy_tpl(ht, all_tpls, t, ht->o->read_len);

		if (t->alive) {
			//p_tpl_reads(t);
			//p_tpl(t);
			show_debug_msg(
					__func__,
					"==== End of tpl %d with length: %d; reads: %d; Alive: %d ==== \n\n",
					t->id, t->len, t->reads->len, t->alive);
			correct_tpl_base(ht->seqs, t, ht->o->read_len);
			while (changed) {
				if (to_branching) {
					branching(ht, all_tpls, t, LESS_MISMATCH, 0);
					if (t->m_juncs)
						pre_n_m_juncs = t->m_juncs->len;
					branching(ht, all_tpls, t, LESS_MISMATCH, 1);
					// In case some junctions are removed because of no pairs in the first round
					if (t->m_juncs && t->m_juncs->len > pre_n_m_juncs)
						branching(ht, all_tpls, t, LESS_MISMATCH, 0);
				}
				changed = prune_tpl_tails(ht, all_tpls, t);
				set_jun_reads(ht, t);
				strip_branches(ht, all_tpls, t);
			}
			tpl_jumping(ht, all_tpls, t);
		}
	}
	rm_global_tpl(all_tpls, t, FRESH);
}

/**
 * For branches without enough junction reads,
 * 1. try to extend to another direction
 * 2. otherwise, destroy it and free the reads to FRESH
 */
void strip_branches(hash_table *ht, tpl_hash *all_tpls, tpl *t) {
	GPtrArray *branches = t->m_juncs, *near_tpls = NULL;
	junction *jun = NULL;
	int i = 0, changed = 0, opp_ori = 0, ori_len = 0, to_connect = 0,
			pre_n_reads = 0;
	tpl *branch = NULL;
	bwa_seq_t *query = NULL, *l_tail = NULL, *r_tail = NULL;
	pool *p = NULL;
	if (!branches || !t->alive)
		return;
	near_tpls = nearby_tpls(t, 1);
	g_ptr_array_sort(branches, (GCompareFunc) cmp_junc_by_branch_id);
	show_debug_msg(__func__, "There are %d branches on [%d, %d]\n",
			branches->len, t->id, t->len);
	//p_junctions(branches);
	for (i = 0; i < branches->len; i++) {
		jun = (junction*) g_ptr_array_index(branches, i);
		if (jun->status != 0)
			continue;
		branch = jun->branch_tpl;
		if (!branch->b_juncs || branch->b_juncs->len != 1)
			continue;
		if (!branch->alive) {
			jun->status = 1;
			continue;
		}
		ori_len = branch->len;
		pre_n_reads = branch->reads->len;
		//p_junction(jun);
		if (jun->weight <= MIN_JUNCTION_READS || branch->cov
				< jun->main_tpl->cov * BRANCHING_THRE || !has_nearby_pairs(ht,
				near_tpls, branch, MIN_PAIRS)) {
			//p_tpl(branch);
			//p_tpl(jun->main_tpl);
			opp_ori = jun->ori ? 0 : 1;
			if (branch->len < kmer_len)
				continue;
			show_debug_msg(__func__, "Striping template [%d, %d] to %s\n",
					branch->id, branch->len, opp_ori ? "left" : "right");
			//p_junction(jun);
			query = get_pure_tail(branch, kmer_len, opp_ori);
			if (!query)
				continue;
			clear_tpl_tails(branch);
			p = new_pool();
			init_pool(ht, p, branch, kmer_len, N_MISMATCHES, opp_ori);
			p_query(__func__, query);
			p_pool("INITIAL_POOL", p, NULL);

			to_connect = kmer_ext_tpl(ht, all_tpls, p, branch, query, opp_ori);
			if (to_connect && branch->len > ori_len) {
				connect_by_full_reads(ht, all_tpls, branch, opp_ori);
			}
			set_rev_com(branch->ctg);
			destroy_pool(p);
			bwa_free_read_seq(1, query);
			//p_tpl_junctions(all_tpls, 93);
			//p_tpl_junctions(all_tpls, 97);
			show_debug_msg(__func__, "This junction is striped\n");
			destroy_junction(jun);
			if (opp_ori) {
				upd_locus_on_tpl(branch, ori_len, pre_n_reads);
			}
		}
		//		if (branch->alive)
		//			finalize_tpl(ht, all_tpls, branch, 1, 0, 0);
	}
	g_ptr_array_free(near_tpls, TRUE);
}

/**
 * Extend a kmer
 */
void *kmer_ext_thread(gpointer data, gpointer thread_params) {
	tpl *t = NULL;
	int pre_len = 0, to_con_left = 0, to_con_right = 0, ori = 1, flag = 0;
	int invalid = 0;
	int pre_n_reads = 0, iter = 0;
	uint64_t read_id = 0;
	kmer_counter *counter = NULL;
	bwa_seq_t *read = NULL, *query = NULL, *seqs = NULL;
	pool *p = NULL;

	kmer_t_meta *params = (kmer_t_meta*) thread_params;
	tpl_hash *all_tpls = params->all_tpls;
	hash_table *ht = params->ht;

	seqs = ht->seqs;
	counter = (kmer_counter*) data;
	read_id = counter->kmer;
	read = &seqs[read_id];
	//read = &seqs[4374716];

	//show_debug_msg(__func__, "============= %s: %" ID64 " ============ \n",
	//		read->name, counter->count);

	if (is_biased_q(read) || has_n(read, 1) || counter->count < 1
			|| is_repetitive_q(read) || is_biased_q(read) || read->status
			!= FRESH) {
		return NULL;
	}

//	if (fresh_trial == 0)
//		read = &seqs[242415];
//	if (fresh_trial == 1)
//		read = &seqs[68550];10897

	printf("\n");
	show_debug_msg(__func__,
			"============= FRESH %s: %" ID64 " ============ \n", read->name,
			counter->count);
	fresh_trial++;
	p_query(__func__, read);
	t = add_global_tpl(all_tpls, read, read->len, 0);
	t->kmer_freq = counter->count;

	if (counter->count > 1)
		mark_init_reads_used(ht, t, read, N_MISMATCHES);
	else
		add2tpl(t, read, 0);
	// Right->left->right->left...until not extendable
	// If it is connected to somewhere, simply stop
	while (iter++ < 2 && t->len > pre_len && (!t->b_juncs || t->b_juncs->len
			== 0)) {
		ori = 1;
		// Extend to the left first
		p = new_pool();
		init_pool(ht, p, t, kmer_len, N_MISMATCHES, ori);
		correct_init_tpl_base(p, t, ori);
		//p_pool("INITIAL_POOL", p, NULL);
		if (iter == 1 && p->reads->len == 0) {
			t->alive = 0;
			flag = 1;
			destroy_pool(p);
			break;
		}

		query = get_tail(t, kmer_len, ori);
		pre_len = t->len;
		pre_n_reads = t->reads->len;
		//p_query(__func__, query);
		to_con_left = kmer_ext_tpl(ht, all_tpls, p, t, query, ori);
		bwa_free_read_seq(1, query);
		destroy_pool(p);
		set_rev_com(t->ctg);
		show_debug_msg(__func__, "tpl %d with length: %d \n", t->id, t->len);

		// Maybe marked as not alive in last extension
		if (!t->alive)
			break;
		upd_locus_on_tpl(t, pre_len, pre_n_reads);

		// Then extend to the left
		ori = 0;
		query = get_tail(t, kmer_len, ori);
		p = new_pool();
		init_pool(ht, p, t, kmer_len, N_MISMATCHES, ori);
		p_query(__func__, query);
		//p_pool("INITIAL_POOL", p, NULL);
		to_con_right = kmer_ext_tpl(ht, all_tpls, p, t, query, ori);
		set_rev_com(t->ctg);
		destroy_pool(p);
		bwa_free_read_seq(1, query);

		// Maybe marked as not alive in last extension
		if (!t->alive)
			break;
	}
	finalize_tpl(ht, all_tpls, t, 1, to_con_left, to_con_right);
	// If the read cannot be even extend one base, reset the read to fresh
	if (flag) {
		reset_to_fresh(read);
	}
	//	if (t->alive)
	//p_tpl_reads(t);
	return NULL;
}

void kmer_threads(kmer_t_meta *params) {
	GThreadPool *thread_pool = NULL;
	uint64_t i = 0;
	hash_table *ht = params->ht;
	bwa_seq_t *r = NULL, *seqs = ht->seqs;
	kmer_counter *counter = NULL;
	GPtrArray *starting_reads = g_ptr_array_sized_new(ht->n_seqs), *low_reads =
			NULL;

	show_msg(__func__, "Getting read frequencies ... \n");
	for (i = 0; i < ht->n_seqs; i++) {
		r = &seqs[i];
		//show_debug_msg(__func__, "Query %s: %d\n", r->name, ht->n_kmers[i]);
		if (r->status == FRESH) {
			if (ht->n_kmers[i] > 1) {
				counter = (kmer_counter*) malloc(sizeof(kmer_counter));
				counter->kmer = i;
				counter->count = ht->n_kmers[i];
				g_ptr_array_add(starting_reads, counter);
			}
		}
	}

	TEST = &seqs[396776];

	// shrink_ht(ht);

	show_msg(__func__, "Sorting %d initial reads ... \n", starting_reads->len);
	g_ptr_array_sort(starting_reads, (GCompareFunc) cmp_kmers_by_count);
	show_msg(__func__, "Extending by reads ...\n");
	show_msg(
			__func__,
			"----------- Stage 1: templates longer than insert size with pairs ----------\n");
	params->to_try_connect = 1;
	thread_pool = g_thread_pool_new((GFunc) kmer_ext_thread, (gpointer) params,
			1, TRUE, NULL);
	for (i = 0; i < starting_reads->len; i++) {
		if (i % 100000 == 0) {
			show_msg(__func__, "%d templates are obtained. \n",
					params->all_tpls->size());
			show_msg(__func__, "Extending %" ID64 "-th read ... \n", i);
		}
		counter = (kmer_counter*) g_ptr_array_index(starting_reads, i);
		//g_thread_pool_push(thread_pool, (gpointer) counter, NULL);
		kmer_ext_thread(counter, params);
		free(counter);
//		if (fresh_trial >= 1)
//		if (kmer_ctg_id >= 123062)
//		break;
	}
	g_ptr_array_free(starting_reads, TRUE);
	show_msg(__func__, "%d templates are obtained. \n",
			params->all_tpls->size());

		show_debug_msg(__func__, "Remaining reads: \n");
		for (i = 0; i < ht->n_seqs; i++) {
			r = &ht->seqs[i];
			if (r->status != USED && r->status != FRESH && r->status != DEAD)
				p_query(__func__, r);
		}

	/**
	 show_msg(__func__,
	 "----------- Stage 2: branching the %d templates -----------\n",
	 params->all_tpls->size());
	 stage = 2;
	 //branching_long(ht, params->all_tpls);
	 show_msg(__func__, "%d templates are obtained. \n",
	 params->all_tpls->size());

	 show_msg(__func__,
	 "----------- Stage 3: remaining fragmented reads -----------\n");
	 stage = 3;
	 show_msg(__func__, "Counting 11-mers of remaining reads ...\n");

	 low_reads = g_ptr_array_sized_new(ht->n_seqs / 10);
	 // Reset not USED/DEAD reads to FRESH
	 for (i = 0; i < ht->n_seqs; i++) {
	 r = &seqs[i];
	 //show_debug_msg(__func__, "Query %s: %d\n", r->name, ht->n_kmers[i]);
	 if (r->status != USED && r->status != DEAD) {
	 reset_to_fresh(r);
	 counter = (kmer_counter*) malloc(sizeof(kmer_counter));
	 counter->kmer = i;
	 counter->count = 0;
	 g_ptr_array_add(low_reads, counter);
	 }
	 }

	 sort_by_kmers(ht, low_reads);
	 //show_msg(__func__, "Shrinking the hash table ... \n");
	 //shrink_ht(ht);
	 show_msg(__func__, "Extending the remaining %d reads ...\n", low_reads->len);
	 params->to_try_connect = 1;
	 for (i = 0; i < low_reads->len / 10; i++) {
	 counter = (kmer_counter*) g_ptr_array_index(low_reads, i);
	 // If the read does not even share any 11-mer with others, ignore
	 if (counter->count <= (ht->o->read_len - ht->o->k) * 2) {
	 continue;
	 }
	 if (i % 100000 == 0)
	 show_msg(__func__, "Extending %" ID64 "-th low read ... \n", i);
	 kmer_ext_thread(counter, params);
	 //g_thread_pool_push(thread_pool, (gpointer) counter, NULL);
	 }
	 for (i = 0; i < low_reads->len; i++) {
	 counter = (kmer_counter*) g_ptr_array_index(low_reads, i);
	 free(counter);
	 }
	 g_ptr_array_free(low_reads, TRUE);
	 **/

	g_thread_pool_free(thread_pool, 0, 1);
}

void test_kmer_ext(kmer_t_meta *params) {
}

int merge_paired_tpls(hash_table *ht, tpl_hash *all_tpls) {
	tpl *t = NULL, *mt = NULL;
	int i = 0, id = 0, m_id = 0;
	int ol = 0, n_mis = 0;
	int rev_com = 0, paired = 0, merged = 0;
	bwa_seq_t *m = NULL;
	tpl_hash::iterator it;

	for (tpl_hash::iterator im = all_tpls->begin(); im != all_tpls->end(); ++im) {
		id = im->first;
		t = (tpl*) im->second;

		show_debug_msg(__func__, "Trying to merge template %d ...\n", t->id);
		//p_tpl(t);
		// If merged before, the alive value is 0
		if (!t->alive || is_high_cov(t))
			continue;
		for (i = 0; i < t->tried->len; i++) {
			m = (bwa_seq_t*) g_ptr_array_index(t->tried, i);
			//p_query(__func__, m);
			// If its mate is used by another template
			if (m->status == USED && m->contig_id != id) {
				m_id = m->contig_id;
				it = all_tpls->find(m_id);
				if (it == all_tpls->end()) {
					g_ptr_array_remove_index_fast(t->tried, i--);
					continue;
				}
				mt = (tpl*) it->second;
				if (!mt->alive) {
					g_ptr_array_remove_index_fast(t->tried, i--);
					continue;
				}
				// If they have junctions, just ignore
				if (tpls_have_junction(t, mt)) {
					g_ptr_array_remove_index_fast(t->tried, i--);
					continue;
				}
				// At least 2 pairs spanning them
				paired = paired_by_reads(ht->seqs, t, mt, MIN_PAIRS);
				//show_debug_msg(__func__, "Paired: %d \n", paired);
				if (!paired) {
					g_ptr_array_remove_index_fast(t->tried, i--);
					continue;
				}
				// At least 11 base overlap
				ol = find_fr_ol_within_k(mt->ctg, t->ctg, MORE_MISMATCH,
						ht->o->k, ht->o->read_len, 0, &rev_com, &n_mis);
				//p_tpl(mt);
				show_debug_msg(__func__, "OVERLAP: %d\n", ol);
				// At most 1 mismatch in 11bp
				if (ol >= ht->o->k && ol >= n_mis * ht->o->k) {
					if (merge_tpls(t, mt, ol, rev_com)) {
						// Update the t->tried
						mv_unpaired_to_tried(ht->seqs, t, kmer_ctg_id);
						tpl_jumping(ht, all_tpls, t);
						merged = 1;
						i = 0;
					}
				} else {
					ol = find_fr_ol_within_k(t->ctg, mt->ctg, MORE_MISMATCH,
							ht->o->k, ht->o->read_len, 0, &rev_com, &n_mis);
					show_debug_msg(__func__, "OVERLAP: %d\n", ol);
					if (ol >= ht->o->k && ol >= n_mis * ht->o->k) {
						if (merge_tpls(mt, t, ol, rev_com)) {
							mv_unpaired_to_tried(ht->seqs, mt, kmer_ctg_id);
							merged = 1;
							tpl_jumping(ht, all_tpls, mt);
							break;
						}
					}
					g_ptr_array_remove_index_fast(t->tried, i--);
					rm_from_tried(t, mt->id);
				}
				//g_ptr_array_sort(t->reads, (GCompareFunc) cmp_reads_by_contig_locus);
				//p_readarray(t->reads, 1);
				// End of overlap checking and merging
			} // End of this read
		} // End of checking all reads
	} // End of checking all templates
	return merged;
}

/**
 * Merge all templates by pairs and overlapping
 */
void iter_merge(hash_table *ht, tpl_hash *all_tpls, kmer_hash *tpl_kmer_hash) {
	int merge_iter = 0, merged = 1;
	tpl *t = NULL;

	// Add all mates not on current template to t->tried.
	for (tpl_hash::iterator im = all_tpls->begin(); im != all_tpls->end(); ++im) {
		t = (tpl*) im->second;
		unfrozen_tried(t);
		mv_unpaired_to_tried(ht->seqs, t, kmer_ctg_id);
		g_ptr_array_sort(t->reads, (GCompareFunc) cmp_reads_by_name);
	}

	// Multiple iterations, until not merged anymore
	while (merge_iter++ < 8 && merged) {
		//merged = merge_tpls_by_hash(ht, tpl_kmer_hash, all_tpls);
		merged = merge_paired_tpls(ht, all_tpls);
	}
	if (merge_iter == 8) {
		show_msg(__func__,
				"[WARNING] 8 iterations for merging. May be some error. \n");
	}

	// Clear the reads on t->tried, to save space
	for (tpl_hash::iterator it = all_tpls->begin(); it != all_tpls->end(); ++it) {
		t = (tpl*) it->second;
		if (t->tried) {
			while (t->tried->len > 0)
				g_ptr_array_remove_index_fast(t->tried, 0);
		}
	}
}

/**
 * For those junction at the right tail/head of the main template,
 * 	simply remove the junction and merge the templates.
 */
void merge_together_tpls(tpl_hash *all_tpls) {
	tpl *t = NULL;
	junction *jun = NULL;
	int i = 0;
	GPtrArray *tpls = g_ptr_array_sized_new(all_tpls->size());
	for (tpl_hash::iterator m = all_tpls->begin(); m != all_tpls->end(); ++m) {
		t = (tpl*) m->second;
		if (!t->alive || !t->m_juncs || t->m_juncs->len <= 0)
			continue;
		for (i = 0; i < t->m_juncs->len; i++) {
			jun = (junction*) g_ptr_array_index(t->m_juncs, i);
			if (jun->status != 0)
				continue;
			if (jun->locus <= 2 && jun->ori == 1) {
				destroy_junction(jun);
				merge_tpls(jun->branch_tpl, t, jun->locus, 0);
			} else {
				if (jun->locus >= t->len - 2 && jun->ori == 0) {
					destroy_junction(jun);
					merge_tpls(t, jun->branch_tpl, t->len - jun->locus, 0);
				}
			}
		}
	}
}

/**
 * Remove those junctions, whose templates don't have pairs spanning them
 */
void rm_no_pair_junctions(hash_table *ht, tpl_hash *all_tpls) {
	tpl *t = NULL;
	junction *jun = NULL;
	int i = 0, n_pairs = 0;
	for (tpl_hash::iterator m = all_tpls->begin(); m != all_tpls->end(); ++m) {
		t = (tpl*) m->second;
		if (t->m_juncs) {
			for (i = 0; i < t->m_juncs->len; i++) {
				jun = (junction*) g_ptr_array_index(t->m_juncs, i);
				n_pairs = count_pairs(ht->seqs, jun);
				if (n_pairs <= MIN_PAIRS) {
					show_debug_msg(__func__,
							"Too few pairs spanning the templates: %d\n",
							n_pairs);
					destroy_junction(jun);
					i--;
				}
			}
		}
	}
}

void test_smith_waterman(char *fn) {
	bwa_seq_t *r1 = NULL, *r2 = NULL;
	int similar = 0;
	uint64_t n_seqs = 0;
	bwa_seq_t *seqs = load_reads(fn, &n_seqs);
	r1 = &seqs[0];
	r2 = &seqs[1];
	p_query(__func__, r1);
	p_query(__func__, r2);
	similar = similar_seqs(r1, r2, MORE_MISMATCH, MAX_GAPS, MATCH_SCORE,
			MISMATCH_SCORE, INDEL_SCORE);
	exit(1);
}

/*
 * After merging, may go through some low coverage junctions
 */
void ext_after_merging(hash_table *ht, tpl_hash *all_tpls) {
	tpl *t = NULL;
	pool *p = new_pool();
	bwa_seq_t *query = NULL;
	int to_connect = 0;
	for (tpl_hash::iterator m = all_tpls->begin(); m != all_tpls->end(); ++m) {
		t = (tpl*) m->second;
		if (!t->alive || t->len <= kmer_len)
			continue;
		show_debug_msg(__func__, "Template [%d, %d]\n", t->id, t->len);
		while (p->reads->len > 0)
			g_ptr_array_remove_index_fast(p->reads, 0);
		if (!t->r_tail) {
			query = new_seq(t->ctg, kmer_len, t->len - kmer_len);
			to_connect = kmer_ext_tpl(ht, all_tpls, p, t, query, 0);
		}
		if (!t->l_tail) {
			query = new_seq(t->ctg, kmer_len, 0);
			to_connect = kmer_ext_tpl(ht, all_tpls, p, t, query, 1);
		}
		show_debug_msg(__func__, "Template [%d, %d]---\n\n", t->id, t->len);
	}
	destroy_pool(p);
}

void ext_by_kmers_core(char *lib_file, const char *solid_file) {
	FILE *contigs = NULL;
	kmer_t_meta *params = (kmer_t_meta*) calloc(1, sizeof(kmer_t_meta));
	tpl_hash all_tpls;
	GPtrArray *read_tpls = NULL, *branching_events = NULL;
	hash_table *ht = NULL;
	char *fn = NULL, *lib_fn_copy = NULL;
	hang_reads = g_ptr_array_sized_new(1024);
	kmer_hash tpl_kmer_hash;

	show_msg(__func__, "Library: %s \n", lib_file);
	lib_fn_copy = str_dup(lib_file);
	ht = load_k_hash(lib_file);

	clock_gettime(CLOCK_MONOTONIC, &kmer_finish_time);
	show_msg(__func__, "Done loading read hash table: %.2f sec\n",
			(float) (kmer_finish_time.tv_sec - kmer_start_time.tv_sec));

	branching_events = g_ptr_array_sized_new(BUFSIZ);
	params->ht = ht;
	params->all_tpls = &all_tpls;

	mark_low_qua_reads(ht->seqs, ht->n_seqs);

	//test_kmer_ext(params);
	//exit(1);
	kmer_threads(params);
	//re_valid_tpls(ht, &all_tpls);
	// Start branching after the frequent kmers are consumed already.
	//start_branching(&all_tpls, params);

	//show_msg(__func__, "Building kmer hash for %d templates ...\n",
	//			all_tpls.size());
	//build_tpl_hash(tpl_kmer_hash, &all_tpls, ht->o->k, ht->o->read_len);

	show_msg(__func__, "Merging %d templates by pairs and overlapping ...\n",
			all_tpls.size());
	show_debug_msg(__func__,
			"Merging %d templates by pairs and overlapping ...\n",
			all_tpls.size());
	merge_together_tpls(&all_tpls);
	iter_merge(ht, &all_tpls, &tpl_kmer_hash);
	ext_after_merging(ht, &all_tpls);

	//show_msg(__func__, "Removing junctions without pairs ...\n");
	//rm_no_pair_junctions(ht, &all_tpls);

	show_msg(__func__, "Saving contigs: %.2f sec\n",
			(float) (kmer_finish_time.tv_sec - kmer_start_time.tv_sec));
	fn = get_output_file("paired.fa", kmer_out);
	contigs = xopen(fn, "w");
	read_tpls = hash_to_array(&all_tpls);
	all_tpls.clear();

	save_tpls(read_tpls, contigs, 0, 0, 0);
	fflush(contigs);
	fclose(contigs);
	free(fn);

	get_junction_arr(read_tpls, branching_events);
	g_ptr_array_sort(branching_events, (GCompareFunc) cmp_junc_by_id);
	fn = get_output_file("paired.junctions", kmer_out);
	store_features(fn, branching_events, read_tpls);
	free(fn);

	show_msg(__func__, "Reloading the hash table ... \n");
	reload_table(ht, lib_fn_copy);
	free(lib_fn_copy);
	process_graph(read_tpls, branching_events, ht, kmer_out);
	destroy_ht(ht);
}

void read_juncs_from_file(char *junc_fn, char *pair_fa, GPtrArray *all_tpls,
		GPtrArray *all_junctions) {
	FILE *junc_fp = xopen(junc_fn, "r");
	bwa_seq_t *seqs = NULL, *ctg = NULL;
	uint64_t n_ctgs = 0, i = 0, id = 0;
	seqs = load_reads(pair_fa, &n_ctgs);
	tpl *t = NULL, *main_tpl = NULL, *branch = NULL;
	tpl_hash tpls;
	char buf[BUFSIZ];
	char *attr[18], *idstr[18];
	junction *jun = NULL;
	for (i = 0; i < n_ctgs; i++) {
		t = new_tpl();
		ctg = &seqs[i];
		t->id = atoi(ctg->name);
		t->ctg = new_seq(ctg, ctg->len, 0);
		id = t->id;

		tpls[t->id] = t;
		//show_debug_msg(__func__, "template %d \n", t->id);
		t->len = t->ctg->len;
		t->alive = 1;
		g_ptr_array_add(all_tpls, t);
	}

	int line = 0;
	while (fgets(buf, sizeof(buf), junc_fp)) {
		line++;
		if (line == 1)
			continue;
		i = 0;
		attr[0] = strtok(buf, "\t");
		while (attr[i] != NULL) { //ensure a pointer was found
			//			printf("fields[%d] = %s\n", i, fields[i]);
			attr[++i] = strtok(NULL, "\t"); //continue to tokenize the string
		}
		idstr[0] = strtok(attr[0], ", ");
		idstr[1] = strtok(NULL, ", ");
		i = 0;
		while (idstr[0][++i] != '\0') {
			idstr[0][i - 1] = idstr[0][i];
		}
		idstr[0][i - 1] = '\0';

		idstr[2] = strtok(attr[1], ", ");
		idstr[3] = strtok(NULL, ", ");
		i = 0;
		while (idstr[2][++i] != '\0') {
			idstr[2][i - 1] = idstr[2][i];
		}
		idstr[2][i - 1] = '\0';
		//printf("%d\n", atoi(idstr[0]));
		//printf("%d\n", atoi(idstr[2]));

		main_tpl = (tpl*) tpls[atoi(idstr[0])];
		branch = (tpl*) tpls[atoi(idstr[2])];
		jun
				= new_junction(main_tpl, branch, 0, atoi(attr[2]),
						atoi(attr[4]), 0);
		//atoi(attr[3]));
		//p_tpl(main_tpl);
		//p_tpl(branch);
		if (!main_tpl->m_juncs)
			main_tpl->m_juncs = g_ptr_array_sized_new(4);
		if (!branch->b_juncs)
			branch->b_juncs = g_ptr_array_sized_new(4);
		g_ptr_array_add(main_tpl->m_juncs, jun);
		g_ptr_array_add(branch->b_juncs, jun);
		g_ptr_array_add(all_junctions, jun);
	}
	bwa_free_read_seq(n_ctgs, seqs);
}

GPtrArray *read_blat_hits(char *blat_psl) {
	char buf[1000];
	int line_no = 0, i = 0;
	char *attr[32], *intstr[32];
	FILE *psl_fp = xopen(blat_psl, "r");
	GPtrArray *hits = g_ptr_array_sized_new(32);
	blat_hit *h = NULL;
	while (fgets(buf, sizeof(buf), psl_fp)) {
		line_no++;
		if (line_no <= 5)
			continue;
		i = 0;
		attr[0] = strtok(buf, "\t");
		while (attr[i] != NULL) { //ensure a pointer was found
			attr[++i] = strtok(NULL, "\t"); //continue to tokenize the string
		}
		if (i < 20)
			break;
		h = (blat_hit*) malloc(sizeof(blat_hit));
		h->n_match = atoi(attr[0]);
		h->n_mismatch = atoi(attr[1]);
		h->n_rep = atoi(attr[2]);
		h->n_n = atoi(attr[3]);
		h->n_query_gap = atoi(attr[4]);
		h->n_query_gap_base = atoi(attr[5]);
		h->n_ref_gap = atoi(attr[6]);
		h->n_ref_gap_base = atoi(attr[7]);
		h->strand = attr[8][0];
		h->query = strdup(attr[9]);
		h->q_len = atoi(attr[10]);
		h->q_start = atoi(attr[11]);
		h->q_end = atoi(attr[12]);
		h->ref = strdup(attr[13]);
		h->r_len = atoi(attr[14]);
		h->r_start = atoi(attr[15]);
		h->r_end = atoi(attr[16]);
		h->n_block = atoi(attr[17]);
		h->block_size = (int*) calloc(h->n_block, sizeof(int));
		h->block_size[0] = atoi(strtok(attr[18], ","));
		i = 0;
		while (i < h->n_block - 1) {
			h->block_size[++i] = atoi(strtok(NULL, ",")); //continue to tokenize the string
		}
		h->query_block_start = (int*) calloc(h->n_block, sizeof(int));
		h->query_block_start[0] = atoi(strtok(attr[19], ","));
		i = 0;
		while (i < h->n_block - 1) {
			h->query_block_start[++i] = atoi(strtok(NULL, ",")); //continue to tokenize the string
		}
		h->ref_block_start = (int*) calloc(h->n_block, sizeof(int));
		h->ref_block_start[0] = atoi(strtok(attr[20], ","));
		i = 0;
		while (i < h->n_block - 1) {
			h->ref_block_start[++i] = atoi(strtok(NULL, ",")); //continue to tokenize the string
		}
		g_ptr_array_add(hits, h);
	}
	fclose(psl_fp);
	return hits;
}

void p_blat_hit(blat_hit *h) {
	int i = 0;
	printf("Hit\t");
	printf("%d\t", h->n_match);
	printf("%d\t", h->n_mismatch);
	printf("%d\t", h->n_rep);
	printf("%d\t", h->n_n);
	printf("%d\t", h->n_query_gap);
	printf("%d\t", h->n_query_gap_base);
	printf("%d\t", h->n_ref_gap);
	printf("%d\t", h->n_ref_gap_base);
	printf("%c\t", h->strand);
	printf("%s\t", h->query);
	printf("%d\t", h->q_len);
	printf("%d\t", h->q_start);
	printf("%d\t", h->q_end);
	printf("%s\t", h->ref);
	printf("%d\t", h->r_len);
	printf("%d\t", h->r_start);
	printf("%d\t", h->r_end);
	printf("%d\t", h->n_block);
	for (i = 0; i < h->n_block; i++) {
		printf("%d,", h->block_size[i]);
	}
	printf("\t");
	for (i = 0; i < h->n_block; i++) {
		printf("%d,", h->query_block_start[i]);
	}
	printf("\t");
	for (i = 0; i < h->n_block; i++) {
		printf("%d,", h->ref_block_start[i]);
	}
	printf("\t\n");
}

void validate_junctions(char *junc_fn, char *pair_fa, char *pair_psl,
		char *hash_fn) {
	int i = 0, j = 0, x = 0, has_hit = 0, has_valid_hit = 0;
	int left_len = 0, right_len = 0, half_len = 0, len = 0;
	tpl *t = NULL;
	junction *jun = NULL;
	GPtrArray *all_tpls = g_ptr_array_sized_new(32);
	GPtrArray *all_junctions = g_ptr_array_sized_new(32);
	GPtrArray *junc_seqs = g_ptr_array_sized_new(32);
	GPtrArray *paired_hits = NULL;
	blat_hit *h = NULL;
	bwa_seq_t *branch_part = NULL, *main_part = NULL, *junc_seq = NULL;

	paired_hits = read_blat_hits(pair_psl);

	read_juncs_from_file(junc_fn, pair_fa, all_tpls, all_junctions);
	hash_table *ht = load_k_hash(hash_fn);
	half_len = ht->o->read_len - JUNCTION_BOUNDARY_BASE;
	for (i = 0; i < all_tpls->len; i++) {
		t = (tpl*) g_ptr_array_index(all_tpls, i);
		t->reads = g_ptr_array_sized_new(64);
		refresh_tpl_reads(ht, t, N_MISMATCHES);
		g_ptr_array_sort(t->reads, (GCompareFunc) cmp_reads_by_name);
	}
	for (i = 0; i < all_tpls->len; i++) {
		t = (tpl*) g_ptr_array_index(all_tpls, i);
		printf("Inspecting template [%d, %d] ...\n", t->id, t->len);
		set_jun_reads(ht, t);
		if (t->m_juncs) {
			printf("\nMain template hits: \n");
			for (j = 0; j < paired_hits->len; j++) {
				h = (blat_hit*) g_ptr_array_index(paired_hits, j);
				if (atoi(h->query) == t->id) {
					p_blat_hit(h);
				}
			}
			printf("\nBranch template hits: \n");
			for (x = 0; x < t->m_juncs->len; x++) {
				jun = (junction*) g_ptr_array_index(t->m_juncs, x);
				p_junction(jun);
				has_hit = 0;
				has_valid_hit = 0;
				if (jun->ori) {
					len = min3(jun->branch_tpl->len, half_len, jun->locus);
					branch_part = new_seq(jun->branch_tpl->ctg, len,
							jun->branch_tpl->len - len);
					main_part = new_seq(t->ctg, len, jun->locus - len);
				} else {
					len = min3(jun->branch_tpl->len, half_len, t->len
							- jun->locus);
					branch_part = new_seq(jun->branch_tpl->ctg, len, 0);
					main_part = new_seq(t->ctg, len, jun->locus);
				}
				printf("Junction_pairs: %d \n", count_pairs(ht->seqs, jun));
				p_query("Branch_junc", branch_part);
				p_query("Main_junc", main_part);
				printf("Mismatches: %d\n", seq_ol(branch_part, main_part,
						main_part->len, main_part->len));
				//bwa_free_read_seq(1, branch_part);
				bwa_free_read_seq(1, main_part);

				if (jun->ori) {
					len = min3(jun->branch_tpl->len, half_len, t->len
							- jun->locus);
					main_part = new_seq(t->ctg, len, jun->locus);
					junc_seq = blank_seq(main_part->len + branch_part->len);
					memcpy(junc_seq->seq, branch_part->seq, branch_part->len);
					memcpy(junc_seq->seq + branch_part->len, main_part->seq,
							main_part->len);
				} else {
					len = min3(jun->branch_tpl->len, half_len, jun->locus);
					main_part = new_seq(t->ctg, len, jun->locus - len);
					junc_seq = blank_seq(main_part->len + branch_part->len);
					memcpy(junc_seq->seq, main_part->seq, main_part->len);
					memcpy(junc_seq->seq + main_part->len, branch_part->seq,
							branch_part->len);
				}
				junc_seq->len = branch_part->len + main_part->len;
				set_rev_com(junc_seq);
				junc_seq->name = (char*) malloc(sizeof(char) * 1024);
				sprintf(junc_seq->name, ">%d_%d_%d_%d\n", t->id,
						jun->branch_tpl->id, jun->locus, jun->ori);
				g_ptr_array_add(junc_seqs, junc_seq);

				for (j = 0; j < paired_hits->len; j++) {
					h = (blat_hit*) g_ptr_array_index(paired_hits, j);
					if (atoi(h->query) == jun->branch_tpl->id) {
						has_hit = 1;
						if (h->n_query_gap == 0 && h->n_ref_gap == 0
								&& h->n_mismatch < 5) {
							if (jun->ori == 0) {
								if (h->q_start <= 2) {
									has_valid_hit = 1;
								}
							} else {
								if (h->q_end >= h->q_len - 2) {
									has_valid_hit = 1;
								}
							}
						}
						p_blat_hit(h);
					}
				}
				if (!has_hit) {
					printf("VALID: No hit for branch template [%d, %d] \n",
							jun->branch_tpl->id, jun->branch_tpl->len);
				} else {
					if (has_valid_hit) {
						printf("VALID: The branch [%d, %d] is valid. \n",
								jun->branch_tpl->id, jun->branch_tpl->len);
					} else {
						printf("VALID: Wrong branch: [%d, %d] \n",
								jun->branch_tpl->id, jun->branch_tpl->len);
					}
				}
			}
			printf("\n+++\n");
		}
	}
	FILE *fp = xopen("../SRR097897_branch/junc.fa", "w");
	for (i = 0; i < junc_seqs->len; i++) {
		junc_seq = (bwa_seq_t*) g_ptr_array_index(junc_seqs, i);
		save_con(junc_seq->name, junc_seq, fp);
	}
	fclose(fp);
}

void blat_ref(char *joint_fa, char *joint_psl) {
	GPtrArray *paired_hits = NULL;
	int i = 0, j = 0, x = 0, has_hit = 0, has_valid_hit = 0;
	uint64_t n_joint = 0;
	int locus = 0, ori = 0;
	char *attr[18];
	bwa_seq_t *joints = load_reads(joint_fa, &n_joint), *s = NULL;
	blat_hit *h = NULL;
	char *good = NULL, *junc_name = NULL, *name_copy = NULL;
	paired_hits = read_blat_hits(joint_psl);
	for (i = 0; i < n_joint; i++) {
		s = &joints[i];
		p_query(__func__, s);
		has_valid_hit = 0;
		has_hit = 0;
		for (j = 0; j < paired_hits->len; j++) {
			h = (blat_hit*) g_ptr_array_index(paired_hits, j);
			if (strcmp(s->name, h->query) == 0) {
				has_hit = 1;
				x = 0;
				name_copy = strdup(s->name);
				attr[0] = strtok(name_copy, "_");
				while (attr[x] != NULL) {
					attr[++x] = strtok(NULL, "_");
				}
				locus = atoi(attr[2]);
				ori = atoi(attr[3]);
				printf("Locus: %d; ori: %d\n", locus, ori);
				if (h->q_start < locus && h->q_end > locus) {
					if (h->n_mismatch <= 4 && h->n_ref_gap == 0
							&& h->n_query_gap == 0) {
						has_valid_hit = 1;
						good = strdup(h->ref);
						junc_name = strdup(h->query);
					}
				}
				p_blat_hit(h);
			}
		}
		if (!has_hit) {
			printf("REPORT No hit for %s. \n", s->name);
		} else {
			if (has_valid_hit) {
				printf("REPORT valid hit for %s on %s \n", junc_name, good);
			} else {
				printf("REPORT invalid hit for %s \n", s->name);
			}
		}
	}

}

void process_only(char *junc_fn, char *pair_fa, char *hash_fn) {
	GPtrArray *all_tpls = g_ptr_array_sized_new(32);
	GPtrArray *all_junctions = g_ptr_array_sized_new(32);
	read_juncs_from_file(junc_fn, pair_fa, all_tpls, all_junctions);
	hash_table *ht = load_k_hash(hash_fn);
	filter_junctions(all_junctions, all_tpls, ht);
	process_graph(all_tpls, all_junctions, ht, kmer_out);
}

int pe_kmer(int argc, char *argv[]) {
	//	validate_junctions("../SRR097897_branch/paired.junctions",
	//			"../SRR097897_branch/paired.fa", "../SRR097897_branch/paired.fa.psl",
	//			"/home/carl/Projects/peta/rnaseq/Spombe/SRR097897/SRR097897.part.fa");
	//	blat_ref("../SRR097897_part/paired.joint.fa", "../SRR097897_part/paired.joint.fa.psl");
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
		g_thread_init(NULL);
	kmer_id_mutex = g_mutex_new();

	ext_by_kmers_core(argv[optind], argv[optind + 1]);

	clock_gettime(CLOCK_MONOTONIC, &kmer_finish_time);
	show_msg(__func__, "Done: %.2f sec\n", (float) (kmer_finish_time.tv_sec
			- kmer_start_time.tv_sec));
	return 0;
}
