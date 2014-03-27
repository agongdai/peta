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

int TESTING = 0;
int DETAIL_ID = -1;

int test_suffix = 0;
int kmer_ctg_id = 1;
int ins_size = 0;
int sd_ins_size = 0;
int kmer_n_threads = 0;
int kmer_len = 0;
int stage = 1;
int fresh_trial = 0;
uint32_t n_used_reads = 0;
char *kmer_out = NULL;
GPtrArray *hang_reads = NULL, *tpls_await_branching = NULL;

bwa_seq_t *TEST = NULL;
bwa_seq_t *CLUSTER_START_READ = NULL;

GMutex *kmer_id_mutex;
struct timespec kmer_start_time, kmer_finish_time;

/**
 * Initialize a template
 */
tpl *blank_tpl(bwa_seq_t *start_read, int len, int ori, char *step) {
	if (kmer_ctg_id > 68200) exit(1);
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
	t->start_read = start_read;//CLUSTER_START_READ;
	t->len = t->ctg->len;
	t->ctg->rev_com = 0;
	t->tinfo = (testing_info*) malloc(sizeof(testing_info));
	t->tinfo->step = (char*) malloc(BUFSIZE * sizeof(char));
	t->tinfo->ref = NULL;
	t->tinfo->starting_read = start_read;
	sprintf(t->tinfo->step, "%s:", step);
	return t;
}

GPtrArray *hash_to_array(tpl_hash *all_tpls) {
	tpl *t = NULL;
	GPtrArray *tpls = g_ptr_array_sized_new(all_tpls->size());
	show_msg(__func__, "Putting hashed templates to array ... \n");
	for (tpl_hash::iterator m = all_tpls->begin(); m != all_tpls->end(); ++m) {
		t = (tpl*) m->second;
		p_tpl(t);
		if (t->alive) {
			//show_debug_msg(__func__, "LAST READ [%d, %d] \n", t->id, t->len);
			//p_query("LAST READ", t->last_read);
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

	if (branch->id == -2) {
		p_tpl(branch);
		show_debug_msg(__func__, "To the %s \n", ori ? "left" : "right");
		p_query(__func__, tail);
		p_readarray(hits, 1);
	}
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

			n_mis = seq_ol(branch_seq, r, ol_len, MORE_MISMATCH);
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
	if (branch->len >= LONG_TPL_LEN || !branch->alive || !branch->b_juncs
			|| branch->b_juncs->len <= 0)
		return 0;
	if (branch->b_juncs->len == 1) {
		exist_junc = (junction*) g_ptr_array_index(branch->b_juncs, 0);
		exist_ori = exist_junc->ori;
		main_tpl = exist_junc->main_tpl;
		on_main = branch_on_main(main_tpl, branch, exist_junc->locus,
				branch->len * (1 - BRANCH_SIMILARITY) + 2, exist_ori);
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
				<= IGNORE_DIFF * 1.5)) {
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
			sub = new_seq(main_seq, right->locus - left->locus,
					left->locus + l_len);
			similar = similar_seqs(sub, branch->ctg,
					branch->len * (1 - BRANCH_SIMILARITY) + 2, MAX_GAPS,
					MATCH_SCORE, MISMATCH_SCORE, INDEL_SCORE);
			//p_ctg_seq("MAIN", sub);
			//p_ctg_seq("BRANCH", branch->ctg);
			//show_debug_msg(__func__, "Similar: %d\n", similar);
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
	if (t->not_covered)
		status = DEAD;
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

tpl *add_global_tpl(tpl_hash *all_tpls, bwa_seq_t *branch_read, char *step,
		int len, int ori) {
	tpl *branch = blank_tpl(branch_read, len, ori, step);
	g_mutex_lock(kmer_id_mutex);
	all_tpls->insert(make_pair<int, tpl*> ((int) branch->id, (tpl*) branch));
	g_mutex_unlock(kmer_id_mutex);
	return branch;
}

void p_test_read() {
	if (TESTING)
	p_query("TEST", TEST);
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
			reset_to_hang(r);
			g_ptr_array_add(hang_reads, r);
		}
	}
	if (t->tried) {
		for (i = 0; i < t->tried->len; i++) {
			r = (bwa_seq_t*) g_ptr_array_index(t->tried, i);
			reset_to_hang(r);
			g_ptr_array_add(hang_reads, r);
		}
	}
}

void copy_to_hang(GPtrArray *reads) {
	int i = 0;
	bwa_seq_t *r = NULL;
	for (i = 0; i < reads->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(reads, i);
		g_ptr_array_add(hang_reads, r);
		reset_to_hang(r);
	}
}

void set_hang_status() {
	int i = 0;
	bwa_seq_t *r = NULL;
	for (i = 0; i < hang_reads->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(hang_reads, i);
		reset_to_hang(r);
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
	g_ptr_array_add(near_tpls, main_tpl);
	int i = 0;
	for (i = 0; i < near_tpls->len; i++) {
		tpl *t = (tpl*) g_ptr_array_index(near_tpls, i);
		//p_tpl(t);
		show_debug_msg(__func__, "Template [%d, %d] is included\n", t->id,
				t->len);
	}
	has_pairs = has_nearby_pairs(ht, near_tpls, branch_tpl, MIN_PAIRS);
	g_ptr_array_free(near_tpls, TRUE);
	if (!has_pairs) {
		//branch_tpl->alive = 0;
		show_debug_msg(__func__, "No pairs found on template [%d, %d] \n",
				branch_tpl->id, branch_tpl->len);
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
	int i = 0, j = 0, tmp = 0, too_far = 0;
	bwa_seq_t *r = NULL, *tail = NULL;
	tpl *main_tpl = NULL, *branch_tpl = NULL;
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
	GPtrArray *near_tpls = NULL;

	if (!branch || !branch->alive || !branch->reads || branch->reads->len == 0)
		return 0;

	// If extending to the left, and it's not connected to any template, mark it 'dead'
	if (ori && (!branch->b_juncs || branch->b_juncs->len == 0) && branch->len
			<= ht->o->read_len) {
		show_debug_msg(__func__, "[WARNING] Bad branch: [%d, %d] \n",
				branch->id, branch->len);
		mark_tpl_dead(branch);
		return 0;
	}

	printf("\n");
	show_debug_msg(__func__, "Trying to connect template [%d, %d] to %s ...\n",
			branch->id, branch->len, ori ? "left" : "right");
	//p_ctg_seq(__func__, branch->ctg);
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

		//if (branch->id == 75) {
		//	p_tpl(branch);
		//	p_tpl_reads(branch);
		//	printf("*****\n");
		//	p_tpl(main_tpl);
		//	p_tpl_reads(main_tpl);
		//	printf("*****\n");
		//}

		too_far = 1;
		near_tpls = nearby_tpls(main_tpl, 1);
		for (j = 0; j < near_tpls->len; j++) {
			branch_tpl = (tpl*) g_ptr_array_index(near_tpls, j);
			if (branch == branch_tpl) {
				too_far = 0;
				break;
			}
		}
		g_ptr_array_free(near_tpls, TRUE);
		if (too_far)
			continue;

		// If the main template is too short, just ignore
		if (main_tpl->len <= ht->o->k || main_tpl == branch || main_tpl->cov
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
			//p_junction(exist_junc);
			if (exist_junc->status == 0 && exist_junc->main_tpl == main_tpl
					&& exist_junc->ori != exist_ori) {
				if (exist_ori == 1 && con_pos < exist_junc->locus) {
					loop_len = exist_junc->locus - con_pos;
					//show_debug_msg(__func__, "Loop len: %d; branch->len: %d; r->cursor: %d \n", loop_len, branch->len, r->cursor);
					//loop_len = (branch->len - r->cursor + loop_len < 0) ? 0
					//		: loop_len;
				}
				if (exist_ori == 0 && con_pos > exist_junc->locus) {
					loop_len = con_pos - exist_junc->locus;
					//loop_len = (branch->len - (r->len - r->cursor - 1)
					//		- loop_len) < 0 ? 0 : loop_len;
				}
				// In case too long loop length
				loop_len = loop_len <= 30 ? loop_len : 0;
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

		if (exist_ori == 0) {
			if (main_tpl->r_tail && con_pos > main_tpl->len - 10) {
				branch->alive = 0;
				break;
			}
		} else {
			if (main_tpl->l_tail && con_pos < 10) {
				branch->alive = 0;
				break;
			}
		}

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
	if (branch->id == -1)
		p_tpl(branch);
	return connected;
}

/**
 * Extend a template until no next kmer
 */
int kmer_ext_tpl(hash_table *ht, tpl_hash *all_tpls, pool *p, tpl *from,
		tpl *t, bwa_seq_t *query, const int ori) {
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
	if (from)
		g_ptr_array_add(near_tpls, from);
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
		if (no_read_len >= ht->o->read_len * 1.5 && !has_n_in_pool(p)) {
			show_debug_msg(__func__,
					"[%d, %d] %d bases not covered by a read. \n", t->id,
					t->len, no_read_len);
			to_connect = 1;
			t->not_covered = 1;
			break;
		}

		if (p->reads->len > 0)
			last_read = (bwa_seq_t*) g_ptr_array_index(p->reads,
					p->reads->len - 1);

		if (t->id == DETAIL_ID || DETAIL_ID == 0) {
			p_test_read();
			p_query("QUERY", tail);
			p_ctg_seq("TEMPLATE", t->ctg);
			p_pool("CURRENT POOL", p, NULL);
		}

		max_c = get_next_char(ht, p, t, ori);

		if (t->id == DETAIL_ID || DETAIL_ID == 0)
			show_debug_msg(__func__,
					"Ori: %d, Template [%d, %d], Next char: %c \n", ori, t->id,
					t->len, "ZACGTN"[max_c + 1]);

		// If cannot extend, try to add mates into the pool
		if (max_c == -1) {
			//p_ctg_seq("TEMPLATE", t->ctg);
			//show_debug_msg(__func__, "Looking for mates on [%d, %d] ...\n",
			//		t->id, t->len);
			//p_test_read();
			//p_tpl_reads(t);
			find_ol_mates(ht, p, near_tpls, t, tail, N_MISMATCHES, ori);
			//p_pool("MATE_POOL", p, NULL);
			max_c = get_next_char(ht, p, t, ori);
			if (max_c == -1) {
				if (adj_tail) {
					bwa_free_read_seq(1, tail);
					tail = adj_tail;
					max_c = get_next_char(ht, p, t, ori);
				}
				if (max_c == -1) {
					show_debug_msg(__func__,
							"No hits, stop ori %d: [%d, %d] \n", ori, t->id,
							t->len);
					//p_query("LAST READ", last_read);
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
		//if (t->len % 1 == 0 || p->reads->len <= 20)
		next_pool(ht, p, near_tpls, t, tail, N_MISMATCHES, ori);

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
	r = get_mate(read, ht->seqs);
	read->rev_com = r->rev_com;
	to = add_global_tpl(all_tpls, read, "JUMPING", read->len, 0);
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
 * This is a unit to extend a template to left/right.
 * @to_connect: if 1, try to connect to existing templates after extension
 * @from: if doing jumping, we jump from some template 'from', whose reads are used
 * 			for mate pool determination
 * @init_p: initial pool; if NULL, initialise the pool
 * @init_q: initial query; if NULL, get template tail as query
 */
int ext_unit(hash_table *ht, tpl_hash *all_tpls, pool *init_p, tpl *from,
		tpl *t, bwa_seq_t *init_q, int to_connect, int ori) {
	pool *p = init_p;
	int to_con = 0, connected = 0;
	bwa_seq_t *query = NULL;

	if (!t || !t->alive)
		return -1;

	if (!init_p) {
		p = new_pool();
		init_pool(ht, p, t, kmer_len, N_MISMATCHES, ori);
		query = get_tail(t, kmer_len, ori);
	} else {
		query = init_q;
	}
	if (p->reads->len == 0 && t->len <= ht->o->read_len) {
		bwa_free_read_seq(1, query);
		destroy_pool(p);
		return -1;
	}
	// Correct the initial bases if the starting read is not good
	if (t->len == ht->o->read_len)
		correct_init_tpl_base(p, t, ori);
	// If extending to left, reverse the locus first
	if (ori)
		reverse_locus(t);
	to_con = kmer_ext_tpl(ht, all_tpls, p, from, t, query, ori);
	bwa_free_read_seq(1, query);
	destroy_pool(p);

	if (!t->alive)
		return -1;

	set_rev_com(t->ctg);
	if (ori)
		upd_locus_on_tpl(t, 0, 0);
	// Some cases never do connection:
	//	1. The end is repetitive, like 'AAAAAAAAA'
	if (to_connect > 0 && to_con) {
		connected = connect_by_full_reads(ht, all_tpls, t, ori);
	}
	return connected;
}

/**
 * Extend a read to a template;
 * No branching and validation
 */
tpl *ext_a_read(hash_table *ht, tpl_hash *all_tpls, bwa_seq_t *read,
		index64 count) {
	tpl *t = NULL;
	int iter = 0, after_unit = -1, pre_len = 0;
	if (is_biased_q(read) || has_n(read, 1) || is_repetitive_q(read)
			|| read->status != FRESH) {
		return NULL;
	}
	printf("\n");
	show_debug_msg(__func__,
			"============= FRESH %s: %" ID64 " ============ \n", read->name,
			count);

	fresh_trial++;
	p_query(__func__, read);
	t = add_global_tpl(all_tpls, read, "FRESH", read->len, 0);
	t->kmer_freq = count;
	if (count > 1)
		mark_init_reads_used(ht, t, read, N_MISMATCHES);
	else
		add2tpl(t, read, 0);
	// Right->left->right->left...until not extendable
	// If it is connected to somewhere, simply stop
	while (iter++ <= 2 && t->len > pre_len && (!t->b_juncs || t->b_juncs->len
			== 0)) {
		// Extend to the left first
		after_unit = ext_unit(ht, all_tpls, NULL, NULL, t, NULL, 0, 1);
		//p_tpl_reads(t);
		show_debug_msg(__func__, "tpl %d with length: %d \n", t->id, t->len);
		//p_ctg_seq(__func__, t->ctg);
		if (iter > 1 && after_unit == -1)
			break;
		pre_len = t->len;
		after_unit = ext_unit(ht, all_tpls, NULL, NULL, t, NULL, 0, 0);
		show_debug_msg(__func__, "tpl %d with length: %d \n", t->id, t->len);
		if (after_unit == -1)
			break;
	}
	// If the read cannot be even extend one base, reset the read to fresh
	if (t->len == read->len) {
		reset_to_fresh(read);
	}
	return t;
}

/**
 * Reaching this step, the initial pool cannot be empty
 */
void do_jumping(hash_table *ht, tpl_hash *all_tpls, tpl *from, tpl *t,
		bwa_seq_t *r) {
	show_debug_msg(__func__,
			"Jumping from template [%d, %d] to read %s as template %d ...\n",
			from->id, from->len, r->name, t->id);
	ext_unit(ht, all_tpls, NULL, from, t, NULL, 0, 1);
	// Maybe marked as not alive in last extension
	if (!t->alive)
		return;
	ext_unit(ht, all_tpls, NULL, from, t, NULL, 0, 1);
	ext_unit(ht, all_tpls, NULL, from, t, NULL, 0, 0);
	ext_unit(ht, all_tpls, NULL, from, t, NULL, 0, 1);
	ext_unit(ht, all_tpls, NULL, from, t, NULL, 0, 0);
	show_debug_msg(__func__, "Jumping tpl %d with length: %d \n", t->id, t->len);
	//p_tpl_reads(t);
}

void tpl_jumping(hash_table *ht, tpl_hash *all_tpls, tpl *from) {
	return;
	int i = 0, merged = 0;
	bwa_seq_t *r = NULL, *m = NULL;
	tpl *to = NULL;
	pool *p = NULL;
	GPtrArray *hanged = NULL;
	if (!from || !from->alive)
		return;
	//	g_ptr_array_sort(from->reads, (GCompareFunc) cmp_reads_by_contig_locus);
	//p_tpl_reads(from);
	show_debug_msg(__func__, "Started jumping from template [%d, %d] \n",
			from->id, from->len);
	//p_tpl_reads(from);
	p_test_read();
	hanged = g_ptr_array_sized_new(4);
	for (i = 0; i < from->reads->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(from->reads, i);
		m = get_mate(r, ht->seqs);
		//show_debug_msg(__func__, "i = %d\n", i);
		//p_query(__func__, r);
		//p_query(__func__, m);
		// If the read is in the middle of the template, must do not jump from its mate
		if (m->status != FRESH || (r->contig_locus > ins_size + 150
				&& r->contig_locus < from->len - ins_size - 150))
			continue;
		if (binary_exist(hanged, m))
			continue;
		//p_query(__func__, r);
		//p_query(__func__, m);
		//continue;
		p_test_read();
		to = please_jump(ht, all_tpls, from, m);
		//if (to->id > 250)
		//	exit(1);
		if (!to)
			continue;
		mark_init_reads_used(ht, to, m, N_MISMATCHES);
		do_jumping(ht, all_tpls, from, to, m);
		//p_tpl_reads(to);
		//refresh_tpl_reads(ht, to, N_MISMATCHES);
		if (!to->alive) {
			concat_arrays(hanged, to->reads);
			g_ptr_array_sort(hanged, (GCompareFunc) cmp_reads_by_name);
			rm_global_tpl(all_tpls, to, FRESH);
			continue;
		}

		refresh_tpl_reads(ht, to, N_MISMATCHES);
		if (!to->alive) {
			concat_arrays(hanged, to->reads);
			g_ptr_array_sort(hanged, (GCompareFunc) cmp_reads_by_name);
			rm_global_tpl(all_tpls, to, FRESH);
			continue;
		}
		correct_tpl_base(ht->seqs, to, ht->o->read_len, 0, to->len);
		if (to->len >= 100) {
			merged = merged_jumped(ht, from, to, m, MORE_MISMATCH);
			//p_test_read();
			if (merged) {
				show_debug_msg(__func__,
						"Jumped to read %s [%d, %d] as [%d, %d]...\n\n",
						m->name, to->id, to->len, from->id, from->len);
				printf("JUMPED, %d, %d, %s\n", from->id, to->id, m->name);
				unfrozen_tried(to);
				g_ptr_array_sort(from->reads,
						(GCompareFunc) cmp_reads_by_contig_locus);
				i = 0;
				while (hanged->len > 0) {
					r = (bwa_seq_t*) g_ptr_array_index(hanged, 0);
					reset_to_fresh(r);
					g_ptr_array_remove_index_fast(hanged, 0);
				}
			} else {
				show_debug_msg(__func__, "Template [%d, %d] not merged. \n",
						to->id, to->len);
			}
		}
		// If it is alive, long but not merged; put it wait for branching
		if (to->alive && to->len >= LONG_TPL_LEN) {
			printf("JUMPED, %d, %d, %s\n", from->id, to->id, m->name);
			unfrozen_tried(to);
			g_ptr_array_add(tpls_await_branching, to);
		} else {
			to->alive = 0;
			concat_arrays(hanged, to->reads);
			g_ptr_array_sort(hanged, (GCompareFunc) cmp_reads_by_name);
			rm_global_tpl(all_tpls, to, FRESH);
		}
		show_debug_msg(__func__, "--- End of jumping %s\n\n", m->name);
	}
	g_ptr_array_free(hanged, TRUE);
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
	//if (t->len >= LONG_TPL_LEN)
	//	return 0;
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
	if (t->cov > HIHG_COV_THRE) {
		is_valid = 1;
	}
	//p_tpl_reads(t);
	//if (is_valid && t->cov <= MIN_BRANCH_COV)
	//	is_valid = 0;
	//with_pairs = has_pairs_on_tpl(ht, t, MIN_PAIRS);
	//show_debug_msg(__func__, "Has pairs on template [%d, %d]: %d\n", t->id, t->len, with_pairs);
	// At any stage, if the template is longer than insert size, require some pairs
	if (is_valid && t->len >= ins_size + 100 && !has_pairs_on_tpl(ht, t,
			MIN_PAIRS)) {
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
	g_ptr_array_sort(t->reads, (GCompareFunc) cmp_reads_by_contig_locus);

	//p_tpl(t);
	//p_junctions(branches);
	//p_tpl_reads(t);
	while (branches->len > 0 && freed && !(added_len_to_left && pruned_right)) {
		freed = 0;
		branches = t->m_juncs;
		for (i = 0; i < branches->len; i++) {
			jun = (junction*) g_ptr_array_index(branches, i);
			branch = jun->branch_tpl;
			//printf("\n");
			//show_debug_msg(__func__, "Branch [%d, %d] %s \n", branch->id,
			//		branch->len, branch->alive ? "Alive" : "Dead");
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
			// If only one end
			if (!(branch->b_juncs && branch->b_juncs->len == 2)) {
				//p_tpl(branch);
				//p_tpl_reads(branch);
				if (added_len_to_left == 0 && jun->locus < ht->o->read_len
						&& jun->ori == 1) {
					main_cov = calc_tpl_cov(t, 0, jun->locus, ht->o->read_len);
					show_debug_msg(__func__,
							"Branch coverage: %.2f; main coverage: %.2f\n",
							branch->cov, main_cov);
					//p_tpl_reads(t);
					if (branch->len > jun->locus && main_cov < LOW_PART_COV
							&& branch->cov > main_cov) {
						p_junction(jun);
						mv_reads_bt_tpls(branch, t, jun->locus, 1);
						new_seq = (ubyte_t*) calloc(branch->len + t->len,
								sizeof(ubyte_t));
						memcpy(new_seq, branch->ctg->seq,
								sizeof(ubyte_t) * branch->len);
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
								if (jun2->locus < jun->locus
										|| jun2->branch_tpl == branch)
									jun2->status = 1;
								else
									jun2->locus += added_len_to_left;
							}
						}
					}
				} else if (!pruned_right && jun->locus > (t->len
						- ht->o->read_len) && jun->ori == 0) {
					p_junction(jun);
					main_cov = calc_tpl_cov(t, jun->locus, t->len,
							ht->o->read_len);
					show_debug_msg(__func__,
							"Branch coverage: %.2f; main coverage: %.2f\n",
							branch->cov, main_cov);
					if (branch->len > (t->len - jun->locus) && main_cov
							< LOW_PART_COV && branch->cov > main_cov) {
						mv_reads_bt_tpls(branch, t, t->len - jun->locus, 0);
						new_seq = (ubyte_t*) calloc(branch->len + t->len,
								sizeof(ubyte_t));
						memcpy(new_seq, t->ctg->seq,
								sizeof(ubyte_t) * jun->locus);
						memcpy(new_seq + jun->locus, branch->ctg->seq,
								branch->len);
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
								if (jun2->locus > jun->locus
										|| jun2->branch_tpl == branch)
									jun2->status = 1;
							}
						}
					}
				}
			}
			merge_branch_to_main(ht, branch);
			//show_debug_msg(__func__, "Branch [%d, %d] %s \n", branch->id,
			//		branch->len, branch->alive ? "Alive" : "Dead");
			try_destroy_tpl(ht, all_tpls, branch, ht->o->read_len);
			//show_debug_msg(__func__, "Branch [%d, %d] %s \n", branch->id,
			//		branch->len, branch->alive ? "Alive" : "Dead");
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
	bwa_seq_t *tail = NULL, *branch_read = NULL, *jun_read = NULL,
			*start_read = NULL, *mate = NULL;
	int i = 0, j = 0, x = 0, shift = 0, cursor = 0, pos = 0,
			read_status = HANG;
	int con_pos = 0, n_junc_reads = 0;
	int dead = 0, to_connect = 1, connected = 1, n_mis = 0;
	int b_s = 0, b_e = 0, m_s = 0, m_e = 0, score = 0;
	tpl *branch = NULL, *branch_tpl = NULL;
	int least_ol_len = kmer_len;
	GPtrArray *b_reads = NULL, *copy_b_reads = NULL, *wait_for_val = NULL;
	junction *jun = NULL;

	if (!t || !t->alive || !t->ctg || t->len <= least_ol_len)
		return;
	unfrozen_tried(t);
	wait_for_val = g_ptr_array_sized_new(4);
	printf("\n");
	show_debug_msg(__func__,
			"===== Branching template [%d, %d] to %s ===== \n", t->id, t->len,
			ori ? "left" : "right");
	tail = ori ? new_seq(t->ctg, least_ol_len, t->len - least_ol_len)
			: new_seq(t->ctg, least_ol_len, 0);
	p_test_read();
	for (i = 0; i <= t->len - least_ol_len; i++) {
		if (i > 0) {
			if (ori) ext_que(tail, t->ctg->seq[t->len - i - least_ol_len], 1);
			else ext_que(tail, t->ctg->seq[i + least_ol_len - 1], 0);
		}

		shift = ori ? t->len - least_ol_len - i : i;
		//show_debug_msg(__func__, "Template [%d, %d] at %d (i: %d) \n", t->id, t->len, shift, i);
		//p_query(__func__, tail);
		b_reads = check_branch_tail(ht, t, tail, shift, mismatches, FRESH, ori);
		if (b_reads->len == 0) {
			g_ptr_array_free(b_reads, TRUE);
			continue;
		}
		//p_readarray(b_reads, 1);

		// The original reads would be used for extension, so the pos and cursor
		// 	values would be changed. Copy them to copy_b_reads.
		copy_b_reads = g_ptr_array_sized_new(b_reads->len);
		for (x = 0; x < b_reads->len; x++) {
			branch_read = (bwa_seq_t*) g_ptr_array_index(b_reads, x);
			jun_read = new_seq(branch_read, branch_read->len, 0);
			reset_to_fresh(branch_read);
			g_ptr_array_add(copy_b_reads, jun_read);
		}

		for (j = 0; j < b_reads->len; j++) {
			// For later truncate the branch template
			branch_read = (bwa_seq_t*) g_ptr_array_index(b_reads, j);
			jun_read = (bwa_seq_t*) g_ptr_array_index(copy_b_reads, j);
			cursor = jun_read->cursor;
			pos = jun_read->pos;

			// con_pos: the locus to connect the two templates
			con_pos = ori ? shift - (pos - cursor - 1) : cursor - pos + shift;
			// In case the junctions create some small loop
			if (has_nearby_junc(t, con_pos)) break;

			// Create a branch template with length 1bp
			branch = add_global_tpl(all_tpls, branch_read, "BRANCHING", 1, ori);
			branch->ctg->seq[0] = jun_read->rev_com ? branch_read->rseq[cursor] : branch_read->seq[cursor];
			branch->ctg->rseq[0] = 3 - branch->ctg->seq[0];
			set_tail(branch, t, con_pos, ht->o->read_len - 1, ori);
			jun = add_a_junction(t, branch, NULL, con_pos, ori, n_junc_reads);

			// Mark the reads similar to the starting reads as USED
			start_read = get_tail(t, ht->o->read_len, ori);
			mark_init_reads_used(ht, branch, start_read, mismatches);
			bwa_free_read_seq(1, start_read);

			// Perform extension
			printf("\n");
			show_debug_msg(__func__,
					"===== Branching template [%d, %d] to %s at %d; Locus %d, Started %d ===== \n",
					t->id, t->len, ori ? "left" : "right", shift, con_pos, branch->id);
			p_query("BRANCH_QUERY", branch_read);
			p_ctg_seq("BRANCH", branch->ctg);

			connected = ext_unit(ht, all_tpls, NULL, NULL, branch, NULL, to_connect, ori);
			branch->cov = calc_tpl_cov(branch, 0, branch->len, ht->o->read_len);
			dead = 0;
			if (!branch->alive || (!connected && (branch->len
					<= branch_read->len) && branch->cov < HIHG_COV_THRE))
				dead = 1;
			show_debug_msg(__func__, "Dead: %d\n", dead);

			if (!dead) {
				refresh_tpl_reads(ht, branch, mismatches);
				if (!branch->alive) dead = 1;
			}

			if (!dead) {
				g_ptr_array_add(wait_for_val, branch);
				unfrozen_tried(branch);
				correct_tpl_base(ht->seqs, branch, ht->o->read_len, 0,
						branch->len);
				p_tpl(branch);
			} else {
				reset_to_hang(branch_read);
				mark_as_hang_tmp(branch);
				branch->alive = 0;
				rm_global_tpl(all_tpls, branch, read_status);
			}
			break;
		}
		for (j = 0; j < copy_b_reads->len; j++) {
			jun_read = (bwa_seq_t*) g_ptr_array_index(copy_b_reads, j);
			bwa_free_read_seq(1, jun_read);
		}
		g_ptr_array_free(copy_b_reads, TRUE);
		g_ptr_array_free(b_reads, TRUE);
	}
	unfrozen_hang_reads();
	// For cases where two transcripts are highly similar
	// There may not be enough spanning pairs during branching
	// Two branches may validate each other
	show_debug_msg(__func__, "Validating branches of template [%d, %d] ... \n",
			t->id, t->len);
	//p_tpl(t);
	for (i = wait_for_val->len - 1; i >= 0; i--) {
		branch = (tpl*) g_ptr_array_index(wait_for_val, i);
		dead = 1;
		// For the reads at the junctions, must be at least one pair
		for (j = 0; j < branch->reads->len; j++) {
			branch_read = (bwa_seq_t*) g_ptr_array_index(branch->reads, j);
			// Check the junction reads only
			if (branch_read->contig_locus < 0
					|| branch_read->contig_locus > branch->len - branch_read->len) {
				mate = get_mate(branch_read, ht->seqs);
				if (mate->status == USED) {
					// If the mate is on the main template or branch template, it is valid
					if (mate->contig_id == t->id || mate->contig_id == branch->id) {
						dead = 0; break;
					} else { // In case it is a complicated component
						for (x = 0; x < wait_for_val->len; x++) {
							branch_tpl = (tpl*) g_ptr_array_index(wait_for_val, x);
							if (mate->contig_id == branch_tpl->id) {
								dead = 0; break;
							}
						}
					}
				}
			}
		}

		if (!dead) {
			if (branch->b_juncs->len == 1 && !val_branch_by_pairs(ht, t, branch)) dead = 1;
		}

		// Either case validates the branch:
		//	1. One pair at the junction
		//	2. MIN_PAIRS pairs spanning
		if (dead) {
			branch->alive = 0;
			rm_global_tpl(all_tpls, branch, FRESH);
		} else {
			if (branch->b_juncs->len == 1 && branch->cov < HIHG_COV_THRE) {
				tpl_jumping(ht, all_tpls, branch);
			}
		}
	}
	g_ptr_array_free(wait_for_val, TRUE);
	p_test_read();
	bwa_free_read_seq(1, tail);
}

void try_connect(hash_table *ht, tpl_hash *all_tpls, int to_con_left,
		int to_con_right, tpl *t) {
	int connected = 0;
	int ori = 0;
	if (!t->alive)
		return;
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

void iter_branching(hash_table *ht, tpl_hash *all_tpls, tpl *t,
		int to_branching) {
	int changed = 1, pre_n_m_juncs = 0;
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
		//strip_branches(ht, all_tpls, t);
		//break;
	}
}

/**
 * Validate a template;
 * Branch the template
 */
void finalize_tpl(hash_table *ht, tpl_hash *all_tpls, tpl *t, int to_branching,
		int to_con_left, int to_con_right) {
	int changed = 1, pre_n_m_juncs = 0, ori_len = 0;
	int i = 0;
	bwa_seq_t *r = NULL;
	if (t->alive) {
		set_rev_com(t->ctg);
		refresh_tpl_reads(ht, t, N_MISMATCHES);
	}
	// TRIED reads are not going to be used by JUMPING, but not BRANCHING
	unfrozen_tried(t);
	try_connect(ht, all_tpls, to_con_left, to_con_right, t);
	// Reactive the TRIED reads to FRESH, for other starting reads
	//show_debug_msg(__func__, "Finalizing template [%d, %d] \n", t->id, t->len);
	//p_tpl_reads(t);
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
			correct_tpl_base(ht->seqs, t, ht->o->read_len, 0, t->len);
			iter_branching(ht, all_tpls, t, to_branching);
			//p_tpl(t);
			if (t->cov < HIHG_COV_THRE) {
				ori_len = t->len;
				tpl_jumping(ht, all_tpls, t);
				if (t->len > ori_len && to_branching) {
					iter_branching(ht, all_tpls, t, to_branching);
				}
			} else
				show_debug_msg(__func__,
						"Template [%d, %d] no branching, coverage %.2f \n",
						t->id, t->len, t->cov);
		}
	}
	rm_global_tpl(all_tpls, t, FRESH);
}

/**
 * Extend a kmer
 */
void *kmer_ext_thread(gpointer data, gpointer thread_params) {
	uint64_t read_id = 0;
	kmer_counter *counter = NULL;
	bwa_seq_t *read = NULL, *seqs = NULL;
	tpl *t = NULL;
	int i = 0, to_con_left = 0, to_con_right = 0;

	kmer_t_meta *params = (kmer_t_meta*) thread_params;
	tpl_hash *all_tpls = params->all_tpls;
	hash_table *ht = params->ht;

	seqs = ht->seqs;
	counter = (kmer_counter*) data;
	read_id = counter->kmer;
	read = &seqs[read_id];

	if (counter->count < 1) {
		return NULL;
	}
	if (TESTING) {
		if (fresh_trial == 0)
			read = &seqs[TESTING];
		if (fresh_trial == 1)
			read = &seqs[2904319];
		//	if (fresh_trial == 2)
		//		read = &seqs[7299565];
	}

	t = ext_a_read(ht, all_tpls, read, counter->count);
	if (t)
		finalize_tpl(ht, all_tpls, t, 1, 0, 0);
	//ext_unit(ht, all_tpls, NULL, NULL, t, NULL, 0, 1);
	//ext_unit(ht, all_tpls, NULL, NULL, t, NULL, 0, 0);
	while (tpls_await_branching->len > 0) {
		t = (tpl*) g_ptr_array_index(tpls_await_branching,
				tpls_await_branching->len - 1);
		show_debug_msg(__func__, "Await branching template [%d, %d] ... \n",
				t->id, t->len);
		to_con_right = ext_unit(ht, all_tpls, NULL, NULL, t, NULL, 0, 1);
		to_con_left = ext_unit(ht, all_tpls, NULL, NULL, t, NULL, 0, 0);
		g_ptr_array_remove_index_fast(tpls_await_branching,
				tpls_await_branching->len - 1);
		finalize_tpl(ht, all_tpls, t, 1, 0, 0);
		//ext_unit(ht, all_tpls, NULL, NULL, t, NULL, 0, 1);
		//ext_unit(ht, all_tpls, NULL, NULL, t, NULL, 0, 0);
	}

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

	TEST = &seqs[630743];

	// shrink_ht(ht);

	if (!TESTING) {
		show_msg(__func__, "Sorting %d initial reads ... \n",
				starting_reads->len);
		g_ptr_array_sort(starting_reads, (GCompareFunc) cmp_kmers_by_count);
	}
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
show_msg		(__func__, "Extending %" ID64 "-th read ... \n", i);
	}
	counter = (kmer_counter*) g_ptr_array_index(starting_reads, i);
	//g_thread_pool_push(thread_pool, (gpointer) counter, NULL);
	kmer_ext_thread(counter, params);
	free(counter);
	//if (fresh_trial >= 2)
	//if (kmer_ctg_id >= 1200)
	//if (params->all_tpls->size() > 1000)
	if (TESTING && fresh_trial >= 1)
	break;
}
g_ptr_array_free(starting_reads, TRUE);
show_msg(__func__, "%d templates are obtained. \n",
		params->all_tpls->size());

//	show_debug_msg(__func__, "Remaining reads: \n");
//	for (i = 0; i < ht->n_seqs; i++) {
//		r = &ht->seqs[i];
//		if (r->status != USED && r->status != DEAD)
//			p_query(__func__, r);
//	}

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
						//tpl_jumping(ht, all_tpls, t);
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
							//tpl_jumping(ht, all_tpls, mt);
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
	show_debug_msg(__func__, "Merging together templates ...\n");
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
		show_debug_msg(__func__,
				"Extending Template [%d, %d] after merging ... \n", t->id,
				t->len);
		while (p->reads->len > 0)
			g_ptr_array_remove_index_fast(p->reads, 0);
		if (!t->r_tail) {
			query = new_seq(t->ctg, kmer_len, t->len - kmer_len);
			to_connect = kmer_ext_tpl(ht, all_tpls, p, NULL, t, query, 0);
		}
		if (!t->l_tail) {
			query = new_seq(t->ctg, kmer_len, 0);
			to_connect = kmer_ext_tpl(ht, all_tpls, p, NULL, t, query, 1);
		}
		t->is_root = 0;
		show_debug_msg(__func__, "Template [%d, %d]---\n\n", t->id, t->len);
	}
	destroy_pool(p);
}

void save_read_status(hash_table *ht) {
	bwa_seq_t *r = NULL;
	int i = 0;
	char buf[BUFSIZE];
	char *fn = get_output_file("reads.status", kmer_out);
	FILE *f = xopen(fn, "w");
	for (i = 0; i < ht->n_seqs; i++) {
		r = &ht->seqs[i];
		sprintf(buf, "%s\t%d\t%d\n", r->name, r->status, r->contig_id);
		fputs(buf, f);
	}
	free(fn);
	fclose(f);
}

void ext_by_kmers_core(char *lib_file, const char *solid_file) {
	FILE *contigs = NULL;
	kmer_t_meta *params = (kmer_t_meta*) calloc(1, sizeof(kmer_t_meta));
	tpl_hash all_tpls;
	GPtrArray *read_tpls = NULL, *branching_events = NULL;
	hash_table *ht = NULL;
	char *fn = NULL, *lib_fn_copy = NULL;
	hang_reads = g_ptr_array_sized_new(1024);
	tpls_await_branching = g_ptr_array_sized_new(32);
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

	g_ptr_array_free(tpls_await_branching, TRUE);
	show_msg(__func__, "Merging %d templates by pairs and overlapping ...\n",
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
	save_read_status(ht);
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

void process_only(char *junc_fn, char *pair_fa, char *hash_fn) {
	GPtrArray *all_tpls = g_ptr_array_sized_new(32);
	GPtrArray *all_junctions = g_ptr_array_sized_new(32);
	read_juncs_from_file(junc_fn, pair_fa, all_tpls, all_junctions);
	hash_table *ht = load_k_hash(hash_fn);
	filter_junctions(all_junctions, all_tpls, ht);
	process_graph(all_tpls, all_junctions, ht, kmer_out);
}

void cluster_set_start_reads(GPtrArray *junctions) {
	junction *j = NULL;
	int i = 0;
	for (i = 0; i < junctions->len; i++) {
		j = (junction*) g_ptr_array_index(junctions, i);
		j->branch_tpl->start_read = j->main_tpl->start_read;
	}
	for (i = 0; i < junctions->len; i++) {
		j = (junction*) g_ptr_array_index(junctions, i);
		j->branch_tpl->start_read = j->main_tpl->start_read;
	}
	for (i = 0; i < junctions->len; i++) {
		j = (junction*) g_ptr_array_index(junctions, i);
		j->branch_tpl->start_read = j->main_tpl->start_read;
	}
	for (i = 0; i < junctions->len; i++) {
		j = (junction*) g_ptr_array_index(junctions, i);
		j->branch_tpl->start_read = j->main_tpl->start_read;
	}
	for (i = 0; i < junctions->len; i++) {
		j = (junction*) g_ptr_array_index(junctions, i);
		j->branch_tpl->start_read = j->main_tpl->start_read;
	}
}

void cluster_ass(hash_table *ht, tpl_hash *all_tpls, bwa_seq_t *r1,
		bwa_seq_t *r2, int *cluster_ids_delimiter, int index, int cluster_id) {
	int i = 0;
	bwa_seq_t *r = NULL;
	pool *p = NULL;
	show_debug_msg(__func__, "Reads range: %d~%d \n",
			cluster_ids_delimiter[index - 1], cluster_ids_delimiter[index]);
	for (i = 0; i < ht->n_seqs; i++) {
		r = &ht->seqs[i];
		if (r->status != USED) {
			r->status = HANG;
			if (i >= cluster_ids_delimiter[index] && i
					< cluster_ids_delimiter[index + 1]) {
				r->status = FRESH;
			}
		}
	}

	CLUSTER_START_READ = r1;
	tpl *t = add_global_tpl(all_tpls, r1, "TEMPLATE", r1->len, 0);
	bwa_free_read_seq(1, t->ctg);
	t->ctg = new_seq(r1, r1->len, 0);
	t->len = t->ctg->len;
	t->start_read = r1;
	p = new_pool();
	init_pool(ht, p, t, kmer_len, N_MISMATCHES, 0);
	p_pool(__func__, p, NULL);
	if (p->reads->len > 0) {
		g_ptr_array_sort(p->reads, (GCompareFunc) cmp_reads_by_cursor);
		r = (bwa_seq_t*) g_ptr_array_index(p->reads, p->reads->len - 1);
		truncate_tpl(t, t->len - r->cursor, 1);
	}
	destroy_pool(p);

	p_tpl(t);
	refresh_tpl_reads(ht, t, N_MISMATCHES);

	// Assemble the cluster
	int to_con_left = ext_unit(ht, all_tpls, NULL, NULL, t, NULL, 0, 0);
	finalize_tpl(ht, all_tpls, t, 1, 0, 0);

	//	return;

	CLUSTER_START_READ = r2;
	t = add_global_tpl(all_tpls, r2, "TEMPLATE", r2->len, 0);
	bwa_free_read_seq(1, t->ctg);
	t->ctg = new_seq(r2, r2->len, 0);
	t->len = t->ctg->len;
	t->start_read = r2;

	p = new_pool();
	init_pool(ht, p, t, kmer_len, N_MISMATCHES, 1);
	p_pool(__func__, p, NULL);
	if (p->reads->len > 0) {
		g_ptr_array_sort(p->reads, (GCompareFunc) cmp_reads_by_cursor);
		r = (bwa_seq_t*) g_ptr_array_index(p->reads, 0);
		truncate_tpl(t, t->len - (r->len - r->cursor), 0);
	}
	destroy_pool(p);

	p_tpl(t);
	refresh_tpl_reads(ht, t, N_MISMATCHES);
	int to_con_right = ext_unit(ht, all_tpls, NULL, NULL, t, NULL, 0, 1);
	finalize_tpl(ht, all_tpls, t, 1, 0, 0);
}

int pe_cluster(int argc, char *argv[]) {
	if (!g_thread_supported())
		g_thread_init( NULL);
	kmer_id_mutex = g_mutex_new();
	hang_reads = g_ptr_array_sized_new(1024);
	tpls_await_branching = g_ptr_array_sized_new(32);

	int c = 0;
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

	// Read template FASTA
	uint64_t n_clusters = 0;
	bwa_seq_t *clusters = load_reads(argv[optind + 2], &n_clusters);
	bwa_seq_t *r = &clusters[0];

	// Read the cluster definition
	char buf[1000];
	char *attr[32];
	int *starts = (int*) calloc(n_clusters / 2, sizeof(int));
	int *ends = (int*) calloc(n_clusters / 2, sizeof(int));
	int *ids = (int*) calloc(n_clusters / 2, sizeof(int));

	attr[0] = strtok(r->name, "_");
	int cluster_index = 0, cluster_id = atoi(attr[0]);
	int start_index = 0, end_index = 0;
	FILE *defination = xopen(argv[optind + 1], "r");
	show_msg(__func__, "Reading template reads range... \n");
	while (fgets(buf, sizeof(buf), defination)) {
		attr[0] = strtok(buf, "\t");
		attr[1] = strtok(NULL, "\t");
		if (atoi(attr[0]) == cluster_id) {
			end_index = atoi(attr[1]);
			//show_debug_msg(__func__, "Cluster %d reads are: [%d, %d). \n",
			//		cluster_id, start_index, end_index, cluster_index);
			ids[cluster_index] = cluster_id;
			starts[cluster_index] = start_index;
			ends[cluster_index++] = end_index;
			if (cluster_index >= n_clusters / 2)
				break;
			r = &clusters[2 * cluster_index];
			attr[0] = strtok(r->name, "_");
			cluster_id = atoi(attr[0]);
		}
		start_index = atoi(attr[1]);
	}
	fclose(defination);

	// Pick raw reads and hash
	FILE *cluster_fa = xopen("tmp.fa", "w");
	FILE *all_fa = xopen(argv[optind], "r");
	int *cluster_ids_delimiter = (int*) calloc(n_clusters / 2 + 1, sizeof(int));
	int s = starts[0], e = ends[0], line_no = 0, id = 0, index = 0;
	while (fgets(buf, sizeof(buf), all_fa)) {
		if (line_no >= s && line_no < e) {
			if (line_no % 2 == 0) {
				fprintf(cluster_fa, ">%d cluster=%d\n", id++, ids[index]);
			} else {
				fprintf(cluster_fa, buf);
			}
			if (line_no % 1000000 == 0) {
				show_msg(__func__, "Looping to line %d...\n", line_no);
			}
			if (index >= n_clusters / 2)
				break;
		}
		if (line_no >= end_index)
			break;
		line_no++;
		if (line_no == e) {
			s = starts[++index];
			e = ends[index];
			cluster_ids_delimiter[index] = id;
		}
	}
	show_msg(__func__, "Saved %d reads to tmp.fa \n", id);
	fclose(all_fa);
	fclose(cluster_fa);

	char cmd[1000];
	sprintf(cmd, "./peta k_hash -k 11 -l 100 -i 2 -b 19 -s 4 tmp.fa");
	system(cmd);

	// Load hash table
	tpl_hash all_tpls;
	hash_table *ht = load_k_hash("tmp.fa");

	// Do the assembly
	bwa_seq_t *r1 = NULL, *r2 = NULL;
	int i = 0;
	for (i = 0; i < n_clusters / 2; i++) {
		r1 = &clusters[2 * i];
		r2 = &clusters[2 * i + 1];
		cluster_ass(ht, &all_tpls, r1, r2, cluster_ids_delimiter, i, ids[i]);
		if (i % 100 == 0)
			show_msg(__func__, "%d/%d: Cluster %d done \n", i, n_clusters / 2,
					ids[i]);
	}

	free(starts);
	free(ends);

	for (i = 0; i < n_clusters; i++) {
		r1 = &clusters[i];
		if (r1->status != USED)
			r1->status = FRESH;
	}

	// Merging
	kmer_hash tpl_kmer_hash;
	show_msg(__func__, "Merging %d templates by pairs and overlapping ...\n",
			all_tpls.size());
	merge_together_tpls(&all_tpls);
	iter_merge(ht, &all_tpls, &tpl_kmer_hash);

	char *fn = get_output_file("paired.fa", kmer_out);
	FILE *contigs = xopen(fn, "w");
	GPtrArray *read_tpls = hash_to_array(&all_tpls);
	all_tpls.clear();

	save_tpls(read_tpls, contigs, 0, 0, 0);
	save_read_status(ht);
	fflush(contigs);
	fclose(contigs);
	free(fn);

	// Build graph and run EM
	GPtrArray *branching_events = g_ptr_array_sized_new(BUFSIZ);
	get_junction_arr(read_tpls, branching_events);
	cluster_set_start_reads(branching_events);
	g_ptr_array_sort(branching_events, (GCompareFunc) cmp_junc_by_id);
	fn = get_output_file("paired.junctions", kmer_out);
	store_features(fn, branching_events, read_tpls);
	free(fn);

	show_msg(__func__, "Reloading the hash table ... \n");
	reload_table(ht, "tmp.fa");
	process_graph(read_tpls, branching_events, ht, kmer_out);
	destroy_ht(ht);

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
	if (optind + 1 > argc) {
		show_msg(__func__, "Parameters error! \n");
		return 1;
	}
	if (!g_thread_supported())
		g_thread_init( NULL);
	kmer_id_mutex = g_mutex_new();

	ext_by_kmers_core(argv[optind], argv[optind + 1]);

	clock_gettime(CLOCK_MONOTONIC, &kmer_finish_time);
	show_msg(__func__, "Done: %.2f sec\n",
			(float) (kmer_finish_time.tv_sec - kmer_start_time.tv_sec));
	return 0;
}
