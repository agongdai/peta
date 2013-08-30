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
GPtrArray *branching_events = NULL;

bwa_seq_t *TEST = NULL;

GMutex *kmer_id_mutex;
struct timespec kmer_start_time, kmer_finish_time;

void add_a_junction(tpl *main_tpl, tpl *branch_tpl, bwa_seq_t *connector,
		int locus, int ori, int weight) {
	// Indicating the templates are in-connect, cannot be reverse-complement
	main_tpl->in_connect = 1;
	branch_tpl->in_connect = 1;
	junction *new_j = new_junction(main_tpl, branch_tpl, connector, locus, ori,
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
	uint64_t id = 0;
	tpl *t = NULL;
	GPtrArray *tpls = g_ptr_array_sized_new(all_tpls->size());
	for (tpl_hash::iterator m = all_tpls->begin(); m != all_tpls->end(); ++m) {
		id = m->first;
		t = (tpl*) m->second;
		if (t->alive) {
			//p_tpl_reads(t);
			g_ptr_array_add(tpls, t);
		} else {
			destroy_tpl(t);
		}
		/**
		 show_debug_msg(__func__, "Tails of template %d\n", t->id);
		 p_ctg_seq(__func__, t->l_tail);
		 p_ctg_seq(__func__, t->r_tail);
		 **/
	}
	g_ptr_array_sort(tpls, (GCompareFunc) cmp_tpl_by_id);
	return tpls;
}

/**
 * Use read-length tail to search,
 * 	find those templates could be connected to current branch
 */
GPtrArray *find_connected_reads(hash_table *ht, tpl_hash *all_tpls,
		tpl *branch, const int ori) {
	bwa_seq_t *tail = NULL, *r = NULL, *branch_seq = NULL, *m_seq = NULL;
	index64 main_id = 0;
	int read_len = ht->o->read_len, m_seq_len = 0, ol_len = 0, n_mis = 0;
	int i = 0, j = 0, branch_main_same_ori = 0;
	ubyte_t x = 0, y = 0, *s = NULL;
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

	//	if (branch->id == 4) {
	//		p_query("TEST", TEST);
	//	}
	//char *test = (char*) malloc(sizeof(char) * 16);
	//sprintf(test, "SOME");
	//tail->name = test;
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
		if (r->contig_id != branch->id) {
			tpl_hash::iterator it = all_tpls->find(r->contig_id);
			if (it != all_tpls->end()) {
				m_tpl = (tpl*) it->second;
				for (j = 0; j < r->len; j++) {
					if ((j + r->contig_locus) <= m_tpl->len - 1 && r->seq[j]
							!= m_tpl->ctg->seq[j + r->contig_locus])
						n_mis++;
				}
				if (n_mis >= 8)
					branch_main_same_ori = r->rev_com ? 1 : 0;
				else
					branch_main_same_ori = r->rev_com ? 0 : 1;
				break;
			}
		}
	}
	//show_debug_msg(__func__, "branch_main_same_ori: %d \n",
	//		branch_main_same_ori);
	for (i = 0; i < hits->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(hits, i);
		if (r->contig_id == branch->id) {
			r->cursor = -1;
			r->pos = IMPOSSIBLE_NEGATIVE;
			continue;
		}
		n_mis = -1;
		if (ori) {
			ol_len = r->len - r->pos;
			if (branch->len < ol_len) {
				r->cursor = -1;
				r->pos = IMPOSSIBLE_NEGATIVE;
				continue;
			}
			branch_seq = new_seq(branch->ctg, ol_len, 0);
			n_mis = seq_ol(r, branch_seq, ol_len, N_MISMATCHES);
			if (n_mis >= 0)
				r->cursor = r->pos - 1;
		} else {
			ol_len = r->pos + tail->len;
			if (branch->len < ol_len) {
				r->cursor = -1;
				r->pos = IMPOSSIBLE_NEGATIVE;
				continue;
			}
			branch_seq = new_seq(branch->ctg, ol_len, branch->len - ol_len);

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
			// Reset the rev_com flag on used reads
			r->cursor = -1;
			r->pos = IMPOSSIBLE_NEGATIVE;
			if (!branch_main_same_ori)
				r->rev_com = r->rev_com ? 0 : 1;
		}
		bwa_free_read_seq(1, branch_seq);
	}

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
 * To connect to templates with length smaller than k (25 by default)
 * The input template short_tpl is supposed to be shorter than k
 */

tpl *connect_to_small_tpl(hash_table *ht, uint64_t query_int, tpl *branch_tpl,
		tpl *short_tpl, int *parent_locus, int *borrow_bases, int ori) {
	uint32_t i = 0;
	uint64_t kmer_int = 0;
	int max_len = ht->o->k - 1;
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

	for (i = 0; i <= junc_seq->len - ht->o->k; i++) {
		kmer_int = get_kmer_int(junc_seq->seq, i, 1, ht->o->k);

		/**
		 show_debug_msg(__func__, "Debug %d: left_len: %d; right_len: %d \n", i,
		 left_seq->len, right_seq->len);
		 bwa_seq_t *debug = get_kmer_seq(kmer_int, 25);
		 p_query(__func__, debug);
		 bwa_free_read_seq(1, debug);
		 **/

		if (kmer_int == query_int) {
			bwa_free_read_seq(1, junc_seq);
			junc_seq = NULL;
			if (ori) {
				if (i >= left_seq->len) {
					*parent_locus = i - left_seq->len;
					parent = short_tpl;
				} else {
					*borrow_bases = ht->o->k - 1 - i;
					*parent_locus = left_j->locus - *borrow_bases;
					parent = left_j->main_tpl;
				}
			} else {
				if (i > left_seq->len) {
					*parent_locus = i - left_seq->len;
					parent = short_tpl;
				} else {
					*borrow_bases = ht->o->k - short_tpl->len - (left_seq->len
							- i) - 1;
					*parent_locus = right_j->locus + *borrow_bases;
					parent = right_j->main_tpl;
				}
			}
			break;
		}
	}
	if (junc_seq)
		bwa_free_read_seq(1, junc_seq);
	bwa_free_read_seq(1, left_seq);
	bwa_free_read_seq(1, right_seq);
	return parent;
}

/**
 * Check whether the branch can be merged to the main template
 */
int merge_branch_to_main(hash_table *ht, tpl *branch) {
	junction *exist_junc = NULL, *right = NULL, *left = NULL;
	int i = 0;
	int loop_len = 0, l_len = 0, r_len = 0, t_len = 0;
	int on_main = 0, merged = 0, similar = 0;
	int exist_ori = 0;
	tpl *main_tpl = NULL;
	bwa_seq_t *main_seq = NULL, *sub = NULL, *r = NULL;
	// If right and left connections are too close, just ignore.
	if (!branch->b_juncs || branch->b_juncs->len <= 0)
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
				while (branch->reads->len)
					g_ptr_array_remove_index_fast(branch->reads, 0);
				branch->alive = 0;
				show_debug_msg(__func__,
						"Ignored the template [%d, %d] because too short\n",
						branch->id, branch->len);
			}
			bwa_free_read_seq(1, sub);
		}
		bwa_free_read_seq(1, main_seq);
	}
	return merged;
}

/**
 * During extension, if it reaches some read which is used already, try to connect to it.
 */
int connect_by_full_reads(hash_table *ht, tpl_hash *all_tpls, tpl *branch,
		int ext_len, const int ori) {
	GPtrArray *con_reads = NULL, *junc_reads = NULL;
	index64 i = 0, j = 0;
	bwa_seq_t *r = NULL, *tail = NULL;
	tpl *main_tpl = NULL;
	// Positions, lengths, etc.
	int locus = 0, con_pos = 0, exist_ori = 0, loop_len = 0;
	int copy_start = 0, copy_end = 0;
	// Indicators, etc
	int on_main = 0, valid = 0, weight = 0;
	int connected = 0, is_rev = 0;
	// If the branch is reverse complement connected, the direction needs to be switch
	int adj_ori = 0;
	int max_trial = 0;
	junction *exist_junc = NULL;

	// If extending to the left, and it's not connected to any template, mark it 'dead'
	if (ori && (!branch->b_juncs || branch->b_juncs->len == 0) && branch->len
			<= ht->o->read_len) {
		branch->alive = 0;
		return 0;
	}
	//p_query(__func__, TEST);
	con_reads = find_connected_reads(ht, all_tpls, branch, ori);
	//show_debug_msg(__func__, "Connecting reads: \n");
	//p_readarray(con_reads, 1);
	max_trial = con_reads->len > 4 ? 4 : con_reads->len;
	//p_query(__func__, TEST);
	p_tpl(branch);
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
		if (main_tpl->len <= ht->o->k || main_tpl == branch) {
			continue;
		}

		// If the branch is reversed complement in last round,
		//	but it is not successfully connected, need to reverse it back for this round.
		if (is_rev) {
			switch_fr(branch->ctg);
			is_rev = 0;
		}

		//p_query("CONNECTOR", r);
		//p_tpl(main_tpl);
		//p_tpl(branch);
		/**
		 * When connecting, need to check whether it's reverse-complement or not
		 * So need to get the subseq on branch and main to compare base-by-base
		 */
		tail = ori ? get_tail(branch, r->len - r->cursor - 1, 1) : get_tail(
				branch, r->cursor, 0);
		//p_query("BRANCH TAIL", tail);
		p_query("CONNECTOR", r);
		con_pos = ori ? r->contig_locus + 1 + r->cursor : r->contig_locus;
		//show_debug_msg(__func__, "CONPOS: %d\n", con_pos);

		if (!similar_bytes(tail->seq, main_tpl->ctg->seq + con_pos, kmer_len,
				N_MISMATCHES + 3)) {
			//p_tpl(branch);
			//p_tpl(main_tpl);
			show_debug_msg(__func__,
					"Branch [%d, %d] is reverse complemented \n", branch->id,
					branch->len);
			// Reverse the branch and the direction
			set_rev_com(branch->ctg);
			adj_ori = ori ? 0 : 1;
			switch_fr(branch->ctg);
			is_rev = 1;
		}
		//p_query("TESTING", TEST);
		bwa_free_read_seq(1, tail);
		//p_tpl(branch);

		locus = r->contig_locus;
		con_pos = adj_ori ? (locus + r->len - r->cursor) : (locus);
		exist_ori = adj_ori ? 0 : 1;

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

		//p_ctg_seq("MAIN", main_tpl->ctg);
		//p_ctg_seq("BRAN", branch->ctg);
		find_reads_ahead(branch, ht->o->read_len, ext_len, &weight, ori);
		//p_query("TESTING", TEST);
		valid = weight >= MIN_JUNCTION_READS ? 1 : 0;
		//show_debug_msg(__func__, "BRANCH READS: %d \n", branch->reads->len);
		//p_readarray(branch->reads, 1);
		//show_debug_msg(__func__, "Junction reads: %d \n", weight);

		if (!valid) {
			show_debug_msg(__func__, "No enough junction reads. Ignore. \n");
			continue;
			//			}
		} // End of checking junction reads and pairs
		// If junction added, the junction reads should be added to branch later
		// Not only free them, but reset the them to FRESH
		//unhold_reads_array(junc_reads);

		// Trim the branch
		if (exist_ori) {
			branch->len -= r->cursor;
		} else {
			memmove(branch->ctg->seq, branch->ctg->seq + r->cursor + 1,
					sizeof(ubyte_t) * (r->len - r->cursor - 1));
			branch->len -= (r->len - r->cursor - 1);
		}
		branch->ctg->len = branch->len;
		//p_ctg_seq("TRUNCATED", branch->ctg);
		set_rev_com(branch->ctg);

		// If there is a small loop, erase it.
		if (branch->b_juncs && branch->b_juncs->len == 1) {
			exist_junc = (junction*) g_ptr_array_index(branch->b_juncs, 0);
			if (exist_junc->main_tpl == main_tpl && exist_junc->ori
					!= exist_ori) {
				loop_len = get_abs(con_pos - exist_junc->locus);
				show_debug_msg(__func__, "Erasing small loop... \n");
				p_junction(exist_junc);
				if (loop_len < ht->o->read_len) {
					// Copy the small loop to the head of the branch, adjust the connecting position
					copy_start = min(con_pos, exist_junc->locus);
					copy_end = max(con_pos, exist_junc->locus);
					show_debug_msg(__func__, "Copy start~end: %d, %d\n",
							copy_start, copy_end);
					if (adj_ori) {
						for (j = copy_end - 1; j >= copy_start; j--) {
							ext_con(branch->ctg, main_tpl->ctg->seq[j], 1);
						}
					} else {
						for (j = copy_start; j < copy_end; j++) {
							ext_con(branch->ctg, main_tpl->ctg->seq[j], 0);
						}
					}
					branch->len = branch->ctg->len;
					set_rev_com(branch->ctg);
					con_pos = exist_junc->locus;
				}
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
		break;
	} // End of connecting all probable templates
	//p_query(__func__, TEST);
	// If it is not connected and is reverse complemented, just get it back
	if (is_rev && !connected) {
		switch_fr(branch->ctg);
	}
	for (i = 0; i < con_reads->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(con_reads, i);
		r->pos = IMPOSSIBLE_NEGATIVE;
		r->cursor = -1;
	}
	//p_query(__func__, TEST);
	g_ptr_array_free(con_reads, TRUE);
	return connected;
}

/**
 * Extend a template until no next kmer
 */
int kmer_ext_tpl(hash_table *ht, tpl_hash *all_tpls, pool *p, tpl *t,
		bwa_seq_t *query, int to_try_connect, const int ori) {
	int max_c = -1, pre_c = -1, *counters = NULL, weight = 0, con_existing = 0;
	int max_c_all = 0, *counters_all = NULL;
	int connected = 0, consuming_pool = 0;
	int ext_len = 0, no_read_len = 0;
	GPtrArray *near_tpls = g_ptr_array_sized_new(4);
	bwa_seq_t *tail = new_seq(query, query->len, 0);

	show_debug_msg(__func__,
			"------ Started extending tpl %d to ori %d ... ------\n", t->id,
			ori);
	p_query(__func__, tail);
	get_nearby_tpls(t, near_tpls);
	//p_ctg_seq("TEMPLATE", t->ctg);
	//p_query(__func__, TEST);
	while (1) {
		// If the query is bad, consume the reads in the pool first,
		// 	some cases would be able to go through the bad query region and continue
		//	otherwise, need to stop current extension.
		if (is_bad_query(tail)) {
			if (!consuming_pool) {
				keep_paired_reads(ht, p, t);
				consuming_pool = 1;
				p_query(__func__, tail);
				p_ctg_seq(__func__, t->ctg);
				show_debug_msg(__func__, "Consuming the pool ...\n");
				p_pool("CURRENT POOL", p, NULL);
			}
		}

		// If the reads in the pool are already consumed
		if (consuming_pool && p->reads->len == 0) {
			if (is_bad_query(tail)) {
				show_debug_msg(__func__, "Repetitive tail, stop.\n");
				p_query(__func__, tail);
				mark_pool_reads_tried(p, t);
				break;
			} else {
				// If current pool is empty, but the tail is not bad,
				// then continue to extend
				consuming_pool = 0;
				// Align the tail and add reads to the pool
				next_pool(ht, p, t, tail, N_MISMATCHES, ori);
			}
		}

		// If any long-enough subsequence is not covered by any read, try to connect
		if (no_read_len >= ht->o->read_len) {
			show_debug_msg(__func__,
					"[%d, %d] %d bases not covered by a read. \n", t->id,
					t->len, no_read_len);
			if (to_try_connect)
				con_existing = connect_by_full_reads(ht, all_tpls, t, ext_len,
						ori);
			show_debug_msg(__func__, "No hits, stop ori %d: [%d, %d] \n", ori,
					t->id, t->len);
			break;
		}

		max_c = get_next_char(ht, p, near_tpls, t, ori);
		// If cannot extend, try to add mates into the pool
		if (max_c == -1) {
			//p_ctg_seq("TEMPLATE", t->ctg);
			find_hashed_mates(ht, p, t, tail->len + 1, N_MISMATCHES, ori);
			max_c = get_next_char(ht, p, near_tpls, t, ori);
			if (max_c == -1) {
				if (to_try_connect)
					con_existing = connect_by_full_reads(ht, all_tpls, t,
							ext_len, ori);
				show_debug_msg(__func__, "No hits, stop ori %d: [%d, %d] \n",
						ori, t->id, t->len);
				break;
			} else {
				p_query("TAIL", tail);
				show_debug_msg(__func__, "Added mates: ori %d \n", ori);
				p_pool("MATE_POOL", p, NULL);
			}
		}

		//p_query(__func__, TEST);
		//if (t->id == 19999 && ori == 0 && t->len > 1340 && t->len < 1400) {
			//p_query("TESTING", TEST);
		//	show_debug_msg(__func__,
		//			"Ori: %d, Template [%d, %d], Next char: %c \n", ori, t->id,
		//			t->len, "ACGTN"[max_c]);
		//	p_query(__func__, tail);
		//	p_ctg_seq("TEMPLATE", t->ctg);
		//	p_pool("CURRENT POOL", p, NULL);
		//}

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
		if (!consuming_pool && (t->len % 8 == 0 || p->reads->len <= 4))
			next_pool(ht, p, t, tail, LESS_MISMATCH, ori);

		if (t->len % 100 == 0)
			show_debug_msg(__func__, "Ori %d, tpl %d, length %d \n", ori,
					t->id, t->len);
		if (t->len > 100000) {
			p_tpl(t);
			err_fatal(__func__,
					"[ERROR] Too long contig [%d, %d]. Maybe some bug.\n",
					t->id, t->len);
		}
	}
	//p_query(__func__, TEST);
	g_ptr_array_free(near_tpls, TRUE);
	bwa_free_read_seq(1, tail);
	return con_existing;
}

int try_destroy_tpl(hash_table *ht, tpl_hash *all_tpls, tpl *t, int read_len) {
	int i = 0, j = 0;
	int is_valid = 1, with_pairs = 0;
	float branch_cov = 0.0, main_cov = 0.0;
	junction *jun = NULL, *jun_pointer_on_main = NULL;
	tpl *main_tpl = NULL;
	GPtrArray *near_tpls = NULL;
	bwa_seq_t *r = NULL;
	if (!t)
		return 1;
	unfrozen_tried(t);
	if (!t->alive)
		is_valid = 0;
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
					if (branch_cov >= main_cov * BRANCHING_COV)
						is_valid = 1;
				}
			}
		}
	}
	// At stage 1, only get templates longer than insert size
	//if (stage == 1 && t->len < ins_size) {
	//	is_valid = 0;
	//}
	with_pairs = has_pairs_on_tpl(ht, t, MIN_PAIRS);
	// At any stage, if the template is longer than insert size, require some pairs
	if (is_valid && t->len >= ins_size && !with_pairs)
		is_valid = 0;
	// At stage 3, obtain transcripts shorter than insert size; otherwise, must has pairs
	if (is_valid && t->len < ins_size && stage != 3 && !with_pairs) {
		if (!t->b_juncs || t->b_juncs->len == 0)
			is_valid = 0;
		else {
			reset_is_root(t);
			near_tpls = g_ptr_array_sized_new(4);
			near_tpls = get_nearby_tpls(t, near_tpls);
			//if (t->id == 2) {
			//	p_tpl_reads(t);
			//}
			if (!has_nearby_pairs(ht, near_tpls, t, MIN_PAIRS)) {
				is_valid = 0;
				show_debug_msg(__func__, "No enough pairs \n");
			}
			reset_is_root(t);
			g_ptr_array_free(near_tpls, TRUE);
		}
	}
	if (!is_valid) {
		//show_debug_msg(__func__, "Destroying template [%d, %d] \n", t->id,
		//		t->len);
		g_mutex_lock(kmer_id_mutex);
		all_tpls->erase(t->id);
		g_mutex_unlock(kmer_id_mutex);
		if (t->b_juncs) {
			for (i = 0; i < t->b_juncs->len; i++) {
				jun = (junction*) g_ptr_array_index(t->b_juncs, i);
				destroy_junction(jun);
				//jun->status = 1;
				g_ptr_array_remove_index_fast(branching_events,
						branching_events->len - 1);
			}
		}
		// The reads on it marked as TRIED
		destroy_tpl(t);
		return 1;
	}
	return 0;
}

/**
 * Perform recursive branching for a template
 */
void branching(hash_table *ht, tpl_hash *all_tpls, tpl *t, int mismatches,
		int ori) {
	bwa_seq_t *tail = NULL, *r = NULL, *branch_read = NULL, *query = NULL;
	int i = 0, j = 0, shift = 0, cursor = 0, pos = 0;
	int near_s = 0, near_e = 0;
	int con_pos = 0, n_junc_reads = 0;
	int exist_ori = ori, dead = 0;
	tpl *branch = NULL;
	pool *p = NULL;
	int *pos_reads = NULL;
	int near_cov = 0, least_ol_len = ht->o->read_len - 2 * UNIQUE_LEN;

	if (!t || !t->alive || !t->ctg)
		return;
	pos_reads = (int*) calloc(t->len, sizeof(int));
	for (i = 0; i < t->reads->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(t->reads, i);
		if (r->contig_locus >= 0 && r->contig_locus < t->len - 1)
			pos_reads[r->contig_locus]++;
	}
	printf("\n");
	show_debug_msg(__func__,
			"===== Branching template [%d, %d] to %s ===== \n", t->id, t->len,
			ori ? "left" : "right");
	if (least_ol_len < kmer_len)
		least_ol_len = kmer_len;
	for (i = 0; i <= t->len - least_ol_len; i++) {
		// If extending to right, get tail from the end
		shift = i;
		if (!t->ctg) {
			show_debug_msg(__func__,
					"[WARNING] Template [%d, %d] suddenly empty.", t->id,
					t->len);
			break;
		}
		tail = new_seq(t->ctg, least_ol_len, shift);
		//if (i == 308) {
		//	p_ctg_seq(__func__, t->ctg);
		//	p_query(__func__, tail);
		//}
		branch_read = check_branch_tail(ht, t, tail, shift, mismatches, FRESH,
				ori);

		if (branch_read) {
			// For later truncate the branch template
			cursor = branch_read->cursor;
			pos = branch_read->pos;
			//branch_read->pos = IMPOSSIBLE_NEGATIVE;

			//p_query("BRANCH_QUERY", branch_read);

			// Create a new template
			branch = blank_tpl(branch_read, branch_read->len, ori);

			g_mutex_lock(kmer_id_mutex);
			all_tpls->insert(make_pair<int, tpl*> ((int) branch->id,
					(tpl*) branch));
			g_mutex_unlock(kmer_id_mutex);

			// Check the junction reads
			p = new_pool();
			init_pool(ht, p, branch, kmer_len, mismatches, ori);
			//g_ptr_array_sort(p->reads, (GCompareFunc) cmp_reads_by_name);
			//keep_paired_reads(ht, p, branch);
			keep_good_cursors(p);
			//p_query(__func__, tail);
			//p_query("BRANCH_QUERY", branch_read);
			//show_debug_msg(__func__, "i: %d; CURSOR: %d\n", i, cursor);
			//p_pool(__func__, p, NULL);

			// If the branch coverage is too low, ignore
			//p_query(__func__, TEST);
			if (p->reads->len >= MIN_JUNCTION_READS) {
				near_s = ori ? i : i - ht->o->read_len - 1;
				near_e = ori ? i + ht->o->read_len - 1 : i;
				near_s = near_s < 0 ? 0 : near_s;
				near_e = near_e >= t->len ? t->len - 1 : near_e;
				near_cov = 0;
				for (j = near_s; j < near_e; j++) {
					near_cov += pos_reads[j];
				}
				//if (i == 695)
				//	p_tpl_reads(t);
				//show_debug_msg(__func__, "Start~End: [%d, %d] \n", near_s,
				//		near_e);
				//show_debug_msg(__func__, "near_cov: %d \n", near_cov);

				if (p->reads->len > near_cov * BRANCHING_THRE) {
					//p_query(__func__, tail);
					//p_ctg_seq(__func__, t->ctg);

					n_junc_reads = p->reads->len;
					mark_init_reads_used(ht, branch, branch_read, mismatches);
					if (ori)
						add2tpl(branch, branch_read, cursor);
					else
						add2tpl(branch, branch_read, 0 - cursor);

					// Add the branching junction first;
					// Later may add connection junction
					con_pos = ori ? i : cursor - pos + i;
					//show_debug_msg(__func__, "Connect position: %d\n", con_pos);
					if (ori) {
						branch->len = cursor + 1;
					} else {
						memmove(branch->ctg->seq, branch->ctg->seq + cursor,
								(branch->len - cursor) * sizeof(ubyte_t));
						branch->len -= cursor;
					}

					branch->ctg->len = branch->len;
					set_rev_com(branch->ctg);
					add_a_junction(t, branch, NULL, con_pos, exist_ori,
							n_junc_reads);
					set_tail(branch, t, con_pos, ht->o->read_len - 1, exist_ori);

					p_query(__func__, tail);
					p_query("BRANCH_QUERY", branch_read);
					show_debug_msg(__func__, "i: %d; CURSOR: %d\n", i, cursor);
					show_debug_msg(__func__, "Branching at %d \n", con_pos);
					p_pool(__func__, p, NULL);
					p_ctg_seq("AFTER TRUNCATE", branch->ctg);

					// Perform extension
					if (branch_read->rev_com)
						switch_fr(branch_read);
					query = ori ? new_seq(branch_read, kmer_len, 0) : new_seq(
							branch_read, kmer_len, branch_read->len - kmer_len);
					query->rev_com = 0;
					if (branch_read->rev_com)
						switch_fr(branch_read);
					printf("\n");
					show_debug_msg(
							__func__,
							"===== Branching template [%d, %d] to %s; Started %d ===== \n",
							t->id, t->len, ori ? "left" : "right", branch->id);
					kmer_ext_tpl(ht, all_tpls, p, branch, query, 1, ori);
					bwa_free_read_seq(1, query);
					merge_branch_to_main(ht, branch);
				}
			} else {
				reset_to_fresh(branch_read);
			}
			destroy_pool(p);
			unfrozen_tried(branch);
			set_rev_com(branch->ctg);
			if (ori)
				upd_locus_on_tpl(branch, 0, 0);
			refresh_tpl_reads(ht, branch, mismatches);
			p_tpl(branch);
			//p_query(__func__, TEST);

			dead = try_destroy_tpl(ht, all_tpls, branch, ht->o->read_len);
			if (!dead) {
				//p_tpl_reads(branch);
			} else {
				//show_debug_msg(__func__, "Template is destroyed.\n");
			}
		}
		bwa_free_read_seq(1, tail);
	}
	free(pos_reads);
}

void branching_long(hash_table *ht, tpl_hash *all_tpls) {
	uint64_t id = 0;
	tpl *t = NULL;
	for (tpl_hash::iterator m = all_tpls->begin(); m != all_tpls->end(); ++m) {
		id = m->first;
		t = (tpl*) m->second;
		branching(ht, all_tpls, t, N_MISMATCHES, 0);
		branching(ht, all_tpls, t, N_MISMATCHES, 1);
	}
}

/**
 * Extend a kmer
 */
void *kmer_ext_thread(gpointer data, gpointer thread_params) {
	tpl *t = NULL;
	int pre_len = 0, round_2_len = 0, connected = 0, ori = 1, flag = 0;
	int invalid = 0;
	int pre_n_reads = 0, iter = 0;
	uint64_t read_id = 0;
	uint32_t i = 0;
	kmer_counter *counter = NULL;
	bwa_seq_t *read = NULL, *query = NULL, *seqs = NULL;
	junction *jun = NULL;
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

	//if (kmer_ctg_id == 1)
	//	read = &seqs[181];
	//	if (kmer_ctg_id == 2)
	//		read = &seqs[3901683];

	printf("\n");
	show_debug_msg(__func__,
			"============= FRESH %s: %" ID64 " ============ \n", read->name,
			counter->count);
	//p_query(__func__, TEST);
	fresh_trial++;
	p_query(__func__, read);
	t = blank_tpl(read, read->len, 0);
	t->kmer_freq = counter->count;
	// Insert first, in case it connects to itself during extension
	g_mutex_lock(kmer_id_mutex);
	all_tpls->insert(make_pair<int, tpl*> ((int) t->id, (tpl*) t));
	g_mutex_unlock(kmer_id_mutex);

	if (counter->count > 1)
		mark_init_reads_used(ht, t, read, N_MISMATCHES);
	else
		add2tpl(t, read, 0);
	// Right->left->right->left...until not extendable
	// If it is connected to somewhere, simply stop
	while (iter++ < 2 && t->len > pre_len && (!t->b_juncs || t->b_juncs->len
			== 0)) {
		// Extend to the right first
		// Make a clone of the original starting read, which is global
		p = new_pool();
		init_pool(ht, p, t, kmer_len, N_MISMATCHES, 0);
		//p_query(__func__, query);
		//g_ptr_array_sort(p->reads, (GCompareFunc) cmp_reads_by_name);
		//p_pool("INITIAL_POOL", p, NULL);
		//exit(1);
		if (iter == 1 && p->reads->len == 0) {
			t->alive = 0;
			flag = 1;
			destroy_pool(p);
			break;
		}

		// The correction is done only once
		if (pre_len == 0)
			correct_init_tpl_base(p, t, kmer_len);

		query = get_tail(t, kmer_len, 0);
		connected = kmer_ext_tpl(ht, all_tpls, p, t, query,
				params->to_try_connect, 0);
		destroy_pool(p);
		pre_len = t->len;
		set_rev_com(t->ctg);
		//refresh_tpl_reads(ht, t, N_MISMATCHES);
		pre_n_reads = t->reads->len;
		bwa_free_read_seq(1, query);
		show_debug_msg(__func__, "tpl %d with length: %d \n", t->id, t->len);

		// Maybe marked as not alive in last extension
		if (!t->alive)
			break;

		// Its reverse complement is Connected to an existing template
		ori = 1;
		if (connected == 2) {
			ori = 0;
		}

		// Then extend to the left
		query = get_tail(t, kmer_len, ori);
		p = new_pool();
		init_pool(ht, p, t, kmer_len, N_MISMATCHES, ori);
		//p_query(__func__, query);
		//p_pool("INITIAL_POOL", p, NULL);

		connected |= kmer_ext_tpl(ht, all_tpls, p, t, query,
				params->to_try_connect, ori);
		set_rev_com(t->ctg);
		destroy_pool(p);
		bwa_free_read_seq(1, query);
		//upd_locus_on_tpl(t, pre_len, pre_n_reads);

		//p_tpl(t);
		// Still necessary because the hashing may not get all reads
		upd_locus_on_tpl(t, pre_len, pre_n_reads);

		// Maybe marked as not alive in last extension
		if (!t->alive)
			break;
	}
	unfrozen_tried(t);
	set_rev_com(t->ctg);
	// Reactive the TRIED reads to FRESH, for other starting reads
	invalid = try_destroy_tpl(ht, all_tpls, t, ht->o->read_len);
	if (!invalid) {
		//show_debug_msg(__func__, "Reads before refresh: %d \n", t->reads->len);
		refresh_tpl_reads(ht, t, N_MISMATCHES);
		merge_branch_to_main(ht, t);
		invalid = try_destroy_tpl(ht, all_tpls, t, ht->o->read_len);

		if (!invalid) {
			g_ptr_array_sort(t->reads, (GCompareFunc) cmp_reads_by_contig_locus);
			//if (t->id == 2)
			//p_readarray(t->reads, 1);
			//p_tpl_reads(t);
			p_tpl(t);
			show_debug_msg(
					__func__,
					"==== End of tpl %d with length: %d; reads: %d; Alive: %d ==== \n\n",
					t->id, t->len, t->reads->len, t->alive);
			branching(ht, all_tpls, t, LESS_MISMATCH, 0);
			branching(ht, all_tpls, t, LESS_MISMATCH, 1);
		}
	} else {
		show_debug_msg(__func__, "Template is destroyed.\n");
	}

	// If the read cannot be even extend one base, reset the read to fresh
	if (flag) {
		reset_to_fresh(read);
	}
	return NULL;
}

void kmer_threads(kmer_t_meta *params) {
	GThreadPool *thread_pool = NULL;
	uint64_t i = 0;
	hash_table *ht = params->ht;
	bwa_seq_t *r = NULL, *seqs = ht->seqs;
	kmer_counter *counter = NULL;
	GPtrArray *starting_reads = g_ptr_array_sized_new(ht->n_seqs);
	GPtrArray *low_reads = NULL;

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

	TEST = &seqs[607];

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
		if (i % 100000 == 0)
			show_msg(__func__, "Extending %" ID64 "-th read ... \n", i);
		counter = (kmer_counter*) g_ptr_array_index(starting_reads, i);
		//g_thread_pool_push(thread_pool, (gpointer) counter, NULL);
		kmer_ext_thread(counter, params);
		free(counter);
		//if (fresh_trial >= 1)
		//	break;
	}
	g_ptr_array_free(starting_reads, TRUE);
	show_msg(__func__, "%d templates are obtained. \n",
			params->all_tpls->size());

	show_debug_msg(__func__, "Remaining reads: \n");
	for (i = 0; i < ht->n_seqs; i++) {
		r = &ht->seqs[i];
		if (r->status != USED)
			p_query(__func__, r);
	}

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

	g_thread_pool_free(thread_pool, 0, 1);
}

void test_kmer_ext(kmer_t_meta *params) {
}

int merge_paired_tpls(hash_table *ht, tpl_hash *all_tpls) {
	tpl *left = NULL, *right = NULL, *t = NULL, *mt = NULL;
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
				p_tpl(mt);
				show_debug_msg(__func__, "OVERLAP: %d\n", ol);
				// At most 1 mismatch in 11bp
				if (ol >= ht->o->k && ol >= n_mis * ht->o->k) {
					if (merge_tpls(t, mt, ol, rev_com)) {
						// Update the t->tried
						mv_unpaired_to_tried(ht->seqs, t, kmer_ctg_id);
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
	int merge_iter = 0, id = 0, merged = 1;
	tpl *t = NULL;

	// Add all mates not on current template to t->tried.
	for (tpl_hash::iterator im = all_tpls->begin(); im != all_tpls->end(); ++im) {
		id = im->first;
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
		id = it->first;
		t = (tpl*) it->second;
		if (t->tried) {
			while (t->tried->len > 0)
				g_ptr_array_remove_index_fast(t->tried, 0);
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

void ext_by_kmers_core(char *lib_file, const char *solid_file) {
	FILE *contigs = NULL;
	kmer_t_meta *params = (kmer_t_meta*) calloc(1, sizeof(kmer_t_meta));
	tpl_hash all_tpls;
	GPtrArray *read_tpls = NULL;
	hash_table *ht = NULL;
	char *fn = NULL, *lib_fn_copy = NULL;
	kmer_hash tpl_kmer_hash;

	show_msg(__func__, "Library: %s \n", lib_file);
	lib_fn_copy = str_dup(lib_file);

	//test_smith_waterman(lib_file);
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
	iter_merge(ht, &all_tpls, &tpl_kmer_hash);

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

	g_ptr_array_sort(branching_events, (GCompareFunc) cmp_junc_by_id);
	clean_junctions(read_tpls, branching_events);
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
	char *attr[18];
	junction *jun = NULL;
	for (i = 0; i < n_ctgs; i++) {
		ctg = &seqs[i];
		t = new_tpl();
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
	hash_table *ht = load_k_hash(hash_fn);
	filter_junctions(all_junctions, all_tpls, ht);
	process_graph(all_tpls, all_junctions, ht, kmer_out);
}

int pe_kmer(int argc, char *argv[]) {

	//	process_only("../SRR097897_out/paired.junctions.nolen",
	//			"../SRR097897_out/paired.fa",
	//			"/home/carl/Projects/peta/rnaseq/Spombe/SRR097897/SRR097897.fa");
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
