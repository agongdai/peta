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
#include "read.h"
#include "pool.hpp"

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

/**
 * Initialize a template
 */
tpl *blank_tpl(bwa_seq_t *start_read, int len, int ori) {
	tpl *t = new_tpl();
	g_mutex_lock(kmer_id_mutex);
	t->id = kmer_ctg_id++;
	g_mutex_unlock(kmer_id_mutex);
	if (len < 0 || len > start_read->len) {
		t->ctg = new_seq(start_read, start_read->len, 0);
	} else {
		t->ctg = ori ? new_seq(start_read, len, 0) : new_seq(start_read, len,
				start_read->len - len);
	}
	t->len = t->ctg->len;
	t->start_read = start_read;
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

/**
 * Validate the junction by checking mate pairs.
 * Depending on ori and con_pos, only partial reads on the main template are counted
 */
int vld_junc_by_mates(tpl *main_tpl, tpl *branch_tpl, GPtrArray *junc_reads,
		hash_table *ht, const int con_pos, const int ori) {
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
	is_valid = vld_tpl_mates(branch_tpl, main_tpl, start, end, MIN_PAIRS);
	if (!is_valid) {
		g_ptr_array_sort(junc_reads, (GCompareFunc) cmp_reads_by_name);
		/**
		 show_debug_msg(__func__, "Tag 1 \n");
		 p_readarray(main_tpl->reads, 1);
		 show_debug_msg(__func__, "Tag 2\n");
		 p_readarray(junc_reads, 1);
		 **/
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
			junc_seq = NULL;
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
	if (junc_seq)
		bwa_free_read_seq(1, junc_seq);
	bwa_free_read_seq(1, left_seq);
	bwa_free_read_seq(1, right_seq);
	return parent;
}

/**
 * Use read-length tail to search,
 * 	find those templates could be connected to current branch
 */
GPtrArray *find_connected_reads(hash_table *ht, tpl_hash *all_tpls,
		tpl *branch, const int ori) {
	bwa_seq_t *tail = NULL, *r = NULL, *tail_shift = NULL;
	index64 main_id = 0;
	int read_len = ht->o->read_len;
	int i = 0;
	ubyte_t x = 0;
	GPtrArray *mains = g_ptr_array_sized_new(0);
	GPtrArray *hits = NULL;

	if (branch->len < ht->o->read_len)
		return mains;

	tail = ori ? new_seq(branch->ctg, read_len, 0) : new_seq(branch->ctg,
			read_len, branch->len - read_len);
	// Try ACGT four directions
	for (x = 0; x < 4; x++) {
		tail_shift = new_seq(tail, tail->len, 0);
		ext_que(tail_shift, x, ori);

		hits = find_both_fr_full_reads(ht, tail_shift, hits, N_MISMATCHES);
		for (i = 0; i < hits->len; i++) {
			r = (bwa_seq_t*) g_ptr_array_index(hits, i);
			if (r->status != USED)
				continue;
			main_id = r->contig_id;
			tpl_hash::iterator it = all_tpls->find(main_id);
			if (it != all_tpls->end()) {
				g_ptr_array_add(mains, r);
			}
		}
		bwa_free_read_seq(1, tail_shift);
	}
	mains = rm_dup_reads_on_tpl(mains);
	for (i = 0; i < mains->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(mains, i);
		p_query(__func__, r);
	}
	bwa_free_read_seq(1, tail);
	return mains;
}

/**
 * During extension, if it reaches some read which is used already, try to connect to it.
 */
int connect_by_full_reads(hash_table *ht, tpl_hash *all_tpls, tpl *branch,
		const int ori) {
	GPtrArray *con_reads = NULL, *junc_reads = NULL;
	index64 i = 0, j = 0;
	bwa_seq_t *r = NULL;
	tpl *main_tpl = NULL;
	// Positions, lengths, etc.
	int locus = 0, con_pos = 0, exist_ori = 0, loop_len = 0;
	// Indicators, etc
	int on_main = 0, valid = 0, weight = 0, borrow_bases = 0;
	int connected = 0;
	junction *r_junc = NULL;

	// If extending to the left, and it's not connected to any template, mark it 'dead'
	if (ori && (!branch->b_juncs || branch->b_juncs->len == 0) && branch->len
			<= 3 * ht->o->k) {
		branch->alive = 0;
		return 0;
	}

	con_reads = find_connected_reads(ht, all_tpls, branch, ori);
	for (i = 0; i < con_reads->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(con_reads, i);
		//p_query(__func__, r);
		tpl_hash::iterator it = all_tpls->find(r->contig_id);
		if (it == all_tpls->end()) {
			continue;
		}
		// The candidate template to connect
		main_tpl = (tpl*) it->second;
		locus = r->contig_locus;
		con_pos = ori ? (locus + 1) : (locus + ht->o->read_len - 1);
		exist_ori = ori ? 0 : 1;

		//p_tpl(main_tpl);

		// If right and left connections are too close, just ignore.
		if (branch->b_juncs && branch->b_juncs->len > 0) {
			r_junc = (junction*) g_ptr_array_index(branch->b_juncs, 0);
			if (r_junc->main_tpl == main_tpl) {
				// If all of them simply too short
				if (get_abs(r_junc->locus - con_pos) <= IGNORE_DIFF
						&& branch->len <= IGNORE_DIFF) {
					branch->alive = 0;
					break;
				}
				// If the branch is likely to be merged to the main template, try it
				loop_len = get_abs(branch->len - (r_junc->locus - con_pos));
				if (loop_len <= IGNORE_DIFF) {
					on_main = branch_on_main(main_tpl->ctg, branch->ctg,
							con_pos, (branch->len / ht->o->k + 2) * 3,
							exist_ori);
					valid = on_main ? 0 : 1;
					if (!valid) {
						// Mark it as 'dead', will be destroyed in kmer_ext_thread.
						branch->alive = 0;
						break;
					}
				}
			} // Check the bubble only if left and right junction connect to the same main_tpl
		} // End of checking the short 'bubble'
		junc_reads = find_junc_reads_w_tails(ht, main_tpl, branch, con_pos,
				(ht->o->read_len - SHORT_BRANCH_SHIFT) * 2, ori, &weight);
		valid = (junc_reads->len >= MIN_JUNCTION_READS) ? 1 : 0;
		//p_readarray(junc_reads, 0);
		if (!valid) {
			show_debug_msg(__func__,
					"No enough junction reads. Check mates now...\n");
			valid = vld_junc_by_mates(main_tpl, branch, junc_reads, ht,
					con_pos, ori);
			if (!valid) {
				show_debug_msg(__func__,
						"Not passed the pair validation. Free it later.\n");
				g_ptr_array_free(junc_reads, TRUE);
				continue;
			}
		} // End of checking junction reads and pairs
		g_ptr_array_free(junc_reads, TRUE);

		// Trim the branch
		if (branch->len < ht->o->read_len) {
			if (exist_ori)
				con_pos -= branch->len;
			else
				con_pos = locus + branch->len + 1;
			branch->len = 0;
			branch->ctg->len = 0;
			set_rev_com(branch->ctg);
		} else {
			// Make the branch not sharing a read length subseq with the main
			if (borrow_bases) {
				if (exist_ori) {
					con_pos -= borrow_bases;
					branch->len -= borrow_bases;
				} else {
					con_pos += borrow_bases;
					memmove(branch->ctg->seq, branch->ctg->seq + borrow_bases,
							sizeof(ubyte_t) * (branch->len - borrow_bases));
					branch->len -= borrow_bases;
				}
				//show_debug_msg(__func__, "Connecting position: %d\n",
				//		con_pos);
			} else {
				if (exist_ori) {
					con_pos -= (ht->o->read_len - 1);
				} else {
					con_pos += (ht->o->read_len - 1);
					memmove(branch->ctg->seq, branch->ctg->seq
							+ (ht->o->read_len - 1), sizeof(ubyte_t)
							* (branch->len - (ht->o->read_len - 1)));
				}
				branch->len -= (ht->o->read_len - 1);
			}
			branch->ctg->len = branch->len;
			//p_ctg_seq("TRUNCATED", branch->ctg);
			set_rev_com(branch->ctg);
		} // End of trimming the branch

		// If there is a small loop, erase it.
		if (ori && branch->b_juncs && branch->b_juncs->len == 1) {
			r_junc = (junction*) g_ptr_array_index(branch->b_juncs, 0);
			if (r_junc->main_tpl == main_tpl && r_junc->ori == 1) {
				loop_len = con_pos - r_junc->locus;
				show_debug_msg(__func__, "Erasing small loop... \n");
				p_junction(r_junc);
				if (loop_len > 0 && loop_len < ht->o->read_len) {
					// Copy the small loop to the head of the branch, adjust the connecting position
					for (j = r_junc->locus + loop_len - 1; j >= r_junc->locus; j--) {
						ext_con(branch->ctg, main_tpl->ctg->seq[j], 1);
					}
					branch->len = branch->ctg->len;
					con_pos = r_junc->locus;
				}
			}
		}

		// Finally! Go to add the junction!
		//p_tpl(branch);
		show_debug_msg(__func__,
				"Connect existing [%d, %d] to [%d, %d] at %d. \n", branch->id,
				branch->len, main_tpl->id, main_tpl->len, con_pos);
		set_tail(branch, main_tpl, con_pos, ht->o->read_len
				- SHORT_BRANCH_SHIFT, exist_ori);
		//p_ctg_seq("Right tail", branch->r_tail);
		//p_ctg_seq("Left  tail", branch->l_tail);
		add_a_junction(main_tpl, branch, 0, con_pos, exist_ori, weight);
		connected = 1;

	} // End of connecting all probable templates
	g_ptr_array_free(con_reads, TRUE);
	return connected;
}

/**
 * Extend a template until no next kmer
 */
int kmer_ext_tpl(hash_table *ht, tpl_hash *all_tpls, pool *p, tpl *t,
		bwa_seq_t *query, const int ori) {
	int max_c = -1, *counters = NULL, weight = 0, con_existing = 0;
	int max_c_all = 0, *counters_all = NULL;
	int connected = 0;
	bwa_seq_t *tail = new_seq(query, query->len, 0);

	show_debug_msg(__func__,
			"------ Started extending tpl %d to ori %d... ------\n", t->id, ori);
	//p_query(__func__, tail);
	//p_ctg_seq("TEMPLATE", t->ctg);

	while (1) {
		max_c = get_next_char(p, t, ori);
		// If cannot extend, try to add mates into the pool
		if (max_c == -1) {
			//p_ctg_seq("TEMPLATE", t->ctg);
			find_match_mates(ht, p, t, tail->len, LESS_MISMATCH, ori);
			max_c = get_next_char(p, t, ori);
			if (max_c == -1) {
				con_existing = connect_by_full_reads(ht, all_tpls, t, ori);
				show_debug_msg(__func__, "No hits, stop ori %d: [%d, %d] \n",
						ori, t->id, t->len);
				break;
			} else {
				show_debug_msg(__func__, "Added mates: ori %d \n", ori);
				p_query("TAIL", tail);
				p_pool("MATE_POOL", p, NULL);
			}
		}
		/**
		 p_query(__func__, tail);
		 show_debug_msg(__func__, "Next char: %c \n", "ACGTN"[max_c]);
		 p_ctg_seq("TEMPLATE", t->ctg);
		 p_pool(__func__, p, NULL);
		 **/
		ext_con(t->ctg, max_c, ori);
		t->len = t->ctg->len;
		ext_que(tail, max_c, ori);
		// If the overlapped region between t and r has too many mismatches, remove from pool
		rm_half_clip_reads(p, t, max_c, N_MISMATCHES, ori);
		forward(p, t, ori);
		next_pool(ht, p, t, tail, N_MISMATCHES, ori);

		if (t->len % 100 == 0)
			show_debug_msg(__func__, "Ori %d, tpl %d, length %d \n", ori,
					t->id, t->len);
	}
	bwa_free_read_seq(1, tail);
	return con_existing;
}

/**
 * Extend a kmer
 */
void *kmer_ext_thread(gpointer data, gpointer thread_params) {
	tpl *t = NULL;
	int pre_len = 0, round_2_len = 0, connected = 0, ori = 1, flag = 0;
	int pre_n_reads = 0;
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

	//show_debug_msg(__func__, "============= %s: %" ID64 " ============ \n",
	//		read->name, counter->count);

	if (is_biased_q(read) || has_n(read, 1) || counter->count < 2
			|| is_repetitive_q(read) || is_biased_q(read) || read->status
			!= FRESH) {
		return NULL;
	}

	show_debug_msg(__func__, "============= %s: %" ID64 " ============ \n",
			read->name, counter->count);
	p_query(__func__, read);
	t = blank_tpl(read, read->len, 0);
	t->kmer_freq = counter->count;
	// Insert first, in case it connects to itself during extension
	g_mutex_lock(kmer_id_mutex);
	all_tpls->insert(make_pair<int, tpl*> ((int) t->id, (tpl*) t));
	g_mutex_unlock(kmer_id_mutex);

	mark_init_reads_used(ht, t, read, N_MISMATCHES);
	// Right->left->right->left...until not extendable
	while (t->len > pre_len) {
		// Extend to the right first
		// Make a clone of the original starting read, which is global
		p = new_pool();
		init_pool(ht, p, t, kmer_len, N_MISMATCHES, 0);
		//p_pool("INITIAL_POOL", p, NULL);
		// The correction is done only once
		if (pre_len == 0)
			correct_tpl_base(p, t, kmer_len);

		query = get_tail(t, kmer_len, 0);
		connected = kmer_ext_tpl(ht, all_tpls, p, t, query, 0);
		destroy_pool(p);
		pre_len = t->len;
		pre_n_reads = t->reads->len;
		show_debug_msg(__func__, "tpl %d with length: %d \n", t->id, t->len);

		// Then extend to the left
		bwa_free_read_seq(1, query);
		query = get_tail(t, kmer_len, 1);
		// Its reverse complement is Connected to an existing template
		ori = 1;
		if (connected == 2) {
			ori = 0;
			switch_fr(query);
		}

		p = new_pool();
		init_pool(ht, p, t, kmer_len, N_MISMATCHES, ori);
		//p_pool("INITIAL_POOL", p, NULL);

		connected |= kmer_ext_tpl(ht, all_tpls, p, t, query, ori);
		destroy_pool(p);
		bwa_free_read_seq(1, query);
		upd_locus_on_tpl(t, pre_len, pre_n_reads);

		show_debug_msg(__func__, "tpl %d with length: %d \n", t->id, t->len);
		p_ctg_seq("TEMPLATE", t->ctg);
	}

	//g_ptr_array_sort(t->reads, (GCompareFunc) cmp_reads_by_contig_locus);
	//p_readarray(t->reads, 1);

	if (!t->alive || (t->len <= read->len && (!t->b_juncs || t->b_juncs->len
			< 2))) {
		g_mutex_lock(kmer_id_mutex);
		all_tpls->erase(t->id);
		g_mutex_unlock(kmer_id_mutex);
		if (t->b_juncs) {
			for (i = 0; i < t->b_juncs->len; i++) {
				jun = (junction*) g_ptr_array_index(t->b_juncs, i);
				destroy_junction(jun);
				g_ptr_array_remove_index_fast(branching_events,
						branching_events->len - 1);
			}
		}
		// The reads on it marked as TRIED
		destroy_tpl(t);
	} else {
		unfrozen_tried(t);
		//upd_tpl_jun_locus(t, branching_events, opt->k);
	}
	return NULL;
}

void kmer_threads(kmer_t_meta *params) {
	GThreadPool *thread_pool = NULL;
	uint64_t i = 0;
	hash_table *ht = params->ht;
	read_hash *rh = params->rh;
	GPtrArray *starting_reads = g_ptr_array_sized_new(rh->n_seqs);
	bwa_seq_t *r = NULL, *seqs = ht->seqs;
	kmer_counter *counter = NULL;

	for (i = 0; i < rh->n_seqs; i++) {
		r = &seqs[i];
		counter = (kmer_counter*) malloc(sizeof(kmer_counter));
		counter->kmer = i;
		counter->count = rh->similar_reads_count[i];
		g_ptr_array_add(starting_reads, counter);
	}

	destroy_rh(rh);
	show_msg(__func__, "Sorting %d initial reads... \n", starting_reads->len);
	g_ptr_array_sort(starting_reads, (GCompareFunc) cmp_kmers_by_count);
	show_msg(__func__, "Extending by reads...\n");
	thread_pool = g_thread_pool_new((GFunc) kmer_ext_thread, params, 1, TRUE,
			NULL);
	for (i = 0; i < starting_reads->len; i++) {
		if (i % 100000 == 0)
			show_debug_msg(__func__, "Extending %" ID64 "-th read... \n", i);
		counter = (kmer_counter*) g_ptr_array_index(starting_reads, i);
		kmer_ext_thread(counter, params);
		free(counter);
		//g_thread_pool_push(thread_pool, (gpointer) counter, NULL);
		if (kmer_ctg_id >= 100)
			break;
	}
	g_thread_pool_free(thread_pool, 0, 1);
	g_ptr_array_free(starting_reads, TRUE);
}

void test_kmer_ext(kmer_t_meta *params) {
}

void ext_by_kmers_core(char *lib_file, const char *solid_file) {
	read_hash *rh = NULL;
	FILE *contigs = NULL;
	kmer_t_meta *params = (kmer_t_meta*) calloc(1, sizeof(kmer_t_meta));
	tpl_hash all_tpls;
	GPtrArray *read_tpls = NULL;
	hash_table *ht = NULL;

	show_msg(__func__, "Library: %s \n", lib_file);

	ht = load_k_hash(lib_file);
	rh = load_read_hash(lib_file);

	clock_gettime(CLOCK_MONOTONIC, &kmer_finish_time);
	show_msg(__func__, "Done loading read hash table: %.2f sec\n",
			(float) (kmer_finish_time.tv_sec - kmer_start_time.tv_sec));

	branching_events = g_ptr_array_sized_new(BUFSIZ);
	params->ht = ht;
	params->rh = rh;
	params->all_tpls = &all_tpls;

	//test_kmer_ext(params);
	//exit(1);
	kmer_threads(params);
	// Start branching after the frequent kmers are consumed already.
	//start_branching(&all_tpls, params);

	show_msg(__func__, "Saving contigs: %.2f sec\n",
			(float) (kmer_finish_time.tv_sec - kmer_start_time.tv_sec));
	contigs = xopen(get_output_file("paired.fa", kmer_out), "w");
	read_tpls = hash_to_array(&all_tpls);
	save_tpls(read_tpls, contigs, 0, 0, 0);
	fflush(contigs);
	fclose(contigs);

	g_ptr_array_sort(branching_events, (GCompareFunc) cmp_junc_by_id);
	clean_junctions(branching_events);
	store_features(get_output_file("paired.junctions", kmer_out),
			branching_events, read_tpls);
	process_graph(read_tpls, branching_events, ht);

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
	process_graph(all_tpls, all_junctions, ht);
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
