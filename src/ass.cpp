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
	t->start_kmer = atol(start_read->name);
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
/**
 int vld_junc_by_mates(tpl *main_tpl, tpl *branch_tpl, GPtrArray *junc_reads,
 hash_map *hm, const int con_pos, const int ori) {
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
/**
 is_valid = vld_tpl_mates(branch_tpl, main_tpl, start, end, MIN_PAIRS);
 if (!is_valid) {
 g_ptr_array_sort(junc_reads, (GCompareFunc) cmp_reads_by_name);
 /**
 show_debug_msg(__func__, "Tag 1 \n");
 p_readarray(main_tpl->reads, 1);
 show_debug_msg(__func__, "Tag 2\n");
 p_readarray(junc_reads, 1);
 **/
/**
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
 **/
/**
 * To connect to templates with length smaller than k (25 by default)
 * The input template short_tpl is supposed to be shorter than k
 */
/**
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
/**
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
 **/

/**
 * If the template reaches some kmer which is used, stop the extension and add a branching event
 * The branch template make be marked 'dead' when ori is 1, then its junctions will be removed in kmer_ext_thread
 */
/**
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
/**
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
/**
 //if (locus > existing->len - hm->o->k)
 //	continue;
 con_pos = ori ? (locus + 1) : (locus + hm->o->k - 1);
 exist_ori = ori ? 0 : 1;

 // If extending to the left, while its right is not connected and it is too short, ignore
 if (ori && (!branch->b_juncs || branch->b_juncs->len == 0)
 && branch->len <= 3 * hm->o->k) {
 branch->alive = 0;
 break;
 }

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
/**
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
/**
 existing = parent_existing;
 con_pos = parent_locus;
 }
 }
 /**
 // If the branch is can be merged to main template, skip
 // Only check when it is extending to the left and
 //		its right end is connected.
 if (branch->len < 100 && branch->r_tail) {
 /**
 show_debug_msg(__func__, "Main template: %d; pos: %d \n",
 existing->id, con_pos);
 show_debug_msg(__func__, "Branch template: %d; ori: %d \n",
 branch->id, exist_ori);
 p_ctg_seq("Right Tail", branch->r_tail);
 **/
/**
 right_junc = (junction*) g_ptr_array_index(branch->b_juncs, 0);
 // If right and left connections are too close, just ignore.
 if (right_junc->main_tpl == existing) {
 if (get_abs(right_junc->locus - con_pos) <= N_MISMATCHES
 && branch->len <= N_MISMATCHES) {
 branch->alive = 0;
 break;
 }
 // If the branch is likely to be merged to the main template, try it
 loop_len = get_abs(branch->len - (right_junc->locus
 - con_pos));
 if (loop_len <= N_MISMATCHES) {
 on_main = branch_on_main(existing->ctg, branch->ctg,
 con_pos, (branch->len / hm->o->k + 2) * 3,
 exist_ori);
 valid = on_main ? 0 : 1;
 if (!valid) {
 // Mark it as 'dead', will be destroyed in kmer_ext_thread.
 branch->alive = 0;
 break;
 }
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
/**
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
 **/

/**
 * If next kmer is used by another template, try to connect to it.
 * Return: 0: cannot be connected
 *         1: forward connected
 *         2: reverse complement connected
 **/
/**
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
 **/

/**
 * Extend a template until no next kmer
 */
int kmer_ext_tpl(hash_table *ht, tpl_hash *all_tpls, tpl *t, bwa_seq_t *query,
		const int ori) {
	int max_c = -1, *counters = NULL, weight = 0, con_existing = 0;
	int max_c_all = 0, *counters_all = NULL;
	pool *p = NULL;
	bwa_seq_t *tail = new_seq(query, query->len, 0);

	show_debug_msg(__func__,
			"------ Started extending tpl %d to ori %d... ------\n", t->id, ori);
	p_query(__func__, tail);
	p_ctg_seq("TEMPLATE", t->ctg);

	p = new_pool();
	next_pool(ht, p, t, tail, N_MISMATCHES, ori);
	if (!ori)
		correct_tpl_base(p, t, tail->len);
	//p_pool("INITIAL_POOL", p, NULL);
	while (1) {
		max_c = get_next_char(p, t, ori);
		if (max_c == -1) {
			show_debug_msg(__func__, "No hits, stop ori %d: [%d, %d] \n", ori,
					t->id, t->len);
			break;
		}

		/**
		if (!ori) {
			p_query(__func__, tail);
			p_pool(__func__, p, NULL);
			show_debug_msg(__func__, "Next char: %c \n", "ACGTN"[max_c]);
			p_ctg_seq("TEMPLATE", t->ctg);
		}
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
	destroy_pool(p);
	return con_existing;
}

/**
 * Extend a kmer
 */
void *kmer_ext_thread(gpointer data, gpointer thread_params) {
	tpl *t = NULL;
	int round_1_len = 0, round_2_len = 0, connected = 0, ori = 1, flag = 0;
	int round_1_n_reads = 0;
	uint64_t read_id = 0;
	uint32_t i = 0;
	kmer_counter *counter = NULL;
	bwa_seq_t *read = NULL, *query = NULL, *seqs = NULL;
	junction *jun = NULL;

	kmer_t_meta *params = (kmer_t_meta*) thread_params;
	tpl_hash *all_tpls = params->all_tpls;
	hash_table *ht = params->ht;

	seqs = ht->seqs;
	counter = (kmer_counter*) data;
	read_id = counter->kmer;
	read = &seqs[read_id];

	//show_debug_msg(__func__, "============= %s: %" ID64 " ============ \n",
	//		read->name, counter->count);

	if (is_biased_q(read) || counter->count < 2 || is_repetitive_q(read)
			|| is_biased_q(read) || read->status != FRESH) {
		return NULL;
	}
	// Make a clone of the original starting read, which is global
	query = new_seq(read, kmer_len, read->len - kmer_len);

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

	// Extend to the right first
	connected = kmer_ext_tpl(ht, all_tpls, t, query, 0);
	round_1_len = t->len;
	round_1_n_reads = t->reads->len;
	show_debug_msg(__func__, "tpl %d with length: %d \n", t->id, t->len);

	bwa_free_read_seq(1, query);
	query = new_seq(read, kmer_len, 0);
	// Its reverse complement is Connected to an existing template
	if (connected == 2) {
		ori = 0;
		switch_fr(query);
	}

	// Then extend to the left
	connected |= kmer_ext_tpl(ht, all_tpls, t, query, ori);
	upd_locus_on_tpl(t, round_1_len, round_1_n_reads);
	g_ptr_array_sort(t->reads, (GCompareFunc) cmp_reads_by_contig_locus);
	//p_readarray(t->reads, 1);
	show_debug_msg(__func__, "tpl %d with length: %d \n", t->id, t->len);
	p_ctg_seq("TEMPLATE", t->ctg);

	if (!t->alive || (t->len <= query->len && (!t->b_juncs || t->b_juncs->len
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
		//upd_tpl_jun_locus(t, branching_events, opt->k);
	}
	bwa_free_read_seq(1, query);
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
		//if (i >= 100100)
		//break;
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
