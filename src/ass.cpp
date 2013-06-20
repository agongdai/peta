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
#include "edge.h"
#include "kmers.hpp"
#include "ass.hpp"
#include "junction.hpp"
#include "rnaseq.h"
#include "pechar.h"
#include "utils.h"
#include "peseq.h"

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

void add_a_junction(edge *main_tpl, edge *branch_tpl, uint64_t kmer, int locus,
		int ori, int weight) {
	// Indicating the templates are in-connect, cannot be reverse-complement
	main_tpl->in_connect = 1;
	branch_tpl->in_connect = 1;
	junction *new_j = new_junction(main_tpl, branch_tpl, kmer, locus, ori,
			weight);
	g_mutex_lock(kmer_id_mutex);
	g_ptr_array_add(branching_events, new_j);
	g_mutex_unlock(kmer_id_mutex);
}

edge *blank_edge(uint64_t query_int, int kmer_len, int init_len, int ori) {
	bwa_seq_t *kmer = NULL;
	edge *eg = new_eg();
	g_mutex_lock(kmer_id_mutex);
	eg->id = kmer_ctg_id++;
	g_mutex_unlock(kmer_id_mutex);
	// Get a copy of the kmer
	kmer = get_kmer_seq(query_int, kmer_len);
	if (ori)
		eg->ctg = new_seq(kmer, init_len, 0);
	else
		eg->ctg = new_seq(kmer, init_len, kmer->len - init_len);
	eg->len = eg->ctg->len;
	eg->tid = atoi(kmer->name);
	eg->start_kmer = query_int;
	bwa_free_read_seq(1, kmer);
	return eg;
}

GPtrArray *hash_to_array(tpl_hash *all_tpls) {
	uint64_t id = 0;
	edge *eg = NULL;
	GPtrArray *tpls = g_ptr_array_sized_new(all_tpls->size());
	for (tpl_hash::iterator m = all_tpls->begin(); m != all_tpls->end(); ++m) {
		id = m->first;
		eg = (edge*) m->second;
		g_ptr_array_add(tpls, eg);
	}
	return tpls;
}

gint cmp_kmers_by_count(gpointer a, gpointer b) {
	kmer_counter *c_a = *((kmer_counter**) a);
	kmer_counter *c_b = *((kmer_counter**) b);
	return ((c_b->count) - c_a->count);
}

void mark_tpl_kmers_used(edge *eg, hash_map *hm, const int kmer_len,
		const int shift) {
	int i = 0;
	uint64_t query_int = 0;
	if (eg->len < kmer_len)
		return;
	for (i = 0; i <= eg->len - kmer_len; i++) {
		query_int = get_kmer_int(eg->ctg->seq, i, 1, kmer_len);
		mark_kmer_used(query_int, hm, eg->id, i + shift);
	}
}

void mark_tpl_kmers_fresh(edge *eg, hash_map *hm, const int kmer_len) {
	int i = 0;
	uint64_t query_int = 0;
	for (i = 0; i < eg->len - kmer_len; i++) {
		query_int = get_kmer_int(eg->ctg->seq, i, 1, kmer_len);
		mark_kmer_not_used(query_int, hm);
	}
}

void cal_coverage(edge *eg, hash_map *hm) {
	int i = 0, *counters = NULL;
	uint64_t query_int = 0;
	float sum = 0;
	if (!eg || eg->len <= 0)
		return;
	for (i = 0; i < eg->len - kmer_len; i++) {
		query_int = get_kmer_int(eg->ctg->seq, i, 1, kmer_len);
		counters = count_next_kmers(hm, query_int, 1, 0);
		sum += counters[get_max_index(counters)];
	}
	eg->coverage = sum / eg->len;
}

/**
 * Called after try_short_tpl_byte_ext.
 * If the edge after short extension is too short but right connected,
 * 	borrow some sequence from the right connected edge.
 */
void borrow_existing_seq(bwa_seq_t *short_seq, edge *connected_tpl,
		const int shift, const int kmer_len, const int max_len, const int ori) {
	bwa_seq_t *tail = NULL;
	int borrow_len = 0, existing_ori = 0, existing_shift = 0;
	if (short_seq->len >= max_len)
		return;
	existing_ori = ori ? 0 : 1;
	existing_shift = existing_ori ? (shift + kmer_len - 1) : shift;
	tail = cut_edge_tail(connected_tpl, existing_shift, max_len, existing_ori);
	borrow_len = max_len - short_seq->len;
	borrow_len = borrow_len > tail->len ? tail->len : borrow_len;
	// The space has be allocated to max_len in try_short_tpl_byte_ext already.
	if (existing_ori)
		memcpy(short_seq->seq + short_seq->len, tail->seq, borrow_len
				* sizeof(ubyte_t));
	else {
		memmove(short_seq->seq + borrow_len, short_seq->seq, short_seq->len);
		memcpy(short_seq->seq, tail->seq, borrow_len * sizeof(ubyte_t));
	}

	short_seq->len += borrow_len;
	free_read_seq(tail);
}

/**
 * Try to get the right connected
 */
edge *get_connected_tpl(hash_map *hm, tpl_hash *all_tpls, uint64_t query_int,
		int *locus, const int ori) {
	int next_c = 0, eg_id = 0;
	uint64_t value = 0;
	edge *right_tpl = NULL;

	next_c = next_char_by_kmers(hm, query_int, 0, ori);
	query_int = shift_bit(query_int, next_c, hm->o->k, ori);
	read_tpl_using_kmer(query_int, hm->hash, &eg_id, locus, &value);

	tpl_hash::iterator it = all_tpls->find(eg_id);
	if (it != all_tpls->end()) {
		right_tpl = (edge*) it->second;
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
int val_branching(hash_map *hm, edge *main_tpl, tpl_hash *all_tpls,
		const int shift, uint64_t query_int, int c, const int ori, int *weight,
		const int max_len) {
	uint64_t query = 0;
	int n_reads = 0, existing_shift = 0;
	bwa_seq_t *branch_seq = NULL, *main_seq = NULL;
	edge *connected_tpl = NULL;

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
		main_seq = cut_edge_tail(main_tpl, shift, max_len / 2, ori);
		n_reads = find_junc_reads(hm, branch_seq, main_seq, max_len, weight);
	} else {
		main_seq = cut_edge_tail(main_tpl, shift + hm->o->k, max_len / 2, ori);
		n_reads = find_junc_reads(hm, main_seq, branch_seq, max_len, weight);
	}
	//p_ctg_seq("MAIN", main_tpl->ctg);
	//p_ctg_seq("SUB", main_seq);
	//p_ctg_seq("BRANCH", branch_seq);
	//show_debug_msg(__func__, "N_READS: %d \n", n_reads);
	bwa_free_read_seq(1, main_seq);
	bwa_free_read_seq(1, branch_seq);
	return n_reads;
}

/**
 * If the template reaches some kmer which is used, stop the extension and add a branching event
 */
int connect(edge *branch, hash_map *hm, tpl_hash *all_tpls, uint64_t query_int,
		const int ori) {
	int *counters = NULL, locus = 0, i = 0, connected = 0, eg_id = 0,
			valid = 0, weight = 0, con_pos = 0, exist_ori = 0;
	uint64_t value = 0, query_copy = query_int;
	edge *existing = NULL;

	//show_debug_msg(__func__, "---------- Connecting to existing, ori %d ----------\n", ori);
	//p_ctg_seq(__func__, branch->ctg);
	counters = count_next_kmers(hm, query_int, 0, ori);
	//bwa_seq_t *debug = get_kmer_seq(query_int, 25);
	//p_query(__func__, debug);
	//bwa_free_read_seq(1, debug);
	//show_debug_msg(__func__, "Counters: %d,%d,%d,%d\n", counters[0],
	//		counters[1], counters[2], counters[3]);

	// In case connecting to the template itself
	mark_tpl_kmers_used(branch, hm, hm->o->k, 0);
	for (i = 0; i < 4; i++) {
		if (counters[i] < MIN_WEIGHT)
			continue;
		query_int = shift_bit(query_copy, i, hm->o->k, ori);
		read_tpl_using_kmer(query_int, hm->hash, &eg_id, &locus, &value);

		tpl_hash::iterator it = all_tpls->find(eg_id);
		if (it != all_tpls->end()) {
			existing = (edge*) it->second;
			// It happens when 'existing' and 'branch' are the same.
			// And the same template has been trimmed before.
			if (locus > existing->len - hm->o->k)
				continue;
			con_pos = ori ? (locus + 1) : (locus + hm->o->k - 1);

			/**
			 bwa_seq_t *debug = get_kmer_seq(query_int, 25);
			 p_query(__func__, debug);
			 bwa_free_read_seq(1, debug);
			 show_debug_msg(__func__, "connect pos: %d; locus: %d \n", con_pos,
			 locus);
			 **/
			valid = find_junc_reads_w_tails(hm, existing, branch, con_pos,
					(hm->o->read_len - SHORT_BRANCH_SHIFT) * 2, ori, &weight);
			if (valid) {
				exist_ori = ori ? 0 : 1;
				if (branch->len < hm->o->k - 1) {
					if (exist_ori)
						con_pos -= branch->len;
					else
						con_pos = locus + branch->len + 1;
					branch->len = 0;
					branch->ctg->len = 0;
					set_rev_com(branch->ctg);
				} else {
					// Make the branch not sharing a 24-mer with the main
					if (exist_ori) {
						con_pos -= (hm->o->k - 1);
					} else {
						con_pos += (hm->o->k - 1);
						memmove(branch->ctg->seq, branch->ctg->seq + (hm->o->k
								- 1), sizeof(ubyte_t) * (branch->len
								- (hm->o->k - 1)));
					}
					branch->len -= (hm->o->k - 1);
					branch->ctg->len = branch->len;
					set_rev_com(branch->ctg);
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
	}
	free(counters);
	return connected;
}

int existing_connect(edge *branch, hash_map *hm, tpl_hash *all_tpls,
		uint64_t query_int, int ori) {
	int connected = 0, rev_ori = 0;
	// During extension, the sequence is actually reversed, here reverse back temp
	if (ori)
		seq_reverse(branch->len, branch->ctg->seq, 0);
	set_rev_com(branch->ctg);
	connected = connect(branch, hm, all_tpls, query_int, ori);
	// Try the reverse complement of the branch and connect
	//	, if there is no other template connecting to it currently
	if (!connected && !branch->in_connect) {
		switch_fr(branch->ctg);
		rev_ori = ori ? 0 : 1;
		query_int = rev_comp_kmer(query_int, hm->o->k);
		connected = connect(branch, hm, all_tpls, query_int, rev_ori);
		// If connected, no need to reverse back, because the extending will be always stopped
		if (!connected)
			switch_fr(branch->ctg);
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
void trim_weak_tails(edge *eg, hash_map *hm, const int ori) {
	uint64_t query_int = 0;
	int i = 0, support = 0, kmer_len = hm->o->k, shift = 0, sum_sup = 0;
	if (ori) {
		for (i = 0; i <= eg->len - kmer_len; i++) {
			query_int = get_kmer_int(eg->ctg->seq, i, 1, kmer_len);
			support = get_kmer_count(query_int, hm, 0);
			if (support <= 1) {
				shift++;
				sum_sup += support;
			} else
				break;
		}
		if (shift > 0) {
			//p_ctg_seq(__func__, eg->ctg);
			memcpy(eg->ctg->seq, eg->ctg->seq + shift, sizeof(ubyte_t)
					* (eg->len - shift));
			eg->ctg->len = eg->len - shift;
			eg->len = eg->ctg->len;
			eg->kmer_freq -= sum_sup;
			set_rev_com(eg->ctg);
			//p_ctg_seq(__func__, eg->ctg);
		}
	} else {
		//p_ctg_seq(__func__, eg->ctg);
		for (i = eg->len - kmer_len; i >= 0; i--) {
			query_int = get_kmer_int(eg->ctg->seq, i, 1, kmer_len);
			support = get_kmer_count(query_int, hm, 0);
			if (support <= 1) {
				sum_sup += support;
				eg->ctg->len--;
			} else
				break;
		}
		eg->kmer_freq -= sum_sup;
		eg->len = eg->ctg->len;
		set_rev_com(eg->ctg);
		//p_ctg_seq(__func__, eg->ctg);
	}
}

/**
 * For an existing edge, try to assemble the branches
 */
void kmer_ext_branch(edge *eg, hash_map *hm, tpl_hash *all_tpls, const int ori) {
	int *counters = NULL, weight = 0, c = 0, i = 0, j = 0, con_pos = 0;
	int kmer_len = hm->o->k, branch_is_valid = 0, max_freq = 0, con_existing =
			0;
	uint64_t query_int = 0, branch_query = 0;
	edge *branch = NULL;
	if (eg->len <= kmer_len)
		return;
	//show_debug_msg(__func__, "^^^^^ Branching [%d, %d] to ori %d ^^^^^^\n",
	//		eg->id, eg->len, ori);
	for (i = 0; i < eg->len - kmer_len; i++) {
		query_int = get_kmer_int(eg->ctg->seq, i, 1, kmer_len);
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
			 eg->id, eg->len, i, ori, counters[0], counters[1],
			 counters[2], counters[3]);
			 show_debug_msg(__func__, "Counters[]: %d, Max frequency: %d \n",
			 counters[j], max_freq);
			 **/

			weight = 0;
			branch_query = shift_bit(query_int, j, kmer_len, ori);
			branch_is_valid = val_branching(hm, eg, all_tpls, i, query_int, j,
					ori, &weight, (hm->o->read_len - SHORT_BRANCH_SHIFT) * 2);
			if (!branch_is_valid || weight < MIN_WEIGHT) {
				continue;
			}
			show_debug_msg(__func__, "Weight: %d \n", weight);
			branch = blank_edge(branch_query, kmer_len, 1, ori);
			branch->kmer_freq = get_kmer_count(branch_query, hm, 1);

			con_pos = ori ? i : i + kmer_len;
			set_tail(branch, eg, con_pos, hm->o->read_len - SHORT_BRANCH_SHIFT,
					ori);
			// Insert first, in case it connects to itself during extension
			g_mutex_lock(kmer_id_mutex);
			all_tpls->insert(make_pair<int, edge*> ((int) branch->id,
					(edge*) branch));
			g_mutex_unlock(kmer_id_mutex);
			con_existing = kmer_ext_edge(branch, branch_query, hm, all_tpls,
					ori);
			cal_coverage(branch, hm);
			// If the branch can be merged into main template, erase the branch
			if (branch_on_main(eg->ctg, branch->ctg, con_pos, (branch->len
					/ hm->o->k + 2) * 3, ori)) {
				show_debug_msg(__func__, "Branch is not valid: [%d, %d] \n",
						branch->id, branch->len);
				// Remove the branch from the global hash table
				g_mutex_lock(kmer_id_mutex);
				all_tpls->erase(branch->id);
				g_mutex_unlock(kmer_id_mutex);
				// Although the branch will be removed, its kmers are marked as used by the main template.
				branch->id = eg->id;
				// The kmer frequencies are added to the main tpl now.
				eg->kmer_freq += branch->kmer_freq;
				if (ori) {
					mark_tpl_kmers_used(branch, hm, kmer_len, eg->len
							- branch->len);
				} else
					mark_tpl_kmers_used(branch, hm, kmer_len, eg->len);
				// Do not destroy first, because it may be right connected to some edge. Mark it as not alive
				// But free most of the memory it occupies, since there are many invalid branches.
				// Destory not-alive junctions and edges in function clean_junctions.
				free_eg_seq(branch);
				branch->alive = 0;
			} else {
				add_a_junction(eg, branch, query_int, con_pos, ori, weight);
				mark_tpl_kmers_used(branch, hm, kmer_len, 0);
				upd_tpl_jun_locus(branch, branching_events, kmer_len);
				// Try to extend branches of current branch
				kmer_ext_branch(branch, hm, all_tpls, 0);
				kmer_ext_branch(branch, hm, all_tpls, 1);
			}
		}
		free(counters);
		counters = NULL;
	}
	free(counters);
}

int kmer_ext_edge(edge *eg, uint64_t query_int, hash_map *hm,
		tpl_hash *all_tpls, const int ori) {
	int max_c = 0, *counters = NULL, weight = 0, con_existing = 0;
	int max_c_all = 0, *counters_all = NULL;

	show_debug_msg(__func__, "------ Started extending edge %d to ori %d...\n",
			eg->id, ori);
	bwa_seq_t *debug = get_kmer_seq(query_int, 25);
	p_query(__func__, debug);
	bwa_free_read_seq(1, debug);
	if (ori)
		seq_reverse(eg->len, eg->ctg->seq, 0);
	while (1) {
		weight = 0;
		counters = count_next_kmers(hm, query_int, 1, ori);
		max_c = get_max_index(counters);
		counters_all = count_next_kmers(hm, query_int, 0, ori);
		max_c_all = get_max_index(counters_all);

		/**
		 if (ori)
		 seq_reverse(eg->len, eg->ctg->seq, 0);
		 bwa_seq_t *debug = get_kmer_seq(query_int, 25);
		 p_query(__func__, debug);
		 bwa_free_read_seq(1, debug);
		 show_debug_msg(__func__, "Next chars: [%d, %d, %d, %d] \n",
		 counters[0], counters[1], counters[2], counters[3]);
		 show_debug_msg(__func__,
		 "Ori %d, Edge %d, length %d, Next char: %d \n", ori, eg->id,
		 eg->len, max_c);
		 p_ctg_seq("Contig", eg->ctg);
		 if (ori)
		 seq_reverse(eg->len, eg->ctg->seq, 0);
		 **/

		// If the max direction is taken already, connect to it.
		if (max_c_all != max_c) {
			if (existing_connect(eg, hm, all_tpls, query_int, ori)) {
				con_existing = 1;
				free(counters_all);
				break;
			}
		}

		if (max_c == -1) {
			if (!existing_connect(eg, hm, all_tpls, query_int, ori))
				show_debug_msg(__func__, "[%d, %d] No hits, stop here. \n",
						eg->id, eg->len);
			else {
				con_existing = 1;
			}
			free(counters);
			break;
		}
		eg->kmer_freq += counters[max_c];
		free(counters);
		free(counters_all);
		ext_con(eg->ctg, max_c, 0);
		eg->len = eg->ctg->len;
		// In case the template has repeats, so mark kmers used TEMPERARIALY. The locus is incorrect.
		mark_kmer_used(query_int, hm, eg->id, eg->len);
		query_int = shift_bit(query_int, max_c, hm->o->k, ori);
		if (eg->len % 100 == 0)
			show_debug_msg(__func__, "Ori %d, Edge %d, length %d \n", ori,
					eg->id, eg->len);
	}
	if (ori)
		seq_reverse(eg->len, eg->ctg->seq, 0);
	return con_existing;
}

/**
 * Extend a kmer
 */
void *kmer_ext_thread(gpointer data, gpointer thread_params) {
	edge *eg = NULL;
	int round_1_len = 0, round_2_len = 0, connected = 0;
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

	if (kmer_int < 64 || rev_kmer_int < 64 || c <= 1 || is_repetitive_q(query)
			|| is_biased_q(query) || kmer_is_used(kmer_int, params->hm)) {
		bwa_free_read_seq(1, query);
		return NULL;
	}
	bwa_free_read_seq(1, query);
	show_debug_msg(__func__,
			"============= %" ID64 ": %" ID64 " ============ \n", kmer_int, c);
	eg = blank_edge(kmer_int, opt->k, opt->k, 0);
	eg->kmer_freq = get_kmer_count(kmer_int, params->hm, 1);
	// Insert first, in case it connects to itself during extension
	g_mutex_lock(kmer_id_mutex);
	all_tpls->insert(make_pair<int, edge*> ((int) eg->id, (edge*) eg));
	g_mutex_unlock(kmer_id_mutex);

	connected = kmer_ext_edge(eg, kmer_int, params->hm, all_tpls, 0);
	round_1_len = eg->len;
	show_debug_msg(__func__, "Edge %d with length: %d \n", eg->id, eg->len);

	// The query int is not changed, just the direction is changed.
	connected |= kmer_ext_edge(eg, kmer_int, params->hm, all_tpls, 1);
	show_debug_msg(__func__, "Edge %d with length: %d \n", eg->id, eg->len);

	if (!connected && eg->len - round_1_len > 2) {
		kmer_int = get_kmer_int(eg->ctg->seq, eg->len - opt->k, 1, opt->k);
		round_2_len = eg->len;
		connected |= kmer_ext_edge(eg, kmer_int, params->hm, all_tpls, 0);
		show_debug_msg(__func__, "Edge %d with length: %d \n", eg->id, eg->len);
		if (eg->len - round_2_len > 2) {
			kmer_int = get_kmer_int(eg->ctg->seq, 0, 1, opt->k);
			connected |= kmer_ext_edge(eg, kmer_int, params->hm, all_tpls, 1);
			show_debug_msg(__func__, "Edge %d with length: %d \n", eg->id,
					eg->len);
		}
	}

	if (!connected && eg->len <= opt->k) {
		g_mutex_lock(kmer_id_mutex);
		all_tpls->erase(eg->id);
		g_mutex_unlock(kmer_id_mutex);
		destroy_eg(eg);
	} else {
		mark_tpl_kmers_used(eg, params->hm, opt->k, 0);
		upd_tpl_jun_locus(eg, branching_events, opt->k);
		cal_coverage(eg, params->hm);
		kmer_ext_branch(eg, params->hm, all_tpls, 0);
		kmer_ext_branch(eg, params->hm, all_tpls, 1);
		eg->start_kmer = *((uint64_t*) data);
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
		if (i % 10000 == 0)
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
	edge *eg = new_eg();
	FILE *contigs = NULL;
	GPtrArray *all_edges = NULL;
	map_opt *opt = params->hm->o;
	tpl_hash *all_tpls = params->all_tpls;

	eg = blank_edge(kmer_int, 25, 25, 0);
	all_tpls->insert(make_pair<int, edge*> ((int) eg->id, (edge*) eg));
	kmer_ext_edge(eg, kmer_int, params->hm, params->all_tpls, 0);
	round_1_len = eg->len;
	show_debug_msg(__func__, "Extending to the right: [%d, %d]... \n", eg->id,
			eg->len);
	kmer_ext_edge(eg, kmer_int, params->hm, params->all_tpls, 1);
	show_debug_msg(__func__, "Done Extending to the right: [%d, %d]... \n",
			eg->id, eg->len);

	mark_tpl_kmers_used(eg, params->hm, opt->k, 0);
	upd_tpl_jun_locus(eg, branching_events, opt->k);
	//kmer_ext_branch(eg, params->hm, all_tpls, 0);
	//kmer_ext_branch(eg, params->hm, all_tpls, 1);

	kmer_int = 502458011902035;
	eg = blank_edge(kmer_int, 25, 25, 0);
	all_tpls->insert(make_pair<int, edge*> ((int) eg->id, (edge*) eg));
	kmer_ext_edge(eg, kmer_int, params->hm, params->all_tpls, 0);
	round_1_len = eg->len;
	show_debug_msg(__func__, "Extending to the right: [%d, %d]... \n", eg->id,
			eg->len);
	kmer_ext_edge(eg, kmer_int, params->hm, params->all_tpls, 1);
	show_debug_msg(__func__, "Done Extending to the right: [%d, %d]... \n",
			eg->id, eg->len);

	mark_tpl_kmers_used(eg, params->hm, opt->k, 0);
	upd_tpl_jun_locus(eg, branching_events, opt->k);
	//	kmer_ext_branch(eg, params->hm, all_tpls, 0);
	//	kmer_ext_branch(eg, params->hm, all_tpls, 1);

	all_edges = hash_to_array(all_tpls);
	contigs = xopen(get_output_file("single.fa", kmer_out), "w");
	save_edges(all_edges, contigs, 0, 0, 0);
	fflush(contigs);
	fclose(contigs);

	clean_junctions(branching_events);
	g_ptr_array_sort(branching_events, (GCompareFunc) cmp_junctions_by_id);
	store_junctions(get_output_file("single.junctions", kmer_out),
			branching_events);
}

/**
 * Start branching after the frequent kmers are consumed already.
 */
void start_branching(tpl_hash *all_tpls, kmer_t_meta *params) {
	uint64_t id = 0;
	edge *eg = NULL;
	for (tpl_hash::iterator m = all_tpls->begin(); m != all_tpls->end(); ++m) {
		id = m->first;
		eg = (edge*) m->second;
		kmer_ext_branch(eg, params->hm, all_tpls, 0);
		kmer_ext_branch(eg, params->hm, all_tpls, 1);
	}
}

void ext_by_kmers_core(char *lib_file, const char *solid_file) {
	hash_map *hm = NULL;
	mer_hash map;
	FILE *contigs = NULL;
	kmer_t_meta *params = (kmer_t_meta*) calloc(1, sizeof(kmer_t_meta));
	tpl_hash all_tpls;
	GPtrArray *kmer_edges = NULL;

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
	kmer_edges = hash_to_array(&all_tpls);
	save_edges(kmer_edges, contigs, 0, 0, 0);
	fflush(contigs);
	fclose(contigs);

	clean_junctions(branching_events);
	g_ptr_array_sort(branching_events, (GCompareFunc) cmp_junctions_by_id);
	store_junctions(get_output_file("paired.junctions", kmer_out),
			branching_events);
	//	reset_edge_ids(all_kmer_edges);
	//	merge_ol_edges(all_kmer_edges, ins_size, sd_ins_size, hm->seqs,
	//			kmer_n_threads);
	//	free(params);
	//	destroy_hm(hm);
	//
	//	contigs = xopen(get_output_file("merged.fa", kmer_out), "w");
	//	save_edges(all_kmer_edges, contigs, 0, 0, 100);
	//	fflush(contigs);
	//	fclose(contigs);
	//	contigs = xopen(get_output_file("merged_all.fa", kmer_out), "w");
	//	save_edges(all_kmer_edges, contigs, 0, 0, 0);
	//	fflush(contigs);
	//	fclose(contigs);
}

int pe_kmer(int argc, char *argv[]) {
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
