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
#include "pool.hpp"
#include "merge.hpp"

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
		if (t->alive) {
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
	bwa_seq_t *tail = NULL, *r = NULL, *tail_shift = NULL;
	index64 main_id = 0;
	int read_len = ht->o->read_len;
	int i = 0;
	ubyte_t x = 0;
	GPtrArray *mains = g_ptr_array_sized_new(0);
	GPtrArray *hits = NULL;

	tail = get_tail(branch, ht->o->read_len, ori);

	// If the tail is like 'AAAAAAATAAAA', ignore
	if (is_biased_q(tail) || tail->len < ht->o->read_len) {
		bwa_free_read_seq(1, tail);
		return mains;
	}

	// Try ACGT four directions
	for (x = 0; x < 4; x++) {
		tail_shift = new_seq(tail, tail->len, 0);
		ext_que(tail_shift, x, ori);

		//p_query(__func__, tail_shift);

		hits = find_both_fr_full_reads(ht, tail_shift, hits, N_MISMATCHES);
		for (i = 0; i < hits->len; i++) {
			r = (bwa_seq_t*) g_ptr_array_index(hits, i);
			//p_query(__func__, r);
			//if (r->status != USED)
			//	continue;
			main_id = r->contig_id;
			tpl_hash::iterator it = all_tpls->find(main_id);
			if (it != all_tpls->end()) {
				if (r->status == USED)
					g_ptr_array_add(mains, r);
			}
		}
		bwa_free_read_seq(1, tail_shift);
	}
	//mains = rm_duplicates(mains);
	mains = rm_dup_connectors(mains);
	for (i = 0; i < mains->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(mains, i);
		//p_query(__func__, r);
	}

	bwa_free_read_seq(1, tail);
	return mains;
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
 * During extension, if it reaches some read which is used already, try to connect to it.
 */
int connect_by_full_reads(hash_table *ht, tpl_hash *all_tpls, tpl *branch,
		const int ori) {
	GPtrArray *con_reads = NULL, *junc_reads = NULL;
	index64 i = 0, j = 0;
	bwa_seq_t *r = NULL, *tail = NULL;
	tpl *main_tpl = NULL;
	// Positions, lengths, etc.
	int locus = 0, con_pos = 0, exist_ori = 0, loop_len = 0;
	int copy_start = 0, copy_end = 0;
	// Indicators, etc
	int on_main = 0, valid = 0, weight = 0, borrow_bases = 0;
	int connected = 0, is_rev = 0;
	// If the branch is reverse complement connected, the direction needs to be switch
	int adj_ori = 0;
	int max_trial = 0;
	junction *exist_junc = NULL;

	// If extending to the left, and it's not connected to any template, mark it 'dead'
	if (ori && (!branch->b_juncs || branch->b_juncs->len == 0) && branch->len
			<= 3 * ht->o->k) {
		branch->alive = 0;
		return 0;
	}

	con_reads = find_connected_reads(ht, all_tpls, branch, ori);
	//show_debug_msg(__func__, "Connecting reads: \n");
	//p_readarray(con_reads, 1);
	max_trial = con_reads->len > 8 ? 8 : con_reads->len;
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
		if (main_tpl->len <= ht->o->k) {
			continue;
		}

		// If the extension reaches some where immediately when starting,
		// 	most likely it has more mismatches which are not capture before.
		// Simply add the starting read to the main template, destory the branch.
		if (branch->len == ht->o->read_len && !branch->b_juncs) {
			locus = ori ? r->contig_locus + 1 : r->contig_locus - 1;
			add2tpl(main_tpl, branch->start_read, locus);
			branch->alive = 0;
			p_query("CONNECTOR", r);
			show_debug_msg(__func__, "Read %s added to template %d at %d \n",
					branch->start_read->name, main_tpl->id, locus);
			break;
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
		 * Its length is read_len - 1.
		 */
		tail = get_tail(branch, ht->o->read_len - 1, ori);
		//p_query("BRANCH TAIL", tail);
		//p_query("CONNECTOR", r);
		con_pos = ori ? r->contig_locus + 1 : r->contig_locus;
		if (!similar_bytes(tail->seq, main_tpl->ctg->seq + con_pos, r->len - 1,
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
		bwa_free_read_seq(1, tail);
		//p_tpl(branch);

		// Here we are sure they are going to connect,
		show_debug_msg(__func__,
				"Trying to connect [%d, %d] and [%d, %d] at %d ori %d...\n",
				main_tpl->id, main_tpl->len, branch->id, branch->len, con_pos,
				adj_ori);

		locus = r->contig_locus;
		con_pos = adj_ori ? (locus + 1) : (locus + ht->o->read_len - 1);
		exist_ori = adj_ori ? 0 : 1;

		// If right and left connections are too close, just ignore.
		if (branch->b_juncs && branch->b_juncs->len > 0) {
			exist_junc = (junction*) g_ptr_array_index(branch->b_juncs, 0);
			p_junction(exist_junc);
			if (exist_junc->main_tpl == main_tpl) {
				// If all of them simply too short, or two short same-direction junctions
				if ((get_abs(exist_junc->locus - con_pos) <= IGNORE_DIFF
						&& branch->len <= IGNORE_DIFF + ht->o->k * 3)
						|| exist_junc->ori == exist_ori) {
					branch->alive = 0;
					show_debug_msg(
							__func__,
							"Ignored the template [%d, %d] because too short\n",
							branch->id, branch->len);
					break;
				}
				// If the branch is likely to be merged to the main template, try it
				loop_len = get_abs(branch->len - (exist_junc->locus - con_pos));
				if (loop_len <= IGNORE_DIFF) {
					on_main = branch_on_main(main_tpl->ctg, branch->ctg,
							con_pos, (branch->len / ht->o->k + 2) * 3,
							exist_ori);
					valid = on_main ? 0 : 1;
					if (!valid) {
						// Mark it as 'dead', will be destroyed in kmer_ext_thread.
						show_debug_msg(
								__func__,
								"Ignore the template [%d, %d] because merged to main template \n",
								branch->id, branch->len);
						branch->alive = 0;
						break;
					}
				}
			} // Check the bubble only if left and right junction connect to the same main_tpl
		} // End of checking the short 'bubble'
		junc_reads = find_junc_reads_w_tails(ht, main_tpl, branch, con_pos,
				(ht->o->read_len - SHORT_BRANCH_SHIFT) * 2, adj_ori, &weight);
		valid = (junc_reads->len >= MIN_JUNCTION_READS) ? 1 : 0;

		//p_readarray(junc_reads, 1);

		if (!valid) {
			show_debug_msg(__func__,
					"No enough junction reads. Check mates now...\n");
			//			valid = vld_junc_by_mates(main_tpl, branch, junc_reads, ht,
			//					con_pos, adj_ori);
			//			if (!valid) {
			//				show_debug_msg(__func__,
			//						"Not passed the pair validation. Free it later.\n");
			unhold_reads_array(junc_reads);
			continue;
			//			}
		} // End of checking junction reads and pairs
		// If junction added, the junction reads should be added to branch later
		// Not only free them, but reset the them to FRESH
		unhold_reads_array(junc_reads);

		// Trim the branch
		if (branch->len < ht->o->read_len) {
			if (exist_ori)
				con_pos -= branch->len;
			else
				con_pos = locus + branch->len + 1;
			branch->len = 0;
			branch->ctg->len = 0;

		} else {
			// Make the branch not sharing a (read_length - 1) subseq with the main
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
		} // End of trimming the branch
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

		//p_tpl(branch);

		// Finally! Go to add the junction!
		show_debug_msg(__func__,
				"Connect existing [%d, %d] to [%d, %d] at %d with ori %d. \n",
				branch->id, branch->len, main_tpl->id, main_tpl->len, con_pos,
				exist_ori);
		set_tail(branch, main_tpl, con_pos, ht->o->read_len - 1, exist_ori);
		//p_ctg_seq("Right tail", branch->r_tail);
		//p_ctg_seq("Left  tail", branch->l_tail);
		add_a_junction(main_tpl, branch, 0, con_pos, exist_ori, weight);
		refresh_reads_on_tail(ht, branch, N_MISMATCHES);
		if (is_rev) {
			// Means reverse complement sequence connected
			connected = 2;
		} else {
			// Means forward sequence connected
			connected = 1;
		}
		// This branch could be connected to single another template
		break;
	} // End of connecting all probable templates
	// If it is not connected and is reverse complemented, just get it back
	if (is_rev && !connected) {
		switch_fr(branch->ctg);
	}
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
	int ext_len = 0;
	bwa_seq_t *tail = new_seq(query, query->len, 0);

	show_debug_msg(__func__,
			"------ Started extending tpl %d to ori %d... ------\n", t->id, ori);
	p_query(__func__, tail);
	//p_ctg_seq("TEMPLATE", t->ctg);

	while (1) {
		if (same_bytes(tail->seq, tail->len) || is_repetitive_q(tail)) {
			show_debug_msg(__func__, "Repetitive tail, stop.\n");
			p_query(__func__, tail);
			break;
		}

		max_c = get_next_char(ht, p, t, ori);
		// If cannot extend, try to add mates into the pool
		if (max_c == -1) {
			//p_ctg_seq("TEMPLATE", t->ctg);
			find_hashed_mates(ht, p, t, tail->len + 1, MORE_MISMATCH, ori);
			max_c = get_next_char(ht, p, t, ori);
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

		//if (t->id == 5) {
		//show_debug_msg(__func__, "Template [%d, %d], Next char: %c \n", t->id, t->len, "ACGTN"[max_c]);
		//p_query(__func__, tail);
		//p_ctg_seq("TEMPLATE", t->ctg);
		//p_pool("CURRENT POOL", p, NULL);
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
		forward(p, t, ori);
		next_pool(ht, p, t, tail, N_MISMATCHES, ori);

		if (t->len % 100 == 0)
			show_debug_msg(__func__, "Ori %d, tpl %d, length %d \n", ori,
					t->id, t->len);
		if (t->len > 50000) {
			p_tpl(t);
			err_fatal(__func__,
					"[ERROR] Too long contig [%d, %d]. Maybe some bug.\n",
					t->id, t->len);
		}
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

//	if (kmer_ctg_id == 1)
//		read = &seqs[22133];
//	if (kmer_ctg_id == 2)
//		read = &seqs[17452];

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
	// If it is connected to somewhere, simply stop
	while (iter++ <= 4 && t->len > pre_len && (!t->b_juncs || t->b_juncs->len
			== 0)) {
		// Extend to the right first
		// Make a clone of the original starting read, which is global
		p = new_pool();
		init_pool(ht, p, t, kmer_len, N_MISMATCHES, 0);
		//p_query(__func__, query);
		//g_ptr_array_sort(p->reads, (GCompareFunc) cmp_reads_by_name);
		//p_pool("INITIAL_POOL", p, NULL);
		//exit(1);

		// The correction is done only once
		if (pre_len == 0)
			correct_init_tpl_base(p, t, kmer_len);

		query = get_tail(t, kmer_len, 0);

		connected = kmer_ext_tpl(ht, all_tpls, p, t, query, 0);
		destroy_pool(p);
		pre_len = t->len;
		set_rev_com(t->ctg);
		refresh_tpl_reads(ht, t, N_MISMATCHES);
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

		connected |= kmer_ext_tpl(ht, all_tpls, p, t, query, ori);
		set_rev_com(t->ctg);
		destroy_pool(p);
		bwa_free_read_seq(1, query);
		//upd_locus_on_tpl(t, pre_len, pre_n_reads);

		//p_tpl(t);
		// Still necessary because the hashing may not get all reads
		upd_locus_on_tpl(t, pre_len, pre_n_reads);
		refresh_tpl_reads(ht, t, N_MISMATCHES);
		//g_ptr_array_sort(t->reads, (GCompareFunc) cmp_reads_by_contig_locus);

		//correct_tpl_base(t, ht->o->read_len);
		//p_readarray(t->reads, 1);
	}
	p_tpl(t);
	show_debug_msg(__func__,
			"==== End of tpl %d with length: %d; reads: %d ==== \n\n", t->id,
			t->len, t->reads->len);
	// Reactive the TRIED reads to FRESH, for other starting reads
	unfrozen_tried(t);
	// Used by find_pairs
	//g_ptr_array_sort(t->reads, (GCompareFunc) cmp_reads_by_contig_locus);
	//p_readarray(t->reads, 1);
	set_rev_com(t->ctg);

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
		// Remove all reads, realign from the hash table
		// It is important because some reads may not be picked during extension
		//	due to not full sensitive hashing.
		//	If not redo aligning, there would be many false junctions
		//refresh_tpl_reads(ht, t, N_MISMATCHES);
		//g_ptr_array_sort(t->reads, (GCompareFunc) cmp_reads_by_name);
		//p_readarray(t->reads, 1);
		//upd_tpl_jun_locus(t, branching_events, opt->k);
		//correct_tpl_base(t, ht->o->read_len);
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

	show_msg(__func__, "Sorting %d initial reads... \n", starting_reads->len);
	g_ptr_array_sort(starting_reads, (GCompareFunc) cmp_kmers_by_count);
	show_msg(__func__, "Extending by reads...\n");
	thread_pool = g_thread_pool_new((GFunc) kmer_ext_thread, params, 1, TRUE,
			NULL);
	for (i = 0; i < starting_reads->len; i++) {
		if (i % 100000 == 0)
		show_msg(__func__, "Extending %" ID64 "-th read... \n", i);
		counter = (kmer_counter*) g_ptr_array_index(starting_reads, i);
		kmer_ext_thread(counter, params);
		free(counter);
		//g_thread_pool_push(thread_pool, (gpointer) counter, NULL);
		//if (kmer_ctg_id >= 3)
		//	break;
	}

	show_msg(__func__, "Counting 11-mers of remaining reads ...\n");

	low_reads = g_ptr_array_sized_new(ht->n_seqs / 10);
	// Reset not USED reads to FRESH
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
	show_msg(__func__, "Extending the remaining %d reads ...\n", low_reads->len);
	for (i = 0; i < low_reads->len / 4; i++) {
		counter = (kmer_counter*) g_ptr_array_index(low_reads, i);
		// If the read does not even share any 11-mer with others, ignore
		if (counter->count <= (ht->o->read_len - ht->o->k) * 2) {
			free(counter);
			continue;
		}
		if (i % 100000 == 0)
			show_msg(__func__, "Extending %" ID64 "-th low read... \n", i);
		kmer_ext_thread(counter, params);
		free(counter);
	}

	g_thread_pool_free(thread_pool, 0, 1);
	g_ptr_array_free(starting_reads, TRUE);
	g_ptr_array_free(low_reads, TRUE);
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
		p_tpl(t);
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
				paired = find_pairs(t->reads, mt->reads, t->id, mt->id, 0,
						mt->len, MIN_PAIRS);
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
void iter_merge(hash_table *ht, tpl_hash *all_tpls) {
	int merge_iter = 0, id = 0, merged = 1;
	tpl *t = NULL;

	// Add all mates not on current template to t->tried.
	for (tpl_hash::iterator im = all_tpls->begin(); im != all_tpls->end(); ++im) {
		id = im->first;
		t = (tpl*) im->second;
		mv_unpaired_to_tried(ht->seqs, t, kmer_ctg_id);
        g_ptr_array_sort(t->reads, (GCompareFunc) cmp_reads_by_name);
	}

	// Multiple iterations, until not merged anymore
	while (merge_iter++ < 8 && merged) {
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

void ext_by_kmers_core(char *lib_file, const char *solid_file) {
	FILE *contigs = NULL;
	kmer_t_meta *params = (kmer_t_meta*) calloc(1, sizeof(kmer_t_meta));
	tpl_hash all_tpls;
	GPtrArray *read_tpls = NULL;
	hash_table *ht = NULL;
    char *fn = NULL;

	show_msg(__func__, "Library: %s \n", lib_file);

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

	show_msg(__func__, "Merging %d templates by pairs and overlapping...\n",
			all_tpls.size());
	iter_merge(ht, &all_tpls);

	show_msg(__func__, "Saving contigs: %.2f sec\n",
			(float) (kmer_finish_time.tv_sec - kmer_start_time.tv_sec));
    fn = get_output_file("paired.fa", kmer_out);
	contigs = xopen(fn, "w");
	read_tpls = hash_to_array(&all_tpls);
	save_tpls(read_tpls, contigs, 0, 0, 0);
	fflush(contigs);
	fclose(contigs);
    free(fn);

	g_ptr_array_sort(branching_events, (GCompareFunc) cmp_junc_by_id);
	clean_junctions(branching_events);
    fn = get_output_file("paired.junctions", kmer_out);
	store_features(fn,branching_events, read_tpls);
    free(fn);

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
		g_thread_init ( NULL);
	kmer_id_mutex = g_mutex_new();

	ext_by_kmers_core(argv[optind], argv[optind + 1]);

	clock_gettime(CLOCK_MONOTONIC, &kmer_finish_time);
	show_msg(__func__, "Done: %.2f sec\n", (float) (kmer_finish_time.tv_sec
			- kmer_start_time.tv_sec));
	return 0;
}
